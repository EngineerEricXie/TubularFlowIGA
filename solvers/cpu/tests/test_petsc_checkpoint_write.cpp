#include "PetscCheckpointWrite.hpp"
#include "PetscGather.hpp"
#include <fstream>
#include <iostream>
#include <iterator>
#include <csignal>
#include <sys/resource.h>

namespace fs = std::filesystem;
namespace {

void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

std::vector<char> Bytes(const fs::path& path)
{
	std::ifstream input(path, std::ios::binary);
	Check(static_cast<bool>(input), "cannot read test state");
	return std::vector<char>((std::istreambuf_iterator<char>(input)), {});
}

void Fill(Vec state, double base, bool nonfinite)
{
	PetscInt begin = 0, end = 0;
	VecGetOwnershipRange(state, &begin, &end);
	PetscScalar* values = nullptr;
	VecGetArray(state, &values);
	for (PetscInt i = begin; i < end; ++i) values[i-begin] = base+i;
	if (nonfinite && end > begin) values[0] = std::numeric_limits<double>::quiet_NaN();
	VecRestoreArray(state, &values);
}

std::vector<char> Snapshot(Vec state)
{
	PetscInt rows = 0; VecGetLocalSize(state, &rows);
	std::vector<char> bytes(static_cast<std::size_t>(rows)*sizeof(PetscScalar));
	iga::PetscReadArray view; view.Acquire(state);
	if (!bytes.empty()) std::memcpy(bytes.data(), view.Data(), bytes.size());
	return bytes;
}

// Scope always restores both the process limit and signal disposition before
// the test reports its result. The hard resource limit is never changed.
struct FileLimit {
	struct rlimit saved{};
	struct sigaction previous{};
	bool active = false;
	void Enable(bool selected)
	{
		if (!selected) return;
		Check(::getrlimit(RLIMIT_FSIZE, &saved) == 0, "cannot get test file limit");
		struct sigaction ignore{}; ignore.sa_handler = SIG_IGN; sigemptyset(&ignore.sa_mask);
		Check(::sigaction(SIGXFSZ, &ignore, &previous) == 0, "cannot ignore test SIGXFSZ");
		auto limited = saved; limited.rlim_cur = 0;
		if (::setrlimit(RLIMIT_FSIZE, &limited)) {
			::sigaction(SIGXFSZ, &previous, nullptr);
			throw std::runtime_error("cannot set test file limit");
		}
		active = true;
	}
	~FileLimit()
	{
		if (active) { ::setrlimit(RLIMIT_FSIZE, &saved); ::sigaction(SIGXFSZ, &previous, nullptr); }
	}
};

void Run(MPI_Comm comm, const fs::path& root, PetscInt rows)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::CollectiveLocalStage(comm, "writer fixtures", [&] {
		if (rank == 0) {
			Check(fs::create_directories(root), "test directory already exists");
			fs::create_directory(root/"directory.state");
			std::ofstream(root/"parent-file") << "not a directory";
			Check(::mkfifo((root/"fifo.state").c_str(),0600) == 0, "cannot make test FIFO");
			fs::create_symlink("state",root/"symlink.state");
		}
	});
	iga::petsc_checkpoint_detail::CandidateVector source;
	iga::RequireCollectivePetscSuccess(comm,"writer fixture vector",VecCreateMPI(comm,PETSC_DECIDE,rows,&source.value));
	Fill(source.value,1.25,false);
	iga::WritePetscCheckpointVector(source.value,comm,root/"state");
	int cases = 0;
	for (int mode = 0; mode < 10; ++mode) {
		if ((mode == 2 && ranks == 1) || (mode == 8 && rows == 0)) continue;
		// Use owning rank 0 when the last rank is empty. That also exercises
		// root file preparation failure; otherwise fail a non-root data write.
		const int fault_rank = rows < ranks ? 0 : ranks-1;
		Fill(source.value,-10.5-mode,mode == 8 && rank == fault_rank);
		const auto before = Snapshot(source.value);
		std::vector<char> old_file;
		iga::CollectiveLocalStage(comm,"writer original file",[&] { if (rank == 0) old_file = Bytes(root/"state"); });
		fs::path path = root/"state";
		if (mode <= 2 && rank == ranks-1) {
			if (mode == 0) path.clear();
			if (mode == 1) path = path.string()+std::string(1,'\0');
			if (mode == 2) path = root/"different.state";
		}
		if (mode == 3) path = root/"missing"/"state";
		if (mode == 4) path = root/"parent-file"/"state";
		if (mode == 5) path = root/"directory.state";
		if (mode == 6) path = root/"fifo.state";
		if (mode == 7) path = root/"symlink.state";
		std::string message;
		{
			FileLimit limit;
			iga::CollectiveLocalStage(comm,"test file limit",[&] { limit.Enable(mode == 9 && rank == fault_rank); });
			try { iga::WritePetscCheckpointVector(source.value,comm,path); }
			catch (const std::exception& error) { message = error.what(); }
		}
		const std::string stage = mode <= 1 || mode == 8 ? "checkpoint write preparation"
			: mode == 2 ? "checkpoint write agreement"
			: mode == 9 && fault_rank != 0 ? "checkpoint local write" : "checkpoint write file preparation";
		iga::CollectiveLocalStage(comm,"writer expected rejection",[&] {
			Check(message.find(stage) != std::string::npos,"wrong or missing writer diagnostic");
			Check(Snapshot(source.value) == before,"writer mutated source state");
			if (rank == 0) {
				Check(Bytes(root/"state") == old_file,"failed writer replaced previous state");
				for (const auto& entry : fs::directory_iterator(root))
					Check(entry.path().filename().string().find(".tmp.") == std::string::npos,"writer leaked temporary file");
			}
		});
		iga::RequireCollectiveSameText(comm,"writer common diagnostic",message);
		Fill(source.value,-10.5-mode,false);
		const auto finite = Snapshot(source.value);
		iga::WritePetscCheckpointVector(source.value,comm,root/"state");
		iga::CollectiveLocalStage(comm,"writer retry state",[&] { Check(Snapshot(source.value) == finite,"successful writer mutated source"); });
		iga::petsc_checkpoint_detail::CandidateVector loaded;
		iga::RequireCollectivePetscSuccess(comm,"writer retry vector",VecDuplicate(source.value,&loaded.value));
		iga::ReadPetscCheckpointVector(loaded.value,comm,root/"state");
		iga::CollectiveLocalStage(comm,"writer retry values",[&] { Check(Snapshot(loaded.value) == finite,"writer retry has wrong coefficients"); });
		++cases;
	}
	// Compare exact output bytes against the original PETSc viewer, and load
	// our output with that viewer as an independent format compatibility check.
	std::exception_ptr compatibility_error;
	if (rank == 0) try {
		Vec serial = nullptr; PetscViewer viewer = nullptr;
		using iga::petsc_checkpoint_detail::Check;
		Check(VecCreateSeq(PETSC_COMM_SELF,rows,&serial),"legacy fixture vector");
		Fill(serial,-19.5,false);
		Check(PetscViewerBinaryOpen(PETSC_COMM_SELF,(root/"legacy").c_str(),FILE_MODE_WRITE,&viewer),"legacy fixture viewer");
		Check(VecView(serial,viewer),"legacy fixture write");
		Check(PetscViewerDestroy(&viewer),"legacy fixture close");
		if (Bytes(root/"state") != Bytes(root/"legacy")) throw std::runtime_error("writer changed binary format");
		const auto expected = Snapshot(serial);
		Check(VecSet(serial,100.0),"legacy clear");
		Check(PetscViewerBinaryOpen(PETSC_COMM_SELF,(root/"state").c_str(),FILE_MODE_READ,&viewer),"legacy load viewer");
		Check(VecLoad(serial,viewer),"legacy load");
		Check(PetscViewerDestroy(&viewer),"legacy load close");
		if (Snapshot(serial) != expected) throw std::runtime_error("legacy VecLoad coefficients differ");
		Check(VecDestroy(&serial),"legacy destroy");
	} catch (...) { compatibility_error = std::current_exception(); }
	iga::CollectiveLocalStage(comm,"writer legacy compatibility",[&] { if (compatibility_error) std::rethrow_exception(compatibility_error); });
	if (rank == 0) std::cout << "petsc_checkpoint_write ranks=" << ranks << " rows=" << rows << " cases=" << cases << " passed\n";
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Check(argc == 2 && ranks == 3,"usage: mpiexec -np 3 petsc_checkpoint_write_test FRESH_OUTPUT");
		const fs::path root(argv[1]);
		Run(PETSC_COMM_WORLD,root/"world",4);
		Run(PETSC_COMM_WORLD,root/"world-sparse",2);
		Run(PETSC_COMM_WORLD,root/"world-empty",0);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
		Run(group,root/(rank == 0 ? "self" : "pair"),3);
		Run(group,root/(rank == 0 ? "self-empty" : "pair-empty"),0);
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "petsc_checkpoint_write all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}

#include "PetscCheckpointRead.hpp"
#include "PetscGather.hpp"

#include <fstream>
#include <iostream>
#include <iterator>

namespace fs = std::filesystem;

namespace {

void Check(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void LegacyFile(const fs::path& path, PetscInt rows, bool nonfinite)
{
	using iga::petsc_checkpoint_detail::Check;
	Vec state = nullptr; PetscViewer viewer = nullptr;
	Check(VecCreateSeq(PETSC_COMM_SELF,rows,&state),"fixture vector");
	PetscScalar* values = nullptr;
	Check(VecGetArray(state,&values),"fixture array");
	for (PetscInt row = 0; row < rows; ++row) values[row] = 1.25+row;
	if (nonfinite && rows) values[rows-1] = std::numeric_limits<double>::quiet_NaN();
	Check(VecRestoreArray(state,&values),"fixture restore");
	Check(PetscViewerBinaryOpen(PETSC_COMM_SELF,path.c_str(),FILE_MODE_WRITE,&viewer),"fixture viewer");
	Check(VecView(state,viewer),"fixture VecView");
	Check(PetscViewerDestroy(&viewer),"fixture viewer destroy");
	Check(VecDestroy(&state),"fixture vector destroy");
}

void Files(const fs::path& root, PetscInt rows)
{
	Check(fs::create_directories(root),"fixture directory already exists");
	LegacyFile(root/"good.state",rows,false);
	LegacyFile(root/"nonfinite.state",rows,true);
	std::ifstream input(root/"good.state",std::ios::binary);
	const std::vector<char> good((std::istreambuf_iterator<char>(input)),{});
	Check(good.size() >= 2*sizeof(PetscInt),"fixture Vec header missing");
	for (int mode = 0; mode < 5; ++mode) {
		auto bytes = good;
		if (mode == 0) bytes.clear();
		if (mode == 1) bytes[0] ^= 1;
		if (mode == 2) bytes[2*sizeof(PetscInt)-1] ^= 1;
		if (mode == 3) bytes.pop_back();
		if (mode == 4) bytes.push_back('x');
		std::ofstream output(root/("bad"+std::to_string(mode)+".state"),std::ios::binary);
		output.write(bytes.data(),static_cast<std::streamsize>(bytes.size())); output.close();
		Check(static_cast<bool>(output),"cannot write invalid fixture");
	}
	fs::copy_file(root/"good.state",root/"alias.state");
	Check(::mkfifo((root/"fifo.state").c_str(),0600) == 0,"cannot create FIFO fixture");
	// The reader has a fixed vector layout and must not load sidecar options.
	std::ofstream info(root/"good.state.info");
	info << "-vecload_block_size 999\n-ksp_type checkpoint_sidecar_must_not_load\n";
	info.close(); Check(static_cast<bool>(info),"cannot write sidecar fixture");
}

void Run(MPI_Comm communicator, const fs::path& root, PetscInt rows)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator,&rank); MPI_Comm_size(communicator,&ranks);
	std::exception_ptr fixture_error;
	if (rank == 0) try { Files(root,rows); } catch (...) { fixture_error = std::current_exception(); }
	iga::CollectiveLocalStage(communicator,"checkpoint fixture",[&] { if (fixture_error) std::rethrow_exception(fixture_error); });
	iga::petsc_checkpoint_detail::CandidateVector target;
	iga::RequireCollectivePetscSuccess(communicator,"test create",VecCreateMPI(communicator,PETSC_DECIDE,rows,&target.value));
	iga::RequireCollectivePetscSuccess(communicator,"test fill",VecSet(target.value,-10.0));
	const auto handle = target.value;
	char option_before[1024]{}, option_after[1024]{};
	PetscOptionsGetString(nullptr,nullptr,"-ksp_type",option_before,sizeof(option_before),nullptr);
	iga::ReadPetscCheckpointVector(target.value,communicator,root/"good.state");
	PetscOptionsGetString(nullptr,nullptr,"-ksp_type",option_after,sizeof(option_after),nullptr);
	const auto expected = iga::GatherAllPetscReal(target.value,communicator,rows);
	iga::CollectiveLocalStage(communicator,"checkpoint fixture values",[&] {
		Check(std::string(option_before) == option_after,"checkpoint loaded sidecar options");
		for (PetscInt row = 0; row < rows; ++row) Check(expected[row] == 1.25+row,"wrong checkpoint coefficient");
	});
	int cases = 0;
	for (int mode = 0; mode < 12; ++mode) {
		if ((mode == 6 && rows == 0) || (mode == 9 && ranks == 1)) continue;
		fs::path path;
		if (mode == 0) path = root/"missing.state";
		if (mode >= 1 && mode <= 5) path = root/("bad"+std::to_string(mode-1)+".state");
		if (mode == 6) path = root/"nonfinite.state";
		if (mode == 7) path = root;
		if (mode == 8) path = root/"fifo.state";
		if (mode >= 9) {
			path = root/"good.state";
			if (rank == ranks-1) {
				if (mode == 9) path = root/"alias.state";
				if (mode == 10) path = path.string()+std::string(1,'\0')+"suffix";
				if (mode == 11) path.clear();
			}
		}
		std::string message;
		try { iga::ReadPetscCheckpointVector(target.value,communicator,path); }
		catch (const std::exception& error) { message = error.what(); }
		iga::CollectiveLocalStage(communicator,"checkpoint expected diagnostic",[&] {
			const auto stage = mode < 9 ? "checkpoint local read" : mode == 9 ? "checkpoint read agreement" : "checkpoint read preparation";
			Check(message.find(stage) != std::string::npos,"missing/wrong checkpoint diagnostic");
		});
		iga::RequireCollectiveSameText(communicator,"checkpoint common diagnostic",message);
		const auto after = iga::GatherAllPetscReal(target.value,communicator,rows);
		iga::CollectiveLocalStage(communicator,"checkpoint unchanged state",[&] {
			Check(after == expected && target.value == handle,"failed load changed state or handle");
		});
		iga::ReadPetscCheckpointVector(target.value,communicator,root/"good.state");
		++cases;
	}
	const auto replica = root/(rank == ranks-1 ? "alias.state" : "good.state");
	iga::ReadReplicatedPetscCheckpointVector(target.value, communicator, replica);
	int replica_cases = 0;
	for (int mode = 0; mode < 7; ++mode) {
		if (mode == 3 && rows == 0) continue;
		auto path = replica;
		const char* stage = "asset content read";
		if (mode == 0 && rank == ranks-1) path = root/"missing.state";
		if (mode == 1) {
			if (rank == ranks-1) path = root/"bad0.state";
			stage = ranks == 1 ? "checkpoint local read" : "checkpoint asset state";
		}
		if (mode == 2) { path = root/"bad1.state"; stage = "checkpoint local read"; }
		if (mode == 3) { path = root/"nonfinite.state"; stage = "checkpoint local read"; }
		if (mode == 4 && rank == ranks-1) path = root/"fifo.state";
		if (mode == 5 && rank == ranks-1) path = root;
		if (mode == 6) { path = root/"bad4.state"; stage = "checkpoint local read"; }
		std::string message;
		try { iga::ReadReplicatedPetscCheckpointVector(target.value, communicator, path); }
		catch (const std::exception& error) { message = error.what(); }
		iga::CollectiveLocalStage(communicator, "replica checkpoint diagnostic", [&] {
			Check(message.find(stage) != std::string::npos, "wrong replica failure stage");
		});
		iga::RequireCollectiveSameText(communicator, "replica common diagnostic", message);
		const auto after = iga::GatherAllPetscReal(target.value, communicator, rows);
		iga::CollectiveLocalStage(communicator, "replica unchanged state", [&] {
			Check(after == expected && target.value == handle, "failed replica read changed target");
		});
		iga::ReadReplicatedPetscCheckpointVector(target.value, communicator, replica);
		++replica_cases;
	}
	iga::RequireCollectivePetscSuccess(communicator,"test destroy",VecDestroy(&target.value));
	if (rank == 0) std::cout << "petsc_checkpoint_replica ranks=" << ranks << " rows=" << rows
		<< " cases=" << replica_cases << " passed\n";
	if (rank == 0) std::cout << "petsc_checkpoint_read ranks=" << ranks << " rows=" << rows << " cases=" << cases << " passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		Check(argc == 2 && ranks == 3,"usage: mpiexec -np 3 petsc_checkpoint_read_test FRESH_OUTPUT");
		const fs::path root(argv[1]);
		Run(PETSC_COMM_WORLD,root/"world",2);
		Run(PETSC_COMM_WORLD,root/"world-empty",0);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
		Run(group,root/(rank == 0 ? "self" : "pair"),rank == 0 ? 3 : 1);
		Run(group,root/(rank == 0 ? "self-empty" : "pair-empty"),0);
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "petsc_checkpoint_read all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}

// Exercise the exact native CLI writer, including its PETSc gather lifetime.
#define main FlowCliMain
#include "../src/iga_navier_stokes.cpp"
#undef main

#include <sys/stat.h>
#include <sys/resource.h>
#include <csignal>

namespace {

void CheckOutput(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void TestOutput(MPI_Comm communicator, const fs::path& directory, double seed)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	const auto mesh = directory/"mesh.vtk";
	iga::CollectiveLocalStage(communicator, "output test fixture", [&] {
		if (rank != 0) return;
		CheckOutput(fs::create_directories(directory), "fixture directory must be new");
		std::ofstream output(mesh);
		output << "# vtk DataFile Version 3.0\ntest\nASCII\n"
			<< "DATASET UNSTRUCTURED_GRID\nPOINTS 1 double\n0 0 0\n"
			<< "CELLS 1 2\n1 0\nCELL_TYPES 1\n1\n";
		output.close();
		CheckOutput(bool(output), "fixture write failed");
	});
	Vec state = nullptr, original = nullptr;
	// Non-root ranks own no rows, including both ranks of the world group.
	iga::RequireCollectivePetscSuccess(communicator, "output test state",
		VecCreateMPI(communicator, rank == 0 ? 4 : 0, 4, &state));
	iga::RequireCollectivePetscSuccess(communicator, "output test initialize", VecSet(state, seed));
	iga::RequireCollectivePetscSuccess(communicator, "output test copy create", VecDuplicate(state, &original));
	iga::RequireCollectivePetscSuccess(communicator, "output test copy", VecCopy(state, original));
	const auto initial_handle = state;
	for (const std::string mode : {"layout", "overflow", "velocity-directory", "pressure-directory",
		"vtu-directory", "velocity-fifo", "pressure-fifo", "vtu-fifo", "missing-mesh", "missing-hdf-writer", "file-limit"}) {
		const auto path = directory/(mode+".txt");
		const auto vtk = directory/(mode+".vtu");
		fs::path damaged;
		iga::CollectiveLocalStage(communicator, "output test damage", [&] {
			if (rank != 0) return;
			if (mode.find("directory") == std::string::npos && mode.find("fifo") == std::string::npos) return;
			damaged = mode.rfind("velocity", 0) == 0 ? path
				: mode.rfind("pressure", 0) == 0 ? fs::path(path.string()+".pressure") : vtk;
			if (mode.find("directory") != std::string::npos) fs::create_directory(damaged);
			else CheckOutput(mkfifo(damaged.c_str(), 0600) == 0, "cannot create output FIFO");
		});
		std::uint64_t nodes = 1;
		if (rank == ranks-1 && mode == "layout") nodes = 2;
		if (rank == ranks-1 && mode == "overflow") nodes = std::numeric_limits<std::uint64_t>::max();
		struct rlimit original_limit{};
		using SignalHandler = void (*)(int);
		SignalHandler original_handler = SIG_DFL;
		// Fault only the writer, after MPI and the communicator are initialized.
		iga::CollectiveLocalStage(communicator, "output test file limit", [&] {
			if (rank != 0 || mode != "file-limit") return;
			CheckOutput(getrlimit(RLIMIT_FSIZE, &original_limit) == 0, "cannot read file limit");
			original_handler = std::signal(SIGXFSZ, SIG_IGN);
			CheckOutput(original_handler != SIG_ERR, "cannot ignore file limit signal");
			auto limited = original_limit;
			limited.rlim_cur = 1;
			CheckOutput(setrlimit(RLIMIT_FSIZE, &limited) == 0, "cannot set file limit");
		});
		std::string diagnostic;
		try {
			WriteFlowOutput(state, nodes, path, mode == "missing-mesh" ? directory/"absent.vtk" : mesh,
				vtk, 0.25, rank, mode == "missing-hdf-writer" ? iga::VisualizationFormat::BezierVtkHdf
					: iga::VisualizationFormat::Vtu, nullptr);
		} catch (const std::exception& error) { diagnostic = error.what(); }
		iga::CollectiveLocalStage(communicator, "output test file limit restore", [&] {
			if (rank != 0 || mode != "file-limit") return;
			CheckOutput(setrlimit(RLIMIT_FSIZE, &original_limit) == 0, "cannot restore file limit");
			CheckOutput(std::signal(SIGXFSZ, original_handler) != SIG_ERR, "cannot restore signal handler");
			CheckOutput(fs::file_size(path) == 1, "writer did not encounter file limit");
		});
		iga::CollectiveLocalStage(communicator, "output test failure assertion", [&] {
			const std::string stage = mode == "layout" || mode == "overflow"
				? "flow output layout" : "flow field output";
			CheckOutput(diagnostic.find(stage+": rank ") == 0, "missing coordinated output failure");
			CheckOutput(state == initial_handle, "output replaced source handle");
		});
		iga::RequireCollectiveSameText(communicator, "output test common diagnostic", diagnostic);
		PetscBool same = PETSC_FALSE;
		iga::RequireCollectivePetscSuccess(communicator, "output test unchanged state", VecEqual(state, original, &same));
		iga::CollectiveLocalStage(communicator, "output test recovery preparation", [&] {
			CheckOutput(same, "output failure changed source state");
			if (rank == 0 && !damaged.empty()) fs::remove(damaged);
		});
		WriteFlowOutput(state, 1, path, mesh, vtk, 0.25, rank, iga::VisualizationFormat::Vtu, nullptr);
		iga::CollectiveLocalStage(communicator, "output test recovery values", [&] {
			if (rank != 0) return;
			std::ifstream velocity(path), pressure(path.string()+".pressure");
			double x = 0, y = 0, z = 0, p = 0;
			CheckOutput(bool(velocity >> x >> y >> z) && bool(pressure >> p), "incomplete recovery fields");
			CheckOutput(x == seed && y == seed && z == seed && p == seed, "wrong recovery values");
			CheckOutput(fs::file_size(vtk) > 0, "missing recovery VTU");
		});
	}
	iga::RequireCollectivePetscSuccess(communicator, "output test cleanup", VecDestroy(&original));
	iga::RequireCollectivePetscSuccess(communicator, "output test cleanup", VecDestroy(&state));
	if (rank == 0) std::cout << "flow_output ranks=" << ranks << " failures=11 retries=11 passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		CheckOutput(argc == 2 && ranks == 3, "require 3 ranks and a new output directory");
		TestOutput(PETSC_COMM_WORLD, fs::path(argv[1])/"world", 10);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		TestOutput(group, fs::path(argv[1])/(rank == 0 ? "self" : "pair"), rank == 0 ? 20 : 30);
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

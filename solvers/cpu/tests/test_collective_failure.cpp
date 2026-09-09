#include "CollectiveFailure.hpp"
#include "PetscReadArray.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"

// Reuse the packed fixture with a partition owner that has zero matrix rows.
#define main OwnershipFixtureMain
#include "test_parallel_ownership.cpp"
#undef main

namespace {

template <class Function>
void RequireFailure(Function&& function, MPI_Comm communicator,
	const std::string& diagnostic)
{
	int matched = 0, all_matched = 0;
	try { function(); }
	catch (const std::runtime_error& error) {
		matched = std::string(error.what()).find(diagnostic) != std::string::npos;
	}
	MPI_Allreduce(&matched, &all_matched, 1, MPI_INT, MPI_MIN, communicator);
	Require(all_matched == 1, "failure diagnostic did not reach every rank");
}

void CheckProtocol(MPI_Comm communicator)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	const int offender = ranks-1;
	const auto prefix = "fault: rank "+std::to_string(offender)+": ";
	RequireFailure([&] {
		iga::CollectiveLocalStage(communicator, "fault", [&] {
			if (rank == offender) throw std::runtime_error("local input failure");
		});
	}, communicator, prefix+"local input failure");
	RequireFailure([&] {
		iga::CollectiveLocalStage(communicator, "fault", [&] {
			if (rank == offender) throw std::bad_alloc();
		});
	}, communicator, prefix);
	RequireFailure([&] {
		iga::CollectiveLocalStage(communicator, "fault", [&] {
			if (rank == offender) throw 7;
		});
	}, communicator, prefix+"unknown exception");
	RequireFailure([&] {
		iga::CollectiveLocalStage(communicator, "multiple", [&] {
			throw std::runtime_error(std::string(4096, 'x'));
		});
	}, communicator, "multiple: rank 0: ");
	iga::CollectiveLocalStage(communicator, "healthy", [] {});
	// Exercise exact comparison across chunk boundaries and embedded NULs.
	std::string text(10003, 'a');
	text[4095] = '\0';
	iga::RequireCollectiveSameText(communicator, "same text", text);
	iga::RequireCollectiveSameText(communicator, "empty text", "");
	if (ranks > 1) {
		for (int variant = 0; variant < 3; ++variant) {
			auto changed = text;
			if (rank == offender) {
				if (variant == 0) changed[8193] = 'b';
				else if (variant == 1) changed.clear();
				else changed.push_back('b');
			}
			RequireFailure([&] {
				iga::RequireCollectiveSameText(communicator, "text agreement", changed);
			}, communicator, "text agreement: rank "+std::to_string(offender)+": input differs");
		}
	}
}

iga::SimulationConfiguration TransportConfiguration()
{
	iga::SimulationConfiguration configuration;
	configuration.fields = {{"scalar", iga::FieldKind::Scalar, 0.0}};
	configuration.time = {0.1, 1};
	iga::EquationSystemDefinition system;
	system.name = "transport";
	system.kind = iga::EquationKind::LinearTransport;
	system.unknowns = {"scalar"};
	system.terms = {{iga::TermKind::TimeDerivative, "scalar", "scalar", 1.0, ""}};
	configuration.equation_systems.push_back(system);
	iga::FieldBoundaryCondition boundary;
	boundary.field = "scalar";
	boundary.kind = iga::FieldBoundaryKind::Dirichlet;
	boundary.value = {2.0};
	configuration.boundaries.push_back({1, "fixed", {boundary}});
	return configuration;
}

void CheckProductionStages(const std::filesystem::path& path, int rank)
{
	iga::Database database(path.string());
	RequireFailure([&] {
		iga::OwnedRowAssembler invalid(database, PETSC_COMM_WORLD, rank == 1 ? 0 : 1);
	}, PETSC_COMM_WORLD, "owned-row database input: rank 1: fields must be positive");
	RequireFailure([&] {
		iga::OwnedRowAssembler invalid(database, PETSC_COMM_WORLD, rank == 1 ? 2 : 1);
	}, PETSC_COMM_WORLD, "owned-row global shape: rank 1:");
	iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 1);
	RequireFailure([&] {
		const auto pattern = iga::FieldCouplingPattern::Dense(rank == 1 ? 2 : 1);
		Mat matrix = assembler.CreateMatrix(pattern);
		MatDestroy(&matrix);
	}, PETSC_COMM_WORLD, "owned-row matrix preallocation: rank 1:");
	Vec vector = assembler.CreateVector();
	VecSet(vector, 2.0);
	RequireFailure([&] {
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "view", [&] {
			iga::PetscReadArray view;
			view.Acquire(vector);
			if (rank == 2) throw std::runtime_error("empty-rank view failure");
		});
	}, PETSC_COMM_WORLD, "view: rank 2: empty-rank view failure");
	PetscReal norm = 0.0;
	VecNorm(vector, NORM_2, &norm);
	Require(std::abs(norm-std::sqrt(8.0)) < 1e-12, "view unwinding corrupted vector");
	VecDestroy(&vector);
	const auto configuration = TransportConfiguration();
	auto system = iga::CompileLinearSystem(configuration, "transport");
	system.velocity_source = "prescribed";
	RequireFailure([&] {
		const std::vector<int> labels = rank == 1 ? std::vector<int>{1} : std::vector<int>{1, 1};
		iga::TransientTransportRuntime invalid(database, PETSC_COMM_WORLD,
			configuration, system, labels);
	}, PETSC_COMM_WORLD, "transport runtime input: rank 1: transport boundary labels");
	RequireFailure([&] {
		auto invalid_system = system;
		if (rank == 1) invalid_system.terms.at(0).trial = system.fields.size();
		iga::TransientTransportRuntime invalid(database, PETSC_COMM_WORLD,
			configuration, std::move(invalid_system), {1, 1});
	}, PETSC_COMM_WORLD, "transport runtime input: rank 1: field-coupling index");
	iga::TransientTransportRuntime runtime(database, PETSC_COMM_WORLD,
		configuration, system, {1, 1});
	const auto initial = runtime.GatherState();
	std::vector<std::array<double, 3>> velocity(runtime.RequiredNodes().size(), {0.0, 0.0, 0.0});
	auto invalid_velocity = velocity;
	if (rank == 1) invalid_velocity.clear();
	runtime.BeginStep();
	RequireFailure([&] {
		runtime.SolveTrial(configuration, runtime.RequiredNodes(), invalid_velocity);
	}, PETSC_COMM_WORLD, "transport trial input: rank 1:");
	runtime.AbortStep();
	Require(runtime.Steps() == 0 && runtime.GatherState() == initial,
		"failed input changed committed transport state");
	runtime.Advance(configuration, runtime.RequiredNodes(), velocity);
	const auto retried = runtime.GatherState();
	double error = 0.0, scale = 0.0;
	for (std::size_t i = 0; i < initial.size(); ++i) {
		error = std::hypot(error, retried[i]-initial[i]);
		scale = std::hypot(scale, initial[i]);
	}
	// The established CPU comparison gate applies to a new solve; rollback
	// above is a VecCopy and must preserve the committed values exactly.
	Require(runtime.Steps() == 1 && error/scale <= 1e-6,
		"retry did not preserve constrained scalar solution");
	if (rank == 0) std::cout << "transport retry relative_l2=" << error/scale << '\n';
	if (rank == 1) system.dt = 0.0;
	iga::TransientTransportRuntime bad_assembly(database, PETSC_COMM_WORLD,
		configuration, system, {1, 1});
	bad_assembly.BeginStep();
	RequireFailure([&] {
		bad_assembly.SolveTrial(configuration, bad_assembly.RequiredNodes(), velocity);
	}, PETSC_COMM_WORLD, "transport element assembly: rank 1: linear transport time step must be positive");
	RequireFailure([&] { bad_assembly.PrepareCommitStep(); }, PETSC_COMM_WORLD,
		"requires a successful trial solve");
	bad_assembly.AbortStep();
	Require(bad_assembly.Steps() == 0 && bad_assembly.GatherState() == initial,
		"failed assembly changed committed transport state");
	iga::ResolvedBoundaryConditions boundaries;
	boundaries.velocity_constrained.assign(2, 1);
	boundaries.pressure_constrained.assign(2, 1);
	boundaries.velocity.resize(2, {0.0, 0.0, 0.0});
	boundaries.pressure.assign(2, 0.0);
	RequireFailure([&] {
		auto malformed = boundaries;
		if (rank == 1) malformed.pressure_constrained.clear();
		iga::TransientFlowRuntime invalid(database, PETSC_COMM_WORLD, false, true,
			{1.0, 0.01, 0.1}, std::move(malformed), {1, 1},
			std::vector<std::array<double, 3>>(2, {0.0, 0.0, 0.0}), {}, {});
	}, PETSC_COMM_WORLD, "flow runtime input: rank 1: flow boundary data");
	iga::TransientFlowRuntime bad_flow(database, PETSC_COMM_WORLD, false, true,
		{rank == 1 ? -1.0 : 1.0, 0.01, 0.1}, boundaries, {1, 1},
		std::vector<std::array<double, 3>>(2, {0.0, 0.0, 0.0}), {}, {});
	bad_flow.BeginStep(0, 0.1, 2, 1e-8, 1e-12, 1e-6);
	RequireFailure([&] { bad_flow.SolveTrial(); }, PETSC_COMM_WORLD,
		"flow element assembly: rank 1: invalid Navier-Stokes density");
	RequireFailure([&] { bad_flow.PrepareCommitStep(); }, PETSC_COMM_WORLD,
		"requires a successful 3D flow trial solve");
	bad_flow.AbortStep();
	VecNorm(bad_flow.State(), NORM_2, &norm);
	Require(norm == 0.0, "failed flow assembly changed committed state");
	// Exercise both history and current array acquisition again after rollback.
	bad_flow.BeginStep(0, 0.1, 2, 1e-8, 1e-12, 1e-6);
	RequireFailure([&] { bad_flow.SolveTrial(); }, PETSC_COMM_WORLD,
		"flow element assembly: rank 1: invalid Navier-Stokes density");
	bad_flow.AbortStep();
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(ranks == 3 && argc == 2, "usage: mpiexec -np 3 collective_failure_test OUTPUT_DIRECTORY");
		CheckProtocol(PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckProtocol(group);
		MPI_Comm_free(&group);
		const auto directory = std::filesystem::path(argv[1]);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "test fixture", [&] {
			if (rank != 0) return;
			Require(std::filesystem::create_directories(directory), "test directory already exists");
			WriteAuditDatabase(directory/"audit.ntiga", false);
		});
		CheckProductionStages(directory/"audit.ntiga", rank);
		std::cout << "collective failure rank=" << rank << " passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	return 0;
}

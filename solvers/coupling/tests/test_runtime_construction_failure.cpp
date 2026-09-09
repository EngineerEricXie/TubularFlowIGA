#include "BoundarySupport.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"
#include "CouplingFixture.hpp"
#include <memory>
#include <type_traits>

static_assert(!std::is_copy_constructible<iga::TransientFlowRuntime>::value,
	"flow PETSc ownership cannot be copied");
static_assert(!std::is_copy_assignable<iga::TransientFlowRuntime>::value,
	"flow PETSc ownership cannot be assigned");
static_assert(!std::is_copy_constructible<iga::TransientTransportRuntime>::value,
	"transport PETSc ownership cannot be copied");
static_assert(!std::is_copy_assignable<iga::TransientTransportRuntime>::value,
	"transport PETSc ownership cannot be assigned");

#define main ConstructionFixtureMain
#include "test_bifurcation_coupling_smoke.cpp"
#undef main

namespace {
void RequireConstruction(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void PetscConstructionCheck(PetscErrorCode code)
{
	if (code) throw std::runtime_error("construction test PETSc operation failed");
}

struct RetainedObjects {
	std::array<PetscObject, 32> objects{};
	std::size_t size = 0;

	void Retain(PetscObject object)
	{
		RequireConstruction(object && size < objects.size(), "invalid construction object observation");
		for (std::size_t i = 0; i < size; ++i)
			RequireConstruction(objects[i] != object, "construction object observed twice");
		PetscConstructionCheck(PetscObjectReference(object));
		objects[size++] = object;
	}

	void VerifyAndRelease(MPI_Comm comm)
	{
		iga::RequireCollectiveSameInt(comm, "construction retained object count", static_cast<int>(size));
		// Release retained parents (KSP/scatter) before their referenced data.
		while (size) {
			const auto object = objects[size - 1];
			iga::CollectiveLocalStage(comm, "construction object ownership", [&] {
				PetscInt references = 0;
				PetscConstructionCheck(PetscObjectGetReference(object, &references));
				RequireConstruction(references == 1, "failed constructor retained a PETSc ownership reference");
			});
			iga::RequireCollectivePetscSuccess(comm, "construction retained object release", PetscObjectDereference(object));
			objects[--size] = nullptr;
		}
	}

	~RetainedObjects()
	{
		while (size) PetscObjectDereference(objects[--size]);
	}
};

struct ProbeReset {
	~ProbeReset() { iga::RuntimeConstructionHooksForTesting() = {}; }
};

template<class Work>
void RejectConstruction(MPI_Comm comm, const char* expected, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm, "construction failure diagnostic", [&] {
		if (message.find(expected) == std::string::npos)
			throw std::runtime_error("missing/wrong construction failure: " + message);
	});
	iga::RequireCollectiveSameText(comm, "construction common diagnostic", message);
}

std::vector<double> GatherConstructionState(Vec vector)
{
	Vec all = nullptr;
	VecScatter scatter = nullptr;
	PetscConstructionCheck(VecScatterCreateToAll(vector, &scatter, &all));
	PetscConstructionCheck(VecScatterBegin(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
	PetscConstructionCheck(VecScatterEnd(scatter, vector, all, INSERT_VALUES, SCATTER_FORWARD));
	PetscInt count = 0;
	PetscConstructionCheck(VecGetSize(all, &count));
	iga::PetscReadArray view;
	view.Acquire(all);
	std::vector<double> values(static_cast<std::size_t>(count));
	for (PetscInt i = 0; i < count; ++i) values[i] = PetscRealPart(view.Data()[i]);
	view.Restore();
	PetscConstructionCheck(VecScatterDestroy(&scatter));
	PetscConstructionCheck(VecDestroy(&all));
	return values;
}

void CompareConstructionState(MPI_Comm comm, const std::vector<double>& reference,
	const std::vector<double>& actual, std::size_t stride, std::size_t component, bool invert = false)
{
	iga::CollectiveLocalStage(comm, "construction retry field", [&] {
		RequireConstruction(reference.size() == actual.size(), "construction retry field size");
		double error = 0, norm = 0;
		for (std::size_t i = 0; i < actual.size(); ++i) if ((i % stride == component) != invert) {
			error = std::hypot(error, actual[i] - reference[i]);
			norm = std::hypot(norm, reference[i]);
		}
		RequireConstruction(std::isfinite(error) && std::isfinite(norm)
			&& (norm > 0 ? error / norm <= 1e-6 : error <= 1e-12), "construction retry field differs");
	});
}

struct ConstructionInput {
	iga::SimulationConfiguration flow, scalar;
	std::vector<int> labels;
	std::vector<std::array<double, 3>> velocity;
	iga::ResolvedBoundaryConditions boundaries;
	std::set<std::int32_t> wall;
	iga::CompiledLinearSystem compiled;
};

std::unique_ptr<iga::TransientFlowRuntime> MakeConstructionFlow(iga::Database& database,
	MPI_Comm comm, const ConstructionInput& input)
{
	return iga::AllocateCollectiveRuntime<iga::TransientFlowRuntime>(comm, database, comm, true, true,
		iga::NavierStokesParameters{kDensity, kViscosity, kDt}, input.boundaries,
		input.labels, input.velocity, input.wall, std::vector<iga::OutletModelState>{});
}

std::unique_ptr<iga::TransientTransportRuntime> MakeConstructionTransport(iga::Database& database,
	MPI_Comm comm, const ConstructionInput& input)
{
	return iga::AllocateCollectiveRuntime<iga::TransientTransportRuntime>(comm, database, comm, input.scalar,
		input.compiled, input.labels);
}

std::vector<double> HealthyConstructionFlow(iga::Database& database, MPI_Comm comm, const ConstructionInput& input)
{
	auto flow = MakeConstructionFlow(database, comm, input);
	flow->InitializeState(input.flow);
	flow->Advance(input.flow, 0, kDt, 12, 1e-8, 1e-12, 1e-6);
	return GatherConstructionState(flow->State());
}

std::vector<double> HealthyConstructionTransport(iga::Database& database, MPI_Comm comm, const ConstructionInput& input)
{
	auto transport = MakeConstructionTransport(database, comm, input);
	std::vector<std::array<double, 3>> zero_velocity(transport->RequiredNodes().size());
	transport->Advance(input.scalar, transport->RequiredNodes(), zero_velocity);
	const auto values = transport->GatherState();
	iga::CollectiveLocalStage(comm, "construction scalar analytic state", [&] {
		for (double value : values)
			RequireConstruction(std::abs(value - (2.0 + kDt)) <= 1e-10, "constant scalar source solution differs");
	});
	return values;
}

struct ConstructionFailureCase {
	const char* stage;
	int occurrence = 1;
	const char* diagnostic = nullptr;
	bool allocation_failure = false;
};

// One constant field basis on a regular Bezier cube. This is a construction
// fixture, not a Navier--Stokes discretization: trailing ranks own zero rows
// and have an empty required-element/halo list. Numerical retries below use
// the separate 64-basis fixture.
void WriteEmptyRankDatabase(const fs::path& path, std::uint32_t ranks)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	RequireConstruction(static_cast<bool>(output), "cannot create empty-rank database");
	constexpr std::uint64_t header_size = 88;
	constexpr std::uint64_t element_offset = header_size + 2*sizeof(std::uint64_t) + sizeof(std::int32_t);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, ranks);
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto rank_index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	for (int axis = 0; axis < 3; ++axis) iga::Write(output, 0.0);
	iga::Write(output, 1.0);
	iga::Write(output, 1.0);
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint32_t{1});
	for (int face = 0; face < 6; ++face) iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::uint8_t{64});
	for (std::uint8_t column = 0; column < 64; ++column) {
		iga::Write(output, column);
		iga::Write(output, 1.0);
	}
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i)
				for (double value : {i/3.0, j/3.0, k/3.0}) iga::Write(output, value);
	const auto rank_index_offset = static_cast<std::uint64_t>(output.tellp());
	iga::Write(output, std::uint64_t{0});
	for (std::uint32_t rank = 0; rank < ranks; ++rank) iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{0});
	output.seekp(header_size + sizeof(std::uint64_t));
	iga::Write(output, rank_index_offset);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	RequireConstruction(static_cast<bool>(output), "cannot finalize empty-rank database");
}

template<class Build, class Retry>
int ConstructionFailures(MPI_Comm comm, const std::vector<ConstructionFailureCase>& cases,
	Build&& build, Retry&& retry)
{
	int rank = 0, ranks = 1, count = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	for (const auto& failure : cases) {
		RetainedObjects objects;
		ProbeReset reset;
		int occurrences = 0;
		auto& hooks = iga::RuntimeConstructionHooksForTesting();
		hooks.object_created = [&](PetscObject object) { objects.Retain(object); };
		hooks.after_local_stage = [&](const char* stage) {
			if (std::strcmp(stage, failure.stage) == 0 && ++occurrences == failure.occurrence && rank == ranks - 1) {
				if (failure.allocation_failure) throw std::bad_alloc();
				throw std::runtime_error(std::string("injected construction failure: ") + stage);
			}
		};
		RejectConstruction(comm, failure.diagnostic ? failure.diagnostic : failure.stage, build);
		if (failure.allocation_failure)
			iga::CollectiveLocalStage(comm, "construction early failure owns no PETSc objects", [&] {
				RequireConstruction(objects.size == 0, "early allocation failure reached PETSc object creation");
			});
		hooks = {};
		objects.VerifyAndRelease(comm);
		retry();
		++count;
		if (rank == 0) std::cout << "construction_failure ranks=" << ranks << " stage=" << failure.stage
			<< " occurrence=" << failure.occurrence << " cleanup_retry=passed\n";
	}
	for (bool mismatched : {false, true}) {
		if (mismatched && ranks == 1) continue;
		RetainedObjects objects;
		ProbeReset reset;
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { objects.Retain(object); };
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_type",
			mismatched ? (rank == ranks - 1 ? "gmres" : "fgmres") : "unavailable_construction_solver"));
		RejectConstruction(comm, mismatched ? "PETSc option agreement" : "solver options", build);
		iga::RuntimeConstructionHooksForTesting() = {};
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_type", "fgmres"));
		objects.VerifyAndRelease(comm);
		retry();
		++count;
	}
	return count;
}

void CheckConstructionOption(const char* name, const char* expected)
{
	char value[64]{};
	PetscBool present = PETSC_FALSE;
	PetscConstructionCheck(PetscOptionsGetString(nullptr, nullptr, name, value, sizeof(value), &present));
	RequireConstruction(present && std::strcmp(value, expected) == 0, "construction changed a default or user option");
}

void RunScalableConstruction(MPI_Comm comm, const fs::path& root)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::CollectiveLocalStage(comm, "construction scalable fixture", [&] {
		if (rank == 0) {
			iga::test::WriteC0SquareDuctDatabase(root / "scalable.ntiga", 3, ranks, 3);
			iga::test::WriteC0SquareDuctControlMesh(root / "scalable.vtk", 3, 3);
		}
	});
	iga::Database database((root / "scalable.ntiga").string());
	ConstructionInput input;
	iga::CollectiveLocalStage(comm, "construction scalable input", [&] {
		RequireConstruction(database.header().nodes == 1000, "fixture does not exercise the scalable default");
		const auto mesh = iga::ReadLabeledHexMesh((root / "scalable.vtk").string(), 1000, 27);
		input.labels = mesh.labels;
		input.velocity.resize(1000);
		input.boundaries.velocity_constrained.resize(1000);
		input.boundaries.pressure_constrained.resize(1000);
		input.boundaries.velocity.resize(1000);
		input.boundaries.pressure.resize(1000);
		input.wall = iga::WallTraceBasis(database, mesh);
	});
	PetscConstructionCheck(PetscOptionsClearValue(nullptr, "-pc_type"));
	// This entry belongs to the user, so rollback must leave it untouched.
	PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-fieldsplit_0_pc_type", "jacobi"));
	const char* inserted[] = {"-fieldsplit_0_ksp_type", "-fieldsplit_1_ksp_type", "-fieldsplit_1_pc_type"};
	for (const char* name : inserted) PetscConstructionCheck(PetscOptionsClearValue(nullptr, name));
	int failures = 0;
	for (const char* stage : {"-fieldsplit_0_ksp_type", "-fieldsplit_1_ksp_type",
		"flow default solver options", "flow runtime ready", "flow solver options"}) {
		const bool invalid_solver = std::strcmp(stage, "flow solver options") == 0;
		if (invalid_solver) PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_type", "unavailable_construction_solver"));
		const auto before = iga::CapturePetscOptions(nullptr);
		RetainedObjects objects;
		ProbeReset reset;
		auto& hooks = iga::RuntimeConstructionHooksForTesting();
		hooks.object_created = [&](PetscObject object) { objects.Retain(object); };
		hooks.after_local_stage = [&](const char* current) {
			if (!invalid_solver && rank == ranks - 1 && std::strcmp(stage, current) == 0)
				throw std::runtime_error(std::string("injected scalable construction failure: ") + stage);
		};
		RejectConstruction(comm, stage, [&] { auto flow = MakeConstructionFlow(database, comm, input); });
		hooks = {};
		objects.VerifyAndRelease(comm);
		iga::CollectiveLocalStage(comm, "construction default rollback", [&] {
			RequireConstruction(iga::CapturePetscOptions(nullptr) == before, "failed construction leaked PETSc defaults");
			CheckConstructionOption("-fieldsplit_0_pc_type", "jacobi");
		});
		if (invalid_solver) PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_type", "fgmres"));
		{
			auto flow = MakeConstructionFlow(database, comm, input);
			flow->Assembler().ValidateOwnership();
			CompareConstructionState(comm, std::vector<double>(4000, 0.0), GatherConstructionState(flow->State()), 1, 0);
		}
		iga::CollectiveLocalStage(comm, "construction successful defaults", [&] {
			CheckConstructionOption("-fieldsplit_0_ksp_type", "preonly");
			CheckConstructionOption("-fieldsplit_0_pc_type", "jacobi");
			CheckConstructionOption("-fieldsplit_1_ksp_type", "preonly");
			CheckConstructionOption("-fieldsplit_1_pc_type", "gamg");
		});
		for (const char* name : inserted) PetscConstructionCheck(PetscOptionsClearValue(nullptr, name));
		++failures;
		if (rank == 0) std::cout << "scalable_construction ranks=" << ranks << " stage=" << stage
			<< " default_rollback=passed user_override=preserved retry=passed\n";
	}
	PetscConstructionCheck(PetscOptionsClearValue(nullptr, "-fieldsplit_0_pc_type"));
	PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-pc_type", "lu"));
	if (rank == 0) std::cout << "scalable_construction_summary ranks=" << ranks << " cases=" << failures << '\n';
}

void RunConstructionGroup(MPI_Comm comm, const fs::path& root)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::CollectiveLocalStage(comm, "construction fixture", [&] {
		if (rank == 0) {
			RequireConstruction(fs::create_directories(root), "construction fixture already exists");
			WriteThreeDCase(root); WriteDatabase(root / "group.ntiga", ranks); WriteDatabase(root / "serial.ntiga", 1);
			WriteEmptyRankDatabase(root / "empty.ntiga", ranks);
		}
	});
	iga::Database database((root / "group.ntiga").string()), serial((root / "serial.ntiga").string());
	ConstructionInput input;
	iga::CollectiveLocalStage(comm, "construction fixture input", [&] {
		input.flow = iga::ReadSimulationConfiguration((root / "simulation_config.json").string());
		const auto mesh = iga::ReadLabeledHexMesh((root / "controlmesh.vtk").string(), database.header().nodes, database.header().elements);
		input.labels = mesh.labels;
		input.velocity = iga::ReadVelocity((root / "initial_velocityfield.txt").string(), database.header().nodes);
		input.boundaries = iga::ResolveFlowBoundaries(input.flow, iga::FirstNavierStokesSystem(input.flow), input.labels, input.velocity);
		input.wall = iga::WallTraceBasis(database, mesh);
		input.scalar.fields = {{"tracer", iga::FieldKind::Scalar, 2.0}};
		input.scalar.time = {kDt, 1};
		iga::EquationSystemDefinition system;
		system.name = "transport"; system.kind = iga::EquationKind::LinearTransport; system.unknowns = {"tracer"};
		system.terms = {{iga::TermKind::TimeDerivative, "tracer", "tracer", 1.0, ""},
			{iga::TermKind::VolumeSource, "tracer", "tracer", 1.0, ""}};
		input.scalar.equation_systems.push_back(system);
		input.scalar.velocity_sources.push_back({"prescribed", "prescribed", "", "", "error"});
		for (int label = 0; label < 4; ++label) {
			iga::FieldBoundaryCondition condition;
			condition.field = "tracer"; condition.kind = iga::FieldBoundaryKind::NoFlux;
			input.scalar.boundaries.push_back({label, "boundary" + std::to_string(label), {condition}});
		}
		input.compiled = iga::CompileLinearSystem(input.scalar, "transport");
	});
	const auto reference_flow = HealthyConstructionFlow(serial, PETSC_COMM_SELF, input);
	const auto reference_scalar = HealthyConstructionTransport(serial, PETSC_COMM_SELF, input);
	const auto retry_flow = [&] {
		const auto state = HealthyConstructionFlow(database, comm, input);
		CompareConstructionState(comm, reference_flow, state, 4, 3);
		CompareConstructionState(comm, reference_flow, state, 4, 3, true);
	};
	const auto retry_scalar = [&] {
		CompareConstructionState(comm, reference_scalar, HealthyConstructionTransport(database, comm, input), 1, 0);
	};
	if (ranks > 1) {
		iga::OwnedRowAssembler assembler(database, comm, 4);
		RejectConstruction(comm, "owned-row matrix pattern policy", [&] {
			Mat matrix = assembler.CreateMatrix(rank == ranks - 1);
			PetscConstructionCheck(MatDestroy(&matrix));
		});
		retry_flow();
		if (rank == 0) std::cout << "construction_matrix_policy ranks=" << ranks << " common_rejection=passed retry=passed\n";
	}
	const std::vector<ConstructionFailureCase> flow_failures{
		{"runtime object allocation", 1, "runtime object storage ready", true},
		{"runtime object storage ready", 1, nullptr, true},
		{"flow boundary input copied", 1, "flow runtime input", true},
		{"flow trial configuration preparation", 1, nullptr, true},
		{"flow runtime input"}, {"flow boundary catalog preparation"}, {"flow boundary catalog layout"},
		{"flow boundary catalog publication"}, {"owned-row matrix created"}, {"owned-row vector created", 1},
		{"owned-row vector created", 3}, {"owned-row vector created", 5}, {"flow boundary row preparation"},
		{"flow halo preparation"}, {"flow halo source index created"}, {"flow halo state created"},
		{"flow halo previous created"}, {"flow halo destination index created"}, {"flow halo scatter created"},
		{"flow solver created"}, {"flow runtime ready"}};
	const std::vector<ConstructionFailureCase> scalar_failures{
		{"runtime object allocation", 1, "runtime object storage ready", true},
		{"runtime object storage ready", 1, nullptr, true},
		{"transport compiled input preparation", 1, nullptr, true},
		{"transport configuration preparation", 1, nullptr, true},
		{"transport runtime input"}, {"owned-row matrix created", 1}, {"owned-row matrix created", 2},
		{"owned-row vector created", 1}, {"owned-row vector created", 3}, {"owned-row vector created", 5},
		{"transport initial array active"}, {"transport initial values"}, {"transport halo preparation"},
		{"transport halo source index created"}, {"transport halo state created"},
		{"transport halo destination index created"}, {"transport halo scatter created"},
		{"transport solver created"}, {"transport runtime ready"}};
	const int flow_cases = ConstructionFailures(comm, flow_failures,
		[&] { auto flow = MakeConstructionFlow(database, comm, input); }, retry_flow);
	const int scalar_cases = ConstructionFailures(comm, scalar_failures,
		[&] { auto transport = MakeConstructionTransport(database, comm, input); }, retry_scalar);
	iga::Database empty_database((root / "empty.ntiga").string());
	ConstructionInput empty_input;
	iga::CollectiveLocalStage(comm, "construction empty-rank input", [&] {
		empty_input = input;
		empty_input.labels = {0};
		empty_input.scalar.boundaries.erase(std::remove_if(empty_input.scalar.boundaries.begin(),
			empty_input.scalar.boundaries.end(), [](const auto& boundary) {
				return boundary.label != 0;
			}), empty_input.scalar.boundaries.end());
		empty_input.velocity.resize(1);
		empty_input.wall.clear();
		empty_input.boundaries = {};
		empty_input.boundaries.velocity_constrained.resize(1);
		empty_input.boundaries.pressure_constrained.resize(1);
		empty_input.boundaries.velocity.resize(1);
		empty_input.boundaries.pressure.resize(1);
		const auto range = empty_database.NodeRange(rank);
		RequireConstruction(range.second - range.first == (rank == 0 ? 1u : 0u), "incorrect empty-rank fixture ownership");
		RequireConstruction(empty_database.RequiredElementIndices(rank).size() == (rank == 0 ? 1u : 0u), "incorrect empty-rank fixture halo");
	});
	const auto retry_empty_flow = [&] {
		auto flow = MakeConstructionFlow(empty_database, comm, empty_input);
		CompareConstructionState(comm, std::vector<double>(4, 0.0), GatherConstructionState(flow->State()), 1, 0);
	};
	const auto retry_empty_scalar = [&] {
		auto transport = MakeConstructionTransport(empty_database, comm, empty_input);
		CompareConstructionState(comm, {2.0}, transport->GatherState(), 1, 0);
	};
	const int empty_flow_cases = ConstructionFailures(comm, flow_failures,
		[&] { auto flow = MakeConstructionFlow(empty_database, comm, empty_input); }, retry_empty_flow);
	const int empty_scalar_cases = ConstructionFailures(comm, scalar_failures,
		[&] { auto transport = MakeConstructionTransport(empty_database, comm, empty_input); }, retry_empty_scalar);
	if (rank == 0) std::cout << "runtime_construction ranks=" << ranks << " flow_cases=" << flow_cases
		<< " scalar_cases=" << scalar_cases << " empty_flow_cases=" << empty_flow_cases
		<< " empty_scalar_cases=" << empty_scalar_cases
		<< " retained_references_released=1 retry_fields=passed empty_row_construction=passed\n";
	RunScalableConstruction(comm, root);
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
	int rank = 0, ranks = 1, status = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	std::cout << std::unitbuf;
	try {
		RequireConstruction(argc == 2 && ranks == 3, "require three ranks and fresh output parent");
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_type", "fgmres"));
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-pc_type", "lu"));
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", "mumps"));
		PetscConstructionCheck(PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-12"));
		RunConstructionGroup(PETSC_COMM_WORLD, fs::path(argv[1]) / "world");
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		RunConstructionGroup(group, fs::path(argv[1]) / (rank == 0 ? "single" : "pair"));
		MPI_Comm_free(&group);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	PetscPopErrorHandler(); PetscFinalize();
	return status;
}

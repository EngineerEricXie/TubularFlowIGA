#include "CheckedText.hpp"
#include "ExecutionResources.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include "PetscGather.hpp"
#include "PetscReadArray.hpp"
#include <memory>
#include <optional>
#include "CaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"
#include "TransportElement.hpp"

#include <petscksp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

namespace {
class ReturnErrors {
public:
	void Enable(MPI_Comm communicator)
	{
		iga::CollectiveLocalStage(communicator, "legacy transport error handler", [&] {
			if (PetscPushErrorHandler(PetscReturnErrorHandler, nullptr))
				throw std::runtime_error("cannot install returning PETSc error handler");
			active_ = true;
		});
	}
	~ReturnErrors() { if (active_) PetscPopErrorHandler(); }
private:
	bool active_ = false;
};

void CheckPetsc(MPI_Comm communicator, const char* operation, PetscErrorCode status)
{
	iga::RequireCollectivePetscSuccess(communicator, operation, status);
}

struct TransportObjects {
	TransportObjects() = default;
	TransportObjects(const TransportObjects&) = delete;
	TransportObjects& operator=(const TransportObjects&) = delete;
	Mat left = nullptr, previous = nullptr;
	Vec current = nullptr, next = nullptr, rhs = nullptr;
	KSP solver = nullptr;
	~TransportObjects()
	{
		if (solver) KSPDestroy(&solver);
		if (rhs) VecDestroy(&rhs);
		if (next) VecDestroy(&next);
		if (current) VecDestroy(&current);
		if (previous) MatDestroy(&previous);
		if (left) MatDestroy(&left);
	}
	void Close(MPI_Comm communicator)
	{
		CheckPetsc(communicator, "legacy transport KSPDestroy", KSPDestroy(&solver));
		CheckPetsc(communicator, "legacy transport VecDestroy rhs", VecDestroy(&rhs));
		CheckPetsc(communicator, "legacy transport VecDestroy next", VecDestroy(&next));
		CheckPetsc(communicator, "legacy transport VecDestroy current", VecDestroy(&current));
		CheckPetsc(communicator, "legacy transport MatDestroy previous", MatDestroy(&previous));
		CheckPetsc(communicator, "legacy transport MatDestroy left", MatDestroy(&left));
	}
};

struct TransportInput {
	iga::TransportParameters parameters;
	iga::SimulationConfiguration converted;
	iga::CompiledLinearSystem system;
	std::vector<int> labels;
	std::vector<std::array<double, 3>> velocity;
	iga::CaseConfiguration case_configuration;
	iga::ResolvedBoundaryConditions boundaries;
	iga::TransportCouplingPatterns coupling_patterns;
	fs::path case_dir, velocity_path;
};

void RequireRegularOutput(const fs::path& path)
{
	if (fs::exists(path) && !fs::is_regular_file(path))
		throw std::runtime_error("legacy transport output is not a regular file: "+path.string());
}

void RunLegacyTransport(int argc, char** argv, MPI_Comm communicator)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	{
		iga::RequireExecutionResources(communicator, &std::cout);
		ReturnErrors errors;
		errors.Enable(communicator);
		std::unique_ptr<iga::Database> database_owner;
		std::string database_fingerprint;
		iga::CollectiveLocalStage(communicator, "legacy transport database input", [&] {
			if (argc < 3) throw std::runtime_error("usage: iga_transport DATABASE.ntiga CASE_DIR [STEPS] [OUTPUT] [VELOCITY]");
			database_fingerprint = iga::ReadAssetFingerprint(argv[1]);
			database_owner = std::make_unique<iga::Database>(argv[1]);
			iga::ValidatePackedExecution(database_owner->header().ranks, database_owner->header().nodes, 2, ranks);
		});
		iga::RequireCollectiveSameText(communicator, "legacy transport asset database", database_fingerprint);
		iga::RequireCollectivePetscOptions(communicator);
		auto& database = *database_owner;
		std::optional<TransportInput> input;
		iga::AssetFileCatalog assets;
		iga::CollectiveLocalStage(communicator, "legacy transport input catalog", [&] {
			input.emplace();
			input->case_dir = argv[2];
			input->velocity_path = argc >= 6 ? fs::path(argv[5]) : input->case_dir/"initial_velocityfield.txt";
			assets.emplace("legacy transport asset parameters", input->case_dir/"simulation_parameter.txt");
			assets.emplace("legacy transport asset mesh", input->case_dir/"controlmesh.vtk");
			assets.emplace("legacy transport asset velocity", input->velocity_path);
			const auto configuration = input->case_dir/"case_config.json";
			if (fs::exists(configuration)) assets.emplace("legacy transport asset configuration", configuration);
		});
		iga::RequireCollectiveAssetFiles(communicator, assets);
		auto& parameters = input->parameters;
		auto& converted = input->converted;
		auto& system = input->system;
		auto& labels = input->labels;
		auto& velocity = input->velocity;
		auto& case_configuration = input->case_configuration;
		auto& boundaries = input->boundaries;
		auto& coupling_patterns = input->coupling_patterns;
		std::string execution_controls;
		iga::CollectiveLocalStage(communicator, "legacy transport case input", [&] {
			const auto& case_dir = input->case_dir;
			parameters = iga::ReadTransportParameters((case_dir / "simulation_parameter.txt").string());
			if (argc >= 4) {
				const std::string text(argv[3]); std::size_t used = 0;
				const auto steps = std::stoll(text, &used);
				if (used != text.size() || steps < 0 || steps > std::numeric_limits<int>::max())
					throw std::runtime_error("STEPS must be a nonnegative integer within int capacity");
				parameters.steps = static_cast<int>(steps);
			}
			converted = iga::ConvertLegacyNeuronTransport(parameters);
			system = iga::CompileLinearSystem(converted, "neuron_transport");
			labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
			velocity = iga::ReadVelocity(input->velocity_path.string(), database.header().nodes);
			case_configuration = iga::ReadCaseConfiguration((case_dir / "case_config.json").string());
			boundaries = iga::ResolveBoundaryConditions(case_configuration, labels, velocity, parameters);
			coupling_patterns = iga::BuildTransportCouplingPatterns(system, converted);
			execution_controls = std::to_string(parameters.steps)+" "+std::to_string(argc >= 5);
		});
		iga::RequireCollectiveSameText(communicator, "legacy transport execution agreement", execution_controls);
		iga::CollectiveLocalStage(communicator, "legacy transport input logging", [&] {
			if (rank == 0) std::cout << "boundary_config=" << (case_configuration.present ? "case_config.json" : "legacy-defaults")
			<< " transport_nodes=" << boundaries.transport_nodes << '\n';
			iga::FlushCheckedText(std::cout);
		});
		iga::OwnedRowAssembler assembler(database, communicator, 2);
		iga::RequireValidGeometry(assembler.elements(),
			[&assembler](const iga::Element& element) {
				return assembler.OwnsElementByMinimumNode(element);
			}, communicator);
		TransportObjects objects;
		auto& left = objects.left;
		auto& previous = objects.previous;
		auto& current = objects.current;
		auto& next = objects.next;
		auto& rhs = objects.rhs;
		auto& solver = objects.solver;
		left = assembler.CreateMatrix(coupling_patterns.left);
		previous = assembler.CreateMatrix(coupling_patterns.previous);
		CheckPetsc(communicator, "legacy transport MatSetOption", MatSetOption(previous, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
		const auto assembly_start = std::chrono::steady_clock::now();
		iga::CollectiveLocalStage(communicator, "legacy transport element assembly", [&] {
			iga::GenericTransportMatrices element_matrices(coupling_patterns);
			for (const auto& element : assembler.elements()) {
				iga::FullCell4x4x4VolumeQuadratureProvider volume_quadrature(element);
				iga::BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
				iga::BuildGenericTransportElement(element, velocity, system, converted,
					element_matrices, volume_quadrature.Rule(), surface_quadrature.Rule());
				assembler.AddElementMatrix(left, element, element_matrices.left);
				assembler.AddElementMatrix(previous, element, element_matrices.previous);
			}
		});
		iga::OwnedRowAssembler::Assemble(left, communicator);
		iga::OwnedRowAssembler::Assemble(previous, communicator);

		std::vector<PetscInt> boundary_rows;
		std::vector<PetscScalar> boundary_values;
		iga::CollectiveLocalStage(communicator, "legacy transport boundary preparation", [&] {
			for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
				if (boundaries.transport_constrained[static_cast<std::size_t>(node)]) {
					boundary_rows.push_back(static_cast<PetscInt>(2*node));
					boundary_rows.push_back(static_cast<PetscInt>(2*node+1));
					boundary_values.push_back(boundaries.n0[static_cast<std::size_t>(node)]);
					boundary_values.push_back(boundaries.nplus[static_cast<std::size_t>(node)]);
				}
		});
		CheckPetsc(communicator, "legacy transport MatZeroRows left", MatZeroRows(left,
			static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 1.0, nullptr, nullptr));
		CheckPetsc(communicator, "legacy transport MatZeroRows previous", MatZeroRows(previous,
			static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 0.0, nullptr, nullptr));

		current = assembler.CreateVector();
		next = assembler.CreateVector();
		rhs = assembler.CreateVector();
		CheckPetsc(communicator, "legacy transport VecSet", VecSet(current, 0.0));
		CheckPetsc(communicator, "legacy transport VecSet", VecSet(next, 0.0));
		CheckPetsc(communicator, "legacy transport VecSet", VecSet(rhs, 0.0));
		// Only owned rows are inserted. The boundary values are time independent
		// in this legacy contract and are reused without per-step allocation.
		CheckPetsc(communicator, "legacy transport initial values", VecSetValues(current,
			static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES));
		iga::OwnedRowAssembler::Assemble(current, communicator);

		CheckPetsc(communicator, "legacy transport KSPCreate", KSPCreate(communicator, &solver));
		CheckPetsc(communicator, "legacy transport KSPSetOperators", KSPSetOperators(solver, left, left));
		CheckPetsc(communicator, "legacy transport KSPSetType", KSPSetType(solver, KSPGMRES));
		CheckPetsc(communicator, "legacy transport KSPGMRESSetRestart", KSPGMRESSetRestart(solver, 50));
		CheckPetsc(communicator, "legacy transport KSPSetTolerances", KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000));
		PC preconditioner = nullptr;
		CheckPetsc(communicator, "legacy transport KSPGetPC", KSPGetPC(solver, &preconditioner));
		CheckPetsc(communicator, "legacy transport PCSetType", PCSetType(preconditioner, PCBJACOBI));
		CheckPetsc(communicator, "legacy transport KSPSetFromOptions", KSPSetFromOptions(solver));
		iga::RequireKspFactorBackend(solver, left, communicator);
		CheckPetsc(communicator, "legacy transport solver setup", KSPSetUp(solver));
		const auto assembly_end = std::chrono::steady_clock::now();

		std::uint64_t total_iterations = 0;
		const auto solve_start = std::chrono::steady_clock::now();
		for (int step = 0; step < parameters.steps; ++step) {
			CheckPetsc(communicator, "legacy transport MatMult", MatMult(previous, current, rhs));
			CheckPetsc(communicator, "legacy transport boundary values", VecSetValues(rhs,
				static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES));
			iga::OwnedRowAssembler::Assemble(rhs, communicator);
			if (step > 0) CheckPetsc(communicator, "legacy transport VecCopy", VecCopy(current, next));
			CheckPetsc(communicator, "legacy transport KSPSetInitialGuessNonzero", KSPSetInitialGuessNonzero(solver, step > 0 ? PETSC_TRUE : PETSC_FALSE));
			CheckPetsc(communicator, "legacy transport linear solve", KSPSolve(solver, rhs, next));
			KSPConvergedReason reason;
			PetscInt iterations = 0;
			iga::CollectiveLocalStage(communicator, "legacy transport convergence", [&] {
				if (KSPGetConvergedReason(solver, &reason) || KSPGetIterationNumber(solver, &iterations))
					throw std::runtime_error("cannot query transport KSP status");
				if (iterations < 0 || static_cast<std::uint64_t>(iterations)
					> std::numeric_limits<std::uint64_t>::max()-total_iterations)
					throw std::runtime_error("transport KSP iteration count is out of range");
				if (reason <= 0) throw std::runtime_error("transport KSP did not converge at step " + std::to_string(step));
			});
			total_iterations += static_cast<std::uint64_t>(iterations);
			CheckPetsc(communicator, "legacy transport VecSwap", VecSwap(current, next));
		}
		const auto solve_end = std::chrono::steady_clock::now();
		PetscReal norm = 0.0;
		CheckPetsc(communicator, "legacy transport VecNorm", VecNorm(current, NORM_2, &norm));
		if (argc >= 5) {
			iga::PetscGatherObjects output_objects;
			iga::RequireCollectivePetscSuccess(communicator, "legacy output scatter create",
				VecScatterCreateToZero(current, &output_objects.scatter, &output_objects.all));
			iga::RequireCollectivePetscSuccess(communicator, "legacy output scatter begin",
				VecScatterBegin(output_objects.scatter, current, output_objects.all, INSERT_VALUES, SCATTER_FORWARD));
			iga::RequireCollectivePetscSuccess(communicator, "legacy output scatter end",
				VecScatterEnd(output_objects.scatter, current, output_objects.all, INSERT_VALUES, SCATTER_FORWARD));
			iga::CollectiveLocalStage(communicator, "legacy transport output", [&] {
				if (rank != 0) return;
				iga::PetscReadArray view; view.Acquire(output_objects.all);
				PetscInt size = 0;
				if (VecGetSize(output_objects.all, &size)
					|| size != assembler.global_rows()) throw std::runtime_error("legacy transport output size mismatch");
				RequireRegularOutput(argv[4]);
				std::ofstream output(argv[4]);
				if (!output) throw std::runtime_error("cannot create final state output");
				output.precision(17);
				for (std::uint64_t node = 0; node < database.header().nodes; ++node)
					output << node << ' ' << PetscRealPart(view.Data()[2*node]) << ' ' << PetscRealPart(view.Data()[2*node+1]) << '\n';
				output.close();
				if (!output) throw std::runtime_error("cannot write final state output");
				view.Restore();
			});
			output_objects.Close(communicator);
		}
		objects.Close(communicator);
		iga::CollectiveLocalStage(communicator, "legacy transport final summary", [&] {
			if (!std::isfinite(norm)) throw std::runtime_error("legacy transport state norm is not finite");
			if (rank != 0) return;
			const auto assembly_seconds = std::chrono::duration<double>(assembly_end-assembly_start).count();
			const auto solve_seconds = std::chrono::duration<double>(solve_end-solve_start).count();
			std::cout << "transport_v2 nodes=" << database.header().nodes << " elements=" << database.header().elements
				<< " steps=" << parameters.steps << " assembly_s=" << assembly_seconds << " solve_s=" << solve_seconds
				<< " total_iterations=" << total_iterations << " final_l2=" << norm << '\n';
			iga::FlushCheckedText(std::cout);
		});
	}
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA transport solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		RunLegacyTransport(argc, argv, PETSC_COMM_WORLD);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

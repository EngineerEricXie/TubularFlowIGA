#include "BoundaryFlow.hpp"
#include "CaseInput.hpp"
#include "FlowCheckpoint.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "NavierStokesElement.hpp"
#include "OwnedRowAssembler.hpp"
#include "OutletModel.hpp"
#include "OutletCheckpoint.hpp"
#include "PressureTraction.hpp"
#include "TemporalFunction.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct FlowOptions {
	fs::path database;
	fs::path case_dir;
	int max_newton = 12;
	fs::path output;
	fs::path checkpoint;
	fs::path restart;
	int output_every = 0;
	int checkpoint_every = 0;
	int stop_after_step = 0;
};

int ParsePositiveInteger(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	int value = 0;
	try {
		value = std::stoi(text, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(option + " requires a positive integer");
	}
	if (used != text.size() || value <= 0)
		throw std::runtime_error(option + " requires a positive integer");
	return value;
}

FlowOptions ParseOptions(int argc, char** argv)
{
	if (argc < 3) throw std::runtime_error(
		"usage: iga_navier_stokes DATABASE.ntiga CASE_DIR [MAX_NEWTON] [OUTPUT] "
		"[--max-newton N] [--output PATH] [--output-every N] "
		"[--checkpoint PREFIX] [--checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N]");
	FlowOptions options;
	options.database = argv[1];
	options.case_dir = argv[2];
	int positional = 0;
	for (int i = 3; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument.rfind("--", 0) != 0) {
			if (positional == 0) options.max_newton = ParsePositiveInteger(argument, "MAX_NEWTON");
			else if (positional == 1) options.output = argument;
			else throw std::runtime_error("too many positional arguments");
			++positional;
			continue;
		}
		if (i+1 >= argc) throw std::runtime_error(argument + " requires a value");
		const std::string value(argv[++i]);
		if (argument == "--max-newton") options.max_newton = ParsePositiveInteger(value, argument);
		else if (argument == "--output") options.output = value;
		else if (argument == "--output-every") options.output_every = ParsePositiveInteger(value, argument);
		else if (argument == "--checkpoint") options.checkpoint = value;
		else if (argument == "--checkpoint-every") options.checkpoint_every = ParsePositiveInteger(value, argument);
		else if (argument == "--restart") options.restart = value;
		else if (argument == "--stop-after-step")
			options.stop_after_step = ParsePositiveInteger(value, argument);
		else throw std::runtime_error("unknown option: " + argument);
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output or legacy OUTPUT");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

void WriteFlowOutput(Vec state, std::uint64_t nodes, const fs::path& path, int rank)
{
	Vec root = nullptr;
	VecScatter root_scatter = nullptr;
	VecScatterCreateToZero(state, &root_scatter, &root);
	VecScatterBegin(root_scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(root_scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
	int write_failed = 0;
	if (rank == 0) {
		try {
			const PetscScalar* values = nullptr;
			VecGetArrayRead(root, &values);
			std::ofstream output(path);
			std::ofstream pressure(path.string() + ".pressure");
			if (!output || !pressure) throw std::runtime_error("cannot create Navier-Stokes output");
			output.precision(17);
			pressure.precision(17);
			for (std::uint64_t node = 0; node < nodes; ++node) {
				output << PetscRealPart(values[4*node]) << ' ' << PetscRealPart(values[4*node+1]) << ' '
					<< PetscRealPart(values[4*node+2]) << '\n';
				pressure << PetscRealPart(values[4*node+3]) << '\n';
			}
			VecRestoreArrayRead(root, &values);
			if (!output || !pressure) throw std::runtime_error("cannot write Navier-Stokes output");
		} catch (const std::exception&) {
			write_failed = 1;
		}
	}
	MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	VecScatterDestroy(&root_scatter);
	VecDestroy(&root);
	if (write_failed) throw std::runtime_error("cannot write Navier-Stokes output: " + path.string());
}

void WriteCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata, int rank)
{
	const auto state_path = iga::FlowCheckpointStatePath(prefix);
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, state_path.string().c_str(), FILE_MODE_WRITE, &viewer);
	VecView(state, viewer);
	PetscViewerDestroy(&viewer);
	int write_failed = 0;
	if (rank == 0) {
		try {
			iga::WriteFlowCheckpointMetadata(prefix, metadata);
		} catch (const std::exception&) {
			write_failed = 1;
		}
	}
	MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	if (write_failed) throw std::runtime_error(
		"cannot write flow checkpoint metadata: " + iga::FlowCheckpointMetadataPath(prefix).string());
}

void ReadCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata)
{
	if (metadata.state_format != "petsc_binary")
		throw std::runtime_error("CPU flow restart requires petsc_binary checkpoint state");
	fs::path state_path(metadata.state_file);
	if (state_path.is_relative()) state_path = iga::FlowCheckpointMetadataPath(prefix).parent_path()/state_path;
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, state_path.string().c_str(), FILE_MODE_READ, &viewer);
	VecLoad(state, viewer);
	PetscViewerDestroy(&viewer);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA stabilized steady/transient Navier-Stokes solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto options = ParseOptions(argc, argv);
		iga::Database database(options.database.string());
		const auto& case_dir = options.case_dir;
		const int max_newton = options.max_newton;
		const auto labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
		const auto boundary_velocity = iga::ReadVelocity((case_dir / "initial_velocityfield.txt").string(), database.header().nodes);
		iga::ResolvedBoundaryConditions boundaries;
		iga::SimulationConfiguration simulation_configuration;
		std::vector<iga::OutletModelState> outlet_models;
		std::map<int, double> pressure_tractions;
		iga::NavierStokesParameters flow_parameters{1.0, 0.1, 0.0};
		bool configured = false;
		bool transient = false;
		std::string boundary_config;
		if (fs::exists(case_dir / "simulation_config.json")) {
			simulation_configuration = iga::ReadSimulationConfiguration(
				(case_dir / "simulation_config.json").string());
			const auto& flow = iga::FirstNavierStokesSystem(simulation_configuration);
			configured = true;
			transient = flow.time_integration == "backward_euler";
			flow_parameters = {flow.density, flow.viscosity,
				transient ? simulation_configuration.time.dt : 0.0};
			const auto waveform_configuration = transient
				? iga::MaterializeBoundaryWaveforms(simulation_configuration, case_dir.string(), 0.0)
				: simulation_configuration;
			outlet_models = iga::InitializeOutletModels(simulation_configuration, flow);
			const auto initial_configuration = iga::MaterializeOutletPressures(
				waveform_configuration, outlet_models);
			pressure_tractions = iga::ExtractPressureTractions(
				initial_configuration, iga::FirstNavierStokesSystem(initial_configuration));
			boundaries = iga::ResolveFlowBoundaries(initial_configuration,
				iga::FirstNavierStokesSystem(initial_configuration), labels, boundary_velocity);
			boundary_config = "simulation_config.json";
		} else {
			const auto parameters = iga::ReadTransportParameters((case_dir / "simulation_parameter.txt").string());
			const auto case_configuration = iga::ReadCaseConfiguration((case_dir / "case_config.json").string());
			boundaries = iga::ResolveBoundaryConditions(case_configuration, labels, boundary_velocity, parameters);
			boundary_config = case_configuration.present ? "case_config.json" : "legacy-defaults";
		}
		const auto physical_steps = transient ? simulation_configuration.time.steps : 1;
		if (options.stop_after_step > physical_steps)
			throw std::runtime_error("--stop-after-step exceeds configured physical steps");
		const auto run_end_step = options.stop_after_step > 0
			? options.stop_after_step : physical_steps;
		if (!transient)
			for (const auto& model : outlet_models)
				if (model.kind != iga::FieldBoundaryKind::Resistance)
					throw std::runtime_error("RC/RCR outlets require backward_euler flow");
		if (!transient && (!options.restart.empty() || !options.checkpoint.empty()
			|| options.output_every > 0 || options.checkpoint_every > 0))
			throw std::runtime_error("restart, checkpoint, and time-indexed output require transient flow");
		if (rank == 0) std::cout << "boundary_config=" << boundary_config
			<< " viscosity=" << flow_parameters.dynamic_viscosity
			<< " density=" << flow_parameters.density
			<< " time_integration=" << (transient ? "backward_euler" : "steady")
			<< " dt=" << flow_parameters.dt << " steps=" << physical_steps
			<< " run_end_step=" << run_end_step
			<< " velocity_nodes=" << boundaries.velocity_nodes
			<< " pressure_nodes=" << boundaries.pressure_nodes << '\n';
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 4);
		iga::RequireValidGeometry(assembler.elements(), rank, PETSC_COMM_WORLD);
		for (const auto& model : outlet_models) {
			long long local_faces = 0;
			for (const auto& element : assembler.elements()) {
				if (element.owner != rank) continue;
				for (const auto label : element.boundary_labels)
					if (label == model.label) ++local_faces;
			}
			long long global_faces = 0;
			MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD);
			if (global_faces == 0)
				throw std::runtime_error("outlet model label " + std::to_string(model.label)
					+ " has no boundary faces in the .ntiga database; repack with iga_pack");
		}
		for (const auto& traction : pressure_tractions) {
			long long local_faces = 0;
			for (const auto& element : assembler.elements()) {
				if (element.owner != rank) continue;
				for (const auto label : element.boundary_labels)
					if (label == traction.first) ++local_faces;
			}
			long long global_faces = 0;
			MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM,
				PETSC_COMM_WORLD);
			if (global_faces == 0)
				throw std::runtime_error("pressure traction label "
					+std::to_string(traction.first)
					+" has no boundary faces in the .ntiga database; repack with iga_pack");
		}
		Mat jacobian = assembler.CreateMatrix(true);
		Vec state = assembler.CreateVector(), previous = assembler.CreateVector();
		Vec update = assembler.CreateVector(), rhs = assembler.CreateVector();
		VecSet(state, 0.0);
		VecSet(previous, 0.0);
		VecSet(update, 0.0);
		VecSet(rhs, 0.0);

		std::vector<std::int32_t> ghost_nodes;
		for (const auto& element : assembler.elements()) ghost_nodes.insert(ghost_nodes.end(), element.connectivity.begin(), element.connectivity.end());
		std::sort(ghost_nodes.begin(), ghost_nodes.end());
		ghost_nodes.erase(std::unique(ghost_nodes.begin(), ghost_nodes.end()), ghost_nodes.end());
		std::unordered_map<std::int32_t, std::size_t> ghost_position;
		std::vector<PetscInt> ghost_rows;
		ghost_rows.reserve(4*ghost_nodes.size());
		for (std::size_t i = 0; i < ghost_nodes.size(); ++i) {
			ghost_position[ghost_nodes[i]] = i;
			for (int field = 0; field < 4; ++field) ghost_rows.push_back(4*ghost_nodes[i]+field);
		}
		IS source_rows = nullptr;
		ISCreateGeneral(PETSC_COMM_WORLD, static_cast<PetscInt>(ghost_rows.size()), ghost_rows.data(), PETSC_COPY_VALUES, &source_rows);
		Vec ghost_state = nullptr;
		VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(ghost_rows.size()), &ghost_state);
		Vec ghost_previous = nullptr;
		VecDuplicate(ghost_state, &ghost_previous);
		IS destination_rows = nullptr;
		ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(ghost_rows.size()), 0, 1, &destination_rows);
		VecScatter scatter = nullptr;
		VecScatterCreate(state, source_rows, ghost_state, destination_rows, &scatter);
		auto compute_outlet_flows = [&]() {
			std::vector<double> local(outlet_models.size(), 0.0), global(outlet_models.size(), 0.0);
			if (outlet_models.empty()) return global;
			std::unordered_map<int, std::size_t> model_index;
			for (std::size_t i = 0; i < outlet_models.size(); ++i)
				model_index.emplace(outlet_models[i].label, i);
			VecScatterBegin(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			const PetscScalar* values = nullptr;
			VecGetArrayRead(ghost_state, &values);
			for (const auto& element : assembler.elements()) {
				if (element.owner != rank) continue;
				std::vector<std::array<double, 4>> nodal(element.connectivity.size());
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto position = ghost_position.at(element.connectivity[a]);
					for (int field = 0; field < 4; ++field)
						nodal[a][field] = PetscRealPart(values[4*position+field]);
				}
				for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
					const auto found = model_index.find(element.boundary_labels[face]);
					if (found != model_index.end())
						local[found->second] += iga::IntegrateBoundaryFlow(element, face, nodal);
				}
			}
			VecRestoreArrayRead(ghost_state, &values);
			MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()),
				MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
			return global;
		};

		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto index = static_cast<std::size_t>(node);
			if (boundaries.velocity_constrained[index])
				for (int field = 0; field < 3; ++field) boundary_rows.push_back(static_cast<PetscInt>(4*node+field));
			if (boundaries.pressure_constrained[index]) boundary_rows.push_back(static_cast<PetscInt>(4*node+3));
		}
		int start_step = 0;
		if (transient && !options.restart.empty()) {
			const auto metadata = iga::ReadFlowCheckpointMetadata(options.restart);
			iga::ValidateFlowCheckpoint(metadata, database.header().nodes, physical_steps,
				flow_parameters.dt, flow_parameters.density, flow_parameters.dynamic_viscosity);
			ReadCheckpoint(state, options.restart, metadata);
			iga::RestoreOutletCheckpoint(metadata, outlet_models);
			start_step = metadata.completed_step;
			if (rank == 0) std::cout << "restart=" << options.restart.string()
				<< " completed_step=" << start_step
				<< " physical_time=" << metadata.physical_time << '\n';
		} else if (transient) {
			std::vector<PetscScalar> initial_values;
			initial_values.reserve(boundary_rows.size());
			for (auto row : boundary_rows) {
				const auto node = static_cast<std::uint64_t>(row/4);
				const auto field = row%4;
				const auto index = static_cast<std::size_t>(node);
				initial_values.push_back(field < 3
					? boundaries.velocity[index][static_cast<std::size_t>(field)]
					: boundaries.pressure[index]);
			}
			VecSetValues(state, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), initial_values.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(state);
		}
		VecCopy(state, previous);

		KSP solver = nullptr;
		KSPCreate(PETSC_COMM_WORLD, &solver);
		KSPSetType(solver, KSPFGMRES);
		KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
		PC pc = nullptr;
		KSPGetPC(solver, &pc);
		PCSetType(pc, PCBJACOBI);
		KSPSetFromOptions(solver);

		PetscInt total_linear_iterations = 0;
		const auto start = std::chrono::steady_clock::now();
		for (int step = start_step; step < run_end_step; ++step) {
			if (step > 0) VecCopy(state, previous);
			const auto physical_time = transient ? (step+1)*flow_parameters.dt : 0.0;
			iga::SimulationConfiguration step_configuration;
			if (configured)
				step_configuration = transient
					? iga::MaterializeBoundaryWaveforms(
						simulation_configuration, case_dir.string(), physical_time)
					: simulation_configuration;
			std::vector<double> previous_capacitor_pressure(outlet_models.size());
			for (std::size_t i = 0; i < outlet_models.size(); ++i)
				previous_capacitor_pressure[i] = outlet_models[i].capacitor_pressure;
			bool outlet_converged = false;
			const int maximum_outlet_iterations = outlet_models.empty() ? 1 : 12;
			for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
				if (configured) {
					const auto boundary_configuration = iga::MaterializeOutletPressures(
						step_configuration, outlet_models);
					pressure_tractions = iga::ExtractPressureTractions(
						boundary_configuration,
						iga::FirstNavierStokesSystem(boundary_configuration));
					boundaries = iga::ResolveFlowBoundaries(boundary_configuration,
						iga::FirstNavierStokesSystem(boundary_configuration), labels, boundary_velocity);
				}
			PetscReal initial_residual = -1.0;
			bool converged = false;
			for (int nonlinear = 0; nonlinear < max_newton; ++nonlinear) {
			const auto iteration_start = std::chrono::steady_clock::now();
			MatZeroEntries(jacobian);
			VecSet(rhs, 0.0);
			VecScatterBegin(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			if (transient) {
				VecScatterBegin(scatter, previous, ghost_previous, INSERT_VALUES, SCATTER_FORWARD);
				VecScatterEnd(scatter, previous, ghost_previous, INSERT_VALUES, SCATTER_FORWARD);
			}
			const PetscScalar* ghost_values = nullptr;
			VecGetArrayRead(ghost_state, &ghost_values);
			const PetscScalar* previous_values = nullptr;
			if (transient) VecGetArrayRead(ghost_previous, &previous_values);
			for (const auto& element : assembler.elements()) {
				std::vector<std::array<double,4>> nodal(element.connectivity.size());
				std::vector<std::array<double,4>> previous_nodal;
				if (transient) previous_nodal.resize(element.connectivity.size());
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto position = ghost_position.at(element.connectivity[a]);
					for (int field = 0; field < 4; ++field) {
						nodal[a][field] = PetscRealPart(ghost_values[4*position+field]);
						if (transient) previous_nodal[a][field] = PetscRealPart(previous_values[4*position+field]);
					}
				}
				auto local = iga::BuildNavierStokesElement(
					element, nodal, previous_nodal, flow_parameters);
				for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
					const auto traction = pressure_tractions.find(
						element.boundary_labels[face]);
					if (traction == pressure_tractions.end()) continue;
					const auto surface = iga::IntegrateBoundaryPressureTraction(
						element, face, traction->second);
					for (std::size_t row = 0; row < surface.size(); ++row)
						local.negative_residual[row] += surface[row];
				}
				assembler.AddElementMatrix(jacobian, element, local.jacobian);
				assembler.AddElementVector(rhs, element, local.negative_residual);
			}
			if (transient) VecRestoreArrayRead(ghost_previous, &previous_values);
			VecRestoreArrayRead(ghost_state, &ghost_values);
			iga::OwnedRowAssembler::Assemble(jacobian);
			iga::OwnedRowAssembler::Assemble(rhs);

			const PetscScalar* owned_state = nullptr;
			VecGetArrayRead(state, &owned_state);
			std::vector<PetscScalar> boundary_update;
			boundary_update.reserve(boundary_rows.size());
			for (auto row : boundary_rows) {
				const auto node = static_cast<std::uint64_t>(row / 4);
				const auto field = row % 4;
				const auto index = static_cast<std::size_t>(node);
				const double target = field < 3
					? boundaries.velocity[index][static_cast<std::size_t>(field)]
					: boundaries.pressure[index];
				const auto local_row = static_cast<std::size_t>(row - 4*assembler.node_begin());
				boundary_update.push_back(target - PetscRealPart(owned_state[local_row]));
			}
			VecRestoreArrayRead(state, &owned_state);
			MatZeroRows(jacobian, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 1.0, nullptr, nullptr);
			VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_update.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(rhs);
			PetscReal residual_norm = 0.0;
			VecNorm(rhs, NORM_2, &residual_norm);
			if (initial_residual < 0.0) initial_residual = residual_norm;
			const auto nonlinear_tolerance = std::max<PetscReal>(1e-10, 1e-5*initial_residual);
			if (residual_norm <= nonlinear_tolerance) {
				if (rank == 0) std::cout << "step=" << step+1 << " time=" << physical_time
					<< " converged newton=" << nonlinear << " residual_l2=" << residual_norm
					<< " tolerance=" << nonlinear_tolerance << " assembly_s="
					<< std::chrono::duration<double>(std::chrono::steady_clock::now()-iteration_start).count() << '\n';
				converged = true;
				break;
			}
			const auto linear_start = std::chrono::steady_clock::now();
			KSPSetOperators(solver, jacobian, jacobian);
			KSPSolve(solver, rhs, update);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver, &reason);
			if (reason <= 0) throw std::runtime_error("Navier-Stokes linear solve failed at nonlinear iteration " + std::to_string(nonlinear));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver, &iterations);
			total_linear_iterations += iterations;
			PetscReal linear_residual = 0.0;
			KSPGetResidualNorm(solver, &linear_residual);
			PetscReal update_norm = 0.0;
			VecNorm(update, NORM_2, &update_norm);
			VecAXPY(state, 1.0, update);
			if (rank == 0) std::cout << "step=" << step+1 << " time=" << physical_time
				<< " newton=" << nonlinear << " residual_l2=" << residual_norm << " update_l2=" << update_norm
				<< " linear_iterations=" << iterations << " linear_residual=" << linear_residual
				<< " assembly_s=" << std::chrono::duration<double>(linear_start-iteration_start).count()
				<< " linear_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-linear_start).count() << '\n';
			if (update_norm < 1e-10) { converged = true; break; }
			}
			if (!converged) throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON at physical step "
				+ std::to_string(step+1));
				if (outlet_models.empty()) {
					outlet_converged = true;
					break;
				}
				const auto flows = compute_outlet_flows();
				const auto evaluated = iga::EvaluateOutletCoupling(outlet_models,
					previous_capacitor_pressure, flows, flow_parameters.dt);
				const auto outlet_tolerance = iga::OutletCouplingTolerance(evaluated);
				if (rank == 0) std::cout << "step=" << step+1 << " time=" << physical_time
					<< " outlet_iteration=" << coupling
					<< " pressure_change=" << evaluated.maximum_pressure_change
					<< " tolerance=" << outlet_tolerance << '\n';
				if (evaluated.maximum_pressure_change <= outlet_tolerance) {
					iga::CommitOutletCoupling(outlet_models, evaluated);
					for (std::size_t i = 0; i < outlet_models.size(); ++i) {
						if (rank == 0) std::cout << "outlet label=" << outlet_models[i].label
							<< " flow=" << outlet_models[i].flow
							<< " pressure=" << outlet_models[i].pressure
							<< " capacitor_pressure=" << outlet_models[i].capacitor_pressure << '\n';
					}
					outlet_converged = true;
					break;
				}
				iga::RelaxOutletCoupling(outlet_models, evaluated);
			}
			if (!outlet_converged)
				throw std::runtime_error("outlet fixed-point iteration did not converge at physical step "
					+ std::to_string(step+1));
			const auto completed_step = step+1;
			if (options.output_every > 0 && completed_step%options.output_every == 0)
				WriteFlowOutput(state, database.header().nodes,
					iga::TimeIndexedPath(options.output, completed_step), rank);
			if (!options.checkpoint.empty()
				&& (completed_step == run_end_step
					|| (options.checkpoint_every > 0 && completed_step%options.checkpoint_every == 0))) {
				iga::FlowCheckpointMetadata metadata;
				metadata.nodes = database.header().nodes;
				metadata.completed_step = completed_step;
				metadata.physical_time = completed_step*flow_parameters.dt;
				metadata.dt = flow_parameters.dt;
				metadata.density = flow_parameters.density;
				metadata.viscosity = flow_parameters.dynamic_viscosity;
				metadata.state_file = iga::FlowCheckpointStatePath(options.checkpoint).filename().string();
				metadata.state_format = "petsc_binary";
				iga::AppendOutletCheckpoint(outlet_models, metadata);
				WriteCheckpoint(state, options.checkpoint, metadata, rank);
				if (rank == 0) std::cout << "checkpoint=" << options.checkpoint.string()
					<< " completed_step=" << completed_step << '\n';
			}
		}
		const auto end = std::chrono::steady_clock::now();
		PetscReal state_norm = 0.0;
		VecNorm(state, NORM_2, &state_norm);
		const PetscScalar* owned_values = nullptr;
		VecGetArrayRead(state, &owned_values);
		double local_velocity_squared = 0.0, local_pressure_squared = 0.0;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto local = static_cast<std::size_t>(node - assembler.node_begin()) * 4;
			for (int field = 0; field < 3; ++field) {
				const auto value = PetscRealPart(owned_values[local+field]);
				local_velocity_squared += value*value;
			}
			const auto pressure = PetscRealPart(owned_values[local+3]);
			local_pressure_squared += pressure*pressure;
		}
		VecRestoreArrayRead(state, &owned_values);
		double velocity_squared = 0.0, pressure_squared = 0.0;
		MPI_Allreduce(&local_velocity_squared, &velocity_squared, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce(&local_pressure_squared, &pressure_squared, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		if (rank == 0) std::cout << "navier_stokes_v2 seconds=" << std::chrono::duration<double>(end-start).count()
			<< " total_linear_iterations=" << total_linear_iterations << " state_l2=" << state_norm
			<< " velocity_l2=" << std::sqrt(velocity_squared) << " pressure_l2=" << std::sqrt(pressure_squared) << '\n';

		if (!options.output.empty())
			WriteFlowOutput(state, database.header().nodes, options.output, rank);
		KSPDestroy(&solver); VecScatterDestroy(&scatter); ISDestroy(&destination_rows);
		VecDestroy(&ghost_previous); VecDestroy(&ghost_state); ISDestroy(&source_rows);
		VecDestroy(&rhs); VecDestroy(&update); VecDestroy(&previous); VecDestroy(&state); MatDestroy(&jacobian);
	} catch (const std::exception& e) {
		std::cerr << "rank " << rank << ": " << e.what() << '\n'; status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

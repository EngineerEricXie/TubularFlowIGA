#include "CaseInput.hpp"
#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"
#include "SimulationConfig.hpp"
#include "TemporalFunction.hpp"
#include "FlowCheckpoint.hpp"
#include "TransportCheckpoint.hpp"
#include "VelocitySeries.hpp"

#include <petscksp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct TransportOptions {
	fs::path database;
	fs::path case_dir;
	std::string system;
	fs::path output;
	fs::path velocity;
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
		throw std::runtime_error(option+" requires a positive integer");
	}
	if (used != text.size() || value <= 0)
		throw std::runtime_error(option+" requires a positive integer");
	return value;
}

TransportOptions ParseOptions(int argc, char** argv)
{
	if (argc < 3) throw std::runtime_error(
		"usage: iga_solve DATABASE.ntiga CASE_DIR [SYSTEM] [OUTPUT] [VELOCITY] "
		"[--system NAME] [--output PATH] [--velocity PATH] [--output-every N] "
		"[--checkpoint PREFIX] [--checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N]");
	TransportOptions options;
	options.database = argv[1];
	options.case_dir = argv[2];
	int positional = 0;
	for (int i = 3; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument.rfind("--", 0) != 0) {
			if (positional == 0) options.system = argument;
			else if (positional == 1) options.output = argument;
			else if (positional == 2) options.velocity = argument;
			else throw std::runtime_error("too many positional arguments");
			++positional;
			continue;
		}
		if (i+1 >= argc) throw std::runtime_error(argument+" requires a value");
		const std::string value(argv[++i]);
		if (argument == "--system") options.system = value;
		else if (argument == "--output") options.output = value;
		else if (argument == "--velocity") options.velocity = value;
		else if (argument == "--output-every")
			options.output_every = ParsePositiveInteger(value, argument);
		else if (argument == "--checkpoint") options.checkpoint = value;
		else if (argument == "--checkpoint-every")
			options.checkpoint_every = ParsePositiveInteger(value, argument);
		else if (argument == "--restart") options.restart = value;
		else if (argument == "--stop-after-step")
			options.stop_after_step = ParsePositiveInteger(value, argument);
		else throw std::runtime_error("unknown option: "+argument);
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

void WriteTransportOutput(Vec state, std::uint64_t nodes,
	const std::vector<std::string>& fields, const fs::path& path, int rank)
{
	Vec root = nullptr;
	VecScatter scatter = nullptr;
	VecScatterCreateToZero(state, &scatter, &root);
	VecScatterBegin(scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
	int write_failed = 0;
	if (rank == 0) {
		try {
			const PetscScalar* values = nullptr;
			VecGetArrayRead(root, &values);
			std::ofstream output(path);
			std::ofstream metadata(path.string()+".fields");
			if (!output || !metadata)
				throw std::runtime_error("cannot create configured PDE output");
			output.precision(17);
			for (const auto& name : fields) metadata << name << '\n';
			for (std::uint64_t node = 0; node < nodes; ++node) {
				output << node;
				for (std::size_t field = 0; field < fields.size(); ++field)
					output << ' ' << PetscRealPart(values[node*fields.size()+field]);
				output << '\n';
			}
			VecRestoreArrayRead(root, &values);
			if (!output || !metadata)
				throw std::runtime_error("cannot write configured PDE output");
		} catch (const std::exception&) {
			write_failed = 1;
		}
	}
	MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	VecScatterDestroy(&scatter);
	VecDestroy(&root);
	if (write_failed)
		throw std::runtime_error("cannot write configured PDE output: "+path.string());
}

void WriteTransportCheckpoint(Vec state, const fs::path& prefix,
	const iga::TransportCheckpointMetadata& metadata, int rank)
{
	const auto state_path = iga::TransportCheckpointStatePath(prefix);
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(
		PETSC_COMM_WORLD, state_path.string().c_str(), FILE_MODE_WRITE, &viewer);
	VecView(state, viewer);
	PetscViewerDestroy(&viewer);
	int write_failed = 0;
	if (rank == 0) {
		try {
			iga::WriteTransportCheckpointMetadata(prefix, metadata);
		} catch (const std::exception&) {
			write_failed = 1;
		}
	}
	MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	if (write_failed)
		throw std::runtime_error("cannot write transport checkpoint metadata: "
			+iga::TransportCheckpointMetadataPath(prefix).string());
}

void ReadTransportCheckpoint(Vec state, const fs::path& prefix,
	const iga::TransportCheckpointMetadata& metadata)
{
	if (metadata.state_format != "petsc_binary")
		throw std::runtime_error("CPU transport restart requires petsc_binary state");
	fs::path state_path(metadata.state_file);
	if (state_path.is_relative())
		state_path = iga::TransportCheckpointMetadataPath(prefix).parent_path()/state_path;
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(
		PETSC_COMM_WORLD, state_path.string().c_str(), FILE_MODE_READ, &viewer);
	VecLoad(state, viewer);
	PetscViewerDestroy(&viewer);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA configured PDE solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto options = ParseOptions(argc, argv);
		iga::Database database(options.database.string());
		const auto& case_dir = options.case_dir;
		const auto configuration = iga::ReadSimulationConfiguration(
			(case_dir / "simulation_config.json").string());
		const auto system_name = !options.system.empty() ? options.system
			: iga::FirstLinearTransportSystem(configuration);
		const auto system = iga::CompileLinearSystem(configuration, system_name);
		if (options.stop_after_step > system.steps)
			throw std::runtime_error("--stop-after-step exceeds configured transport steps");
		const auto run_end_step = options.stop_after_step > 0
			? options.stop_after_step : system.steps;
		const auto labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
		std::vector<std::array<double, 3>> prescribed_velocity;
		std::vector<iga::VelocitySnapshot> velocity_snapshots;
		const iga::VelocitySourceDefinition* velocity_source = nullptr;
		if (system.velocity_source == "prescribed") {
			const auto velocity_path = !options.velocity.empty()
				? options.velocity : case_dir / "initial_velocityfield.txt";
			prescribed_velocity = iga::ReadVelocity(velocity_path.string(), database.header().nodes);
		} else {
			if (!options.velocity.empty())
				throw std::runtime_error("a positional VELOCITY file cannot override a configured snapshot_series");
			velocity_source = &iga::FindVelocitySource(configuration, system.velocity_source);
			velocity_snapshots = iga::ReadVelocityManifest(case_dir / velocity_source->manifest);
		}
		const auto initial_configuration = iga::MaterializeBoundaryWaveforms(configuration, case_dir, 0.0);
		const auto boundaries = iga::ResolveScalarBoundaries(initial_configuration, system, labels);
		const auto initial = iga::InitialScalarValues(configuration, system);
		const auto fields = system.fields.size();
		if (rank == 0) {
			std::cout << "configuration=simulation_config.json system=" << system.name
				<< " fields=" << fields << " velocity_source=" << system.velocity_source
				<< " dirichlet_dofs=" << boundaries.constrained_dofs << '\n';
			for (std::size_t i = 0; i < system.fields.size(); ++i)
				std::cout << "field[" << i << "]=" << system.fields[i] << '\n';
		}

		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, fields);
		iga::RequireValidGeometry(assembler.elements(), rank, PETSC_COMM_WORLD);
		std::map<int, long long> surface_faces;
		for (const auto& boundary : configuration.boundaries)
			for (const auto& condition : boundary.conditions)
				if (system.field_index.count(condition.field)
					&& (condition.kind == iga::FieldBoundaryKind::Flux
						|| condition.kind == iga::FieldBoundaryKind::Robin))
					surface_faces.emplace(boundary.label, 0);
		for (const auto& element : assembler.elements()) {
			if (element.owner != rank) continue;
			for (const auto label : element.boundary_labels) {
				auto found = surface_faces.find(label);
				if (found != surface_faces.end()) ++found->second;
			}
		}
		for (auto& entry : surface_faces) {
			long long global_faces = 0;
			MPI_Allreduce(&entry.second, &global_faces, 1, MPI_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD);
			entry.second = global_faces;
			if (global_faces == 0)
				throw std::runtime_error("configured scalar surface boundary label "
					+ std::to_string(entry.first) + " has no boundary faces in the .ntiga database; repack with iga_pack");
		}
		Mat left = assembler.CreateMatrix();
		Mat previous = assembler.CreateMatrix();
		MatSetOption(left, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
		MatSetOption(previous, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
		MatSetOption(previous, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
		Vec forcing = assembler.CreateVector();
		VecSet(forcing, 0.0);

		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
			for (std::size_t field = 0; field < fields; ++field)
				if (boundaries.constrained[static_cast<std::size_t>(node)*fields+field])
					boundary_rows.push_back(static_cast<PetscInt>(node*fields+field));
		double assembly_seconds = 0.0;
		auto assemble_operators = [&](const std::vector<std::array<double, 3>>& velocity) {
			const auto start = std::chrono::steady_clock::now();
			MatZeroEntries(left);
			MatZeroEntries(previous);
			VecSet(forcing, 0.0);
			for (const auto& element : assembler.elements()) {
				const auto matrices = iga::BuildGenericTransportElement(
					element, velocity, system, configuration);
				assembler.AddElementMatrix(left, element, matrices.left);
				assembler.AddElementMatrix(previous, element, matrices.previous);
				assembler.AddElementVector(forcing, element, matrices.source);
			}
			iga::OwnedRowAssembler::Assemble(left);
			iga::OwnedRowAssembler::Assemble(previous);
			iga::OwnedRowAssembler::Assemble(forcing);
			MatZeroRows(left, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 1.0, nullptr, nullptr);
			MatZeroRows(previous, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 0.0, nullptr, nullptr);
			assembly_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-start).count();
		};
		auto velocity_at = [&](double time) {
			const auto interpolation = iga::ResolveVelocityInterpolation(
				velocity_snapshots, time, velocity_source->out_of_range);
			auto lower = iga::ReadVelocity(
				(case_dir/velocity_snapshots[interpolation.lower].file).string(), database.header().nodes);
			if (interpolation.lower == interpolation.upper) return lower;
			const auto upper = iga::ReadVelocity(
				(case_dir/velocity_snapshots[interpolation.upper].file).string(), database.header().nodes);
			for (std::size_t node = 0; node < lower.size(); ++node)
				for (int component = 0; component < 3; ++component)
					lower[node][component] = (1.0-interpolation.upper_weight)*lower[node][component]
						+ interpolation.upper_weight*upper[node][component];
			return lower;
		};
		if (!velocity_source) assemble_operators(prescribed_velocity);

		Vec current = assembler.CreateVector();
		Vec next = assembler.CreateVector();
		Vec rhs = assembler.CreateVector();
		VecSet(current, 0.0);
		VecSet(next, 0.0);
		VecSet(rhs, 0.0);
		PetscScalar* current_array = nullptr;
		VecGetArray(current, &current_array);
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
			for (std::size_t field = 0; field < fields; ++field) {
				const auto local = static_cast<std::size_t>(node-assembler.node_begin())*fields+field;
				const auto global = static_cast<std::size_t>(node)*fields+field;
				current_array[local] = boundaries.constrained[global] ? boundaries.value[global] : initial[field];
			}
		VecRestoreArray(current, &current_array);
		int start_step = 0;
		if (!options.restart.empty()) {
			const auto metadata = iga::ReadTransportCheckpointMetadata(options.restart);
			iga::ValidateTransportCheckpoint(metadata, database.header().nodes,
				system.fields, system.name, system.velocity_source, system.steps, system.dt);
			ReadTransportCheckpoint(current, options.restart, metadata);
			start_step = metadata.completed_step;
			if (rank == 0) std::cout << "restart=" << options.restart.string()
				<< " completed_step=" << start_step
				<< " physical_time=" << metadata.physical_time << '\n';
		}

		KSP solver = nullptr;
		KSPCreate(PETSC_COMM_WORLD, &solver);
		if (!velocity_source) KSPSetOperators(solver, left, left);
		KSPSetType(solver, KSPGMRES);
		KSPGMRESSetRestart(solver, 50);
		KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
		PC preconditioner = nullptr;
		KSPGetPC(solver, &preconditioner);
		PCSetType(preconditioner, PCBJACOBI);
		KSPSetFromOptions(solver);
		if (!velocity_source) KSPSetUp(solver);

		PetscInt total_iterations = 0;
		const auto solve_start = std::chrono::steady_clock::now();
		double linear_seconds = 0.0;
		for (int step = start_step; step < run_end_step; ++step) {
			if (velocity_source) {
				const auto velocity = velocity_at((step+1)*system.dt);
				assemble_operators(velocity);
				KSPSetOperators(solver, left, left);
				KSPSetUp(solver);
			}
			const auto step_configuration = iga::MaterializeBoundaryWaveforms(
				configuration, case_dir, (step+1)*system.dt);
			const auto step_boundaries = iga::ResolveScalarBoundaries(step_configuration, system, labels);
			MatMult(previous, current, rhs);
			VecAXPY(rhs, 1.0, forcing);
			std::vector<PetscScalar> boundary_values;
			boundary_values.reserve(boundary_rows.size());
			for (auto row : boundary_rows)
				boundary_values.push_back(step_boundaries.value[static_cast<std::size_t>(row)]);
			VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(rhs);
			if (step > 0) VecCopy(current, next);
			KSPSetInitialGuessNonzero(solver, step > 0 ? PETSC_TRUE : PETSC_FALSE);
			const auto linear_start = std::chrono::steady_clock::now();
			KSPSolve(solver, rhs, next);
			linear_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-linear_start).count();
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver, &reason);
			if (reason <= 0) throw std::runtime_error("KSP did not converge at step " + std::to_string(step));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver, &iterations);
			total_iterations += iterations;
			VecSwap(current, next);
			const auto completed_step = step+1;
			if (options.output_every > 0
				&& completed_step%options.output_every == 0)
				WriteTransportOutput(current, database.header().nodes, system.fields,
					iga::TimeIndexedPath(options.output, completed_step), rank);
			if (!options.checkpoint.empty()
				&& (completed_step == run_end_step
					|| (options.checkpoint_every > 0
						&& completed_step%options.checkpoint_every == 0))) {
				iga::TransportCheckpointMetadata metadata;
				metadata.nodes = database.header().nodes;
				metadata.fields = system.fields;
				metadata.system = system.name;
				metadata.velocity_source = system.velocity_source;
				metadata.completed_step = completed_step;
				metadata.physical_time = completed_step*system.dt;
				metadata.dt = system.dt;
				metadata.state_file = iga::TransportCheckpointStatePath(
					options.checkpoint).filename().string();
				metadata.state_format = "petsc_binary";
				WriteTransportCheckpoint(current, options.checkpoint, metadata, rank);
				if (rank == 0) std::cout << "checkpoint=" << options.checkpoint.string()
					<< " completed_step=" << completed_step << '\n';
			}
		}
		const auto solve_end = std::chrono::steady_clock::now();

		if (!options.output.empty())
			WriteTransportOutput(current, database.header().nodes,
				system.fields, options.output, rank);
		PetscReal norm = 0.0;
		VecNorm(current, NORM_2, &norm);
		if (rank == 0) std::cout << "iga_solve system=" << system.name
			<< " steps=" << system.steps << " run_end_step=" << run_end_step
			<< " assembly_s=" << assembly_seconds
			<< " solve_s=" << linear_seconds
			<< " time_loop_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
			<< " total_iterations=" << total_iterations << " final_l2=" << norm << '\n';
		KSPDestroy(&solver);
		VecDestroy(&rhs); VecDestroy(&next); VecDestroy(&current); VecDestroy(&forcing);
		MatDestroy(&previous); MatDestroy(&left);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

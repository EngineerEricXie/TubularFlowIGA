#include "CheckedText.hpp"
#include "ExecutionResources.hpp"
#include "CaseInput.hpp"
#include "CollectiveAssetInput.hpp"
#include "CollectivePetscOptions.hpp"
#include "PetscCheckpointRead.hpp"
#include "PetscCheckpointWrite.hpp"
#include "PetscGather.hpp"
#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "IgaDatabase.hpp"
#include "MemoryReport.hpp"
#include "PhaseProfile.hpp"
#include "PetscPhaseProfile.hpp"
#include "OwnedRowAssembler.hpp"
#include "SimulationConfig.hpp"
#include "TemporalFunction.hpp"
#include "TransportBoundaryPreflight.hpp"
#include "FlowCheckpoint.hpp"
#include "TransportCheckpoint.hpp"
#include "VelocitySeries.hpp"
#include "PhysiologyOutput.hpp"
#include "TemporalVtkHdf.hpp"
#include "VtkOutput.hpp"

#include <petscksp.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

namespace {

void CheckPetsc(const char* stage, PetscErrorCode code)
{
	iga::RequireCollectivePetscSuccess(PETSC_COMM_WORLD, stage, code);
}

class ReturnErrors {
public:
	void Enable()
	{
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport error handler", [&] {
			if (PetscPushErrorHandler(PetscReturnErrorHandler, nullptr))
				throw std::runtime_error("cannot install returning PETSc error handler");
			active_ = true;
		});
	}
	~ReturnErrors() { if (active_) PetscPopErrorHandler(); }
private:
	bool active_ = false;
};

struct TransportPetscObjects {
	Mat left = nullptr, previous = nullptr;
	Vec forcing = nullptr, current = nullptr, next = nullptr, rhs = nullptr;
	KSP solver = nullptr;
	~TransportPetscObjects()
	{
		if (solver) KSPDestroy(&solver);
		if (rhs) VecDestroy(&rhs);
		if (next) VecDestroy(&next);
		if (current) VecDestroy(&current);
		if (forcing) VecDestroy(&forcing);
		if (previous) MatDestroy(&previous);
		if (left) MatDestroy(&left);
	}
};

class WriteArray {
public:
	void Acquire(Vec vector)
	{
		if (VecGetArray(vector, &values_)) throw std::runtime_error("VecGetArray failed");
		vector_ = vector;
	}
	~WriteArray() { if (vector_) VecRestoreArray(vector_, &values_); }
	PetscScalar* Data() const { return values_; }
	void Restore()
	{
		if (VecRestoreArray(vector_, &values_)) throw std::runtime_error("VecRestoreArray failed");
		vector_ = nullptr;
	}
private:
	Vec vector_ = nullptr;
	PetscScalar* values_ = nullptr;
};

void RequireRegularOutput(const fs::path& path)
{
	if (fs::exists(path) && !fs::is_regular_file(path))
		throw std::runtime_error("output must be a regular file: "+path.string());
}

void SetupTransportKsp(KSP solver, Mat matrix)
{
	iga::RequireKspFactorBackend(solver, matrix, PETSC_COMM_WORLD);
	iga::PhaseScope phase(iga::ProfilePhase::SolverSetup);
	CheckPetsc("transport solver setup", KSPSetUp(solver));
	CheckPetsc("transport block solver setup", KSPSetUpOnBlocks(solver));
}

void EnableMemoryTracking(int argc, char** argv)
{
	bool requested = false;
	for (int i = 1; i < argc; ++i)
		if (std::string(argv[i]) == "--memory-report") requested = true;
	if (!requested) return;
	std::string petsc_options = std::getenv("PETSC_OPTIONS")
		? std::getenv("PETSC_OPTIONS") : "";
	if (petsc_options.find("-malloc_debug") == std::string::npos) {
		if (!petsc_options.empty()) petsc_options += ' ';
		petsc_options += "-malloc_debug";
		if (setenv("PETSC_OPTIONS", petsc_options.c_str(), 1) != 0)
			throw std::runtime_error("cannot enable PETSc allocation tracking");
	}
}

struct TransportOptions {
	fs::path database;
	fs::path case_dir;
	std::string system;
	fs::path output;
	fs::path velocity;
	fs::path checkpoint;
	fs::path restart;
	fs::path memory_report;
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
	int output_every = 0;
	int checkpoint_every = 0;
	int stop_after_step = 0;
};

struct TransportInput {
	iga::SimulationConfiguration configuration;
	iga::CompiledLinearSystem system;
	std::vector<int> labels;
	std::vector<std::array<double, 3>> prescribed_velocity;
	std::vector<iga::VelocitySnapshot> velocity_snapshots;
	const iga::VelocitySourceDefinition* velocity_source = nullptr;
	iga::ResolvedScalarBoundaries boundaries;
	std::vector<double> initial;
	iga::TransportCouplingPatterns coupling_patterns;
	int run_end_step = 0;
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
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
		"[--stop-after-step N] [--memory-report PATH] "
		"[--visualization-format auto|vtu|vtkhdf]");
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
		else if (argument == "--memory-report") options.memory_report = value;
		else if (argument == "--visualization-format")
			options.visualization_format = iga::ParseVisualizationFormat(value);
		else if (argument == "--stop-after-step")
			options.stop_after_step = ParsePositiveInteger(value, argument);
		else throw std::runtime_error("unknown option: "+argument);
		if (PetscOptionsClearValue(nullptr, argument.c_str()))
			throw std::runtime_error("cannot consume application option: "+argument);
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

void WriteTransportOutput(Vec state, std::uint64_t nodes,
	const std::vector<std::string>& fields, const fs::path& path,
	const fs::path& mesh_path, const fs::path& vtk_path, double physical_time,
	const iga::PhysiologyDefinition& physiology, int rank,
	iga::VisualizationFormat visualization_format,
	iga::TemporalVtkHdfWriter* vtkhdf)
{
	iga::PhaseScope output_phase(iga::ProfilePhase::Output);
	iga::PetscGatherObjects objects;
	auto& root = objects.all;
	auto& scatter = objects.scatter;
	iga::PhaseScope gather_phase(iga::ProfilePhase::Communication);
	CheckPetsc("transport output scatter create", VecScatterCreateToZero(state, &scatter, &root));
	CheckPetsc("transport output scatter begin", VecScatterBegin(scatter, state, root, INSERT_VALUES, SCATTER_FORWARD));
	CheckPetsc("transport output scatter end", VecScatterEnd(scatter, state, root, INSERT_VALUES, SCATTER_FORWARD));
	gather_phase.Stop();
	iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport field output", [&] {
		if (rank == 0) {
			iga::PetscReadArray view;
			view.Acquire(root);
			const auto* values = view.Data();
			std::vector<double> interleaved(static_cast<std::size_t>(nodes)*fields.size());
			RequireRegularOutput(path);
			RequireRegularOutput(path.string()+".fields");
			std::ofstream output(path);
			std::ofstream metadata(path.string()+".fields");
			if (!output || !metadata)
				throw std::runtime_error("cannot create configured PDE output");
			output.precision(17);
			for (const auto& name : fields) metadata << name << '\n';
			for (std::uint64_t node = 0; node < nodes; ++node) {
				output << node;
				for (std::size_t field = 0; field < fields.size(); ++field) {
					const auto value = PetscRealPart(values[node*fields.size()+field]);
					interleaved[static_cast<std::size_t>(node)*fields.size()+field] = value;
					output << ' ' << value;
				}
				output << '\n';
			}
			view.Restore();
			output.close();
			metadata.close();
			if (!output || !metadata)
				throw std::runtime_error("cannot write configured PDE output");
			std::vector<iga::VtkPointArray> arrays;
			for (std::size_t field = 0; field < fields.size(); ++field) {
				std::vector<double> field_values(static_cast<std::size_t>(nodes));
				for (std::size_t node = 0; node < static_cast<std::size_t>(nodes); ++node)
					field_values[node] = interleaved[node*fields.size()+field];
				arrays.push_back({fields[field], 1, std::move(field_values)});
			}
			auto derived = iga::ComputePhysiologyPointArrays(physiology, fields, interleaved);
			arrays.insert(arrays.end(), std::make_move_iterator(derived.begin()),
				std::make_move_iterator(derived.end()));
			if (visualization_format == iga::VisualizationFormat::Vtu) {
				RequireRegularOutput(vtk_path);
				iga::WriteVtu(mesh_path, vtk_path, arrays, physical_time);
			} else {
				if (!vtkhdf) throw std::runtime_error("VTKHDF writer is unavailable");
				vtkhdf->Append(physical_time, arrays);
			}
		}
	});
	objects.Close(PETSC_COMM_WORLD);
}

void WriteTransportCheckpoint(Vec state, const fs::path& prefix,
	const iga::TransportCheckpointMetadata& metadata, int rank)
{
	iga::PhaseScope output_phase(iga::ProfilePhase::Output);
	fs::path state_path;
	std::string text;
	iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint write preparation", [&] {
		state_path = iga::TransportCheckpointStatePath(prefix);
		text = iga::SerializeTransportCheckpointMetadata(metadata);
	});
	iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "transport checkpoint write agreement", text);
	iga::WritePetscCheckpointVector(state, PETSC_COMM_WORLD, state_path);
	iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint metadata write", [&] {
		if (rank != 0) return;
		const auto path = iga::TransportCheckpointMetadataPath(prefix);
		RequireRegularOutput(path);
		std::ofstream output(path);
		output << text;
		output.close();
		if (!output) throw std::runtime_error("cannot write transport checkpoint metadata: "+path.string());
	});
}

void ReadTransportCheckpoint(Vec state, const fs::path& prefix,
	const iga::TransportCheckpointMetadata& metadata)
{
	fs::path state_path;
	iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint path", [&] {
		if (metadata.state_format != "petsc_binary")
			throw std::runtime_error("CPU transport restart requires petsc_binary state");
		state_path = metadata.state_file;
		if (state_path.is_relative())
			state_path = iga::TransportCheckpointMetadataPath(prefix).parent_path()/state_path;
	});
	iga::ReadReplicatedPetscCheckpointVector(state, PETSC_COMM_WORLD, state_path);
}

} // namespace

int main(int argc, char** argv)
{
	std::exception_ptr tracking_error;
	try { EnableMemoryTracking(argc, argv); }
	catch (...) { tracking_error = std::current_exception(); }
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA configured PDE solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int ranks = 1;
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	int status = 0;
	try {
		iga::PhaseScope input_phase(iga::ProfilePhase::Input);
		iga::RequireExecutionResources(PETSC_COMM_WORLD, &std::cout);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport allocation tracking", [&] {
			if (tracking_error) std::rethrow_exception(tracking_error);
		});
		ReturnErrors errors;
		errors.Enable();
		TransportOptions options;
		std::string controls, database_fingerprint;
		std::unique_ptr<iga::Database> database_owner;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport arguments", [&] {
			options = ParseOptions(argc, argv);
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << options.system.size() << ':' << options.system << ' '
				<< options.output_every << ' ' << options.checkpoint_every << ' ' << options.stop_after_step << ' '
				<< static_cast<int>(options.visualization_format) << ' ' << !options.output.empty() << ' '
				<< !options.velocity.empty() << ' ' << !options.checkpoint.empty() << ' '
				<< !options.restart.empty() << ' ' << !options.memory_report.empty();
			controls = text.str();
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "transport execution controls", controls);
		iga::RequireCollectivePetscOptions(PETSC_COMM_WORLD);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport database preflight", [&] {
			database_fingerprint = iga::ReadAssetFingerprint(options.database);
			database_owner = std::make_unique<iga::Database>(options.database.string());
			iga::ValidatePackedExecution(database_owner->header().ranks, database_owner->header().nodes, 1, ranks);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "transport asset database", database_fingerprint);
		std::optional<iga::DistributedMemoryRecorder> memory_owner;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport memory report setup", [&] {
			// The recorder constructor only uses local MPI inquiries.
			memory_owner.emplace(PETSC_COMM_WORLD, options.memory_report);
		});
		auto& memory = *memory_owner;
		auto& database = *database_owner;
		memory.Record("database_open");
		std::unique_ptr<iga::BezierVisualizationMesh> bezier_mesh;
		std::unique_ptr<iga::TemporalVtkHdfWriter> vtkhdf;
		const auto& case_dir = options.case_dir;
		std::optional<TransportInput> input;
		iga::AssetFileCatalog case_assets, forcing_assets;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport input catalog", [&] {
			input.emplace();
			case_assets.emplace("transport asset configuration", case_dir/"simulation_config.json");
			case_assets.emplace("transport asset mesh", case_dir/"controlmesh.vtk");
		});
		iga::RequireCollectiveAssetFiles(PETSC_COMM_WORLD, case_assets);
		auto& configuration = input->configuration;
		auto& system = input->system;
		auto& labels = input->labels;
		auto& prescribed_velocity = input->prescribed_velocity;
		auto& velocity_snapshots = input->velocity_snapshots;
		auto& velocity_source = input->velocity_source;
		auto& boundaries = input->boundaries;
		auto& initial = input->initial;
		auto& coupling_patterns = input->coupling_patterns;
		auto& run_end_step = input->run_end_step;
		auto& visualization_format = input->visualization_format;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport case input", [&] {
			configuration = iga::ReadSimulationConfiguration((case_dir/"simulation_config.json").string());
			iga::RequireFlowOnlyCoupling(configuration.coupling, "CPU configured transport");
			const auto system_name = !options.system.empty() ? options.system
				: iga::FirstLinearTransportSystem(configuration);
			system = iga::CompileLinearSystem(configuration, system_name);
			iga::ValidatePackedExecution(database.header().ranks, database.header().nodes, system.fields.size(), ranks);
			if (options.stop_after_step > system.steps)
				throw std::runtime_error("--stop-after-step exceeds configured transport steps");
			run_end_step = options.stop_after_step > 0 ? options.stop_after_step : system.steps;
			visualization_format = iga::ResolveVisualizationFormat(options.visualization_format, true);
			labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
			if (system.velocity_source == "prescribed")
				forcing_assets.emplace("transport asset velocity", !options.velocity.empty()
					? options.velocity : case_dir/"initial_velocityfield.txt");
			else {
				if (!options.velocity.empty())
					throw std::runtime_error("a positional VELOCITY file cannot override a configured snapshot_series");
				velocity_source = &iga::FindVelocitySource(configuration, system.velocity_source);
				forcing_assets.emplace("transport asset velocity manifest", case_dir/velocity_source->manifest);
			}
			for (const auto& boundary : configuration.boundaries)
				for (const auto& condition : boundary.conditions) {
					if (condition.waveform.empty()) continue;
					const auto& function = iga::FindTemporalFunction(configuration, condition.waveform);
					if (function.kind == iga::TemporalFunctionKind::PeriodicTable)
						forcing_assets.emplace("transport asset temporal " + function.name, case_dir/function.file);
				}
		});
		iga::RequireCollectiveAssetFiles(PETSC_COMM_WORLD, forcing_assets);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport forcing input", [&] {
			if (velocity_source) velocity_snapshots = iga::ReadVelocityManifest(case_dir/velocity_source->manifest);
			else prescribed_velocity = iga::ReadVelocity(forcing_assets.at("transport asset velocity").string(), database.header().nodes);
			const auto initial_configuration = iga::MaterializeBoundaryWaveforms(configuration, case_dir, 0.0);
			boundaries = iga::ResolveScalarBoundaries(initial_configuration, system, labels);
			initial = iga::InitialScalarValues(configuration, system);
			coupling_patterns = iga::BuildTransportCouplingPatterns(system, configuration);
		});
		const auto fields = system.fields.size();
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport input logging", [&] {
			if (rank == 0) {
				std::cout << "configuration=simulation_config.json system=" << system.name
					<< " fields=" << fields << " velocity_source=" << system.velocity_source
					<< " dirichlet_dofs=" << boundaries.constrained_dofs << '\n';
				iga::FlushCheckedText(std::cout);
				for (std::size_t i = 0; i < system.fields.size(); ++i)
					std::cout << "field[" << i << "]=" << system.fields[i] << '\n';
				iga::FlushCheckedText(std::cout);
			}
		});

		input_phase.Stop();
		iga::PhaseScope geometry_phase(iga::ProfilePhase::Geometry);
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, fields);
		memory.Record("required_elements_loaded");
		iga::RequireValidGeometry(assembler.elements(),
			[&assembler](const iga::Element& element) {
				return assembler.OwnsElementByMinimumNode(element);
			}, PETSC_COMM_WORLD);
		iga::RequireConfiguredScalarSurfaceFaces(
			configuration, system, assembler, PETSC_COMM_WORLD);
		geometry_phase.Stop();
		TransportPetscObjects objects;
		auto& left = objects.left;
		auto& previous = objects.previous;
		auto& forcing = objects.forcing;
		auto& current = objects.current;
		auto& next = objects.next;
		auto& rhs = objects.rhs;
		auto& solver = objects.solver;
		left = assembler.CreateMatrix(coupling_patterns.left);
		previous = assembler.CreateMatrix(coupling_patterns.previous);
		memory.Record("matrix_preallocation");
		CheckPetsc("transport MatSetOption", MatSetOption(left, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
		CheckPetsc("transport MatSetOption", MatSetOption(previous, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
		CheckPetsc("transport MatSetOption", MatSetOption(previous, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
		forcing = assembler.CreateVector();
		CheckPetsc("transport VecSet", VecSet(forcing, 0.0));

		std::vector<PetscInt> boundary_rows;
		std::optional<iga::GenericTransportMatrices> element_matrices_owner;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport assembly preparation", [&] {
			for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
				for (std::size_t field = 0; field < fields; ++field)
					if (boundaries.constrained[static_cast<std::size_t>(node)*fields+field])
						boundary_rows.push_back(static_cast<PetscInt>(node*fields+field));
			element_matrices_owner.emplace(coupling_patterns);
		});
		double assembly_seconds = 0.0;
		auto& element_matrices = *element_matrices_owner;
		bool assembly_memory_recorded = false;
		auto assemble_operators = [&](const std::vector<std::array<double, 3>>& velocity) {
			iga::PhaseScope assembly_phase(iga::ProfilePhase::Assembly);
			const auto start = std::chrono::steady_clock::now();
			CheckPetsc("transport MatZeroEntries", MatZeroEntries(left));
			CheckPetsc("transport MatZeroEntries", MatZeroEntries(previous));
			CheckPetsc("transport VecSet", VecSet(forcing, 0.0));
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport element assembly", [&] {
				for (const auto& element : assembler.elements()) {
					iga::FullCell4x4x4VolumeQuadratureProvider volume_quadrature(element);
					iga::BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
					iga::BuildGenericTransportElement(
						element, velocity, system, configuration, element_matrices,
						volume_quadrature.Rule(), surface_quadrature.Rule());
					assembler.AddElementMatrix(left, element, element_matrices.left);
					assembler.AddElementMatrix(previous, element, element_matrices.previous);
					assembler.AddElementVector(forcing, element, element_matrices.source);
				}
			});
			iga::PhaseScope finalize_phase(iga::ProfilePhase::Communication);
			iga::OwnedRowAssembler::Assemble(left, PETSC_COMM_WORLD);
			iga::OwnedRowAssembler::Assemble(previous, PETSC_COMM_WORLD);
			iga::OwnedRowAssembler::Assemble(forcing, PETSC_COMM_WORLD);
			CheckPetsc("transport MatZeroRows", MatZeroRows(left, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 1.0, nullptr, nullptr));
			CheckPetsc("transport MatZeroRows", MatZeroRows(previous, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 0.0, nullptr, nullptr));
			finalize_phase.Stop();
			if (!assembly_memory_recorded) {
				memory.Record("operator_assembly");
				assembly_memory_recorded = true;
			}
			assembly_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-start).count();
		};
		auto velocity_at = [&](double time) {
			iga::PhaseScope input_phase(iga::ProfilePhase::Input);
			iga::VelocityInterpolation interpolation;
			iga::AssetFileCatalog snapshots;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport velocity selection", [&] {
				interpolation = iga::ResolveVelocityInterpolation(velocity_snapshots, time, velocity_source->out_of_range);
				snapshots.emplace("transport asset lower snapshot", case_dir/velocity_snapshots[interpolation.lower].file);
				if (interpolation.lower != interpolation.upper)
					snapshots.emplace("transport asset upper snapshot", case_dir/velocity_snapshots[interpolation.upper].file);
			});
			iga::RequireCollectiveAssetFiles(PETSC_COMM_WORLD, snapshots);
			std::vector<std::array<double, 3>> lower;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport velocity input", [&] {
				lower = iga::ReadVelocity(
					(case_dir/velocity_snapshots[interpolation.lower].file).string(), database.header().nodes);
				if (interpolation.lower == interpolation.upper) return;
				const auto upper = iga::ReadVelocity(
					(case_dir/velocity_snapshots[interpolation.upper].file).string(), database.header().nodes);
				for (std::size_t node = 0; node < lower.size(); ++node)
					for (int component = 0; component < 3; ++component)
						lower[node][component] = (1.0-interpolation.upper_weight)*lower[node][component]
							+ interpolation.upper_weight*upper[node][component];
			});
			return lower;
		};
		if (!velocity_source) assemble_operators(prescribed_velocity);

		current = assembler.CreateVector();
		next = assembler.CreateVector();
		rhs = assembler.CreateVector();
		CheckPetsc("transport VecSet", VecSet(current, 0.0));
		CheckPetsc("transport VecSet", VecSet(next, 0.0));
		CheckPetsc("transport VecSet", VecSet(rhs, 0.0));
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport initial state", [&] {
			WriteArray view;
			view.Acquire(current);
			auto* current_array = view.Data();
			for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
				for (std::size_t field = 0; field < fields; ++field) {
					const auto local = static_cast<std::size_t>(node-assembler.node_begin())*fields+field;
					const auto global = static_cast<std::size_t>(node)*fields+field;
					current_array[local] = boundaries.constrained[global] ? boundaries.value[global] : initial[field];
				}
			view.Restore();
		});
		int start_step = 0;
		if (!options.restart.empty()) {
			iga::TransportCheckpointMetadata metadata;
			std::string metadata_text;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint metadata", [&] {
				(void)iga::ReadAssetFingerprint(iga::TransportCheckpointMetadataPath(options.restart));
				metadata = iga::ReadTransportCheckpointMetadata(options.restart);
				iga::ValidateTransportCheckpoint(metadata, database.header().nodes,
					system.fields, system.name, system.velocity_source, system.steps, system.dt);
				metadata_text = iga::SerializeTransportCheckpointMetadata(metadata);
			});
			iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "transport checkpoint metadata agreement", metadata_text);
			ReadTransportCheckpoint(current, options.restart, metadata);
			start_step = metadata.completed_step;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport restart logging", [&] {
				if (rank == 0) std::cout << "restart=" << options.restart.string()
					<< " completed_step=" << start_step
					<< " physical_time=" << metadata.physical_time << '\n';
				iga::FlushCheckedText(std::cout);
			});
		}
		if (!options.output.empty()
			&& visualization_format == iga::VisualizationFormat::BezierVtkHdf) {
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport visualization initialization", [&] {
				if (rank != 0) return;
				bezier_mesh = std::make_unique<iga::BezierVisualizationMesh>(
					iga::BuildBezierVisualizationMesh(database, false));
				const auto report = iga::BezierGeometryReportPath(options.output);
				RequireRegularOutput(report);
				RequireRegularOutput(iga::VtkHdfPath(options.output));
				iga::WriteBezierGeometryReport(report, bezier_mesh->validation);
				iga::RequireValidBezierGeometry(bezier_mesh->validation);
				vtkhdf = std::make_unique<iga::TemporalVtkHdfWriter>(
					iga::VtkHdfPath(options.output), *bezier_mesh,
					!options.restart.empty());
				std::cout << "bezier_geometry_points=" << bezier_mesh->points.size()
					<< " local_point_references=" << bezier_mesh->validation.local_points
					<< " geometry_report=" << report.string()
					<< " vtkhdf=" << vtkhdf->path().string() << '\n';
				iga::FlushCheckedText(std::cout);
			});
		}

		CheckPetsc("transport KSPCreate", KSPCreate(PETSC_COMM_WORLD, &solver));
		if (!velocity_source) CheckPetsc("transport KSPSetOperators", KSPSetOperators(solver, left, left));
		CheckPetsc("transport KSPSetType", KSPSetType(solver, KSPGMRES));
		CheckPetsc("transport KSPGMRESSetRestart", KSPGMRESSetRestart(solver, 50));
		CheckPetsc("transport KSPSetTolerances", KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000));
		PC preconditioner = nullptr;
		CheckPetsc("transport KSPGetPC", KSPGetPC(solver, &preconditioner));
		CheckPetsc("transport PCSetType", PCSetType(preconditioner, PCBJACOBI));
		CheckPetsc("transport KSPSetFromOptions", KSPSetFromOptions(solver));
		if (!velocity_source) {
			SetupTransportKsp(solver, left);
			memory.Record("ksp_setup");
		}

		PetscInt total_iterations = 0;
		bool setup_memory_recorded = !velocity_source;
		const auto solve_start = std::chrono::steady_clock::now();
		double linear_seconds = 0.0;
		std::vector<std::pair<double, fs::path>> vtk_snapshots;
		auto write_output = [&](int step, bool final) {
			fs::path text_path, vtk_path, mesh_path;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport output paths", [&] {
				text_path = final ? options.output : iga::TimeIndexedPath(options.output, step);
				vtk_path = final ? iga::VtuFinalPath(options.output) : iga::VtuStepPath(options.output, step);
				mesh_path = case_dir/"controlmesh.vtk";
			});
			WriteTransportOutput(current, database.header().nodes, system.fields,
				text_path, mesh_path, vtk_path, step*system.dt,
				configuration.physiology, rank, visualization_format, vtkhdf.get());
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport output bookkeeping", [&] {
				if (!final && visualization_format == iga::VisualizationFormat::Vtu)
					vtk_snapshots.push_back({step*system.dt, vtk_path});
			});
		};
		if (options.output_every > 0) {
			write_output(start_step, false);
		}
		for (int step = start_step; step < run_end_step; ++step) {
			if (velocity_source) {
				const auto velocity = velocity_at((step+1)*system.dt);
				assemble_operators(velocity);
				CheckPetsc("transport KSPSetOperators", KSPSetOperators(solver, left, left));
				SetupTransportKsp(solver, left);
				if (!setup_memory_recorded) {
					memory.Record("ksp_setup");
					setup_memory_recorded = true;
				}
			}
			iga::ResolvedScalarBoundaries step_boundaries;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport step input", [&] {
				const auto step_configuration = iga::MaterializeBoundaryWaveforms(
					configuration, case_dir, (step+1)*system.dt);
				step_boundaries = iga::ResolveScalarBoundaries(step_configuration, system, labels);
				for (const auto value : step_boundaries.value)
					if (!std::isfinite(value)) throw std::runtime_error("transport boundary value is not finite");
			});
			CheckPetsc("transport MatMult", MatMult(previous, current, rhs));
			CheckPetsc("transport VecAXPY", VecAXPY(rhs, 1.0, forcing));
			std::vector<PetscScalar> boundary_values;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport boundary values", [&] {
				boundary_values.reserve(boundary_rows.size());
				for (auto row : boundary_rows)
					boundary_values.push_back(step_boundaries.value[static_cast<std::size_t>(row)]);
				if (VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES))
					throw std::runtime_error("cannot insert transport boundary values");
			});
			iga::OwnedRowAssembler::Assemble(rhs, PETSC_COMM_WORLD);
			if (step > 0) CheckPetsc("transport VecCopy", VecCopy(current, next));
			CheckPetsc("transport KSPSetInitialGuessNonzero", KSPSetInitialGuessNonzero(solver, step > 0 ? PETSC_TRUE : PETSC_FALSE));
			const auto linear_start = std::chrono::steady_clock::now();
			SetupTransportKsp(solver, left);
			{
				iga::PhaseScope linear_phase(iga::ProfilePhase::LinearSolve);
				CheckPetsc("transport linear solve", KSPSolve(solver, rhs, next));
			}
			linear_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-linear_start).count();
			KSPConvergedReason reason;
			PetscInt iterations = 0;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport convergence", [&] {
				if (KSPGetConvergedReason(solver, &reason) || KSPGetIterationNumber(solver, &iterations))
					throw std::runtime_error("cannot query transport KSP status");
				if (reason <= 0) throw std::runtime_error("KSP did not converge at step " + std::to_string(step));
			});
			total_iterations += iterations;
			CheckPetsc("transport VecSwap", VecSwap(current, next));
			const auto completed_step = step+1;
			if (options.output_every > 0
				&& completed_step%options.output_every == 0) {
				write_output(completed_step, false);
			}
			if (!options.checkpoint.empty()
				&& (completed_step == run_end_step
					|| (options.checkpoint_every > 0
						&& completed_step%options.checkpoint_every == 0))) {
				iga::TransportCheckpointMetadata metadata;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint metadata preparation", [&] {
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
				});
				WriteTransportCheckpoint(current, options.checkpoint, metadata, rank);
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport checkpoint logging", [&] {
					if (rank == 0) std::cout << "checkpoint=" << options.checkpoint.string()
						<< " completed_step=" << completed_step << '\n';
					iga::FlushCheckedText(std::cout);
				});
			}
		}
		memory.Record("linear_solve");
		const auto solve_end = std::chrono::steady_clock::now();

		if (!options.output.empty()) {
			write_output(run_end_step, true);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport output index", [&] {
				if (rank == 0) {
					if (visualization_format == iga::VisualizationFormat::Vtu) {
						if (vtk_snapshots.empty())
							vtk_snapshots.push_back({run_end_step*system.dt,
								iga::VtuFinalPath(options.output)});
						RequireRegularOutput(iga::PvdPath(options.output));
						iga::WritePvd(iga::PvdPath(options.output), vtk_snapshots);
					}
					RequireRegularOutput(options.output.parent_path()/"physiology_fields.json");
					iga::WritePhysiologyManifest(
						options.output.parent_path()/"physiology_fields.json",
						configuration.physiology, system.fields);
				}
			});
		}
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport visualization close", [&] {
			if (rank == 0 && vtkhdf) vtkhdf->Close();
		});
		PetscReal norm = 0.0;
		CheckPetsc("transport VecNorm", VecNorm(current, NORM_2, &norm));
		CheckPetsc("transport KSPDestroy", KSPDestroy(&solver));
		CheckPetsc("transport VecDestroy", VecDestroy(&rhs));
		CheckPetsc("transport VecDestroy", VecDestroy(&next));
		CheckPetsc("transport VecDestroy", VecDestroy(&current));
		CheckPetsc("transport VecDestroy", VecDestroy(&forcing));
		CheckPetsc("transport MatDestroy", MatDestroy(&previous));
		CheckPetsc("transport MatDestroy", MatDestroy(&left));
		memory.Close();
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport final summary", [&] {
			if (!std::isfinite(norm)) throw std::runtime_error("transport state norm is not finite");
			if (rank == 0) std::cout << "iga_solve system=" << system.name
				<< " steps=" << system.steps << " run_end_step=" << run_end_step
				<< " assembly_s=" << assembly_seconds
				<< " solve_s=" << linear_seconds
				<< " time_loop_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
				<< " total_iterations=" << total_iterations << " final_l2=" << norm << '\n';
			iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	try {
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "transport profile logging", [&] {
			iga::CurrentPhaseProfile().Write(std::cout, rank, ranks, global_status);
			iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		global_status = 1;
	}
	PetscFinalize();
	return global_status;
}

#include <cctype>
#include "CheckedText.hpp"
#include "ExecutionResources.hpp"
#include "BoundarySupport.hpp"
#include "CaseInput.hpp"
#include "CollectiveAssetInput.hpp"
#include "CouplingHistory.hpp"
#include "FlowCheckpoint.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "OutletCheckpoint.hpp"
#include "PetscCheckpointRead.hpp"
#include "PetscCheckpointWrite.hpp"
#include "PetscGather.hpp"
#include "PetscBezierVisualization.hpp"
#include "TemporalFunction.hpp"
#include "ThreeDVcaCoupling.hpp"
#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"
#include "VcaCheckpoint.hpp"
#include "VelocitySeries.hpp"
#include "TemporalVtkHdf.hpp"
#include "VtkOutput.hpp"

#include <petscksp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <type_traits>
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
	bool parallel_output = false;
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
	double nonlinear_relative_tolerance = 1e-5;
	double nonlinear_absolute_tolerance = 1e-10;
	double mass_relative_tolerance = 1e-3;
};

// Construct defaults and read local data within coordinated stages, before
// any runtime constructor enters PETSc collectives.
struct FlowCaseInput {
	iga::LabeledHexMesh mesh;
	std::set<std::int32_t> wall_trace_basis;
	std::vector<std::array<double, 3>> boundary_velocity;
	iga::SimulationConfiguration configuration;
	iga::ResolvedBoundaryConditions boundaries;
	std::vector<iga::OutletModelState> outlet_models;
	iga::NavierStokesParameters parameters{1.0, 0.1, 0.0};
	bool configured = false;
	bool transient = false;
	std::string boundary_config;
	std::unique_ptr<iga::VcaExternalCircuit> vca_circuit;
	iga::CompiledLinearSystem vca_transport_system;
	bool vca_has_transport = false;
	int physical_steps = 1;
	int run_end_step = 1;
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
};

struct FlowStepInput {
	iga::SimulationConfiguration configuration;
	iga::VascularInletState inlet;
};

FlowStepInput PrepareFlowStepInput(MPI_Comm communicator, const char* stage,
	bool configured, bool transient, const iga::SimulationConfiguration& configuration,
	const fs::path& case_directory, double physical_time, double inlet_time,
	iga::VcaExternalCircuit* circuit, const iga::CompiledLinearSystem* transport,
	double reference_inlet_flow)
{
	static_assert(std::is_nothrow_move_constructible<FlowStepInput>::value,
		"step input must leave its coordinated stage without another allocation");
	std::optional<FlowStepInput> input;
	iga::CollectiveLocalStage(communicator, stage, [&] {
		input.emplace();
		if (configured)
			input->configuration = transient
				? iga::MaterializeBoundaryWaveforms(configuration, case_directory, physical_time)
				: configuration;
		if (circuit) {
			input->inlet = circuit->InletState(inlet_time);
			iga::ApplyThreeDVascularInlet(input->configuration,
				iga::FirstNavierStokesSystem(input->configuration), input->inlet, reference_inlet_flow);
			if (transport)
				iga::ApplyThreeDVascularSpeciesInlet(input->configuration, *transport, input->inlet);
		}
	});
	return std::move(*input);
}

iga::VascularStepResult PrepareFlowVcaResult(MPI_Comm communicator, double physical_time,
	double dt, const iga::VascularInletState& inlet, const iga::ThreeDVascularPortDefinition& definition,
	const iga::FlowPortMeasurements& ports, double flow_epsilon)
{
	static_assert(std::is_nothrow_move_constructible<iga::VascularStepResult>::value,
		"VCA result must leave its coordinated stage without another allocation");
	std::optional<iga::VascularStepResult> result;
	iga::CollectiveLocalStage(communicator, "flow VCA result preparation", [&] {
		result.emplace(iga::BuildThreeDFlowPortResult(physical_time, dt, inlet,
			definition, ports.flows, ports.pressures));
		for (auto& outlet : result->outlets) {
			const auto flux = ports.species_fluxes.find(outlet.outlet_id);
			if (flux == ports.species_fluxes.end()) continue;
			outlet.species_flux = flux->second;
			outlet.average_valid = std::abs(outlet.flow_m3_s) > flow_epsilon;
			const auto concentration = ports.species_concentrations.find(outlet.outlet_id);
			if (concentration != ports.species_concentrations.end())
				outlet.flux_weighted_concentration = concentration->second;
		}
	});
	return std::move(*result);
}

void AdvanceFlowVcaCircuit(MPI_Comm communicator, iga::VcaExternalCircuit& circuit,
	iga::CouplingHistoryWriter& history, const iga::VascularStepResult& result, double flow_epsilon)
{
	iga::CollectiveLocalStage(communicator, "flow VCA circuit advance", [&] {
		const auto venous = iga::AggregateVascularOutlets(result.outlets, flow_epsilon);
		const auto report = circuit.Advance(venous, result.dt_s, result.time_s);
		history.Add(result, venous, report);
	});
}

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

double ParsePositiveFiniteDouble(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	double value = 0.0;
	try {
		value = std::stod(text, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(option+" requires a finite positive value");
	}
	if (used != text.size() || !std::isfinite(value) || value <= 0.0)
		throw std::runtime_error(option+" requires a finite positive value");
	return value;
}

FlowOptions ParseOptions(int argc, char** argv)
{
	if (argc < 3) throw std::runtime_error(
		"usage: iga_navier_stokes DATABASE.ntiga CASE_DIR [MAX_NEWTON] [OUTPUT] "
		"[--max-newton N] [--output PATH] [--output-every N] "
		"[--checkpoint PREFIX] [--checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N] [--nonlinear-rtol R] [--nonlinear-atol A] [--mass-rtol R] "
		"[--visualization-format auto|vtu|vtkhdf|pvtu]");
	FlowOptions options;
	options.database = argv[1];
	options.case_dir = argv[2];
	int positional = 0;
	for (int i = 3; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument.size() > 1 && argument[0] == '-' && std::isalpha(static_cast<unsigned char>(argument[1]))) {
			// PETSc parsed these keys during initialization; keep their values
			// out of the legacy positional application interface.
			if (i+1 < argc) {
				const std::string next(argv[i+1]);
				const bool next_key = next.size() > 1 && next[0] == '-'
					&& (next[1] == '-' || std::isalpha(static_cast<unsigned char>(next[1])));
				if (!next_key) ++i;
			}
			continue;
		}
		if (argument.rfind("--", 0) != 0) {
			if (positional == 0) options.max_newton = ParsePositiveInteger(argument, "MAX_NEWTON");
			else if (positional == 1) options.output = argument;
			else throw std::runtime_error("too many positional arguments");
			++positional;
			continue;
		}
		if (i+1 >= argc) throw std::runtime_error(argument+" requires a value");
		const std::string value(argv[++i]);
		if (argument == "--max-newton") options.max_newton = ParsePositiveInteger(value, argument);
		else if (argument == "--output") options.output = value;
		else if (argument == "--output-every") options.output_every = ParsePositiveInteger(value, argument);
		else if (argument == "--checkpoint") options.checkpoint = value;
		else if (argument == "--checkpoint-every") options.checkpoint_every = ParsePositiveInteger(value, argument);
		else if (argument == "--restart") options.restart = value;
		else if (argument == "--stop-after-step") options.stop_after_step = ParsePositiveInteger(value, argument);
		else if (argument == "--nonlinear-rtol")
			options.nonlinear_relative_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--nonlinear-atol")
			options.nonlinear_absolute_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--mass-rtol")
			options.mass_relative_tolerance = ParsePositiveFiniteDouble(value, argument);
		else if (argument == "--visualization-format") {
			options.parallel_output = value == "pvtu";
			options.visualization_format = options.parallel_output ? iga::VisualizationFormat::Vtu : iga::ParseVisualizationFormat(value);
		}
		else throw std::runtime_error("unknown option: "+argument);
		if (PetscOptionsClearValue(nullptr, argument.c_str()))
			throw std::runtime_error("cannot consume application option: "+argument);
	}
	if (options.parallel_output && options.output.empty())
		throw std::runtime_error("pvtu requires --output or legacy OUTPUT");
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output or legacy OUTPUT");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

void RequireRegularOutput(const fs::path& path)
{
	if (fs::exists(path) && !fs::is_regular_file(path))
		throw std::runtime_error("output target is not a regular file: "+path.string());
}

void WriteFlowOutput(Vec state, std::uint64_t nodes, const fs::path& path,
	const fs::path& mesh_path, const fs::path& vtk_path, double physical_time, int rank,
	iga::VisualizationFormat visualization_format,
	iga::TemporalVtkHdfWriter* vtkhdf)
{
	iga::PhaseScope output_phase(iga::ProfilePhase::Output);
	const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(state));
	iga::CollectiveLocalStage(communicator, "flow output layout", [&] {
		PetscInt rows = 0;
		if (VecGetSize(state, &rows)) throw std::runtime_error("cannot query flow output rows");
		if (rows < 0 || nodes > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max())/4
			|| static_cast<std::uint64_t>(rows) != 4*nodes)
			throw std::runtime_error("flow output state size differs from database nodes");
	});
	iga::PetscGatherObjects objects;
	iga::PhaseScope gather_phase(iga::ProfilePhase::Communication);
	iga::RequireCollectivePetscSuccess(communicator, "flow output gather create",
		VecScatterCreateToZero(state, &objects.scatter, &objects.all));
	iga::RequireCollectivePetscSuccess(communicator, "flow output gather begin",
		VecScatterBegin(objects.scatter, state, objects.all, INSERT_VALUES, SCATTER_FORWARD));
	iga::RequireCollectivePetscSuccess(communicator, "flow output gather end",
		VecScatterEnd(objects.scatter, state, objects.all, INSERT_VALUES, SCATTER_FORWARD));
	gather_phase.Stop();
	iga::CollectiveLocalStage(communicator, "flow field output", [&] {
		if (rank == 0) {
			iga::PetscReadArray view;
			view.Acquire(objects.all);
			const auto* values = view.Data();
			std::vector<double> velocity(3*static_cast<std::size_t>(nodes));
			std::vector<double> pressure(static_cast<std::size_t>(nodes));
			RequireRegularOutput(path);
			RequireRegularOutput(path.string()+".pressure");
			std::ofstream output(path);
			std::ofstream pressure_output(path.string()+".pressure");
			if (!output || !pressure_output) throw std::runtime_error("cannot create Navier-Stokes output");
			output.precision(17);
			pressure_output.precision(17);
			for (std::uint64_t node = 0; node < nodes; ++node) {
				for (int component = 0; component < 3; ++component)
					velocity[3*static_cast<std::size_t>(node)+component]
						= PetscRealPart(values[4*node+component]);
				pressure[static_cast<std::size_t>(node)] = PetscRealPart(values[4*node+3]);
				output << velocity[3*static_cast<std::size_t>(node)] << ' '
					<< velocity[3*static_cast<std::size_t>(node)+1] << ' '
					<< velocity[3*static_cast<std::size_t>(node)+2] << '\n';
				pressure_output << pressure[static_cast<std::size_t>(node)] << '\n';
			}
			view.Restore();
			output.close();
			pressure_output.close();
			if (!output || !pressure_output) throw std::runtime_error("cannot write Navier-Stokes output");
			std::vector<iga::VtkPointArray> arrays;
			arrays.reserve(2);
			arrays.push_back({"velocity", 3, std::move(velocity)});
			arrays.push_back({"pressure", 1, std::move(pressure)});
			if (visualization_format == iga::VisualizationFormat::Vtu) {
				RequireRegularOutput(vtk_path);
				iga::WriteVtu(mesh_path, vtk_path, arrays, physical_time);
			} else {
				if (!vtkhdf) throw std::runtime_error("VTKHDF writer is unavailable");
				vtkhdf->Append(physical_time, arrays);
			}
		}
	});
	objects.Close(communicator);
}

void WriteCheckpointMetadataText(const fs::path& path, const std::string& text)
{
	std::ofstream output(path);
	output << text;
	output.close();
	if (!output) throw std::runtime_error("cannot write checkpoint metadata: "+path.string());
}

void WriteCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata, int rank)
{
	iga::PhaseScope output_phase(iga::ProfilePhase::Output);
	const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(state));
	fs::path state_path;
	std::string metadata_text;
	iga::CollectiveLocalStage(communicator, "flow checkpoint write preparation", [&] {
		state_path = iga::FlowCheckpointStatePath(prefix);
		metadata_text = iga::SerializeFlowCheckpointMetadata(metadata);
	});
	iga::RequireCollectiveSameText(communicator, "flow checkpoint write agreement", metadata_text);
	iga::WritePetscCheckpointVector(state, communicator, state_path);
	iga::CollectiveLocalStage(communicator, "flow checkpoint metadata write", [&] {
		if (rank == 0) WriteCheckpointMetadataText(iga::FlowCheckpointMetadataPath(prefix), metadata_text);
	});
}

void ReadCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata)
{
	const auto communicator = PetscObjectComm(reinterpret_cast<PetscObject>(state));
	fs::path path;
	iga::CollectiveLocalStage(communicator, "flow checkpoint path", [&] {
		if (metadata.state_format != "petsc_binary")
			throw std::runtime_error("CPU flow restart requires petsc_binary checkpoint state");
		path = metadata.state_file;
		if (path.is_relative()) path = iga::FlowCheckpointMetadataPath(prefix).parent_path()/path;
	});
	iga::ReadPetscCheckpointVector(state, communicator, path);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr,
		"TubularFlowIGA stabilized steady/transient Navier-Stokes solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int ranks = 1;
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	iga::CurrentPhaseProfile().EnableFromEnvironment();
	int status = 0;
	try {
		iga::PhaseScope input_phase(iga::ProfilePhase::Input);
		iga::RequireExecutionResources(PETSC_COMM_WORLD, &std::cout);
		FlowOptions options;
		std::string controls, database_fingerprint;
		std::unique_ptr<iga::Database> database_owner;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow arguments", [&] {
			options = ParseOptions(argc, argv);
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10)
				<< options.max_newton << ' ' << options.output_every << ' '
				<< options.checkpoint_every << ' ' << options.stop_after_step << ' '
				<< static_cast<int>(options.visualization_format) << ' ' << options.parallel_output << ' '
				<< options.nonlinear_relative_tolerance << ' ' << options.nonlinear_absolute_tolerance << ' '
				<< options.mass_relative_tolerance << ' ' << !options.output.empty() << ' '
				<< !options.checkpoint.empty() << ' ' << !options.restart.empty();
			controls = text.str();
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "flow execution controls", controls);
		iga::RequireCollectivePetscOptions(PETSC_COMM_WORLD);
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow database preflight", [&] {
			// Reject special files before Database can open a blocking stream.
			database_fingerprint = iga::ReadAssetFingerprint(options.database);
			database_owner = std::make_unique<iga::Database>(options.database.string());
			iga::ValidatePackedExecution(database_owner->header().ranks, database_owner->header().nodes, 4, ranks);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "flow asset database", database_fingerprint);
		auto& database = *database_owner;
		std::unique_ptr<iga::BezierVisualizationMesh> bezier_mesh;
		std::unique_ptr<iga::TemporalVtkHdfWriter> vtkhdf;
		std::optional<FlowCaseInput> input;
		iga::AssetFileCatalog case_assets;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow input catalog", [&] {
			input.emplace();
			input->configured = fs::exists(options.case_dir/"simulation_config.json");
			case_assets.emplace("flow asset mesh", options.case_dir/"controlmesh.vtk");
			case_assets.emplace("flow asset velocity", options.case_dir/"initial_velocityfield.txt");
			if (input->configured)
				case_assets.emplace("flow asset configuration", options.case_dir/"simulation_config.json");
			else {
				case_assets.emplace("flow asset parameters", options.case_dir/"simulation_parameter.txt");
				if (fs::exists(options.case_dir/"case_config.json"))
					case_assets.emplace("flow asset legacy configuration", options.case_dir/"case_config.json");
			}
		});
		iga::RequireCollectiveAssetFiles(PETSC_COMM_WORLD, case_assets);
		auto& mesh = input->mesh;
		auto& wall_trace_basis = input->wall_trace_basis;
		auto& boundary_velocity = input->boundary_velocity;
		auto& configuration = input->configuration;
		auto& boundaries = input->boundaries;
		auto& outlet_models = input->outlet_models;
		auto& parameters = input->parameters;
		auto& configured = input->configured;
		auto& transient = input->transient;
		auto& boundary_config = input->boundary_config;
		auto& vca_circuit = input->vca_circuit;
		auto& vca_transport_system = input->vca_transport_system;
		auto& vca_has_transport = input->vca_has_transport;
		auto& physical_steps = input->physical_steps;
		auto& run_end_step = input->run_end_step;
		auto& visualization_format = input->visualization_format;
		const auto& labels = mesh.labels;
		iga::AssetFileCatalog waveform_assets;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow case input", [&] {
			mesh = iga::ReadLabeledHexMesh((options.case_dir/"controlmesh.vtk").string(),
				database.header().nodes, database.header().elements);
			wall_trace_basis = iga::WallTraceBasis(database, mesh);
			boundary_velocity = iga::ReadVelocity(
				(options.case_dir/"initial_velocityfield.txt").string(), database.header().nodes);
			if (configured) {
				configuration = iga::ReadSimulationConfiguration(
					(options.case_dir/"simulation_config.json").string());
				transient = iga::FirstNavierStokesSystem(configuration).time_integration == "backward_euler";
				// MaterializeBoundaryWaveforms evaluates referenced boundary functions.
				if (transient)
					for (const auto& boundary : configuration.boundaries)
						for (const auto& condition : boundary.conditions) {
							if (condition.waveform.empty()) continue;
							const auto& function = iga::FindTemporalFunction(configuration, condition.waveform);
							if (function.kind == iga::TemporalFunctionKind::PeriodicTable)
								waveform_assets.emplace("flow asset temporal " + function.name,
									options.case_dir/function.file);
						}
			}
		});
		iga::RequireCollectiveAssetFiles(PETSC_COMM_WORLD, waveform_assets);
		std::unique_ptr<iga::TransientTransportRuntime> vca_transport;
		iga::VcaCheckpointIdentity vca_checkpoint_identity;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow boundary input", [&] {
			if (configured) {
				const auto& flow = iga::FirstNavierStokesSystem(configuration);
				configured = true;
				transient = flow.time_integration == "backward_euler";
				if (configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly) {
					iga::RequireThreeDVascularPorts(configuration.coupling, "CPU 3D VCA flow bridge");
					if (configuration.coupling.external_circuit.reservoir.species.empty())
						iga::RequireThreeDFlowOnlyCircuit(configuration.coupling);
					else {
						vca_transport_system = iga::RequireThreeDVcaTransportSystem(configuration);
						vca_has_transport = true;
					}
					if (!transient)
						throw std::runtime_error("CPU 3D VCA flow bridge requires backward_euler Navier-Stokes");
					vca_circuit = std::make_unique<iga::VcaExternalCircuit>(configuration.coupling);
				}
				parameters = {flow.density, flow.viscosity, transient ? configuration.time.dt : 0.0};
				const auto waveform = transient
					? iga::MaterializeBoundaryWaveforms(configuration, options.case_dir.string(), 0.0)
					: configuration;
				outlet_models = iga::InitializeOutletModels(configuration, flow);
				const auto initial = iga::MaterializeOutletPressures(waveform, outlet_models);
				boundaries = iga::ResolveFlowBoundaries(initial, iga::FirstNavierStokesSystem(initial),
					labels, boundary_velocity);
				boundary_config = "simulation_config.json";
			} else {
				const auto transport = iga::ReadTransportParameters(
					(options.case_dir/"simulation_parameter.txt").string());
				const auto case_config = iga::ReadCaseConfiguration((options.case_dir/"case_config.json").string());
				boundaries = iga::ResolveBoundaryConditions(case_config, labels, boundary_velocity, transport);
				boundary_config = case_config.present ? "case_config.json" : "legacy-defaults";
			}
			physical_steps = transient ? configuration.time.steps : 1;
			if (options.stop_after_step > physical_steps)
				throw std::runtime_error("--stop-after-step exceeds configured physical steps");
			run_end_step = options.stop_after_step > 0 ? options.stop_after_step : physical_steps;
			visualization_format = iga::ResolveVisualizationFormat(
				options.visualization_format, transient);
			if (vca_circuit && !vca_has_transport
				&& (!options.restart.empty() || !options.checkpoint.empty()))
				throw std::runtime_error("VCA flow-only checkpoint/restart requires a transport state and is unavailable");
			if (!transient)
				for (const auto& model : outlet_models)
					if (model.kind != iga::FieldBoundaryKind::Resistance)
						throw std::runtime_error("RC/RCR outlets require backward_euler flow");
			if (!transient && (!options.restart.empty() || !options.checkpoint.empty()
				|| options.output_every > 0 || options.checkpoint_every > 0))
				throw std::runtime_error("restart, checkpoint, and time-indexed output require transient flow");
			if (rank == 0) std::cout << "boundary_config=" << boundary_config
				<< " viscosity=" << parameters.dynamic_viscosity << " density=" << parameters.density
				<< " time_integration=" << (transient ? "backward_euler" : "steady")
				<< " dt=" << parameters.dt << " steps=" << physical_steps
				<< " run_end_step=" << run_end_step << " velocity_nodes=" << boundaries.velocity_nodes
				<< " pressure_nodes=" << boundaries.pressure_nodes
				<< " wall_trace_velocity_nodes=" << wall_trace_basis.size() << '\n';
			iga::FlushCheckedText(std::cout);
		});

		std::string flow_solver_prefix, transport_solver_prefix;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow solver prefixes", [&] {
			flow_solver_prefix = iga::PetscDomainOptionsPrefix(configured ? iga::FirstNavierStokesSystem(configuration).name : "flow", "flow");
			if (vca_has_transport) transport_solver_prefix = iga::PetscDomainOptionsPrefix(vca_transport_system.name, "transport");
		});
		input_phase.Stop();
		iga::PhaseScope geometry_phase(iga::ProfilePhase::Geometry);
		iga::TransientFlowRuntime flow(database, PETSC_COMM_WORLD, configured, transient,
			parameters, boundaries, labels, boundary_velocity, wall_trace_basis,
			std::move(outlet_models), {}, {}, 0.0, flow_solver_prefix);
		iga::RequireValidGeometry(flow.Elements(), rank, PETSC_COMM_WORLD);
		geometry_phase.Stop();
		if (vca_circuit) {
			std::map<int, long long> port_faces;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow VCA face catalog", [&] {
				port_faces.emplace(configuration.coupling.three_d_ports.inlet_label, 0);
				for (const auto label : configuration.coupling.three_d_ports.outlet_labels)
					port_faces.emplace(label, 0);
				for (const auto& element : flow.OwnedElements())
					for (const auto label : element.boundary_labels) {
						auto found = port_faces.find(label);
						if (found != port_faces.end()) ++found->second;
					}
			});
			for (auto& item : port_faces) {
				long long global_faces = 0;
				MPI_Allreduce(&item.second, &global_faces, 1, MPI_LONG_LONG, MPI_SUM,
					PETSC_COMM_WORLD);
				if (global_faces == 0) throw std::runtime_error("VCA port label "
					+std::to_string(item.first)+" has no boundary faces in the .ntiga database; repack with iga_pack");
			}
		}
		for (const auto& model : flow.OutletModels()) {
			long long local_faces = 0;
			for (const auto& element : flow.OwnedElements())
				for (const auto label : element.boundary_labels)
					if (label == model.label) ++local_faces;
			long long global_faces = 0;
			MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD);
			if (global_faces == 0) throw std::runtime_error("outlet model label "
				+std::to_string(model.label)+" has no boundary faces in the .ntiga database; repack with iga_pack");
		}
		if (configured) {
			std::map<int, double> tractions;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow traction catalog", [&] {
				const auto traction_configuration = iga::MaterializeOutletPressures(
					configuration, flow.OutletModels());
				tractions = iga::ExtractPressureTractions(traction_configuration,
					iga::FirstNavierStokesSystem(traction_configuration));
			});
			for (const auto& traction : tractions) {
				long long local_faces = 0;
				for (const auto& element : flow.OwnedElements())
					for (const auto label : element.boundary_labels)
						if (label == traction.first) ++local_faces;
				long long global_faces = 0;
				MPI_Allreduce(&local_faces, &global_faces, 1, MPI_LONG_LONG, MPI_SUM,
					PETSC_COMM_WORLD);
				if (global_faces == 0) throw std::runtime_error("pressure traction label "
					+std::to_string(traction.first)+" has no boundary faces in the .ntiga database; repack with iga_pack");
			}
		}
		double vca_reference_inlet_flow = 0.0;
		if (vca_circuit) {
			vca_reference_inlet_flow = flow.ReferenceBoundaryFlow(
				configuration.coupling.three_d_ports.inlet_label);
			if (!(std::abs(vca_reference_inlet_flow) > configuration.coupling.flow_epsilon_m3_s))
				throw std::runtime_error("VCA inlet initial_velocityfield.txt profile has zero integrated flow");
		}
		if (vca_has_transport)
			vca_transport = iga::AllocateCollectiveRuntime<iga::TransientTransportRuntime>(PETSC_COMM_WORLD, database,
				PETSC_COMM_WORLD, configuration, vca_transport_system, labels,
				std::map<std::uint64_t, iga::VolumeQuadratureRule>{}, std::set<std::string>{}, std::string{}, transport_solver_prefix);
		std::unique_ptr<iga::CouplingHistoryWriter> vca_history;
		fs::path vca_history_path;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow VCA identity and history", [&] {
			if (vca_transport && vca_transport->RequiredNodes() != flow.RequiredNodes())
				throw std::runtime_error("VCA flow and transport required-node layouts differ");
			if (!vca_circuit) return;
			vca_checkpoint_identity.configuration_fingerprint = iga::VcaConfigurationFingerprint(
				options.case_dir/"simulation_config.json");
			vca_checkpoint_identity.transport_system = vca_transport_system.name;
			vca_checkpoint_identity.inlet_label = configuration.coupling.three_d_ports.inlet_label;
			vca_checkpoint_identity.outlet_labels = configuration.coupling.three_d_ports.outlet_labels;
			vca_checkpoint_identity.device_model = iga::VcaDeviceModelIdentity(configuration.coupling);
			fs::path directory = options.output.empty() ? options.case_dir/"results"/"vca_flow"
				: options.output.parent_path();
			if (directory.empty()) directory = ".";
			vca_history_path = directory/"coupling_manifest.json";
			vca_history = std::make_unique<iga::CouplingHistoryWriter>(directory, configuration.coupling.mode);
		});

		int start_step = 0;
		std::vector<double> vca_species_state;
		std::map<std::string, double> vca_previous_mass;
		if (transient && !options.restart.empty()) {
			iga::FlowCheckpointMetadata metadata;
			std::string metadata_text;
			std::vector<iga::OutletModelState> restored_outlets;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow checkpoint metadata", [&] {
				metadata = iga::ReadFlowCheckpointMetadata(options.restart);
				iga::ValidateFlowCheckpoint(metadata, database.header().nodes, physical_steps,
					parameters.dt, parameters.density, parameters.dynamic_viscosity);
				metadata_text = iga::SerializeFlowCheckpointMetadata(metadata);
				restored_outlets = flow.OutletModels();
				iga::RestoreOutletCheckpoint(metadata, restored_outlets);
			});
			iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "flow checkpoint metadata agreement", metadata_text);
			ReadCheckpoint(flow.State(), options.restart, metadata);
			flow.OutletModels().swap(restored_outlets);
			start_step = metadata.completed_step;
			if (vca_circuit) {
				iga::VcaCheckpointMetadata vca_metadata;
				fs::path transport_path;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "VCA checkpoint metadata", [&] {
					if (!vca_transport) throw std::runtime_error("VCA checkpoint requires in-process transport state");
					vca_metadata = iga::ReadVcaCheckpointMetadata(options.restart);
					iga::ValidateVcaCheckpoint(vca_metadata, start_step, parameters.dt,
						vca_transport->System().fields, vca_checkpoint_identity);
					metadata_text = iga::SerializeVcaCheckpointMetadata(vca_metadata);
					transport_path = vca_metadata.transport_state_file;
					if (transport_path.is_relative())
						transport_path = iga::VcaCheckpointMetadataPath(options.restart).parent_path()/transport_path;
				});
				iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "VCA checkpoint metadata agreement", metadata_text);
				vca_transport->ReadState(transport_path);
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "VCA checkpoint reservoir", [&] {
					vca_circuit->RestoreState(vca_metadata.reservoir);
				});
				vca_species_state = vca_transport->GatherRequiredState();
				vca_previous_mass = vca_transport->TotalMass();
			}
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "checkpoint restart logging", [&] {
				if (rank == 0) std::cout << "restart=" << options.restart.string()
					<< " completed_step=" << start_step << " physical_time=" << metadata.physical_time << '\n';
				iga::FlushCheckedText(std::cout);
			});
		} else {
			if (vca_circuit) {
				const auto initial = PrepareFlowStepInput(PETSC_COMM_WORLD, "flow initial VCA input",
					true, true, configuration, options.case_dir, 0.0, 0.0, vca_circuit.get(),
					vca_transport ? &vca_transport->System() : nullptr, vca_reference_inlet_flow);
				flow.InitializeState(initial.configuration);
			} else {
				flow.InitializeState();
			}
		}
		flow.CopyStateToPrevious();
		if (!options.output.empty()
			&& (visualization_format == iga::VisualizationFormat::BezierVtkHdf || options.parallel_output)) {
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow visualization initialization", [&] {
				if (rank != 0) return;
				if (options.parallel_output && fs::exists(iga::PvdPath(options.output)))
					throw std::runtime_error("parallel output requires a new PVD path");
				bezier_mesh = std::make_unique<iga::BezierVisualizationMesh>(
					iga::BuildBezierVisualizationMesh(database, false));
				const auto report = iga::BezierGeometryReportPath(options.output);
				RequireRegularOutput(report);
				if (!options.parallel_output) RequireRegularOutput(iga::VtkHdfPath(options.output));
				iga::WriteBezierGeometryReport(report, bezier_mesh->validation);
				iga::RequireValidBezierGeometry(bezier_mesh->validation);
				if (options.parallel_output) {
					std::cout << "parallel_bezier_geometry_points=" << bezier_mesh->points.size()
						<< " geometry_report=" << report.string() << '\n';
					iga::FlushCheckedText(std::cout);bezier_mesh.reset();return;
				}
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

		const auto start = std::chrono::steady_clock::now();
		std::vector<std::pair<double, fs::path>> vtk_snapshots;
		std::vector<iga::VelocitySnapshot> velocity_snapshots;
		int last_parallel_output_step = -1;
		auto write_output = [&](int step, double physical_time, bool final_output) {
			fs::path text_path, vtk_path, mesh_path;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow output preparation", [&] {
				text_path = final_output ? options.output : iga::TimeIndexedPath(options.output, step);
				vtk_path = final_output ? iga::VtuFinalPath(options.output) : iga::VtuStepPath(options.output, step);
				mesh_path = options.case_dir/"controlmesh.vtk";
			});
			if (options.parallel_output) {
				if (last_parallel_output_step == step) return;
				iga::PhaseScope output_phase(iga::ProfilePhase::Output);
				fs::path directory;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"parallel flow output path",[&] {
					directory=iga::VtuStepPath(options.output,step);directory.replace_extension();
				});
				const auto extracted=iga::BuildPetscBezierPartition(flow.State(),flow.OwnedElements(),
					database.header().elements,database.header().nodes,{{"velocity",3},{"pressure",1}});
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,directory,extracted.piece,physical_time);
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"parallel flow output bookkeeping",[&] {
					vtk_snapshots.push_back({physical_time,directory/"snapshot.pvtu"});
					std::cout<<"parallel_flow_output rank="<<rank<<" step="<<step<<" selected_rows="<<extracted.requested_rows
						<<" global_rows="<<extracted.global_rows<<'\n';iga::FlushCheckedText(std::cout);
				});
				iga::WriteParallelVtkSeries(PETSC_COMM_WORLD,iga::PvdPath(options.output),vtk_snapshots);
				last_parallel_output_step=step;return;
			}
			WriteFlowOutput(flow.State(), database.header().nodes, text_path,
				mesh_path, vtk_path, physical_time, rank, visualization_format, vtkhdf.get());
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow output bookkeeping", [&] {
				if (!final_output) {
					if (visualization_format == iga::VisualizationFormat::Vtu)
						vtk_snapshots.push_back({physical_time, vtk_path});
					velocity_snapshots.push_back({physical_time, text_path});
				}
			});
		};
		if (transient && options.output_every > 0)
			write_output(start_step, start_step*parameters.dt, false);
		for (int step = start_step; step < run_end_step; ++step) {
			const auto physical_time = transient ? (step+1)*parameters.dt : 0.0;
			const auto step_input = PrepareFlowStepInput(PETSC_COMM_WORLD, "flow step input",
				configured, transient, configuration, options.case_dir, physical_time, step*parameters.dt,
				vca_circuit.get(), vca_transport ? &vca_transport->System() : nullptr, vca_reference_inlet_flow);
			const auto& step_configuration = step_input.configuration;
			const auto& inlet = step_input.inlet;
			flow.BeginStep(step, physical_time, options.max_newton,
				options.nonlinear_relative_tolerance, options.nonlinear_absolute_tolerance,
				options.mass_relative_tolerance);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow step boundaries", [&] {
				if (configured) flow.SetTrialBoundaryConfiguration(step_configuration);
			});
			flow.SolveTrial();
			flow.CommitStep();
			if (vca_transport) {
				const auto velocity = flow.GatherRequiredVelocity();
				vca_transport->Advance(step_configuration, flow.RequiredNodes(), velocity);
				vca_species_state = vca_transport->GatherRequiredState();
			}
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow solver diagnostics", [&] {
				if (rank != 0) return;
				const auto write = [&](const iga::PetscKspConfiguration& solver) {
					std::cout << "solver_configuration prefix=" << solver.prefix << " ksp=" << solver.ksp
						<< " pc=" << solver.pc << " factor_backend=" << solver.factor_backend
						<< " step=" << step << " iterations=" << solver.last_iterations
						<< " reason=" << static_cast<int>(solver.last_reason) << '\n';
				};
				write(flow.SolverConfiguration());
				if (vca_transport) write(vca_transport->SolverConfiguration());
				iga::FlushCheckedText(std::cout);
			});
			if (vca_circuit) {
				const std::vector<std::string> empty_fields;
				const auto& species_fields = vca_transport ? vca_transport->System().fields : empty_fields;
				const auto ports = flow.MeasurePorts(configuration.coupling.three_d_ports,
					species_fields, vca_species_state, vca_transport ? &vca_transport->System() : nullptr);
				auto result = PrepareFlowVcaResult(PETSC_COMM_WORLD, physical_time, parameters.dt, inlet,
					configuration.coupling.three_d_ports, ports, configuration.coupling.flow_epsilon_m3_s);
				if (vca_transport) {
					auto total_mass = vca_transport->TotalMass();
					auto source_integrals = vca_transport->SourceIntegrals();
					iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow VCA transport budget", [&] {
						result.total_mass = std::move(total_mass);
						result.source_integrals = std::move(source_integrals);
						for (const auto& mass : result.total_mass) {
							const auto old = vca_previous_mass.find(mass.first);
							if (old == vca_previous_mass.end()) continue;
							double boundary_flux = 0.0;
							const auto inlet_flux = ports.species_fluxes.find(
								configuration.coupling.three_d_ports.inlet_label);
							if (inlet_flux != ports.species_fluxes.end()) {
								const auto species = inlet_flux->second.find(mass.first);
								if (species != inlet_flux->second.end()) boundary_flux += species->second;
							}
							for (const auto& outlet : result.outlets) {
								const auto flux = outlet.species_flux.find(mass.first);
								if (flux != outlet.species_flux.end()) boundary_flux += flux->second;
							}
							const auto source = result.source_integrals.find(mass.first);
							result.balance_residuals[mass.first]
								= (mass.second-old->second)/parameters.dt+boundary_flux
								-(source == result.source_integrals.end() ? 0.0 : source->second);
						}
						vca_previous_mass = result.total_mass;
					});
				}
				AdvanceFlowVcaCircuit(PETSC_COMM_WORLD, *vca_circuit, *vca_history,
					result, configuration.coupling.flow_epsilon_m3_s);
			}
			const auto completed_step = step+1;
			if (options.output_every > 0 && completed_step%options.output_every == 0)
				write_output(completed_step, physical_time, false);
			if (!options.checkpoint.empty() && (completed_step == run_end_step
				|| (options.checkpoint_every > 0 && completed_step%options.checkpoint_every == 0))) {
				iga::FlowCheckpointMetadata metadata;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow checkpoint metadata preparation", [&] {
					metadata.nodes = database.header().nodes;
					metadata.completed_step = completed_step;
					metadata.physical_time = completed_step*parameters.dt;
					metadata.dt = parameters.dt;
					metadata.density = parameters.density;
					metadata.viscosity = parameters.dynamic_viscosity;
					metadata.state_file = iga::FlowCheckpointStatePath(options.checkpoint).filename().string();
					metadata.state_format = "petsc_binary";
					iga::AppendOutletCheckpoint(flow.OutletModels(), metadata);
				});
				WriteCheckpoint(flow.State(), options.checkpoint, metadata, rank);
				if (vca_circuit) {
					iga::VcaCheckpointMetadata vca_metadata;
					fs::path transport_path;
					std::string metadata_text;
					iga::CollectiveLocalStage(PETSC_COMM_WORLD, "VCA checkpoint write preparation", [&] {
						if (!vca_transport) throw std::runtime_error("VCA checkpoint requires in-process transport state");
						transport_path = iga::VcaCheckpointTransportStatePath(options.checkpoint);
						vca_metadata.completed_step = completed_step;
						vca_metadata.physical_time = completed_step*parameters.dt;
						vca_metadata.dt = parameters.dt;
						vca_metadata.fields = vca_transport->System().fields;
						vca_metadata.identity = vca_checkpoint_identity;
						vca_metadata.transport_state_file = transport_path.filename().string();
						vca_metadata.reservoir = vca_circuit->State();
						metadata_text = iga::SerializeVcaCheckpointMetadata(vca_metadata);
					});
					iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "VCA checkpoint write agreement", metadata_text);
					vca_transport->WriteState(transport_path);
					iga::CollectiveLocalStage(PETSC_COMM_WORLD, "VCA checkpoint metadata write", [&] {
						if (rank == 0) WriteCheckpointMetadataText(iga::VcaCheckpointMetadataPath(options.checkpoint), metadata_text);
					});
				}
				iga::CollectiveLocalStage(PETSC_COMM_WORLD, "checkpoint write logging", [&] {
					if (rank == 0) std::cout << "checkpoint=" << options.checkpoint.string()
						<< " completed_step=" << completed_step << '\n';
					iga::FlushCheckedText(std::cout);
				});
			}
		}
		if (vca_circuit) {
			iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow VCA history output", [&] {
				if (rank != 0) return;
				RequireRegularOutput(vca_history_path);
				vca_history->Write("cpu_3d_navier_stokes");
			});
		}
		const auto summary = flow.Summary();
		// Preserve the existing timing interval while publishing success only
		// after the final field and index writers have returned successfully.
		const auto solve_seconds = std::chrono::duration<double>(
			std::chrono::steady_clock::now()-start).count();
		if (!options.output.empty()) {
			const auto final_time = transient ? run_end_step*parameters.dt : 0.0;
			write_output(run_end_step, final_time, true);
			if (!options.parallel_output) iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow output index", [&] {
				if (rank != 0) return;
				if (visualization_format == iga::VisualizationFormat::Vtu) {
					if (vtk_snapshots.empty())
						vtk_snapshots.push_back({final_time, iga::VtuFinalPath(options.output)});
					RequireRegularOutput(iga::PvdPath(options.output));
					iga::WritePvd(iga::PvdPath(options.output), vtk_snapshots);
				}
				if (!velocity_snapshots.empty()) {
					RequireRegularOutput(iga::VelocityManifestPath(options.output));
					iga::WriteVelocityManifest(iga::VelocityManifestPath(options.output), velocity_snapshots);
				}
			});
		}
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow visualization close", [&] {
			if (rank == 0 && vtkhdf) vtkhdf->Close();
		});
		if (vca_transport) vca_transport->Close();
		flow.Close();
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow completion logging", [&] {
			if (rank == 0) std::cout << "navier_stokes_v2 seconds=" << solve_seconds
				<< " total_linear_iterations=" << summary.linear_iterations
				<< " state_l2=" << summary.state_l2 << " velocity_l2=" << summary.velocity_l2
				<< " pressure_l2=" << summary.pressure_l2 << '\n';
			iga::FlushCheckedText(std::cout);
		});
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	try {
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "flow profile logging", [&] {
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

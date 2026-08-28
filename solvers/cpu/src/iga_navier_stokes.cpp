#include "BoundarySupport.hpp"
#include "CaseInput.hpp"
#include "CouplingHistory.hpp"
#include "FlowCheckpoint.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "OutletCheckpoint.hpp"
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
#include <map>
#include <memory>
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
	iga::VisualizationFormat visualization_format = iga::VisualizationFormat::Automatic;
	double nonlinear_relative_tolerance = 1e-5;
	double nonlinear_absolute_tolerance = 1e-10;
	double mass_relative_tolerance = 1e-3;
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
		"[--visualization-format auto|vtu|vtkhdf]");
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
		else if (argument == "--visualization-format")
			options.visualization_format = iga::ParseVisualizationFormat(value);
		else throw std::runtime_error("unknown option: "+argument);
		PetscOptionsClearValue(nullptr, argument.c_str());
	}
	if (options.output_every > 0 && options.output.empty())
		throw std::runtime_error("--output-every requires --output or legacy OUTPUT");
	if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

void WriteFlowOutput(Vec state, std::uint64_t nodes, const fs::path& path,
	const fs::path& mesh_path, const fs::path& vtk_path, double physical_time, int rank,
	iga::VisualizationFormat visualization_format,
	iga::TemporalVtkHdfWriter* vtkhdf)
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
			std::vector<double> velocity(3*static_cast<std::size_t>(nodes));
			std::vector<double> pressure(static_cast<std::size_t>(nodes));
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
			VecRestoreArrayRead(root, &values);
			if (!output || !pressure_output) throw std::runtime_error("cannot write Navier-Stokes output");
			std::vector<iga::VtkPointArray> arrays{
				{"velocity", 3, std::move(velocity)},
				{"pressure", 1, std::move(pressure)}};
			if (visualization_format == iga::VisualizationFormat::Vtu)
				iga::WriteVtu(mesh_path, vtk_path, arrays, physical_time);
			else {
				if (!vtkhdf) throw std::runtime_error("VTKHDF writer is unavailable");
				vtkhdf->Append(physical_time, arrays);
			}
		} catch (const std::exception& error) {
			std::cerr << "rank 0: " << error.what() << '\n';
			write_failed = 1;
		}
	}
	MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	VecScatterDestroy(&scatter);
	VecDestroy(&root);
	if (write_failed) throw std::runtime_error("cannot write Navier-Stokes output: "+path.string());
}

void WriteCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata, int rank)
{
	PetscViewer viewer = nullptr;
	const auto state_path = iga::FlowCheckpointStatePath(prefix);
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
		"cannot write flow checkpoint metadata: "+iga::FlowCheckpointMetadataPath(prefix).string());
}

void ReadCheckpoint(Vec state, const fs::path& prefix,
	const iga::FlowCheckpointMetadata& metadata)
{
	if (metadata.state_format != "petsc_binary")
		throw std::runtime_error("CPU flow restart requires petsc_binary checkpoint state");
	fs::path path(metadata.state_file);
	if (path.is_relative()) path = iga::FlowCheckpointMetadataPath(prefix).parent_path()/path;
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, path.string().c_str(), FILE_MODE_READ, &viewer);
	VecLoad(state, viewer);
	PetscViewerDestroy(&viewer);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr,
		"TubularFlowIGA stabilized steady/transient Navier-Stokes solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto options = ParseOptions(argc, argv);
		iga::Database database(options.database.string());
		std::unique_ptr<iga::BezierVisualizationMesh> bezier_mesh;
		std::unique_ptr<iga::TemporalVtkHdfWriter> vtkhdf;
		const auto mesh = iga::ReadLabeledHexMesh((options.case_dir/"controlmesh.vtk").string(),
			database.header().nodes, database.header().elements);
		const auto& labels = mesh.labels;
		const auto wall_trace_basis = iga::WallTraceBasis(database, mesh);
		const auto boundary_velocity = iga::ReadVelocity(
			(options.case_dir/"initial_velocityfield.txt").string(), database.header().nodes);
		iga::SimulationConfiguration configuration;
		iga::ResolvedBoundaryConditions boundaries;
		std::vector<iga::OutletModelState> outlet_models;
		iga::NavierStokesParameters parameters{1.0, 0.1, 0.0};
		bool configured = false;
		bool transient = false;
		std::string boundary_config;
		std::unique_ptr<iga::VcaExternalCircuit> vca_circuit;
		std::unique_ptr<iga::TransientTransportRuntime> vca_transport;
		iga::CompiledLinearSystem vca_transport_system;
		iga::VcaCheckpointIdentity vca_checkpoint_identity;
		bool vca_has_transport = false;
		if (fs::exists(options.case_dir/"simulation_config.json")) {
			configuration = iga::ReadSimulationConfiguration(
				(options.case_dir/"simulation_config.json").string());
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
		const auto physical_steps = transient ? configuration.time.steps : 1;
		if (options.stop_after_step > physical_steps)
			throw std::runtime_error("--stop-after-step exceeds configured physical steps");
		const auto run_end_step = options.stop_after_step > 0 ? options.stop_after_step : physical_steps;
		const auto visualization_format = iga::ResolveVisualizationFormat(
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

		iga::TransientFlowRuntime flow(database, PETSC_COMM_WORLD, configured, transient,
			parameters, boundaries, labels, boundary_velocity, wall_trace_basis,
			std::move(outlet_models));
		iga::RequireValidGeometry(flow.Elements(), rank, PETSC_COMM_WORLD);
		if (vca_circuit) {
			std::map<int, long long> port_faces;
			port_faces.emplace(configuration.coupling.three_d_ports.inlet_label, 0);
			for (const auto label : configuration.coupling.three_d_ports.outlet_labels)
				port_faces.emplace(label, 0);
			for (const auto& element : flow.OwnedElements())
				for (const auto label : element.boundary_labels) {
					auto found = port_faces.find(label);
					if (found != port_faces.end()) ++found->second;
				}
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
			const auto traction_configuration = iga::MaterializeOutletPressures(
				configuration, flow.OutletModels());
			const auto tractions = iga::ExtractPressureTractions(traction_configuration,
				iga::FirstNavierStokesSystem(traction_configuration));
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
			vca_transport = std::make_unique<iga::TransientTransportRuntime>(database,
				PETSC_COMM_WORLD, configuration, vca_transport_system, labels);
		if (vca_transport && vca_transport->RequiredNodes() != flow.RequiredNodes())
			throw std::runtime_error("VCA flow and transport required-node layouts differ");
		if (vca_circuit) {
			vca_checkpoint_identity.configuration_fingerprint = iga::VcaConfigurationFingerprint(
				options.case_dir/"simulation_config.json");
			vca_checkpoint_identity.transport_system = vca_transport_system.name;
			vca_checkpoint_identity.inlet_label = configuration.coupling.three_d_ports.inlet_label;
			vca_checkpoint_identity.outlet_labels = configuration.coupling.three_d_ports.outlet_labels;
			vca_checkpoint_identity.device_model = iga::VcaDeviceModelIdentity(configuration.coupling);
		}
		std::unique_ptr<iga::CouplingHistoryWriter> vca_history;
		if (vca_circuit) {
			fs::path directory = options.output.empty() ? options.case_dir/"results"/"vca_flow"
				: options.output.parent_path();
			if (directory.empty()) directory = ".";
			vca_history = std::make_unique<iga::CouplingHistoryWriter>(directory, configuration.coupling.mode);
		}

		int start_step = 0;
		std::vector<double> vca_species_state;
		std::map<std::string, double> vca_previous_mass;
		if (transient && !options.restart.empty()) {
			const auto metadata = iga::ReadFlowCheckpointMetadata(options.restart);
			iga::ValidateFlowCheckpoint(metadata, database.header().nodes, physical_steps,
				parameters.dt, parameters.density, parameters.dynamic_viscosity);
			ReadCheckpoint(flow.State(), options.restart, metadata);
			iga::RestoreOutletCheckpoint(metadata, flow.OutletModels());
			start_step = metadata.completed_step;
			if (vca_circuit) {
				if (!vca_transport) throw std::runtime_error("VCA checkpoint requires in-process transport state");
				const auto vca_metadata = iga::ReadVcaCheckpointMetadata(options.restart);
				iga::ValidateVcaCheckpoint(vca_metadata, start_step, parameters.dt,
					vca_transport->System().fields, vca_checkpoint_identity);
				fs::path transport_path(vca_metadata.transport_state_file);
				if (transport_path.is_relative())
					transport_path = iga::VcaCheckpointMetadataPath(options.restart).parent_path()/transport_path;
				vca_transport->ReadState(transport_path);
				vca_circuit->RestoreState(vca_metadata.reservoir);
				vca_species_state = vca_transport->GatherRequiredState();
				vca_previous_mass = vca_transport->TotalMass();
			}
			if (rank == 0) std::cout << "restart=" << options.restart.string()
				<< " completed_step=" << start_step << " physical_time=" << metadata.physical_time << '\n';
		} else {
			if (vca_circuit) {
				auto initial_configuration = iga::MaterializeBoundaryWaveforms(configuration,
					options.case_dir.string(), 0.0);
				const auto initial_inlet = vca_circuit->InletState(0.0);
				iga::ApplyThreeDVascularInlet(initial_configuration,
					iga::FirstNavierStokesSystem(initial_configuration), initial_inlet,
					vca_reference_inlet_flow);
				if (vca_transport)
					iga::ApplyThreeDVascularSpeciesInlet(initial_configuration,
						vca_transport->System(), initial_inlet);
				flow.InitializeState(initial_configuration);
			} else {
				flow.InitializeState();
			}
		}
		flow.CopyStateToPrevious();
		if (!options.output.empty()
			&& visualization_format == iga::VisualizationFormat::BezierVtkHdf) {
			int initialization_failed = 0;
			if (rank == 0) {
				try {
					bezier_mesh = std::make_unique<iga::BezierVisualizationMesh>(
						iga::BuildBezierVisualizationMesh(database, false));
					const auto report = iga::BezierGeometryReportPath(options.output);
					iga::WriteBezierGeometryReport(report, bezier_mesh->validation);
					iga::RequireValidBezierGeometry(bezier_mesh->validation);
					vtkhdf = std::make_unique<iga::TemporalVtkHdfWriter>(
						iga::VtkHdfPath(options.output), *bezier_mesh,
						!options.restart.empty());
					std::cout << "bezier_geometry_points=" << bezier_mesh->points.size()
						<< " local_point_references=" << bezier_mesh->validation.local_points
						<< " geometry_report=" << report.string()
						<< " vtkhdf=" << vtkhdf->path().string() << '\n';
				} catch (const std::exception& error) {
					std::cerr << "rank 0: " << error.what() << '\n';
					initialization_failed = 1;
				}
			}
			MPI_Bcast(&initialization_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			if (initialization_failed)
				throw std::runtime_error("cannot initialize Bezier VTKHDF output");
		}

		const auto start = std::chrono::steady_clock::now();
		std::vector<std::pair<double, fs::path>> vtk_snapshots;
		std::vector<iga::VelocitySnapshot> velocity_snapshots;
		if (transient && options.output_every > 0) {
			const auto text_path = iga::TimeIndexedPath(options.output, start_step);
			const auto vtk_path = iga::VtuStepPath(options.output, start_step);
			WriteFlowOutput(flow.State(), database.header().nodes, text_path,
				options.case_dir/"controlmesh.vtk", vtk_path, start_step*parameters.dt, rank,
				visualization_format, vtkhdf.get());
			if (visualization_format == iga::VisualizationFormat::Vtu)
				vtk_snapshots.push_back({start_step*parameters.dt, vtk_path});
			velocity_snapshots.push_back({start_step*parameters.dt, text_path});
		}
		for (int step = start_step; step < run_end_step; ++step) {
			const auto physical_time = transient ? (step+1)*parameters.dt : 0.0;
			iga::SimulationConfiguration step_configuration;
			if (configured)
				step_configuration = transient
					? iga::MaterializeBoundaryWaveforms(configuration, options.case_dir.string(), physical_time)
					: configuration;
			iga::VascularInletState inlet;
			if (vca_circuit) {
				inlet = vca_circuit->InletState(step*parameters.dt);
				iga::ApplyThreeDVascularInlet(step_configuration,
					iga::FirstNavierStokesSystem(step_configuration), inlet, vca_reference_inlet_flow);
				if (vca_transport)
					iga::ApplyThreeDVascularSpeciesInlet(step_configuration,
						vca_transport->System(), inlet);
			}
			flow.Advance(step_configuration, step, physical_time, options.max_newton,
				options.nonlinear_relative_tolerance, options.nonlinear_absolute_tolerance,
				options.mass_relative_tolerance);
			if (vca_transport) {
				vca_transport->Advance(step_configuration, flow.RequiredNodes(),
					flow.GatherRequiredVelocity());
				vca_species_state = vca_transport->GatherRequiredState();
			}
			if (vca_circuit) {
				const auto ports = flow.MeasurePorts(configuration.coupling.three_d_ports,
					vca_transport ? vca_transport->System().fields : std::vector<std::string>{},
					vca_species_state);
				auto result = iga::BuildThreeDFlowPortResult(physical_time, parameters.dt, inlet,
					configuration.coupling.three_d_ports, ports.flows, ports.pressures);
				for (auto& outlet : result.outlets) {
					const auto flux = ports.species_fluxes.find(outlet.outlet_id);
					if (flux == ports.species_fluxes.end()) continue;
					outlet.species_flux = flux->second;
					outlet.average_valid = std::abs(outlet.flow_m3_s)
						> configuration.coupling.flow_epsilon_m3_s;
					if (outlet.average_valid)
						for (const auto& species : outlet.species_flux)
							outlet.flux_weighted_concentration[species.first]
								= species.second/outlet.flow_m3_s;
				}
				if (vca_transport) {
					result.total_mass = vca_transport->TotalMass();
					result.source_integrals = vca_transport->SourceIntegrals();
					for (const auto& mass : result.total_mass) {
						const auto old = vca_previous_mass.find(mass.first);
						if (old == vca_previous_mass.end()) continue;
						const auto inlet_flux = inlet.species.count(mass.first)
							? inlet.flow_m3_s*inlet.species.at(mass.first) : 0.0;
						double outlet_flux = 0.0;
						for (const auto& outlet : result.outlets) {
							const auto flux = outlet.species_flux.find(mass.first);
							if (flux != outlet.species_flux.end()) outlet_flux += flux->second;
						}
						const auto source = result.source_integrals.find(mass.first);
						result.balance_residuals[mass.first]
							= (mass.second-old->second)/parameters.dt-inlet_flux+outlet_flux
							-(source == result.source_integrals.end() ? 0.0 : source->second);
					}
					vca_previous_mass = result.total_mass;
				}
				const auto venous = iga::AggregateVascularOutlets(result.outlets,
					configuration.coupling.flow_epsilon_m3_s);
				const auto report = vca_circuit->Advance(venous, parameters.dt, physical_time);
				vca_history->Add(result, venous, report);
			}
			const auto completed_step = step+1;
			if (options.output_every > 0 && completed_step%options.output_every == 0) {
				const auto text_path = iga::TimeIndexedPath(options.output, completed_step);
				const auto vtk_path = iga::VtuStepPath(options.output, completed_step);
				WriteFlowOutput(flow.State(), database.header().nodes, text_path,
					options.case_dir/"controlmesh.vtk", vtk_path, physical_time, rank,
					visualization_format, vtkhdf.get());
				if (visualization_format == iga::VisualizationFormat::Vtu)
					vtk_snapshots.push_back({physical_time, vtk_path});
				velocity_snapshots.push_back({physical_time, text_path});
			}
			if (!options.checkpoint.empty() && (completed_step == run_end_step
				|| (options.checkpoint_every > 0 && completed_step%options.checkpoint_every == 0))) {
				iga::FlowCheckpointMetadata metadata;
				metadata.nodes = database.header().nodes;
				metadata.completed_step = completed_step;
				metadata.physical_time = completed_step*parameters.dt;
				metadata.dt = parameters.dt;
				metadata.density = parameters.density;
				metadata.viscosity = parameters.dynamic_viscosity;
				metadata.state_file = iga::FlowCheckpointStatePath(options.checkpoint).filename().string();
				metadata.state_format = "petsc_binary";
				iga::AppendOutletCheckpoint(flow.OutletModels(), metadata);
				WriteCheckpoint(flow.State(), options.checkpoint, metadata, rank);
				if (vca_circuit) {
					if (!vca_transport)
						throw std::runtime_error("VCA checkpoint requires in-process transport state");
					vca_transport->WriteState(iga::VcaCheckpointTransportStatePath(options.checkpoint));
					int write_failed = 0;
					if (rank == 0) {
						try {
							iga::VcaCheckpointMetadata vca_metadata;
							vca_metadata.completed_step = completed_step;
							vca_metadata.physical_time = completed_step*parameters.dt;
							vca_metadata.dt = parameters.dt;
							vca_metadata.fields = vca_transport->System().fields;
							vca_metadata.identity = vca_checkpoint_identity;
							vca_metadata.transport_state_file
								= iga::VcaCheckpointTransportStatePath(options.checkpoint).filename().string();
							vca_metadata.reservoir = vca_circuit->State();
							iga::WriteVcaCheckpointMetadata(options.checkpoint, vca_metadata);
						} catch (const std::exception&) {
							write_failed = 1;
						}
					}
					MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
					if (write_failed) throw std::runtime_error("cannot write VCA checkpoint metadata");
				}
				if (rank == 0) std::cout << "checkpoint=" << options.checkpoint.string()
					<< " completed_step=" << completed_step << '\n';
			}
		}
		if (vca_circuit) {
			int write_failed = 0;
			if (rank == 0) {
				try {
					vca_history->Write("cpu_3d_navier_stokes");
				} catch (const std::exception&) {
					write_failed = 1;
				}
			}
			MPI_Bcast(&write_failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
			if (write_failed) throw std::runtime_error("cannot write VCA coupling manifest");
		}
		const auto summary = flow.Summary();
		if (rank == 0) std::cout << "navier_stokes_v2 seconds="
			<< std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
			<< " total_linear_iterations=" << summary.linear_iterations
			<< " state_l2=" << summary.state_l2 << " velocity_l2=" << summary.velocity_l2
			<< " pressure_l2=" << summary.pressure_l2 << '\n';
		if (!options.output.empty()) {
			const auto final_time = transient ? run_end_step*parameters.dt : 0.0;
			WriteFlowOutput(flow.State(), database.header().nodes, options.output,
				options.case_dir/"controlmesh.vtk", iga::VtuFinalPath(options.output), final_time, rank,
				visualization_format, vtkhdf.get());
			if (rank == 0) {
				if (visualization_format == iga::VisualizationFormat::Vtu) {
					if (vtk_snapshots.empty())
						vtk_snapshots.push_back({final_time, iga::VtuFinalPath(options.output)});
					iga::WritePvd(iga::PvdPath(options.output), vtk_snapshots);
				}
				if (!velocity_snapshots.empty())
					iga::WriteVelocityManifest(iga::VelocityManifestPath(options.output), velocity_snapshots);
			}
		}
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

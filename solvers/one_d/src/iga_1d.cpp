#include "OneDCheckpoint.hpp"
#include "OneDCoupling.hpp"
#include "OneDImplicit.hpp"
#include "OneDOutput.hpp"
#include "CouplingReplay.hpp"

#include <petscsys.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Options {
	fs::path case_directory;
	std::string system;
	fs::path output_directory;
	bool check = false;
	fs::path checkpoint;
	int checkpoint_every = 0;
	fs::path restart;
	int stop_after_step = 0;
};

int PositiveInteger(const std::string& text, const std::string& option)
{
	std::size_t used = 0;
	int value = 0;
	try { value = std::stoi(text, &used); }
	catch (const std::exception&) { throw std::runtime_error(option+" requires a positive integer"); }
	if (used != text.size() || value < 1) throw std::runtime_error(option+" requires a positive integer");
	return value;
}

Options ParseOptions(int argc, char** argv)
{
	if (argc < 2) throw std::runtime_error(
		"usage: iga_1d CASE_DIR [--system NAME] [--output-dir DIR] [--check] "
		"[--checkpoint PREFIX --checkpoint-every N] [--restart PREFIX] "
		"[--stop-after-step N] [PETSc options]");
	Options options;
	options.case_directory = argv[1];
	for (int i = 2; i < argc; ++i) {
		const std::string argument(argv[i]);
		if (argument == "--check") { options.check = true; continue; }
		if (argument == "--system" || argument == "--output-dir"
			|| argument == "--checkpoint" || argument == "--checkpoint-every"
			|| argument == "--restart" || argument == "--stop-after-step") {
			if (++i >= argc) throw std::runtime_error(argument+" requires a value");
			const std::string value(argv[i]);
			if (argument == "--system") options.system = value;
			else if (argument == "--output-dir") options.output_directory = value;
			else if (argument == "--checkpoint") options.checkpoint = value;
			else if (argument == "--checkpoint-every") options.checkpoint_every = PositiveInteger(value, argument);
			else if (argument == "--restart") options.restart = value;
			else options.stop_after_step = PositiveInteger(value, argument);
			continue;
		}
		if (!argument.empty() && argument[0] == '-') continue;
		throw std::runtime_error("unexpected argument: "+argument);
	}
		if (options.checkpoint_every > 0 && options.checkpoint.empty())
		throw std::runtime_error("--checkpoint-every requires --checkpoint");
	return options;
}

std::string ReadText(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open file: "+path.string());
	std::ostringstream text;
	text << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error("cannot read file: "+path.string());
	return text.str();
}

const iga::OneDFlowSystemDefinition& SelectFlow(const iga::OneDConfiguration& configuration,
	const std::string& name)
{
	if (name.empty()) {
		if (configuration.flow_systems.size() != 1)
			throw std::runtime_error("--system is required when a 1d case contains multiple flow systems");
		return configuration.flow_systems.front();
	}
	for (const auto& flow : configuration.flow_systems)
		if (flow.name == name) return flow;
	throw std::runtime_error("unknown 1d flow system '"+name+"'");
}

std::vector<iga::OneDTransportState> InitializeTransports(
	const iga::OneDConfiguration& configuration,
	const iga::OneDFlowSystemDefinition& flow, const iga::OneDNetwork& network)
{
	std::vector<iga::OneDTransportState> result;
	for (const auto& transport : configuration.transport_systems)
		if (transport.flow_system == flow.name)
			result.push_back(iga::InitializeOneDTransport(configuration, transport, network));
	return result;
}

iga::OneDCheckpointMetadata CheckpointMetadata(
	const iga::OneDConfiguration& configuration, const iga::OneDNetwork& network,
	const iga::OneDFlowState& state, const std::vector<iga::OneDTransportState>& transports,
	std::uint64_t config_fingerprint, const fs::path& prefix)
{
	iga::OneDCheckpointMetadata metadata;
	metadata.completed_step = state.completed_step;
	metadata.internal_substeps = state.internal_substeps;
	metadata.physical_time = state.physical_time;
	metadata.dt = configuration.time.dt;
	metadata.inlet_flow = state.inlet_flow;
	metadata.cells = network.cells;
	metadata.nodes = static_cast<int>(network.nodes.size());
	metadata.segments = static_cast<int>(network.segments.size());
	metadata.outlets = static_cast<int>(state.outlets.size());
	metadata.species = iga::OneDCheckpointSpecies(transports);
	metadata.config_fingerprint = config_fingerprint;
	metadata.network_fingerprint = iga::OneDNetworkFingerprint(network);
	metadata.state_file = iga::OneDCheckpointStatePath(prefix).filename().string();
	return metadata;
}

void ValidateRestart(const iga::OneDCheckpointMetadata& metadata,
	const iga::OneDConfiguration& configuration, const iga::OneDNetwork& network,
	const std::vector<iga::OneDTransportState>& transports, std::uint64_t config_fingerprint)
{
	auto close = [](double first, double second) {
		return std::abs(first-second) <= 1.0e-12*std::max({1.0, std::abs(first), std::abs(second)});
	};
	if (metadata.cells != network.cells || metadata.nodes != static_cast<int>(network.nodes.size())
		|| metadata.segments != static_cast<int>(network.segments.size())
		|| metadata.network_fingerprint != iga::OneDNetworkFingerprint(network))
		throw std::runtime_error("1d checkpoint network does not match the case");
	if (metadata.config_fingerprint != config_fingerprint || !close(metadata.dt, configuration.time.dt))
		throw std::runtime_error("1d checkpoint configuration does not match the case");
	if (metadata.completed_step > configuration.time.steps
		|| !close(metadata.physical_time, metadata.completed_step*configuration.time.dt))
		throw std::runtime_error("1d checkpoint time is inconsistent with the configured run");
	if (metadata.species != iga::OneDCheckpointSpecies(transports))
		throw std::runtime_error("1d checkpoint species do not match the case");
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA native one-dimensional solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		const auto setup_start = std::chrono::steady_clock::now();
		auto options = ParseOptions(argc, argv);
		const auto config_path = options.case_directory/"simulation_config.json";
		const auto config_text = ReadText(config_path);
		auto configuration = iga::ParseOneDConfiguration(config_text);
		const auto& flow = SelectFlow(configuration, options.system);
		auto network = iga::ReadOneDNetwork(options.case_directory/configuration.geometry.file,
			configuration.geometry.length_scale_to_m, flow.discretization.cells_per_segment,
			flow.dynamic_viscosity, configuration.geometry.root_node_id);
		iga::ValidateOneDTopologyReferences(configuration, network);
		const auto inlet = iga::ResolveOneDInlet(configuration);
		auto flow_state = iga::OneDFlowState{};
		flow_state.outlets = iga::ResolveOneDOutlets(configuration, network);
		auto transports = InitializeTransports(configuration, flow, network);
		std::unique_ptr<iga::ReplayInletProvider> replay;
		std::unique_ptr<iga::VcaExternalCircuit> circuit;
		if (configuration.coupling.mode == iga::SimulationScopeMode::VcaReplay)
			replay = std::make_unique<iga::ReplayInletProvider>(
				iga::ReplayInletProvider::Read(options.case_directory
					/configuration.coupling.replay_file));
		if (configuration.coupling.mode == iga::SimulationScopeMode::VcaClosedLoop)
			circuit = std::make_unique<iga::VcaExternalCircuit>(configuration.coupling);
		if (circuit && (!options.restart.empty() || !options.checkpoint.empty()))
			throw std::runtime_error(
				"closed-loop checkpoint/restart requires coupled reservoir state and is not yet enabled");
		const auto config_fingerprint = iga::OneDFingerprint(config_text);
		if (options.check) {
			if (rank == 0) std::cout << "schema_version=3 dimension=1d system=" << flow.name
				<< " nodes=" << network.nodes.size() << " segments=" << network.segments.size()
				<< " root_id=" << network.nodes[static_cast<std::size_t>(network.root)].id
				<< " cells=" << network.cells << " outlets=" << network.outlet_nodes.size()
				<< " transport_systems=" << transports.size() << '\n';
			PetscFinalize();
			return 0;
		}
		const double inlet_area = network.segments.front().area0;
		auto open_loop_inlet = [&](double time, double flow_value) {
			iga::VascularInletState state;
			state.time_s = time;
			state.has_flow = true;
			state.flow_m3_s = flow_value;
			for (const auto& transport : transports)
				for (const auto& species : transport.species)
					state.species[species.definition.field]
						= iga::EvaluateOneDSpeciesInlet(configuration, species,
							options.case_directory, time);
			if (configuration.physiology.enabled) {
				state.has_hematocrit = true;
				state.hematocrit_percent = configuration.physiology.hematocrit_percent;
			}
			return state;
		};
		double initial_inlet = iga::EvaluateOneDInlet(configuration, inlet,
			options.case_directory, 0.0, inlet_area);
		iga::VascularInletState initial_port;
		if (replay) {
			initial_port = replay->StateAt(0.0);
			initial_inlet = iga::ApplyOneDCoupledInlet(configuration, transports, initial_port);
		} else if (circuit) {
			initial_port = circuit->InletState(0.0);
			initial_inlet = iga::ApplyOneDCoupledInlet(configuration, transports, initial_port);
		} else initial_port = open_loop_inlet(0.0, initial_inlet);
		if (flow.model == iga::OneDFlowModel::Rigid)
			iga::SolveRigidOneD(network, flow, flow_state, initial_inlet, configuration.time.dt);
		else iga::InitializeCompliantOneDFromRigid(network, flow, flow_state,
			initial_inlet, configuration.time.dt);
		if (!options.restart.empty()) {
			const auto metadata = iga::ReadOneDCheckpoint(options.restart, flow_state, transports,
				network, flow.dynamic_viscosity);
			ValidateRestart(metadata, configuration, network, transports, config_fingerprint);
		}
		if (options.stop_after_step > configuration.time.steps)
			throw std::runtime_error("--stop-after-step exceeds configured steps");
		if (options.output_directory.empty())
			options.output_directory = options.case_directory/"results"/"one_d"/flow.name;
		std::unique_ptr<iga::OneDOutputWriter> writer;
		std::unique_ptr<iga::CouplingHistoryWriter> coupling_writer;
		if (rank == 0) {
			iga::WriteOneDSkeletonFiles(options.output_directory, network,
				configuration.geometry.length_scale_to_m);
			writer = std::make_unique<iga::OneDOutputWriter>(
				options.output_directory, network, flow);
		}
		if (rank == 0 && configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly)
			coupling_writer = std::make_unique<iga::CouplingHistoryWriter>(
				options.output_directory, configuration.coupling.mode);
		auto derived = iga::ComputeOneDDerivedFields(configuration, transports);
		auto initial_result = iga::BuildOneDStepResult(configuration, network,
			flow_state, transports, initial_port, 0.0);
		auto previous_mass = initial_result.total_mass;
		if (coupling_writer) coupling_writer->Add(initial_result,
			iga::AggregateVascularOutlets(initial_result.outlets,
				configuration.coupling.flow_epsilon_m3_s));
		double output_seconds = 0.0;
		double solve_output_seconds = 0.0;
		if (rank == 0) {
			const auto before = std::chrono::steady_clock::now();
			writer->Write(flow_state.completed_step, flow_state.physical_time,
				flow_state, transports, derived);
			output_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
		}
		const auto setup_end = std::chrono::steady_clock::now();
		const auto solve_start = std::chrono::steady_clock::now();
		const int final_step = options.stop_after_step > 0
			? options.stop_after_step : configuration.time.steps;
		for (int step = flow_state.completed_step+1; step <= final_step; ++step) {
			const double time = step*configuration.time.dt;
			double inlet_flow = iga::EvaluateOneDInlet(configuration, inlet,
				options.case_directory, time, inlet_area);
			iga::VascularInletState port_inlet;
			if (replay) {
				port_inlet = replay->StateAt(time);
				inlet_flow = iga::ApplyOneDCoupledInlet(
					configuration, transports, port_inlet);
			} else if (circuit) {
				port_inlet = circuit->InletState(time);
				inlet_flow = iga::ApplyOneDCoupledInlet(
					configuration, transports, port_inlet);
			} else port_inlet = open_loop_inlet(time, inlet_flow);
			if (flow.scheme == iga::OneDFlowScheme::SteadyPoiseuille)
				iga::SolveRigidOneD(network, flow, flow_state, inlet_flow, configuration.time.dt);
			else if (flow.scheme == iga::OneDFlowScheme::ExplicitRusanov)
				iga::AdvanceExplicitOneD(network, flow, flow_state, inlet_flow, configuration.time.dt);
			else iga::AdvanceImplicitOneD(network, flow, flow_state, inlet_flow, configuration.time.dt);
			for (auto& transport : transports)
				iga::AdvanceOneDTransport(configuration, network, flow_state, transport,
					options.case_directory, (step-1)*configuration.time.dt, configuration.time.dt);
			iga::ApplyOneDVasodilation(configuration, network, transports,
				configuration.time.dt, flow.dynamic_viscosity);
			flow_state.completed_step = step;
			flow_state.physical_time = time;
			derived = iga::ComputeOneDDerivedFields(configuration, transports);
			if (configuration.coupling.mode != iga::SimulationScopeMode::FlowOnly) {
				auto coupled_result = iga::BuildOneDStepResult(configuration, network,
					flow_state, transports, port_inlet, configuration.time.dt,
					previous_mass);
				previous_mass = coupled_result.total_mass;
				auto venous = iga::AggregateVascularOutlets(coupled_result.outlets,
					configuration.coupling.flow_epsilon_m3_s);
				iga::CircuitAdvanceReport circuit_report;
				if (circuit) circuit_report = circuit->Advance(
					venous, configuration.time.dt, time);
				if (coupling_writer) coupling_writer->Add(
					std::move(coupled_result), std::move(venous),
					std::move(circuit_report));
			}
			if (rank == 0 && (step%configuration.time.output_every == 0
				|| step == final_step)) {
				const auto before = std::chrono::steady_clock::now();
				writer->Write(step, time, flow_state, transports, derived);
				const double elapsed = std::chrono::duration<double>(
					std::chrono::steady_clock::now()-before).count();
				output_seconds += elapsed;
				solve_output_seconds += elapsed;
			}
			if (!options.checkpoint.empty() && options.checkpoint_every > 0
				&& (step%options.checkpoint_every == 0 || step == final_step)) {
				fs::create_directories(options.checkpoint.parent_path().empty()
					? fs::path(".") : options.checkpoint.parent_path());
				iga::WriteOneDCheckpoint(options.checkpoint,
					CheckpointMetadata(configuration, network, flow_state, transports,
						config_fingerprint, options.checkpoint), flow_state, transports, network, rank);
			}
		}
		const auto solve_end = std::chrono::steady_clock::now();
		if (rank == 0) {
			if (coupling_writer) coupling_writer->Write();
			iga::WriteOneDPhysiologyManifest(options.output_directory,
				configuration, transports, derived);
			writer->Finish(flow_state,
				std::chrono::duration<double>(setup_end-setup_start).count(),
				std::chrono::duration<double>(solve_end-solve_start).count()-solve_output_seconds,
				output_seconds);
			std::cout << "completed 1d system=" << flow.name << " steps="
				<< flow_state.completed_step << " output=" << options.output_directory << '\n';
		}
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}

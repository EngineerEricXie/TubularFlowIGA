#include "OneDCheckpoint.hpp"
#include "OneDImplicit.hpp"
#include "OneDOutput.hpp"
#include "OneDRuntime.hpp"
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
		const auto& selected_flow = SelectFlow(configuration, options.system);
		auto network = iga::ReadOneDNetwork(options.case_directory/configuration.geometry.file,
			configuration.geometry.length_scale_to_m, selected_flow.discretization.cells_per_segment,
			selected_flow.dynamic_viscosity, configuration.geometry.root_node_id);
		iga::ValidateOneDTopologyReferences(configuration, network);
		const auto inlet = iga::ResolveOneDInlet(configuration);
		iga::OneDFlowRuntime runtime(configuration, selected_flow, std::move(network), inlet,
			options.case_directory, [](const iga::OneDNetwork& runtime_network,
				const iga::OneDFlowSystemDefinition& runtime_flow, iga::OneDFlowState& state,
				double inlet_flow, double dt) {
				iga::AdvanceImplicitOneD(runtime_network, runtime_flow, state, inlet_flow, dt);
			});
		const auto& flow = runtime.FlowSystem();
		std::unique_ptr<iga::ReplayInletProvider> replay;
		std::unique_ptr<iga::VcaExternalCircuit> circuit;
		if (runtime.Configuration().coupling.mode == iga::SimulationScopeMode::VcaReplay)
			replay = std::make_unique<iga::ReplayInletProvider>(
				iga::ReplayInletProvider::Read(options.case_directory
					/runtime.Configuration().coupling.replay_file));
		if (runtime.Configuration().coupling.mode == iga::SimulationScopeMode::VcaClosedLoop)
			circuit = std::make_unique<iga::VcaExternalCircuit>(runtime.Configuration().coupling);
		if (circuit && (!options.restart.empty() || !options.checkpoint.empty()))
			throw std::runtime_error(
				"closed-loop checkpoint/restart requires coupled reservoir state and is not yet enabled");
		const auto config_fingerprint = iga::OneDFingerprint(config_text);
		if (options.check) {
			if (rank == 0) std::cout << "schema_version=3 dimension=1d system=" << flow.name
				<< " nodes=" << runtime.Network().nodes.size()
				<< " segments=" << runtime.Network().segments.size()
				<< " root_id=" << runtime.Network().nodes[
					static_cast<std::size_t>(runtime.Network().root)].id
				<< " cells=" << runtime.Network().cells
				<< " outlets=" << runtime.Network().outlet_nodes.size()
				<< " transport_systems=" << runtime.Transports().size() << '\n';
			PetscFinalize();
			return 0;
		}
		const double inlet_area = runtime.Network().segments.front().area0;
		double initial_inlet = iga::EvaluateOneDInlet(runtime.Configuration(), inlet,
			options.case_directory, 0.0, inlet_area);
		iga::VascularInletState initial_port;
		if (replay) {
			initial_port = replay->StateAt(0.0);
			runtime.InitializeCoupled(initial_port);
		} else if (circuit) {
			initial_port = circuit->InletState(0.0);
			runtime.InitializeCoupled(initial_port);
		} else {
			initial_port = runtime.OpenLoopInlet(0.0, initial_inlet);
			runtime.InitializeOpenLoop(initial_port);
		}
		if (!options.restart.empty()) {
			auto restart_flow = runtime.FlowState();
			auto restart_transports = runtime.Transports();
			auto restart_network = runtime.Network();
			const auto metadata = iga::ReadOneDCheckpoint(options.restart, restart_flow,
				restart_transports, restart_network, flow.dynamic_viscosity);
			ValidateRestart(metadata, runtime.Configuration(), restart_network,
				restart_transports, config_fingerprint);
			runtime.RestoreCommittedState(std::move(restart_flow), std::move(restart_transports),
				std::move(restart_network));
		}
		if (options.stop_after_step > runtime.Configuration().time.steps)
			throw std::runtime_error("--stop-after-step exceeds configured steps");
		if (options.output_directory.empty())
			options.output_directory = options.case_directory/"results"/"one_d"/flow.name;
		std::unique_ptr<iga::OneDOutputWriter> writer;
		std::unique_ptr<iga::CouplingHistoryWriter> coupling_writer;
		if (rank == 0) {
			iga::WriteOneDSkeletonFiles(options.output_directory, runtime.Network(),
				runtime.Configuration().geometry.length_scale_to_m);
			writer = std::make_unique<iga::OneDOutputWriter>(
				options.output_directory, runtime.Network(), flow);
		}
		if (rank == 0 && runtime.Configuration().coupling.mode != iga::SimulationScopeMode::FlowOnly)
			coupling_writer = std::make_unique<iga::CouplingHistoryWriter>(
				options.output_directory, runtime.Configuration().coupling.mode);
		auto derived = iga::ComputeOneDDerivedFields(runtime.Configuration(), runtime.Transports());
		auto initial_result = iga::BuildOneDStepResult(runtime.Configuration(), runtime.Network(),
			runtime.FlowState(), runtime.Transports(), initial_port, 0.0);
		auto previous_mass = initial_result.total_mass;
		if (coupling_writer) coupling_writer->Add(initial_result,
			iga::AggregateVascularOutlets(initial_result.outlets,
				runtime.Configuration().coupling.flow_epsilon_m3_s));
		double output_seconds = 0.0;
		double solve_output_seconds = 0.0;
		if (rank == 0) {
			const auto before = std::chrono::steady_clock::now();
			writer->Write(runtime.FlowState().completed_step, runtime.FlowState().physical_time,
				runtime.FlowState(), runtime.Transports(), derived);
			output_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
		}
		const auto setup_end = std::chrono::steady_clock::now();
		const auto solve_start = std::chrono::steady_clock::now();
		const int final_step = options.stop_after_step > 0
			? options.stop_after_step : runtime.Configuration().time.steps;
		for (int step = runtime.FlowState().completed_step+1; step <= final_step; ++step) {
			const double time = step*runtime.Configuration().time.dt;
			double inlet_flow = iga::EvaluateOneDInlet(runtime.Configuration(), inlet,
				options.case_directory, time, inlet_area);
			iga::VascularInletState port_inlet;
			runtime.BeginStep(runtime.FlowState().physical_time, runtime.Configuration().time.dt);
			if (replay) {
				port_inlet = replay->StateAt(time);
				runtime.SetCoupledInlet(port_inlet);
			} else if (circuit) {
				port_inlet = circuit->InletState(time);
				runtime.SetCoupledInlet(port_inlet);
			} else {
				port_inlet = runtime.OpenLoopInlet(time, inlet_flow);
				runtime.SetOpenLoopInlet(port_inlet);
			}
			runtime.SolveTrial();
			runtime.CommitStep();
			derived = iga::ComputeOneDDerivedFields(runtime.Configuration(), runtime.Transports());
			if (runtime.Configuration().coupling.mode != iga::SimulationScopeMode::FlowOnly) {
				auto coupled_result = iga::BuildOneDStepResult(runtime.Configuration(),
					runtime.Network(), runtime.FlowState(), runtime.Transports(), port_inlet,
					runtime.Configuration().time.dt, previous_mass);
				previous_mass = coupled_result.total_mass;
				auto venous = iga::AggregateVascularOutlets(coupled_result.outlets,
					runtime.Configuration().coupling.flow_epsilon_m3_s);
				iga::CircuitAdvanceReport circuit_report;
				if (circuit) circuit_report = circuit->Advance(
					venous, runtime.Configuration().time.dt, time);
				if (coupling_writer) coupling_writer->Add(
					std::move(coupled_result), std::move(venous),
					std::move(circuit_report));
			}
			if (rank == 0 && (step%runtime.Configuration().time.output_every == 0
				|| step == final_step)) {
				const auto before = std::chrono::steady_clock::now();
				writer->Write(step, time, runtime.FlowState(), runtime.Transports(), derived);
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
					CheckpointMetadata(runtime.Configuration(), runtime.Network(),
						runtime.FlowState(), runtime.Transports(), config_fingerprint,
						options.checkpoint), runtime.FlowState(), runtime.Transports(),
					runtime.Network(), rank);
			}
		}
		const auto solve_end = std::chrono::steady_clock::now();
		if (rank == 0) {
			if (coupling_writer) coupling_writer->Write();
			iga::WriteOneDPhysiologyManifest(options.output_directory,
				runtime.Configuration(), runtime.Transports(), derived);
			writer->Finish(runtime.FlowState(),
				std::chrono::duration<double>(setup_end-setup_start).count(),
				std::chrono::duration<double>(solve_end-solve_start).count()-solve_output_seconds,
				output_seconds);
			std::cout << "completed 1d system=" << flow.name << " steps="
				<< runtime.FlowState().completed_step << " output=" << options.output_directory << '\n';
		}
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}

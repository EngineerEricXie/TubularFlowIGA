#define main FlowCliMain
#include "../src/iga_navier_stokes.cpp"
#undef main

namespace {

void CheckStep(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

template <class Work>
void StepFailure(MPI_Comm communicator, const char* stage, Work&& work)
{
	std::string diagnostic;
	try { work(); } catch (const std::exception& error) { diagnostic = error.what(); }
	iga::CollectiveLocalStage(communicator, "step test assertion", [&] {
		CheckStep(diagnostic.find(std::string(stage)+": rank ") == 0, "wrong step failure diagnostic");
	});
	iga::RequireCollectiveSameText(communicator, "step test common diagnostic", diagnostic);
}

void TestFlowStep(MPI_Comm communicator, const fs::path& fixture, const fs::path& directory, double concentration)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	std::optional<iga::SimulationConfiguration> base;
	const auto local = directory/("rank-"+std::to_string(rank));
	iga::CollectiveLocalStage(communicator, "step test fixture", [&] {
		CheckStep(fs::create_directories(local), "test output must be new");
		base.emplace(iga::ReadSimulationConfiguration((fixture/"simulation_config.json").string()));
		base->coupling.external_circuit.reservoir.species["oxygen"] = concentration;
		iga::TemporalFunctionDefinition table;
		table.name = "step_table";
		table.kind = iga::TemporalFunctionKind::PeriodicTable;
		table.period = 1;
		table.file = "table.csv";
		base->temporal_functions = {table};
		base->boundaries.front().conditions.front().waveform = table.name;
		base->boundaries.front().conditions.front().scale = 0.001;
		std::ofstream output(local/table.file);
		output << "time,value\n0,1\n0.04,1.2\n";
		output.close();
		CheckStep(bool(output), "cannot write test table");
	});
	const auto system = iga::RequireThreeDVcaTransportSystem(*base);
	auto healthy_input = [&] {
		iga::VcaExternalCircuit circuit(base->coupling);
		const auto input = PrepareFlowStepInput(communicator, "flow step input", true, true,
			*base, local, 0.01, 0.0, &circuit, &system, -1.0);
		iga::CollectiveLocalStage(communicator, "step test inlet values", [&] {
			CheckStep(input.inlet.time_s == 0.0 && input.inlet.species.at("oxygen") == concentration,
				"VCA inlet lag or species changed");
			const auto& condition = input.configuration.boundaries.front().conditions.front();
			CheckStep(condition.waveform.empty() && condition.scale == input.inlet.flow_m3_s,
				"VCA reference profile scale changed");
		});
		return input.inlet;
	};
	const auto inlet = healthy_input();
	int cases = 0;
	for (const std::string mode : {"unknown-table", "missing-table", "bad-table", "absolute-table",
		"missing-flow", "zero-reference", "nan-reference", "missing-species", "missing-species-boundary", "nan-inlet-time"}) {
		auto configuration = *base;
		auto transport = system;
		double reference = -1, inlet_time = 0;
		if (rank == ranks-1) {
			if (mode == "unknown-table") configuration.boundaries.front().conditions.front().waveform = "absent";
			if (mode == "missing-table") configuration.temporal_functions.front().file = "absent.csv";
			if (mode == "bad-table") {
				configuration.temporal_functions.front().file = "bad.csv";
				std::ofstream output(local/"bad.csv"); output << "malformed\n";
			}
			if (mode == "absolute-table") configuration.temporal_functions.front().file = (local/"table.csv").string();
			if (mode == "missing-flow") configuration.equation_systems.erase(configuration.equation_systems.begin());
			if (mode == "zero-reference") reference = 0;
			if (mode == "nan-reference") reference = std::numeric_limits<double>::quiet_NaN();
			if (mode == "missing-species") transport.field_index.erase("oxygen");
			if (mode == "missing-species-boundary") configuration.boundaries.front().conditions.pop_back();
			if (mode == "nan-inlet-time") inlet_time = std::numeric_limits<double>::quiet_NaN();
		}
		iga::VcaExternalCircuit circuit(base->coupling);
		const char* stage = cases%2 == 0 ? "flow initial VCA input" : "flow step input";
		StepFailure(communicator, stage, [&] {
			PrepareFlowStepInput(communicator, stage, true, true, configuration, local,
				0.01, inlet_time, &circuit, &transport, reference);
		});
		healthy_input();
		++cases;
	}
	// Standalone flow materializes the current physical time; steady and legacy
	// paths retain their existing copy/default behavior.
	const auto standalone = PrepareFlowStepInput(communicator, "flow step input", true, true,
		*base, local, .01, 0, nullptr, nullptr, 0);
	const auto steady = PrepareFlowStepInput(communicator, "flow step input", true, false,
		*base, local, .01, 0, nullptr, nullptr, 0);
	const auto legacy = PrepareFlowStepInput(communicator, "flow step input", false, false,
		*base, local, .01, 0, nullptr, nullptr, 0);
	iga::CollectiveLocalStage(communicator, "step test standalone modes", [&] {
		CheckStep(std::abs(standalone.configuration.boundaries.front().conditions.front().scale-.00105) < 1e-15,
			"waveform used the wrong physical time");
		CheckStep(steady.configuration.boundaries.front().conditions.front().waveform == "step_table",
			"steady flow unexpectedly materialized a waveform");
		CheckStep(legacy.configuration.boundaries.empty(), "legacy input unexpectedly copied configuration");
	});
	iga::FlowPortMeasurements ports;
	for (int label : base->coupling.three_d_ports.outlet_labels) {
		const double flow = inlet.flow_m3_s/base->coupling.three_d_ports.outlet_labels.size();
		ports.flows[label] = flow;
		ports.pressures[label] = 1.0;
		ports.species_fluxes[label]["oxygen"] = concentration*flow;
		ports.species_concentrations[label]["oxygen"] = concentration;
	}
	const auto epsilon = base->coupling.flow_epsilon_m3_s;
	auto healthy_result = [&] {
		return PrepareFlowVcaResult(communicator, .01, .01, inlet, base->coupling.three_d_ports, ports, epsilon);
	};
	for (const std::string mode : {"nan-time", "negative-dt", "invalid-inlet", "missing-flow", "missing-pressure"}) {
		auto values = ports;
		auto input = inlet;
		double time = .01, dt = .01;
		if (rank == ranks-1) {
			if (mode == "nan-time") time = std::numeric_limits<double>::quiet_NaN();
			if (mode == "negative-dt") dt = -1;
			if (mode == "invalid-inlet") input.has_flow = false;
			if (mode == "missing-flow") values.flows.clear();
			if (mode == "missing-pressure") values.pressures.clear();
		}
		StepFailure(communicator, "flow VCA result preparation", [&] {
			PrepareFlowVcaResult(communicator, time, dt, input, base->coupling.three_d_ports, values, epsilon);
		});
		healthy_result();
		++cases;
	}
	for (const std::string mode : {"empty-outlets", "duplicate-outlets", "nan-flow", "negative-epsilon",
		"negative-dt", "empty-reservoir", "negative-species-mass"}) {
		auto result = healthy_result();
		double threshold = epsilon;
		if (rank == ranks-1) {
			if (mode == "empty-outlets") result.outlets.clear();
			if (mode == "duplicate-outlets") result.outlets.push_back(result.outlets.front());
			if (mode == "nan-flow") result.outlets.front().flow_m3_s = std::numeric_limits<double>::quiet_NaN();
			if (mode == "negative-epsilon") threshold = -1;
			if (mode == "negative-dt") result.dt_s = -1;
			if (mode == "empty-reservoir") result.outlets.front().flow_m3_s = -1e6;
			if (mode == "negative-species-mass") result.outlets.front().species_flux["oxygen"] = -1e6;
		}
		iga::VcaExternalCircuit circuit(base->coupling);
		iga::CouplingHistoryWriter history(local, base->coupling.mode);
		StepFailure(communicator, "flow VCA circuit advance", [&] {
			AdvanceFlowVcaCircuit(communicator, circuit, history, result, threshold);
		});
		// The native CLI exits on a failed coupling step. Reconstruct here; this
		// test does not claim an in-place rollback of partially advanced circuits.
		iga::VcaExternalCircuit retry(base->coupling);
		iga::CouplingHistoryWriter retry_history(local, base->coupling.mode);
		const auto accepted = healthy_result();
		AdvanceFlowVcaCircuit(communicator, retry, retry_history, accepted, epsilon);
		iga::CollectiveLocalStage(communicator, "step test circuit values", [&] {
			CheckStep(std::abs(retry.State().species.at("oxygen")-concentration) < 1e-12,
				"closed-loop species conservation failed");
			CheckStep(std::abs(retry.State().volume_m3-base->coupling.external_circuit.reservoir.volume_m3) < 1e-12,
				"closed-loop volume conservation failed");
		});
		++cases;
	}
	// Unequal work counts make an accidental WORLD collective in a helper
	// observable while the independent pair finishes its own communicator.
	if (ranks == 1) healthy_input();
	if (rank == 0) std::cout << "flow_step ranks=" << ranks << " failures=" << cases
		<< " healthy_retries=" << cases << " passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		CheckStep(argc == 3 && ranks == 3, "require 3 ranks, VCA fixture, and new output directory");
		const fs::path output = fs::absolute(argv[2]);
		TestFlowStep(PETSC_COMM_WORLD, argv[1], output/"world", .2);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		TestFlowStep(group, argv[1], output/(rank == 0 ? "self" : "pair"), rank == 0 ? .3 : .4);
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

#include "CouplingReplay.hpp"
#include "OneDCoupling.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;

namespace {

bool Close(double first, double second, double tolerance = 1.0e-12)
{
	return std::abs(first-second) <= tolerance*std::max({1.0, std::abs(first), std::abs(second)});
}

iga::CouplingDefinition PfcClosedLoopDefinition()
{
	iga::CouplingDefinition definition;
	definition.mode = iga::SimulationScopeMode::VcaClosedLoop;
	definition.perfusate.type = "pfc";
	definition.perfusate.oxygen.hematocrit_percent = 0.0;
	definition.perfusate.oxygen.hemoglobin_g_dl = 0.0;
	definition.perfusate.oxygen.oxygen_solubility_ml_dl_mmhg = 0.034;
	definition.external_circuit.reservoir.enabled = true;
	definition.external_circuit.reservoir.volume_m3 = 1.0;
	definition.external_circuit.reservoir.temperature_c = 37.0;
	definition.external_circuit.pump.mode = "flow_control";
	definition.external_circuit.pump.flow_m3_s = 0.1;
	return definition;
}

} // namespace

int main()
{
	iga::OxygenCapacityParameters oxygen;
	oxygen.hematocrit_percent = 40.0;
	oxygen.hemoglobin_g_dl = 13.6;
	const auto equilibrium = iga::OxygenFromPo2(100.0, oxygen);
	assert(equilibrium.total_oxygen_mol_m3 > equilibrium.dissolved_oxygen_mol_m3);
	assert(Close(iga::Po2FromTotalOxygen(equilibrium.total_oxygen_mol_m3, oxygen), 100.0, 1.0e-9));

	const auto venous = iga::AggregateVascularOutlets({
		{2, 2.0, 10.0, {{"glucose", 8.0}}, {}, {}, true},
		{3, 1.0, 4.0, {{"glucose", 1.0}}, {}, {}, true},
	});
	assert(Close(venous.flow_m3_s, 3.0));
	assert(Close(venous.pressure_pa, 8.0));
	assert(Close(venous.flux_weighted_concentration.at("glucose"), 3.0));

	{
		iga::CouplingDefinition definition;
		definition.mode = iga::SimulationScopeMode::VcaReplay;
		bool rejected = false;
		try { iga::RequireFlowOnlyCoupling(definition, "test backend"); }
		catch (const std::runtime_error&) { rejected = true; }
		assert(rejected);
		definition.mode = iga::SimulationScopeMode::FlowOnly;
		iga::RequireFlowOnlyCoupling(definition, "test backend");
	}

	{
		auto definition = PfcClosedLoopDefinition();
		definition.external_circuit.reservoir.species = {{"oxygen", 1.0}};
		definition.external_circuit.oxygenator.enabled = true;
		definition.external_circuit.oxygenator.po2_mmhg = 100.0;
		iga::VcaExternalCircuit circuit(definition);
		const auto arterial = circuit.InletState(0.0);
		assert(arterial.species.at("oxygen") > 1.0);
		const auto returned = iga::AggregateVascularOutlets({
			{2, 0.1, 0.0, {{"oxygen", 0.1}}, {}, {}, true},
		});
		const auto report = circuit.Advance(returned, 1.0, 1.0);
		assert(Close(circuit.State().species.at("oxygen"), 1.0));
		assert(report.device_source_rate.at("oxygen") > 0.0);
	}

	{
		iga::OneDConfiguration configuration;
		configuration.physiology.hematocrit_percent = 40.0;
		configuration.physiology.hemoglobin_g_dl = 13.6;
		configuration.coupling.perfusate.oxygen.hematocrit_percent = 40.0;
		configuration.coupling.perfusate.oxygen.hemoglobin_g_dl = 13.6;
		iga::VascularInletState inlet;
		inlet.time_s = 0.0;
		inlet.has_flow = true;
		inlet.flow_m3_s = 1.0e-9;
		inlet.has_hematocrit = true;
		inlet.hematocrit_percent = 20.0;
		std::vector<iga::OneDTransportState> transports;
		iga::ApplyOneDCoupledInlet(configuration, transports, inlet);
		assert(Close(configuration.physiology.hematocrit_percent, 20.0));
		assert(Close(configuration.physiology.hemoglobin_g_dl, 6.8));
	}

	{
		iga::OneDConfiguration configuration;
		configuration.physiology.enabled = true;
		configuration.physiology.metabolism_rates["glucose"] = -1.0;
		iga::OneDNetwork network;
		network.cells = 1;
		network.segments.push_back({0, 0, 1, 1.0, 0.0, 0.0, 0.0, 0.0, {}, 0, 1});
		iga::OneDFlowState flow;
		flow.area = {1.0};
		iga::OneDSpeciesState glucose;
		glucose.definition.field = "glucose";
		glucose.definition.volume_source = 3.0;
		glucose.definition.reaction_rate = 0.5;
		glucose.concentration = {2.0};
		glucose.wall_kind = iga::OneDWallBoundaryKind::ConstantFlux;
		glucose.wall_value = 0.25;
		const double expected = 1.0-0.5*std::sqrt(iga::OneDPi);
		assert(Close(iga::OneDSpeciesSourceIntegral(configuration, network, flow, glucose), expected));
	}

	{
		const auto path = fs::temp_directory_path()/"tubularflowiga-vca-replay.csv";
		{
			std::ofstream output(path);
			output << "time_s,flow_m3_s,species.glucose\n0,1,2\n1,2,3\n";
		}
		const auto replay = iga::ReplayInletProvider::Read(path);
		const auto state = replay.StateAt(0.5);
		assert(Close(state.flow_m3_s, 1.0));
		assert(Close(state.species.at("glucose"), 2.0));
		fs::remove(path);

		const auto json_path = fs::temp_directory_path()/"tubularflowiga-vca-replay.json";
		{
			std::ofstream output(json_path);
			output << R"json({"states":[{"time_s":0,"flow_m3_s":1,"species":{"glucose":2}},{"time_s":1,"flow_m3_s":2,"species":{"glucose":3}}]})json";
		}
		const auto json_replay = iga::ReplayInletProvider::Read(json_path);
		const auto json_state = json_replay.StateAt(1.0);
		assert(Close(json_state.flow_m3_s, 2.0));
		assert(Close(json_state.species.at("glucose"), 3.0));
		fs::remove(json_path);
	}

	{
		const auto directory = fs::temp_directory_path()/"tubularflowiga-vca-manifest";
		fs::create_directories(directory);
		iga::VascularStepResult step;
		step.time_s = 1.0;
		step.dt_s = 1.0;
		step.inlet.time_s = 1.0;
		step.inlet.has_flow = true;
		step.inlet.flow_m3_s = 0.1;
		step.inlet.species = {{"glucose", 2.0}};
		step.total_mass = {{"glucose", 1.0}};
		step.source_integrals = {{"glucose", -0.2}};
		step.balance_residuals = {{"glucose", 0.0}};
		iga::AggregatedVascularReturn return_state;
		return_state.flow_m3_s = 0.1;
		return_state.pressure_valid = true;
		return_state.pressure_pa = 4.0;
		return_state.species_flux = {{"glucose", 0.1}};
		return_state.flux_weighted_concentration = {{"glucose", 1.0}};
		return_state.outlet_ids = {2};
		iga::CircuitAdvanceReport report;
		report.time_s = 1.0;
		report.species_mass_change = {{"glucose", -0.1}};
		report.device_source_rate = {{"glucose", 0.0}};
		iga::CouplingHistoryWriter writer(directory, iga::SimulationScopeMode::VcaClosedLoop);
		writer.Add(step, return_state, report);
		writer.Write();
		std::ifstream input(directory/"coupling_manifest.json");
		std::ostringstream text;
		text << input.rdbuf();
		assert(text.str().find("\"units\"") != std::string::npos);
		assert(text.str().find("aggregated_venous_return_history") != std::string::npos);
		assert(text.str().find("vascular_source_integrals_mol_s") != std::string::npos);
		fs::remove(directory/"coupling_manifest.json");
		fs::remove(directory);
	}

	std::cout << "one-dimensional coupling tests passed\n";
}

#include "OneDRuntime.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

namespace {

bool Close(double first, double second)
{
	return std::abs(first-second) <= 1.0e-12*std::max({1.0, std::abs(first), std::abs(second)});
}

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
}

std::string Configuration()
{
	return R"json({
 "schema_version":3,"dimension":"1d",
 "geometry":{"kind":"swc_network","file":"tree.swc","length_scale_to_m":1.0},
 "fields":[{"name":"area","kind":"scalar"},{"name":"flow_rate","kind":"scalar"},{"name":"pressure","kind":"pressure"},{"name":"signal","kind":"scalar","initial_value":1.0}],
 "time":{"dt":0.001,"steps":2,"output_every":1},
 "equation_systems":[
  {"name":"flow","kind":"network_flow_1d","unknowns":["area","flow_rate","pressure"],"model":"rigid","scheme":"steady_poiseuille","dynamic_viscosity":0.004,"density":1060.0,"discretization":{"cells_per_segment":1}},
  {"name":"transport","kind":"network_transport_1d","unknowns":["signal"],"flow_system":"flow","species":[{"field":"signal","diffusivity":0.0}]}
 ],
 "boundaries":[
  {"name":"inlet","role":"inlet","conditions":[{"field":"flow_rate","type":"dirichlet","value":1e-9},{"field":"signal","type":"dirichlet","value":1.0}]},
  {"name":"outlet","role":"outlet","conditions":[{"field":"pressure","type":"pressure","value":0.0}]}
 ],
 "physiology":{"enabled":true,"oxygen_capacity":{"enabled":true,"hematocrit_percent":40.0,"hemoglobin_g_dl":13.6},"vasodilation":{"enabled":true,"field":"signal","emax_radius_fraction":0.1,"ec50":1.0,"relaxation_tau":0.01}}
})json";
}

} // namespace

int main()
{
	const auto directory = fs::temp_directory_path()/"tubularflowiga-one-d-runtime-test";
	fs::create_directories(directory);
	{
		std::ofstream tree(directory/"tree.swc");
		tree << "1 2 0 0 0 0.001 -1\n2 2 0.01 0 0 0.001 1\n";
	}
	auto configuration = iga::ParseOneDConfiguration(Configuration());
	const auto flow = configuration.flow_systems.front();
	auto network = iga::ReadOneDNetwork(directory/"tree.swc", 1.0, 1, flow.dynamic_viscosity);
	iga::OneDFlowRuntime runtime(configuration, flow, network,
		iga::ResolveOneDInlet(configuration), directory);
	runtime.InitializeOpenLoop(1.0e-9);
	RequireRejected([&runtime] { runtime.CommitStep(); });
	const double initial_radius = runtime.Network().segments.front().radius0;
	const double initial_perfusate_hematocrit
		= runtime.Configuration().coupling.perfusate.oxygen.hematocrit_percent;
	const double initial_transport_inlet
		= runtime.Transports().front().species.front().inlet_value;
	runtime.BeginStep(0.0, configuration.time.dt);
	iga::VascularInletState inlet;
	inlet.time_s = configuration.time.dt;
	inlet.has_flow = true;
	inlet.flow_m3_s = 2.0e-9;
	inlet.has_hematocrit = true;
	inlet.hematocrit_percent = 20.0;
	inlet.species.emplace("signal", 2.0);
	runtime.SetCoupledInlet(inlet);
	runtime.SolveTrial();
	const auto first_root = runtime.GetPortState("root");
	const auto first_outlet = runtime.GetPortState("outlet:2");
	assert(*first_root.outward_flow_m3_s < 0.0);
	assert(*first_outlet.outward_flow_m3_s > 0.0);
	assert(Close(first_root.concentration.at("signal"), 2.0));
	assert(Close(first_root.outward_species_flux.at("signal"), -4.0e-9));
	assert(runtime.FlowState().completed_step == 1);
	assert(runtime.FlowState().internal_substeps == 0);
	assert(!Close(runtime.Network().segments.front().radius0, initial_radius));
	const double first_radius = runtime.Network().segments.front().radius0;
	const auto first_area = runtime.FlowState().area;
	const auto first_flow = runtime.FlowState().flow;
	const auto first_pressure = runtime.FlowState().pressure;
	assert(Close(runtime.Configuration().physiology.hematocrit_percent, 20.0));
	runtime.RollbackTrial();
	assert(runtime.FlowState().completed_step == 0);
	assert(Close(runtime.Network().segments.front().radius0, initial_radius));
	assert(Close(runtime.Configuration().physiology.hematocrit_percent, 40.0));
	assert(Close(runtime.Configuration().coupling.perfusate.oxygen.hematocrit_percent,
		initial_perfusate_hematocrit));
	assert(Close(runtime.Transports().front().species.front().inlet_value,
		initial_transport_inlet));
	assert(Close(runtime.LastInlet().flow_m3_s, 1.0e-9));
	runtime.SolveTrial();
	const auto repeated_root = runtime.GetPortState("root");
	const auto repeated_outlet = runtime.GetPortState("outlet:2");
	assert(Close(*first_root.outward_flow_m3_s, *repeated_root.outward_flow_m3_s));
	assert(Close(*first_outlet.outward_flow_m3_s, *repeated_outlet.outward_flow_m3_s));
	assert(runtime.FlowState().area == first_area);
	assert(runtime.FlowState().flow == first_flow);
	assert(runtime.FlowState().pressure == first_pressure);
	assert(runtime.Network().segments.front().radius0 == first_radius);
	runtime.CommitStep();
	assert(runtime.FlowState().completed_step == 1);
	runtime.BeginStep(configuration.time.dt, configuration.time.dt);
	iga::PortBoundaryData generic_input;
	generic_input.time_s = 2.0*configuration.time.dt;
	generic_input.outward_flow_m3_s = -3.0e-9;
	generic_input.concentration.emplace("signal", 3.0);
	runtime.SetPortInput("root", generic_input);
	iga::PortBoundaryData outlet_pressure;
	outlet_pressure.time_s = 2.0*configuration.time.dt;
	outlet_pressure.mean_pressure_pa = 5.0;
	runtime.SetPortInput("outlet:2", outlet_pressure);
	runtime.SolveTrial();
	assert(Close(*runtime.GetPortState("root").outward_flow_m3_s, -3.0e-9));
	assert(Close(*runtime.GetPortState("outlet:2").mean_pressure_pa, 5.0));
	RequireRejected([&runtime] { runtime.GetPortState("outlet:2x"); });
	runtime.CommitStep();
	assert(runtime.FlowState().completed_step == 2);

	auto failing_flow = flow;
	failing_flow.scheme = iga::OneDFlowScheme::ImplicitPetsc;
	iga::OneDFlowRuntime failing(configuration, failing_flow, network,
		iga::ResolveOneDInlet(configuration), directory,
		[](const iga::OneDNetwork&, const iga::OneDFlowSystemDefinition&,
			iga::OneDFlowState& state, double, double) {
			state.completed_step = 99;
			throw std::runtime_error("injected implicit failure");
		});
	failing.InitializeOpenLoop(1.0e-9);
	failing.BeginStep(0.0, configuration.time.dt);
	failing.SetOpenLoopInlet(failing.OpenLoopInlet(configuration.time.dt, 1.0e-9));
	RequireRejected([&failing] { failing.SolveTrial(); });
	assert(failing.CurrentPhase() == iga::OneDFlowRuntime::Phase::TrialSolved);
	RequireRejected([&failing] { failing.CommitStep(); });
	failing.RollbackTrial();
	assert(failing.FlowState().completed_step == 0);
	assert(failing.CurrentPhase() == iga::OneDFlowRuntime::Phase::TrialOpen);

	auto resistance_configuration = configuration;
	auto& resistance = resistance_configuration.boundaries.back().conditions.front();
	resistance.type = "resistance";
	resistance.resistance = 1.0e12;
	resistance.reference_pressure = 0.0;
	iga::OneDFlowRuntime resistance_runtime(resistance_configuration, flow, network,
		iga::ResolveOneDInlet(resistance_configuration), directory);
	resistance_runtime.InitializeOpenLoop(1.0e-9);
	resistance_runtime.BeginStep(0.0, configuration.time.dt);
	resistance_runtime.SetOpenLoopInlet(
		resistance_runtime.OpenLoopInlet(configuration.time.dt, 1.0e-9));
	iga::PortBoundaryData forbidden_pressure;
	forbidden_pressure.time_s = configuration.time.dt;
	forbidden_pressure.mean_pressure_pa = 5.0;
	RequireRejected([&resistance_runtime, &forbidden_pressure] {
		resistance_runtime.SetPortInput("outlet:2", forbidden_pressure);
	});
	fs::remove_all(directory);
	std::cout << "one-dimensional runtime tests passed\n";
}

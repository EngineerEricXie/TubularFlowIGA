#include "OneDRuntime.hpp"
#include "OneDFlowDomainAdapter.hpp"

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
	{
		std::ofstream tree(directory/"root-bifurcation.swc");
		tree << "1 2 0 0 0 0.001 -1\n"
			<< "2 2 0.01 0.01 0 0.001 1\n"
			<< "3 2 0.01 -0.01 0 0.001 1\n";
	}
	auto configuration = iga::ParseOneDConfiguration(Configuration());
	const auto flow = configuration.flow_systems.front();
	auto network = iga::ReadOneDNetwork(directory/"tree.swc", 1.0, 1, flow.dynamic_viscosity);
	{
		auto bifurcated_network = iga::ReadOneDNetwork(directory/"root-bifurcation.swc",
			1.0, 1, flow.dynamic_viscosity);
		assert(iga::OneDSegmentsOutOfNode(bifurcated_network,
			bifurcated_network.root).size() == 2);
		iga::OneDFlowRuntime bifurcated(configuration, flow, bifurcated_network,
			iga::ResolveOneDInlet(configuration), directory);
		bifurcated.InitializeOpenLoop(1.0e-9);
		bifurcated.BeginStep(0.0, configuration.time.dt);
		iga::VascularInletState inlet;
		inlet.time_s = configuration.time.dt;
		inlet.has_flow = true;
		inlet.flow_m3_s = 2.0e-9;
		inlet.species.emplace("signal", 2.0);
		bifurcated.SetCoupledInlet(inlet);
		bifurcated.SolveTrial();
		const auto& species = bifurcated.Transports().front().species.front();
		double expected_native_flux = 0.0;
		for (const int segment_index : iga::OneDSegmentsOutOfNode(
			bifurcated.Network(), bifurcated.Network().root)) {
			const auto& segment = bifurcated.Network().segments.at(
				static_cast<std::size_t>(segment_index));
			expected_native_flux += bifurcated.FlowState().flow.at(
				static_cast<std::size_t>(segment.cell_offset))*2.0;
		}
		assert(Close(species.root_native_flux, expected_native_flux));
		assert(Close(bifurcated.GetPortState("root").outward_species_flux.at("signal"),
			-expected_native_flux));
	}
	{
		auto macro_configuration = configuration;
		macro_configuration.time.steps = 4;
		iga::OneDFlowRuntime subcycled(macro_configuration, flow, network,
			iga::ResolveOneDInlet(macro_configuration), directory);
		iga::OneDFlowRuntime sequential(macro_configuration, flow, network,
			iga::ResolveOneDInlet(macro_configuration), directory);
		subcycled.InitializeOpenLoop(1.0e-9);
		sequential.InitializeOpenLoop(1.0e-9);
		subcycled.BeginStep(0.0, 4.0*macro_configuration.time.dt);
		subcycled.SetOpenLoopInlet(subcycled.OpenLoopInlet(4.0*macro_configuration.time.dt, 1.0e-9));
		subcycled.SolveTrial();
		assert(subcycled.Diagnostics().planned_configured_substeps == 4);
		assert(subcycled.Diagnostics().attempted_configured_substeps == 4);
		assert(subcycled.Diagnostics().completed_configured_substeps == 4);
		for (int step = 0; step < 4; ++step) {
			sequential.BeginStep(sequential.FlowState().physical_time, macro_configuration.time.dt);
			sequential.SetOpenLoopInlet(sequential.OpenLoopInlet((step+1)*macro_configuration.time.dt, 1.0e-9));
			sequential.SolveTrial();
			sequential.CommitStep();
		}
		assert(subcycled.FlowState().area == sequential.FlowState().area);
		assert(subcycled.FlowState().flow == sequential.FlowState().flow);
		assert(subcycled.FlowState().pressure == sequential.FlowState().pressure);
		assert(subcycled.FlowState().node_pressure == sequential.FlowState().node_pressure);
		assert(subcycled.FlowState().segment_flow == sequential.FlowState().segment_flow);
		assert(subcycled.FlowState().outlets.size() == sequential.FlowState().outlets.size());
		for (std::size_t i = 0; i < subcycled.FlowState().outlets.size(); ++i) {
			const auto& macro_outlet = subcycled.FlowState().outlets[i];
			const auto& sequential_outlet = sequential.FlowState().outlets[i];
			assert(Close(macro_outlet.pressure, sequential_outlet.pressure)
				&& Close(macro_outlet.flow, sequential_outlet.flow)
				&& Close(macro_outlet.capacitor_pressure, sequential_outlet.capacitor_pressure));
		}
		assert(subcycled.Transports().size() == sequential.Transports().size());
		for (std::size_t i = 0; i < subcycled.Transports().size(); ++i) {
			assert(subcycled.Transports()[i].species.size() == sequential.Transports()[i].species.size());
			for (std::size_t species = 0; species < subcycled.Transports()[i].species.size(); ++species) {
				const auto& macro_species = subcycled.Transports()[i].species[species];
				const auto& sequential_species = sequential.Transports()[i].species[species];
				assert(macro_species.concentration == sequential_species.concentration
					&& Close(macro_species.inlet_value, sequential_species.inlet_value));
			}
		}
		assert(subcycled.Network().segments.size() == sequential.Network().segments.size());
		for (std::size_t i = 0; i < subcycled.Network().segments.size(); ++i) {
			const auto& macro_segment = subcycled.Network().segments[i];
			const auto& sequential_segment = sequential.Network().segments[i];
			assert(Close(macro_segment.radius0, sequential_segment.radius0)
				&& Close(macro_segment.area0, sequential_segment.area0)
				&& Close(macro_segment.resistance, sequential_segment.resistance));
		}
		assert(Close(subcycled.Configuration().physiology.hematocrit_percent,
			sequential.Configuration().physiology.hematocrit_percent));
		assert(subcycled.FlowState().completed_step == 4);
		RequireRejected([&sequential] { sequential.BeginStep(sequential.FlowState().physical_time, sequential.Configuration().time.dt); });
		subcycled.RollbackTrial();
		subcycled.SolveTrial();
		assert(subcycled.FlowState().area == sequential.FlowState().area);
		assert(subcycled.FlowState().flow == sequential.FlowState().flow);
		assert(subcycled.Transports().front().species.front().concentration
			== sequential.Transports().front().species.front().concentration);
		assert(Close(subcycled.Network().segments.front().radius0,
			sequential.Network().segments.front().radius0));
		auto configured_configuration = configuration;
		configured_configuration.time.dt = 0.0025;
		configured_configuration.time.steps = 4;
		iga::OneDFlowRuntime configured(configured_configuration, flow, network,
			iga::ResolveOneDInlet(configured_configuration), directory);
		configured.InitializeOpenLoop(1.0e-9);
		configured.BeginStep(0.0, 0.01);
		configured.SetConfiguredOpenLoopInlet();
		configured.SolveTrial();
		const auto& endpoints = configured.Diagnostics().configured_open_loop_endpoint_times_s;
		assert(endpoints.size() == 4);
		assert(Close(endpoints[0], 0.0025) && Close(endpoints[1], 0.005)
			&& Close(endpoints[2], 0.0075) && Close(endpoints[3], 0.01));
		auto waveform_configuration = configured_configuration;
		iga::TemporalFunctionDefinition waveform;
		waveform.name = "pulse";
		waveform.kind = iga::TemporalFunctionKind::Sinusoid;
		waveform.units = "m3/s";
		waveform.mean = 1.0e-9;
		waveform.amplitude = 0.5e-9;
		waveform.period = 0.01;
		waveform_configuration.temporal_functions.push_back(waveform);
		waveform_configuration.boundaries.front().conditions.front().waveform = "pulse";
		iga::OneDFlowRuntime waveform_runtime(waveform_configuration, flow, network,
			iga::ResolveOneDInlet(waveform_configuration), directory);
		waveform_runtime.InitializeOpenLoop(1.0e-9);
		waveform_runtime.BeginStep(0.0, 0.01);
		waveform_runtime.SetConfiguredOpenLoopInlet();
		waveform_runtime.SolveTrial();
		const auto& waveform_flows = waveform_runtime.Diagnostics().configured_open_loop_endpoint_flows_m3_s;
		assert(waveform_flows.size() == 4);
		for (std::size_t i = 0; i < waveform_flows.size(); ++i) {
			const double time_s = (i+1)*0.0025;
			assert(Close(waveform_flows[i], 1.0e-9+0.5e-9*std::sin(2.0*std::acos(-1.0)*time_s/0.01)));
		}
		const auto first_waveform_flows = waveform_flows;
		waveform_runtime.RollbackTrial();
		waveform_runtime.SolveTrial();
		assert(waveform_runtime.Diagnostics().configured_open_loop_endpoint_flows_m3_s == first_waveform_flows);
		iga::OneDFlowRuntime held(configured_configuration, flow, network,
			iga::ResolveOneDInlet(configured_configuration), directory);
		held.InitializeOpenLoop(1.0e-9);
		held.BeginStep(0.0, 0.01);
		held.SetOpenLoopInlet(held.OpenLoopInlet(0.01, 2.0e-9));
		held.SolveTrial();
		assert(held.Diagnostics().configured_open_loop_endpoint_flows_m3_s.empty());
		assert(held.Diagnostics().substep_endpoint_flows_m3_s.size() == 4);
		for (const double inlet_flow : held.Diagnostics().substep_endpoint_flows_m3_s)
			assert(Close(inlet_flow, 2.0e-9));
	}
	{
		auto inertance_configuration = configuration;
		auto inertance_flow = flow;
		inertance_flow.scheme = iga::OneDFlowScheme::RigidInertance;
		iga::OneDFlowRuntime inertance_runtime(inertance_configuration, inertance_flow, network,
			iga::ResolveOneDInlet(inertance_configuration), directory);
		inertance_runtime.InitializeOpenLoop(1.0e-9);
		const auto committed = inertance_runtime.FlowState();
		inertance_runtime.BeginStep(0.0, inertance_configuration.time.dt);
		inertance_runtime.SetOpenLoopInlet(inertance_runtime.OpenLoopInlet(
			inertance_configuration.time.dt, 1.2e-9));
		inertance_runtime.SolveTrial();
		const auto trial = inertance_runtime.FlowState();
		assert(trial.segment_flow.front() > committed.segment_flow.front());
		assert(trial.node_pressure.front() > committed.node_pressure.front());
		inertance_runtime.RollbackTrial();
		assert(inertance_runtime.FlowState().segment_flow == committed.segment_flow);
		assert(inertance_runtime.FlowState().node_pressure == committed.node_pressure);
		inertance_runtime.SolveTrial();
		assert(inertance_runtime.FlowState().segment_flow == trial.segment_flow);
		assert(inertance_runtime.FlowState().node_pressure == trial.node_pressure);
		inertance_runtime.CommitStep();
		assert(inertance_runtime.FlowState().completed_step == 1);
	}
	{
		auto explicit_configuration = configuration;
		explicit_configuration.time.dt = 0.0025;
		explicit_configuration.time.steps = 4;
		auto explicit_flow = flow;
		explicit_flow.model = iga::OneDFlowModel::Compliant;
		explicit_flow.scheme = iga::OneDFlowScheme::ExplicitRusanov;
		iga::OneDFlowRuntime explicit_runtime(explicit_configuration, explicit_flow, network,
			iga::ResolveOneDInlet(explicit_configuration), directory);
		explicit_runtime.InitializeOpenLoop(1.0e-9);
		const auto committed = explicit_runtime.FlowState();
		explicit_runtime.BeginStep(0.0, 0.01);
		explicit_runtime.SetOpenLoopInlet(explicit_runtime.OpenLoopInlet(0.01, 1.0e-9));
		explicit_runtime.SolveTrial();
		const auto trial = explicit_runtime.FlowState();
		const auto diagnostics = explicit_runtime.Diagnostics();
		assert(diagnostics.planned_configured_substeps == 4 && diagnostics.attempted_configured_substeps == 4
			&& diagnostics.completed_configured_substeps == 4 && diagnostics.explicit_cfl_substep_delta >= 4);
		assert(trial.internal_substeps >= committed.internal_substeps+4);
		explicit_runtime.RollbackTrial();
		assert(explicit_runtime.FlowState().area == committed.area && explicit_runtime.FlowState().flow == committed.flow
			&& explicit_runtime.FlowState().pressure == committed.pressure
			&& explicit_runtime.FlowState().internal_substeps == committed.internal_substeps);
		explicit_runtime.SolveTrial();
		assert(explicit_runtime.FlowState().area == trial.area && explicit_runtime.FlowState().flow == trial.flow
			&& explicit_runtime.FlowState().pressure == trial.pressure
			&& explicit_runtime.Diagnostics().explicit_cfl_substep_delta == diagnostics.explicit_cfl_substep_delta);
		explicit_runtime.CommitStep();
		assert(explicit_runtime.FlowState().completed_step == 4 && explicit_runtime.FlowState().internal_substeps == trial.internal_substeps);
		auto failing_explicit_configuration = explicit_configuration;
		iga::TemporalFunctionDefinition jump;
		jump.name = "late_extreme_jump";
		jump.kind = iga::TemporalFunctionKind::Sinusoid;
		jump.units = "m3/s";
		jump.mean = 1.0e5;
		jump.amplitude = -1.0e5;
		jump.period = 0.01;
		jump.phase = 0.0;
		failing_explicit_configuration.temporal_functions.push_back(jump);
		failing_explicit_configuration.boundaries.front().conditions.front().waveform = jump.name;
		iga::OneDFlowRuntime failing_explicit(failing_explicit_configuration, explicit_flow, network,
			iga::ResolveOneDInlet(failing_explicit_configuration), directory);
		failing_explicit.InitializeOpenLoop(1.0e-9);
		const auto explicit_committed = failing_explicit.FlowState();
		failing_explicit.BeginStep(0.0, 0.01);
		failing_explicit.SetConfiguredOpenLoopInlet();
		RequireRejected([&failing_explicit] { failing_explicit.SolveTrial(); });
		const auto partial = failing_explicit.Diagnostics();
		assert(failing_explicit.CurrentPhase() == iga::OneDFlowRuntime::Phase::TrialSolved);
		assert(partial.attempted_configured_substeps == 2 && partial.completed_configured_substeps == 1
			&& partial.explicit_cfl_substep_delta > 0
			&& failing_explicit.FlowState().internal_substeps-explicit_committed.internal_substeps
			== partial.explicit_cfl_substep_delta);
		const auto frozen_schedule = partial.configured_open_loop_endpoint_flows_m3_s;
		failing_explicit.RollbackTrial();
		assert(failing_explicit.FlowState().area == explicit_committed.area
			&& failing_explicit.FlowState().flow == explicit_committed.flow
			&& failing_explicit.FlowState().internal_substeps == explicit_committed.internal_substeps);
		RequireRejected([&failing_explicit] { failing_explicit.SolveTrial(); });
		assert(failing_explicit.Diagnostics().configured_open_loop_endpoint_flows_m3_s == frozen_schedule
			&& failing_explicit.Diagnostics().attempted_configured_substeps == partial.attempted_configured_substeps
			&& failing_explicit.Diagnostics().completed_configured_substeps == partial.completed_configured_substeps);
		failing_explicit.RollbackTrial();
		failing_explicit.SetOpenLoopInlet(failing_explicit.OpenLoopInlet(0.01, 1.0e-9));
		failing_explicit.SolveTrial();
		assert(failing_explicit.Diagnostics().completed_configured_substeps == 4);
	}
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
	assert(runtime.Transports().front().species.front().boundary_flux_valid);
	assert(Close(first_root.outward_species_flux.at("signal"),
		-runtime.Transports().front().species.front().root_native_flux));
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
	runtime.PrepareCommitStep();
	assert(runtime.CurrentPhase() == iga::OneDFlowRuntime::Phase::CommitPrepared);
	runtime.FinalizeCommitStep();
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
	failing.AbortStep();
	assert(failing.FlowState().completed_step == 0);
	assert(failing.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready);
	assert(failing.Diagnostics().attempted_configured_substeps == 0);
	failing.BeginStep(0.0, configuration.time.dt);
	failing.AbortStep();
	assert(failing.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready);
	failing.AbortStep();
	{
		iga::OneDFlowRuntime prepared_abort(configuration, flow, network,
			iga::ResolveOneDInlet(configuration), directory);
		prepared_abort.InitializeOpenLoop(1.0e-9);
		const auto committed = prepared_abort.FlowState();
		prepared_abort.BeginStep(0.0, configuration.time.dt);
		prepared_abort.SetOpenLoopInlet(
			prepared_abort.OpenLoopInlet(configuration.time.dt, 2.0e-9));
		prepared_abort.SolveTrial();
		prepared_abort.AbortStep();
		assert(prepared_abort.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready);
		assert(prepared_abort.FlowState().area == committed.area);
		assert(prepared_abort.Diagnostics().attempted_configured_substeps == 0);
		prepared_abort.BeginStep(0.0, configuration.time.dt);
		prepared_abort.SetOpenLoopInlet(
			prepared_abort.OpenLoopInlet(configuration.time.dt, 2.0e-9));
		prepared_abort.SolveTrial();
		prepared_abort.PrepareCommitStep();
		prepared_abort.AbortStep();
		assert(prepared_abort.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready);
		assert(prepared_abort.FlowState().area == committed.area);
		assert(prepared_abort.FlowState().flow == committed.flow);
		assert(prepared_abort.FlowState().node_pressure == committed.node_pressure);
		assert(prepared_abort.Diagnostics().attempted_configured_substeps == 0);
	}

	auto sub_configuration = configuration;
	sub_configuration.time.dt = 0.0025;
	sub_configuration.time.steps = 4;
	std::vector<double> implicit_dts;
	int implicit_calls = 0;
	bool fail_implicit_substep = true;
	iga::OneDFlowRuntime implicit_subcycled(sub_configuration, failing_flow, network,
		iga::ResolveOneDInlet(sub_configuration), directory,
		[&implicit_dts, &implicit_calls, &fail_implicit_substep](const iga::OneDNetwork&, const iga::OneDFlowSystemDefinition&,
			iga::OneDFlowState& state, double inlet_flow, double dt_s) {
			implicit_dts.push_back(dt_s);
			++implicit_calls;
			if (fail_implicit_substep && implicit_calls == 3) throw std::runtime_error("substep three failure");
			state.inlet_flow = inlet_flow;
		});
	implicit_subcycled.InitializeOpenLoop(1.0e-9);
	implicit_subcycled.BeginStep(0.0, 0.01);
	implicit_subcycled.SetOpenLoopInlet(implicit_subcycled.OpenLoopInlet(0.01, 1.0e-9));
	RequireRejected([&implicit_subcycled] { implicit_subcycled.SolveTrial(); });
	assert(implicit_subcycled.Diagnostics().planned_configured_substeps == 4);
	assert(implicit_subcycled.Diagnostics().attempted_configured_substeps == 3);
	assert(implicit_subcycled.Diagnostics().completed_configured_substeps == 2);
	assert(implicit_subcycled.Diagnostics().explicit_cfl_substep_delta == 0);
	implicit_subcycled.RollbackTrial();
	assert(implicit_subcycled.FlowState().completed_step == 0);
	implicit_calls = 0;
	fail_implicit_substep = false;
	implicit_dts.clear();
	implicit_subcycled.SolveTrial();
	assert(implicit_dts.size() == 4);
	for (const double dt_s : implicit_dts) assert(Close(dt_s, 0.0025));
	assert(implicit_subcycled.FlowState().completed_step == 4);

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

	auto rcr_configuration = configuration;
	auto& rcr = rcr_configuration.boundaries.back().conditions.front();
	rcr.type = "windkessel_rcr";
	rcr.proximal_resistance = 1.0e8;
	rcr.distal_resistance = 1.0e9;
	rcr.capacitance = 1.0e-10;
	rcr.reference_pressure = 0.0;
	rcr.initial_pressure = 3.0;
	iga::OneDFlowRuntime rcr_runtime(rcr_configuration, flow, network,
		iga::ResolveOneDInlet(rcr_configuration), directory);
	rcr_runtime.InitializeOpenLoop(1.0e-9);
	const double committed_capacitor = rcr_runtime.FlowState().outlets.front().capacitor_pressure;
	rcr_runtime.BeginStep(0.0, configuration.time.dt);
	rcr_runtime.SetOpenLoopInlet(rcr_runtime.OpenLoopInlet(configuration.time.dt, 1.0e-9));
	RequireRejected([&rcr_runtime, &forbidden_pressure] {
		rcr_runtime.SetPortInput("outlet:2", forbidden_pressure);
	});
	rcr_runtime.SolveTrial();
	const double trial_capacitor = rcr_runtime.FlowState().outlets.front().capacitor_pressure;
	assert(!Close(trial_capacitor, committed_capacitor));
	rcr_runtime.RollbackTrial();
	assert(rcr_runtime.FlowState().outlets.front().capacitor_pressure == committed_capacitor);
	rcr_runtime.SolveTrial();
	assert(rcr_runtime.FlowState().outlets.front().capacitor_pressure == trial_capacitor);
	rcr_runtime.CommitStep();
	assert(rcr_runtime.FlowState().outlets.front().capacitor_pressure == trial_capacitor);
	RequireRejected([&rcr_runtime] { rcr_runtime.CommitStep(); });
	{
		auto rcr_sub_configuration = rcr_configuration;
		rcr_sub_configuration.time.dt = 0.0025;
		rcr_sub_configuration.time.steps = 4;
		iga::OneDFlowRuntime macro(rcr_sub_configuration, flow, network,
			iga::ResolveOneDInlet(rcr_sub_configuration), directory);
		iga::OneDFlowRuntime sequential(rcr_sub_configuration, flow, network,
			iga::ResolveOneDInlet(rcr_sub_configuration), directory);
		macro.InitializeOpenLoop(1.0e-9);
		sequential.InitializeOpenLoop(1.0e-9);
		macro.BeginStep(0.0, 0.01);
		macro.SetOpenLoopInlet(macro.OpenLoopInlet(0.01, 1.0e-9));
		macro.SolveTrial();
		for (int step = 0; step < 4; ++step) {
			sequential.BeginStep(sequential.FlowState().physical_time, 0.0025);
			sequential.SetOpenLoopInlet(sequential.OpenLoopInlet((step+1)*0.0025, 1.0e-9));
			sequential.SolveTrial();
			sequential.CommitStep();
		}
		assert(Close(macro.FlowState().outlets.front().capacitor_pressure,
			sequential.FlowState().outlets.front().capacitor_pressure));
		assert(Close(macro.FlowState().outlets.front().pressure, sequential.FlowState().outlets.front().pressure));
		assert(macro.FlowState().completed_step == 4 && Close(macro.FlowState().physical_time, 0.01));
		macro.RollbackTrial();
		macro.SolveTrial();
		assert(Close(macro.FlowState().outlets.front().capacitor_pressure,
			sequential.FlowState().outlets.front().capacitor_pressure));
	}
	{
		iga::CouplingPort terminal;
		terminal.id = "terminal";
		terminal.subsystem_id = "configured";
		terminal.locator_kind = "runtime_port";
		terminal.locator = "outlet:2";
		terminal.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
			iga::PortQuantity::MeanPressure};
		terminal.requires = {iga::PortQuantity::MeanPressure};
		iga::OneDFlowRuntime native(configuration, flow, network,
			iga::ResolveOneDInlet(configuration), directory);
		native.InitializeOpenLoop(1.0e-9);
		iga::OneDFlowDomainAdapter adapter("configured", native, {terminal},
			iga::OneDInletPolicy::ConfiguredOpenLoop);
		adapter.BeginStep({0, 0.0, configuration.time.dt});
		iga::PortBoundaryData pressure;
		pressure.time_s = configuration.time.dt;
		pressure.mean_pressure_pa = 5.0;
		adapter.SetPortInput("terminal", pressure);
		adapter.SolveTrial();
		assert(adapter.GetPortState("terminal").mean_pressure_pa.has_value());
		adapter.PrepareCommitStep();
		adapter.FinalizeCommitStep();
		assert(native.FlowState().completed_step == 1);

		auto root = terminal;
		root.id = "inlet";
		root.subsystem_id = "coupled";
		root.locator = "root";
		root.requires = {iga::PortQuantity::FlowRate};
		iga::OneDFlowRuntime coupled_native(configuration, flow, network,
			iga::ResolveOneDInlet(configuration), directory);
		coupled_native.InitializeOpenLoop(1.0e-9);
		iga::OneDFlowDomainAdapter coupled("coupled", coupled_native, {root},
			iga::OneDInletPolicy::CoupledRoot);
		coupled.BeginStep({0, 0.0, configuration.time.dt});
		iga::PortBoundaryData root_flow;
		root_flow.time_s = configuration.time.dt;
		root_flow.outward_flow_m3_s = -2.0e-9;
		coupled.SetPortInput("inlet", root_flow);
		coupled.SolveTrial();
		assert(Close(*coupled.GetPortState("inlet").outward_flow_m3_s, -2.0e-9));
		coupled.AbortStep();
		assert(coupled_native.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready);
		RequireRejected([&native, &root] {
			iga::OneDFlowDomainAdapter invalid("configured", native, {root},
				iga::OneDInletPolicy::ConfiguredOpenLoop);
		});
		auto reversed = terminal;
		reversed.orientation.native_to_outward_sign = -1;
		RequireRejected([&native, &reversed] {
			iga::OneDFlowDomainAdapter invalid("configured", native, {reversed},
				iga::OneDInletPolicy::ConfiguredOpenLoop);
		});
	}
	fs::remove_all(directory);
	std::cout << "one-dimensional runtime tests passed\n";
}

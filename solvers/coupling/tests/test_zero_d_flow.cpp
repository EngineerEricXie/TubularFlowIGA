#include "ZeroDFlowDomain.hpp"
#include "DomainRuntimeRegistry.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace {

void RequireRejected(const std::function<void()>& operation)
{
	bool rejected = false;
	try { operation(); }
	catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

bool Near(double first, double second, double tolerance = 1.0e-12)
{
	return std::abs(first-second) <= tolerance*std::max(1.0, std::max(std::abs(first), std::abs(second)));
}

iga::ZeroDFlowModel Source()
{
	iga::ZeroDFlowModel model;
	model.role = iga::ZeroDFlowRole::SourceReservoir;
	model.source = {2.0, 4.0, 3.0};
	return model;
}

iga::ZeroDFlowModel Terminal()
{
	iga::ZeroDFlowModel model;
	model.role = iga::ZeroDFlowRole::TerminalRcr;
	model.terminal = {1.5, 5.0, 2.0, 7.0};
	return model;
}

iga::PortBoundaryData PressureInput(double time_s, double pressure_pa)
{
	iga::PortBoundaryData input;
	input.time_s = time_s;
	input.mean_pressure_pa = pressure_pa;
	return input;
}

iga::PortBoundaryData FlowInput(double time_s, double outward_flow_m3_s)
{
	iga::PortBoundaryData input;
	input.time_s = time_s;
	input.outward_flow_m3_s = outward_flow_m3_s;
	return input;
}

} // namespace

int main()
{
	{
		const auto model = Source();
		const iga::ZeroDFlowState state{11.0};
		const double gamma = 5.0;
		const double dt = 0.25;
		const auto trial = iga::EvaluateZeroDFlowTrial(model, state, gamma, dt);
		const double expected_pressure = ((2.0/dt)*11.0+3.0+gamma/4.0)/(2.0/dt+1.0/4.0);
		assert(Near(trial.state.stored_pressure_pa, expected_pressure));
		assert(Near(*trial.port.mean_pressure_pa, gamma));
		assert(Near(*trial.port.outward_flow_m3_s, (expected_pressure-gamma)/4.0));
		assert(Near(2.0*(expected_pressure-11.0)/dt,
			3.0-*trial.port.outward_flow_m3_s));
		assert(std::abs(trial.storage.residual_m3) < 1.0e-14);
	}
	{
		const auto model = Terminal();
		const iga::ZeroDFlowState state{9.0};
		const double outward_port_flow = -2.0; // Qin=2 m^3/s into the terminal model.
		const double dt = 0.1;
		const auto trial = iga::EvaluateZeroDFlowTrial(model, state, outward_port_flow, dt);
		const double inlet = -outward_port_flow;
		const double expected = ((2.0/dt)*9.0+inlet+7.0/5.0)/(2.0/dt+1.0/5.0);
		assert(Near(trial.state.stored_pressure_pa, expected));
		assert(Near(*trial.port.mean_pressure_pa, expected+1.5*inlet));
		assert(Near(2.0*(expected-9.0)/dt, inlet-(expected-7.0)/5.0));
		assert(std::abs(trial.storage.residual_m3) < 1.0e-14);
		const auto reverse = iga::EvaluateZeroDFlowTrial(model, state, 0.75, dt);
		assert(*reverse.port.outward_flow_m3_s > 0.0); // Reverse Qin is intentionally not clamped.
		const double reverse_inlet = -0.75;
		const double expected_reverse = ((2.0/dt)*9.0+reverse_inlet+7.0/5.0)
			/(2.0/dt+1.0/5.0);
		assert(Near(reverse.state.stored_pressure_pa, expected_reverse));
		assert(Near(*reverse.port.mean_pressure_pa, expected_reverse+1.5*reverse_inlet));
		assert(Near(2.0*(expected_reverse-9.0)/dt,
			reverse_inlet-(expected_reverse-7.0)/5.0));
		assert(std::abs(reverse.storage.residual_m3) < 1.0e-14);
	}
	{
		// Exact equilibrium for source: Ps=pGamma+Rs*Qp.
		const auto model = Source();
		const iga::ZeroDFlowState equilibrium{17.0};
		const auto trial = iga::EvaluateZeroDFlowTrial(model, equilibrium, 5.0, 0.4);
		assert(Near(trial.state.stored_pressure_pa, equilibrium.stored_pressure_pa));
		assert(Near(*trial.port.outward_flow_m3_s, 3.0));
	}
	{
		const auto source_port = iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir);
		const auto terminal_port = iga::MakeZeroDFlowPort("terminal", iga::ZeroDFlowRole::TerminalRcr);
		assert(!source_port.provides.count(iga::PortQuantity::Area));
		iga::ValidateZeroDFlowDomainMetadata("source", {source_port}, iga::ZeroDFlowRole::SourceReservoir);
		iga::ValidateZeroDFlowDomainMetadata("terminal", {terminal_port}, iga::ZeroDFlowRole::TerminalRcr);
	}
	{
		const auto model = Source();
		const iga::ZeroDFlowState state{1.0};
		const auto source_identity = iga::BuildZeroDFlowModelIdentitySha256(model);
		auto changed_source = model;
		changed_source.source.capacitance_m3_pa += 1.0;
		assert(source_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_source));
		changed_source = model;
		changed_source.source.resistance_pa_s_m3 += 1.0;
		assert(source_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_source));
		changed_source = model;
		changed_source.source.prescribed_flow_m3_s += 1.0;
		assert(source_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_source));
		const auto terminal = Terminal();
		const auto terminal_identity = iga::BuildZeroDFlowModelIdentitySha256(terminal);
		auto changed_terminal = terminal;
		changed_terminal.terminal.proximal_resistance_pa_s_m3 += 1.0;
		assert(terminal_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_terminal));
		changed_terminal = terminal;
		changed_terminal.terminal.distal_resistance_pa_s_m3 += 1.0;
		assert(terminal_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_terminal));
		changed_terminal = terminal;
		changed_terminal.terminal.capacitance_m3_pa += 1.0;
		assert(terminal_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_terminal));
		changed_terminal = terminal;
		changed_terminal.terminal.distal_pressure_pa += 1.0;
		assert(terminal_identity != iga::BuildZeroDFlowModelIdentitySha256(changed_terminal));
		assert(iga::BuildZeroDFlowStateIdentitySha256(model, state)
			!= iga::BuildZeroDFlowStateIdentitySha256(model, iga::ZeroDFlowState{2.0}));
		const iga::ZeroDFlowStepAccounting accounting{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
		const auto accounting_identity = iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			accounting);
		auto changed_accounting = accounting;
		changed_accounting.initial_stored_volume_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		changed_accounting = accounting;
		changed_accounting.final_stored_volume_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		changed_accounting = accounting;
		changed_accounting.prescribed_source_amount_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		changed_accounting = accounting;
		changed_accounting.distal_sink_amount_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		changed_accounting = accounting;
		changed_accounting.outward_graph_port_amount_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		changed_accounting = accounting;
		changed_accounting.residual_m3 += 1.0;
		assert(accounting_identity != iga::BuildZeroDFlowStepAccountingIdentitySha256(model,
			changed_accounting));
		RequireRejected([] { auto bad = Source(); bad.source.capacitance_m3_pa = 0.0;
			iga::ValidateZeroDFlowModel(bad); });
		RequireRejected([] { auto bad = Terminal(); bad.terminal.distal_resistance_pa_s_m3
			= std::numeric_limits<double>::infinity(); iga::ValidateZeroDFlowModel(bad); });
		RequireRejected([&] { (void)iga::EvaluateZeroDFlowTrial(model, state,
			0.0, 0.0); });
		RequireRejected([&] { (void)iga::EvaluateZeroDFlowTrial(model, state,
			std::numeric_limits<double>::quiet_NaN(), 0.1); });
	}
	{
		// Analytic source RC reference: C dP/dt=Qp-P/R, P(0)=P0, pGamma=0.
		// P(t)=R Qp+(P0-R Qp) exp(-t/(R C)); backward Euler is first order.
		const auto model = Source();
		const iga::ZeroDFlowState initial{1.0};
		const double final_time = 1.0;
		auto error = [&](int steps) {
			auto state = initial;
			const double dt = final_time/steps;
			for (int i = 0; i < steps; ++i)
				state = iga::EvaluateZeroDFlowTrial(model, state, 0.0, dt).state;
			const double exact = 12.0+(1.0-12.0)*std::exp(-final_time/8.0);
			return std::abs(state.stored_pressure_pa-exact);
		};
		const double coarse = error(40);
		const double fine = error(80);
		assert(coarse/fine > 1.9 && coarse/fine < 2.1);
	}
	{
		// The production wrapper keeps its committed pressure and publication
		// separate from the trial produced at t_{n+1}.
		iga::ZeroDFlowDomainRuntime runtime("source", Source(), {11.0},
			{iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir)});
		const auto before = runtime.CommittedStateIdentitySha256();
		RequireRejected([&] { runtime.SolveTrial(); });
		RequireRejected([&] { runtime.SetPortInput("port", PressureInput(0.25, 5.0)); });
		runtime.BeginStep({0, 0.0, 0.25});
		RequireRejected([&] { runtime.GetPortState("port"); });
		RequireRejected([&] { runtime.PrepareCommitStep(); });
		RequireRejected([&] { runtime.SetPortInput("wrong", PressureInput(0.25, 5.0)); });
		RequireRejected([&] { runtime.SetPortInput("port", FlowInput(0.25, -1.0)); });
		RequireRejected([&] { runtime.SetPortInput("port", PressureInput(0.2, 5.0)); });
		RequireRejected([&] { runtime.SetPortInput("port", PressureInput(0.25,
			std::numeric_limits<double>::quiet_NaN())); });
		runtime.SetPortInput("port", PressureInput(0.25, 5.0));
		RequireRejected([&] { runtime.SetPortInput("port", PressureInput(0.25, 5.0)); });
		runtime.SolveTrial();
		const auto expected = iga::EvaluateZeroDFlowTrial(Source(), {11.0}, 5.0, 0.25, 0.0);
		assert(runtime.CommittedStateIdentitySha256() == before);
		assert(!runtime.CommittedPortState());
		const auto state = runtime.GetPortState("port");
		assert(state.time_s == 0.25);
		assert(Near(*state.outward_flow_m3_s, *expected.port.outward_flow_m3_s));
		const auto accounting = runtime.TrialStepAccounting();
		assert(Near(accounting.initial_stored_volume_m3, 22.0));
		assert(Near(accounting.final_stored_volume_m3,
			2.0*expected.state.stored_pressure_pa));
		assert(std::abs(accounting.residual_m3) < 1.0e-14);
		const auto accounting_identity = iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(),
			accounting);
		runtime.PrepareCommitStep();
		assert(!runtime.CommittedPortState()); // Preparation does not publish.
		runtime.FinalizeCommitStep();
		runtime.FinalizeCommitStep(); // Idempotent, single promotion.
		RequireRejected([&] { (void)runtime.TrialStepAccounting(); });
		assert(iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(), accounting)
			== accounting_identity);
		assert(runtime.CommittedStepCount() == 1);
		assert(runtime.CommittedStepIndex() == 0);
		assert(runtime.CommittedTime() == 0.25);
		assert(Near(runtime.CommittedState().stored_pressure_pa,
			expected.state.stored_pressure_pa));
		assert(runtime.CommittedPortState());
		assert(runtime.CommittedPortState()->time_s == 0.25);
		assert(runtime.CommittedStepAccounting());
		RequireRejected([&] { runtime.BeginStep({2, 0.25, 0.25}); });
		RequireRejected([&] { runtime.BeginStep({1, 0.0, 0.25}); });
		iga::ZeroDFlowDomainRuntime rerun("source", Source(), {11.0},
			{iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir)});
		rerun.BeginStep({0, 0.0, 0.25});
		rerun.SetPortInput("port", PressureInput(0.25, 5.0));
		rerun.SolveTrial(); rerun.PrepareCommitStep(); rerun.FinalizeCommitStep();
		assert(rerun.CommittedStateIdentitySha256() == runtime.CommittedStateIdentitySha256());
		assert(*rerun.CommittedPortState()->outward_flow_m3_s
			== *runtime.CommittedPortState()->outward_flow_m3_s);
	}
	{
		// A changed strong-iteration input is replayed from the committed t_n
		// state, and neither rollback nor abort can leak a trial result.
		iga::ZeroDFlowDomainRuntime runtime("source", Source(), {11.0},
			{iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir)});
		runtime.BeginStep({0, 0.0, 0.25});
		runtime.SetPortInput("port", PressureInput(0.25, 2.0));
		runtime.SolveTrial();
		const auto first = runtime.GetPortState("port");
		const auto rolled_back_accounting = runtime.TrialStepAccounting();
		const auto rolled_back_accounting_identity =
			iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(), rolled_back_accounting);
		const auto committed = runtime.CommittedStateIdentitySha256();
		runtime.RollbackTrial();
		assert(runtime.CommittedStateIdentitySha256() == committed);
		RequireRejected([&] { (void)runtime.TrialStepAccounting(); });
		assert(iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(), rolled_back_accounting)
			== rolled_back_accounting_identity);
		runtime.SetPortInput("port", PressureInput(0.25, 5.0));
		runtime.SolveTrial();
		const auto replay = runtime.GetPortState("port");
		const auto expected = iga::EvaluateZeroDFlowTrial(Source(), {11.0}, 5.0, 0.25);
		assert(!Near(*first.outward_flow_m3_s, *replay.outward_flow_m3_s));
		assert(Near(*replay.outward_flow_m3_s, *expected.port.outward_flow_m3_s));
		const auto aborted_accounting = runtime.TrialStepAccounting();
		const auto aborted_accounting_identity =
			iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(), aborted_accounting);
		runtime.PrepareCommitStep();
		runtime.AbortStep();
		assert(runtime.CommittedStateIdentitySha256() == committed);
		assert(runtime.CommittedStepCount() == 0);
		RequireRejected([&] { runtime.GetPortState("port"); });
		RequireRejected([&] { (void)runtime.TrialStepAccounting(); });
		assert(iga::BuildZeroDFlowStepAccountingIdentitySha256(Source(), aborted_accounting)
			== aborted_accounting_identity);
		RequireRejected([&] { runtime.AbortStep(); });
	}
	{
		iga::ZeroDFlowDomainRuntime runtime("terminal", Terminal(), {9.0},
			{iga::MakeZeroDFlowPort("terminal", iga::ZeroDFlowRole::TerminalRcr)}, 1.5);
		runtime.BeginStep({0, 1.5, 0.1});
		RequireRejected([&] { runtime.SetPortInput("port", PressureInput(1.6, 3.0)); });
		runtime.SetPortInput("port", FlowInput(1.6, -2.0));
		runtime.SolveTrial();
		const auto expected = iga::EvaluateZeroDFlowTrial(Terminal(), {9.0}, -2.0, 0.1, 1.5);
		const auto state = runtime.GetPortState("port");
		assert(state.time_s == 1.6);
		assert(Near(*state.mean_pressure_pa, *expected.port.mean_pressure_pa));
		assert(Near(*state.outward_flow_m3_s, -2.0));
		assert(std::abs(runtime.TrialStepAccounting().residual_m3) < 1.0e-14);
		runtime.PrepareCommitStep();
		runtime.FinalizeCommitStep();
		assert(runtime.CommittedTime() == 1.6 && runtime.CommittedStepIndex() == 0);
	}
	{
		// Metadata/model-role binding is strict, and a normal registry sees the
		// runtime as a ZeroDFlow node without any runner integration.
		RequireRejected([] {
			iga::ZeroDFlowDomainRuntime("source", Source(), {1.0},
				{iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::TerminalRcr)});
		});
		auto port = iga::MakeZeroDFlowPort("source", iga::ZeroDFlowRole::SourceReservoir);
		const iga::SimulationGraph graph({{"source", iga::DomainKind::ZeroDFlow, {port}}}, {});
		assert((iga::MakeSequentialPlan(graph, "source").domain_ids
			== std::vector<std::string>{"source"}));
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.emplace_back(std::make_unique<iga::ZeroDFlowDomainRuntime>(
			"source", Source(), iga::ZeroDFlowState{1.0}, std::vector<iga::CouplingPort>{port}));
		iga::DomainRuntimeRegistry registry(graph, std::move(runtimes));
		assert(registry.Runtime("source").Kind() == iga::DomainKind::ZeroDFlow);
		const iga::SimulationGraph wrong_kind({{"source", iga::DomainKind::OneDFlow, {port}}}, {});
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> mismatched;
		mismatched.emplace_back(std::make_unique<iga::ZeroDFlowDomainRuntime>(
			"source", Source(), iga::ZeroDFlowState{1.0}, std::vector<iga::CouplingPort>{port}));
		RequireRejected([&] { iga::DomainRuntimeRegistry(wrong_kind, std::move(mismatched)); });
	}
	return 0;
}

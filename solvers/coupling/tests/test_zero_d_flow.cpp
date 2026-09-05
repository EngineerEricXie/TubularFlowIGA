#include "ZeroDFlowDomain.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

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
		assert(iga::BuildZeroDFlowModelIdentitySha256(model)
			== iga::BuildZeroDFlowModelIdentitySha256(model));
		assert(iga::BuildZeroDFlowStateIdentitySha256(model, state)
			== iga::BuildZeroDFlowStateIdentitySha256(model, state));
		auto changed = model;
		changed.source.prescribed_flow_m3_s += 1.0;
		assert(iga::BuildZeroDFlowModelIdentitySha256(model)
			!= iga::BuildZeroDFlowModelIdentitySha256(changed));
		assert(iga::BuildZeroDFlowStateIdentitySha256(model, state)
			!= iga::BuildZeroDFlowStateIdentitySha256(model, iga::ZeroDFlowState{2.0}));
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
	return 0;
}

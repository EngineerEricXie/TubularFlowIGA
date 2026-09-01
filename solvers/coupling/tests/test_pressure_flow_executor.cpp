#include "PressureFlowComponentExecutor.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void RequireRejected(const std::function<void()>& operation)
{
	bool rejected = false;
	try { operation(); }
	catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

iga::CouplingPort Port(const std::string& domain, const std::string& id,
	std::set<iga::PortQuantity> requires)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = domain;
	port.locator_kind = "fake";
	port.locator = id;
	port.provides = {iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure};
	port.requires = std::move(requires);
	iga::ValidateCouplingPort(port);
	return port;
}

iga::SimulationGraph Chain()
{
	auto upstream = iga::DomainNode{"up", iga::DomainKind::OneDFlow,
		{Port("up", "terminal", {iga::PortQuantity::MeanPressure})}};
	auto middle = iga::DomainNode{"mid", iga::DomainKind::ThreeDBodyFittedFlow,
		{Port("mid", "inlet", {iga::PortQuantity::FlowRate}),
		 Port("mid", "outlet", {iga::PortQuantity::MeanPressure})}};
	auto downstream = iga::DomainNode{"down", iga::DomainKind::OneDFlow,
		{Port("down", "root", {iga::PortQuantity::FlowRate})}};
	return iga::SimulationGraph({std::move(upstream), std::move(middle), std::move(downstream)}, {
		{"left", {"up", "terminal"}, {"mid", "inlet"}, iga::CouplingLaw::PressureFlow},
		{"right", {"mid", "outlet"}, {"down", "root"}, iga::CouplingLaw::PressureFlow}});
}

iga::SimulationGraph Branch(bool permuted = false)
{
	auto source = iga::DomainNode{"up", iga::DomainKind::OneDFlow,
		{Port("up", "terminal", {iga::PortQuantity::MeanPressure})}};
	auto junction = iga::DomainNode{"junction", iga::DomainKind::ThreeDBodyFittedFlow,
		{Port("junction", "inlet", {iga::PortQuantity::FlowRate}),
		 Port("junction", "left", {iga::PortQuantity::MeanPressure}),
		 Port("junction", "right", {iga::PortQuantity::MeanPressure})}};
	auto left = iga::DomainNode{"branch_a", iga::DomainKind::OneDFlow,
		{Port("branch_a", "root", {iga::PortQuantity::FlowRate})}};
	auto right = iga::DomainNode{"branch_b", iga::DomainKind::OneDFlow,
		{Port("branch_b", "root", {iga::PortQuantity::FlowRate})}};
	std::vector<iga::CouplingEdge> edges{
		{"c_right", {"junction", "right"}, {"branch_b", "root"},
			iga::CouplingLaw::PressureFlow},
		{"a_source", {"up", "terminal"}, {"junction", "inlet"},
			iga::CouplingLaw::PressureFlow},
		{"b_left", {"junction", "left"}, {"branch_a", "root"},
			iga::CouplingLaw::PressureFlow}};
	if (permuted) {
		std::reverse(edges.begin(), edges.end());
		for (auto& edge : edges) std::swap(edge.first, edge.second);
	}
	return iga::SimulationGraph({std::move(right), std::move(source), std::move(left),
		std::move(junction)}, std::move(edges));
}

class FakeRuntime : public iga::CoupledDomainRuntime {
public:
	FakeRuntime(std::string id, iga::DomainKind kind, std::vector<iga::CouplingPort> ports)
		: id_(std::move(id)), kind_(kind), ports_(std::move(ports)) {}

	const std::string& DomainId() const noexcept override { return id_; }
	iga::DomainKind Kind() const noexcept override { return kind_; }
	const std::vector<iga::CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const iga::DomainStepContext& step) override
	{
		step.Validate();
		if (phase_ != Phase::Ready) throw std::runtime_error("fake BeginStep phase");
		end_time_ = step.EndTime();
		inputs_.clear();
		phase_ = Phase::Open;
		++begins;
	}

	void SetPortInput(const std::string& port_id,
		const iga::PortBoundaryData& input) override
	{
		if (phase_ != Phase::Open) throw std::runtime_error("fake SetPortInput phase");
		iga::ValidatePortBoundaryData(input);
		inputs_[port_id] = input;
	}

	void SolveTrial() override
	{
		if (phase_ != Phase::Open) throw std::runtime_error("fake SolveTrial phase");
		if (fail_solve) {
			phase_ = Phase::Solved;
			throw std::runtime_error("injected fake solve failure");
		}
		states_.clear();
		if (id_ == "up") {
			RequirePressure("terminal");
			states_["terminal"] = State(nonfinite_state
				? std::numeric_limits<double>::quiet_NaN() : 2.0,
				*inputs_.at("terminal").mean_pressure_pa);
		} else if (id_ == "mid") {
			RequireFlow("inlet");
			RequirePressure("outlet");
			states_["inlet"] = State(*inputs_.at("inlet").outward_flow_m3_s, 10.0);
			states_["outlet"] = State(3.0, *inputs_.at("outlet").mean_pressure_pa);
		} else if (id_ == "down") {
			RequireFlow("root");
			states_["root"] = State(*inputs_.at("root").outward_flow_m3_s, 20.0);
		} else if (id_ == "junction") {
			RequireFlow("inlet");
			RequirePressure("left");
			RequirePressure("right");
			states_["inlet"] = State(*inputs_.at("inlet").outward_flow_m3_s, 10.0);
			states_["left"] = State(1.0, *inputs_.at("left").mean_pressure_pa);
			states_["right"] = State(2.0, *inputs_.at("right").mean_pressure_pa);
		} else if (id_ == "branch_a" || id_ == "branch_b") {
			RequireFlow("root");
			states_["root"] = State(*inputs_.at("root").outward_flow_m3_s,
				id_ == "branch_a" ? 20.0 : 30.0);
		} else {
			throw std::runtime_error("unknown fake runtime domain");
		}
		phase_ = Phase::Solved;
		++solves;
	}

	iga::PortState GetPortState(const std::string& port_id) const override
	{
		if (phase_ != Phase::Solved) throw std::runtime_error("fake GetPortState phase");
		return states_.at(port_id);
	}

	void RollbackTrial() override
	{
		if (phase_ != Phase::Solved) throw std::runtime_error("fake RollbackTrial phase");
		inputs_.clear();
		states_.clear();
		phase_ = Phase::Open;
		++rollbacks;
	}

	void AbortStep() override
	{
		if (phase_ == Phase::Ready) return;
		inputs_.clear();
		states_.clear();
		phase_ = Phase::Ready;
		++aborts;
		if (fail_abort) throw std::runtime_error("injected fake abort failure");
	}

	void PrepareCommitStep() override
	{
		if (phase_ != Phase::Solved) throw std::runtime_error("fake prepare phase");
		if (fail_prepare) throw std::runtime_error("injected fake prepare failure");
		phase_ = Phase::Prepared;
	}

	void FinalizeCommitStep() noexcept override
	{
		if (phase_ != Phase::Prepared) std::terminate();
		phase_ = Phase::Ready;
		++commits;
	}

	bool fail_solve = false;
	bool fail_prepare = false;
	bool fail_abort = false;
	bool nonfinite_state = false;
	int begins = 0;
	int solves = 0;
	int rollbacks = 0;
	int aborts = 0;
	int commits = 0;

private:
	enum class Phase { Ready, Open, Solved, Prepared };

	iga::PortState State(double flow, double pressure) const
	{
		iga::PortState state;
		state.time_s = end_time_;
		state.outward_flow_m3_s = flow;
		state.mean_pressure_pa = pressure;
		return state;
	}

	void RequireFlow(const std::string& port) const
	{
		if (!inputs_.count(port) || !inputs_.at(port).outward_flow_m3_s)
			throw std::runtime_error("fake missing flow input");
	}

	void RequirePressure(const std::string& port) const
	{
		if (!inputs_.count(port) || !inputs_.at(port).mean_pressure_pa)
			throw std::runtime_error("fake missing pressure input");
	}

	std::string id_;
	iga::DomainKind kind_;
	std::vector<iga::CouplingPort> ports_;
	Phase phase_ = Phase::Ready;
	double end_time_ = 0.0;
	std::map<std::string, iga::PortBoundaryData> inputs_;
	std::map<std::string, iga::PortState> states_;
};

struct Fixture {
	explicit Fixture(const iga::SimulationGraph& graph)
	{
		auto upstream = std::make_unique<FakeRuntime>("up", iga::DomainKind::OneDFlow,
			graph.Domain("up").ports);
		auto middle = std::make_unique<FakeRuntime>("mid",
			iga::DomainKind::ThreeDBodyFittedFlow, graph.Domain("mid").ports);
		auto downstream = std::make_unique<FakeRuntime>("down", iga::DomainKind::OneDFlow,
			graph.Domain("down").ports);
		up = upstream.get();
		mid = middle.get();
		down = downstream.get();
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.push_back(std::move(middle));
		runtimes.push_back(std::move(downstream));
		runtimes.push_back(std::move(upstream));
		registry = std::make_unique<iga::DomainRuntimeRegistry>(graph, std::move(runtimes));
	}

	FakeRuntime* up = nullptr;
	FakeRuntime* mid = nullptr;
	FakeRuntime* down = nullptr;
	std::unique_ptr<iga::DomainRuntimeRegistry> registry;
};

struct BranchFixture {
	explicit BranchFixture(const iga::SimulationGraph& graph)
	{
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		for (const auto& id : {"branch_b", "junction", "up", "branch_a"}) {
			auto runtime = std::make_unique<FakeRuntime>(id, graph.Domain(id).kind,
				graph.Domain(id).ports);
			by_id[id] = runtime.get();
			runtimes.push_back(std::move(runtime));
		}
		registry = std::make_unique<iga::DomainRuntimeRegistry>(graph, std::move(runtimes));
	}

	std::map<std::string, FakeRuntime*> by_id;
	std::unique_ptr<iga::DomainRuntimeRegistry> registry;
};

std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> RuntimeSet(
	const iga::SimulationGraph& graph, iga::DomainKind upstream_kind,
	std::vector<iga::CouplingPort> upstream_ports)
{
	std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
	runtimes.push_back(std::make_unique<FakeRuntime>("up", upstream_kind,
		std::move(upstream_ports)));
	runtimes.push_back(std::make_unique<FakeRuntime>("mid",
		iga::DomainKind::ThreeDBodyFittedFlow, graph.Domain("mid").ports));
	runtimes.push_back(std::make_unique<FakeRuntime>("down", iga::DomainKind::OneDFlow,
		graph.Domain("down").ports));
	return runtimes;
}

} // namespace

int main()
{
	const auto graph = Chain();
	const iga::DomainStepContext step{0, 0.0, 0.1};
	{
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step, {{"left", 0.0}, {"right", 0.0}});
		assert(result.iterations.size() == 1);
		assert(result.iterations.front().edges.size() == 2);
		assert(result.accepted_ports.size() == 4);
		assert(fixture.up->commits == 1 && fixture.mid->commits == 1
			&& fixture.down->commits == 1);
		for (const auto& edge : result.iterations.front().edges)
			assert(edge.flow_residual_m3_s == 0.0);
	}
	{
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Fixed;
		controls.maximum_iterations = 40;
		controls.pressure_relative_tolerance = 1.0e-7;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step, {{"left", 0.0}, {"right", 0.0}});
		assert(result.iterations.size() > 1);
		assert(result.iterations.back().converged);
		assert(fixture.up->rollbacks+1 == fixture.up->solves);
		assert(fixture.mid->rollbacks+1 == fixture.mid->solves);
		assert(fixture.down->rollbacks+1 == fixture.down->solves);
	}
	{
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Aitken;
		controls.maximum_iterations = 8;
		controls.pressure_relative_tolerance = 1.0e-12;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step, {{"left", 0.0}, {"right", 0.0}});
		assert(result.iterations.size() == 3);
		assert(result.iterations.back().converged);
	}
	{
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Fixed;
		controls.maximum_iterations = 2;
		controls.pressure_relative_tolerance = 1.0e-14;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		bool rejected = false;
		try { executor.Advance(step, {{"left", 0.0}, {"right", 0.0}}); }
		catch (const iga::PressureFlowConvergenceError& error) {
			assert(error.Result().iterations.size() == 2);
			assert(error.Result().iterations.back().edges.size() == 2);
			assert(error.Result().iterations.back().edges.front().measured_pressure_pa == 10.0);
			rejected = true;
		}
		assert(rejected);
		assert(fixture.up->commits == 0 && fixture.up->aborts == 1);
	}
	{
		Fixture fixture(graph);
		fixture.up->nonfinite_state = true;
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step, {{"left", 10.0}, {"right", 20.0}}); });
		assert(fixture.up->commits == 0 && fixture.mid->commits == 0
			&& fixture.down->commits == 0);
		assert(fixture.up->aborts == 1 && fixture.mid->aborts == 1
			&& fixture.down->aborts == 1);
	}
	{
		Fixture fixture(graph);
		fixture.mid->fail_prepare = true;
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step, {{"left", 10.0}, {"right", 20.0}}); });
		assert(fixture.up->commits == 0 && fixture.mid->commits == 0
			&& fixture.down->commits == 0);
		assert(fixture.up->aborts == 1 && fixture.mid->aborts == 1
			&& fixture.down->aborts == 1);
	}
	{
		Fixture fixture(graph);
		fixture.down->fail_solve = true;
		fixture.mid->fail_abort = true;
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		bool rejected = false;
		try { executor.Advance(step, {{"left", 10.0}, {"right", 20.0}}); }
		catch (const std::runtime_error& error) {
			const std::string message = error.what();
			assert(message.find("injected fake solve failure") != std::string::npos);
			assert(message.find("injected fake abort failure") != std::string::npos);
			rejected = true;
		}
		assert(rejected);
		assert(fixture.up->commits == 0 && fixture.mid->commits == 0
			&& fixture.down->commits == 0);
		assert(fixture.up->aborts == 1 && fixture.mid->aborts == 1
			&& fixture.down->aborts == 1);
	}
	{
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		RequireRejected([&] {
			iga::PressureFlowComponentExecutor reverse(*fixture.registry, "down", controls);
		});
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step,
			{{"c_right", 0.0}, {"a_source", 0.0}, {"b_left", 0.0}});
		assert((executor.Plan().domain_order == std::vector<std::string>{"up", "junction",
			"branch_a", "branch_b"}));
		assert(result.iterations.size() == 1);
		assert(result.iterations.front().edges.size() == 3);
		assert(result.iterations.front().edges[0].edge_id == "a_source");
		assert(result.iterations.front().edges[1].edge_id == "b_left");
		assert(result.iterations.front().edges[2].edge_id == "c_right");
		assert(*result.accepted_ports.at({"junction", "left"}).outward_flow_m3_s == 1.0);
		assert(*result.accepted_ports.at({"branch_a", "root"}).outward_flow_m3_s == -1.0);
		assert(*result.accepted_ports.at({"junction", "right"}).outward_flow_m3_s == 2.0);
		assert(*result.accepted_ports.at({"branch_b", "root"}).outward_flow_m3_s == -2.0);
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->solves == 1);
			assert(runtime.second->commits == 1);
		}
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Fixed;
		controls.maximum_iterations = 40;
		controls.pressure_relative_tolerance = 1.0e-7;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}, {"c_right", 0.0}});
		assert(result.iterations.back().converged);
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->solves == static_cast<int>(result.iterations.size()));
			assert(runtime.second->rollbacks+1 == runtime.second->solves);
		}
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Aitken;
		controls.maximum_iterations = 8;
		controls.pressure_relative_tolerance = 1.0e-12;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		const auto result = executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}, {"c_right", 0.0}});
		assert(result.iterations.size() == 3 && result.iterations.back().converged);
		assert(fixture.by_id.at("junction")->solves == 3);
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		fixture.by_id.at("branch_b")->fail_solve = true;
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}, {"c_right", 0.0}}); });
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->commits == 0);
			assert(runtime.second->aborts == 1);
		}
	}
	{
		const auto branch = Branch(true);
		BranchFixture fixture(branch);
		fixture.by_id.at("branch_a")->fail_prepare = true;
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}, {"c_right", 0.0}}); });
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->commits == 0);
			assert(runtime.second->aborts == 1);
		}
	}
	{
		const auto canonical_graph = Branch();
		const auto permuted_graph = Branch(true);
		BranchFixture canonical(canonical_graph);
		BranchFixture permuted(permuted_graph);
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor canonical_executor(*canonical.registry, "up", controls);
		iga::PressureFlowComponentExecutor permuted_executor(*permuted.registry, "up", controls);
		const std::map<std::string, double> pressure{{"a_source", 4.0},
			{"b_left", 5.0}, {"c_right", 6.0}};
		const auto canonical_result = canonical_executor.Advance(step, pressure);
		const auto permuted_result = permuted_executor.Advance(step, pressure);
		for (const auto& state : canonical_result.accepted_ports) {
			const auto& other = permuted_result.accepted_ports.at(state.first);
			assert(state.second.outward_flow_m3_s == other.outward_flow_m3_s);
			assert(state.second.mean_pressure_pa == other.mean_pressure_pa);
		}
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}}); });
		for (const auto& runtime : fixture.by_id)
			assert(runtime.second->begins == 0);
	}
	{
		const auto branch = Branch();
		BranchFixture fixture(branch);
		iga::PressureFlowExecutionControls controls;
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls);
		RequireRejected([&] { executor.Advance(step,
			{{"a_source", 0.0}, {"b_left", 0.0}, {"c_right", 0.0}},
			[](const iga::PressureFlowStepResult&) {
				throw std::runtime_error("injected branch callback failure");
			}); });
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->commits == 0);
			assert(runtime.second->aborts == 1);
		}
	}
	{
		auto runtimes = RuntimeSet(graph, iga::DomainKind::ThreeDBodyFittedFlow,
			graph.Domain("up").ports);
		RequireRejected([&] { iga::DomainRuntimeRegistry registry(graph, std::move(runtimes)); });
	}
	{
		auto ports = graph.Domain("up").ports;
		ports.front().locator = "different";
		auto runtimes = RuntimeSet(graph, iga::DomainKind::OneDFlow, std::move(ports));
		RequireRejected([&] { iga::DomainRuntimeRegistry registry(graph, std::move(runtimes)); });
	}
	{
		auto ports = graph.Domain("up").ports;
		ports.front().orientation.native_to_outward_sign = -1;
		auto runtimes = RuntimeSet(graph, iga::DomainKind::OneDFlow, std::move(ports));
		RequireRejected([&] { iga::DomainRuntimeRegistry registry(graph, std::move(runtimes)); });
	}
	std::cout << "pressure-flow component executor tests passed\n";
	return 0;
}

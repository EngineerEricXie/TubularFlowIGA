#include "SimulationGraph.hpp"

#include <cassert>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
}

iga::CouplingPort LogicalPort(const std::string& domain_id, const std::string& id,
	const std::string& locator, iga::PortQuantity accepted)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = domain_id;
	port.locator_kind = "runtime_port";
	port.locator = locator;
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires = {accepted};
	iga::ValidateCouplingPort(port);
	return port;
}

iga::DomainNode Domain(const std::string& id, iga::DomainKind kind,
	std::vector<iga::CouplingPort> ports)
{
	return {id, kind, std::move(ports)};
}

iga::SimulationGraph StraightChain()
{
	auto upstream = Domain("upstream", iga::DomainKind::OneDFlow,
		{LogicalPort("upstream", "terminal", "outlet:2", iga::PortQuantity::MeanPressure)});
	auto three_d = Domain("roi", iga::DomainKind::ThreeDBodyFittedFlow,
		{LogicalPort("roi", "inlet", "boundary:1", iga::PortQuantity::FlowRate),
		 LogicalPort("roi", "outlet", "boundary:2", iga::PortQuantity::MeanPressure)});
	auto downstream = Domain("downstream", iga::DomainKind::OneDFlow,
		{LogicalPort("downstream", "root", "root", iga::PortQuantity::FlowRate)});
	iga::CouplingEdge left{"upstream_to_roi", {"upstream", "terminal"}, {"roi", "inlet"},
		iga::CouplingLaw::PressureFlow};
	iga::CouplingEdge right{"roi_to_downstream", {"roi", "outlet"}, {"downstream", "root"},
		iga::CouplingLaw::PressureFlow};
	return iga::SimulationGraph({std::move(downstream), std::move(three_d), std::move(upstream)},
		{std::move(right), std::move(left)});
}

} // namespace

int main()
{
	const auto graph = StraightChain();
	assert(graph.Domains().size() == 3);
	assert(graph.Edges().size() == 2);
	assert(graph.Domain("roi").kind == iga::DomainKind::ThreeDBodyFittedFlow);
	assert(graph.Port({"roi", "inlet"}).requires.count(iga::PortQuantity::FlowRate));
	assert(graph.Port({"roi", "inlet"}).provides.count(iga::PortQuantity::FlowRate));
	assert(graph.Edge("upstream_to_roi").law == iga::CouplingLaw::PressureFlow);
	const auto plan = iga::MakeSequentialPlan(graph, "upstream");
	assert((plan.domain_ids == std::vector<std::string>{"upstream", "roi", "downstream"}));
	assert((plan.edge_ids == std::vector<std::string>{"upstream_to_roi", "roi_to_downstream"}));
	const auto reverse = iga::MakeSequentialPlan(graph, "downstream");
	assert((reverse.domain_ids == std::vector<std::string>{"downstream", "roi", "upstream"}));
	{
		std::vector<iga::DomainNode> copied_domains;
		for (const auto& domain : graph.Domains()) copied_domains.push_back(domain.second);
		auto copied_edges = graph.Edges();
		const iga::SimulationGraph immutable(copied_domains, copied_edges);
		copied_domains.front().id = "mutated";
		copied_edges.front().id = "mutated";
		assert(immutable.Domains().count("mutated") == 0);
		assert(immutable.Edge("roi_to_downstream").id == "roi_to_downstream");
	}

	RequireRejected([&graph] { (void)graph.Domain("missing"); });
	RequireRejected([&graph] { (void)graph.Port({"roi", "missing"}); });
	RequireRejected([&graph] { (void)graph.Edge("missing"); });
	RequireRejected([&graph] { (void)iga::MakeSequentialPlan(graph, "missing"); });
	RequireRejected([&graph] { (void)iga::MakeSequentialPlan(graph, "roi"); });
	RequireRejected([] { iga::SimulationGraph({}, {}); });
	RequireRejected([] {
		iga::SimulationGraph({{"", iga::DomainKind::OneDFlow, {}}}, {});
	});
	RequireRejected([] {
		iga::SimulationGraph({{"empty", iga::DomainKind::OneDFlow, {}}}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("unknown", "port", "port", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("unknown", static_cast<iga::DomainKind>(99), {port})}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("one", "terminal", "outlet:2", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {port}),
			Domain("one", iga::DomainKind::OneDFlow, {port})}, {});
	});
	RequireRejected([] {
		auto port = LogicalPort("wrong", "terminal", "outlet:2", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {port})}, {});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "same", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("one", "b", "same", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first, second})}, {});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"bad", {"one", "missing"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"empty_endpoint", {"", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"ambiguous", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		second.provides.erase(iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"missing_diagnostic", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"same", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow},
			 {"same", {"one", "a"}, {"two", "b"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first})},
			{{"self", {"one", "a"}, {"one", "a"}, iga::CouplingLaw::PressureFlow}});
	});
	RequireRejected([] {
		auto first = LogicalPort("one", "a", "a", iga::PortQuantity::MeanPressure);
		auto second = LogicalPort("two", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph({Domain("one", iga::DomainKind::OneDFlow, {first}),
			Domain("two", iga::DomainKind::ThreeDBodyFittedFlow, {second})},
			{{"unknown_law", {"one", "a"}, {"two", "b"},
				static_cast<iga::CouplingLaw>(99)}});
	});
	RequireRejected([] {
		auto source = Domain("source", iga::DomainKind::OneDFlow,
			{LogicalPort("source", "out", "out", iga::PortQuantity::MeanPressure)});
		auto first = Domain("first", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("first", "in", "in", iga::PortQuantity::FlowRate)});
		auto second = Domain("second", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("second", "in", "in", iga::PortQuantity::FlowRate)});
		iga::SimulationGraph({std::move(source), std::move(first), std::move(second)},
			{{"first_edge", {"source", "out"}, {"first", "in"}, iga::CouplingLaw::PressureFlow},
			 {"second_edge", {"source", "out"}, {"second", "in"}, iga::CouplingLaw::PressureFlow}});
	});

	{
		auto center = Domain("center", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("center", "a", "a", iga::PortQuantity::FlowRate),
			 LogicalPort("center", "b", "b", iga::PortQuantity::FlowRate),
			 LogicalPort("center", "c", "c", iga::PortQuantity::FlowRate)});
		std::vector<iga::DomainNode> domains{std::move(center)};
		std::vector<iga::CouplingEdge> edges;
		for (const auto& id : {"a", "b", "c"}) {
			domains.push_back(Domain(id, iga::DomainKind::OneDFlow,
				{LogicalPort(id, "terminal", "terminal", iga::PortQuantity::MeanPressure)}));
			edges.push_back({std::string("edge_")+id, {"center", id}, {id, "terminal"},
				iga::CouplingLaw::PressureFlow});
		}
		const iga::SimulationGraph branch(std::move(domains), std::move(edges));
		assert(branch.Domains().size() == 4);
		RequireRejected([&branch] { (void)iga::MakeSequentialPlan(branch, "a"); });
	}
	{
		auto graph_with_island = StraightChain();
		std::vector<iga::DomainNode> domains;
		for (const auto& domain : graph_with_island.Domains()) domains.push_back(domain.second);
		domains.push_back(Domain("island", iga::DomainKind::ThreeDBodyFittedFlow,
			{LogicalPort("island", "unused", "unused", iga::PortQuantity::FlowRate)}));
		iga::SimulationGraph disconnected(std::move(domains), graph_with_island.Edges());
		RequireRejected([&disconnected] {
			(void)iga::MakeSequentialPlan(disconnected, "upstream");
		});
	}
	{
		auto left = LogicalPort("left", "a", "a", iga::PortQuantity::MeanPressure);
		auto right = LogicalPort("right", "b", "b", iga::PortQuantity::FlowRate);
		iga::SimulationGraph same_kind({Domain("left", iga::DomainKind::OneDFlow, {left}),
			Domain("right", iga::DomainKind::OneDFlow, {right})},
			{{"same_kind", {"left", "a"}, {"right", "b"}, iga::CouplingLaw::PressureFlow}});
		RequireRejected([&same_kind] { (void)iga::MakeSequentialPlan(same_kind, "left"); });
	}
	{
		std::vector<iga::DomainNode> domains;
		for (int index = 0; index < 4; ++index) {
			const std::string id = "cycle_"+std::to_string(index);
			const bool one_d = index%2 == 0;
			domains.push_back(Domain(id, one_d ? iga::DomainKind::OneDFlow
				: iga::DomainKind::ThreeDBodyFittedFlow,
				{LogicalPort(id, "left", "left", one_d ? iga::PortQuantity::MeanPressure
					: iga::PortQuantity::FlowRate),
				 LogicalPort(id, "right", "right", one_d ? iga::PortQuantity::MeanPressure
					: iga::PortQuantity::FlowRate)}));
		}
		const iga::SimulationGraph cycle(std::move(domains), {
			{"cycle_edge_0", {"cycle_0", "right"}, {"cycle_1", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_1", {"cycle_1", "right"}, {"cycle_2", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_2", {"cycle_2", "right"}, {"cycle_3", "left"}, iga::CouplingLaw::PressureFlow},
			{"cycle_edge_3", {"cycle_3", "right"}, {"cycle_0", "left"}, iga::CouplingLaw::PressureFlow}});
		RequireRejected([&cycle] { (void)iga::MakeSequentialPlan(cycle, "cycle_0"); });
	}
	return 0;
}

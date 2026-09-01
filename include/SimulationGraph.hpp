#ifndef IGA_SIMULATION_GRAPH_HPP
#define IGA_SIMULATION_GRAPH_HPP

#include "CouplingEdge.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct DomainNode {
	std::string id;
	DomainKind kind = DomainKind::OneDFlow;
	std::vector<CouplingPort> ports;
};

struct SequentialCouplingPlan {
	std::vector<std::string> domain_ids;
	std::vector<std::string> edge_ids;
};

struct PressureFlowInterfacePlan {
	std::string edge_id;
	PortRef flow_provider;
	PortRef flow_receiver;
	PortRef pressure_provider;
	PortRef pressure_receiver;
};

struct PressureFlowComponentPlan {
	std::vector<std::string> domain_order;
	std::vector<PressureFlowInterfacePlan> interfaces;
};

class SimulationGraph {
public:
	SimulationGraph(std::vector<DomainNode> domains, std::vector<CouplingEdge> edges)
		: edges_(std::move(edges))
	{
		if (domains.empty()) throw std::runtime_error("simulation graph requires at least one domain");
		for (auto& domain : domains) {
			ValidateDomain(domain);
			if (!domains_.emplace(domain.id, std::move(domain)).second)
				throw std::runtime_error("simulation graph domain ids must be unique");
		}
		std::set<std::string> edge_ids;
		std::set<PortRef> attached_ports;
		std::set<std::pair<PortRef, PortRef>> endpoint_pairs;
		for (const auto& edge : edges_) {
			ValidateCouplingEdge(edge);
			if (!edge_ids.insert(edge.id).second)
				throw std::runtime_error("simulation graph coupling edge ids must be unique");
			const auto& first = Port(edge.first);
			const auto& second = Port(edge.second);
			if (!attached_ports.insert(edge.first).second
				|| !attached_ports.insert(edge.second).second)
				throw std::runtime_error("simulation graph port is attached to multiple edges");
			auto pair = std::make_pair(edge.first, edge.second);
			if (pair.second < pair.first) std::swap(pair.first, pair.second);
			if (!endpoint_pairs.insert(pair).second)
				throw std::runtime_error("simulation graph has a duplicate endpoint pair");
			if (edge.law == CouplingLaw::PressureFlow)
				ValidatePressureFlowEdgePorts(first, second);
		}
	}

	const DomainNode& Domain(const std::string& id) const
	{
		const auto found = domains_.find(id);
		if (found == domains_.end())
			throw std::runtime_error("simulation graph references unknown domain '"+id+"'");
		return found->second;
	}

	const CouplingPort& Port(const PortRef& reference) const
	{
		ValidatePortRef(reference);
		const auto& domain = Domain(reference.domain_id);
		const auto found = std::find_if(domain.ports.begin(), domain.ports.end(),
			[&](const CouplingPort& port) { return port.id == reference.port_id; });
		if (found == domain.ports.end())
			throw std::runtime_error("simulation graph references unknown port '"
				+reference.domain_id+"."+reference.port_id+"'");
		return *found;
	}

	const CouplingEdge& Edge(const std::string& id) const
	{
		const auto found = std::find_if(edges_.begin(), edges_.end(),
			[&](const CouplingEdge& edge) { return edge.id == id; });
		if (found == edges_.end())
			throw std::runtime_error("simulation graph references unknown edge '"+id+"'");
		return *found;
	}

	const std::map<std::string, DomainNode>& Domains() const noexcept { return domains_; }
	const std::vector<CouplingEdge>& Edges() const noexcept { return edges_; }

private:
	static void ValidateDomain(const DomainNode& domain)
	{
		if (domain.id.empty()) throw std::runtime_error("simulation graph domain id must be nonempty");
		if (!IsKnownDomainKind(domain.kind))
			throw std::runtime_error("simulation graph domain has an unsupported kind");
		if (domain.ports.empty())
			throw std::runtime_error("simulation graph domain requires at least one logical port");
		ValidateCouplingPorts(domain.ports);
		std::set<std::pair<std::string, std::string>> locators;
		for (const auto& port : domain.ports) {
			if (port.subsystem_id != domain.id)
				throw std::runtime_error("simulation graph port subsystem_id does not match its domain");
			if (!locators.emplace(port.locator_kind, port.locator).second)
				throw std::runtime_error("simulation graph logical ports must have unique physical locators within a domain");
		}
	}

	std::map<std::string, DomainNode> domains_;
	std::vector<CouplingEdge> edges_;
};

inline SequentialCouplingPlan MakeSequentialPlan(const SimulationGraph& graph,
	const std::string& start_domain_id)
{
	(void)graph.Domain(start_domain_id);
	if (graph.Domains().size() == 1) return {{start_domain_id}, {}};
	std::map<std::string, std::vector<std::pair<std::string, std::string>>> adjacency;
	for (const auto& domain : graph.Domains()) adjacency.emplace(domain.first,
		std::vector<std::pair<std::string, std::string>>{});
	for (const auto& edge : graph.Edges()) {
		const auto& first = graph.Domain(edge.first.domain_id);
		const auto& second = graph.Domain(edge.second.domain_id);
		if (first.kind == second.kind)
			throw std::runtime_error("sequential coupling plan requires heterogeneous neighboring domains");
		adjacency.at(first.id).emplace_back(second.id, edge.id);
		adjacency.at(second.id).emplace_back(first.id, edge.id);
	}
	for (auto& neighbors : adjacency) {
		std::sort(neighbors.second.begin(), neighbors.second.end());
		if (neighbors.second.size() > 2)
			throw std::runtime_error("sequential coupling plan does not support graph branches");
	}
	if (adjacency.at(start_domain_id).size() != 1)
		throw std::runtime_error("sequential coupling plan must start at a chain endpoint");
	SequentialCouplingPlan plan;
	std::set<std::string> visited;
	std::string current = start_domain_id;
	std::string previous;
	while (true) {
		if (!visited.insert(current).second)
			throw std::runtime_error("sequential coupling plan does not support graph cycles");
		plan.domain_ids.push_back(current);
		const auto& neighbors = adjacency.at(current);
		auto next = neighbors.end();
		for (auto candidate = neighbors.begin(); candidate != neighbors.end(); ++candidate)
			if (candidate->first != previous) {
				if (next != neighbors.end())
					throw std::runtime_error("sequential coupling plan has an ambiguous next domain");
				next = candidate;
			}
		if (next == neighbors.end()) break;
		plan.edge_ids.push_back(next->second);
		previous = current;
		current = next->first;
	}
	if (visited.size() != graph.Domains().size())
		throw std::runtime_error("sequential coupling plan requires one connected graph component");
	return plan;
}

inline void ValidateConnectedGraph(const SimulationGraph& graph,
	const std::string& start_domain_id)
{
	(void)graph.Domain(start_domain_id);
	std::map<std::string, std::vector<std::string>> adjacency;
	for (const auto& domain : graph.Domains()) adjacency.emplace(domain.first,
		std::vector<std::string>{});
	for (const auto& edge : graph.Edges()) {
		adjacency.at(edge.first.domain_id).push_back(edge.second.domain_id);
		adjacency.at(edge.second.domain_id).push_back(edge.first.domain_id);
	}
	std::vector<std::string> pending{start_domain_id};
	std::set<std::string> visited;
	while (!pending.empty()) {
		const auto domain = pending.back();
		pending.pop_back();
		if (!visited.insert(domain).second) continue;
		for (const auto& neighbor : adjacency.at(domain)) pending.push_back(neighbor);
	}
	if (visited.size() != graph.Domains().size())
		throw std::runtime_error("simulation graph must contain one connected component");
}

inline SequentialCouplingPlan MakeSequentialPressureFlowPlan(const SimulationGraph& graph,
	const std::string& start_domain_id)
{
	auto plan = MakeSequentialPlan(graph, start_domain_id);
	for (std::size_t i = 0; i < plan.edge_ids.size(); ++i) {
		const auto& edge = graph.Edge(plan.edge_ids[i]);
		const auto& earlier = plan.domain_ids[i];
		const auto& later = plan.domain_ids[i+1];
		const auto earlier_ref = edge.first.domain_id == earlier ? edge.first : edge.second;
		const auto later_ref = edge.first.domain_id == later ? edge.first : edge.second;
		if (earlier_ref.domain_id != earlier || later_ref.domain_id != later)
			throw std::runtime_error("sequential edge does not join adjacent plan domains");
		const auto& earlier_port = graph.Port(earlier_ref);
		const auto& later_port = graph.Port(later_ref);
		if (!earlier_port.requires.count(PortQuantity::MeanPressure)
			|| !later_port.requires.count(PortQuantity::FlowRate))
			throw std::runtime_error(
				"sequential pressure-flow plan requires pressure receivers before flow receivers");
	}
	return plan;
}

inline PressureFlowComponentPlan MakeAcyclicPressureFlowPlan(
	const SimulationGraph& graph, const std::string& start_domain_id)
{
	ValidateConnectedGraph(graph, start_domain_id);
	if (graph.Edges().empty())
		throw std::runtime_error("pressure-flow component requires at least one interface");
	PressureFlowComponentPlan plan;
	std::map<std::string, int> indegree;
	std::map<std::string, std::vector<std::string>> downstream;
	for (const auto& domain : graph.Domains()) {
		indegree.emplace(domain.first, 0);
		downstream.emplace(domain.first, std::vector<std::string>{});
	}
	for (const auto& edge : graph.Edges()) {
		if (graph.Domain(edge.first.domain_id).kind
			== graph.Domain(edge.second.domain_id).kind)
			throw std::runtime_error(
				"pressure-flow component requires heterogeneous neighboring domains");
		const auto& first = graph.Port(edge.first);
		const bool first_receives_pressure
			= first.requires.count(PortQuantity::MeanPressure) != 0;
		const auto pressure_receiver = first_receives_pressure ? edge.first : edge.second;
		const auto flow_receiver = first_receives_pressure ? edge.second : edge.first;
		plan.interfaces.push_back({edge.id, pressure_receiver, flow_receiver,
			flow_receiver, pressure_receiver});
		downstream.at(pressure_receiver.domain_id).push_back(flow_receiver.domain_id);
		++indegree.at(flow_receiver.domain_id);
	}
	std::sort(plan.interfaces.begin(), plan.interfaces.end(),
		[](const PressureFlowInterfacePlan& first,
			const PressureFlowInterfacePlan& second) {
			return first.edge_id < second.edge_id;
		});
	std::set<std::string> ready;
	for (const auto& domain : indegree)
		if (domain.second == 0) ready.insert(domain.first);
	if (ready.size() != 1 || *ready.begin() != start_domain_id)
		throw std::runtime_error(
			"pressure-flow component requires the declared start domain to be its unique source");
	while (!ready.empty()) {
		const auto domain = *ready.begin();
		ready.erase(ready.begin());
		plan.domain_order.push_back(domain);
		auto children = downstream.at(domain);
		std::sort(children.begin(), children.end());
		for (const auto& child : children) {
			auto& remaining = indegree.at(child);
			--remaining;
			if (remaining == 0) ready.insert(child);
		}
	}
	if (plan.domain_order.size() != graph.Domains().size())
		throw std::runtime_error("pressure-flow component does not support graph cycles");
	return plan;
}

} // namespace iga

#endif

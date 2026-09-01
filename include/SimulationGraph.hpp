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

} // namespace iga

#endif

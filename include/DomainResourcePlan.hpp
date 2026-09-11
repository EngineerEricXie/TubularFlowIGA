#ifndef IGA_DOMAIN_RESOURCE_PLAN_HPP
#define IGA_DOMAIN_RESOURCE_PLAN_HPP

#include "SimulationGraph.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class DomainResourceMode { Shared, DomainGroups };

struct DomainResourceGroupDefinition {
	std::string id;
	int ranks = 0;
	std::vector<std::string> domains;
};

struct DomainResourceDefinition {
	DomainResourceMode mode = DomainResourceMode::Shared;
	std::vector<DomainResourceGroupDefinition> groups;
};

struct DomainResourceGroup {
	std::string id;
	std::vector<int> world_ranks;
	std::vector<std::string> domains;
};

struct DomainResourcePlan {
	DomainResourceMode mode = DomainResourceMode::Shared;
	std::vector<DomainResourceGroup> groups;
	std::map<std::string, std::size_t> domain_group;
	std::vector<std::vector<std::string>> hydraulic_batches;

	const DomainResourceGroup& GroupFor(const std::string& domain_id) const
	{
		const auto found = domain_group.find(domain_id);
		if (found == domain_group.end())
			throw std::runtime_error("resource plan has no group for domain '"+domain_id+"'");
		return groups.at(found->second);
	}
};

inline std::vector<std::vector<std::string>> MakeDomainSolveBatches(
	const std::vector<std::string>& order,
	const std::vector<std::pair<std::string, std::string>>& dependencies,
	const std::map<std::string, std::size_t>& domain_group)
{
	std::map<std::string, int> depth;
	for (const auto& domain : order) {
		if (!domain_group.count(domain))
			throw std::runtime_error("domain solve schedule has no resource group for '"+domain+"'");
		int value = 0;
		for (const auto& edge : dependencies) if (edge.second == domain) {
			const auto provider = depth.find(edge.first);
			if (provider == depth.end())
				throw std::runtime_error("domain solve order violates a dependency");
			value = std::max(value, provider->second+1);
		}
		if (!depth.emplace(domain, value).second)
			throw std::runtime_error("domain solve order contains a duplicate domain");
	}
	for (const auto& edge : dependencies)
		if (!depth.count(edge.first) || !depth.count(edge.second))
			throw std::runtime_error("domain solve dependency references an unscheduled domain");
	std::map<int, std::vector<std::string>> levels;
	for (const auto& domain : order) levels[depth.at(domain)].push_back(domain);
	std::vector<std::vector<std::string>> result;
	for (const auto& level : levels) {
		std::vector<std::vector<std::string>> waves;
		std::vector<std::set<std::size_t>> used;
		for (const auto& domain : level.second) {
			const auto group = domain_group.at(domain);
			std::size_t wave = 0;
			while (wave < used.size() && used[wave].count(group)) ++wave;
			if (wave == waves.size()) {
				waves.emplace_back();
				used.emplace_back();
			}
			waves[wave].push_back(domain);
			used[wave].insert(group);
		}
		result.insert(result.end(), waves.begin(), waves.end());
	}
	return result;
}

inline DomainResourcePlan MakeDomainResourcePlan(const SimulationGraph& graph,
	const std::string& start_domain_id, const DomainResourceDefinition& definition,
	int world_ranks)
{
	if (world_ranks < 1) throw std::runtime_error("domain resource plan requires at least one rank");
	DomainResourcePlan result;
	result.mode = definition.mode;
	if (definition.mode == DomainResourceMode::Shared) {
		if (!definition.groups.empty())
			throw std::runtime_error("shared resource mode cannot declare domain groups");
		DomainResourceGroup group;
		group.id = "shared";
		for (int rank = 0; rank < world_ranks; ++rank) group.world_ranks.push_back(rank);
		for (const auto& domain : graph.Domains()) {
			group.domains.push_back(domain.first);
			result.domain_group.emplace(domain.first, 0);
		}
		result.groups.push_back(std::move(group));
	} else {
		if (definition.groups.empty())
			throw std::runtime_error("domain_groups resource mode requires groups");
		int next_rank = 0;
		std::set<std::string> group_ids;
		for (const auto& requested : definition.groups) {
			if (requested.id.empty() || !group_ids.insert(requested.id).second)
				throw std::runtime_error("domain resource group ids must be nonempty and unique");
			if (requested.ranks < 1 || next_rank+requested.ranks > world_ranks)
				throw std::runtime_error("domain resource group rank request exceeds allocation");
			if (requested.domains.empty())
				throw std::runtime_error("domain resource group must own at least one domain");
			DomainResourceGroup group;
			group.id = requested.id;
			group.domains = requested.domains;
			for (int rank = 0; rank < requested.ranks; ++rank)
				group.world_ranks.push_back(next_rank++);
			const auto index = result.groups.size();
			for (const auto& domain : group.domains) {
				(void)graph.Domain(domain);
				if (!result.domain_group.emplace(domain, index).second)
					throw std::runtime_error("domain resource groups must own each domain exactly once");
			}
			result.groups.push_back(std::move(group));
		}
		if (result.domain_group.size() != graph.Domains().size())
			throw std::runtime_error("domain resource groups require exact graph-domain coverage");
	}

	const auto flow = MakeAcyclicPressureFlowPlan(graph, start_domain_id);
	std::vector<std::pair<std::string, std::string>> dependencies;
	for (const auto& edge : flow.interfaces)
		dependencies.emplace_back(edge.flow_provider.domain_id, edge.flow_receiver.domain_id);
	result.hydraulic_batches = MakeDomainSolveBatches(flow.domain_order,
		dependencies, result.domain_group);
	return result;
}

} // namespace iga

#endif

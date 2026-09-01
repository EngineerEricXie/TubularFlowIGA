#ifndef IGA_BIFURCATION_MULTIDOMAIN_CASE_HPP
#define IGA_BIFURCATION_MULTIDOMAIN_CASE_HPP

#include "SequentialMultidomainCase.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OneDThreeDBranchDefinition {
	std::string edge_id;
	std::string downstream_domain_id;
	PortRef three_d_outlet;
	PortRef downstream_root;
	PortRef downstream_terminal_observation;
};

struct OneDThreeDBifurcationDefinition {
	std::string upstream_domain_id;
	std::string three_d_domain_id;
	std::string upstream_edge_id;
	PortRef upstream_terminal;
	PortRef upstream_root_observation;
	PortRef three_d_inlet;
	PortRef three_d_wall_observation;
	std::vector<OneDThreeDBranchDefinition> branches;
};

inline PortRef UniqueOneDTerminalObservation(const SimulationGraph& graph,
	const std::string& domain_id)
{
	const CouplingPort* terminal = nullptr;
	for (const auto& port : graph.Domain(domain_id).ports) {
		if (port.locator_kind != "runtime_port" || port.locator == "root"
			|| !port.requires.empty()) continue;
		(void)ParseOneDOutletNodeLocator(port);
		if (terminal)
			throw std::runtime_error(
				"bifurcation runner requires exactly one terminal observation per downstream 1D branch");
		terminal = &port;
	}
	if (!terminal)
		throw std::runtime_error(
			"bifurcation runner requires one terminal observation per downstream 1D branch");
	return {domain_id, terminal->id};
}

inline OneDThreeDBifurcationDefinition ResolveOneDThreeDBifurcation(
	const MultidomainConfiguration& configuration)
{
	const auto plan = MakeAcyclicPressureFlowPlan(configuration.graph,
		configuration.start_domain_id);
	if (plan.domain_order.size() < 4 || plan.interfaces.size() < 3)
		throw std::runtime_error(
			"bifurcation runner requires one source, one 3D junction, and at least two 1D branches");
	OneDThreeDBifurcationDefinition result;
	result.upstream_domain_id = plan.domain_order.front();
	if (configuration.graph.Domain(result.upstream_domain_id).kind != DomainKind::OneDFlow
		|| GraphDomainDefinitionFor(configuration, result.upstream_domain_id).one_d_inlet_policy
			!= OneDInletPolicy::ConfiguredOpenLoop)
		throw std::runtime_error(
			"bifurcation runner requires a configured-open-loop 1D source");
	for (const auto& domain : configuration.graph.Domains())
		if (domain.second.kind == DomainKind::ThreeDBodyFittedFlow) {
			if (!result.three_d_domain_id.empty())
				throw std::runtime_error("bifurcation runner requires exactly one 3D junction");
			result.three_d_domain_id = domain.first;
		}
	if (result.three_d_domain_id.empty())
		throw std::runtime_error("bifurcation runner requires exactly one 3D junction");
	for (const auto& interface : plan.interfaces) {
		if (interface.flow_provider.domain_id == result.upstream_domain_id
			&& interface.flow_receiver.domain_id == result.three_d_domain_id) {
			if (!result.upstream_edge_id.empty())
				throw std::runtime_error("bifurcation runner requires exactly one 3D inlet");
			result.upstream_edge_id = interface.edge_id;
			result.upstream_terminal = interface.flow_provider;
			result.three_d_inlet = interface.flow_receiver;
			continue;
		}
		if (interface.flow_provider.domain_id != result.three_d_domain_id)
			throw std::runtime_error("bifurcation runner supports only one-level 3D branching");
		const auto& downstream = configuration.graph.Domain(interface.flow_receiver.domain_id);
		const auto& definition = GraphDomainDefinitionFor(configuration, downstream.id);
		if (downstream.kind != DomainKind::OneDFlow
			|| definition.one_d_inlet_policy != OneDInletPolicy::CoupledRoot
			|| interface.flow_receiver.port_id.empty())
			throw std::runtime_error(
				"bifurcation runner requires coupled-root 1D branch leaves");
		const auto& root = configuration.graph.Port(interface.flow_receiver);
		if (root.locator_kind != "runtime_port" || root.locator != "root")
			throw std::runtime_error("bifurcation branch input must be a 1D root runtime port");
		result.branches.push_back({interface.edge_id, downstream.id,
			interface.flow_provider, interface.flow_receiver,
			UniqueOneDTerminalObservation(configuration.graph, downstream.id)});
	}
	if (result.upstream_edge_id.empty() || result.branches.size() < 2
		|| plan.domain_order.size() != result.branches.size()+2)
		throw std::runtime_error(
			"bifurcation runner requires one 1D source--3D inlet and at least two 1D branch leaves");
	const auto& upstream_port = configuration.graph.Port(result.upstream_terminal);
	if (upstream_port.locator_kind != "runtime_port" || upstream_port.locator == "root")
		throw std::runtime_error("bifurcation source port must be a 1D outlet runtime port");
	(void)ParseOneDOutletNodeLocator(upstream_port);
	result.upstream_root_observation = UniqueObservationPort(configuration.graph,
		result.upstream_domain_id, "runtime_port", "root");
	result.three_d_wall_observation = UniqueObservationPort(configuration.graph,
		result.three_d_domain_id, "boundary_label", "0");
	std::sort(result.branches.begin(), result.branches.end(),
		[](const OneDThreeDBranchDefinition& first,
			const OneDThreeDBranchDefinition& second) {
			return first.edge_id < second.edge_id;
		});
	return result;
}

} // namespace iga

#endif

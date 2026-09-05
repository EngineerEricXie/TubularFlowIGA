#ifndef IGA_SEQUENTIAL_MULTIDOMAIN_CASE_HPP
#define IGA_SEQUENTIAL_MULTIDOMAIN_CASE_HPP

#include "MultidomainConfig.hpp"
#include "PressureFlowComponentExecutor.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"

#include <algorithm>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>

namespace iga {

struct ResolvedGraphDomainAssets {
	std::filesystem::path case_directory;
	std::filesystem::path database;
};

struct SequentialOneDThreeDOneDDefinition {
	std::string upstream_domain_id;
	std::string three_d_domain_id;
	std::string downstream_domain_id;
	std::string upstream_edge_id;
	std::string downstream_edge_id;
	PortRef upstream_terminal;
	PortRef upstream_root_observation;
	PortRef three_d_inlet;
	PortRef three_d_outlet;
	PortRef three_d_wall_observation;
	PortRef downstream_root;
	PortRef downstream_terminal_observation;
};

inline const GraphDomainDefinition& GraphDomainDefinitionFor(
	const MultidomainConfiguration& configuration, const std::string& id)
{
	const auto found = std::find_if(configuration.domains.begin(), configuration.domains.end(),
		[&](const GraphDomainDefinition& domain) { return domain.id == id; });
	if (found == configuration.domains.end())
		throw std::runtime_error("multidomain configuration has no domain definition for '"
			+id+"'");
	return *found;
}

inline PortRef EdgeEndpointForDomain(const CouplingEdge& edge, const std::string& domain_id)
{
	if (edge.first.domain_id == domain_id) return edge.first;
	if (edge.second.domain_id == domain_id) return edge.second;
	throw std::runtime_error("coupling edge does not contain expected domain '"+domain_id+"'");
}

inline PortRef UniqueObservationPort(const SimulationGraph& graph,
	const std::string& domain_id, const std::string& locator_kind,
	const std::string& locator)
{
	const CouplingPort* match = nullptr;
	for (const auto& port : graph.Domain(domain_id).ports) {
		if (port.locator_kind != locator_kind || port.locator != locator
			|| !port.requires.empty()) continue;
		if (match) throw std::runtime_error("sequential graph has ambiguous observation port '"
			+domain_id+"."+locator+"'");
		match = &port;
	}
	if (!match) throw std::runtime_error("sequential graph requires an observation port '"
		+domain_id+"."+locator+"'");
	return {domain_id, match->id};
}

inline SequentialOneDThreeDOneDDefinition ResolveSequentialOneDThreeDOneD(
	const MultidomainConfiguration& configuration)
{
	const auto plan = MakeSequentialPressureFlowPlan(configuration.graph,
		configuration.start_domain_id);
	if (plan.domain_ids.size() != 3 || plan.edge_ids.size() != 2)
		throw std::runtime_error("current graph runner requires exactly one 1D--3D--1D chain");
	const auto& upstream = configuration.graph.Domain(plan.domain_ids[0]);
	const auto& three_d = configuration.graph.Domain(plan.domain_ids[1]);
	const auto& downstream = configuration.graph.Domain(plan.domain_ids[2]);
	if (upstream.kind != DomainKind::OneDFlow
		|| three_d.kind != DomainKind::ThreeDBodyFittedFlow
		|| downstream.kind != DomainKind::OneDFlow)
		throw std::runtime_error("current graph runner requires domain kinds 1D--3D--1D");
	const auto& upstream_definition = GraphDomainDefinitionFor(configuration, upstream.id);
	const auto& downstream_definition = GraphDomainDefinitionFor(configuration, downstream.id);
	if (upstream_definition.one_d_inlet_policy != OneDInletPolicy::ConfiguredOpenLoop
		|| downstream_definition.one_d_inlet_policy != OneDInletPolicy::CoupledRoot)
		throw std::runtime_error(
			"current graph runner requires configured-open-loop upstream and coupled-root downstream 1D domains");
	const auto& left_edge = configuration.graph.Edge(plan.edge_ids[0]);
	const auto& right_edge = configuration.graph.Edge(plan.edge_ids[1]);
	SequentialOneDThreeDOneDDefinition result;
	result.upstream_domain_id = upstream.id;
	result.three_d_domain_id = three_d.id;
	result.downstream_domain_id = downstream.id;
	result.upstream_edge_id = left_edge.id;
	result.downstream_edge_id = right_edge.id;
	result.upstream_terminal = EdgeEndpointForDomain(left_edge, upstream.id);
	result.three_d_inlet = EdgeEndpointForDomain(left_edge, three_d.id);
	result.three_d_outlet = EdgeEndpointForDomain(right_edge, three_d.id);
	result.downstream_root = EdgeEndpointForDomain(right_edge, downstream.id);
	const auto& upstream_port = configuration.graph.Port(result.upstream_terminal);
	const auto& downstream_port = configuration.graph.Port(result.downstream_root);
	if (upstream_port.locator == "root"
		|| upstream_port.locator_kind != "runtime_port")
		throw std::runtime_error("upstream coupled port must be a 1D outlet runtime port");
	(void)ParseOneDOutletNodeLocator(upstream_port);
	if (downstream_port.locator_kind != "runtime_port" || downstream_port.locator != "root")
		throw std::runtime_error("downstream coupled port must be the 1D root runtime port");
	result.upstream_root_observation = UniqueObservationPort(configuration.graph,
		upstream.id, "runtime_port", "root");
	result.three_d_wall_observation = UniqueObservationPort(configuration.graph,
		three_d.id, "boundary_label", "0");
	const CouplingPort* downstream_terminal = nullptr;
	for (const auto& port : downstream.ports) {
		if (port.locator_kind != "runtime_port" || port.locator == "root"
			|| !port.requires.empty()) continue;
		(void)ParseOneDOutletNodeLocator(port);
		if (downstream_terminal)
			throw std::runtime_error(
				"current graph runner requires exactly one downstream terminal observation");
		downstream_terminal = &port;
	}
	if (!downstream_terminal)
		throw std::runtime_error(
			"current graph runner requires one downstream terminal observation");
	result.downstream_terminal_observation = {downstream.id, downstream_terminal->id};
	return result;
}

// Retain the public containment helpers while sharing their implementation
// with configuration-time 0D model loading.
inline bool PathIsWithin(const std::filesystem::path& root,
	const std::filesystem::path& candidate)
{
	return multidomain_detail::PathIsWithin(root, candidate);
}

inline std::filesystem::path CanonicalGraphCaseRoot(
	const std::filesystem::path& graph_case_root)
{
	return multidomain_detail::CanonicalGraphCaseRoot(graph_case_root);
}

inline std::filesystem::path ResolveContainedGraphAsset(
	const std::filesystem::path& canonical_root, const std::filesystem::path& relative,
	bool require_directory, const std::string& context)
{
	return multidomain_detail::ResolveContainedGraphAsset(canonical_root, relative,
		require_directory, context);
}

inline std::filesystem::path ResolveContainedCaseFile(
	const std::filesystem::path& canonical_case_directory,
	const std::filesystem::path& relative, const std::string& context)
{
	return multidomain_detail::ResolveContainedCaseFile(canonical_case_directory, relative,
		context);
}

inline std::map<std::string, ResolvedGraphDomainAssets> ResolveGraphDomainAssets(
	const MultidomainConfiguration& configuration, const std::filesystem::path& graph_case_root)
{
	const auto root = multidomain_detail::CanonicalGraphCaseRoot(graph_case_root);
	std::map<std::string, ResolvedGraphDomainAssets> result;
	for (const auto& domain : configuration.domains) {
		ResolvedGraphDomainAssets assets;
		assets.case_directory = multidomain_detail::ResolveContainedGraphAsset(root, domain.case_directory,
			true, "domain '"+domain.id+"' case directory");
		if (domain.kind == DomainKind::ThreeDBodyFittedFlow)
			assets.database = multidomain_detail::ResolveContainedGraphAsset(root, domain.database, false,
				"domain '"+domain.id+"' database");
		result.emplace(domain.id, std::move(assets));
	}
	return result;
}

inline PressureFlowExecutionControls PressureFlowControlsFor(
	const GraphExecutionDefinition& definition)
{
	PressureFlowExecutionControls controls;
	if (definition.kind == GraphExecutionKind::Explicit)
		controls.method = PressureFlowIterationMethod::Explicit;
	else if (definition.kind == GraphExecutionKind::Fixed)
		controls.method = PressureFlowIterationMethod::Fixed;
	else if (definition.kind == GraphExecutionKind::Aitken)
		controls.method = PressureFlowIterationMethod::Aitken;
	else throw std::runtime_error("graph execution has an unknown method");
	controls.maximum_iterations = definition.maximum_iterations;
	controls.pressure_relative_tolerance = definition.pressure_relative_tolerance;
	controls.pressure_reference_pa = definition.pressure_reference_pa;
	controls.flow_relative_tolerance = definition.flow_relative_tolerance;
	controls.relaxation_factor = definition.relaxation_factor;
	controls.minimum_relaxation = definition.minimum_relaxation;
	controls.maximum_relaxation = definition.maximum_relaxation;
	ValidatePressureFlowExecutionControls(controls);
	return controls;
}

inline SpeciesPressureFlowExecutionControls SpeciesPressureFlowControlsFor(
	const GraphExecutionDefinition& definition)
{
	if (!definition.species_routing || definition.species_amount_tolerances.empty())
		throw std::runtime_error(
			"graph execution has no complete schema-v6 species execution controls");
	SpeciesPressureFlowExecutionControls controls;
	controls.hydraulic = PressureFlowControlsFor(definition);
	controls.routing = *definition.species_routing;
	controls.amount_tolerances = definition.species_amount_tolerances;
	controls.Validate();
	return controls;
}

inline SpeciesPressureFlowExecutionControls SpeciesPressureFlowControlsFor(
	const MultidomainConfiguration& configuration)
{
	auto controls = SpeciesPressureFlowControlsFor(configuration.execution);
	if (controls.amount_tolerances.size() != configuration.graph.Species().size())
		throw std::runtime_error("graph execution species amount tolerances have incomplete coverage");
	for (const auto& species : configuration.graph.Species())
		if (!controls.amount_tolerances.count(species.first))
			throw std::runtime_error("graph execution species amount tolerances have incomplete coverage");
	return controls;
}

} // namespace iga

#endif

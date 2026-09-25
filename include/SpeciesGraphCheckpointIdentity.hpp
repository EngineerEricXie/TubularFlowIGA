#ifndef IGA_SPECIES_GRAPH_CHECKPOINT_IDENTITY_HPP
#define IGA_SPECIES_GRAPH_CHECKPOINT_IDENTITY_HPP

#include "CoupledCheckpointManifest.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"

#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>

namespace iga {

namespace species_graph_checkpoint_detail {

inline void Text(Sha256& hash,const std::string& value)
{
	hash.AppendLittleEndian64(value.size());hash.Append(value.data(),value.size());
}

inline void Port(Sha256& hash,const CouplingPort& port)
{
	Text(hash,port.id);Text(hash,port.subsystem_id);
	Text(hash,port.locator_kind);Text(hash,port.locator);
	hash.AppendLittleEndian32(static_cast<std::uint32_t>(port.orientation.native_to_outward_sign));
	hash.AppendLittleEndian64(port.provides.size());
	for(const auto value:port.provides)
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(value));
	hash.AppendLittleEndian64(port.requires.size());
	for(const auto value:port.requires)
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(value));
	hash.AppendLittleEndian64(port.species.size());
	for(const auto& value:port.species)Text(hash,value);
}

} // namespace species_graph_checkpoint_detail

inline std::string BuildSpeciesGraphIdentitySha256(const SimulationGraph& graph)
{
	if(!graph.FsiEdges().empty())
		throw std::invalid_argument("species graph checkpoint does not support FSI edges");
	Sha256 hash;
	using species_graph_checkpoint_detail::Text;
	Text(hash,"SpeciesPressureFlowGraph/v1");
	hash.AppendLittleEndian64(graph.Domains().size());
	for(const auto& entry:graph.Domains()){
		const auto& domain=entry.second;
		if(!domain.surface_interfaces.empty()||!domain.surface_layouts.empty())
			throw std::invalid_argument("species graph checkpoint does not support surfaces");
		Text(hash,entry.first);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(domain.kind));
		hash.AppendLittleEndian64(domain.ports.size());
		for(const auto& port:domain.ports)
			species_graph_checkpoint_detail::Port(hash,port);
		hash.AppendLittleEndian64(domain.species_bindings.size());
		for(const auto& binding:domain.species_bindings){
			Text(hash,binding.first);Text(hash,binding.second);
		}
	}
	hash.AppendLittleEndian64(graph.Edges().size());
	for(const auto& edge:graph.Edges()){
		Text(hash,edge.id);
		Text(hash,edge.first.domain_id);Text(hash,edge.first.port_id);
		Text(hash,edge.second.domain_id);Text(hash,edge.second.port_id);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(edge.law));
		hash.AppendLittleEndian64(edge.species.size());
		for(const auto& species:edge.species)Text(hash,species);
	}
	hash.AppendLittleEndian64(graph.Species().size());
	for(const auto& species:graph.Species()){
		Text(hash,species.first);Text(hash,species.second.concentration_unit);
	}
	return hash.Hex();
}

inline CoupledCheckpointCompatibility BuildSpeciesGraphCheckpointCompatibility(
	const SimulationGraph& graph,
	const std::map<std::string,std::string>& domain_models,
	const SpeciesPressureFlowExecutionControls& controls,
	std::uint32_t ranks)
{
	controls.Validate();
	if(domain_models.size()!=graph.Domains().size())
		throw std::invalid_argument("species graph checkpoint domain model catalog differs");
	CoupledCheckpointCompatibility result;
	result.case_sha256=BuildSpeciesGraphIdentitySha256(graph);
	Sha256 configuration;
	species_graph_checkpoint_detail::Text(configuration,"SpeciesGraphModels/v1");
	species_graph_checkpoint_detail::Text(configuration,result.case_sha256);
	for(const auto& domain:graph.Domains()){
		const auto found=domain_models.find(domain.first);
		if(found==domain_models.end()||!coupled_checkpoint_detail::Digest(found->second))
			throw std::invalid_argument("species graph checkpoint model identity is missing");
		species_graph_checkpoint_detail::Text(configuration,domain.first);
		species_graph_checkpoint_detail::Text(configuration,found->second);
	}
	result.configuration_sha256=configuration.Hex();
	Sha256 execution;
	species_graph_checkpoint_detail::Text(execution,"SpeciesGraphExecution/v1");
	const auto& hydraulic=controls.hydraulic;
	execution.AppendLittleEndian32(static_cast<std::uint32_t>(hydraulic.method));
	execution.AppendLittleEndian32(static_cast<std::uint32_t>(hydraulic.maximum_iterations));
	for(const double value:{hydraulic.pressure_relative_tolerance,
		hydraulic.pressure_reference_pa,hydraulic.flow_relative_tolerance,
		hydraulic.relaxation_factor,hydraulic.minimum_relaxation,
		hydraulic.maximum_relaxation,controls.routing.flow_switch_m3_s,
		controls.routing.flow_absolute_tolerance_m3_s,
		controls.routing.flow_relative_tolerance})
		execution.AppendNormalizedDouble(value);
	execution.AppendLittleEndian64(controls.amount_tolerances.size());
	for(const auto& entry:controls.amount_tolerances){
		species_graph_checkpoint_detail::Text(execution,entry.first);
		execution.AppendNormalizedDouble(entry.second.absolute_tolerance);
		execution.AppendNormalizedDouble(entry.second.reference_amount);
		execution.AppendNormalizedDouble(entry.second.relative_tolerance);
	}
	result.execution_sha256=execution.Hex();
	result.ranks=ranks;
	coupled_checkpoint_detail::Validate(result);
	return result;
}

} // namespace iga

#endif

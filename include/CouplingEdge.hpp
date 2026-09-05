#ifndef IGA_COUPLING_EDGE_HPP
#define IGA_COUPLING_EDGE_HPP

#include "CouplingPort.hpp"

#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace iga {

enum class DomainKind {
	OneDFlow,
	ThreeDBodyFittedFlow,
	ThreeDImmersedFlow,
	SurfaceMembraneStructure,
	// Appended to preserve the serialized identities of the pre-existing kinds.
	ZeroDFlow
};

inline const char* DomainKindName(DomainKind kind)
{
	if (kind == DomainKind::OneDFlow) return "one_d_flow";
	if (kind == DomainKind::ThreeDBodyFittedFlow) return "three_d_body_fitted_flow";
	if (kind == DomainKind::ThreeDImmersedFlow) return "three_d_immersed_flow";
	if (kind == DomainKind::SurfaceMembraneStructure) return "surface_membrane_structure";
	if (kind == DomainKind::ZeroDFlow) return "zero_d_flow";
	return "unknown";
}

inline bool IsKnownDomainKind(DomainKind kind)
{
	return kind == DomainKind::OneDFlow || kind == DomainKind::ThreeDBodyFittedFlow
		|| kind == DomainKind::ThreeDImmersedFlow
		|| kind == DomainKind::SurfaceMembraneStructure || kind == DomainKind::ZeroDFlow;
}

inline bool IsFlowDomainKind(DomainKind kind)
{
	return kind == DomainKind::OneDFlow || kind == DomainKind::ThreeDBodyFittedFlow
		|| kind == DomainKind::ThreeDImmersedFlow || kind == DomainKind::ZeroDFlow;
}

// The membrane is a two-dimensional material topology.  Its embedding is
// three-dimensional, but it must never be treated as a volumetric flow domain.
inline int DomainTopologyDimensionOf(DomainKind kind)
{
	if (kind == DomainKind::OneDFlow) return 1;
	if (kind == DomainKind::ZeroDFlow) return 0;
	if (kind == DomainKind::ThreeDBodyFittedFlow
		|| kind == DomainKind::ThreeDImmersedFlow) return 3;
	if (kind == DomainKind::SurfaceMembraneStructure) return 2;
	throw std::runtime_error("domain kind has no topology dimension");
}

inline int DomainEmbeddingDimensionOf(DomainKind kind)
{
	if (kind == DomainKind::SurfaceMembraneStructure) return 3;
	return DomainTopologyDimensionOf(kind);
}

// Retained for existing scalar-coupling callers.  New code must select either
// topology or embedding dimension deliberately.
inline int DomainDimensionOf(DomainKind kind)
{
	return DomainTopologyDimensionOf(kind);
}

struct PortRef {
	std::string domain_id;
	std::string port_id;
};

inline bool operator==(const PortRef& first, const PortRef& second)
{
	return first.domain_id == second.domain_id && first.port_id == second.port_id;
}

inline bool operator<(const PortRef& first, const PortRef& second)
{
	return std::tie(first.domain_id, first.port_id) < std::tie(second.domain_id, second.port_id);
}

inline void ValidatePortRef(const PortRef& reference)
{
	if (reference.domain_id.empty())
		throw std::runtime_error("coupling endpoint domain_id must be nonempty");
	if (reference.port_id.empty())
		throw std::runtime_error("coupling endpoint port_id must be nonempty");
}

enum class CouplingLaw {
	PressureFlow
};

inline const char* CouplingLawName(CouplingLaw law)
{
	if (law == CouplingLaw::PressureFlow) return "pressure_flow";
	return "unknown";
}

struct CouplingEdge {
	std::string id;
	PortRef first;
	PortRef second;
	CouplingLaw law = CouplingLaw::PressureFlow;
	std::set<std::string> species;

	CouplingEdge() = default;
	CouplingEdge(std::string edge_id, PortRef first_endpoint, PortRef second_endpoint,
		CouplingLaw edge_law, std::set<std::string> coupled_species = {})
		: id(std::move(edge_id)), first(std::move(first_endpoint)),
		  second(std::move(second_endpoint)), law(edge_law),
		  species(std::move(coupled_species)) {}
};

inline void ValidateCouplingEdge(const CouplingEdge& edge)
{
	if (edge.id.empty()) throw std::runtime_error("coupling edge id must be nonempty");
	ValidatePortRef(edge.first);
	ValidatePortRef(edge.second);
	if (edge.first == edge.second || edge.first.domain_id == edge.second.domain_id)
		throw std::runtime_error("coupling edge endpoints must belong to distinct domains");
	if (edge.law != CouplingLaw::PressureFlow)
		throw std::runtime_error("coupling edge has an unsupported law");
}

inline void ValidatePressureFlowEdgePorts(const CouplingPort& first,
	const CouplingPort& second)
{
	for (const auto* port : {&first, &second}) {
		if (!port->provides.count(PortQuantity::FlowRate)
			|| !port->provides.count(PortQuantity::MeanPressure))
			throw std::runtime_error(
				"pressure_flow endpoints must report flow_rate and mean_pressure");
		for (const auto quantity : port->requires)
			if (quantity != PortQuantity::FlowRate && quantity != PortQuantity::MeanPressure
				&& quantity != PortQuantity::SpeciesConcentration
				&& quantity != PortQuantity::SpeciesFlux)
				throw std::runtime_error(
					"pressure_flow endpoints declare an unsupported input quantity");
	}
	const bool first_accepts_flow = first.requires.count(PortQuantity::FlowRate) != 0;
	const bool second_accepts_flow = second.requires.count(PortQuantity::FlowRate) != 0;
	const bool first_accepts_pressure = first.requires.count(PortQuantity::MeanPressure) != 0;
	const bool second_accepts_pressure = second.requires.count(PortQuantity::MeanPressure) != 0;
	if (first_accepts_flow == second_accepts_flow
		|| first_accepts_pressure == second_accepts_pressure
		|| first_accepts_flow == first_accepts_pressure)
		throw std::runtime_error(
			"pressure_flow requires one flow receiver and the opposite pressure receiver");
	if ((first_accepts_flow && !second.provides.count(PortQuantity::FlowRate))
		|| (second_accepts_flow && !first.provides.count(PortQuantity::FlowRate))
		|| (first_accepts_pressure && !second.provides.count(PortQuantity::MeanPressure))
		|| (second_accepts_pressure && !first.provides.count(PortQuantity::MeanPressure)))
		throw std::runtime_error("pressure_flow endpoint does not provide its peer's input");
}

inline void ValidatePressureFlowSpeciesPorts(const CouplingEdge& edge,
	const CouplingPort& first, const CouplingPort& second)
{
	if (edge.species.empty()) {
		if (!first.species.empty() || !second.species.empty())
			throw std::runtime_error("uncoupled species cannot be declared on an edge endpoint");
		return;
	}
	if (first.species != edge.species || second.species != edge.species)
		throw std::runtime_error("coupled species sets must match both edge endpoints");
	const std::set<PortQuantity> capabilities{PortQuantity::SpeciesConcentration,
		PortQuantity::SpeciesFlux};
	for (const auto* port : {&first, &second})
		for (const auto quantity : capabilities)
			if (!port->provides.count(quantity) || !port->requires.count(quantity))
				throw std::runtime_error(
					"species edge endpoints must provide and accept concentration and total flux");
}

} // namespace iga

#endif

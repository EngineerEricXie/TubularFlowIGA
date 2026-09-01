#ifndef IGA_COUPLING_EDGE_HPP
#define IGA_COUPLING_EDGE_HPP

#include "CouplingPort.hpp"

#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

namespace iga {

enum class DomainKind {
	OneDFlow,
	ThreeDBodyFittedFlow
};

inline const char* DomainKindName(DomainKind kind)
{
	if (kind == DomainKind::OneDFlow) return "one_d_flow";
	if (kind == DomainKind::ThreeDBodyFittedFlow) return "three_d_body_fitted_flow";
	return "unknown";
}

inline bool IsKnownDomainKind(DomainKind kind)
{
	return kind == DomainKind::OneDFlow || kind == DomainKind::ThreeDBodyFittedFlow;
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
			if (quantity != PortQuantity::FlowRate && quantity != PortQuantity::MeanPressure)
				throw std::runtime_error(
					"pressure_flow endpoints may accept only flow_rate or mean_pressure");
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

} // namespace iga

#endif

#ifndef IGA_FLOW_DOMAIN_PORT_METADATA_HPP
#define IGA_FLOW_DOMAIN_PORT_METADATA_HPP

#include "CouplingPort.hpp"

#include <charconv>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class OneDInletPolicy { ConfiguredOpenLoop, CoupledRoot };

inline int ParseOneDOutletNodeLocator(const CouplingPort& port)
{
	const std::string prefix = "outlet:";
	if (port.locator.compare(0, prefix.size(), prefix) != 0)
		throw std::runtime_error("1D runtime port must be root or outlet:<node-id>");
	int node = 0;
	const auto first = port.locator.data()+prefix.size();
	const auto parsed = std::from_chars(first, port.locator.data()+port.locator.size(), node);
	if (first == port.locator.data()+port.locator.size() || parsed.ec != std::errc{}
		|| parsed.ptr != port.locator.data()+port.locator.size() || node <= 0)
		throw std::runtime_error(
			"1D outlet runtime-port locator must end in a positive integer node id");
	return node;
}

inline void ValidateOneDFlowDomainMetadata(const std::string& domain_id,
	const std::vector<CouplingPort>& ports, OneDInletPolicy inlet_policy)
{
	if (domain_id.empty()) throw std::runtime_error("1D flow domain id must be nonempty");
	if (inlet_policy != OneDInletPolicy::ConfiguredOpenLoop
		&& inlet_policy != OneDInletPolicy::CoupledRoot)
		throw std::runtime_error("1D flow domain has an unknown inlet policy");
	ValidateCouplingPorts(ports);
	const std::set<PortQuantity> supported_provides = {PortQuantity::Area,
		PortQuantity::FlowRate, PortQuantity::MeanPressure,
		PortQuantity::SpeciesConcentration, PortQuantity::SpeciesFlux};
	int root_flow_receivers = 0;
	for (const auto& port : ports) {
		if (port.subsystem_id != domain_id || port.locator_kind != "runtime_port")
			throw std::runtime_error("1D flow domain requires matching runtime_port metadata");
		if (port.orientation.native_to_outward_sign != 1)
			throw std::runtime_error(
				"1D runtime_port metadata must use canonical outward orientation +1");
		for (const auto quantity : port.provides)
			if (!supported_provides.count(quantity))
				throw std::runtime_error("1D flow runtime port declares an unsupported output quantity");
		if (port.locator == "root") {
			if (!port.requires.empty()
				&& port.requires != std::set<PortQuantity>{PortQuantity::FlowRate})
				throw std::runtime_error("1D root runtime port may receive only flow_rate");
			if (port.requires.count(PortQuantity::FlowRate)) ++root_flow_receivers;
		} else {
			(void)ParseOneDOutletNodeLocator(port);
			if (!port.requires.empty()
				&& port.requires != std::set<PortQuantity>{PortQuantity::MeanPressure})
				throw std::runtime_error("1D outlet runtime port may receive only mean_pressure");
		}
	}
	if (inlet_policy == OneDInletPolicy::ConfiguredOpenLoop && root_flow_receivers != 0)
		throw std::runtime_error(
			"configured-open-loop 1D domain cannot expose a coupled root flow receiver");
	if (inlet_policy == OneDInletPolicy::CoupledRoot && root_flow_receivers != 1)
		throw std::runtime_error(
			"coupled-root 1D domain requires exactly one root flow receiver");
}

inline int ParseThreeDFlowBoundaryLabel(const CouplingPort& port)
{
	ValidateCouplingPort(port);
	if (port.locator_kind != "boundary_label")
		throw std::runtime_error("3D flow port requires a boundary_label locator");
	int label = -1;
	const auto first = port.locator.data();
	const auto parsed = std::from_chars(first, first+port.locator.size(), label);
	if (first == first+port.locator.size() || parsed.ec != std::errc{}
		|| parsed.ptr != first+port.locator.size() || label < 0)
		throw std::runtime_error(
			"3D flow boundary_label locator must be a nonnegative integer");
	return label;
}

inline void ValidateThreeDBodyFittedFlowDomainMetadata(const std::string& domain_id,
	const std::vector<CouplingPort>& ports)
{
	if (domain_id.empty()) throw std::runtime_error("3D flow domain id must be nonempty");
	ValidateCouplingPorts(ports);
	const std::set<PortQuantity> supported_provides = {PortQuantity::Area,
		PortQuantity::FlowRate, PortQuantity::MeanPressure};
	const std::set<PortQuantity> supported_requires = {PortQuantity::FlowRate,
		PortQuantity::MeanPressure, PortQuantity::MeanNormalTraction};
	for (const auto& port : ports) {
		if (port.subsystem_id != domain_id)
			throw std::runtime_error("3D flow domain port metadata has a mismatched domain id");
		(void)ParseThreeDFlowBoundaryLabel(port);
		if (port.requires.size() > 1)
			throw std::runtime_error("3D flow port may receive at most one input quantity");
		for (const auto quantity : port.provides)
			if (!supported_provides.count(quantity))
				throw std::runtime_error("3D flow port declares an unsupported output quantity");
		for (const auto quantity : port.requires)
			if (!supported_requires.count(quantity))
				throw std::runtime_error("3D flow port declares an unsupported input quantity");
	}
}

} // namespace iga

#endif

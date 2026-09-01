#ifndef IGA_COUPLING_PORT_HPP
#define IGA_COUPLING_PORT_HPP

#include "SpeciesCoupling.hpp"

#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// All flow and species-flux signs below are positive outward from a subsystem.
enum class PortQuantity {
	Area,
	FlowRate,
	MeanPressure,
	MeanNormalTraction,
	TotalPressure,
	SpeciesConcentration,
	SpeciesFlux
};

inline const char* PortQuantityName(PortQuantity quantity)
{
	switch (quantity) {
	case PortQuantity::Area: return "area";
	case PortQuantity::FlowRate: return "flow_rate";
	case PortQuantity::MeanPressure: return "mean_pressure";
	case PortQuantity::MeanNormalTraction: return "mean_normal_traction";
	case PortQuantity::TotalPressure: return "total_pressure";
	case PortQuantity::SpeciesConcentration: return "species_concentration";
	case PortQuantity::SpeciesFlux: return "species_flux";
	}
	return "unknown";
}

inline bool IsKnownPortQuantity(PortQuantity quantity)
{
	return quantity == PortQuantity::Area || quantity == PortQuantity::FlowRate
		|| quantity == PortQuantity::MeanPressure
		|| quantity == PortQuantity::MeanNormalTraction
		|| quantity == PortQuantity::TotalPressure
		|| quantity == PortQuantity::SpeciesConcentration
		|| quantity == PortQuantity::SpeciesFlux;
}

inline void RequireFinitePortValue(const std::string& name, double value)
{
	if (!std::isfinite(value)) throw std::runtime_error(name+" must be finite");
}

inline void ValidatePortSpeciesValues(const std::map<std::string, double>& values,
	const std::string& name)
{
	for (const auto& value : values) {
		if (value.first.empty()) throw std::runtime_error(name+" species name must be nonempty");
		RequireFinitePortValue(name+" species "+value.first, value.second);
	}
}

struct PortOrientation {
	// Converts a backend-native positive direction to positive subsystem outflow.
	int native_to_outward_sign = 1;

	void Validate() const
	{
		if (native_to_outward_sign != 1 && native_to_outward_sign != -1)
			throw std::runtime_error("port native_to_outward_sign must be +1 or -1");
	}

	double ToOutward(double native_value) const
	{
		Validate();
		RequireFinitePortValue("native port value", native_value);
		return native_to_outward_sign*native_value;
	}

	double ToNative(double outward_value) const
	{
		Validate();
		RequireFinitePortValue("outward port value", outward_value);
		return native_to_outward_sign*outward_value;
	}
};

struct PortState {
	double time_s = 0.0;
	std::optional<double> area_m2;
	std::optional<double> outward_flow_m3_s;
	std::optional<double> mean_pressure_pa;
	std::optional<double> mean_normal_traction_pa;
	std::optional<double> total_pressure_pa;
	std::map<std::string, double> concentration;
	std::map<std::string, double> outward_species_flux;
};

struct PortBoundaryData {
	double time_s = 0.0;
	std::optional<double> outward_flow_m3_s;
	std::optional<double> mean_pressure_pa;
	std::optional<double> mean_normal_traction_pa;
	std::optional<double> total_pressure_pa;
	std::map<std::string, double> concentration;
	std::map<std::string, double> outward_species_flux;
};

struct CouplingPort {
	std::string id;
	std::string subsystem_id;
	std::string locator_kind;
	std::string locator;
	PortOrientation orientation;
	std::set<PortQuantity> provides;
	std::set<PortQuantity> requires;
	std::set<std::string> species;
};

struct CouplingResidual {
	double outward_flow_m3_s = 0.0;
	std::map<std::string, double> outward_species_flux;
};

inline void ValidatePortState(const PortState& state)
{
	RequireFinitePortValue("port state time_s", state.time_s);
	if (state.area_m2) {
		RequireFinitePortValue("port state area_m2", *state.area_m2);
		if (!(*state.area_m2 > 0.0))
			throw std::runtime_error("port state area_m2 must be positive when provided");
	}
	if (state.outward_flow_m3_s)
		RequireFinitePortValue("port state outward_flow_m3_s", *state.outward_flow_m3_s);
	if (state.mean_pressure_pa)
		RequireFinitePortValue("port state mean_pressure_pa", *state.mean_pressure_pa);
	if (state.mean_normal_traction_pa)
		RequireFinitePortValue("port state mean_normal_traction_pa",
			*state.mean_normal_traction_pa);
	if (state.total_pressure_pa)
		RequireFinitePortValue("port state total_pressure_pa", *state.total_pressure_pa);
	ValidatePortSpeciesValues(state.concentration, "port state concentration");
	ValidatePortSpeciesValues(state.outward_species_flux, "port state outward_species_flux");
}

inline void ValidatePortBoundaryData(const PortBoundaryData& data)
{
	RequireFinitePortValue("port boundary time_s", data.time_s);
	if (data.outward_flow_m3_s)
		RequireFinitePortValue("port boundary outward_flow_m3_s", *data.outward_flow_m3_s);
	if (data.mean_pressure_pa)
		RequireFinitePortValue("port boundary mean_pressure_pa", *data.mean_pressure_pa);
	if (data.mean_normal_traction_pa)
		RequireFinitePortValue("port boundary mean_normal_traction_pa",
			*data.mean_normal_traction_pa);
	if (data.total_pressure_pa)
		RequireFinitePortValue("port boundary total_pressure_pa", *data.total_pressure_pa);
	ValidatePortSpeciesValues(data.concentration, "port boundary concentration");
	ValidatePortSpeciesValues(data.outward_species_flux,
		"port boundary outward_species_flux");
}

inline void ValidatePortQuantityDeclarations(const std::vector<PortQuantity>& provides,
	const std::vector<PortQuantity>& requires)
{
	std::set<PortQuantity> provided;
	for (const auto quantity : provides) {
		if (!IsKnownPortQuantity(quantity))
			throw std::runtime_error("port provides an unknown quantity");
		if (!provided.insert(quantity).second)
			throw std::runtime_error(std::string("port provides duplicate quantity ")
				+PortQuantityName(quantity));
	}
	std::set<PortQuantity> required;
	for (const auto quantity : requires) {
		if (!IsKnownPortQuantity(quantity))
			throw std::runtime_error("port requires an unknown quantity");
		if (!required.insert(quantity).second)
			throw std::runtime_error(std::string("port requires duplicate quantity ")
				+PortQuantityName(quantity));
	}
	if (provided.empty() && required.empty())
		throw std::runtime_error("port must provide or require at least one quantity");
}

inline void ValidateCouplingPort(const CouplingPort& port)
{
	if (port.id.empty()) throw std::runtime_error("coupling port id must be nonempty");
	if (port.subsystem_id.empty())
		throw std::runtime_error("coupling port subsystem_id must be nonempty");
	if (port.locator_kind.empty())
		throw std::runtime_error("coupling port locator_kind must be nonempty");
	if (port.locator.empty()) throw std::runtime_error("coupling port locator must be nonempty");
	port.orientation.Validate();
	ValidatePortQuantityDeclarations(
		std::vector<PortQuantity>(port.provides.begin(), port.provides.end()),
		std::vector<PortQuantity>(port.requires.begin(), port.requires.end()));
	for (const auto& species : port.species)
		if (species.empty()) throw std::runtime_error("port species id must be nonempty");
	const bool has_species_quantity
		= port.provides.count(PortQuantity::SpeciesConcentration)
		|| port.provides.count(PortQuantity::SpeciesFlux)
		|| port.requires.count(PortQuantity::SpeciesConcentration)
		|| port.requires.count(PortQuantity::SpeciesFlux);
	if (!port.species.empty() && !has_species_quantity)
		throw std::runtime_error(
			"port species declarations require species quantities");
}

inline void ValidateCouplingPorts(const std::vector<CouplingPort>& ports)
{
	std::set<std::pair<std::string, std::string>> ids;
	for (const auto& port : ports) {
		ValidateCouplingPort(port);
		if (!ids.emplace(port.subsystem_id, port.id).second)
			throw std::runtime_error(
				"coupling port ids must be unique within each subsystem");
	}
}

inline CouplingResidual ConservativeEdgeResidual(const PortState& first,
	const PortState& second)
{
	ValidatePortState(first);
	ValidatePortState(second);
	if (!first.outward_flow_m3_s || !second.outward_flow_m3_s)
		throw std::runtime_error("conservative edge residual requires outward flow at both ports");
	CouplingResidual result;
	result.outward_flow_m3_s = *first.outward_flow_m3_s+*second.outward_flow_m3_s;
	if (first.outward_species_flux.size() != second.outward_species_flux.size())
		throw std::runtime_error("conservative edge residual requires matching species sets");
	for (const auto& species : first.outward_species_flux) {
		const auto other = second.outward_species_flux.find(species.first);
		if (other == second.outward_species_flux.end())
			throw std::runtime_error("conservative edge residual requires matching species sets");
		result.outward_species_flux.emplace(species.first, species.second+other->second);
	}
	return result;
}

inline CouplingResidual ConservativeEdgeResidual(const PortState& first,
	const PortState& second, const std::set<std::string>& species)
{
	const auto result = ConservativeEdgeResidual(first, second);
	std::set<std::string> actual;
	for (const auto& value : result.outward_species_flux) actual.insert(value.first);
	if (actual != species)
		throw std::runtime_error(
			"conservative edge residual does not match its declared species set");
	return result;
}

} // namespace iga

#endif

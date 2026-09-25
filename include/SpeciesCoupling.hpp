#ifndef IGA_SPECIES_COUPLING_HPP
#define IGA_SPECIES_COUPLING_HPP

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct SpeciesDefinition {
	std::string id;
	std::string concentration_unit;
};

inline bool operator==(const SpeciesDefinition& first, const SpeciesDefinition& second)
{
	return first.id == second.id
		&& first.concentration_unit == second.concentration_unit;
}

inline std::string SpeciesFluxUnit(const SpeciesDefinition& species)
{
	if (species.concentration_unit.empty())
		throw std::runtime_error("species concentration unit must be nonempty");
	return "("+species.concentration_unit+")*m^3/s";
}

inline std::map<std::string, SpeciesDefinition> MakeSpeciesRegistry(
	const std::vector<SpeciesDefinition>& definitions)
{
	std::map<std::string, SpeciesDefinition> result;
	for (const auto& definition : definitions) {
		if (definition.id.empty() || definition.concentration_unit.empty())
			throw std::runtime_error("species id and concentration unit must be nonempty");
		if (!result.emplace(definition.id, definition).second)
			throw std::runtime_error("species ids must be unique");
	}
	return result;
}

struct SpeciesRoutingControls {
	double flow_switch_m3_s = 0.0;
	double flow_absolute_tolerance_m3_s = 0.0;
	double flow_relative_tolerance = 1.0e-10;

	void Validate() const
	{
		if (flow_switch_m3_s < 0.0 || !std::isfinite(flow_switch_m3_s)
			|| flow_absolute_tolerance_m3_s < 0.0
			|| !std::isfinite(flow_absolute_tolerance_m3_s)
			|| !(flow_relative_tolerance > 0.0)
			|| !std::isfinite(flow_relative_tolerance))
			throw std::runtime_error("species routing controls are invalid");
	}
};

struct SpeciesAmountTolerance {
	double absolute_tolerance = 0.0;
	double reference_amount = 0.0;
	double relative_tolerance = 1.0e-10;

	void Validate() const
	{
		if (!(absolute_tolerance >= 0.0) || !std::isfinite(absolute_tolerance)
			|| !(reference_amount >= 0.0) || !std::isfinite(reference_amount)
			|| !(relative_tolerance > 0.0) || !std::isfinite(relative_tolerance))
			throw std::runtime_error("species amount tolerance is invalid");
	}
};

enum class SpeciesDonor { First, Second };

inline SpeciesDonor ResolveSpeciesDonor(double first_outward_flow_m3_s,
	double second_outward_flow_m3_s, const SpeciesRoutingControls& controls,
	std::optional<SpeciesDonor> last_committed_donor = std::nullopt)
{
	controls.Validate();
	if (!std::isfinite(first_outward_flow_m3_s)
		|| !std::isfinite(second_outward_flow_m3_s))
		throw std::runtime_error("species routing requires finite outward flows");
	const double scale = std::max(std::abs(first_outward_flow_m3_s),
		std::abs(second_outward_flow_m3_s));
	const double tolerance = std::max(controls.flow_absolute_tolerance_m3_s,
		controls.flow_relative_tolerance*scale);
	if (std::abs(first_outward_flow_m3_s+second_outward_flow_m3_s) > tolerance) {
		std::ostringstream message;
		message.precision(17);
		message << "species routing requires a conservative flow pair: first="
			<< first_outward_flow_m3_s << " second=" << second_outward_flow_m3_s
			<< " tolerance=" << tolerance;
		throw std::runtime_error(message.str());
	}
	if (std::abs(first_outward_flow_m3_s) <= controls.flow_switch_m3_s
		&& std::abs(second_outward_flow_m3_s) <= controls.flow_switch_m3_s) {
		if (!last_committed_donor)
			throw std::runtime_error("near-zero species routing requires committed ownership");
		return *last_committed_donor;
	}
	if (first_outward_flow_m3_s > controls.flow_switch_m3_s
		&& second_outward_flow_m3_s < -controls.flow_switch_m3_s)
		return SpeciesDonor::First;
	if (second_outward_flow_m3_s > controls.flow_switch_m3_s
		&& first_outward_flow_m3_s < -controls.flow_switch_m3_s)
		return SpeciesDonor::Second;
	throw std::runtime_error("species routing requires one outward donor and one inward receiver");
}

} // namespace iga

#endif

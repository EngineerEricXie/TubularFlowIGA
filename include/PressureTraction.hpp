#ifndef IGA_PRESSURE_TRACTION_HPP
#define IGA_PRESSURE_TRACTION_HPP

#include "SimulationConfig.hpp"

#include <map>
#include <stdexcept>

namespace iga {

inline std::map<int, double> ExtractPressureTractions(
	const SimulationConfiguration& configuration,
	const EquationSystemDefinition& system)
{
	if (system.kind != EquationKind::NavierStokes || system.unknowns.size() != 2)
		throw std::runtime_error("pressure traction requires a Navier-Stokes system");
	const auto& pressure_name = system.unknowns[1];
	std::map<int, double> result;
	for (const auto& boundary : configuration.boundaries)
		for (const auto& condition : boundary.conditions) {
			if (condition.kind != FieldBoundaryKind::PressureTraction) continue;
			if (condition.field != pressure_name || condition.value.size() != 1)
				throw std::runtime_error(
					"pressure_traction must provide one value for the pressure field");
			if (!result.emplace(boundary.label, condition.value[0]).second)
				throw std::runtime_error("duplicate pressure traction boundary label");
		}
	return result;
}

} // namespace iga

#endif

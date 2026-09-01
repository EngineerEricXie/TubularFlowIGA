#ifndef IGA_THREE_D_FLOW_COUPLING_HPP
#define IGA_THREE_D_FLOW_COUPLING_HPP

#include "CouplingPort.hpp"
#include "SimulationConfig.hpp"

#include <charconv>
#include <cmath>
#include <stdexcept>
#include <string>

namespace iga {

inline int ParseThreeDFlowBoundaryLabel(const CouplingPort& port)
{
	ValidateCouplingPort(port);
	if (port.locator_kind != "boundary_label")
		throw std::runtime_error("3D reference profile input requires a boundary_label locator");
	int label = -1;
	const auto first = port.locator.data();
	const auto parsed = std::from_chars(first, first+port.locator.size(), label);
	if (first == first+port.locator.size() || parsed.ec != std::errc{}
		|| parsed.ptr != first+port.locator.size() || label < 0)
		throw std::runtime_error("3D reference profile boundary_label locator must be a nonnegative integer");
	return label;
}

// Applies a scalar outward-positive flow to an explicitly configured velocity
// reference profile. The caller supplies the reference profile's outward flow.
inline void ApplyThreeDReferenceProfileInput(SimulationConfiguration& configuration,
	const EquationSystemDefinition& system, const CouplingPort& port,
	const PortBoundaryData& input, double reference_outward_flow_m3_s)
{
	ValidatePortBoundaryData(input);
	const int label = ParseThreeDFlowBoundaryLabel(port);
	if (!port.requires.count(PortQuantity::FlowRate))
		throw std::runtime_error("3D reference profile port must require flow_rate");
	if (!input.outward_flow_m3_s || input.mean_pressure_pa || input.mean_normal_traction_pa
		|| input.total_pressure_pa || !input.concentration.empty()
		|| !input.outward_species_flux.empty())
		throw std::runtime_error(
			"3D reference profile input requires exactly outward_flow_m3_s");
	if (!std::isfinite(reference_outward_flow_m3_s) || reference_outward_flow_m3_s == 0.0)
		throw std::runtime_error("3D reference profile outward flow must be finite and nonzero");
	if (system.kind != EquationKind::NavierStokes || system.unknowns.empty())
		throw std::runtime_error("3D reference profile input requires a Navier-Stokes system");
	const auto& velocity = system.unknowns.front();
	FieldBoundaryCondition* profile = nullptr;
	for (auto& boundary : configuration.boundaries) {
		if (boundary.label != label) continue;
		for (auto& condition : boundary.conditions) {
			if (condition.field != velocity) continue;
			if (profile)
				throw std::runtime_error("3D reference profile input has ambiguous velocity boundary declarations");
			profile = &condition;
		}
	}
	if (!profile)
		throw std::runtime_error("3D reference profile port label has no velocity boundary declaration");
	if (profile->kind != FieldBoundaryKind::Dirichlet
		|| profile->profile != "initial_velocityfield.txt" || !profile->waveform.empty())
		throw std::runtime_error(
			"3D reference profile input requires a materialized initial_velocityfield.txt velocity Dirichlet profile");
	const double scale = *input.outward_flow_m3_s/reference_outward_flow_m3_s;
	if (!std::isfinite(scale) || (*input.outward_flow_m3_s != 0.0 && scale == 0.0))
		throw std::runtime_error("3D reference profile scale is not representable");
	profile->scale = scale;
}

} // namespace iga

#endif

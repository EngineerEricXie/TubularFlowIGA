#ifndef IGA_THREE_D_VCA_COUPLING_HPP
#define IGA_THREE_D_VCA_COUPLING_HPP

#include "SimulationConfig.hpp"
#include "VascularCoupling.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

inline void RequireThreeDVascularPorts(const CouplingDefinition& definition,
	const std::string& backend)
{
	const auto& ports = definition.three_d_ports;
	if (ports.inlet_label < 0 || ports.outlet_labels.empty())
		throw std::runtime_error(backend+" requires coupling.three_d_ports with an inlet_label and outlet_labels");
	for (const int label : ports.outlet_labels)
		if (label < 0 || label == ports.inlet_label)
			throw std::runtime_error(backend+" has invalid coupling.three_d_ports labels");
}

inline void RequireThreeDFlowOnlyCircuit(const CouplingDefinition& definition)
{
	if (definition.mode != SimulationScopeMode::VcaClosedLoop)
		throw std::runtime_error("CPU 3D VCA bridge currently supports vca_closed_loop only");
	const auto& circuit = definition.external_circuit;
	if (!circuit.reservoir.species.empty() || circuit.oxygenator.enabled
		|| circuit.dialyzer.enabled || !circuit.infusion_rates.empty())
		throw std::runtime_error("CPU 3D VCA flow bridge does not yet transport species; reservoir species and external solute devices require iga_1d");
}

inline CompiledLinearSystem RequireThreeDVcaTransportSystem(
	const SimulationConfiguration& configuration)
{
	std::string name;
	for (const auto& system : configuration.equation_systems)
		if (system.kind == EquationKind::LinearTransport) {
			if (!name.empty())
				throw std::runtime_error("CPU 3D VCA species bridge requires exactly one linear_transport system");
			name = system.name;
		}
	if (name.empty())
		throw std::runtime_error("CPU 3D VCA species bridge requires a linear_transport system");
	const auto result = CompileLinearSystem(configuration, name);
	if (result.dt != configuration.time.dt || result.steps != configuration.time.steps)
		throw std::runtime_error("CPU 3D VCA flow and transport must have matching time integration");
	return result;
}

inline void ApplyThreeDVascularSpeciesInlet(SimulationConfiguration& configuration,
	const CompiledLinearSystem& system, const VascularInletState& inlet)
{
	ValidateVascularInletState(inlet);
	RequireThreeDVascularPorts(configuration.coupling, "CPU 3D VCA species bridge");
	for (const auto& species : inlet.species) {
		const auto field = system.field_index.find(species.first);
		if (field == system.field_index.end())
			throw std::runtime_error("CPU 3D VCA inlet species '"+species.first
				+"' is not transported by the configured system");
		bool applied = false;
		for (auto& boundary : configuration.boundaries) {
			if (boundary.label != configuration.coupling.three_d_ports.inlet_label)
				continue;
			for (auto& condition : boundary.conditions) {
				if (condition.field != species.first) continue;
				if (condition.kind != FieldBoundaryKind::Dirichlet)
					throw std::runtime_error("CPU 3D VCA inlet species requires a Dirichlet boundary");
				condition.value = {species.second};
				condition.waveform.clear();
				applied = true;
			}
		}
		if (!applied)
			throw std::runtime_error("CPU 3D VCA inlet label has no Dirichlet boundary for species '"
				+species.first+"'");
	}
}

inline void ApplyThreeDVascularInlet(SimulationConfiguration& configuration,
	const EquationSystemDefinition& system, const VascularInletState& inlet,
	double reference_inlet_flow_m3_s)
{
	ValidateVascularInletState(inlet);
	RequireThreeDVascularPorts(configuration.coupling, "CPU 3D VCA bridge");
	if (!inlet.has_flow || inlet.has_pressure)
		throw std::runtime_error("CPU 3D VCA bridge currently requires a flow-only inlet state");
	if (!std::isfinite(reference_inlet_flow_m3_s)
		|| !(std::abs(reference_inlet_flow_m3_s)
			> configuration.coupling.flow_epsilon_m3_s))
		throw std::runtime_error("CPU 3D VCA bridge reference inlet profile has zero flow");
	if (system.kind != EquationKind::NavierStokes || system.unknowns.empty())
		throw std::runtime_error("CPU 3D VCA bridge requires a Navier-Stokes system");
	const auto& velocity = system.unknowns.front();
	bool applied = false;
	for (auto& boundary : configuration.boundaries) {
		if (boundary.label != configuration.coupling.three_d_ports.inlet_label)
			continue;
		for (auto& condition : boundary.conditions) {
			if (condition.field != velocity) continue;
			if (condition.kind != FieldBoundaryKind::Dirichlet
				|| condition.profile != "initial_velocityfield.txt")
				throw std::runtime_error("CPU 3D VCA inlet requires an initial_velocityfield.txt velocity Dirichlet profile");
			condition.scale = -inlet.flow_m3_s/reference_inlet_flow_m3_s;
			applied = true;
		}
	}
	if (!applied)
		throw std::runtime_error("CPU 3D VCA inlet label has no velocity profile boundary");
}

inline VascularStepResult BuildThreeDFlowPortResult(double time_s, double dt_s,
	const VascularInletState& inlet, const ThreeDVascularPortDefinition& ports,
	const std::map<int, double>& outlet_flows,
	const std::map<int, double>& outlet_pressures)
{
	RequireFiniteCouplingValue("3D VCA result time", time_s);
	if (!(dt_s >= 0.0) || !std::isfinite(dt_s))
		throw std::runtime_error("3D VCA result timestep must be finite and nonnegative");
	ValidateVascularInletState(inlet);
	VascularStepResult result;
	result.time_s = time_s;
	result.dt_s = dt_s;
	result.inlet = inlet;
	for (const int label : ports.outlet_labels) {
		const auto flow = outlet_flows.find(label);
		const auto pressure = outlet_pressures.find(label);
		if (flow == outlet_flows.end() || pressure == outlet_pressures.end())
			throw std::runtime_error("3D VCA result is missing configured outlet label "+std::to_string(label));
		result.outlets.push_back({label, flow->second, pressure->second, {}, {}, {}, true});
	}
	return result;
}

} // namespace iga

#endif

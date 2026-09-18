#ifndef GENERIC_CASE_INPUT_HPP
#define GENERIC_CASE_INPUT_HPP

#include "CaseInput.hpp"

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <vector>

namespace iga {

struct ResolvedScalarBoundaries {
	std::vector<int> constrained;
	std::vector<double> value;
	std::size_t constrained_dofs = 0;
};

inline const EquationSystemDefinition& FindEquationSystem(
	const SimulationConfiguration& configuration, const std::string& name)
{
	for (const auto& system : configuration.equation_systems)
		if (system.name == name) return system;
	throw std::runtime_error("unknown equation system '" + name + "'");
}

inline std::string FirstLinearTransportSystem(const SimulationConfiguration& configuration)
{
	for (const auto& system : configuration.equation_systems)
		if (system.kind == EquationKind::LinearTransport) return system.name;
	throw std::runtime_error("simulation configuration has no linear_transport equation system");
}

inline ResolvedScalarBoundaries ResolveScalarBoundaries(
	const SimulationConfiguration& configuration,
	const CompiledLinearSystem& system,
	const std::vector<int>& labels)
{
	for (const auto& field_name : system.fields) {
		const auto found = std::find_if(configuration.fields.begin(), configuration.fields.end(),
			[&](const FieldDefinition& field) { return field.name == field_name; });
		if (found == configuration.fields.end() || found->kind != FieldKind::Scalar)
			throw std::runtime_error("linear_transport currently requires scalar fields: " + field_name);
	}
	const auto fields = system.fields.size();
	ResolvedScalarBoundaries result;
	result.constrained.assign(labels.size()*fields, 0);
	result.value.assign(labels.size()*fields, 0.0);
	std::set<int> mesh_labels;
	for (int label : labels) if (label >= 0) mesh_labels.insert(label);
	std::set<int> configured_labels;
	for (const auto& boundary : configuration.boundaries) {
		if (!mesh_labels.count(boundary.label))
			throw std::runtime_error("simulation_config.json boundary label " + std::to_string(boundary.label)
				+ " is not present in controlmesh.vtk");
		configured_labels.insert(boundary.label);
		for (const auto& condition : boundary.conditions) {
			const auto field = system.field_index.find(condition.field);
			if (field == system.field_index.end()) continue;
			if (!condition.waveform.empty())
				throw std::runtime_error("time-dependent scalar boundary execution is not implemented");
			if (condition.kind == FieldBoundaryKind::NoFlux || condition.kind == FieldBoundaryKind::AdvectiveOutflow
				|| condition.kind == FieldBoundaryKind::Flux || condition.kind == FieldBoundaryKind::Robin) continue;
			if (condition.kind != FieldBoundaryKind::Dirichlet)
				throw std::runtime_error("unsupported scalar boundary condition for field '" + condition.field + "'");
			if (condition.value.size() != 1)
				throw std::runtime_error("scalar Dirichlet field '" + condition.field + "' requires one value");
			for (std::size_t node = 0; node < labels.size(); ++node)
				if (labels[node] == boundary.label) {
					const auto dof = node*fields + field->second;
					result.constrained[dof] = 1;
					result.value[dof] = condition.value[0];
				}
		}
	}
	for (int label : mesh_labels)
		if (!configured_labels.count(label))
			throw std::runtime_error("mesh boundary label " + std::to_string(label)
				+ " has no simulation_config.json boundary definition");
	result.constrained_dofs = static_cast<std::size_t>(std::count(
		result.constrained.begin(), result.constrained.end(), 1));
	return result;
}

inline std::vector<double> InitialScalarValues(const SimulationConfiguration& configuration,
	const CompiledLinearSystem& system)
{
	std::vector<double> values(system.fields.size(), 0.0);
	for (std::size_t i = 0; i < system.fields.size(); ++i)
		for (const auto& field : configuration.fields)
			if (field.name == system.fields[i]) values[i] = field.initial_value;
	return values;
}

inline const EquationSystemDefinition& FirstNavierStokesSystem(
	const SimulationConfiguration& configuration)
{
	for (const auto& system : configuration.equation_systems)
		if (system.kind == EquationKind::NavierStokes) return system;
	throw std::runtime_error("simulation configuration has no navier_stokes equation system");
}

inline const EquationSystemDefinition& FindNavierStokesSystem(
	const SimulationConfiguration& configuration, const std::string& name)
{
	if (name.empty()) return FirstNavierStokesSystem(configuration);
	for (const auto& system : configuration.equation_systems)
		if (system.name == name) {
			if (system.kind != EquationKind::NavierStokes)
				throw std::runtime_error("equation system '"+name+"' is not navier_stokes");
			return system;
		}
	throw std::runtime_error("simulation configuration has no navier_stokes system named '"+name+"'");
}

inline ResolvedBoundaryConditions ResolveFlowBoundaries(
	const SimulationConfiguration& configuration, const EquationSystemDefinition& system,
	const std::vector<int>& labels,
	const std::vector<std::array<double, 3>>& reference_velocity)
{
	if (system.unknowns.size() != 2)
		throw std::runtime_error("navier_stokes requires [velocity, pressure] unknowns");
	const auto& velocity_name = system.unknowns[0];
	const auto& pressure_name = system.unknowns[1];
	ResolvedBoundaryConditions result;
	const auto nodes = labels.size();
	if (reference_velocity.size() != nodes)
		throw std::runtime_error("velocity profile and boundary labels have different node counts");
	result.velocity_constrained.assign(nodes, 0);
	result.pressure_constrained.assign(nodes, 0);
	result.transport_constrained.assign(nodes, 0);
	result.velocity.assign(nodes, {0.0, 0.0, 0.0});
	result.pressure.assign(nodes, 0.0);
	result.n0.assign(nodes, 0.0);
	result.nplus.assign(nodes, 0.0);
	std::set<int> mesh_labels;
	for (int label : labels) if (label >= 0) mesh_labels.insert(label);
	std::set<int> configured_labels;
	bool has_pressure_traction = false;
	for (const auto& boundary : configuration.boundaries) {
		if (!mesh_labels.count(boundary.label))
			throw std::runtime_error("simulation_config.json boundary label " + std::to_string(boundary.label)
				+ " is not present in controlmesh.vtk");
		configured_labels.insert(boundary.label);
		for (const auto& condition : boundary.conditions) {
			if (condition.field != velocity_name && condition.field != pressure_name) continue;
			if (!condition.waveform.empty())
				throw std::runtime_error("time-dependent Navier-Stokes boundary execution is not implemented");
			if (condition.kind == FieldBoundaryKind::PressureTraction) {
				if (condition.field != pressure_name || condition.value.size() != 1)
					throw std::runtime_error(
						"pressure_traction requires the Navier-Stokes pressure field");
				has_pressure_traction = true;
				if (condition.pressure_gauge)
					throw std::runtime_error(
						"pressure_traction already fixes the pressure reference; remove pressure_gauge "
						"to preserve every continuity equation");
				continue;
			}
			if (condition.kind != FieldBoundaryKind::Dirichlet)
				throw std::runtime_error("navier_stokes boundary field '" + condition.field
					+ "' currently supports only dirichlet conditions");
			if (condition.field == velocity_name) {
				if (condition.profile.empty() && condition.value.size() != 3)
					throw std::runtime_error("velocity Dirichlet condition requires three values or a profile");
				if (!condition.profile.empty() && condition.profile != "initial_velocityfield.txt")
					throw std::runtime_error("unsupported velocity profile '" + condition.profile + "'");
				for (std::size_t node = 0; node < nodes; ++node)
					if (labels[node] == boundary.label) {
						result.velocity_constrained[node] = 1;
						if (condition.profile.empty())
							result.velocity[node] = {condition.value[0], condition.value[1], condition.value[2]};
						else
							for (std::size_t component = 0; component < 3; ++component)
								result.velocity[node][component] = condition.scale*reference_velocity[node][component];
					}
			} else if (condition.field == pressure_name) {
				if (condition.value.size() != 1)
					throw std::runtime_error("pressure Dirichlet condition requires one value");
				for (std::size_t node = 0; node < nodes; ++node)
					if (labels[node] == boundary.label) {
						result.pressure_constrained[node] = 1;
						result.pressure[node] = condition.value[0];
					}
			}
		}
	}
	for (int label : mesh_labels)
		if (!configured_labels.count(label))
			throw std::runtime_error("mesh boundary label " + std::to_string(label)
				+ " has no simulation_config.json boundary definition");
	result.velocity_nodes = static_cast<std::size_t>(std::count(
		result.velocity_constrained.begin(), result.velocity_constrained.end(), 1));
	result.pressure_nodes = static_cast<std::size_t>(std::count(
		result.pressure_constrained.begin(), result.pressure_constrained.end(), 1));
	if (result.velocity_nodes == 0
		|| (result.pressure_nodes == 0 && !has_pressure_traction))
		throw std::runtime_error(
			"navier_stokes requires velocity Dirichlet and pressure gauge/traction boundaries");
	return result;
}
} // namespace iga

#endif

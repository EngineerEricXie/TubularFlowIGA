#ifndef CASE_INPUT_HPP
#define CASE_INPUT_HPP

#include "CaseConfig.hpp"
#include "SimulationConfig.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct TransportParameters {
	double diffusion = 0.0;
	double vplus = 0.0;
	double vminus = 0.0;
	double kplus = 0.0;
	double kminus = 0.0;
	double detach_plus = 0.0;
	double detach_minus = 0.0;
	double dt = 0.0;
	int steps = 0;
	double n0_bc = 0.0;
	double nplus_bc = 0.0;
	double nminus_bc = 0.0;
	double artificial_diffusion = 0.0;
};

struct ResolvedBoundaryConditions {
	std::vector<int> velocity_constrained;
	std::vector<int> pressure_constrained;
	std::vector<int> transport_constrained;
	std::vector<std::array<double, 3>> velocity;
	std::vector<double> pressure;
	std::vector<double> n0;
	std::vector<double> nplus;
	std::size_t velocity_nodes = 0;
	std::size_t pressure_nodes = 0;
	std::size_t transport_nodes = 0;
};

inline std::vector<int> ReadPointLabels(const std::string& path, std::uint64_t expected_nodes)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open mesh: " + path);
	std::string token;
	std::uint64_t point_data = 0;
	while (in >> token) {
		if (token == "POINT_DATA") {
			if (!(in >> point_data)) throw std::runtime_error("invalid POINT_DATA record");
			break;
		}
	}
	if (point_data != expected_nodes) throw std::runtime_error("VTK POINT_DATA count does not match database nodes");
	while (in >> token) if (token == "LOOKUP_TABLE") { in >> token; break; }
	if (!in) throw std::runtime_error("VTK point labels were not found");
	std::vector<int> labels(static_cast<std::size_t>(expected_nodes));
	for (auto& label : labels) {
		double value = 0.0;
		if (!(in >> value)) throw std::runtime_error("VTK point-label array is truncated");
		label = static_cast<int>(value);
	}
	return labels;
}

inline std::vector<std::array<double, 3>> ReadVelocity(const std::string& path, std::uint64_t expected_nodes)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open velocity field: " + path);
	std::vector<std::array<double, 3>> velocity(static_cast<std::size_t>(expected_nodes));
	for (auto& value : velocity)
		if (!(in >> value[0] >> value[1] >> value[2]))
			throw std::runtime_error("velocity field is truncated");
	std::string extra;
	if (in >> extra) throw std::runtime_error("velocity field contains extra records");
	return velocity;
}

inline TransportParameters ReadTransportParameters(const std::string& path)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open simulation parameters: " + path);
	std::map<std::string, double> values;
	std::string key;
	double value = 0.0;
	while (in >> key >> value) values[key] = value;
	auto required = [&](const std::string& name) {
		auto it = values.find(name);
		if (it == values.end()) throw std::runtime_error("missing transport parameter: " + name);
		return it->second;
	};
	TransportParameters p;
	p.diffusion = required("D");
	p.vplus = required("vplus");
	p.vminus = required("vminus");
	p.kplus = required("kplus");
	p.kminus = required("kminus");
	p.detach_plus = required("k'plus");
	p.detach_minus = required("k'minus");
	p.dt = required("dt");
	p.steps = static_cast<int>(required("nstep"));
	p.n0_bc = required("N0bc");
	p.nplus_bc = required("Nplusbc");
	p.nminus_bc = required("Nminusbc");
	auto artificial = values.find("artificial_diffusion");
	if (artificial != values.end()) p.artificial_diffusion = artificial->second;
	else {
		artificial = values.find("artificial_diffusion_weight");
		if (artificial != values.end()) p.artificial_diffusion = artificial->second;
	}
	if (p.dt <= 0.0 || p.steps < 0) throw std::runtime_error("dt and nstep must be positive");
	return p;
}

inline SimulationConfiguration ConvertLegacyNeuronTransport(const TransportParameters& p)
{
	SimulationConfiguration configuration;
	configuration.fields = {
		{"N0", FieldKind::Scalar, 0.0},
		{"Nplus", FieldKind::Scalar, 0.0}
	};
	EquationSystemDefinition system;
	system.name = "neuron_transport";
	system.kind = EquationKind::LinearTransport;
	system.unknowns = {"N0", "Nplus"};
	system.terms = {
		{TermKind::TimeDerivative, "N0", "N0", 1.0, ""},
		{TermKind::Diffusion, "N0", "N0", p.diffusion, ""},
		{TermKind::LinearCoupling, "N0", "N0", p.kplus+p.kminus, ""},
		{TermKind::LinearCoupling, "N0", "Nplus", -p.detach_plus, ""},
		{TermKind::TimeDerivative, "Nplus", "Nplus", 1.0, ""},
		{TermKind::Advection, "Nplus", "Nplus", 1.0, "prescribed"},
		{TermKind::Diffusion, "Nplus", "Nplus", p.vplus*p.artificial_diffusion, ""},
		{TermKind::LinearCoupling, "Nplus", "N0", -p.kplus, ""},
		{TermKind::LinearCoupling, "Nplus", "Nplus", p.detach_plus, ""}
	};
	system.stabilization = {{"Nplus", "supg", "prescribed"}};
	configuration.equation_systems.push_back(std::move(system));
	configuration.time = {p.dt, p.steps};
	configuration.boundaries = {
		{0, "wall", {{"N0", FieldBoundaryKind::NoFlux, {}, 0.0, 0.0, "", 1.0},
			{"Nplus", FieldBoundaryKind::NoFlux, {}, 0.0, 0.0, "", 1.0}}},
		{1, "inlet", {{"N0", FieldBoundaryKind::Dirichlet, {p.n0_bc}, 0.0, 0.0, "", 1.0},
			{"Nplus", FieldBoundaryKind::Dirichlet, {p.nplus_bc}, 0.0, 0.0, "", 1.0}}},
		{2, "outlet", {{"N0", FieldBoundaryKind::AdvectiveOutflow, {}, 0.0, 0.0, "", 1.0},
			{"Nplus", FieldBoundaryKind::AdvectiveOutflow, {}, 0.0, 0.0, "", 1.0}}}
	};
	return configuration;
}

inline ResolvedBoundaryConditions ResolveBoundaryConditions(
	const CaseConfiguration& configuration,
	const std::vector<int>& labels,
	const std::vector<std::array<double, 3>>& reference_velocity,
	const TransportParameters& parameters)
{
	if (labels.size() != reference_velocity.size())
		throw std::runtime_error("boundary labels and velocity field have different node counts");
	ResolvedBoundaryConditions resolved;
	const auto nodes = labels.size();
	resolved.velocity_constrained.assign(nodes, 0);
	resolved.pressure_constrained.assign(nodes, 0);
	resolved.transport_constrained.assign(nodes, 0);
	resolved.velocity.assign(nodes, {0.0, 0.0, 0.0});
	resolved.pressure.assign(nodes, 0.0);
	resolved.n0.assign(nodes, 0.0);
	resolved.nplus.assign(nodes, 0.0);

	auto clear_node = [&](std::size_t node) {
		resolved.velocity_constrained[node] = 0;
		resolved.pressure_constrained[node] = 0;
		resolved.transport_constrained[node] = 0;
		resolved.velocity[node] = {0.0, 0.0, 0.0};
		resolved.pressure[node] = 0.0;
		resolved.n0[node] = 0.0;
		resolved.nplus[node] = 0.0;
	};

	if (!configuration.present || configuration.inherit_legacy) {
		for (std::size_t node = 0; node < nodes; ++node) {
			if (labels[node] == 0) {
				resolved.velocity_constrained[node] = 1;
			} else if (labels[node] == 1) {
				resolved.velocity_constrained[node] = 1;
				for (std::size_t field = 0; field < 3; ++field)
					resolved.velocity[node][field] = parameters.vplus * reference_velocity[node][field];
				resolved.transport_constrained[node] = 1;
				resolved.n0[node] = parameters.n0_bc;
				resolved.nplus[node] = parameters.nplus_bc;
			} else if (labels[node] >= 2) {
				resolved.pressure_constrained[node] = 1;
			}
		}
	}

	std::set<int> configured_labels;
	for (const auto& rule : configuration.boundaries) {
		configured_labels.insert(rule.label);
		bool found = false;
		for (std::size_t node = 0; node < nodes; ++node) {
			if (labels[node] != rule.label) continue;
			found = true;
			clear_node(node);
			if (rule.type == BoundaryType::Wall) {
				resolved.velocity_constrained[node] = 1;
				if (rule.has_velocity) resolved.velocity[node] = rule.velocity;
			} else if (rule.type == BoundaryType::Inlet) {
				resolved.velocity_constrained[node] = 1;
				if (rule.has_velocity) {
					resolved.velocity[node] = rule.velocity;
				} else {
					const double scale = rule.has_velocity_scale ? rule.velocity_scale : parameters.vplus;
					for (std::size_t field = 0; field < 3; ++field)
						resolved.velocity[node][field] = scale * reference_velocity[node][field];
				}
				resolved.transport_constrained[node] = 1;
				resolved.n0[node] = rule.has_transport ? rule.n0 : parameters.n0_bc;
				resolved.nplus[node] = rule.has_transport ? rule.nplus : parameters.nplus_bc;
			} else if (rule.type == BoundaryType::Outlet) {
				resolved.pressure_constrained[node] = 1;
				resolved.pressure[node] = rule.has_pressure ? rule.pressure : 0.0;
			}
		}
		if (!found)
			throw std::runtime_error("case_config.json boundary label " + std::to_string(rule.label)
				+ " is not present in controlmesh.vtk");
	}

	if (configuration.present && !configuration.inherit_legacy) {
		for (int label : labels)
			if (label >= 0 && !configured_labels.count(label))
				throw std::runtime_error("mesh boundary label " + std::to_string(label)
					+ " has no case_config.json rule while inherit_legacy is false");
	}

	resolved.velocity_nodes = static_cast<std::size_t>(std::count(
		resolved.velocity_constrained.begin(), resolved.velocity_constrained.end(), 1));
	resolved.pressure_nodes = static_cast<std::size_t>(std::count(
		resolved.pressure_constrained.begin(), resolved.pressure_constrained.end(), 1));
	resolved.transport_nodes = static_cast<std::size_t>(std::count(
		resolved.transport_constrained.begin(), resolved.transport_constrained.end(), 1));
	if (resolved.velocity_nodes == 0)
		throw std::runtime_error("case has no velocity boundary nodes");
	if (resolved.pressure_nodes == 0)
		throw std::runtime_error("case has no pressure boundary nodes");
	return resolved;
}

} // namespace iga

#endif

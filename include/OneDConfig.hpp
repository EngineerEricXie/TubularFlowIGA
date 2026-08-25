#ifndef IGA_ONE_D_CONFIG_HPP
#define IGA_ONE_D_CONFIG_HPP

#include "SimulationConfig.hpp"
#include "RadiusAnnotatedObj.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class OneDFlowModel { Rigid, Compliant };
enum class OneDFlowScheme { SteadyPoiseuille, ExplicitRusanov, ImplicitPetsc };
enum class OneDImplicitFormulation { PressureNetwork, LinearizedAQ, NonlinearAQ, ImplicitPde };
enum class OneDOutletKind { Pressure, Resistance, WindkesselRcr };
enum class OneDWallModel { Linear, Olufsen };
enum class OneDWallBoundaryKind { NoFlux, ConstantFlux, Robin };

struct OneDGeometryDefinition {
	std::string kind = "swc_network";
	std::string file;
	double length_scale_to_m = 1.0;
	int root_node_id = 0;
};

struct OneDTimeDefinition {
	double dt = 0.0;
	int steps = 0;
	int output_every = 1;
};

struct OneDWallDefinition {
	OneDWallModel model = OneDWallModel::Linear;
	double young_modulus = 1.0e6;
	double thickness_ratio = 0.1;
	double reference_pressure = 0.0;
	double olufsen_k1 = 0.0;
	double olufsen_k2 = -22.5267;
	double olufsen_k3 = 8.65e5;
};

struct OneDDiscretizationDefinition {
	int cells_per_segment = 2;
	double cfl = 0.75;
	double alpha = 4.0/3.0;
	double min_area_fraction = 0.2;
};

struct OneDJunctionDefinition {
	std::string pressure_balance = "static";
	std::string loss_model = "none";
	double coefficient = 0.0;
	std::string reference_velocity = "parent";
	std::map<int, double> node_coefficients;
	std::vector<std::pair<double, double>> angle_table;
};

struct OneDFlowSystemDefinition {
	std::string name;
	std::vector<std::string> unknowns;
	OneDFlowModel model = OneDFlowModel::Rigid;
	OneDFlowScheme scheme = OneDFlowScheme::SteadyPoiseuille;
	OneDImplicitFormulation formulation = OneDImplicitFormulation::PressureNetwork;
	double dynamic_viscosity = 0.0;
	double density = 1.0;
	OneDWallDefinition wall;
	OneDDiscretizationDefinition discretization;
	OneDJunctionDefinition junctions;
};

struct OneDSpeciesDefinition {
	std::string field;
	double diffusivity = 0.0;
	double reaction_rate = 0.0;
	double volume_source = 0.0;
};

struct OneDTransportSystemDefinition {
	std::string name;
	std::string flow_system;
	std::vector<std::string> unknowns;
	std::vector<OneDSpeciesDefinition> species;
};

struct OneDBoundaryCondition {
	std::string field;
	std::string type;
	double value = 0.0;
	std::string waveform;
	std::string quantity;
	double resistance = 0.0;
	double proximal_resistance = 0.0;
	double distal_resistance = 0.0;
	double capacitance = 0.0;
	double reference_pressure = 0.0;
	double initial_pressure = 0.0;
	double coefficient = 0.0;
	double exterior_value = 0.0;
};

struct OneDBoundaryDefinition {
	std::string name;
	std::string role;
	std::vector<int> node_ids;
	std::vector<OneDBoundaryCondition> conditions;
};

struct OneDPhysiologyDefinition {
	bool enabled = false;
	std::map<std::string, double> metabolism_rates;
	bool oxygen_capacity = false;
	double hematocrit_percent = 0.0;
	double oxygen_solubility = 0.0031;
	double hemoglobin_g_dl = 15.0;
	double p50_mmhg = 26.8;
	double hill_exponent = 2.7;
	bool vasodilation = false;
	std::string vasodilator_field = "vasodilator";
	double emax_radius_fraction = 0.0;
	double ec50 = 1.0;
	double relaxation_tau = 1.0;
	std::vector<std::string> derived_fields;
};

struct OneDConfiguration {
	int schema_version = 3;
	std::string dimension = "1d";
	OneDGeometryDefinition geometry;
	OneDTimeDefinition time;
	std::vector<FieldDefinition> fields;
	std::vector<TemporalFunctionDefinition> temporal_functions;
	std::vector<OneDFlowSystemDefinition> flow_systems;
	std::vector<OneDTransportSystemDefinition> transport_systems;
	std::vector<OneDBoundaryDefinition> boundaries;
	OneDPhysiologyDefinition physiology;
};

namespace one_d_config_detail {

using config_detail::Find;
using config_detail::JsonValue;
using config_detail::RequireArray;
using config_detail::RequireBoolean;
using config_detail::RequireInteger;
using config_detail::RequireKnownKeys;
using config_detail::RequireNumber;
using config_detail::RequireObject;
using config_detail::RequireString;

inline const JsonValue& Required(const std::map<std::string, JsonValue>& object,
	const std::string& key, const std::string& context)
{
	const auto* value = Find(object, key);
	if (!value) throw std::runtime_error("simulation_config.json: " + context + " requires '" + key + "'");
	return *value;
}

inline bool OptionalBoolean(const std::map<std::string, JsonValue>& object,
	const std::string& key, bool fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireBoolean(*value, context + "." + key) : fallback;
}

inline double OptionalNumber(const std::map<std::string, JsonValue>& object,
	const std::string& key, double fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireNumber(*value, context + "." + key) : fallback;
}

inline std::string OptionalString(const std::map<std::string, JsonValue>& object,
	const std::string& key, const std::string& fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireString(*value, context + "." + key) : fallback;
}

inline FieldKind ParseOneDFieldKind(const std::string& value)
{
	if (value == "scalar") return FieldKind::Scalar;
	if (value == "pressure") return FieldKind::Pressure;
	throw std::runtime_error("simulation_config.json: 1d field kind must be 'scalar' or 'pressure'");
}

inline OneDWallDefinition ParseWall(const JsonValue* value, const std::string& context)
{
	OneDWallDefinition wall;
	if (!value) return wall;
	const auto& object = RequireObject(*value, context);
	RequireKnownKeys(object, {"model", "young_modulus", "thickness_ratio",
		"reference_pressure", "olufsen_k1", "olufsen_k2", "olufsen_k3"}, context);
	const auto model = OptionalString(object, "model", "linear", context);
	if (model == "linear") wall.model = OneDWallModel::Linear;
	else if (model == "olufsen") wall.model = OneDWallModel::Olufsen;
	else throw std::runtime_error("simulation_config.json: unsupported 1d wall model '" + model + "'");
	wall.young_modulus = OptionalNumber(object, "young_modulus", wall.young_modulus, context);
	wall.thickness_ratio = OptionalNumber(object, "thickness_ratio", wall.thickness_ratio, context);
	wall.reference_pressure = OptionalNumber(object, "reference_pressure", wall.reference_pressure, context);
	wall.olufsen_k1 = OptionalNumber(object, "olufsen_k1", wall.olufsen_k1, context);
	wall.olufsen_k2 = OptionalNumber(object, "olufsen_k2", wall.olufsen_k2, context);
	wall.olufsen_k3 = OptionalNumber(object, "olufsen_k3", wall.olufsen_k3, context);
	if (!(wall.young_modulus > 0.0) || !(wall.thickness_ratio > 0.0))
		throw std::runtime_error("simulation_config.json: 1d wall modulus and thickness ratio must be positive");
	return wall;
}

inline OneDDiscretizationDefinition ParseDiscretization(const JsonValue* value,
	const std::string& context)
{
	OneDDiscretizationDefinition result;
	if (!value) return result;
	const auto& object = RequireObject(*value, context);
	RequireKnownKeys(object, {"cells_per_segment", "cfl", "alpha", "min_area_fraction"}, context);
	if (const auto* cells = Find(object, "cells_per_segment"))
		result.cells_per_segment = RequireInteger(*cells, context + ".cells_per_segment");
	result.cfl = OptionalNumber(object, "cfl", result.cfl, context);
	result.alpha = OptionalNumber(object, "alpha", result.alpha, context);
	result.min_area_fraction = OptionalNumber(object, "min_area_fraction", result.min_area_fraction, context);
	if (result.cells_per_segment < 1 || !(result.cfl > 0.0 && result.cfl <= 1.0)
		|| !(result.alpha >= 1.0) || !(result.min_area_fraction > 0.0 && result.min_area_fraction < 1.0))
		throw std::runtime_error("simulation_config.json: invalid 1d discretization controls");
	return result;
}

inline OneDJunctionDefinition ParseJunctions(const JsonValue* value,
	const std::string& context)
{
	OneDJunctionDefinition result;
	if (!value) return result;
	const auto& object = RequireObject(*value, context);
	RequireKnownKeys(object, {"pressure_balance", "loss_model", "coefficient",
		"reference_velocity", "node_coefficients", "angle_table"}, context);
	result.pressure_balance = OptionalString(object, "pressure_balance", result.pressure_balance, context);
	result.loss_model = OptionalString(object, "loss_model", result.loss_model, context);
	result.coefficient = OptionalNumber(object, "coefficient", result.coefficient, context);
	result.reference_velocity = OptionalString(object, "reference_velocity", result.reference_velocity, context);
	const std::set<std::string> balances{"static", "total"};
	const std::set<std::string> models{"none", "constant", "table", "angle_sin2", "mynard_valen_sendstad"};
	if (!balances.count(result.pressure_balance) || !models.count(result.loss_model)
		|| result.coefficient < 0.0
		|| (result.reference_velocity != "parent" && result.reference_velocity != "child"))
		throw std::runtime_error("simulation_config.json: invalid 1d junction configuration");
	if (const auto* coefficients = Find(object, "node_coefficients")) {
		const auto& items = RequireObject(*coefficients, context + ".node_coefficients");
		for (const auto& item : items) {
			std::size_t used = 0;
			int node = -1;
			try { node = std::stoi(item.first, &used); }
			catch (const std::exception&) { used = 0; }
			if (used != item.first.size() || node < 0)
				throw std::runtime_error("simulation_config.json: junction node coefficient keys must be node ids");
			const double coefficient = RequireNumber(item.second,
				context + ".node_coefficients." + item.first);
			if (coefficient < 0.0)
				throw std::runtime_error("simulation_config.json: junction loss coefficients cannot be negative");
			result.node_coefficients.emplace(node, coefficient);
		}
	}
	if (const auto* table = Find(object, "angle_table")) {
		const auto& rows = RequireArray(*table, context + ".angle_table");
		for (std::size_t i = 0; i < rows.size(); ++i) {
			const auto row_context = context + ".angle_table[" + std::to_string(i) + "]";
			const auto& row = RequireObject(rows[i], row_context);
			RequireKnownKeys(row, {"angle_degrees", "coefficient"}, row_context);
			const double angle = RequireNumber(Required(row, "angle_degrees", row_context),
				row_context + ".angle_degrees");
			const double coefficient = RequireNumber(Required(row, "coefficient", row_context),
				row_context + ".coefficient");
			if (angle < 0.0 || angle > 180.0 || coefficient < 0.0)
				throw std::runtime_error("simulation_config.json: junction angle table values are invalid");
			result.angle_table.emplace_back(angle, coefficient);
		}
		std::sort(result.angle_table.begin(), result.angle_table.end());
		for (std::size_t i = 1; i < result.angle_table.size(); ++i)
			if (result.angle_table[i-1].first == result.angle_table[i].first)
				throw std::runtime_error("simulation_config.json: junction angle table has duplicate angles");
	}
	if (result.loss_model == "table" && result.angle_table.empty()
		&& result.node_coefficients.empty())
		throw std::runtime_error("simulation_config.json: table junction loss requires angle_table or node_coefficients");
	if (result.loss_model != "table" && !result.angle_table.empty())
		throw std::runtime_error("simulation_config.json: angle_table requires table junction loss");
	return result;
}

inline TemporalFunctionDefinition ParseTemporalFunction(const JsonValue& value,
	const std::string& context)
{
	const auto& object = RequireObject(value, context);
	RequireKnownKeys(object, {"name", "kind", "units", "value", "mean", "amplitude",
		"period", "phase", "file", "interpolation", "cosine", "sine"}, context);
	TemporalFunctionDefinition function;
	function.name = RequireString(Required(object, "name", context), context + ".name");
	const auto kind = RequireString(Required(object, "kind", context), context + ".kind");
	if (kind == "constant") function.kind = TemporalFunctionKind::Constant;
	else if (kind == "sinusoid") function.kind = TemporalFunctionKind::Sinusoid;
	else if (kind == "periodic_table") function.kind = TemporalFunctionKind::PeriodicTable;
	else if (kind == "fourier") function.kind = TemporalFunctionKind::Fourier;
	else throw std::runtime_error("simulation_config.json: unsupported temporal function kind '" + kind + "'");
	function.units = RequireString(Required(object, "units", context), context + ".units");
	function.value = OptionalNumber(object, "value", function.value, context);
	function.mean = OptionalNumber(object, "mean", function.mean, context);
	function.amplitude = OptionalNumber(object, "amplitude", function.amplitude, context);
	function.period = OptionalNumber(object, "period", function.period, context);
	function.phase = OptionalNumber(object, "phase", function.phase, context);
	function.file = OptionalString(object, "file", "", context);
	function.interpolation = OptionalString(object, "interpolation", "", context);
	auto parse_coefficients = [&](const char* key, std::vector<double>& target) {
		if (const auto* values = Find(object, key)) {
			const auto& array = RequireArray(*values, context + "." + key);
			for (std::size_t i = 0; i < array.size(); ++i)
				target.push_back(RequireNumber(array[i], context + "." + key + "[" + std::to_string(i) + "]"));
		}
	};
	parse_coefficients("cosine", function.cosine);
	parse_coefficients("sine", function.sine);
	if (function.name.empty() || function.units.empty())
		throw std::runtime_error("simulation_config.json: temporal name and units cannot be empty");
	if (function.kind == TemporalFunctionKind::Sinusoid && !(function.period > 0.0))
		throw std::runtime_error("simulation_config.json: sinusoid requires positive period");
	if (function.kind == TemporalFunctionKind::PeriodicTable
		&& (!(function.period > 0.0) || function.file.empty() || function.interpolation != "linear"))
		throw std::runtime_error("simulation_config.json: periodic_table requires period, file, and linear interpolation");
	if (function.kind == TemporalFunctionKind::Fourier
		&& (!(function.period > 0.0) || function.cosine.empty()
			|| function.cosine.size() != function.sine.size()))
		throw std::runtime_error("simulation_config.json: fourier coefficients must be equal non-empty arrays");
	return function;
}

} // namespace one_d_config_detail

inline OneDConfiguration ParseOneDConfiguration(const std::string& text)
{
	using namespace one_d_config_detail;
	const auto root_value = config_detail::JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "root");
	RequireKnownKeys(root, {"schema_version", "dimension", "geometry", "fields", "time",
		"temporal_functions", "equation_systems", "boundaries", "physiology"}, "root");
	OneDConfiguration result;
	result.schema_version = RequireInteger(Required(root, "schema_version", "root"), "schema_version");
	result.dimension = RequireString(Required(root, "dimension", "root"), "dimension");
	if (result.schema_version != 3 || result.dimension != "1d")
		throw std::runtime_error("simulation_config.json: 1d configuration requires schema_version 3 and dimension '1d'");

	const auto& geometry = RequireObject(Required(root, "geometry", "root"), "geometry");
	RequireKnownKeys(geometry, {"kind", "file", "length_scale_to_m", "root_node_id"}, "geometry");
	result.geometry.kind = RequireString(Required(geometry, "kind", "geometry"), "geometry.kind");
	result.geometry.file = RequireString(Required(geometry, "file", "geometry"), "geometry.file");
	result.geometry.length_scale_to_m = RequireNumber(
		Required(geometry, "length_scale_to_m", "geometry"), "geometry.length_scale_to_m");
	if (const auto* root_id = Find(geometry, "root_node_id"))
		result.geometry.root_node_id = RequireInteger(*root_id, "geometry.root_node_id");
	if ((result.geometry.kind != "swc_network" && result.geometry.kind != "obj_network")
		|| result.geometry.file.empty()
		|| std::filesystem::path(result.geometry.file).is_absolute()
		|| !(result.geometry.length_scale_to_m > 0.0)
		|| result.geometry.root_node_id < 0)
		throw std::runtime_error("simulation_config.json: geometry requires swc_network or obj_network, a relative file, and positive length_scale_to_m");
	if (result.geometry.kind == "swc_network" && result.geometry.root_node_id != 0)
		throw std::runtime_error("simulation_config.json: root_node_id applies only to obj_network");
	if ((result.geometry.kind == "obj_network") != IsRadiusAnnotatedObjPath(result.geometry.file))
		throw std::runtime_error("simulation_config.json: geometry kind and skeleton file extension do not match");

	const auto& time = RequireObject(Required(root, "time", "root"), "time");
	RequireKnownKeys(time, {"dt", "steps", "output_every"}, "time");
	result.time.dt = RequireNumber(Required(time, "dt", "time"), "time.dt");
	result.time.steps = RequireInteger(Required(time, "steps", "time"), "time.steps");
	result.time.output_every = RequireInteger(Required(time, "output_every", "time"), "time.output_every");
	if (!(result.time.dt > 0.0) || result.time.steps < 1 || result.time.output_every < 1)
		throw std::runtime_error("simulation_config.json: time.dt, steps, and output_every must be positive");

	std::set<std::string> field_names;
	const auto& fields = RequireArray(Required(root, "fields", "root"), "fields");
	for (std::size_t i = 0; i < fields.size(); ++i) {
		const auto context = "fields[" + std::to_string(i) + "]";
		const auto& object = RequireObject(fields[i], context);
		RequireKnownKeys(object, {"name", "kind", "initial_value"}, context);
		FieldDefinition field;
		field.name = RequireString(Required(object, "name", context), context + ".name");
		field.kind = ParseOneDFieldKind(RequireString(Required(object, "kind", context), context + ".kind"));
		field.initial_value = OptionalNumber(object, "initial_value", 0.0, context);
		if (field.name.empty() || !field_names.insert(field.name).second)
			throw std::runtime_error("simulation_config.json: 1d field names must be non-empty and unique");
		result.fields.push_back(std::move(field));
	}

	std::set<std::string> temporal_names;
	if (const auto* functions = Find(root, "temporal_functions")) {
		const auto& array = RequireArray(*functions, "temporal_functions");
		for (std::size_t i = 0; i < array.size(); ++i) {
			auto function = ParseTemporalFunction(array[i], "temporal_functions[" + std::to_string(i) + "]");
			if (!temporal_names.insert(function.name).second)
				throw std::runtime_error("simulation_config.json: duplicate temporal function name");
			result.temporal_functions.push_back(std::move(function));
		}
	}

	std::set<std::string> system_names;
	const auto& systems = RequireArray(Required(root, "equation_systems", "root"), "equation_systems");
	for (std::size_t i = 0; i < systems.size(); ++i) {
		const auto context = "equation_systems[" + std::to_string(i) + "]";
		const auto& object = RequireObject(systems[i], context);
		const auto name = RequireString(Required(object, "name", context), context + ".name");
		const auto kind = RequireString(Required(object, "kind", context), context + ".kind");
		if (name.empty() || !system_names.insert(name).second)
			throw std::runtime_error("simulation_config.json: equation-system names must be non-empty and unique");
		const auto& unknown_array = RequireArray(Required(object, "unknowns", context), context + ".unknowns");
		std::vector<std::string> unknowns;
		std::set<std::string> local;
		for (std::size_t j = 0; j < unknown_array.size(); ++j) {
			const auto field = RequireString(unknown_array[j], context + ".unknowns[" + std::to_string(j) + "]");
			if (!field_names.count(field) || !local.insert(field).second)
				throw std::runtime_error("simulation_config.json: system unknowns must name unique declared fields");
			unknowns.push_back(field);
		}
		if (kind == "network_flow_1d") {
			RequireKnownKeys(object, {"name", "kind", "unknowns", "model", "scheme", "formulation",
				"dynamic_viscosity", "density", "wall", "discretization", "junctions"}, context);
			OneDFlowSystemDefinition flow;
			flow.name = name;
			flow.unknowns = unknowns;
			const auto model = RequireString(Required(object, "model", context), context + ".model");
			const auto scheme = RequireString(Required(object, "scheme", context), context + ".scheme");
			if (model == "rigid") flow.model = OneDFlowModel::Rigid;
			else if (model == "compliant") flow.model = OneDFlowModel::Compliant;
			else throw std::runtime_error("simulation_config.json: 1d flow model must be rigid or compliant");
			if (scheme == "steady_poiseuille") flow.scheme = OneDFlowScheme::SteadyPoiseuille;
			else if (scheme == "explicit_rusanov") flow.scheme = OneDFlowScheme::ExplicitRusanov;
			else if (scheme == "implicit_petsc") flow.scheme = OneDFlowScheme::ImplicitPetsc;
			else throw std::runtime_error("simulation_config.json: unsupported 1d flow scheme '" + scheme + "'");
			if ((flow.model == OneDFlowModel::Rigid) != (flow.scheme == OneDFlowScheme::SteadyPoiseuille))
				throw std::runtime_error("simulation_config.json: rigid requires steady_poiseuille and compliant requires explicit_rusanov or implicit_petsc");
			const auto formulation = OptionalString(object, "formulation", "pressure_network", context);
			if (formulation == "pressure_network") flow.formulation = OneDImplicitFormulation::PressureNetwork;
			else if (formulation == "linearized_aq") flow.formulation = OneDImplicitFormulation::LinearizedAQ;
			else if (formulation == "nonlinear_aq") flow.formulation = OneDImplicitFormulation::NonlinearAQ;
			else if (formulation == "implicit_1d_pde") flow.formulation = OneDImplicitFormulation::ImplicitPde;
			else throw std::runtime_error("simulation_config.json: unsupported implicit 1d formulation '" + formulation + "'");
			if (flow.scheme != OneDFlowScheme::ImplicitPetsc && Find(object, "formulation"))
				throw std::runtime_error("simulation_config.json: formulation applies only to implicit_petsc");
			flow.dynamic_viscosity = RequireNumber(Required(object, "dynamic_viscosity", context), context + ".dynamic_viscosity");
			flow.density = RequireNumber(Required(object, "density", context), context + ".density");
			const std::set<std::string> required_flow_fields{"area", "flow_rate", "pressure"};
			if (!(flow.dynamic_viscosity > 0.0) || !(flow.density > 0.0)
				|| std::set<std::string>(flow.unknowns.begin(), flow.unknowns.end()) != required_flow_fields)
				throw std::runtime_error("simulation_config.json: 1d flow requires positive fluid properties and area, flow_rate, pressure unknowns");
			flow.wall = ParseWall(Find(object, "wall"), context + ".wall");
			flow.discretization = ParseDiscretization(Find(object, "discretization"), context + ".discretization");
			flow.junctions = ParseJunctions(Find(object, "junctions"), context + ".junctions");
			result.flow_systems.push_back(std::move(flow));
		} else if (kind == "network_transport_1d") {
			RequireKnownKeys(object, {"name", "kind", "unknowns", "flow_system", "species"}, context);
			OneDTransportSystemDefinition transport;
			transport.name = name;
			transport.unknowns = unknowns;
			transport.flow_system = RequireString(Required(object, "flow_system", context), context + ".flow_system");
			const auto& species = RequireArray(Required(object, "species", context), context + ".species");
			std::set<std::string> configured;
			for (std::size_t j = 0; j < species.size(); ++j) {
				const auto species_context = context + ".species[" + std::to_string(j) + "]";
				const auto& item = RequireObject(species[j], species_context);
				RequireKnownKeys(item, {"field", "diffusivity", "reaction_rate", "volume_source"}, species_context);
				OneDSpeciesDefinition definition;
				definition.field = RequireString(Required(item, "field", species_context), species_context + ".field");
				definition.diffusivity = OptionalNumber(item, "diffusivity", 0.0, species_context);
				definition.reaction_rate = OptionalNumber(item, "reaction_rate", 0.0, species_context);
				definition.volume_source = OptionalNumber(item, "volume_source", 0.0, species_context);
				if (!local.count(definition.field) || !configured.insert(definition.field).second
					|| definition.diffusivity < 0.0 || definition.reaction_rate < 0.0)
					throw std::runtime_error("simulation_config.json: invalid 1d species definition");
				transport.species.push_back(std::move(definition));
			}
			if (transport.species.size() != transport.unknowns.size())
				throw std::runtime_error("simulation_config.json: every 1d transport unknown requires one species definition");
			result.transport_systems.push_back(std::move(transport));
		} else {
			throw std::runtime_error("simulation_config.json: schema v3 1d supports network_flow_1d and network_transport_1d systems");
		}
	}
	if (result.flow_systems.empty()) throw std::runtime_error("simulation_config.json: 1d configuration requires a flow system");
	for (const auto& transport : result.transport_systems) {
		const auto found = std::find_if(result.flow_systems.begin(), result.flow_systems.end(),
			[&](const OneDFlowSystemDefinition& flow) { return flow.name == transport.flow_system; });
		if (found == result.flow_systems.end())
			throw std::runtime_error("simulation_config.json: 1d transport references unknown flow system '" + transport.flow_system + "'");
	}

	const auto& boundaries = RequireArray(Required(root, "boundaries", "root"), "boundaries");
	std::set<std::string> boundary_names;
	for (std::size_t i = 0; i < boundaries.size(); ++i) {
		const auto context = "boundaries[" + std::to_string(i) + "]";
		const auto& object = RequireObject(boundaries[i], context);
		RequireKnownKeys(object, {"name", "role", "node_ids", "conditions"}, context);
		OneDBoundaryDefinition boundary;
		boundary.name = RequireString(Required(object, "name", context), context + ".name");
		if (boundary.name.empty() || !boundary_names.insert(boundary.name).second)
			throw std::runtime_error("simulation_config.json: boundary names must be non-empty and unique");
		boundary.role = RequireString(Required(object, "role", context), context + ".role");
		if (boundary.role != "inlet" && boundary.role != "outlet" && boundary.role != "wall")
			throw std::runtime_error("simulation_config.json: 1d boundary role must be inlet, outlet, or wall");
		if (const auto* ids = Find(object, "node_ids")) {
			const auto& array = RequireArray(*ids, context + ".node_ids");
			std::set<int> unique_ids;
			for (std::size_t j = 0; j < array.size(); ++j) {
				const int id = RequireInteger(array[j], context + ".node_ids[" + std::to_string(j) + "]");
				if (id < 1 || !unique_ids.insert(id).second)
					throw std::runtime_error("simulation_config.json: boundary node_ids must be unique positive skeleton ids");
				boundary.node_ids.push_back(id);
			}
		}
		const auto& conditions = RequireArray(Required(object, "conditions", context), context + ".conditions");
		std::set<std::string> conditioned_fields;
		for (std::size_t j = 0; j < conditions.size(); ++j) {
			const auto item_context = context + ".conditions[" + std::to_string(j) + "]";
			const auto& item = RequireObject(conditions[j], item_context);
			RequireKnownKeys(item, {"field", "type", "value", "waveform", "quantity", "resistance",
				"proximal_resistance", "distal_resistance", "capacitance", "reference_pressure",
				"initial_pressure", "coefficient", "exterior_value"}, item_context);
			OneDBoundaryCondition condition;
			condition.field = RequireString(Required(item, "field", item_context), item_context + ".field");
			condition.type = RequireString(Required(item, "type", item_context), item_context + ".type");
			if (!field_names.count(condition.field))
				throw std::runtime_error("simulation_config.json: 1d boundary condition references unknown field");
			if (!conditioned_fields.insert(condition.field).second)
				throw std::runtime_error("simulation_config.json: a 1d boundary may define only one condition per field");
			condition.value = OptionalNumber(item, "value", 0.0, item_context);
			condition.waveform = OptionalString(item, "waveform", "", item_context);
			condition.quantity = OptionalString(item, "quantity", "", item_context);
			condition.resistance = OptionalNumber(item, "resistance", 0.0, item_context);
			condition.proximal_resistance = OptionalNumber(item, "proximal_resistance", 0.0, item_context);
			condition.distal_resistance = OptionalNumber(item, "distal_resistance", 0.0, item_context);
			condition.capacitance = OptionalNumber(item, "capacitance", 0.0, item_context);
			condition.reference_pressure = OptionalNumber(item, "reference_pressure", 0.0, item_context);
			condition.initial_pressure = OptionalNumber(item, "initial_pressure", condition.reference_pressure, item_context);
			condition.coefficient = OptionalNumber(item, "coefficient", 0.0, item_context);
			condition.exterior_value = OptionalNumber(item, "exterior_value", 0.0, item_context);
			if (!condition.waveform.empty() && !temporal_names.count(condition.waveform))
				throw std::runtime_error("simulation_config.json: boundary waveform names an unknown temporal function");
			const std::set<std::string> types{"dirichlet", "pressure", "resistance", "windkessel_rcr",
				"no_flux", "constant_flux", "robin"};
			if (!types.count(condition.type)) throw std::runtime_error("simulation_config.json: unsupported 1d boundary type '" + condition.type + "'");
			if (condition.type == "resistance" && !(condition.resistance > 0.0))
				throw std::runtime_error("simulation_config.json: 1d resistance outlet requires positive resistance");
			if (condition.type == "windkessel_rcr" && (condition.proximal_resistance < 0.0
				|| !(condition.distal_resistance > 0.0) || !(condition.capacitance > 0.0)))
				throw std::runtime_error("simulation_config.json: 1d RCR requires Rp>=0, Rd>0, and C>0");
			boundary.conditions.push_back(std::move(condition));
		}
		result.boundaries.push_back(std::move(boundary));
	}
	for (const auto& boundary : result.boundaries)
		for (const auto& condition : boundary.conditions) {
			if (condition.type == "robin" && condition.coefficient < 0.0)
				throw std::runtime_error("simulation_config.json: Robin coefficient cannot be negative");
			if (boundary.role != "inlet" || condition.type != "dirichlet"
				|| (condition.field != "flow_rate" && condition.field != "velocity")) continue;
			const auto quantity = condition.quantity.empty()
				? (condition.field == "velocity" ? "centerline_velocity" : "flow_rate")
				: condition.quantity;
			if (quantity != "flow_rate" && quantity != "centerline_velocity")
				throw std::runtime_error("simulation_config.json: inlet quantity must be flow_rate or centerline_velocity");
			if (condition.waveform.empty()) continue;
			const auto function = std::find_if(result.temporal_functions.begin(), result.temporal_functions.end(),
				[&](const TemporalFunctionDefinition& item) { return item.name == condition.waveform; });
			const bool units_ok = quantity == "flow_rate"
				? (function->units == "m3/s" || function->units == "m^3/s")
				: function->units == "m/s";
			if (!units_ok)
				throw std::runtime_error("simulation_config.json: inlet waveform units do not match quantity '"
					+ quantity + "'");
		}

	if (const auto* physiology = Find(root, "physiology")) {
		const auto& object = RequireObject(*physiology, "physiology");
		RequireKnownKeys(object, {"enabled", "metabolism", "oxygen_capacity", "vasodilation", "derived_fields"}, "physiology");
		result.physiology.enabled = OptionalBoolean(object, "enabled", true, "physiology");
		if (const auto* metabolism = Find(object, "metabolism")) {
			const auto& rates = RequireObject(*metabolism, "physiology.metabolism");
			for (const auto& rate : rates) {
				const double value = RequireNumber(rate.second, "physiology.metabolism." + rate.first);
				result.physiology.metabolism_rates.emplace(rate.first, value);
			}
		}
		if (const auto* oxygen = Find(object, "oxygen_capacity")) {
			const auto& item = RequireObject(*oxygen, "physiology.oxygen_capacity");
			RequireKnownKeys(item, {"enabled", "hematocrit_percent", "oxygen_solubility",
				"hemoglobin_g_dl", "p50_mmhg", "hill_exponent"}, "physiology.oxygen_capacity");
			result.physiology.oxygen_capacity = OptionalBoolean(item, "enabled", true, "physiology.oxygen_capacity");
			result.physiology.hematocrit_percent = OptionalNumber(item, "hematocrit_percent", 0.0, "physiology.oxygen_capacity");
			result.physiology.oxygen_solubility = OptionalNumber(item, "oxygen_solubility", 0.0031, "physiology.oxygen_capacity");
			result.physiology.hemoglobin_g_dl = OptionalNumber(item, "hemoglobin_g_dl", 15.0, "physiology.oxygen_capacity");
			result.physiology.p50_mmhg = OptionalNumber(item, "p50_mmhg", 26.8, "physiology.oxygen_capacity");
			result.physiology.hill_exponent = OptionalNumber(item, "hill_exponent", 2.7, "physiology.oxygen_capacity");
			if (result.physiology.hematocrit_percent < 0.0
				|| result.physiology.hematocrit_percent > 100.0
				|| !(result.physiology.oxygen_solubility > 0.0)
				|| !(result.physiology.hemoglobin_g_dl >= 0.0)
				|| !(result.physiology.p50_mmhg > 0.0)
				|| !(result.physiology.hill_exponent > 0.0))
				throw std::runtime_error("simulation_config.json: oxygen-capacity parameters are invalid");
		}
		if (const auto* vasodilation = Find(object, "vasodilation")) {
			const auto& item = RequireObject(*vasodilation, "physiology.vasodilation");
			RequireKnownKeys(item, {"enabled", "field", "emax_radius_fraction", "ec50", "relaxation_tau"}, "physiology.vasodilation");
			result.physiology.vasodilation = OptionalBoolean(item, "enabled", true, "physiology.vasodilation");
			result.physiology.vasodilator_field = OptionalString(item, "field", "vasodilator", "physiology.vasodilation");
			result.physiology.emax_radius_fraction = OptionalNumber(item, "emax_radius_fraction", 0.0, "physiology.vasodilation");
			result.physiology.ec50 = OptionalNumber(item, "ec50", 1.0, "physiology.vasodilation");
			result.physiology.relaxation_tau = OptionalNumber(item, "relaxation_tau", 1.0, "physiology.vasodilation");
			if (result.physiology.emax_radius_fraction < 0.0
				|| !(result.physiology.ec50 > 0.0) || !(result.physiology.relaxation_tau > 0.0))
				throw std::runtime_error("simulation_config.json: vasodilation ec50 and relaxation_tau must be positive");
		}
		if (const auto* derived = Find(object, "derived_fields")) {
			const auto& array = RequireArray(*derived, "physiology.derived_fields");
			std::set<std::string> unique;
			for (std::size_t i = 0; i < array.size(); ++i) {
				const auto name = RequireString(array[i], "physiology.derived_fields[" + std::to_string(i) + "]");
				if (!unique.insert(name).second)
					throw std::runtime_error("simulation_config.json: duplicate physiology derived field '" + name + "'");
				result.physiology.derived_fields.push_back(name);
			}
		}
	}
	if (result.physiology.enabled) {
		std::set<std::string> transported;
		for (const auto& transport : result.transport_systems)
			for (const auto& species : transport.species) transported.insert(species.field);
		for (const auto& rate : result.physiology.metabolism_rates)
			if (!transported.count(rate.first))
				throw std::runtime_error("simulation_config.json: metabolism references missing transported field '"
					+ rate.first + "'");
		if (result.physiology.vasodilation
			&& !transported.count(result.physiology.vasodilator_field))
			throw std::runtime_error("simulation_config.json: vasodilation requires transported field '"
				+ result.physiology.vasodilator_field + "'");
		const std::set<std::string> oxygen_fields{"pO2", "SaO2", "SvO2",
			"dissolved_oxygen", "bound_oxygen", "total_oxygen"};
		const std::set<std::string> allowed{"pO2", "pCO2", "pH", "SaO2", "SvO2",
			"dissolved_oxygen", "bound_oxygen", "total_oxygen", "hematocrit"};
		for (const auto& name : result.physiology.derived_fields) {
			if (!allowed.count(name))
				throw std::runtime_error("simulation_config.json: unsupported physiology derived field '"
					+ name + "'");
			if (oxygen_fields.count(name) && !transported.count("oxygen"))
				throw std::runtime_error("simulation_config.json: derived field '" + name
					+ "' requires transported oxygen");
			if (name == "pCO2" && !transported.count("carbon_dioxide"))
				throw std::runtime_error("simulation_config.json: pCO2 requires transported carbon_dioxide");
			if (name == "pH" && (!transported.count("carbon_dioxide")
				|| !transported.count("bicarbonate")))
				throw std::runtime_error("simulation_config.json: pH requires transported carbon_dioxide and bicarbonate");
			if ((name == "bound_oxygen" || name == "total_oxygen" || name == "hematocrit")
				&& !result.physiology.oxygen_capacity)
				throw std::runtime_error("simulation_config.json: derived field '" + name
					+ "' requires enabled oxygen_capacity");
		}
	}
	return result;
}

inline OneDConfiguration ReadOneDConfiguration(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open 1d simulation configuration: " + path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error("cannot read 1d simulation configuration: " + path.string());
	return ParseOneDConfiguration(contents.str());
}

inline const FieldDefinition& FindOneDField(const OneDConfiguration& configuration,
	const std::string& name)
{
	const auto found = std::find_if(configuration.fields.begin(), configuration.fields.end(),
		[&](const FieldDefinition& field) { return field.name == name; });
	if (found == configuration.fields.end()) throw std::runtime_error("unknown 1d field '" + name + "'");
	return *found;
}

inline const TemporalFunctionDefinition& FindOneDTemporalFunction(
	const OneDConfiguration& configuration, const std::string& name)
{
	const auto found = std::find_if(configuration.temporal_functions.begin(), configuration.temporal_functions.end(),
		[&](const TemporalFunctionDefinition& function) { return function.name == name; });
	if (found == configuration.temporal_functions.end()) throw std::runtime_error("unknown 1d temporal function '" + name + "'");
	return *found;
}

} // namespace iga

#endif

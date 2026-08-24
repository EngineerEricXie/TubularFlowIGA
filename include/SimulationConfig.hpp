#ifndef IGA_SIMULATION_CONFIG_HPP
#define IGA_SIMULATION_CONFIG_HPP

#include "CaseConfig.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class FieldKind { Scalar, Vector3, Pressure };
enum class EquationKind { LinearTransport, NavierStokes };
enum class TermKind { TimeDerivative, Diffusion, Advection, LinearCoupling, VolumeSource };
enum class FieldBoundaryKind { Dirichlet, NoFlux, Flux, Robin, AdvectiveOutflow };
enum class TemporalFunctionKind { Constant, Sinusoid, PeriodicTable, Fourier };

struct FieldDefinition {
	std::string name;
	FieldKind kind = FieldKind::Scalar;
	double initial_value = 0.0;
};

struct WeakFormTerm {
	TermKind kind = TermKind::TimeDerivative;
	std::string equation;
	std::string trial;
	double coefficient = 1.0;
	std::string velocity;
};

struct StabilizationDefinition {
	std::string equation;
	std::string method;
	std::string velocity;
};

struct EquationSystemDefinition {
	std::string name;
	EquationKind kind = EquationKind::LinearTransport;
	std::vector<std::string> unknowns;
	std::vector<WeakFormTerm> terms;
	std::vector<StabilizationDefinition> stabilization;
	double viscosity = 0.0;
	double density = 1.0;
	std::string time_integration = "steady";
};

struct FieldBoundaryCondition {
	std::string field;
	FieldBoundaryKind kind = FieldBoundaryKind::Dirichlet;
	std::vector<double> value;
	double coefficient = 0.0;
	double exterior_value = 0.0;
	std::string profile;
	double scale = 1.0;
	std::string waveform;
};

struct TemporalFunctionDefinition {
	std::string name;
	TemporalFunctionKind kind = TemporalFunctionKind::Constant;
	std::string units;
	double value = 1.0;
	double mean = 0.0;
	double amplitude = 0.0;
	double period = 0.0;
	double phase = 0.0;
	std::string file;
	std::string interpolation;
	std::vector<double> cosine;
	std::vector<double> sine;
};

struct NamedBoundaryDefinition {
	int label = -1;
	std::string name;
	std::vector<FieldBoundaryCondition> conditions;
};

struct TimeDefinition {
	double dt = 0.0;
	int steps = 0;
};

struct SimulationConfiguration {
	int schema_version = 2;
	std::vector<FieldDefinition> fields;
	std::vector<EquationSystemDefinition> equation_systems;
	std::vector<NamedBoundaryDefinition> boundaries;
	std::vector<TemporalFunctionDefinition> temporal_functions;
	TimeDefinition time;
};

struct CompiledTerm {
	TermKind kind = TermKind::TimeDerivative;
	std::size_t equation = 0;
	std::size_t trial = 0;
	double coefficient = 1.0;
	std::string velocity;
};

struct CompiledLinearSystem {
	std::string name;
	std::vector<std::string> fields;
	std::map<std::string, std::size_t> field_index;
	std::vector<CompiledTerm> terms;
	std::vector<StabilizationDefinition> stabilization;
	double dt = 0.0;
	int steps = 0;
};

namespace simulation_detail {

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

inline FieldKind ParseFieldKind(const std::string& value)
{
	if (value == "scalar") return FieldKind::Scalar;
	if (value == "vector3") return FieldKind::Vector3;
	if (value == "pressure") return FieldKind::Pressure;
	throw std::runtime_error("simulation_config.json: unsupported field kind '" + value + "'");
}

inline EquationKind ParseEquationKind(const std::string& value)
{
	if (value == "linear_transport") return EquationKind::LinearTransport;
	if (value == "navier_stokes") return EquationKind::NavierStokes;
	throw std::runtime_error("simulation_config.json: unsupported equation-system kind '" + value + "'");
}

inline TermKind ParseTermKind(const std::string& value)
{
	if (value == "time_derivative") return TermKind::TimeDerivative;
	if (value == "diffusion") return TermKind::Diffusion;
	if (value == "advection") return TermKind::Advection;
	if (value == "linear_coupling") return TermKind::LinearCoupling;
	if (value == "volume_source") return TermKind::VolumeSource;
	throw std::runtime_error("simulation_config.json: unsupported weak-form term '" + value + "'");
}

inline FieldBoundaryKind ParseFieldBoundaryKind(const std::string& value)
{
	if (value == "dirichlet") return FieldBoundaryKind::Dirichlet;
	if (value == "no_flux") return FieldBoundaryKind::NoFlux;
	if (value == "flux") return FieldBoundaryKind::Flux;
	if (value == "robin") return FieldBoundaryKind::Robin;
	if (value == "advective_outflow") return FieldBoundaryKind::AdvectiveOutflow;
	throw std::runtime_error("simulation_config.json: unsupported boundary condition '" + value + "'");
}

inline TemporalFunctionKind ParseTemporalFunctionKind(const std::string& value)
{
	if (value == "constant") return TemporalFunctionKind::Constant;
	if (value == "sinusoid") return TemporalFunctionKind::Sinusoid;
	if (value == "periodic_table") return TemporalFunctionKind::PeriodicTable;
	if (value == "fourier") return TemporalFunctionKind::Fourier;
	throw std::runtime_error("simulation_config.json: unsupported temporal function kind '" + value + "'");
}

inline std::vector<double> ParseValueVector(const JsonValue& value, const std::string& context)
{
	if (value.type == JsonValue::Type::Number) return {RequireNumber(value, context)};
	const auto& array = RequireArray(value, context);
	std::vector<double> result;
	for (std::size_t i = 0; i < array.size(); ++i)
		result.push_back(RequireNumber(array[i], context + "[" + std::to_string(i) + "]"));
	return result;
}

} // namespace simulation_detail

inline SimulationConfiguration ParseSimulationConfiguration(const std::string& text)
{
	using namespace simulation_detail;
	const auto root_value = config_detail::JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "root");
	RequireKnownKeys(root, {"schema_version", "fields", "equation_systems", "boundaries",
		"time", "temporal_functions"}, "root");
	SimulationConfiguration configuration;
	configuration.schema_version = RequireInteger(Required(root, "schema_version", "root"), "schema_version");
	if (configuration.schema_version != 2)
		throw std::runtime_error("simulation_config.json: only schema_version 2 is supported");

	const auto& fields = RequireArray(Required(root, "fields", "root"), "fields");
	std::set<std::string> field_names;
	for (std::size_t i = 0; i < fields.size(); ++i) {
		const auto context = "fields[" + std::to_string(i) + "]";
		const auto& object = RequireObject(fields[i], context);
		RequireKnownKeys(object, {"name", "kind", "initial_value"}, context);
		FieldDefinition field;
		field.name = RequireString(Required(object, "name", context), context + ".name");
		field.kind = ParseFieldKind(RequireString(Required(object, "kind", context), context + ".kind"));
		if (const auto* value = Find(object, "initial_value")) field.initial_value = RequireNumber(*value, context + ".initial_value");
		if (field.name.empty() || !field_names.insert(field.name).second)
			throw std::runtime_error("simulation_config.json: field names must be non-empty and unique");
		configuration.fields.push_back(std::move(field));
	}
	if (configuration.fields.empty()) throw std::runtime_error("simulation_config.json: fields cannot be empty");

	const auto& time = RequireObject(Required(root, "time", "root"), "time");
	RequireKnownKeys(time, {"dt", "steps"}, "time");
	configuration.time.dt = RequireNumber(Required(time, "dt", "time"), "time.dt");
	configuration.time.steps = RequireInteger(Required(time, "steps", "time"), "time.steps");
	if (!(configuration.time.dt > 0.0)) throw std::runtime_error("simulation_config.json: time.dt must be positive");

	std::set<std::string> temporal_names;
	if (const auto* temporal = Find(root, "temporal_functions")) {
		const auto& functions = RequireArray(*temporal, "temporal_functions");
		for (std::size_t i = 0; i < functions.size(); ++i) {
			const auto context = "temporal_functions[" + std::to_string(i) + "]";
			const auto& object = RequireObject(functions[i], context);
			RequireKnownKeys(object, {"name", "kind", "units", "value", "mean", "amplitude",
				"period", "phase", "file", "interpolation", "cosine", "sine"}, context);
			TemporalFunctionDefinition function;
			function.name = RequireString(Required(object, "name", context), context + ".name");
			function.kind = ParseTemporalFunctionKind(
				RequireString(Required(object, "kind", context), context + ".kind"));
			function.units = RequireString(Required(object, "units", context), context + ".units");
			if (function.name.empty() || !temporal_names.insert(function.name).second)
				throw std::runtime_error("simulation_config.json: temporal function names must be non-empty and unique");
			if (function.units.empty())
				throw std::runtime_error("simulation_config.json: temporal function units cannot be empty");
			if (const auto* value = Find(object, "value"))
				function.value = RequireNumber(*value, context + ".value");
			if (const auto* mean = Find(object, "mean"))
				function.mean = RequireNumber(*mean, context + ".mean");
			if (const auto* amplitude = Find(object, "amplitude"))
				function.amplitude = RequireNumber(*amplitude, context + ".amplitude");
			if (const auto* period = Find(object, "period"))
				function.period = RequireNumber(*period, context + ".period");
			if (const auto* phase = Find(object, "phase"))
				function.phase = RequireNumber(*phase, context + ".phase");
			if (const auto* file = Find(object, "file"))
				function.file = RequireString(*file, context + ".file");
			if (const auto* interpolation = Find(object, "interpolation"))
				function.interpolation = RequireString(*interpolation, context + ".interpolation");
			if (const auto* cosine = Find(object, "cosine"))
				function.cosine = ParseValueVector(*cosine, context + ".cosine");
			if (const auto* sine = Find(object, "sine"))
				function.sine = ParseValueVector(*sine, context + ".sine");
			if (function.kind == TemporalFunctionKind::Constant && !Find(object, "value"))
				throw std::runtime_error("simulation_config.json: constant temporal function requires value");
			if (function.kind == TemporalFunctionKind::Sinusoid
				&& (!Find(object, "mean") || !Find(object, "amplitude") || !(function.period > 0.0)))
				throw std::runtime_error("simulation_config.json: sinusoid requires mean, amplitude, and positive period");
			if (function.kind == TemporalFunctionKind::PeriodicTable
				&& (function.file.empty() || !(function.period > 0.0)
					|| function.interpolation != "linear"))
				throw std::runtime_error("simulation_config.json: periodic_table requires file, positive period, and linear interpolation");
			if (function.kind == TemporalFunctionKind::Fourier
				&& (!(function.period > 0.0) || function.cosine.size() != function.sine.size()
					|| function.cosine.empty()))
				throw std::runtime_error("simulation_config.json: fourier requires positive period and equal non-empty cosine/sine arrays");
			configuration.temporal_functions.push_back(std::move(function));
		}
	}

	const auto& systems = RequireArray(Required(root, "equation_systems", "root"), "equation_systems");
	std::set<std::string> system_names;
	for (std::size_t i = 0; i < systems.size(); ++i) {
		const auto context = "equation_systems[" + std::to_string(i) + "]";
		const auto& object = RequireObject(systems[i], context);
		RequireKnownKeys(object, {"name", "kind", "unknowns", "terms", "stabilization", "viscosity",
			"density", "time_integration"}, context);
		EquationSystemDefinition system;
		system.name = RequireString(Required(object, "name", context), context + ".name");
		system.kind = ParseEquationKind(RequireString(Required(object, "kind", context), context + ".kind"));
		if (system.name.empty() || !system_names.insert(system.name).second)
			throw std::runtime_error("simulation_config.json: equation-system names must be non-empty and unique");
		const auto& unknowns = RequireArray(Required(object, "unknowns", context), context + ".unknowns");
		std::set<std::string> local_fields;
		for (std::size_t j = 0; j < unknowns.size(); ++j) {
			const auto name = RequireString(unknowns[j], context + ".unknowns[" + std::to_string(j) + "]");
			if (!field_names.count(name)) throw std::runtime_error("simulation_config.json: unknown field '" + name + "'");
			if (!local_fields.insert(name).second) throw std::runtime_error("simulation_config.json: duplicate system unknown '" + name + "'");
			system.unknowns.push_back(name);
		}
		if (system.unknowns.empty()) throw std::runtime_error("simulation_config.json: system unknowns cannot be empty");
		if (const auto* viscosity = Find(object, "viscosity")) system.viscosity = RequireNumber(*viscosity, context + ".viscosity");
		if (const auto* density = Find(object, "density"))
			system.density = RequireNumber(*density, context + ".density");
		if (const auto* integration = Find(object, "time_integration"))
			system.time_integration = RequireString(*integration, context + ".time_integration");
		if (system.kind == EquationKind::NavierStokes && !(system.viscosity > 0.0))
			throw std::runtime_error("simulation_config.json: navier_stokes viscosity must be positive");
		if (system.kind == EquationKind::NavierStokes && !(system.density > 0.0))
			throw std::runtime_error("simulation_config.json: navier_stokes density must be positive");
		if (system.kind == EquationKind::NavierStokes
			&& system.time_integration != "steady" && system.time_integration != "backward_euler")
			throw std::runtime_error("simulation_config.json: navier_stokes time_integration must be 'steady' or 'backward_euler'");
		if (system.kind != EquationKind::NavierStokes
			&& (Find(object, "density") || Find(object, "time_integration")))
			throw std::runtime_error("simulation_config.json: density and time_integration apply only to navier_stokes");
		if (const auto* terms = Find(object, "terms")) {
			const auto& array = RequireArray(*terms, context + ".terms");
			for (std::size_t j = 0; j < array.size(); ++j) {
				const auto term_context = context + ".terms[" + std::to_string(j) + "]";
				const auto& term_object = RequireObject(array[j], term_context);
				RequireKnownKeys(term_object, {"operator", "equation", "trial", "coefficient", "velocity"}, term_context);
				WeakFormTerm term;
				term.kind = ParseTermKind(RequireString(Required(term_object, "operator", term_context), term_context + ".operator"));
				term.equation = RequireString(Required(term_object, "equation", term_context), term_context + ".equation");
				term.trial = term.equation;
				if (const auto* trial = Find(term_object, "trial")) term.trial = RequireString(*trial, term_context + ".trial");
				if (const auto* coefficient = Find(term_object, "coefficient")) term.coefficient = RequireNumber(*coefficient, term_context + ".coefficient");
				if (const auto* velocity = Find(term_object, "velocity")) term.velocity = RequireString(*velocity, term_context + ".velocity");
				if (!local_fields.count(term.equation) || !local_fields.count(term.trial))
					throw std::runtime_error("simulation_config.json: term equation and trial must be unknowns of their system");
				if (term.kind == TermKind::Advection && term.velocity.empty())
					throw std::runtime_error("simulation_config.json: advection requires a velocity source");
				system.terms.push_back(std::move(term));
			}
		}
		if (system.kind == EquationKind::LinearTransport && system.terms.empty())
			throw std::runtime_error("simulation_config.json: linear_transport requires terms");
		if (const auto* stabilization = Find(object, "stabilization")) {
			const auto& array = RequireArray(*stabilization, context + ".stabilization");
			for (std::size_t j = 0; j < array.size(); ++j) {
				const auto item_context = context + ".stabilization[" + std::to_string(j) + "]";
				const auto& item = RequireObject(array[j], item_context);
				RequireKnownKeys(item, {"equation", "method", "velocity"}, item_context);
				StabilizationDefinition definition;
				definition.equation = RequireString(Required(item, "equation", item_context), item_context + ".equation");
				definition.method = RequireString(Required(item, "method", item_context), item_context + ".method");
				definition.velocity = RequireString(Required(item, "velocity", item_context), item_context + ".velocity");
				if (!local_fields.count(definition.equation) || definition.method != "supg")
					throw std::runtime_error("simulation_config.json: stabilization requires a system equation and supported method 'supg'");
				system.stabilization.push_back(std::move(definition));
			}
		}
		configuration.equation_systems.push_back(std::move(system));
	}
	if (configuration.equation_systems.empty()) throw std::runtime_error("simulation_config.json: equation_systems cannot be empty");

	const auto& boundaries = RequireArray(Required(root, "boundaries", "root"), "boundaries");
	std::set<int> labels;
	for (std::size_t i = 0; i < boundaries.size(); ++i) {
		const auto context = "boundaries[" + std::to_string(i) + "]";
		const auto& object = RequireObject(boundaries[i], context);
		RequireKnownKeys(object, {"label", "name", "conditions"}, context);
		NamedBoundaryDefinition boundary;
		boundary.label = RequireInteger(Required(object, "label", context), context + ".label");
		boundary.name = RequireString(Required(object, "name", context), context + ".name");
		if (!labels.insert(boundary.label).second) throw std::runtime_error("simulation_config.json: duplicate boundary label");
		const auto& conditions = RequireArray(Required(object, "conditions", context), context + ".conditions");
		std::set<std::string> conditioned_fields;
		for (std::size_t j = 0; j < conditions.size(); ++j) {
			const auto condition_context = context + ".conditions[" + std::to_string(j) + "]";
			const auto& condition_object = RequireObject(conditions[j], condition_context);
			RequireKnownKeys(condition_object, {"field", "type", "value", "coefficient",
				"exterior_value", "profile", "scale", "waveform"}, condition_context);
			FieldBoundaryCondition condition;
			condition.field = RequireString(Required(condition_object, "field", condition_context), condition_context + ".field");
			condition.kind = ParseFieldBoundaryKind(RequireString(Required(condition_object, "type", condition_context), condition_context + ".type"));
			if (!field_names.count(condition.field) || !conditioned_fields.insert(condition.field).second)
				throw std::runtime_error("simulation_config.json: boundary fields must exist and be unique per boundary");
			if (const auto* value = Find(condition_object, "value")) condition.value = ParseValueVector(*value, condition_context + ".value");
			if (const auto* coefficient = Find(condition_object, "coefficient")) condition.coefficient = RequireNumber(*coefficient, condition_context + ".coefficient");
			if (const auto* exterior = Find(condition_object, "exterior_value")) condition.exterior_value = RequireNumber(*exterior, condition_context + ".exterior_value");
			if (const auto* profile = Find(condition_object, "profile")) condition.profile = RequireString(*profile, condition_context + ".profile");
			if (const auto* scale = Find(condition_object, "scale")) condition.scale = RequireNumber(*scale, condition_context + ".scale");
			if (const auto* waveform = Find(condition_object, "waveform"))
				condition.waveform = RequireString(*waveform, condition_context + ".waveform");
			if (condition.kind == FieldBoundaryKind::Dirichlet && condition.value.empty() && condition.profile.empty())
				throw std::runtime_error("simulation_config.json: dirichlet requires value or profile");
			if (!condition.value.empty() && !condition.profile.empty())
				throw std::runtime_error("simulation_config.json: boundary condition cannot set both value and profile");
			if (condition.kind == FieldBoundaryKind::Flux && condition.value.size() != 1)
				throw std::runtime_error("simulation_config.json: flux requires one value");
			if (condition.kind == FieldBoundaryKind::Robin && condition.coefficient == 0.0)
				throw std::runtime_error("simulation_config.json: robin requires nonzero coefficient");
			if (!condition.waveform.empty()
				&& (condition.kind != FieldBoundaryKind::Dirichlet
					|| !temporal_names.count(condition.waveform)))
				throw std::runtime_error("simulation_config.json: waveform must name a temporal function on a Dirichlet condition");
			boundary.conditions.push_back(std::move(condition));
		}
		configuration.boundaries.push_back(std::move(boundary));
	}
	return configuration;
}

inline SimulationConfiguration ReadSimulationConfiguration(const std::string& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open simulation configuration: " + path);
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error("cannot read simulation configuration: " + path);
	return ParseSimulationConfiguration(contents.str());
}

inline CompiledLinearSystem CompileLinearSystem(const SimulationConfiguration& configuration,
	const std::string& system_name)
{
	const EquationSystemDefinition* source = nullptr;
	for (const auto& system : configuration.equation_systems)
		if (system.name == system_name) source = &system;
	if (!source) throw std::runtime_error("unknown equation system '" + system_name + "'");
	if (source->kind != EquationKind::LinearTransport)
		throw std::runtime_error("equation system '" + system_name + "' is not linear_transport");
	CompiledLinearSystem result;
	result.name = source->name;
	result.fields = source->unknowns;
	result.dt = configuration.time.dt;
	result.steps = configuration.time.steps;
	result.stabilization = source->stabilization;
	for (std::size_t i = 0; i < result.fields.size(); ++i) result.field_index.emplace(result.fields[i], i);
	for (const auto& term : source->terms) {
		result.terms.push_back({term.kind, result.field_index.at(term.equation),
			result.field_index.at(term.trial), term.coefficient, term.velocity});
	}
	return result;
}

inline const char* TermKindName(TermKind kind)
{
	switch (kind) {
	case TermKind::TimeDerivative: return "time_derivative";
	case TermKind::Diffusion: return "diffusion";
	case TermKind::Advection: return "advection";
	case TermKind::LinearCoupling: return "linear_coupling";
	case TermKind::VolumeSource: return "volume_source";
	}
	return "unknown";
}

} // namespace iga

#endif

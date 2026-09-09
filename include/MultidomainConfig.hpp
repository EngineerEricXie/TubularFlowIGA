#ifndef IGA_MULTIDOMAIN_CONFIG_HPP
#define IGA_MULTIDOMAIN_CONFIG_HPP

#include "CheckedText.hpp"
#include "CaseConfig.hpp"
#include "FlowDomainPortMetadata.hpp"
#include "SimulationGraph.hpp"
#include "SpeciesCoupling.hpp"
#include "ZeroDFlowDomain.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class GraphExecutionKind { Explicit, Fixed, Aitken };
struct GraphTimeDefinition {
	double dt_s = 0.0;
	int steps = 0;
};

struct GraphExecutionDefinition {
	GraphExecutionKind kind = GraphExecutionKind::Explicit;
	int maximum_iterations = 1;
	double pressure_relative_tolerance = 1.0e-6;
	double pressure_reference_pa = 1.0;
	double flow_relative_tolerance = 1.0e-10;
	double relaxation_factor = 0.5;
	double minimum_relaxation = 0.05;
	double maximum_relaxation = 1.0;
	std::optional<SpeciesRoutingControls> species_routing;
	std::map<std::string, SpeciesAmountTolerance> species_amount_tolerances;
};

struct GraphDomainDefinition {
	std::string id;
	DomainKind kind = DomainKind::OneDFlow;
	std::filesystem::path case_directory;
	std::filesystem::path database;
	std::filesystem::path zero_d_model;
	std::optional<ZeroDFlowModelConfiguration> zero_d_flow_model;
	OneDInletPolicy one_d_inlet_policy = OneDInletPolicy::ConfiguredOpenLoop;
	std::vector<CouplingPort> ports;
	std::map<std::string, std::string> species_bindings;
};

struct MultidomainConfiguration {
	int schema_version = 5;
	GraphTimeDefinition time;
	std::string start_domain_id;
	GraphExecutionDefinition execution;
	std::vector<GraphDomainDefinition> domains;
	std::map<std::string, double> initial_pressure_pa;
	SimulationGraph graph;

	MultidomainConfiguration(int schema, GraphTimeDefinition time_value, std::string start_domain,
		GraphExecutionDefinition execution_value, std::vector<GraphDomainDefinition> domain_values,
		std::map<std::string, double> initial_pressures, SimulationGraph graph_value)
		: schema_version(schema), time(time_value), start_domain_id(std::move(start_domain)),
		  execution(execution_value), domains(std::move(domain_values)),
		  initial_pressure_pa(std::move(initial_pressures)), graph(std::move(graph_value)) {}
};

namespace multidomain_detail {

using config_detail::Find;
using config_detail::JsonValue;
using config_detail::RequireArray;
using config_detail::RequireInteger;
using config_detail::RequireKnownKeys;
using config_detail::RequireNumber;
using config_detail::RequireObject;
using config_detail::RequireString;

inline const JsonValue& Required(const std::map<std::string, JsonValue>& object,
	const std::string& key, const std::string& context)
{
	const auto* value = Find(object, key);
	if (!value) throw std::runtime_error("simulation_config.json: "+context
		+" requires '"+key+"'");
	return *value;
}

inline std::filesystem::path RelativePath(const JsonValue& value, const std::string& context)
{
	const auto text = RequireString(value, context);
	const std::filesystem::path path(text);
	if (text.empty() || path.is_absolute())
		throw std::runtime_error("simulation_config.json: "+context
			+" must be a nonempty relative path");
	for (const auto& component : path)
		if (component == "..")
			throw std::runtime_error("simulation_config.json: "+context
				+" cannot traverse outside the graph case");
	return path.lexically_normal();
}

inline bool PathIsWithin(const std::filesystem::path& root,
	const std::filesystem::path& candidate)
{
	const auto mismatch = std::mismatch(root.begin(), root.end(),
		candidate.begin(), candidate.end());
	return mismatch.first == root.end();
}

inline std::filesystem::path CanonicalGraphCaseRoot(
	const std::filesystem::path& graph_case_root)
{
	std::error_code error;
	const auto root = std::filesystem::canonical(
		graph_case_root.empty() ? std::filesystem::path(".") : graph_case_root, error);
	if (error || !std::filesystem::is_directory(root))
		throw std::runtime_error("graph-case root must be an existing directory");
	return root;
}

inline std::filesystem::path ResolveContainedGraphAsset(
	const std::filesystem::path& canonical_root, const std::filesystem::path& relative,
	bool require_directory, const std::string& context)
{
	std::error_code error;
	const auto resolved = std::filesystem::canonical(canonical_root/relative, error);
	if (error || !PathIsWithin(canonical_root, resolved))
		throw std::runtime_error("graph case "+context
			+" does not resolve to an existing asset inside the graph-case root");
	if (require_directory ? !std::filesystem::is_directory(resolved)
		: !std::filesystem::is_regular_file(resolved))
		throw std::runtime_error("graph case "+context+" has the wrong filesystem type");
	return resolved;
}

inline std::filesystem::path ResolveContainedCaseFile(
	const std::filesystem::path& canonical_case_directory,
	const std::filesystem::path& relative, const std::string& context)
{
	return ResolveContainedGraphAsset(canonical_case_directory, relative, false, context);
}

inline PortQuantity ParsePortQuantity(const std::string& value)
{
	if (value == "area") return PortQuantity::Area;
	if (value == "flow_rate") return PortQuantity::FlowRate;
	if (value == "mean_pressure") return PortQuantity::MeanPressure;
	if (value == "mean_normal_traction") return PortQuantity::MeanNormalTraction;
	if (value == "total_pressure") return PortQuantity::TotalPressure;
	if (value == "species_concentration") return PortQuantity::SpeciesConcentration;
	if (value == "species_flux") return PortQuantity::SpeciesFlux;
	throw std::runtime_error("simulation_config.json: unsupported port quantity '"+value+"'");
}

inline std::set<PortQuantity> ParseQuantities(const JsonValue& value,
	const std::string& context)
{
	std::set<PortQuantity> result;
	const auto& values = RequireArray(value, context);
	for (std::size_t i = 0; i < values.size(); ++i) {
		const auto quantity = ParsePortQuantity(
			RequireString(values[i], context+"["+std::to_string(i)+"]"));
		if (!result.insert(quantity).second)
			throw std::runtime_error("simulation_config.json: "+context
				+" contains a duplicate quantity");
	}
	return result;
}

inline std::set<std::string> ParseSpeciesIds(const JsonValue& value,
	const std::string& context)
{
	std::set<std::string> result;
	const auto& values = RequireArray(value, context);
	for (std::size_t i = 0; i < values.size(); ++i) {
		const auto id = RequireString(values[i], context+"["+std::to_string(i)+"]");
		if (id.empty() || !result.insert(id).second)
			throw std::runtime_error("simulation_config.json: "+context
				+" requires unique nonempty species ids");
	}
	return result;
}

inline CouplingPort ParsePort(const JsonValue& value, const std::string& domain_id,
	std::size_t index, int schema_version)
{
	const auto context = "domains."+domain_id+".ports["+std::to_string(index)+"]";
	const auto& object = RequireObject(value, context);
	auto keys = std::set<std::string>{"id", "locator_kind", "locator",
		"native_to_outward_sign", "provides", "requires"};
	if (schema_version == 6) keys.insert("species");
	RequireKnownKeys(object, keys, context);
	CouplingPort port;
	port.id = RequireString(Required(object, "id", context), context+".id");
	port.subsystem_id = domain_id;
	port.locator_kind = RequireString(Required(object, "locator_kind", context),
		context+".locator_kind");
	port.locator = RequireString(Required(object, "locator", context), context+".locator");
	if (const auto* sign = Find(object, "native_to_outward_sign")) {
		const double number = RequireNumber(*sign, context+".native_to_outward_sign");
		if (std::floor(number) != number || number < -2147483648.0
			|| number > 2147483647.0)
			throw std::runtime_error("simulation_config.json: "+context
				+".native_to_outward_sign must be an integer");
		port.orientation.native_to_outward_sign = static_cast<int>(number);
	}
	port.provides = ParseQuantities(Required(object, "provides", context),
		context+".provides");
	port.requires = ParseQuantities(Required(object, "requires", context),
		context+".requires");
	if (schema_version == 6)
		port.species = ParseSpeciesIds(Required(object, "species", context),
			context+".species");
	ValidateCouplingPort(port);
	if (schema_version == 6) {
		const bool has_species_quantity
			= port.provides.count(PortQuantity::SpeciesConcentration)
			|| port.provides.count(PortQuantity::SpeciesFlux)
			|| port.requires.count(PortQuantity::SpeciesConcentration)
			|| port.requires.count(PortQuantity::SpeciesFlux);
		if (has_species_quantity != !port.species.empty())
			throw std::runtime_error("simulation_config.json: "+context
				+" species declarations and quantities must be present together");
	}
	return port;
}

inline GraphDomainDefinition ParseDomain(const JsonValue& value, std::size_t index,
	int schema_version)
{
	const auto context = "domains["+std::to_string(index)+"]";
	const auto& object = RequireObject(value, context);
	auto keys = std::set<std::string>{"id", "dimension", "kind", "case", "database",
		"inlet_policy", "ports"};
	if (schema_version == 5) keys.insert("zero_d_model");
	if (schema_version == 6) {
		keys.insert("species_bindings");
		// Keep this key recognized solely so the explicit schema-v6 0D diagnostic
		// below is reachable instead of being masked as an unknown key.
		keys.insert("zero_d_model");
	}
	RequireKnownKeys(object, keys, context);
	GraphDomainDefinition domain;
	domain.id = RequireString(Required(object, "id", context), context+".id");
	const auto dimension = RequireString(Required(object, "dimension", context),
		context+".dimension");
	const auto kind = RequireString(Required(object, "kind", context), context+".kind");
	if (Find(object, "zero_d_model")
		&& !(schema_version == 5 && dimension == "0d" && kind == "zero_d_flow")) {
		if (schema_version == 6 && dimension == "0d" && kind == "zero_d_flow")
			throw std::runtime_error(
				"simulation_config.json: schema_version 6 does not support zero_d_flow");
		throw std::runtime_error("simulation_config.json: "+context
			+" zero_d_model is only valid for schema_version 5 zero_d_flow domains");
	}
	if (dimension == "1d" && (kind == "network_flow" || kind == "one_d_flow")) {
		domain.kind = DomainKind::OneDFlow;
		if (Find(object, "database"))
			throw std::runtime_error("simulation_config.json: "+context
				+" 1D domain does not accept database");
		const auto policy = RequireString(Required(object, "inlet_policy", context),
			context+".inlet_policy");
		if (policy == "configured_open_loop")
			domain.one_d_inlet_policy = OneDInletPolicy::ConfiguredOpenLoop;
		else if (policy == "coupled_root")
			domain.one_d_inlet_policy = OneDInletPolicy::CoupledRoot;
		else throw std::runtime_error("simulation_config.json: unsupported 1D inlet_policy '"
			+policy+"'");
	} else if (dimension == "3d"
		&& (kind == "body_fitted_iga_flow" || kind == "three_d_body_fitted_flow")) {
		domain.kind = DomainKind::ThreeDBodyFittedFlow;
		if (Find(object, "inlet_policy"))
			throw std::runtime_error("simulation_config.json: "+context
				+" 3D domain does not accept inlet_policy");
		domain.database = RelativePath(Required(object, "database", context),
			context+".database");
	} else if (dimension == "3d" && kind == "three_d_immersed_flow") {
		domain.kind = DomainKind::ThreeDImmersedFlow;
		if (Find(object, "inlet_policy"))
			throw std::runtime_error("simulation_config.json: "+context
				+" immersed 3D domain does not accept inlet_policy");
		if (Find(object, "database"))
			throw std::runtime_error("simulation_config.json: "+context
			+" immersed 3D domain does not accept body-fitted database");
	} else if (dimension == "0d" && kind == "zero_d_flow") {
		if (schema_version != 5)
			throw std::runtime_error("simulation_config.json: schema_version 6 does not support zero_d_flow");
		domain.kind = DomainKind::ZeroDFlow;
		if (Find(object, "database") || Find(object, "inlet_policy"))
			throw std::runtime_error("simulation_config.json: "+context
				+" 0D domain does not accept database or inlet_policy");
		domain.zero_d_model = RelativePath(Required(object, "zero_d_model", context),
			context+".zero_d_model");
		if (domain.zero_d_model.filename() != "zero_d_model.json")
			throw std::runtime_error("simulation_config.json: "+context
				+" zero_d_model must name zero_d_model.json");
	} else throw std::runtime_error("simulation_config.json: "+context
		+" has an inconsistent or unsupported dimension/kind");
	domain.case_directory = RelativePath(Required(object, "case", context), context+".case");
	const auto& ports = RequireArray(Required(object, "ports", context), context+".ports");
	for (std::size_t port = 0; port < ports.size(); ++port)
		domain.ports.push_back(ParsePort(ports[port], domain.id, port, schema_version));
	if (schema_version == 6) {
		const auto& bindings = RequireObject(Required(object, "species_bindings", context),
			context+".species_bindings");
		for (const auto& binding : bindings) {
			const auto field = RequireString(binding.second,
				context+".species_bindings."+binding.first);
			if (binding.first.empty() || field.empty())
				throw std::runtime_error(
					"simulation_config.json: species bindings must be nonempty");
			domain.species_bindings.emplace(binding.first, field);
		}
	}
	if (domain.kind == DomainKind::OneDFlow)
		ValidateOneDFlowDomainMetadata(domain.id, domain.ports, domain.one_d_inlet_policy);
	else if (domain.kind == DomainKind::ThreeDBodyFittedFlow)
		ValidateThreeDBodyFittedFlowDomainMetadata(domain.id, domain.ports);
	else if (domain.kind == DomainKind::ThreeDImmersedFlow)
		ValidateThreeDImmersedFlowDomainMetadata(domain.id, domain.ports);
	else {
		// The referenced model chooses the source or terminal role.  Before the
		// case file is read, accept only the two exact role-shaped port contracts.
		try { ValidateZeroDFlowDomainMetadata(domain.id, domain.ports,
			ZeroDFlowRole::SourceReservoir); }
		catch (const std::exception&) { ValidateZeroDFlowDomainMetadata(domain.id,
			domain.ports, ZeroDFlowRole::TerminalRcr); }
	}
	return domain;
}

inline PortRef ParseEndpoint(const JsonValue& value, const std::string& context)
{
	const auto& object = RequireObject(value, context);
	RequireKnownKeys(object, {"domain", "port"}, context);
	return {RequireString(Required(object, "domain", context), context+".domain"),
		RequireString(Required(object, "port", context), context+".port")};
}

inline SpeciesRoutingControls ParseSpeciesRoutingControls(const JsonValue& value)
{
	const auto& object = RequireObject(value, "execution.species_routing");
	RequireKnownKeys(object, {"flow_switch_m3_s", "flow_absolute_tolerance_m3_s",
		"flow_relative_tolerance"}, "execution.species_routing");
	SpeciesRoutingControls result;
	result.flow_switch_m3_s = RequireNumber(Required(object, "flow_switch_m3_s",
		"execution.species_routing"), "execution.species_routing.flow_switch_m3_s");
	result.flow_absolute_tolerance_m3_s = RequireNumber(Required(object,
		"flow_absolute_tolerance_m3_s", "execution.species_routing"),
		"execution.species_routing.flow_absolute_tolerance_m3_s");
	result.flow_relative_tolerance = RequireNumber(Required(object,
		"flow_relative_tolerance", "execution.species_routing"),
		"execution.species_routing.flow_relative_tolerance");
	try { result.Validate(); }
	catch (const std::exception&) {
		throw std::runtime_error("simulation_config.json: species routing controls are invalid");
	}
	return result;
}

inline std::map<std::string, SpeciesAmountTolerance> ParseSpeciesAmountTolerances(
	const JsonValue& value, const std::map<std::string, SpeciesDefinition>& species)
{
	const auto& object = RequireObject(value, "execution.species_amount_tolerances");
	if (object.empty() || object.size() != species.size())
		throw std::runtime_error("simulation_config.json: species amount tolerances require exact species coverage");
	std::map<std::string, SpeciesAmountTolerance> result;
	for (const auto& definition : species) {
		const auto found = object.find(definition.first);
		if (found == object.end())
			throw std::runtime_error("simulation_config.json: species amount tolerances are missing '"
				+definition.first+"'");
		const auto context = "execution.species_amount_tolerances."+definition.first;
		const auto& tolerance = RequireObject(found->second, context);
		RequireKnownKeys(tolerance, {"absolute_tolerance", "reference_amount",
			"relative_tolerance"}, context);
		SpeciesAmountTolerance parsed;
		parsed.absolute_tolerance = RequireNumber(Required(tolerance, "absolute_tolerance", context),
			context+".absolute_tolerance");
		parsed.reference_amount = RequireNumber(Required(tolerance, "reference_amount", context),
			context+".reference_amount");
		parsed.relative_tolerance = RequireNumber(Required(tolerance, "relative_tolerance", context),
			context+".relative_tolerance");
		try { parsed.Validate(); }
		catch (const std::exception&) {
			throw std::runtime_error("simulation_config.json: species amount tolerance is invalid for '"
				+definition.first+"'");
		}
		result.emplace(definition.first, parsed);
	}
	for (const auto& value : object)
		if (!species.count(value.first))
			throw std::runtime_error("simulation_config.json: species amount tolerance references undeclared species '"
				+value.first+"'");
	return result;
}

inline GraphExecutionDefinition ParseExecution(const JsonValue& value, int schema_version,
	const std::map<std::string, SpeciesDefinition>& species)
{
	const auto& object = RequireObject(value, "execution");
	auto keys = std::set<std::string>{"kind", "maximum_iterations", "pressure_relative_tolerance",
		"pressure_reference_pa", "flow_relative_tolerance", "relaxation_factor",
		"minimum_relaxation", "maximum_relaxation"};
	if (schema_version == 6) {
		keys.insert("species_routing");
		keys.insert("species_amount_tolerances");
	}
	RequireKnownKeys(object, keys, "execution");
	GraphExecutionDefinition result;
	const auto kind = RequireString(Required(object, "kind", "execution"), "execution.kind");
	if (kind == "explicit") result.kind = GraphExecutionKind::Explicit;
	else if (kind == "fixed") result.kind = GraphExecutionKind::Fixed;
	else if (kind == "aitken") result.kind = GraphExecutionKind::Aitken;
	else throw std::runtime_error("simulation_config.json: unsupported graph execution kind '"
		+kind+"'");
	if (const auto* item = Find(object, "maximum_iterations"))
		result.maximum_iterations = RequireInteger(*item, "execution.maximum_iterations");
	if (const auto* item = Find(object, "pressure_relative_tolerance"))
		result.pressure_relative_tolerance = RequireNumber(*item,
			"execution.pressure_relative_tolerance");
	if (const auto* item = Find(object, "pressure_reference_pa"))
		result.pressure_reference_pa = RequireNumber(*item, "execution.pressure_reference_pa");
	if (const auto* item = Find(object, "flow_relative_tolerance"))
		result.flow_relative_tolerance = RequireNumber(*item, "execution.flow_relative_tolerance");
	if (const auto* item = Find(object, "relaxation_factor"))
		result.relaxation_factor = RequireNumber(*item, "execution.relaxation_factor");
	if (const auto* item = Find(object, "minimum_relaxation"))
		result.minimum_relaxation = RequireNumber(*item, "execution.minimum_relaxation");
	if (const auto* item = Find(object, "maximum_relaxation"))
		result.maximum_relaxation = RequireNumber(*item, "execution.maximum_relaxation");
	if (schema_version == 6) {
		result.species_routing = ParseSpeciesRoutingControls(Required(object,
			"species_routing", "execution"));
		result.species_amount_tolerances = ParseSpeciesAmountTolerances(Required(object,
			"species_amount_tolerances", "execution"), species);
	}
	if (result.maximum_iterations < 1 || !(result.pressure_relative_tolerance > 0.0)
		|| !(result.pressure_reference_pa > 0.0) || !(result.flow_relative_tolerance > 0.0)
		|| !(result.relaxation_factor > 0.0) || result.relaxation_factor > 1.0
		|| !(result.minimum_relaxation > 0.0)
		|| result.minimum_relaxation > result.relaxation_factor
		|| result.relaxation_factor > result.maximum_relaxation
		|| result.maximum_relaxation > 1.0
		|| !std::isfinite(result.pressure_relative_tolerance)
		|| !std::isfinite(result.pressure_reference_pa)
		|| !std::isfinite(result.flow_relative_tolerance)
		|| !std::isfinite(result.relaxation_factor)
		|| !std::isfinite(result.minimum_relaxation)
		|| !std::isfinite(result.maximum_relaxation)
		|| (result.kind == GraphExecutionKind::Explicit && result.maximum_iterations != 1))
		throw std::runtime_error("simulation_config.json: graph execution controls are invalid");
	return result;
}

} // namespace multidomain_detail

inline ZeroDFlowModelConfiguration ParseZeroDFlowModelConfiguration(const std::string& text)
{
	using namespace multidomain_detail;
	const auto root_value = config_detail::JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "zero_d_model.json root");
	const auto role = RequireString(Required(root, "role", "zero_d_model.json root"),
		"zero_d_model.json.role");
	ZeroDFlowModelConfiguration result;
	auto& model = result.model;
	if (role == "source_reservoir") {
		RequireKnownKeys(root, {"role", "capacitance_m3_pa", "resistance_pa_s_m3",
			"initial_pressure_pa", "prescribed_flow_m3_s"}, "zero_d_model.json root");
		model.role = ZeroDFlowRole::SourceReservoir;
		model.source.capacitance_m3_pa = RequireNumber(Required(root, "capacitance_m3_pa",
			"zero_d_model.json root"), "zero_d_model.json.capacitance_m3_pa");
		model.source.resistance_pa_s_m3 = RequireNumber(Required(root, "resistance_pa_s_m3",
			"zero_d_model.json root"), "zero_d_model.json.resistance_pa_s_m3");
		model.source.prescribed_flow_m3_s = RequireNumber(Required(root, "prescribed_flow_m3_s",
			"zero_d_model.json root"), "zero_d_model.json.prescribed_flow_m3_s");
		// Initial pressure belongs to the returned committed-state seed.
		result.initial_state.stored_pressure_pa = RequireNumber(
			Required(root, "initial_pressure_pa", "zero_d_model.json root"),
			"zero_d_model.json.initial_pressure_pa");
	} else if (role == "terminal_rcr") {
		RequireKnownKeys(root, {"role", "proximal_resistance_pa_s_m3",
			"distal_resistance_pa_s_m3", "capacitance_m3_pa", "distal_pressure_pa",
			"initial_pressure_pa"}, "zero_d_model.json root");
		model.role = ZeroDFlowRole::TerminalRcr;
		model.terminal.proximal_resistance_pa_s_m3 = RequireNumber(Required(root,
			"proximal_resistance_pa_s_m3", "zero_d_model.json root"),
			"zero_d_model.json.proximal_resistance_pa_s_m3");
		model.terminal.distal_resistance_pa_s_m3 = RequireNumber(Required(root,
			"distal_resistance_pa_s_m3", "zero_d_model.json root"),
			"zero_d_model.json.distal_resistance_pa_s_m3");
		model.terminal.capacitance_m3_pa = RequireNumber(Required(root, "capacitance_m3_pa",
			"zero_d_model.json root"), "zero_d_model.json.capacitance_m3_pa");
		model.terminal.distal_pressure_pa = RequireNumber(Required(root, "distal_pressure_pa",
			"zero_d_model.json root"), "zero_d_model.json.distal_pressure_pa");
		result.initial_state.stored_pressure_pa = RequireNumber(
			Required(root, "initial_pressure_pa", "zero_d_model.json root"),
			"zero_d_model.json.initial_pressure_pa");
	} else throw std::runtime_error("zero_d_model.json: unsupported role '"+role+"'");
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(result.initial_state);
	return result;
}

inline ZeroDFlowModelConfiguration ReadZeroDFlowModelConfiguration(const std::filesystem::path& path,
	std::string& source_text)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open 0D flow model configuration: "+path.string());
	const auto text = iga::ReadCheckedText(input);
	source_text = text;
	return ParseZeroDFlowModelConfiguration(source_text);
}

inline ZeroDFlowModelConfiguration ReadZeroDFlowModelConfiguration(const std::filesystem::path& path)
{
	std::string source_text;
	return ReadZeroDFlowModelConfiguration(path, source_text);
}

inline MultidomainConfiguration ParseMultidomainConfiguration(const std::string& text)
{
	using namespace multidomain_detail;
	const auto root_value = config_detail::JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "root");
	const int schema_version = RequireInteger(Required(root, "schema_version", "root"),
		"schema_version");
	if (schema_version != 5 && schema_version != 6)
		throw std::runtime_error(
			"simulation_config.json: multidomain configuration requires schema_version 5 or 6");
	auto root_keys = std::set<std::string>{"schema_version", "time", "start_domain",
		"execution", "domains", "couplings"};
	if (schema_version == 6) root_keys.insert("species");
	RequireKnownKeys(root, root_keys, "root");
	std::vector<SpeciesDefinition> species;
	if (schema_version == 6) {
		const auto& definitions = RequireArray(Required(root, "species", "root"), "species");
		for (std::size_t i = 0; i < definitions.size(); ++i) {
			const auto context = "species["+std::to_string(i)+"]";
			const auto& object = RequireObject(definitions[i], context);
			RequireKnownKeys(object, {"id", "concentration_unit"}, context);
			species.push_back({RequireString(Required(object, "id", context), context+".id"),
				RequireString(Required(object, "concentration_unit", context),
					context+".concentration_unit")});
		}
		if (species.empty())
			throw std::runtime_error("simulation_config.json: schema_version 6 requires species");
		(void)MakeSpeciesRegistry(species);
	}
	GraphTimeDefinition time;
	const auto& time_object = RequireObject(Required(root, "time", "root"), "time");
	RequireKnownKeys(time_object, {"dt", "steps"}, "time");
	time.dt_s = RequireNumber(Required(time_object, "dt", "time"), "time.dt");
	time.steps = RequireInteger(Required(time_object, "steps", "time"), "time.steps");
	if (!(time.dt_s > 0.0) || !std::isfinite(time.dt_s) || time.steps < 1)
		throw std::runtime_error("simulation_config.json: graph time requires positive dt and steps");
	const auto start_domain = RequireString(Required(root, "start_domain", "root"),
		"start_domain");
	const auto execution = ParseExecution(Required(root, "execution", "root"), schema_version,
		MakeSpeciesRegistry(species));
	std::vector<GraphDomainDefinition> domains;
	std::vector<DomainNode> nodes;
	const auto& domain_values = RequireArray(Required(root, "domains", "root"), "domains");
	for (std::size_t i = 0; i < domain_values.size(); ++i) {
		auto domain = ParseDomain(domain_values[i], i, schema_version);
		nodes.push_back({domain.id, domain.kind, domain.ports, domain.species_bindings});
		domains.push_back(std::move(domain));
	}
	std::vector<CouplingEdge> edges;
	std::map<std::string, double> initial_pressure;
	const auto& coupling_values = RequireArray(Required(root, "couplings", "root"),
		"couplings");
	for (std::size_t i = 0; i < coupling_values.size(); ++i) {
		const auto context = "couplings["+std::to_string(i)+"]";
		const auto& object = RequireObject(coupling_values[i], context);
		auto keys = std::set<std::string>{"id", "a", "b", "mode", "initial_pressure_pa"};
		if (schema_version == 6) keys.insert("species");
		RequireKnownKeys(object, keys, context);
		CouplingEdge edge;
		edge.id = RequireString(Required(object, "id", context), context+".id");
		edge.first = ParseEndpoint(Required(object, "a", context), context+".a");
		edge.second = ParseEndpoint(Required(object, "b", context), context+".b");
		const auto mode = RequireString(Required(object, "mode", context), context+".mode");
		if (mode != "pressure_flow")
			throw std::runtime_error("simulation_config.json: unsupported coupling mode '"+mode+"'");
		edge.law = CouplingLaw::PressureFlow;
		if (schema_version == 6)
			edge.species = ParseSpeciesIds(Required(object, "species", context),
				context+".species");
		const double pressure = RequireNumber(Required(object, "initial_pressure_pa", context),
			context+".initial_pressure_pa");
		if (!initial_pressure.emplace(edge.id, pressure).second)
			throw std::runtime_error("simulation_config.json: coupling ids must be unique");
		edges.push_back(std::move(edge));
	}
	SimulationGraph graph(std::move(nodes), std::move(edges), std::move(species));
	if (graph.Domains().size() < 2 || graph.Edges().empty())
		throw std::runtime_error(
			"simulation_config.json: multidomain schema requires multiple coupled domains");
	(void)graph.Domain(start_domain);
	ValidateConnectedGraph(graph, start_domain);
	std::set<PortRef> attached_ports;
	for (const auto& edge : graph.Edges()) {
		attached_ports.insert(edge.first);
		attached_ports.insert(edge.second);
	}
	for (const auto& domain : graph.Domains())
		for (const auto& port : domain.second.ports)
			if (!port.requires.empty()
				&& !attached_ports.count({domain.first, port.id}))
				throw std::runtime_error("simulation_config.json: required port '"
					+domain.first+"."+port.id+"' is not attached to a coupling");
	return {schema_version, time, start_domain, execution, std::move(domains),
		std::move(initial_pressure), std::move(graph)};
}

// Capture the exact bytes parsed, keyed by logical input identity rather than
// machine-specific paths. The original one-argument overload remains available.
inline MultidomainConfiguration ReadMultidomainConfiguration(const std::string& path,
	std::map<std::string, std::string>& source_texts)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open multidomain simulation configuration: "+path);
	const auto text = iga::ReadCheckedText(input);
	source_texts.clear();
	auto& manifest = source_texts["graph manifest"];
	manifest = text;
	auto result = ParseMultidomainConfiguration(manifest);
	const auto root = multidomain_detail::CanonicalGraphCaseRoot(
		std::filesystem::path(path).parent_path());
	for (auto& domain : result.domains) {
		if (domain.kind != DomainKind::ZeroDFlow) continue;
		const auto case_directory = multidomain_detail::ResolveContainedGraphAsset(root,
			domain.case_directory, true, "domain '"+domain.id+"' case directory");
		const auto model_path = multidomain_detail::ResolveContainedCaseFile(case_directory,
			domain.zero_d_model, "domain '"+domain.id+"' zero_d_model");
		domain.zero_d_flow_model = ReadZeroDFlowModelConfiguration(model_path,
			source_texts["0D model "+domain.id]);
		ValidateZeroDFlowDomainMetadata(domain.id, domain.ports,
			domain.zero_d_flow_model->model.role);
	}
	// A 0D hydraulic component is deliberately narrow in this schema: one
	// compliant source starts the directed tree and RCR models terminate it.
	// Do this after model files are loaded, because port shape alone cannot
	// distinguish those two otherwise compatible one-port roles.
	std::size_t source_count = 0;
	bool has_zero_d = false;
	for (const auto& domain : result.domains) {
		if (domain.kind != DomainKind::ZeroDFlow) continue;
		has_zero_d = true;
		if (!domain.zero_d_flow_model)
			throw std::runtime_error("0D flow domain model was not loaded");
		if (domain.zero_d_flow_model->model.role == ZeroDFlowRole::SourceReservoir)
			++source_count;
	}
	if (has_zero_d) {
		if (source_count != 1)
			throw std::runtime_error("schema-v5 0D flow topology requires exactly one source reservoir");
		const auto plan = MakeAcyclicPressureFlowPlan(result.graph, result.start_domain_id);
		if (result.graph.Edges().size()+1 != result.graph.Domains().size())
			throw std::runtime_error("schema-v5 0D flow topology requires one connected acyclic tree");
		const auto start_it = std::find_if(result.domains.begin(), result.domains.end(),
			[&](const GraphDomainDefinition& domain) {
				return domain.id == result.start_domain_id;
			});
		if (start_it == result.domains.end())
			throw std::runtime_error("0D flow start domain is not declared");
		const auto& start = *start_it;
		if (start.kind != DomainKind::ZeroDFlow || !start.zero_d_flow_model
			|| start.zero_d_flow_model->model.role != ZeroDFlowRole::SourceReservoir)
			throw std::runtime_error("schema-v5 0D flow source reservoir must be the declared start domain");
		for (const auto& domain : result.domains) {
			if (domain.kind != DomainKind::ZeroDFlow) continue;
			const bool pressure_receiver = std::any_of(plan.interfaces.begin(), plan.interfaces.end(),
				[&](const PressureFlowInterfacePlan& edge) {
					return edge.pressure_receiver.domain_id == domain.id;
				});
			const bool flow_receiver = std::any_of(plan.interfaces.begin(), plan.interfaces.end(),
				[&](const PressureFlowInterfacePlan& edge) {
					return edge.flow_receiver.domain_id == domain.id;
				});
			if (domain.zero_d_flow_model->model.role == ZeroDFlowRole::SourceReservoir) {
				if (!pressure_receiver || flow_receiver || domain.id != result.start_domain_id)
					throw std::runtime_error("0D source reservoir must be the unique pressure-receiving start");
			} else if (!flow_receiver || pressure_receiver) {
				throw std::runtime_error("0D terminal RCR must be a flow-receiving terminal leaf");
			}
		}
	}
	return result;
}

inline MultidomainConfiguration ReadMultidomainConfiguration(const std::string& path)
{
	std::map<std::string, std::string> source_texts;
	return ReadMultidomainConfiguration(path, source_texts);
}

} // namespace iga

#endif

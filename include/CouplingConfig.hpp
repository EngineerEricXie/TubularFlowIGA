#ifndef IGA_COUPLING_CONFIG_HPP
#define IGA_COUPLING_CONFIG_HPP

#include "CaseConfig.hpp"
#include "VascularCoupling.hpp"

#include <filesystem>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace iga {

namespace coupling_config_detail {

using config_detail::Find;
using config_detail::JsonValue;
using config_detail::RequireArray;
using config_detail::RequireBoolean;
using config_detail::RequireInteger;
using config_detail::RequireKnownKeys;
using config_detail::RequireNumber;
using config_detail::RequireObject;
using config_detail::RequireString;

inline double OptionalNumber(const std::map<std::string, JsonValue>& object,
	const std::string& key, double fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireNumber(*value, context+"."+key) : fallback;
}

inline std::string OptionalString(const std::map<std::string, JsonValue>& object,
	const std::string& key, const std::string& fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireString(*value, context+"."+key) : fallback;
}

inline bool OptionalBoolean(const std::map<std::string, JsonValue>& object,
	const std::string& key, bool fallback, const std::string& context)
{
	const auto* value = Find(object, key);
	return value ? RequireBoolean(*value, context+"."+key) : fallback;
}

inline std::map<std::string, double> ParseNumberMap(const JsonValue* value,
	const std::string& context, bool nonnegative)
{
	std::map<std::string, double> result;
	if (!value) return result;
	const auto& object = RequireObject(*value, context);
	for (const auto& item : object) {
		const double number = RequireNumber(item.second, context+"."+item.first);
		if (nonnegative && number < 0.0)
			throw std::runtime_error("simulation_config.json: "+context
				+" values must be nonnegative");
		result.emplace(item.first, number);
	}
	return result;
}

} // namespace coupling_config_detail

inline CouplingDefinition ParseCouplingDefinition(
	const std::map<std::string, config_detail::JsonValue>& root)
{
	using namespace coupling_config_detail;
	CouplingDefinition result;
	if (const auto* scope = Find(root, "simulation_scope")) {
		const auto& object = RequireObject(*scope, "simulation_scope");
		RequireKnownKeys(object, {"mode"}, "simulation_scope");
		result.mode = ParseSimulationScopeMode(
			OptionalString(object, "mode", "flow_only", "simulation_scope"));
	}
	if (const auto* perfusate = Find(root, "perfusate")) {
		const auto& object = RequireObject(*perfusate, "perfusate");
		RequireKnownKeys(object, {"type", "oxygen_state", "hematocrit_percent",
			"oxygen_solubility_ml_dl_mmhg", "hemoglobin_g_dl", "p50_mmhg",
			"hill_exponent", "gas_molar_volume_ml_mmol"}, "perfusate");
		result.perfusate.type = OptionalString(object, "type", "rbc", "perfusate");
		if (result.perfusate.type != "rbc" && result.perfusate.type != "pfc")
			throw std::runtime_error("simulation_config.json: perfusate.type must be rbc or pfc");
		result.perfusate.oxygen_state = ParseOxygenTransportState(
			OptionalString(object, "oxygen_state", "dissolved_oxygen", "perfusate"));
		auto& oxygen = result.perfusate.oxygen;
		oxygen.hematocrit_percent = OptionalNumber(object,
			"hematocrit_percent", 0.0, "perfusate");
		oxygen.oxygen_solubility_ml_dl_mmhg = OptionalNumber(object,
			"oxygen_solubility_ml_dl_mmhg", 0.0031, "perfusate");
		oxygen.hemoglobin_g_dl = OptionalNumber(object, "hemoglobin_g_dl",
			0.34*oxygen.hematocrit_percent, "perfusate");
		oxygen.p50_mmhg = OptionalNumber(object, "p50_mmhg", 26.8, "perfusate");
		oxygen.hill_exponent = OptionalNumber(object,
			"hill_exponent", 2.7, "perfusate");
		oxygen.gas_molar_volume_ml_mmol = OptionalNumber(object,
			"gas_molar_volume_ml_mmol", 22.4, "perfusate");
		ValidateOxygenCapacity(oxygen);
		if (result.perfusate.type == "pfc"
			&& (oxygen.hematocrit_percent != 0.0 || oxygen.hemoglobin_g_dl != 0.0))
			throw std::runtime_error(
				"simulation_config.json: pfc perfusate requires zero hematocrit and hemoglobin");
	}
	if (const auto* external = Find(root, "external_circuit")) {
		const auto& object = RequireObject(*external, "external_circuit");
		RequireKnownKeys(object, {"reservoir", "pump", "oxygenator", "dialyzer",
			"infusion_rates"}, "external_circuit");
		if (const auto* reservoir = Find(object, "reservoir")) {
			const auto& item = RequireObject(*reservoir, "external_circuit.reservoir");
			RequireKnownKeys(item, {"volume_m3", "temperature_c",
				"hematocrit_loss_percent_per_hour", "species"},
				"external_circuit.reservoir");
			result.external_circuit.reservoir.enabled = true;
			result.external_circuit.reservoir.volume_m3 = OptionalNumber(item,
				"volume_m3", 0.0, "external_circuit.reservoir");
			result.external_circuit.reservoir.temperature_c = OptionalNumber(item,
				"temperature_c", 37.0, "external_circuit.reservoir");
			result.external_circuit.reservoir.hematocrit_loss_percent_per_hour
				= OptionalNumber(item, "hematocrit_loss_percent_per_hour", 0.0,
					"external_circuit.reservoir");
			result.external_circuit.reservoir.species = ParseNumberMap(
				Find(item, "species"), "external_circuit.reservoir.species", true);
		}
		if (const auto* pump = Find(object, "pump")) {
			const auto& item = RequireObject(*pump, "external_circuit.pump");
			RequireKnownKeys(item, {"mode", "flow_m3_s"}, "external_circuit.pump");
			result.external_circuit.pump.mode = OptionalString(item, "mode",
				"flow_control", "external_circuit.pump");
			result.external_circuit.pump.flow_m3_s = OptionalNumber(item,
				"flow_m3_s", 0.0, "external_circuit.pump");
		}
		if (const auto* oxygenator = Find(object, "oxygenator")) {
			const auto& item = RequireObject(*oxygenator,
				"external_circuit.oxygenator");
			RequireKnownKeys(item, {"enabled", "mode", "po2_mmhg",
				"atmospheric_pressure_mmhg", "co2_to_o2_fraction",
				"co2_solubility_mol_m3_mmhg"}, "external_circuit.oxygenator");
			auto& definition = result.external_circuit.oxygenator;
			definition.enabled = OptionalBoolean(item, "enabled", true,
				"external_circuit.oxygenator");
			definition.mode = OptionalString(item, "mode", definition.mode,
				"external_circuit.oxygenator");
			definition.po2_mmhg = OptionalNumber(item, "po2_mmhg", 0.0,
				"external_circuit.oxygenator");
			definition.atmospheric_pressure_mmhg = OptionalNumber(item,
				"atmospheric_pressure_mmhg", definition.atmospheric_pressure_mmhg,
				"external_circuit.oxygenator");
			definition.co2_to_o2_fraction = OptionalNumber(item,
				"co2_to_o2_fraction", definition.co2_to_o2_fraction,
				"external_circuit.oxygenator");
			definition.co2_solubility_mol_m3_mmhg = OptionalNumber(item,
				"co2_solubility_mol_m3_mmhg",
				definition.co2_solubility_mol_m3_mmhg,
				"external_circuit.oxygenator");
		}
		if (const auto* dialyzer = Find(object, "dialyzer")) {
			const auto& item = RequireObject(*dialyzer, "external_circuit.dialyzer");
			RequireKnownKeys(item, {"enabled", "exchange_flow_m3_s", "dialysate"},
				"external_circuit.dialyzer");
			result.external_circuit.dialyzer.enabled = OptionalBoolean(item,
				"enabled", true, "external_circuit.dialyzer");
			result.external_circuit.dialyzer.exchange_flow_m3_s = OptionalNumber(item,
				"exchange_flow_m3_s", 0.0, "external_circuit.dialyzer");
			result.external_circuit.dialyzer.dialysate = ParseNumberMap(
				Find(item, "dialysate"), "external_circuit.dialyzer.dialysate", true);
		}
		result.external_circuit.infusion_rates = ParseNumberMap(
			Find(object, "infusion_rates"), "external_circuit.infusion_rates", false);
	}
	if (const auto* coupling = Find(root, "coupling")) {
		const auto& object = RequireObject(*coupling, "coupling");
		RequireKnownKeys(object, {"scheme", "replay_file", "flow_epsilon_m3_s",
			"three_d_ports"},
			"coupling");
		result.scheme = OptionalString(object, "scheme", result.scheme, "coupling");
		result.replay_file = OptionalString(object, "replay_file", "", "coupling");
		result.flow_epsilon_m3_s = OptionalNumber(object, "flow_epsilon_m3_s",
			result.flow_epsilon_m3_s, "coupling");
		if (const auto* ports = Find(object, "three_d_ports")) {
			const auto& item = RequireObject(*ports, "coupling.three_d_ports");
			RequireKnownKeys(item, {"inlet_label", "outlet_labels"},
				"coupling.three_d_ports");
			const auto* inlet = Find(item, "inlet_label");
			const auto* outlets = Find(item, "outlet_labels");
			if (!inlet || !outlets)
				throw std::runtime_error("simulation_config.json: coupling.three_d_ports requires inlet_label and outlet_labels");
			result.three_d_ports.inlet_label = RequireInteger(*inlet,
				"coupling.three_d_ports.inlet_label");
			const auto& array = RequireArray(*outlets,
				"coupling.three_d_ports.outlet_labels");
			std::set<int> labels;
			for (std::size_t index = 0; index < array.size(); ++index) {
				const int label = RequireInteger(array[index],
					"coupling.three_d_ports.outlet_labels["+std::to_string(index)+"]");
				if (!labels.insert(label).second)
					throw std::runtime_error("simulation_config.json: coupling.three_d_ports outlet labels must be unique");
				result.three_d_ports.outlet_labels.push_back(label);
			}
			if (result.three_d_ports.outlet_labels.empty()
				|| labels.count(result.three_d_ports.inlet_label))
				throw std::runtime_error("simulation_config.json: coupling.three_d_ports requires distinct inlet and nonempty outlets");
		}
	}
	if (result.scheme != "explicit_staggered")
		throw std::runtime_error(
			"simulation_config.json: native coupling currently supports explicit_staggered");
	if (!(result.flow_epsilon_m3_s >= 0.0))
		throw std::runtime_error(
			"simulation_config.json: coupling flow epsilon must be nonnegative");
	if (result.mode == SimulationScopeMode::VcaReplay) {
		if (result.replay_file.empty()
			|| std::filesystem::path(result.replay_file).is_absolute())
			throw std::runtime_error(
				"simulation_config.json: vca_replay requires a relative coupling.replay_file");
	}
	if (result.mode == SimulationScopeMode::VcaClosedLoop) {
		if (!result.external_circuit.reservoir.enabled)
			throw std::runtime_error(
				"simulation_config.json: vca_closed_loop requires external_circuit.reservoir");
		if (result.external_circuit.pump.mode != "flow_control"
			|| !(result.external_circuit.pump.flow_m3_s > 0.0))
			throw std::runtime_error(
				"simulation_config.json: vca_closed_loop requires a positive flow-control pump");
	}
	return result;
}

} // namespace iga

#endif

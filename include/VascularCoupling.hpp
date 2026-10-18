#ifndef IGA_VASCULAR_COUPLING_HPP
#define IGA_VASCULAR_COUPLING_HPP

#include "OxygenCapacity.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class SimulationScopeMode {
	FlowOnly, VascularOpenLoop, VcaReplay, VcaClosedLoop
};

inline SimulationScopeMode ParseSimulationScopeMode(const std::string& value)
{
	if (value == "flow_only") return SimulationScopeMode::FlowOnly;
	if (value == "vascular_open_loop") return SimulationScopeMode::VascularOpenLoop;
	if (value == "vca_replay") return SimulationScopeMode::VcaReplay;
	if (value == "vca_closed_loop") return SimulationScopeMode::VcaClosedLoop;
	throw std::runtime_error("unsupported simulation_scope.mode '"+value+"'");
}

inline const char* SimulationScopeModeName(SimulationScopeMode value)
{
	if (value == SimulationScopeMode::VascularOpenLoop) return "vascular_open_loop";
	if (value == SimulationScopeMode::VcaReplay) return "vca_replay";
	if (value == SimulationScopeMode::VcaClosedLoop) return "vca_closed_loop";
	return "flow_only";
}

struct VascularInletState {
	double time_s = 0.0;
	bool has_flow = false;
	double flow_m3_s = 0.0;
	bool has_pressure = false;
	double pressure_pa = 0.0;
	std::map<std::string, double> species;
	bool has_temperature = false;
	double temperature_c = 0.0;
	bool has_hematocrit = false;
	double hematocrit_percent = 0.0;
	std::map<std::string, std::string> metadata;
};

struct VascularOutletState {
	int outlet_id = -1;
	double flow_m3_s = 0.0;
	double pressure_pa = 0.0;
	std::map<std::string, double> species_flux;
	std::map<std::string, double> flux_weighted_concentration;
	std::map<std::string, double> derived;
	bool average_valid = true;
};

struct AggregatedVascularReturn {
	double flow_m3_s = 0.0;
	bool pressure_valid = false;
	double pressure_pa = 0.0;
	std::map<std::string, double> species_flux;
	std::map<std::string, double> flux_weighted_concentration;
	bool average_valid = true;
	std::vector<int> outlet_ids;
};

struct VascularStepResult {
	double time_s = 0.0;
	double dt_s = 0.0;
	VascularInletState inlet;
	std::vector<VascularOutletState> outlets;
	std::map<std::string, double> total_mass;
	std::map<std::string, double> source_integrals;
	std::map<std::string, double> balance_residuals;
	std::map<std::string, double> hemodynamics;
};

inline void RequireFiniteCouplingValue(const std::string& name, double value)
{
	if (!std::isfinite(value)) throw std::runtime_error(name+" must be finite");
}

inline void ValidateVascularInletState(const VascularInletState& state)
{
	RequireFiniteCouplingValue("vascular inlet time", state.time_s);
	if (!state.has_flow && !state.has_pressure)
		throw std::runtime_error("vascular inlet requires flow or pressure");
	if (state.has_flow) RequireFiniteCouplingValue("vascular inlet flow", state.flow_m3_s);
	if (state.has_pressure) RequireFiniteCouplingValue("vascular inlet pressure", state.pressure_pa);
	if (state.has_temperature)
		RequireFiniteCouplingValue("vascular inlet temperature", state.temperature_c);
	if (state.has_hematocrit) {
		RequireFiniteCouplingValue("vascular inlet hematocrit", state.hematocrit_percent);
		if (state.hematocrit_percent < 0.0 || state.hematocrit_percent > 100.0)
			throw std::runtime_error("vascular inlet hematocrit must be in [0,100]");
	}
	for (const auto& item : state.species) {
		RequireFiniteCouplingValue("vascular inlet species "+item.first, item.second);
		if (item.second < 0.0)
			throw std::runtime_error("vascular inlet concentrations must be nonnegative");
	}
}

inline void ValidateVascularOutletState(const VascularOutletState& state)
{
	if (state.outlet_id < 0) throw std::runtime_error("vascular outlet id must be nonnegative");
	RequireFiniteCouplingValue("vascular outlet flow", state.flow_m3_s);
	RequireFiniteCouplingValue("vascular outlet pressure", state.pressure_pa);
	for (const auto& item : state.species_flux)
		RequireFiniteCouplingValue("vascular outlet species flux "+item.first, item.second);
	for (const auto& item : state.flux_weighted_concentration) {
		RequireFiniteCouplingValue("vascular outlet concentration "+item.first, item.second);
		if (item.second < 0.0)
			throw std::runtime_error("vascular outlet concentrations must be nonnegative");
	}
}

inline AggregatedVascularReturn AggregateVascularOutlets(
	const std::vector<VascularOutletState>& outlets, double flow_epsilon_m3_s = 1.0e-14)
{
	if (outlets.empty()) throw std::runtime_error("vascular aggregation requires an outlet");
	if (!(flow_epsilon_m3_s >= 0.0) || !std::isfinite(flow_epsilon_m3_s))
		throw std::runtime_error("flow aggregation epsilon must be finite and nonnegative");
	AggregatedVascularReturn result;
	std::map<int, bool> ids;
	std::map<std::string, double> concentration_numerator;
	for (const auto& outlet : outlets) {
		ValidateVascularOutletState(outlet);
		if (!ids.emplace(outlet.outlet_id, true).second)
			throw std::runtime_error("vascular outlet ids must be unique");
		result.outlet_ids.push_back(outlet.outlet_id);
		result.flow_m3_s += outlet.flow_m3_s;
		result.pressure_pa += outlet.flow_m3_s*outlet.pressure_pa;
		for (const auto& item : outlet.species_flux)
			result.species_flux[item.first] += item.second;
		for (const auto& item : outlet.species_flux) {
			const auto concentration = outlet.flux_weighted_concentration.find(item.first);
			concentration_numerator[item.first] += concentration
				== outlet.flux_weighted_concentration.end()
				? item.second : outlet.flow_m3_s*concentration->second;
		}
	}
	result.average_valid = std::abs(result.flow_m3_s) > flow_epsilon_m3_s;
	result.pressure_valid = result.average_valid;
	if (result.average_valid) {
		result.pressure_pa /= result.flow_m3_s;
		for (const auto& item : concentration_numerator)
			result.flux_weighted_concentration[item.first]
				= item.second/result.flow_m3_s;
	} else result.pressure_pa = 0.0;
	return result;
}

struct PerfusateDefinition {
	std::string type = "rbc";
	OxygenTransportState oxygen_state = OxygenTransportState::Dissolved;
	OxygenCapacityParameters oxygen;
};

struct ReservoirDefinition {
	bool enabled = false;
	double volume_m3 = 0.0;
	double temperature_c = 37.0;
	double hematocrit_loss_percent_per_hour = 0.0;
	std::map<std::string, double> species;
};

struct PumpDefinition {
	// flow_control: the pump sets flow_m3_s. pressure_control: the pump holds
	// pressure_pa at the vascular inlet and the vascular model returns the flow;
	// at steady state the pump returns the venous flow to the vascular bed.
	std::string mode = "flow_control";
	double flow_m3_s = 0.0;
	double pressure_pa = 0.0;
};

struct OxygenatorDefinition {
	bool enabled = false;
	std::string mode = "full_equilibration";
	double po2_mmhg = 0.0;
	double atmospheric_pressure_mmhg = 760.0;
	double co2_to_o2_fraction = 0.053;
	double co2_solubility_mol_m3_mmhg = 0.0301;
};

struct DialyzerDefinition {
	bool enabled = false;
	double exchange_flow_m3_s = 0.0;          // volume-neutral exchange
	double inflow_m3_s = 0.0;                 // dialysate added to the reservoir
	double outflow_m3_s = 0.0;                // reservoir perfusate removed
	std::map<std::string, double> dialysate;
};

struct ExternalCircuitDefinition {
	ReservoirDefinition reservoir;
	PumpDefinition pump;
	OxygenatorDefinition oxygenator;
	DialyzerDefinition dialyzer;
	std::map<std::string, double> infusion_rates;
};

struct ThreeDVascularPortDefinition {
	int inlet_label = -1;
	std::vector<int> outlet_labels;
};

struct CouplingDefinition {
	SimulationScopeMode mode = SimulationScopeMode::FlowOnly;
	std::string scheme = "explicit_staggered";
	std::string replay_file;
	double flow_epsilon_m3_s = 1.0e-14;
	ThreeDVascularPortDefinition three_d_ports;
	PerfusateDefinition perfusate;
	ExternalCircuitDefinition external_circuit;
};

inline void RequireFlowOnlyCoupling(const CouplingDefinition& definition,
	const std::string& backend)
{
	if (definition.mode != SimulationScopeMode::FlowOnly)
		throw std::runtime_error(backend+" does not yet implement VCA vascular coupling; "
			+"use the native 1d iga_1d solver or select simulation_scope.mode=flow_only");
}

struct ReservoirState {
	double volume_m3 = 0.0;
	std::map<std::string, double> species;
	double temperature_c = 37.0;
	double hematocrit_percent = 0.0;
};

struct CircuitAdvanceReport {
	double time_s = 0.0;
	double volume_change_m3 = 0.0;
	std::map<std::string, double> species_mass_change;
	std::map<std::string, double> device_source_rate;
};

class VcaExternalCircuit {
public:
	explicit VcaExternalCircuit(const CouplingDefinition& definition)
		: definition_(definition)
	{
		if (definition_.mode != SimulationScopeMode::VcaClosedLoop)
			throw std::runtime_error("VCA external circuit requires vca_closed_loop mode");
		const auto& pump = definition_.external_circuit.pump;
		if (pump.mode == "flow_control") {
			if (!(pump.flow_m3_s > 0.0))
				throw std::runtime_error("VCA flow-control pump requires positive flow_m3_s");
		} else if (pump.mode == "pressure_control") {
			if (!std::isfinite(pump.pressure_pa))
				throw std::runtime_error("VCA pressure-control pump requires a finite pressure_pa");
		} else throw std::runtime_error("VCA pump mode must be flow_control or pressure_control");
		const auto& dialyzer = definition_.external_circuit.dialyzer;
		if (!(dialyzer.exchange_flow_m3_s >= 0.0) || !(dialyzer.inflow_m3_s >= 0.0)
			|| !(dialyzer.outflow_m3_s >= 0.0))
			throw std::runtime_error("VCA dialyzer flows must be nonnegative");
		const auto& reservoir = definition_.external_circuit.reservoir;
		if (!reservoir.enabled || !(reservoir.volume_m3 > 0.0))
			throw std::runtime_error("VCA closed loop requires a positive reservoir volume");
		state_.volume_m3 = reservoir.volume_m3;
		state_.species = reservoir.species;
		state_.temperature_c = reservoir.temperature_c;
		state_.hematocrit_percent = definition_.perfusate.oxygen.hematocrit_percent;
		ValidateState();
	}

	const ReservoirState& State() const { return state_; }

	void RestoreState(const ReservoirState& state)
	{
		state_ = state;
		ValidateState();
	}

	// Mutable controls for a perfusion controller acting between steps.
	ExternalCircuitDefinition& Controls() { return definition_.external_circuit; }

	// Adds solute amounts [mol, or the species' amount unit] and an optional
	// carrier volume to the well-mixed reservoir.
	void ApplyBolus(const std::map<std::string, double>& amount, double volume_m3 = 0.0)
	{
		if (!(volume_m3 >= 0.0) || !std::isfinite(volume_m3))
			throw std::runtime_error("VCA bolus volume must be finite and nonnegative");
		// Validate every amount first so a rejected bolus leaves the state untouched.
		for (const auto& item : amount)
			if (!(item.second >= 0.0) || !std::isfinite(item.second))
				throw std::runtime_error("VCA bolus amounts must be finite and nonnegative");
		const double new_volume = state_.volume_m3+volume_m3;
		for (auto& item : state_.species) item.second *= state_.volume_m3/new_volume;
		for (const auto& item : amount) state_.species[item.first] += item.second/new_volume;
		state_.hematocrit_percent *= state_.volume_m3/new_volume;
		state_.volume_m3 = new_volume;
		ValidateState();
	}

	VascularInletState InletState(double time_s)
	{
		VascularInletState result;
		result.time_s = time_s;
		if (definition_.external_circuit.pump.mode == "pressure_control") {
			result.has_pressure = true;
			result.pressure_pa = definition_.external_circuit.pump.pressure_pa;
		} else {
			result.has_flow = true;
			result.flow_m3_s = definition_.external_circuit.pump.flow_m3_s;
		}
		result.species = state_.species;
		result.has_temperature = true;
		result.temperature_c = state_.temperature_c;
		result.has_hematocrit = true;
		result.hematocrit_percent = state_.hematocrit_percent;
		const auto& oxygenator = definition_.external_circuit.oxygenator;
		if (oxygenator.enabled) {
			if (oxygenator.mode != "full_equilibration")
				throw std::runtime_error("unsupported oxygenator mode '"+oxygenator.mode+"'");
			double po2 = oxygenator.po2_mmhg;
			if (!(po2 > 0.0)) po2 = oxygenator.atmospheric_pressure_mmhg;
			const auto equilibrium = OxygenFromPo2(po2, OxygenParameters());
			if (definition_.perfusate.oxygen_state == OxygenTransportState::Total)
				result.species["total_oxygen"] = equilibrium.total_oxygen_mol_m3;
			else result.species["oxygen"] = equilibrium.dissolved_oxygen_mol_m3;
			if (result.species.count("carbon_dioxide"))
				result.species["carbon_dioxide"] = po2*oxygenator.co2_to_o2_fraction
					*oxygenator.co2_solubility_mol_m3_mmhg;
		}
		last_arterial_species_ = result.species;
		ValidateVascularInletState(result);
		return result;
	}

	CircuitAdvanceReport Advance(const AggregatedVascularReturn& venous,
		double dt_s, double time_s)
	{
		if (!(dt_s >= 0.0) || !std::isfinite(dt_s))
			throw std::runtime_error("VCA timestep must be finite and nonnegative");
		if (definition_.external_circuit.pump.mode != "flow_control")
			throw std::runtime_error("VCA Advance requires a flow-control pump; use AdvanceWithUptake");
		const auto& volume_dialyzer = definition_.external_circuit.dialyzer;
		if (volume_dialyzer.enabled && (volume_dialyzer.inflow_m3_s > 0.0 || volume_dialyzer.outflow_m3_s > 0.0))
			throw std::runtime_error("VCA Advance supports volume-neutral dialysis only; use AdvanceWithUptake");
		const double old_volume = state_.volume_m3;
		const double pump_flow = definition_.external_circuit.pump.flow_m3_s;
		const double new_volume = old_volume+dt_s*(venous.flow_m3_s-pump_flow);
		if (!(new_volume > 0.0) || !std::isfinite(new_volume))
			throw std::runtime_error("VCA reservoir volume became non-positive");
		std::map<std::string, bool> names;
		for (const auto& item : state_.species) names[item.first] = true;
		for (const auto& item : venous.species_flux) names[item.first] = true;
		for (const auto& item : definition_.external_circuit.infusion_rates)
			names[item.first] = true;
		CircuitAdvanceReport report;
		report.time_s = time_s;
		report.volume_change_m3 = new_volume-old_volume;
		std::map<std::string, double> next;
		for (const auto& named : names) {
			const auto& name = named.first;
			const double concentration = state_.species.count(name)
				? state_.species.at(name) : 0.0;
			const double venous_flux = venous.species_flux.count(name)
				? venous.species_flux.at(name) : 0.0;
			const double infusion = definition_.external_circuit.infusion_rates.count(name)
				? definition_.external_circuit.infusion_rates.at(name) : 0.0;
			double mass = old_volume*concentration
				+dt_s*(venous_flux-pump_flow*concentration+infusion);
			if (mass < -1.0e-12 || !std::isfinite(mass))
				throw std::runtime_error("VCA reservoir species mass became negative for '"+name+"'");
			mass = std::max(mass, 0.0);
			next[name] = mass/new_volume;
			report.species_mass_change[name] = mass-old_volume*concentration;
			const double arterial = last_arterial_species_.count(name)
				? last_arterial_species_.at(name) : concentration;
			report.device_source_rate[name] = infusion
				+pump_flow*(arterial-concentration);
		}
		const auto& dialyzer = definition_.external_circuit.dialyzer;
		if (dialyzer.enabled && dialyzer.exchange_flow_m3_s > 0.0 && dt_s > 0.0) {
			for (auto& item : next) {
				const double target = dialyzer.dialysate.count(item.first)
					? dialyzer.dialysate.at(item.first) : item.second;
				const double delta = dialyzer.exchange_flow_m3_s*dt_s*(target-item.second);
				const double mass = new_volume*item.second+delta;
				if (mass < -1.0e-12)
					throw std::runtime_error("dialyzer produced negative mass for '"+item.first+"'");
				item.second = std::max(mass, 0.0)/new_volume;
				report.species_mass_change[item.first] += delta;
				report.device_source_rate[item.first] += delta/dt_s;
			}
		}
		const double rbc_volume = old_volume*state_.hematocrit_percent;
		state_.hematocrit_percent = std::max(0.0,
			rbc_volume/new_volume
			-definition_.external_circuit.reservoir.hematocrit_loss_percent_per_hour
				*dt_s/3600.0);
		state_.volume_m3 = new_volume;
		state_.species = std::move(next);
		ValidateState();
		return report;
	}

	// Advances the reservoir by dt_s with the organ treated as a quasi-steady
	// load. `organ_uptake` is the net removal rate by the vascular bed per
	// species (positive = consumed), evaluated at the current reservoir state;
	// `arterial` holds species that a gas exchanger resets on every pass
	// (O2, CO2), which relax towards the arterial value at the circuit flow.
	// Uptake is linearised as first-order clearance (k = U/c) and production as
	// a constant source; together with dialysis, infusion and the exchanger
	// this gives a linear ODE per species that is integrated exactly, so the
	// update is stable for any step when circuit transit is fast. Dialysate
	// inflow dilutes; perfusate outflow removes volume at reservoir composition.
	CircuitAdvanceReport AdvanceWithUptake(const std::map<std::string, double>& organ_uptake,
		const std::map<std::string, double>& arterial, double circuit_flow_m3_s,
		double dt_s, double time_s)
	{
		if (!(dt_s >= 0.0) || !std::isfinite(dt_s) || !(circuit_flow_m3_s >= 0.0)
			|| !std::isfinite(circuit_flow_m3_s))
			throw std::runtime_error("VCA uptake step requires finite nonnegative dt and flow");
		const auto& circuit = definition_.external_circuit;
		const auto& dialyzer = circuit.dialyzer;
		const double inflow = dialyzer.enabled ? dialyzer.inflow_m3_s : 0.0;
		const double outflow = dialyzer.enabled ? dialyzer.outflow_m3_s : 0.0;
		const double exchange = dialyzer.enabled ? dialyzer.exchange_flow_m3_s : 0.0;
		const double old_volume = state_.volume_m3;
		const double new_volume = old_volume+dt_s*(inflow-outflow);
		if (!(new_volume > 0.0) || !std::isfinite(new_volume))
			throw std::runtime_error("VCA reservoir volume became non-positive");
		std::map<std::string, bool> names;
		for (const auto& item : state_.species) names[item.first] = true;
		for (const auto& item : organ_uptake) names[item.first] = true;
		for (const auto& item : circuit.infusion_rates) names[item.first] = true;
		CircuitAdvanceReport report;
		report.time_s = time_s;
		report.volume_change_m3 = new_volume-old_volume;
		std::map<std::string, double> next;
		for (const auto& named : names) {
			const auto& name = named.first;
			const double c0 = state_.species.count(name) ? state_.species.at(name) : 0.0;
			const double uptake = organ_uptake.count(name) ? organ_uptake.at(name) : 0.0;
			RequireFiniteCouplingValue("VCA organ uptake "+name, uptake);
			const double infusion = circuit.infusion_rates.count(name) ? circuit.infusion_rates.at(name) : 0.0;
			const double dialysate = dialyzer.dialysate.count(name) ? dialyzer.dialysate.at(name) : c0;
			// V dc/dt = -alpha c + beta (outflow at reservoir composition leaves c unchanged).
			double alpha = inflow+exchange, beta = infusion+(inflow+exchange)*dialysate;
			if (uptake > 0.0 && c0 > 0.0) alpha += uptake/c0;
			else beta -= uptake;
			const auto reset = arterial.find(name);
			if (reset != arterial.end()) {
				alpha += circuit_flow_m3_s;
				beta += circuit_flow_m3_s*reset->second;
			}
			// V(t) = V0 + r t; dc/dt = (beta - alpha c)/V(t) has the exact solution
			// c - c_inf = (c0 - c_inf) (V(t)/V0)^(-alpha/r), or exp(-alpha t/V0) at r = 0.
			const double rate = inflow-outflow;
			double c1 = c0;
			if (alpha > 0.0) {
				const double target = beta/alpha;
				const double decay = std::abs(rate) > 1.0e-14*old_volume/std::max(dt_s, 1.0e-300)
					? std::pow(new_volume/old_volume, -alpha/rate)
					: std::exp(-alpha*dt_s/old_volume);
				c1 = target+(c0-target)*decay;
			} else c1 = c0+beta*(std::abs(rate) > 0.0
				? std::log(new_volume/old_volume)/rate : dt_s/old_volume);
			c1 = std::max(c1, 0.0);
			next[name] = c1;
			report.species_mass_change[name] = new_volume*c1-old_volume*c0;
			report.device_source_rate[name] = dt_s > 0.0
				? (new_volume*c1-old_volume*c0)/dt_s+uptake : 0.0;
		}
		const double rbc_volume = old_volume*state_.hematocrit_percent;
		state_.hematocrit_percent = std::max(0.0, rbc_volume/new_volume
			-circuit.reservoir.hematocrit_loss_percent_per_hour*dt_s/3600.0);
		state_.volume_m3 = new_volume;
		state_.species = std::move(next);
		ValidateState();
		return report;
	}

private:
	OxygenCapacityParameters OxygenParameters() const
	{
		auto parameters = definition_.perfusate.oxygen;
		const double baseline_hematocrit = parameters.hematocrit_percent;
		const double hemoglobin_per_hematocrit = baseline_hematocrit > 0.0
			? parameters.hemoglobin_g_dl/baseline_hematocrit : 0.34;
		parameters.hematocrit_percent = state_.hematocrit_percent;
		parameters.hemoglobin_g_dl = state_.hematocrit_percent
			*hemoglobin_per_hematocrit;
		return parameters;
	}

	void ValidateState() const
	{
		if (!(state_.volume_m3 > 0.0) || !std::isfinite(state_.volume_m3))
			throw std::runtime_error("reservoir volume must be positive and finite");
		if (!std::isfinite(state_.temperature_c)
			|| !std::isfinite(state_.hematocrit_percent)
			|| state_.hematocrit_percent < 0.0 || state_.hematocrit_percent > 100.0)
			throw std::runtime_error("reservoir metadata is invalid");
		for (const auto& item : state_.species)
			if (!std::isfinite(item.second) || item.second < 0.0)
				throw std::runtime_error("reservoir concentrations must be finite and nonnegative");
	}

	CouplingDefinition definition_;
	ReservoirState state_;
	std::map<std::string, double> last_arterial_species_;
};

} // namespace iga

#endif

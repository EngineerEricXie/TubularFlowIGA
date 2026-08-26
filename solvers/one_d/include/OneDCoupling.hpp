#ifndef IGA_ONE_D_COUPLING_HPP
#define IGA_ONE_D_COUPLING_HPP

#include "OneDTransport.hpp"
#include "VascularCoupling.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

inline std::string CouplingEscapeJson(const std::string& value)
{
	std::string result;
	for (const char character : value) {
		if (character == '\\' || character == '"') result.push_back('\\');
		result.push_back(character);
	}
	return result;
}

inline void WriteCouplingNumberMap(std::ostream& output,
	const std::map<std::string, double>& values)
{
	output << '{';
	bool first = true;
	for (const auto& item : values) {
		output << (first ? "" : ", ") << '"' << CouplingEscapeJson(item.first)
			<< "\": " << item.second;
		first = false;
	}
	output << '}';
}

inline void WriteCouplingOutletIds(std::ostream& output,
	const std::vector<int>& outlet_ids)
{
	output << '[';
	for (std::size_t i = 0; i < outlet_ids.size(); ++i)
		output << (i == 0 ? "" : ", ") << outlet_ids[i];
	output << ']';
}

inline double ApplyOneDCoupledInlet(OneDConfiguration& configuration,
	std::vector<OneDTransportState>& transports, const VascularInletState& inlet)
{
	ValidateVascularInletState(inlet);
	if (!inlet.has_flow)
		throw std::runtime_error("native 1d coupling currently requires a flow-controlled inlet");
	if (inlet.has_hematocrit) {
		const double prior_hematocrit
			= configuration.coupling.perfusate.oxygen.hematocrit_percent;
		const double hemoglobin_per_hematocrit = prior_hematocrit > 0.0
			? configuration.coupling.perfusate.oxygen.hemoglobin_g_dl/prior_hematocrit
			: 0.34;
		configuration.physiology.hematocrit_percent = inlet.hematocrit_percent;
		configuration.physiology.hemoglobin_g_dl = inlet.hematocrit_percent
			*hemoglobin_per_hematocrit;
		configuration.coupling.perfusate.oxygen.hematocrit_percent
			= inlet.hematocrit_percent;
		configuration.coupling.perfusate.oxygen.hemoglobin_g_dl
			= inlet.hematocrit_percent*hemoglobin_per_hematocrit;
	}
	for (const auto& item : inlet.species) {
		auto* species = FindOneDSpecies(transports, item.first);
		if (!species)
			throw std::runtime_error("coupled inlet species '"+item.first
				+"' is not transported by the selected 1d system");
		species->inlet_value = item.second;
		species->inlet_waveform.clear();
	}
	return inlet.flow_m3_s;
}

inline std::map<std::string, double> OneDTotalSpeciesMass(
	const OneDNetwork& network, const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports)
{
	std::map<std::string, double> result;
	for (const auto& transport : transports)
		for (const auto& species : transport.species) {
			double mass = 0.0;
			for (const auto& segment : network.segments) {
				const double dx = segment.length/segment.cells;
				for (int cell = 0; cell < segment.cells; ++cell) {
					const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
					mass += flow.area[index]*species.concentration[index]*dx;
				}
			}
			result.emplace(species.definition.field, mass);
		}
	return result;
}

inline VascularStepResult BuildOneDStepResult(const OneDConfiguration& configuration,
	const OneDNetwork& network, const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports,
	const VascularInletState& inlet, double dt_s,
	const std::map<std::string, double>& previous_mass = {})
{
	VascularStepResult result;
	result.time_s = inlet.time_s;
	result.dt_s = dt_s;
	result.inlet = inlet;
	result.total_mass = OneDTotalSpeciesMass(network, flow, transports);
	for (const auto& transport : transports)
		for (const auto& species : transport.species)
			result.source_integrals.emplace(species.definition.field,
				OneDSpeciesSourceIntegral(configuration, network, flow, species));
	for (const auto& outlet : flow.outlets) {
		const int incoming = OneDSegmentIntoNode(network, outlet.node);
		if (incoming < 0) throw std::runtime_error("1d outlet does not have an incoming segment");
		const auto& segment = network.segments[static_cast<std::size_t>(incoming)];
		const auto cell = static_cast<std::size_t>(segment.cell_offset+segment.cells-1);
		VascularOutletState port;
		port.outlet_id = network.nodes[static_cast<std::size_t>(outlet.node)].id;
		port.flow_m3_s = outlet.flow;
		port.pressure_pa = outlet.pressure;
		port.average_valid = std::abs(port.flow_m3_s)
			> configuration.coupling.flow_epsilon_m3_s;
		for (const auto& transport : transports)
			for (const auto& species : transport.species) {
				const double concentration = species.concentration[cell];
				port.species_flux[species.definition.field]
					= port.flow_m3_s*concentration;
				if (port.average_valid)
					port.flux_weighted_concentration[species.definition.field]
						= concentration;
			}
		ValidateVascularOutletState(port);
		result.outlets.push_back(std::move(port));
	}
	double outlet_flow = 0.0;
	double outlet_pressure_flow = 0.0;
	for (const auto& outlet : result.outlets) {
		outlet_flow += outlet.flow_m3_s;
		outlet_pressure_flow += outlet.flow_m3_s*outlet.pressure_pa;
	}
	const double outlet_pressure = std::abs(outlet_flow)
		> configuration.coupling.flow_epsilon_m3_s
		? outlet_pressure_flow/outlet_flow : 0.0;
	const double inlet_pressure = flow.node_pressure.empty() ? 0.0
		: flow.node_pressure[static_cast<std::size_t>(network.root)];
	result.hemodynamics["inlet_flow_m3_s"] = flow.inlet_flow;
	result.hemodynamics["outlet_flow_m3_s"] = outlet_flow;
	result.hemodynamics["inlet_pressure_pa"] = inlet_pressure;
	result.hemodynamics["outlet_pressure_pa"] = outlet_pressure;
	if (std::abs(flow.inlet_flow) > configuration.coupling.flow_epsilon_m3_s)
		result.hemodynamics["vascular_resistance_pa_s_m3"]
			= (inlet_pressure-outlet_pressure)/flow.inlet_flow;
	if (dt_s > 0.0) {
		for (const auto& mass : result.total_mass) {
			const auto old = previous_mass.find(mass.first);
			if (old == previous_mass.end()) continue;
			const auto* species = FindOneDSpecies(transports, mass.first);
			const double inlet_flux = species ? flow.inlet_flow*species->inlet_value : 0.0;
			double outlet_flux = 0.0;
			for (const auto& outlet : result.outlets) {
				const auto flux = outlet.species_flux.find(mass.first);
				if (flux != outlet.species_flux.end()) outlet_flux += flux->second;
			}
			const auto source = result.source_integrals.find(mass.first);
			const double source_integral = source == result.source_integrals.end()
				? 0.0 : source->second;
			result.balance_residuals[mass.first]
				= (mass.second-old->second)/dt_s-inlet_flux+outlet_flux-source_integral;
		}
	}
	return result;
}

class CouplingHistoryWriter {
public:
	explicit CouplingHistoryWriter(std::filesystem::path output_directory,
		SimulationScopeMode mode)
		: output_directory_(std::move(output_directory)), mode_(mode) {}

	void Add(VascularStepResult step, AggregatedVascularReturn venous,
		CircuitAdvanceReport circuit = {})
	{
		steps_.push_back(std::move(step));
		venous_.push_back(std::move(venous));
		circuit_.push_back(std::move(circuit));
	}

	void Write() const
	{
		std::filesystem::create_directories(output_directory_);
		std::ofstream output(output_directory_/"coupling_manifest.json");
		if (!output) throw std::runtime_error("cannot create coupling manifest");
		output << std::setprecision(17)
			<< "{\n  \"schema_version\": 1,\n  \"mode\": \""
			<< SimulationScopeModeName(mode_)
			<< "\",\n  \"scheme\": \"explicit_staggered\",\n"
			<< "  \"lag\": \"vascular step uses arterial state at t; closed-loop reservoir advances from that step's aggregated venous return\",\n"
			<< "  \"units\": {\"time\": \"s\", \"flow\": \"m^3/s\", \"pressure\": \"Pa\", \"concentration\": \"mol/m^3\", \"species_flux\": \"mol/s\", \"species_mass\": \"mol\", \"source_rate\": \"mol/s\"},\n"
			<< "  \"arterial_inlet_history\": [";
		for (std::size_t i = 0; i < steps_.size(); ++i) {
			const auto& step = steps_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< step.inlet.time_s << ", \"flow_m3_s\": ";
			if (step.inlet.has_flow) output << step.inlet.flow_m3_s;
			else output << "null";
			output << ", \"pressure_pa\": ";
			if (step.inlet.has_pressure) output << step.inlet.pressure_pa;
			else output << "null";
			output << ", \"species\": ";
			WriteCouplingNumberMap(output, step.inlet.species);
			output << ", \"temperature_c\": ";
			if (step.inlet.has_temperature) output << step.inlet.temperature_c;
			else output << "null";
			output << ", \"hematocrit_percent\": ";
			if (step.inlet.has_hematocrit) output << step.inlet.hematocrit_percent;
			else output << "null";
			output << '}';
		}
		if (!steps_.empty()) output << '\n';
		output << "  ],\n  \"aggregated_venous_return_history\": [";
		for (std::size_t i = 0; i < venous_.size(); ++i) {
			const auto& venous = venous_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< steps_[i].time_s << ", \"flow_m3_s\": " << venous.flow_m3_s
				<< ", \"pressure_pa\": ";
			if (venous.pressure_valid) output << venous.pressure_pa;
			else output << "null";
			output << ", \"species_flux_mol_s\": ";
			WriteCouplingNumberMap(output, venous.species_flux);
			output << ", \"flux_weighted_concentration\": ";
			WriteCouplingNumberMap(output, venous.flux_weighted_concentration);
			output << ", \"average_valid\": "
				<< (venous.average_valid ? "true" : "false") << ", \"outlet_ids\": ";
			WriteCouplingOutletIds(output, venous.outlet_ids);
			output << '}';
		}
		if (!venous_.empty()) output << '\n';
		output << "  ],\n  \"vca_balance_history\": [";
		for (std::size_t i = 0; i < circuit_.size(); ++i) {
			const auto& circuit = circuit_[i];
			const auto& step = steps_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< circuit.time_s << ", \"volume_change_m3\": "
				<< circuit.volume_change_m3 << ", \"species_mass_change_mol\": ";
			WriteCouplingNumberMap(output, circuit.species_mass_change);
			output << ", \"device_source_mol_s\": ";
			WriteCouplingNumberMap(output, circuit.device_source_rate);
			output << ", \"vascular_total_mass_mol\": ";
			WriteCouplingNumberMap(output, step.total_mass);
			output << ", \"vascular_source_integrals_mol_s\": ";
			WriteCouplingNumberMap(output, step.source_integrals);
			output << ", \"vascular_balance_residuals_mol_s\": ";
			WriteCouplingNumberMap(output, step.balance_residuals);
			output << '}';
		}
		if (!circuit_.empty()) output << '\n';
		output << "  ]\n}\n";
	}

private:
	std::filesystem::path output_directory_;
	SimulationScopeMode mode_;
	std::vector<VascularStepResult> steps_;
	std::vector<AggregatedVascularReturn> venous_;
	std::vector<CircuitAdvanceReport> circuit_;
};

} // namespace iga

#endif

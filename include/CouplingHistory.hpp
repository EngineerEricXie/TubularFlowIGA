#ifndef IGA_COUPLING_HISTORY_HPP
#define IGA_COUPLING_HISTORY_HPP

#include "VascularCoupling.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
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

inline void WriteCouplingMap(std::ostream& output,
	const std::map<std::string, double>& values)
{
	output << '{';
	bool first = true;
	for (const auto& value : values) {
		output << (first ? "" : ", ") << '"' << CouplingEscapeJson(value.first)
			<< "\": " << value.second;
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

class CouplingHistoryWriter {
public:
	explicit CouplingHistoryWriter(std::filesystem::path directory,
		SimulationScopeMode mode)
		: directory_(std::move(directory)), mode_(mode) {}

	void Add(VascularStepResult step, AggregatedVascularReturn venous,
		CircuitAdvanceReport report = {})
	{
		steps_.push_back(std::move(step));
		venous_.push_back(std::move(venous));
		reports_.push_back(std::move(report));
	}

	void Write(const std::string& backend = {}) const
	{
		std::filesystem::create_directories(directory_);
		std::ofstream output(directory_/"coupling_manifest.json");
		if (!output) throw std::runtime_error("cannot create coupling manifest");
		output << std::setprecision(17)
			<< "{\n  \"schema_version\": 1,\n  \"mode\": \""
			<< SimulationScopeModeName(mode_) << '\"';
		if (!backend.empty()) output << ",\n  \"backend\": \"" << CouplingEscapeJson(backend) << '\"';
		output << ",\n  \"scheme\": \"explicit_staggered\",\n"
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
			WriteCouplingMap(output, step.inlet.species);
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
			const auto& value = venous_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< steps_[i].time_s << ", \"flow_m3_s\": " << value.flow_m3_s
				<< ", \"pressure_pa\": ";
			if (value.pressure_valid) output << value.pressure_pa;
			else output << "null";
			output << ", \"species_flux_mol_s\": ";
			WriteCouplingMap(output, value.species_flux);
			output << ", \"flux_weighted_concentration\": ";
			WriteCouplingMap(output, value.flux_weighted_concentration);
			output << ", \"average_valid\": "
				<< (value.average_valid ? "true" : "false") << ", \"outlet_ids\": ";
			WriteCouplingOutletIds(output, value.outlet_ids);
			output << '}';
		}
		if (!venous_.empty()) output << '\n';
		output << "  ],\n  \"vca_balance_history\": [";
		for (std::size_t i = 0; i < reports_.size(); ++i) {
			const auto& step = steps_[i];
			const auto& report = reports_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< report.time_s << ", \"volume_change_m3\": " << report.volume_change_m3
				<< ", \"species_mass_change_mol\": ";
			WriteCouplingMap(output, report.species_mass_change);
			output << ", \"device_source_mol_s\": ";
			WriteCouplingMap(output, report.device_source_rate);
			output << ", \"vascular_total_mass_mol\": ";
			WriteCouplingMap(output, step.total_mass);
			output << ", \"vascular_source_integrals_mol_s\": ";
			WriteCouplingMap(output, step.source_integrals);
			output << ", \"vascular_balance_residuals_mol_s\": ";
			WriteCouplingMap(output, step.balance_residuals);
			output << '}';
		}
		if (!reports_.empty()) output << '\n';
		output << "  ]\n}\n";
	}

private:
	std::filesystem::path directory_;
	SimulationScopeMode mode_;
	std::vector<VascularStepResult> steps_;
	std::vector<AggregatedVascularReturn> venous_;
	std::vector<CircuitAdvanceReport> reports_;
};

} // namespace iga

#endif

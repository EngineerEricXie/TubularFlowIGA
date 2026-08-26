#ifndef IGA_COUPLING_HISTORY_HPP
#define IGA_COUPLING_HISTORY_HPP

#include "VascularCoupling.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

inline void WriteCouplingMap(std::ostream& output,
	const std::map<std::string, double>& values)
{
	output << '{';
	bool first = true;
	for (const auto& value : values) {
		output << (first ? "" : ", ") << '"' << value.first << "\": " << value.second;
		first = false;
	}
	output << '}';
}

class CouplingHistoryWriter {
public:
	explicit CouplingHistoryWriter(std::filesystem::path directory,
		SimulationScopeMode mode)
		: directory_(std::move(directory)), mode_(mode) {}

	void Add(VascularStepResult step, AggregatedVascularReturn venous,
		CircuitAdvanceReport report)
	{
		steps_.push_back(std::move(step));
		venous_.push_back(std::move(venous));
		reports_.push_back(std::move(report));
	}

	void Write(const std::string& backend) const
	{
		std::filesystem::create_directories(directory_);
		std::ofstream output(directory_/"coupling_manifest.json");
		if (!output) throw std::runtime_error("cannot create coupling manifest");
		output << std::setprecision(17)
			<< "{\n  \"schema_version\": 1,\n  \"mode\": \""
			<< SimulationScopeModeName(mode_) << "\",\n  \"backend\": \""
			<< backend << "\",\n  \"arterial_inlet_history\": [";
		for (std::size_t i = 0; i < steps_.size(); ++i) {
			const auto& inlet = steps_[i].inlet;
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": " << inlet.time_s
				<< ", \"flow_m3_s\": " << inlet.flow_m3_s << ", \"species\": ";
			WriteCouplingMap(output, inlet.species);
			output << '}';
		}
		if (!steps_.empty()) output << '\n';
		output << "  ],\n  \"aggregated_venous_return_history\": [";
		for (std::size_t i = 0; i < venous_.size(); ++i) {
			const auto& value = venous_[i];
			output << (i == 0 ? "\n" : ",\n") << "    {\"time_s\": "
				<< steps_[i].time_s << ", \"flow_m3_s\": " << value.flow_m3_s
				<< ", \"pressure_pa\": " << value.pressure_pa
				<< ", \"species_flux_mol_s\": ";
			WriteCouplingMap(output, value.species_flux);
			output << ", \"outlet_ids\": [";
			for (std::size_t j = 0; j < value.outlet_ids.size(); ++j)
				output << (j == 0 ? "" : ", ") << value.outlet_ids[j];
			output << "]}";
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

#ifndef IGA_PHYSIOLOGY_OUTPUT_HPP
#define IGA_PHYSIOLOGY_OUTPUT_HPP

#include "SimulationConfig.hpp"
#include "VtkOutput.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

inline std::string EscapePhysiologyJson(const std::string& value)
{
	std::string result;
	for (const char character : value) {
		if (character == '\\' || character == '"') result.push_back('\\');
		result.push_back(character);
	}
	return result;
}

inline double HillSaturation(double pressure_mmhg, double p50, double exponent)
{
	const double pressure = std::max(pressure_mmhg, 0.0);
	const double numerator = std::pow(pressure, exponent);
	return numerator/(numerator+std::pow(p50, exponent));
}

inline std::vector<VtkPointArray> ComputePhysiologyPointArrays(
	const PhysiologyDefinition& physiology, const std::vector<std::string>& fields,
	const std::vector<double>& interleaved_values)
{
	std::vector<VtkPointArray> result;
	if (!physiology.enabled || fields.empty()) return result;
	if (interleaved_values.size()%fields.size() != 0)
		throw std::runtime_error("physiology output has inconsistent field storage");
	const auto nodes = interleaved_values.size()/fields.size();
	std::map<std::string, std::size_t> index;
	for (std::size_t i = 0; i < fields.size(); ++i) index.emplace(fields[i], i);
	auto values = [&](const std::string& field) {
		std::vector<double> output(nodes);
		const auto position = index.at(field);
		for (std::size_t node = 0; node < nodes; ++node)
			output[node] = interleaved_values[node*fields.size()+position];
		return output;
	};
	for (const auto& name : physiology.derived_fields) {
		VtkPointArray array{name, 1, std::vector<double>(nodes, 0.0)};
		if (name == "hematocrit") {
			std::fill(array.values.begin(), array.values.end(), physiology.hematocrit_percent);
		} else if (name == "pO2" || name == "SaO2" || name == "SvO2"
			|| name == "dissolved_oxygen" || name == "bound_oxygen"
			|| name == "total_oxygen") {
			const auto oxygen = values("oxygen");
			for (std::size_t node = 0; node < nodes; ++node) {
				const double po2 = oxygen[node]/std::max(physiology.oxygen_solubility, 1.0e-30);
				const double saturation = HillSaturation(po2,
					physiology.p50_mmhg, physiology.hill_exponent);
				if (name == "pO2") array.values[node] = po2;
				else if (name == "SaO2" || name == "SvO2") array.values[node] = saturation;
				else if (name == "dissolved_oxygen") array.values[node] = oxygen[node];
				else {
					const double bound = 1.34*physiology.hemoglobin_g_dl*saturation;
					array.values[node] = name == "bound_oxygen" ? bound : oxygen[node]+bound;
				}
			}
		} else if (name == "pCO2") {
			const auto carbon_dioxide = values("carbon_dioxide");
			for (std::size_t node = 0; node < nodes; ++node)
				array.values[node] = carbon_dioxide[node]/0.0301;
		} else if (name == "pH") {
			const auto carbon_dioxide = values("carbon_dioxide");
			const auto bicarbonate = values("bicarbonate");
			for (std::size_t node = 0; node < nodes; ++node) {
				const double pco2 = std::max(carbon_dioxide[node]/0.0301, 1.0e-30);
				array.values[node] = 6.1+std::log10(std::max(bicarbonate[node], 1.0e-30)/(0.0301*pco2));
			}
		}
		result.push_back(std::move(array));
	}
	return result;
}

inline void WritePhysiologyManifest(const std::filesystem::path& path,
	const PhysiologyDefinition& physiology, const std::vector<std::string>& fields)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create physiology manifest: "+path.string());
	output << "{\n  \"enabled\": " << (physiology.enabled ? "true" : "false")
		<< ",\n  \"fields\": {";
	bool first = true;
	for (const auto& field : fields) {
		output << (first ? "\n" : ",\n") << "    \"" << EscapePhysiologyJson(field)
			<< "\": {\"status\": \"solved\"}";
		first = false;
	}
	for (const auto& field : physiology.derived_fields) {
		output << (first ? "\n" : ",\n") << "    \"" << EscapePhysiologyJson(field)
			<< "\": {\"status\": \"derived\"}";
		first = false;
	}
	if (!first) output << '\n';
	output << "  }\n}\n";
}

} // namespace iga

#endif

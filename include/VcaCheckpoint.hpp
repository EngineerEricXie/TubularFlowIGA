#ifndef IGA_VCA_CHECKPOINT_HPP
#define IGA_VCA_CHECKPOINT_HPP

#include "CaseConfig.hpp"
#include "VascularCoupling.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct VcaCheckpointMetadata {
	int schema_version = 1;
	int completed_step = 0;
	double physical_time = 0.0;
	double dt = 0.0;
	std::vector<std::string> fields;
	std::string transport_state_file;
	ReservoirState reservoir;
};

inline std::filesystem::path VcaCheckpointMetadataPath(const std::filesystem::path& prefix)
{
	return prefix.string()+".vca.json";
}

inline std::filesystem::path VcaCheckpointTransportStatePath(const std::filesystem::path& prefix)
{
	return prefix.string()+".vca_transport.state";
}

inline std::string SerializeVcaCheckpointMetadata(const VcaCheckpointMetadata& metadata)
{
	std::ostringstream output;
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": " << metadata.schema_version
		<< ",\n  \"completed_step\": " << metadata.completed_step
		<< ",\n  \"physical_time\": " << metadata.physical_time
		<< ",\n  \"dt\": " << metadata.dt
		<< ",\n  \"fields\": [";
	for (std::size_t i = 0; i < metadata.fields.size(); ++i)
		output << (i == 0 ? "" : ", ") << '"' << metadata.fields[i] << '"';
	output << "],\n  \"transport_state_file\": \"" << metadata.transport_state_file
		<< "\",\n  \"reservoir\": {\"volume_m3\": " << metadata.reservoir.volume_m3
		<< ", \"temperature_c\": " << metadata.reservoir.temperature_c
		<< ", \"hematocrit_percent\": " << metadata.reservoir.hematocrit_percent
		<< ", \"species\": {";
	bool first = true;
	for (const auto& species : metadata.reservoir.species) {
		output << (first ? "" : ", ") << '"' << species.first << "\": " << species.second;
		first = false;
	}
	output << "}}\n}\n";
	return output.str();
}

inline VcaCheckpointMetadata ParseVcaCheckpointMetadata(const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "VCA checkpoint");
	RequireKnownKeys(root, {"schema_version", "completed_step", "physical_time", "dt",
		"fields", "transport_state_file", "reservoir"}, "VCA checkpoint");
	auto required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(root, key);
		if (!value) throw std::runtime_error("VCA checkpoint requires '"+key+"'");
		return *value;
	};
	VcaCheckpointMetadata result;
	result.schema_version = RequireInteger(required("schema_version"), "VCA checkpoint.schema_version");
	result.completed_step = RequireInteger(required("completed_step"), "VCA checkpoint.completed_step");
	result.physical_time = RequireNumber(required("physical_time"), "VCA checkpoint.physical_time");
	result.dt = RequireNumber(required("dt"), "VCA checkpoint.dt");
	const auto& fields = RequireArray(required("fields"), "VCA checkpoint.fields");
	for (std::size_t i = 0; i < fields.size(); ++i)
		result.fields.push_back(RequireString(fields[i], "VCA checkpoint.fields["+std::to_string(i)+"]"));
	result.transport_state_file = RequireString(required("transport_state_file"),
		"VCA checkpoint.transport_state_file");
	const auto& reservoir = RequireObject(required("reservoir"), "VCA checkpoint.reservoir");
	RequireKnownKeys(reservoir, {"volume_m3", "temperature_c", "hematocrit_percent", "species"},
		"VCA checkpoint.reservoir");
	auto reservoir_required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(reservoir, key);
		if (!value) throw std::runtime_error("VCA checkpoint.reservoir requires '"+key+"'");
		return *value;
	};
	result.reservoir.volume_m3 = RequireNumber(reservoir_required("volume_m3"), "VCA checkpoint.reservoir.volume_m3");
	result.reservoir.temperature_c = RequireNumber(reservoir_required("temperature_c"), "VCA checkpoint.reservoir.temperature_c");
	result.reservoir.hematocrit_percent = RequireNumber(reservoir_required("hematocrit_percent"), "VCA checkpoint.reservoir.hematocrit_percent");
	const auto& species = RequireObject(reservoir_required("species"), "VCA checkpoint.reservoir.species");
	for (const auto& item : species)
		result.reservoir.species.emplace(item.first, RequireNumber(item.second,
			"VCA checkpoint.reservoir.species."+item.first));
	if (result.schema_version != 1 || result.completed_step < 0 || result.physical_time < 0.0
		|| !(result.dt > 0.0) || result.fields.empty() || result.transport_state_file.empty()
		|| !(result.reservoir.volume_m3 > 0.0) || !std::isfinite(result.reservoir.temperature_c)
		|| result.reservoir.hematocrit_percent < 0.0 || result.reservoir.hematocrit_percent > 100.0)
		throw std::runtime_error("VCA checkpoint metadata contains invalid values");
	return result;
}

inline void WriteVcaCheckpointMetadata(const std::filesystem::path& prefix,
	const VcaCheckpointMetadata& metadata)
{
	std::ofstream output(VcaCheckpointMetadataPath(prefix));
	if (!output) throw std::runtime_error("cannot create VCA checkpoint metadata");
	output << SerializeVcaCheckpointMetadata(metadata);
	if (!output) throw std::runtime_error("cannot write VCA checkpoint metadata");
}

inline VcaCheckpointMetadata ReadVcaCheckpointMetadata(const std::filesystem::path& prefix)
{
	std::ifstream input(VcaCheckpointMetadataPath(prefix));
	if (!input) throw std::runtime_error("cannot open VCA checkpoint metadata");
	std::ostringstream text;
	text << input.rdbuf();
	return ParseVcaCheckpointMetadata(text.str());
}

inline void ValidateVcaCheckpoint(const VcaCheckpointMetadata& metadata, int completed_step,
	double dt, const std::vector<std::string>& fields)
{
	auto close = [](double a, double b) { return std::abs(a-b) <= 1e-12*std::max({1.0, std::abs(a), std::abs(b)}); };
	if (metadata.completed_step != completed_step || !close(metadata.physical_time, completed_step*dt)
		|| !close(metadata.dt, dt) || metadata.fields != fields)
		throw std::runtime_error("VCA checkpoint does not match flow/transport configuration");
}

} // namespace iga

#endif

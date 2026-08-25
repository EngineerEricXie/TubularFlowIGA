#ifndef IGA_TRANSPORT_CHECKPOINT_HPP
#define IGA_TRANSPORT_CHECKPOINT_HPP

#include "CaseConfig.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct TransportCheckpointMetadata {
	int schema_version = 1;
	std::uint64_t nodes = 0;
	std::vector<std::string> fields;
	std::string system;
	std::string velocity_source;
	int completed_step = 0;
	double physical_time = 0.0;
	double dt = 0.0;
	std::string state_file;
	std::string state_format = "petsc_binary";
};

inline std::filesystem::path TransportCheckpointMetadataPath(
	const std::filesystem::path& prefix)
{
	return prefix.string()+".json";
}

inline std::filesystem::path TransportCheckpointStatePath(
	const std::filesystem::path& prefix)
{
	return prefix.string()+".state";
}

inline std::string EscapeTransportCheckpointJson(const std::string& value)
{
	std::string result;
	for (const char character : value) {
		if (character == '\\' || character == '"') result.push_back('\\');
		if (character == '\n') result += "\\n";
		else if (character == '\r') result += "\\r";
		else if (character == '\t') result += "\\t";
		else result.push_back(character);
	}
	return result;
}

inline std::string SerializeTransportCheckpointMetadata(
	const TransportCheckpointMetadata& metadata)
{
	std::ostringstream output;
	output << std::setprecision(17)
		<< "{\n"
		<< "  \"schema_version\": " << metadata.schema_version << ",\n"
		<< "  \"nodes\": " << metadata.nodes << ",\n"
		<< "  \"fields\": [";
	for (std::size_t i = 0; i < metadata.fields.size(); ++i)
		output << (i == 0 ? "" : ", ") << '"'
			<< EscapeTransportCheckpointJson(metadata.fields[i]) << '"';
	output << "],\n"
		<< "  \"system\": \"" << EscapeTransportCheckpointJson(metadata.system) << "\",\n"
		<< "  \"velocity_source\": \""
		<< EscapeTransportCheckpointJson(metadata.velocity_source) << "\",\n"
		<< "  \"completed_step\": " << metadata.completed_step << ",\n"
		<< "  \"physical_time\": " << metadata.physical_time << ",\n"
		<< "  \"dt\": " << metadata.dt << ",\n"
		<< "  \"state_file\": \"" << EscapeTransportCheckpointJson(metadata.state_file) << "\",\n"
		<< "  \"state_format\": \"" << EscapeTransportCheckpointJson(metadata.state_format) << "\"\n"
		<< "}\n";
	return output.str();
}

inline TransportCheckpointMetadata ParseTransportCheckpointMetadata(
	const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "transport checkpoint");
	RequireKnownKeys(root, {"schema_version", "nodes", "fields", "system",
		"velocity_source", "completed_step", "physical_time", "dt",
		"state_file", "state_format"}, "transport checkpoint");
	auto required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(root, key);
		if (!value) throw std::runtime_error("transport checkpoint requires '"+key+"'");
		return *value;
	};
	TransportCheckpointMetadata metadata;
	metadata.schema_version = RequireInteger(
		required("schema_version"), "transport checkpoint.schema_version");
	const auto nodes = RequireNumber(required("nodes"), "transport checkpoint.nodes");
	if (!(nodes > 0.0) || std::floor(nodes) != nodes
		|| nodes > static_cast<double>(std::numeric_limits<std::uint64_t>::max()))
		throw std::runtime_error("transport checkpoint.nodes must be a positive integer");
	metadata.nodes = static_cast<std::uint64_t>(nodes);
	const auto& fields = RequireArray(required("fields"), "transport checkpoint.fields");
	std::set<std::string> unique_fields;
	for (std::size_t i = 0; i < fields.size(); ++i) {
		const auto field = RequireString(
			fields[i], "transport checkpoint.fields["+std::to_string(i)+"]");
		if (field.empty() || !unique_fields.insert(field).second)
			throw std::runtime_error("transport checkpoint fields must be non-empty and unique");
		metadata.fields.push_back(field);
	}
	metadata.system = RequireString(required("system"), "transport checkpoint.system");
	metadata.velocity_source = RequireString(
		required("velocity_source"), "transport checkpoint.velocity_source");
	metadata.completed_step = RequireInteger(
		required("completed_step"), "transport checkpoint.completed_step");
	metadata.physical_time = RequireNumber(
		required("physical_time"), "transport checkpoint.physical_time");
	metadata.dt = RequireNumber(required("dt"), "transport checkpoint.dt");
	metadata.state_file = RequireString(
		required("state_file"), "transport checkpoint.state_file");
	metadata.state_format = RequireString(
		required("state_format"), "transport checkpoint.state_format");
	if (metadata.schema_version != 1 || metadata.fields.empty()
		|| metadata.system.empty() || metadata.velocity_source.empty()
		|| metadata.completed_step < 0 || metadata.physical_time < 0.0
		|| !(metadata.dt > 0.0) || metadata.state_file.empty()
		|| (metadata.state_format != "petsc_binary"
			&& metadata.state_format != "raw_float64"))
		throw std::runtime_error("transport checkpoint metadata contains invalid values");
	return metadata;
}

inline void WriteTransportCheckpointMetadata(const std::filesystem::path& prefix,
	const TransportCheckpointMetadata& metadata)
{
	const auto path = TransportCheckpointMetadataPath(prefix);
	std::ofstream output(path);
	if (!output)
		throw std::runtime_error("cannot create transport checkpoint metadata: "+path.string());
	output << SerializeTransportCheckpointMetadata(metadata);
	if (!output)
		throw std::runtime_error("cannot write transport checkpoint metadata: "+path.string());
}

inline TransportCheckpointMetadata ReadTransportCheckpointMetadata(
	const std::filesystem::path& prefix)
{
	const auto path = TransportCheckpointMetadataPath(prefix);
	std::ifstream input(path);
	if (!input)
		throw std::runtime_error("cannot open transport checkpoint metadata: "+path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof())
		throw std::runtime_error("cannot read transport checkpoint metadata: "+path.string());
	return ParseTransportCheckpointMetadata(contents.str());
}

inline void ValidateTransportCheckpoint(const TransportCheckpointMetadata& metadata,
	std::uint64_t nodes, const std::vector<std::string>& fields,
	const std::string& system, const std::string& velocity_source,
	int configured_steps, double dt)
{
	auto close = [](double first, double second) {
		return std::abs(first-second)
			<= 1e-12*std::max({1.0, std::abs(first), std::abs(second)});
	};
	if (metadata.nodes != nodes)
		throw std::runtime_error("transport checkpoint node count does not match database");
	if (metadata.fields != fields)
		throw std::runtime_error("transport checkpoint fields do not match equation system");
	if (metadata.system != system)
		throw std::runtime_error("transport checkpoint system does not match configuration");
	if (metadata.velocity_source != velocity_source)
		throw std::runtime_error("transport checkpoint velocity source does not match configuration");
	if (metadata.completed_step > configured_steps)
		throw std::runtime_error("transport checkpoint is beyond configured final step");
	if (!close(metadata.physical_time, metadata.completed_step*dt))
		throw std::runtime_error(
			"transport checkpoint physical time is inconsistent with completed step");
	if (!close(metadata.dt, dt))
		throw std::runtime_error("transport checkpoint dt does not match configuration");
}

} // namespace iga

#endif

#ifndef IGA_FLOW_CHECKPOINT_HPP
#define IGA_FLOW_CHECKPOINT_HPP

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

struct FlowCheckpointOutlet {
	int label = -1;
	std::string kind;
	double flow = 0.0;
	double pressure = 0.0;
	double capacitor_pressure = 0.0;
};

struct FlowCheckpointMetadata {
	int schema_version = 1;
	std::uint64_t nodes = 0;
	int fields = 4;
	int completed_step = 0;
	double physical_time = 0.0;
	double dt = 0.0;
	double density = 0.0;
	double viscosity = 0.0;
	std::string state_file;
	std::string state_format = "petsc_binary";
	std::vector<FlowCheckpointOutlet> outlets;
};

inline std::filesystem::path FlowCheckpointMetadataPath(const std::filesystem::path& prefix)
{
	return prefix.string() + ".json";
}

inline std::filesystem::path FlowCheckpointStatePath(const std::filesystem::path& prefix)
{
	return prefix.string() + ".state";
}

inline std::filesystem::path TimeIndexedPath(const std::filesystem::path& base, int step)
{
	if (step < 0) throw std::runtime_error("time-indexed output step cannot be negative");
	std::ostringstream suffix;
	suffix << ".step" << std::setw(6) << std::setfill('0') << step;
	const auto extension = base.extension().string();
	if (extension.empty()) return base.string() + suffix.str();
	return base.parent_path() / (base.stem().string() + suffix.str() + extension);
}

inline std::string SerializeFlowCheckpointMetadata(const FlowCheckpointMetadata& metadata)
{
	auto escape = [](const std::string& value) {
		std::string result;
		for (const char character : value) {
			if (character == '\\' || character == '"') result.push_back('\\');
			if (character == '\n') result += "\\n";
			else if (character == '\r') result += "\\r";
			else if (character == '\t') result += "\\t";
			else result.push_back(character);
		}
		return result;
	};
	std::ostringstream output;
	output << std::setprecision(17)
		<< "{\n"
		<< "  \"schema_version\": " << metadata.schema_version << ",\n"
		<< "  \"nodes\": " << metadata.nodes << ",\n"
		<< "  \"fields\": " << metadata.fields << ",\n"
		<< "  \"completed_step\": " << metadata.completed_step << ",\n"
		<< "  \"physical_time\": " << metadata.physical_time << ",\n"
		<< "  \"dt\": " << metadata.dt << ",\n"
		<< "  \"density\": " << metadata.density << ",\n"
		<< "  \"viscosity\": " << metadata.viscosity << ",\n"
		<< "  \"state_file\": \"" << escape(metadata.state_file) << "\",\n"
		<< "  \"state_format\": \"" << escape(metadata.state_format) << "\",\n"
		<< "  \"outlets\": [";
	for (std::size_t i = 0; i < metadata.outlets.size(); ++i) {
		const auto& outlet = metadata.outlets[i];
		output << (i == 0 ? "\n" : ",\n")
			<< "    {\"label\": " << outlet.label
			<< ", \"kind\": \"" << escape(outlet.kind)
			<< "\", \"flow\": " << outlet.flow
			<< ", \"pressure\": " << outlet.pressure
			<< ", \"capacitor_pressure\": " << outlet.capacitor_pressure << "}";
	}
	if (!metadata.outlets.empty()) output << '\n';
	output << "  ]\n"
		<< "}\n";
	return output.str();
}

inline FlowCheckpointMetadata ParseFlowCheckpointMetadata(const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "flow checkpoint");
	RequireKnownKeys(root, {"schema_version", "nodes", "fields", "completed_step",
		"physical_time", "dt", "density", "viscosity", "state_file", "state_format",
		"outlets"}, "flow checkpoint");
	auto required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(root, key);
		if (!value) throw std::runtime_error("flow checkpoint requires '" + key + "'");
		return *value;
	};
	FlowCheckpointMetadata metadata;
	metadata.schema_version = RequireInteger(required("schema_version"), "flow checkpoint.schema_version");
	const auto nodes = RequireNumber(required("nodes"), "flow checkpoint.nodes");
	if (!(nodes > 0.0) || std::floor(nodes) != nodes
		|| nodes > static_cast<double>(std::numeric_limits<std::uint64_t>::max()))
		throw std::runtime_error("flow checkpoint.nodes must be a positive integer");
	metadata.nodes = static_cast<std::uint64_t>(nodes);
	metadata.fields = RequireInteger(required("fields"), "flow checkpoint.fields");
	metadata.completed_step = RequireInteger(required("completed_step"), "flow checkpoint.completed_step");
	metadata.physical_time = RequireNumber(required("physical_time"), "flow checkpoint.physical_time");
	metadata.dt = RequireNumber(required("dt"), "flow checkpoint.dt");
	metadata.density = RequireNumber(required("density"), "flow checkpoint.density");
	metadata.viscosity = RequireNumber(required("viscosity"), "flow checkpoint.viscosity");
	metadata.state_file = RequireString(required("state_file"), "flow checkpoint.state_file");
	if (const auto* format = Find(root, "state_format"))
		metadata.state_format = RequireString(*format, "flow checkpoint.state_format");
	if (const auto* outlets = Find(root, "outlets")) {
		const auto& array = RequireArray(*outlets, "flow checkpoint.outlets");
		std::set<int> labels;
		for (std::size_t i = 0; i < array.size(); ++i) {
			const auto context = "flow checkpoint.outlets[" + std::to_string(i) + "]";
			const auto& object = RequireObject(array[i], context);
			RequireKnownKeys(object, {"label", "kind", "flow", "pressure",
				"capacitor_pressure"}, context);
			auto item = [&](const std::string& key) -> const JsonValue& {
				const auto* value = Find(object, key);
				if (!value) throw std::runtime_error(context + " requires '" + key + "'");
				return *value;
			};
			FlowCheckpointOutlet outlet;
			outlet.label = RequireInteger(item("label"), context + ".label");
			outlet.kind = RequireString(item("kind"), context + ".kind");
			outlet.flow = RequireNumber(item("flow"), context + ".flow");
			outlet.pressure = RequireNumber(item("pressure"), context + ".pressure");
			outlet.capacitor_pressure = RequireNumber(
				item("capacitor_pressure"), context + ".capacitor_pressure");
			if (!labels.insert(outlet.label).second
				|| (outlet.kind != "resistance" && outlet.kind != "windkessel_rc"
					&& outlet.kind != "windkessel_rcr")
				|| !std::isfinite(outlet.flow) || !std::isfinite(outlet.pressure)
				|| !std::isfinite(outlet.capacitor_pressure))
				throw std::runtime_error(context + " contains invalid outlet state");
			metadata.outlets.push_back(std::move(outlet));
		}
	}
	if (metadata.schema_version != 1) throw std::runtime_error("unsupported flow checkpoint schema_version");
	if (metadata.nodes == 0 || metadata.fields != 4 || metadata.completed_step < 0
		|| metadata.physical_time < 0.0 || !(metadata.dt > 0.0)
		|| !(metadata.density > 0.0) || !(metadata.viscosity > 0.0)
		|| metadata.state_file.empty()
		|| (metadata.state_format != "petsc_binary" && metadata.state_format != "raw_float64"))
		throw std::runtime_error("flow checkpoint metadata contains invalid values");
	return metadata;
}

inline void WriteFlowCheckpointMetadata(const std::filesystem::path& prefix,
	const FlowCheckpointMetadata& metadata)
{
	const auto path = FlowCheckpointMetadataPath(prefix);
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create flow checkpoint metadata: " + path.string());
	output << SerializeFlowCheckpointMetadata(metadata);
	if (!output) throw std::runtime_error("cannot write flow checkpoint metadata: " + path.string());
}

inline FlowCheckpointMetadata ReadFlowCheckpointMetadata(const std::filesystem::path& prefix)
{
	const auto path = FlowCheckpointMetadataPath(prefix);
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open flow checkpoint metadata: " + path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof())
		throw std::runtime_error("cannot read flow checkpoint metadata: " + path.string());
	return ParseFlowCheckpointMetadata(contents.str());
}

inline void ValidateFlowCheckpoint(const FlowCheckpointMetadata& metadata,
	std::uint64_t nodes, int configured_steps, double dt, double density, double viscosity)
{
	auto close = [](double first, double second) {
		return std::abs(first-second) <= 1e-12*std::max({1.0, std::abs(first), std::abs(second)});
	};
	if (metadata.nodes != nodes) throw std::runtime_error("flow checkpoint node count does not match database");
	if (metadata.completed_step > configured_steps)
		throw std::runtime_error("flow checkpoint is beyond configured final step");
	if (!close(metadata.physical_time, metadata.completed_step*dt))
		throw std::runtime_error("flow checkpoint physical time is inconsistent with completed step");
	if (!close(metadata.dt, dt)) throw std::runtime_error("flow checkpoint dt does not match configuration");
	if (!close(metadata.density, density)) throw std::runtime_error("flow checkpoint density does not match configuration");
	if (!close(metadata.viscosity, viscosity)) throw std::runtime_error("flow checkpoint viscosity does not match configuration");
}

} // namespace iga

#endif

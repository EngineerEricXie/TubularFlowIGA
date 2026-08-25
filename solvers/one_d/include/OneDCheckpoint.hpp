#ifndef IGA_ONE_D_CHECKPOINT_HPP
#define IGA_ONE_D_CHECKPOINT_HPP

#include "OneDImplicit.hpp"
#include "OneDOutput.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OneDCheckpointMetadata {
	int schema_version = 2;
	int completed_step = 0;
	int internal_substeps = 0;
	double physical_time = 0.0;
	double dt = 0.0;
	double inlet_flow = 0.0;
	int cells = 0;
	int nodes = 0;
	int segments = 0;
	int outlets = 0;
	std::vector<std::string> species;
	std::uint64_t config_fingerprint = 0;
	std::uint64_t network_fingerprint = 0;
	std::string state_file;
};

inline std::uint64_t OneDFingerprint(const std::string& text)
{
	std::uint64_t hash = 1469598103934665603ull;
	for (const unsigned char character : text) {
		hash ^= character;
		hash *= 1099511628211ull;
	}
	return hash;
}

inline std::uint64_t OneDNetworkFingerprint(const OneDNetwork& network)
{
	std::ostringstream text;
	text << std::setprecision(17);
	for (const auto& node : network.nodes)
		text << node.id << ' ' << node.parent_id << ' ' << node.position[0] << ' '
			<< node.position[1] << ' ' << node.position[2] << ' ' << node.radius << '\n';
	return OneDFingerprint(text.str());
}

inline std::vector<std::string> OneDCheckpointSpecies(
	const std::vector<OneDTransportState>& transports)
{
	std::vector<std::string> names;
	for (const auto& transport : transports)
		for (const auto& species : transport.species) names.push_back(species.definition.field);
	return names;
}

inline std::vector<double> PackOneDCheckpointState(const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports, const OneDNetwork& network)
{
	std::vector<double> values;
	values.insert(values.end(), flow.area.begin(), flow.area.end());
	values.insert(values.end(), flow.flow.begin(), flow.flow.end());
	values.insert(values.end(), flow.pressure.begin(), flow.pressure.end());
	values.insert(values.end(), flow.node_pressure.begin(), flow.node_pressure.end());
	values.insert(values.end(), flow.segment_flow.begin(), flow.segment_flow.end());
	for (const auto& outlet : flow.outlets) {
		values.push_back(outlet.pressure);
		values.push_back(outlet.capacitor_pressure);
		values.push_back(outlet.flow);
	}
	for (const auto& segment : network.segments) values.push_back(segment.radius0);
	for (const auto& transport : transports)
		for (const auto& species : transport.species)
			values.insert(values.end(), species.concentration.begin(), species.concentration.end());
	return values;
}

inline void UnpackOneDCheckpointState(const std::vector<double>& values,
	OneDFlowState& flow, std::vector<OneDTransportState>& transports,
	const OneDCheckpointMetadata& metadata, OneDNetwork& network,
	double dynamic_viscosity)
{
	const std::size_t expected = static_cast<std::size_t>(3*metadata.cells+metadata.nodes
		+2*metadata.segments+3*metadata.outlets+metadata.cells*metadata.species.size());
	if (values.size() != expected) throw std::runtime_error("1d checkpoint state size is invalid");
	std::size_t offset = 0;
	auto assign = [&](std::vector<double>& target, std::size_t count) {
		target.assign(values.begin()+static_cast<std::ptrdiff_t>(offset),
			values.begin()+static_cast<std::ptrdiff_t>(offset+count));
		offset += count;
	};
	assign(flow.area, metadata.cells);
	assign(flow.flow, metadata.cells);
	assign(flow.pressure, metadata.cells);
	assign(flow.node_pressure, metadata.nodes);
	assign(flow.segment_flow, metadata.segments);
	if (flow.outlets.size() != static_cast<std::size_t>(metadata.outlets))
		throw std::runtime_error("1d checkpoint outlet count does not match configuration");
	for (auto& outlet : flow.outlets) {
		outlet.pressure = values[offset++];
		outlet.capacitor_pressure = values[offset++];
		outlet.flow = values[offset++];
	}
	if (network.segments.size() != static_cast<std::size_t>(metadata.segments))
		throw std::runtime_error("1d checkpoint segment count does not match configuration");
	for (auto& segment : network.segments) {
		segment.radius0 = values[offset++];
		if (!(segment.radius0 > 0.0) || !std::isfinite(segment.radius0))
			throw std::runtime_error("1d checkpoint contains an invalid dynamic radius");
		segment.area0 = OneDPi*segment.radius0*segment.radius0;
		segment.resistance = 8.0*dynamic_viscosity*segment.length
			/(OneDPi*std::pow(segment.radius0, 4.0));
	}
	std::size_t species_index = 0;
	for (auto& transport : transports)
		for (auto& species : transport.species) {
			if (species_index >= metadata.species.size()
				|| species.definition.field != metadata.species[species_index])
				throw std::runtime_error("1d checkpoint species ordering does not match configuration");
			assign(species.concentration, metadata.cells);
			++species_index;
		}
}

inline std::filesystem::path OneDCheckpointMetadataPath(const std::filesystem::path& prefix)
{
	return prefix.string()+".json";
}

inline std::filesystem::path OneDCheckpointStatePath(const std::filesystem::path& prefix)
{
	return prefix.string()+".state";
}

inline std::string SerializeOneDCheckpointMetadata(const OneDCheckpointMetadata& metadata)
{
	std::ostringstream output;
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": " << metadata.schema_version << ",\n"
		<< "  \"completed_step\": " << metadata.completed_step << ",\n"
		<< "  \"internal_substeps\": " << metadata.internal_substeps << ",\n"
		<< "  \"physical_time\": " << metadata.physical_time << ",\n"
		<< "  \"dt\": " << metadata.dt << ",\n"
		<< "  \"inlet_flow\": " << metadata.inlet_flow << ",\n"
		<< "  \"cells\": " << metadata.cells << ",\n"
		<< "  \"nodes\": " << metadata.nodes << ",\n"
		<< "  \"segments\": " << metadata.segments << ",\n"
		<< "  \"outlets\": " << metadata.outlets << ",\n"
		<< "  \"config_fingerprint\": \"" << metadata.config_fingerprint << "\",\n"
		<< "  \"network_fingerprint\": \"" << metadata.network_fingerprint << "\",\n"
		<< "  \"state_file\": \"" << OneDEscapeJson(metadata.state_file) << "\",\n"
		<< "  \"species\": [";
	for (std::size_t i = 0; i < metadata.species.size(); ++i)
		output << (i ? ", " : "") << '"' << OneDEscapeJson(metadata.species[i]) << '"';
	output << "]\n}\n";
	return output.str();
}

inline OneDCheckpointMetadata ParseOneDCheckpointMetadata(const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "1d checkpoint");
	RequireKnownKeys(root, {"schema_version", "completed_step", "internal_substeps", "physical_time", "dt", "inlet_flow",
		"cells", "nodes", "segments", "outlets", "config_fingerprint",
		"network_fingerprint", "state_file", "species"}, "1d checkpoint");
	auto required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(root, key);
		if (!value) throw std::runtime_error("1d checkpoint requires '"+key+"'");
		return *value;
	};
	OneDCheckpointMetadata result;
	result.schema_version = RequireInteger(required("schema_version"), "1d checkpoint.schema_version");
	result.completed_step = RequireInteger(required("completed_step"), "1d checkpoint.completed_step");
	result.internal_substeps = RequireInteger(required("internal_substeps"), "1d checkpoint.internal_substeps");
	result.physical_time = RequireNumber(required("physical_time"), "1d checkpoint.physical_time");
	result.dt = RequireNumber(required("dt"), "1d checkpoint.dt");
	result.inlet_flow = RequireNumber(required("inlet_flow"), "1d checkpoint.inlet_flow");
	result.cells = RequireInteger(required("cells"), "1d checkpoint.cells");
	result.nodes = RequireInteger(required("nodes"), "1d checkpoint.nodes");
	result.segments = RequireInteger(required("segments"), "1d checkpoint.segments");
	result.outlets = RequireInteger(required("outlets"), "1d checkpoint.outlets");
	auto integer64 = [&](const char* key) {
		const auto value = RequireString(required(key), std::string("1d checkpoint.")+key);
		std::size_t used = 0;
		std::uint64_t parsed = 0;
		try { parsed = std::stoull(value, &used); }
		catch (const std::exception&) { used = 0; }
		if (used != value.size())
			throw std::runtime_error(std::string("1d checkpoint.")+key+" must be an unsigned integer string");
		return parsed;
	};
	result.config_fingerprint = integer64("config_fingerprint");
	result.network_fingerprint = integer64("network_fingerprint");
	result.state_file = RequireString(required("state_file"), "1d checkpoint.state_file");
	const auto& species = RequireArray(required("species"), "1d checkpoint.species");
	for (std::size_t i = 0; i < species.size(); ++i)
		result.species.push_back(RequireString(species[i], "1d checkpoint.species["+std::to_string(i)+"]"));
	if (result.schema_version != 2 || result.completed_step < 0 || result.internal_substeps < 0
		|| result.physical_time < 0.0 || !std::isfinite(result.inlet_flow)
		|| !(result.dt > 0.0) || result.cells < 1 || result.nodes < 1 || result.segments < 1
		|| result.outlets < 1 || result.state_file.empty())
		throw std::runtime_error("1d checkpoint metadata contains invalid values");
	return result;
}

inline void WriteOneDCheckpoint(const std::filesystem::path& prefix,
	const OneDCheckpointMetadata& metadata, const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports, const OneDNetwork& network,
	int rank)
{
	const auto values = PackOneDCheckpointState(flow, transports, network);
	Vec state = nullptr;
	VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, static_cast<PetscInt>(values.size()), &state);
	OneDSetInitialVector(state, values);
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, OneDCheckpointStatePath(prefix).string().c_str(), FILE_MODE_WRITE, &viewer);
	VecView(state, viewer);
	PetscViewerDestroy(&viewer);
	VecDestroy(&state);
	int failed = 0;
	if (rank == 0) {
		try {
			std::ofstream output(OneDCheckpointMetadataPath(prefix));
			if (!output) throw std::runtime_error("cannot create 1d checkpoint metadata");
			output << SerializeOneDCheckpointMetadata(metadata);
			if (!output) throw std::runtime_error("cannot write 1d checkpoint metadata");
		} catch (const std::exception&) { failed = 1; }
	}
	MPI_Bcast(&failed, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	if (failed) throw std::runtime_error("cannot write 1d checkpoint metadata");
}

inline OneDCheckpointMetadata ReadOneDCheckpoint(const std::filesystem::path& prefix,
	OneDFlowState& flow, std::vector<OneDTransportState>& transports,
	OneDNetwork& network, double dynamic_viscosity)
{
	std::ifstream input(OneDCheckpointMetadataPath(prefix));
	if (!input) throw std::runtime_error("cannot open 1d checkpoint metadata");
	std::ostringstream contents; contents << input.rdbuf();
	auto metadata = ParseOneDCheckpointMetadata(contents.str());
	Vec state = nullptr;
	VecCreate(PETSC_COMM_WORLD, &state);
	PetscViewer viewer = nullptr;
	PetscViewerBinaryOpen(PETSC_COMM_WORLD, OneDCheckpointStatePath(prefix).string().c_str(), FILE_MODE_READ, &viewer);
	VecLoad(state, viewer);
	PetscViewerDestroy(&viewer);
	std::vector<double> values; OneDGetVectorAll(state, values);
	VecDestroy(&state);
	UnpackOneDCheckpointState(values, flow, transports, metadata, network, dynamic_viscosity);
	flow.completed_step = metadata.completed_step;
	flow.internal_substeps = metadata.internal_substeps;
	flow.physical_time = metadata.physical_time;
	flow.inlet_flow = metadata.inlet_flow;
	return metadata;
}

} // namespace iga

#endif

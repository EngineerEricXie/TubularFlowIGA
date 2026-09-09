#ifndef IGA_ONE_D_CHECKPOINT_HPP
#define IGA_ONE_D_CHECKPOINT_HPP

#include "OneDImplicit.hpp"
#include "OneDOutput.hpp"
#include "CollectiveFailure.hpp"
#include "CheckedText.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <exception>
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
	text.exceptions(std::ios::badbit | std::ios::failbit);
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

inline std::size_t OneDCheckpointValueCount(const OneDCheckpointMetadata& metadata)
{
	if (metadata.cells < 1 || metadata.nodes < 1 || metadata.segments < 1 || metadata.outlets < 1)
		throw std::runtime_error("1d checkpoint counts must be positive");
	const auto cells = static_cast<std::uint64_t>(metadata.cells);
	const auto limit = static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max());
	const auto fixed = 3*cells+static_cast<std::uint64_t>(metadata.nodes)
		+2*static_cast<std::uint64_t>(metadata.segments)+3*static_cast<std::uint64_t>(metadata.outlets);
	if (fixed > limit || metadata.species.size() > (limit-fixed)/cells)
		throw std::runtime_error("1d checkpoint state exceeds PETSc index capacity");
	const auto count = fixed+cells*metadata.species.size();
	if (count > std::numeric_limits<std::size_t>::max()/sizeof(double))
		throw std::runtime_error("1d checkpoint state exceeds addressable size");
	return static_cast<std::size_t>(count);
}

inline void UnpackOneDCheckpointState(const std::vector<double>& values,
	OneDFlowState& flow, std::vector<OneDTransportState>& transports,
	const OneDCheckpointMetadata& metadata, OneDNetwork& network,
	double dynamic_viscosity)
{
	const std::size_t expected = OneDCheckpointValueCount(metadata);
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
	output.exceptions(std::ios::badbit | std::ios::failbit);
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

namespace one_d_checkpoint_detail {

// Only used with COMM_SELF. Explicit Close reports errors; the destructor
// releases partially created objects during exception unwinding.
struct LocalPetscFile {
	LocalPetscFile()
	{
		OneDPetscCheck(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr), "checkpoint error handler");
	}
	LocalPetscFile(const LocalPetscFile&) = delete;
	LocalPetscFile& operator=(const LocalPetscFile&) = delete;
	~LocalPetscFile()
	{
		if (array) VecRestoreArrayRead(state, &array);
		if (viewer) PetscViewerDestroy(&viewer);
		if (state) VecDestroy(&state);
		PetscPopErrorHandler();
	}
	void Open(const std::string& path, PetscFileMode mode)
	{
		OneDPetscCheck(PetscViewerCreate(PETSC_COMM_SELF, &viewer), "checkpoint viewer create");
		OneDPetscCheck(PetscViewerSetType(viewer, PETSCVIEWERBINARY), "checkpoint viewer type");
		// A checkpoint has its own fixed filename/header contract. Do not let
		// a viewer option or an adjacent .info file silently redirect it.
		OneDPetscCheck(PetscViewerBinarySetSkipOptions(viewer, PETSC_TRUE), "checkpoint skip options");
		OneDPetscCheck(PetscViewerBinarySetSkipInfo(viewer, PETSC_TRUE), "checkpoint skip info");
		OneDPetscCheck(PetscViewerFileSetMode(viewer, mode), "checkpoint file mode");
		OneDPetscCheck(PetscViewerFileSetName(viewer, path.c_str()), "checkpoint file name");
	}
	void Close()
	{
		OneDPetscCheck(PetscViewerDestroy(&viewer), "checkpoint viewer close");
		OneDPetscCheck(VecDestroy(&state), "checkpoint vector destroy");
	}
	Vec state = nullptr;
	PetscViewer viewer = nullptr;
	const PetscScalar* array = nullptr;
};

inline void ValidateValues(const std::vector<double>& values, const OneDCheckpointMetadata& metadata)
{
	if (values.size() != OneDCheckpointValueCount(metadata))
		throw std::runtime_error("1d checkpoint state size is invalid");
	for (double value : values)
		if (!std::isfinite(value)) throw std::runtime_error("1d checkpoint state contains a nonfinite value");
	for (int cell = 0; cell < metadata.cells; ++cell)
		if (!(values[static_cast<std::size_t>(cell)] > 0.0))
			throw std::runtime_error("1d checkpoint state contains a nonpositive area");
}

inline std::string_view Bytes(const std::vector<double>& values)
{
	return {reinterpret_cast<const char*>(values.data()), values.size()*sizeof(double)};
}

} // namespace one_d_checkpoint_detail

// The native 1D runtime already holds complete state on each rank. File I/O
// uses COMM_SELF so a local PETSc I/O error cannot strand peers in a viewer
// collective. Group agreement is performed only after every local call returns.
inline void WriteOneDCheckpoint(const std::filesystem::path& prefix,
	const OneDCheckpointMetadata& metadata, const OneDFlowState& flow,
	const std::vector<OneDTransportState>& transports, const OneDNetwork& network,
	int rank, MPI_Comm communicator = PETSC_COMM_WORLD)
{
	std::vector<double> values;
	std::string text, state_path, metadata_path;
	CollectiveLocalStage(communicator, "1d checkpoint write preparation", [&] {
		int actual_rank = 0;
		MPI_Comm_rank(communicator, &actual_rank);
		if (rank != actual_rank) throw std::runtime_error("checkpoint writer rank does not match communicator");
		values = PackOneDCheckpointState(flow, transports, network);
		one_d_checkpoint_detail::ValidateValues(values, metadata);
		text = SerializeOneDCheckpointMetadata(metadata);
		ParseOneDCheckpointMetadata(text);
		state_path = OneDCheckpointStatePath(prefix).string();
		metadata_path = OneDCheckpointMetadataPath(prefix).string();
	});
	RequireCollectiveSameText(communicator, "1d checkpoint metadata agreement", text);
	RequireCollectiveSameText(communicator, "1d checkpoint state agreement", one_d_checkpoint_detail::Bytes(values));
	std::exception_ptr error;
	if (rank == 0) {
		try {
			one_d_checkpoint_detail::LocalPetscFile file;
			OneDPetscCheck(VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(values.size()), &file.state),
				"checkpoint vector create");
			OneDSetInitialVector(file.state, values);
			file.Open(state_path, FILE_MODE_WRITE);
			OneDPetscCheck(VecView(file.state, file.viewer), "checkpoint state write");
			file.Close();
		} catch (...) { error = std::current_exception(); }
	}
	CollectiveLocalStage(communicator, "1d checkpoint state write", [&] {
		if (error) std::rethrow_exception(error);
	});
	CollectiveLocalStage(communicator, "1d checkpoint metadata write", [&] {
		if (rank != 0) return;
		std::ofstream output(metadata_path);
		if (!output) throw std::runtime_error("cannot create 1d checkpoint metadata");
		output << text;
		output.close();
		if (!output) throw std::runtime_error("cannot write 1d checkpoint metadata");
	});
}

inline OneDCheckpointMetadata ReadOneDCheckpoint(const std::filesystem::path& prefix,
	OneDFlowState& flow, std::vector<OneDTransportState>& transports,
	OneDNetwork& network, double dynamic_viscosity, MPI_Comm communicator = PETSC_COMM_WORLD)
{
	OneDCheckpointMetadata metadata;
	std::string text, state_path;
	std::size_t count = 0;
	CollectiveLocalStage(communicator, "1d checkpoint metadata read", [&] {
		std::ifstream input(OneDCheckpointMetadataPath(prefix));
		if (!input) throw std::runtime_error("cannot open 1d checkpoint metadata");
		text = ReadCheckedText(input);
		metadata = ParseOneDCheckpointMetadata(text);
		count = OneDCheckpointValueCount(metadata);
		if (metadata.cells != network.cells || static_cast<std::size_t>(metadata.nodes) != network.nodes.size()
			|| static_cast<std::size_t>(metadata.segments) != network.segments.size()
			|| static_cast<std::size_t>(metadata.outlets) != flow.outlets.size()
			|| metadata.species != OneDCheckpointSpecies(transports))
			throw std::runtime_error("1d checkpoint layout does not match configuration");
		state_path = OneDCheckpointStatePath(prefix).string();
	});
	RequireCollectiveSameText(communicator, "1d checkpoint metadata agreement", text);
	std::vector<double> values;
	std::exception_ptr error;
	try {
		one_d_checkpoint_detail::LocalPetscFile file;
		OneDPetscCheck(VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(count), &file.state),
			"checkpoint vector create");
		file.Open(state_path, FILE_MODE_READ);
		OneDPetscCheck(VecLoad(file.state, file.viewer), "checkpoint state load");
		OneDPetscCheck(VecGetArrayRead(file.state, &file.array), "checkpoint vector read");
		values.resize(count);
		for (std::size_t i = 0; i < count; ++i) values[i] = PetscRealPart(file.array[i]);
		OneDPetscCheck(VecRestoreArrayRead(file.state, &file.array), "checkpoint vector restore");
		file.Close();
		one_d_checkpoint_detail::ValidateValues(values, metadata);
	} catch (...) { error = std::current_exception(); }
	CollectiveLocalStage(communicator, "1d checkpoint state read", [&] {
		if (error) std::rethrow_exception(error);
	});
	RequireCollectiveSameText(communicator, "1d checkpoint state agreement", one_d_checkpoint_detail::Bytes(values));
	CollectiveLocalStage(communicator, "1d checkpoint unpack", [&] {
		UnpackOneDCheckpointState(values, flow, transports, metadata, network, dynamic_viscosity);
		flow.completed_step = metadata.completed_step;
		flow.internal_substeps = metadata.internal_substeps;
		flow.physical_time = metadata.physical_time;
		flow.inlet_flow = metadata.inlet_flow;
	});
	return metadata;
}

} // namespace iga

#endif

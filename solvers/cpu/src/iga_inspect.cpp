#include "CheckedText.hpp"
#include "IgaDatabase.hpp"
#include "Sha256.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void AppendString(iga::Sha256& hash, const std::string& value)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(value.size()));
	hash.Append(value.data(), value.size());
}

std::string FileSha256(const std::string& path)
{
	std::ifstream input(path, std::ios::binary);
	if (!input) throw std::runtime_error("cannot open database for hashing: " + path);
	iga::Sha256 hash;
	std::array<char, 64*1024> buffer{};
	while (input) {
		input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
		if (input.gcount() > 0)
			hash.Append(buffer.data(), static_cast<std::size_t>(input.gcount()));
	}
	if (!input.eof()) throw std::runtime_error("cannot hash database: " + path);
	return hash.Hex();
}

bool IsRegionRole(const std::string& value)
{
	static const std::set<std::string> roles{
		"fluid", "solid", "shell", "porous", "network", "lumped"};
	return roles.count(value) != 0;
}

std::string JsonEscape(const std::string& value)
{
	std::ostringstream result;
	for (const unsigned char character : value) {
		if (character == '"' || character == '\\') result << '\\' << character;
		else if (character == '\n') result << "\\n";
		else if (character < 0x20)
			result << "\\u" << std::hex << std::setw(4) << std::setfill('0')
				<< static_cast<unsigned>(character) << std::dec;
		else result << character;
	}
	return result.str();
}

void WriteGeometryManifest(const std::string& path, const std::string& database_path,
	const std::string& role, const iga::Header& header, const std::string& identity,
	const std::string& artifact_sha256, const std::map<std::int32_t, std::uint64_t>& boundaries)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create geometry manifest: " + path);
	output << std::setprecision(17)
		<< "{\n  \"schema_version\": 1,\n"
		<< "  \"route\": \"centerline_to_iga_volume\",\n"
		<< "  \"coordinate_system\": \"cartesian\",\n"
		<< "  \"length_unit\": \"m\",\n"
		<< "  \"reference_geometry\": {\"identity_sha256\": \"" << identity
		<< "\", \"source_artifact\": \"" << JsonEscape(database_path)
		<< "\", \"source_artifact_sha256\": \"" << artifact_sha256 << "\"},\n"
		<< "  \"current_geometry\": {\"kind\": \"same_as_reference\", "
		<< "\"identity_sha256\": \"" << identity << "\"},\n"
		<< "  \"source_to_solver_transform\": {\"source_origin\": ["
		<< header.geometry_transform.source_origin[0] << ", "
		<< header.geometry_transform.source_origin[1] << ", "
		<< header.geometry_transform.source_origin[2] << "], "
		<< "\"source_units_per_normalized_unit\": "
		<< header.geometry_transform.source_units_per_normalized_unit << ", "
		<< "\"source_length_scale_to_m\": "
		<< header.geometry_transform.source_length_scale_to_m << "},\n"
		<< "  \"stable_ids\": {\"volume_elements\": \"ntiga_element_id\", "
		<< "\"control_nodes\": \"ntiga_connectivity_id\", "
		<< "\"boundary_faces\": \"ntiga_element_id_and_local_face\"},\n"
		<< "  \"regions\": [{\"id\": \"volume\", \"role\": \"" << role
		<< "\", \"dimension\": 3}],\n  \"boundaries\": [";
	std::size_t index = 0;
	for (const auto& entry : boundaries) {
		if (index++) output << ",";
		output << "{\"label\": " << entry.first
			<< ", \"faces\": " << entry.second
			<< ", \"semantic_role\": \"explicitly_configured_elsewhere\"}";
	}
	output << "]\n}\n";
	output.close();
	if (!output) throw std::runtime_error("cannot finalize geometry manifest: " + path);
}

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc != 2 && argc != 6) {
			std::cerr << "usage: iga_inspect DATABASE.ntiga "
				"[--geometry-manifest FILE.json --region-role ROLE]\n";
			return 2;
		}
		std::string manifest_path;
		std::string region_role;
		if (argc == 6) {
			if (std::string(argv[2]) != "--geometry-manifest"
				|| std::string(argv[4]) != "--region-role") {
				std::cerr << "usage: iga_inspect DATABASE.ntiga "
					"[--geometry-manifest FILE.json --region-role ROLE]\n";
				return 2;
			}
			manifest_path = argv[3];
			region_role = argv[5];
			if (!IsRegionRole(region_role))
				throw std::invalid_argument("region role is not part of the geometry contract");
		}
		iga::Database db(argv[1]);
		if (!manifest_path.empty() && db.header().version < iga::kVersion)
			throw std::runtime_error(
				"geometry manifest requires a version-5 database with an explicit transform");
		iga::Sha256 geometry_hash;
		if (!manifest_path.empty()) {
			AppendString(geometry_hash, "TubularFlowIGA/body-fitted-geometry/v1");
			geometry_hash.AppendLittleEndian64(db.header().elements);
			geometry_hash.AppendLittleEndian64(db.header().nodes);
			for (const double value : db.header().geometry_transform.source_origin)
				geometry_hash.AppendNormalizedDouble(value);
			geometry_hash.AppendNormalizedDouble(
				db.header().geometry_transform.source_units_per_normalized_unit);
			geometry_hash.AppendNormalizedDouble(
				db.header().geometry_transform.source_length_scale_to_m);
		}
		std::vector<std::uint64_t> counts(db.header().ranks, 0);
		for (auto owner : db.owners()) ++counts.at(static_cast<std::size_t>(owner));
		std::uint64_t total_basis = 0;
		std::size_t max_basis = 0;
		std::uint64_t total_required = 0;
		std::uint64_t min_required = std::numeric_limits<std::uint64_t>::max();
		std::uint64_t max_required = 0;
		std::map<std::int32_t, std::uint64_t> boundary_faces;
		for (std::uint32_t rank = 0; rank < db.header().ranks; ++rank) {
			const auto n = db.RequiredElementIndices(static_cast<std::int32_t>(rank)).size();
			total_required += n;
			min_required = std::min(min_required, static_cast<std::uint64_t>(n));
			max_required = std::max(max_required, static_cast<std::uint64_t>(n));
		}
		for (std::uint64_t i = 0; i < db.header().elements; ++i) {
			auto element = db.Load(i);
			if (!manifest_path.empty()) {
				geometry_hash.AppendLittleEndian64(element.id);
				geometry_hash.AppendLittleEndian32(static_cast<std::uint32_t>(element.type));
				for (const auto label : element.boundary_labels)
					geometry_hash.AppendLittleEndian32(static_cast<std::uint32_t>(label));
				geometry_hash.AppendLittleEndian64(
					static_cast<std::uint64_t>(element.connectivity.size()));
				for (const auto node : element.connectivity)
					geometry_hash.AppendLittleEndian32(static_cast<std::uint32_t>(node));
				for (const auto& row : element.extraction)
					for (const double value : row) geometry_hash.AppendNormalizedDouble(value);
				for (const auto& point : element.bezier_points)
					for (const double value : point) geometry_hash.AppendNormalizedDouble(value);
			}
			total_basis += element.connectivity.size();
			max_basis = std::max(max_basis, element.connectivity.size());
			for (const auto label : element.boundary_labels)
				if (label >= 0) ++boundary_faces[label];
		}
		auto [minimum, maximum] = std::minmax_element(counts.begin(), counts.end());
		std::cout << "elements: " << db.header().elements << '\n'
			<< "nodes: " << db.header().nodes << '\n'
			<< "ranks: " << db.header().ranks << '\n'
			<< "source origin: [" << db.header().geometry_transform.source_origin[0] << ", "
			<< db.header().geometry_transform.source_origin[1] << ", "
			<< db.header().geometry_transform.source_origin[2] << "]\n"
			<< "source units/normalized unit: "
			<< db.header().geometry_transform.source_units_per_normalized_unit << '\n'
			<< "source length scale to m: "
			<< db.header().geometry_transform.source_length_scale_to_m << '\n'
			<< "basis/element mean: " << static_cast<double>(total_basis) / db.header().elements << '\n'
			<< "basis/element max: " << max_basis << '\n'
			<< "elements/rank min: " << *minimum << '\n'
			<< "elements/rank max: " << *maximum << '\n';
		std::cout << "row-touching elements/rank min: " << min_required << '\n'
			<< "row-touching elements/rank max: " << max_required << '\n'
			<< "row-touching duplication: " << static_cast<double>(total_required) / db.header().elements << "x\n";
		for (const auto& entry : boundary_faces)
			std::cout << "boundary_faces[" << entry.first << "]: " << entry.second << '\n';
		if (!manifest_path.empty()) {
			const auto geometry_identity = geometry_hash.Hex();
			std::cout << "geometry identity sha256: " << geometry_identity << '\n';
			WriteGeometryManifest(manifest_path, argv[1], region_role, db.header(),
				geometry_identity, FileSha256(argv[1]), boundary_faces);
		}
		iga::FlushCheckedText(std::cout);
		return 0;
	} catch (const std::exception& e) {
		std::cerr << "iga_inspect: " << e.what() << '\n';
		return 1;
	}
}

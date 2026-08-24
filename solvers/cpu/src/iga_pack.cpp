#include "IgaDatabase.hpp"
#include "IgaPreprocessCache.hpp"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

namespace fs = std::filesystem;

namespace {

struct InputElement
{
	std::uint64_t id = 0;
	std::int32_t type = 0;
	std::vector<std::int32_t> connectivity;
	std::vector<std::array<double, iga::kBezierPointCount>> extraction;
	std::array<std::array<double, 3>, iga::kBezierPointCount> bezier_points{};
};

template <class T>
T Next(std::istream& in, const char* file, std::uint64_t element, const char* field)
{
	T value{};
	if (!(in >> value)) {
		std::ostringstream msg;
		msg << file << " is truncated or malformed at element " << element << " (" << field << ')';
		throw std::runtime_error(msg.str());
	}
	return value;
}

std::vector<std::int32_t> ReadPartition(const fs::path& path, std::uint64_t elements, std::uint32_t ranks)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open partition file: " + path.string());
	std::vector<std::int32_t> owners;
	owners.reserve(static_cast<std::size_t>(elements));
	std::int64_t owner = 0;
	while (in >> owner) {
		if (owner < 0 || owner >= ranks) throw std::runtime_error("partition contains an invalid rank");
		owners.push_back(static_cast<std::int32_t>(owner));
	}
	if (owners.size() != elements)
		throw std::runtime_error("partition element count does not match preprocessing data");
	return owners;
}

void RequireOnlyWhitespace(std::istream& in, const char* file)
{
	std::string extra;
	if (in >> extra) throw std::runtime_error(std::string(file) + " has extra data after the declared elements");
}

void RequireCacheEnd(std::istream& in)
{
	char extra = 0;
	if (in.read(&extra, 1)) throw std::runtime_error("preprocessing cache has extra trailing data");
	if (!in.eof()) throw std::runtime_error("failed while finalizing preprocessing cache read");
}

struct ControlMesh
{
	std::uint64_t nodes = 0;
	std::vector<std::array<std::int32_t, 8>> cells;
	std::vector<std::int32_t> labels;
};

ControlMesh ReadControlMesh(const fs::path& path)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open mesh file: " + path.string());
	ControlMesh mesh;
	std::string token;
	while (in >> token) {
		if (token == "POINTS") {
			std::string type;
			if (!(in >> mesh.nodes >> type) || mesh.nodes == 0)
				throw std::runtime_error("invalid POINTS record in controlmesh.vtk");
			double coordinate = 0.0;
			for (std::uint64_t i = 0; i < 3*mesh.nodes; ++i)
				if (!(in >> coordinate)) throw std::runtime_error("truncated POINTS data in controlmesh.vtk");
		} else if (token == "CELLS") {
			std::uint64_t cells = 0, entries = 0;
			if (!(in >> cells >> entries) || cells == 0 || entries != 9*cells)
				throw std::runtime_error("invalid CELLS record in controlmesh.vtk");
			mesh.cells.resize(static_cast<std::size_t>(cells));
			for (auto& cell : mesh.cells) {
				int count = 0;
				if (!(in >> count) || count != 8) throw std::runtime_error("controlmesh.vtk requires hexahedral cells");
				for (auto& node : cell) {
					std::int64_t value = -1;
					if (!(in >> value) || value < 0 || static_cast<std::uint64_t>(value) >= mesh.nodes)
						throw std::runtime_error("invalid cell connectivity in controlmesh.vtk");
					node = static_cast<std::int32_t>(value);
				}
			}
		} else if (token == "POINT_DATA") {
			std::uint64_t count = 0;
			std::string scalars, name, type, lookup, table;
			int components = 0;
			if (!(in >> count >> scalars >> name >> type >> components >> lookup >> table)
				|| count != mesh.nodes || scalars != "SCALARS" || name != "label"
				|| components != 1 || lookup != "LOOKUP_TABLE")
				throw std::runtime_error("invalid point-label data in controlmesh.vtk");
			mesh.labels.resize(static_cast<std::size_t>(mesh.nodes));
			for (auto& label : mesh.labels)
				if (!(in >> label)) throw std::runtime_error("truncated point labels in controlmesh.vtk");
			break;
		}
	}
	if (mesh.nodes == 0 || mesh.cells.empty() || mesh.labels.size() != mesh.nodes)
		throw std::runtime_error("controlmesh.vtk is missing points, cells, or labels");
	return mesh;
}

std::vector<std::array<std::int32_t, 6>> BoundaryFaceLabels(const ControlMesh& mesh)
{
	constexpr int face_nodes[6][4] = {
		{0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5},
		{2, 3, 7, 6}, {0, 4, 7, 3}, {4, 5, 6, 7}
	};
	using FaceKey = std::array<std::int32_t, 4>;
	using FaceLocation = std::pair<std::size_t, std::size_t>;
	std::map<FaceKey, std::vector<FaceLocation>> incidence;
	for (std::size_t element = 0; element < mesh.cells.size(); ++element)
		for (std::size_t face = 0; face < 6; ++face) {
			FaceKey key{};
			for (std::size_t corner = 0; corner < 4; ++corner)
				key[corner] = mesh.cells[element][face_nodes[face][corner]];
			std::sort(key.begin(), key.end());
			incidence[key].push_back({element, face});
		}
	std::vector<std::array<std::int32_t, 6>> result(mesh.cells.size());
	for (auto& labels : result) labels.fill(-1);
	for (const auto& item : incidence) {
		if (item.second.size() == 2) continue;
		if (item.second.size() != 1) throw std::runtime_error("non-manifold control-mesh face");
		const auto [element, face] = item.second.front();
		std::int32_t positive = -1;
		bool all_wall = true;
		for (std::size_t corner = 0; corner < 4; ++corner) {
			const auto node = mesh.cells[element][face_nodes[face][corner]];
			const auto label = mesh.labels[static_cast<std::size_t>(node)];
			all_wall = all_wall && label == 0;
			if (label > 0) {
				if (positive >= 0 && positive != label)
					throw std::runtime_error("boundary face contains multiple positive labels");
				positive = label;
			}
		}
		if (positive >= 0) result[element][face] = positive;
		else if (all_wall) result[element][face] = 0;
		else throw std::runtime_error("boundary face has no unambiguous mesh label");
	}
	return result;
}

std::uint32_t NodeOwner(std::uint64_t node, std::uint64_t nodes, std::uint32_t ranks)
{
	if (node >= nodes) throw std::runtime_error("connectivity exceeds controlmesh point count");
	const auto q = nodes / ranks;
	const auto rem = nodes % ranks;
	const auto wide_end = (q + 1) * rem;
	if (node < wide_end) return static_cast<std::uint32_t>(node / (q + 1));
	if (q == 0) throw std::runtime_error("more ranks than mesh nodes");
	return static_cast<std::uint32_t>(rem + (node - wide_end) / q);
}

InputElement ReadLegacyElement(std::istream& cmat, std::istream& bzpt, std::uint64_t expected)
{
	InputElement element;
	const auto id = Next<std::uint64_t>(cmat, "cmat.txt", expected, "id");
	const auto nen64 = Next<std::uint64_t>(cmat, "cmat.txt", expected, "basis count");
	const auto type64 = Next<std::int64_t>(cmat, "cmat.txt", expected, "type");
	if (id != expected) throw std::runtime_error("cmat.txt element ids are not contiguous");
	if (nen64 == 0 || nen64 > 4096) throw std::runtime_error("invalid basis count in cmat.txt");
	if (type64 < std::numeric_limits<std::int32_t>::min() ||
		type64 > std::numeric_limits<std::int32_t>::max())
		throw std::runtime_error("element type is out of range");
	element.id = id;
	element.type = static_cast<std::int32_t>(type64);
	const auto nen = static_cast<std::uint32_t>(nen64);
	element.connectivity.resize(nen);
	element.extraction.resize(nen);
	for (auto& node : element.connectivity) {
		const auto input = Next<std::int64_t>(cmat, "cmat.txt", expected, "connectivity");
		if (input < 0 || input > std::numeric_limits<std::int32_t>::max())
			throw std::runtime_error("connectivity index is out of range");
		node = static_cast<std::int32_t>(input);
	}
	for (auto& row : element.extraction)
		for (double& value : row)
			value = Next<double>(cmat, "cmat.txt", expected, "extraction coefficient");
	for (auto& point : element.bezier_points)
		for (double& value : point)
			value = Next<double>(bzpt, "bzpt.txt", expected, "Bezier point");
	return element;
}

InputElement ReadCacheElement(std::istream& cache, std::uint64_t expected)
{
	InputElement element;
	std::uint32_t nen = 0;
	igacache::Read(cache, element.id);
	igacache::Read(cache, element.type);
	igacache::Read(cache, nen);
	if (element.id != expected) throw std::runtime_error("preprocessing cache element ids are not contiguous");
	if (nen == 0 || nen > 4096) throw std::runtime_error("invalid basis count in preprocessing cache");
	element.connectivity.resize(nen);
	element.extraction.resize(nen);
	for (auto& node : element.connectivity) {
		igacache::Read(cache, node);
		if (node < 0) throw std::runtime_error("negative connectivity in preprocessing cache");
	}
	for (auto& row : element.extraction) {
		row.fill(0.0);
		std::array<bool, iga::kBezierPointCount> seen{};
		std::uint8_t nonzeros = 0;
		igacache::Read(cache, nonzeros);
		if (nonzeros > iga::kBezierPointCount)
			throw std::runtime_error("too many extraction entries in preprocessing cache row");
		for (std::uint8_t entry = 0; entry < nonzeros; entry++) {
			std::uint8_t column = 0;
			double value = 0.0;
			igacache::Read(cache, column);
			igacache::Read(cache, value);
			if (column >= iga::kBezierPointCount || seen[column] || value == 0.0)
				throw std::runtime_error("invalid sparse extraction entry in preprocessing cache");
			seen[column] = true;
			row[column] = value;
		}
	}
	for (auto& point : element.bezier_points)
		for (double& value : point)
			igacache::Read(cache, value);
	return element;
}

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc < 4 || argc > 5) {
			std::cerr << "usage: iga_pack CASE_DIR RANKS OUTPUT.ntiga [--legacy-text]\n";
			return 2;
		}
		const bool force_legacy = argc == 5;
		if (force_legacy && std::string(argv[4]) != "--legacy-text") {
			std::cerr << "usage: iga_pack CASE_DIR RANKS OUTPUT.ntiga [--legacy-text]\n";
			return 2;
		}
		const fs::path dir(argv[1]);
		const auto ranks64 = std::stoull(argv[2]);
		if (ranks64 == 0 || ranks64 > std::numeric_limits<std::uint32_t>::max())
			throw std::runtime_error("invalid rank count");
		const auto ranks = static_cast<std::uint32_t>(ranks64);
		const fs::path output(argv[3]);
		const auto control_mesh = ReadControlMesh(dir / "controlmesh.vtk");
		const auto nodes = control_mesh.nodes;

		const fs::path cache_path = dir / "spline_cache.igacache";
		const bool use_cache = !force_legacy && fs::exists(cache_path);
		std::ifstream cache;
		std::ifstream cmat;
		std::ifstream bzpt;
		std::uint64_t elements = 0;
		if (use_cache) {
			cache.open(cache_path, std::ios::binary);
			if (!cache) throw std::runtime_error("cannot open " + cache_path.string());
			const auto header = igacache::ReadHeader(cache);
			if (header.nodes != nodes)
				throw std::runtime_error("preprocessing cache node count does not match controlmesh.vtk");
			if (header.mesh_hash != igacache::HashFile((dir / "controlmesh.vtk").string()))
				throw std::runtime_error("preprocessing cache does not match controlmesh.vtk content");
			elements = header.elements;
			std::cout << "using sparse preprocessing cache " << cache_path << '\n';
		} else {
			cmat.open(dir / "cmat.txt");
			bzpt.open(dir / "bzpt.txt");
			if (!cmat) throw std::runtime_error("cannot open " + (dir / "cmat.txt").string());
			if (!bzpt) throw std::runtime_error("cannot open " + (dir / "bzpt.txt").string());
			elements = Next<std::uint64_t>(cmat, "cmat.txt", 0, "element count");
			const auto point_values = Next<std::uint64_t>(bzpt, "bzpt.txt", 0, "point count");
			if (elements == 0) throw std::runtime_error("cmat.txt declares zero elements");
			if (point_values != elements * iga::kBezierPointCount)
				throw std::runtime_error("bzpt.txt point count is not 64 times the element count");
			std::cout << "using legacy cmat.txt and bzpt.txt preprocessing data\n";
		}

		if (control_mesh.cells.size() != elements)
			throw std::runtime_error("controlmesh cell count does not match preprocessing data");
		const auto boundary_faces = BoundaryFaceLabels(control_mesh);

		const auto partition = ReadPartition(
			dir / ("bzmeshinfo.txt.epart." + std::to_string(ranks)), elements, ranks);
		std::ofstream out(output, std::ios::binary | std::ios::trunc);
		if (!out) throw std::runtime_error("cannot create " + output.string());

		out.write(iga::kMagic.data(), iga::kMagic.size());
		iga::Write(out, iga::kVersion);
		iga::Write(out, ranks);
		iga::Write(out, elements);
		iga::Write(out, nodes);
		iga::Write(out, iga::kBezierPointCount);
		iga::Write(out, std::uint32_t{0});
		const auto rank_index_header_position = out.tellp();
		iga::Write(out, std::uint64_t{0});
		const auto index_position = out.tellp();
		std::vector<std::uint64_t> offsets(static_cast<std::size_t>(elements + 1), 0);
		out.write(reinterpret_cast<const char*>(offsets.data()),
			static_cast<std::streamsize>(offsets.size() * sizeof(std::uint64_t)));
		out.write(reinterpret_cast<const char*>(partition.data()),
			static_cast<std::streamsize>(partition.size() * sizeof(std::int32_t)));
		std::vector<std::vector<std::uint64_t>> required(ranks);

		for (std::uint64_t e = 0; e < elements; ++e) {
			offsets[static_cast<std::size_t>(e)] = static_cast<std::uint64_t>(out.tellp());
			const InputElement element = use_cache ? ReadCacheElement(cache, e) : ReadLegacyElement(cmat, bzpt, e);
			const auto nen = static_cast<std::uint32_t>(element.connectivity.size());
			iga::Write(out, element.id);
			iga::Write(out, element.type);
			iga::Write(out, partition[static_cast<std::size_t>(e)]);
			iga::Write(out, nen);
			out.write(reinterpret_cast<const char*>(boundary_faces[static_cast<std::size_t>(e)].data()),
				static_cast<std::streamsize>(sizeof(boundary_faces[static_cast<std::size_t>(e)])));
			std::set<std::uint32_t> row_owners;
			for (const auto node : element.connectivity) {
				iga::Write(out, node);
				row_owners.insert(NodeOwner(static_cast<std::uint64_t>(node), nodes, ranks));
			}
			for (auto owner : row_owners) required[owner].push_back(e);
			for (const auto& row : element.extraction) {
				std::uint8_t nonzeros = 0;
				for (const double value : row) if (value != 0.0) ++nonzeros;
				iga::Write(out, nonzeros);
				for (std::uint8_t column = 0; column < iga::kBezierPointCount; ++column) {
					if (row[column] == 0.0) continue;
					iga::Write(out, column);
					iga::Write(out, row[column]);
				}
			}
			out.write(reinterpret_cast<const char*>(element.bezier_points.data()),
				static_cast<std::streamsize>(sizeof(element.bezier_points)));
			if (!out) throw std::runtime_error("failed to write element record");
			if ((e + 1) % 10000 == 0) std::cerr << "packed " << (e + 1) << '/' << elements << " elements\n";
		}
		offsets.back() = static_cast<std::uint64_t>(out.tellp());
		if (use_cache) RequireCacheEnd(cache);
		else {
			RequireOnlyWhitespace(cmat, "cmat.txt");
			RequireOnlyWhitespace(bzpt, "bzpt.txt");
		}
		const auto rank_index_offset = static_cast<std::uint64_t>(out.tellp());
		std::vector<std::uint64_t> rank_offsets(static_cast<std::size_t>(ranks) + 1, 0);
		for (std::uint32_t rank = 0; rank < ranks; ++rank)
			rank_offsets[rank + 1] = rank_offsets[rank] + required[rank].size();
		out.write(reinterpret_cast<const char*>(rank_offsets.data()),
			static_cast<std::streamsize>(rank_offsets.size() * sizeof(std::uint64_t)));
		for (const auto& indices : required)
			out.write(reinterpret_cast<const char*>(indices.data()),
				static_cast<std::streamsize>(indices.size() * sizeof(std::uint64_t)));
		out.seekp(index_position);
		out.write(reinterpret_cast<const char*>(offsets.data()),
			static_cast<std::streamsize>(offsets.size() * sizeof(std::uint64_t)));
		if (!out) throw std::runtime_error("failed to finalize database index");
		out.seekp(rank_index_header_position);
		iga::Write(out, rank_index_offset);
		out.close();
		std::cout << "wrote " << output << " with " << elements << " validated elements for " << ranks << " ranks\n";
		return 0;
	} catch (const std::exception& e) {
		std::cerr << "iga_pack: " << e.what() << '\n';
		return 1;
	}
}

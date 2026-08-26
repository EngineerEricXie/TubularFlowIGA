#ifndef BOUNDARY_SUPPORT_HPP
#define BOUNDARY_SUPPORT_HPP

#include "IgaDatabase.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

using HexCell = std::array<std::int32_t, 8>;

struct LabeledHexMesh {
	std::vector<HexCell> cells;
	std::vector<int> labels;
};

inline LabeledHexMesh ReadLabeledHexMesh(const std::string& path,
	std::uint64_t expected_nodes, std::uint64_t expected_elements)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open mesh: " + path);
	std::string token;
	std::uint64_t nodes = 0;
	while (in >> token) if (token == "POINTS") {
		std::string data_type;
		if (!(in >> nodes >> data_type)) throw std::runtime_error("invalid POINTS record");
		break;
	}
	if (nodes != expected_nodes) throw std::runtime_error("VTK point count does not match database nodes");
	for (std::uint64_t value = 0; value < 3*nodes; ++value) {
		double ignored = 0.0;
		if (!(in >> ignored)) throw std::runtime_error("VTK POINTS record is truncated");
	}
	std::uint64_t elements = 0, cell_values = 0;
	while (in >> token) if (token == "CELLS") {
		if (!(in >> elements >> cell_values)) throw std::runtime_error("invalid CELLS record");
		break;
	}
	if (elements != expected_elements) throw std::runtime_error("VTK cell count does not match database elements");
	LabeledHexMesh mesh;
	mesh.cells.resize(static_cast<std::size_t>(elements));
	for (auto& cell : mesh.cells) {
		int count = 0;
		if (!(in >> count) || count != 8) throw std::runtime_error("only eight-node VTK hexahedra are supported");
		for (auto& node : cell)
			if (!(in >> node) || node < 0 || static_cast<std::uint64_t>(node) >= nodes)
				throw std::runtime_error("invalid hexahedron connectivity");
	}
	std::uint64_t point_data = 0;
	while (in >> token) if (token == "POINT_DATA") {
		if (!(in >> point_data)) throw std::runtime_error("invalid POINT_DATA record");
		break;
	}
	if (point_data != nodes) throw std::runtime_error("VTK POINT_DATA count does not match database nodes");
	while (in >> token) if (token == "LOOKUP_TABLE") { in >> token; break; }
	if (!in) throw std::runtime_error("VTK point labels were not found");
	mesh.labels.resize(static_cast<std::size_t>(nodes));
	for (auto& label : mesh.labels) {
		double value = 0.0;
		if (!(in >> value)) throw std::runtime_error("VTK point-label array is truncated");
		label = static_cast<int>(value);
		if (value != static_cast<double>(label)) throw std::runtime_error("VTK point label is not integral");
	}
	return mesh;
}

struct BoundaryFace {
	std::uint64_t element = 0;
	int local_face = -1;
	std::array<std::int32_t, 4> nodes{};
	int uses = 0;
};

inline std::vector<BoundaryFace> ExternalFaces(const LabeledHexMesh& mesh)
{
	constexpr std::array<std::array<int, 4>, 6> corners{{
		{{0, 3, 2, 1}}, {{0, 1, 5, 4}}, {{1, 2, 6, 5}},
		{{2, 3, 7, 6}}, {{0, 4, 7, 3}}, {{4, 5, 6, 7}}
	}};
	std::map<std::array<std::int32_t, 4>, BoundaryFace> faces;
	for (std::size_t element = 0; element < mesh.cells.size(); ++element)
		for (int local_face = 0; local_face < 6; ++local_face) {
			std::array<std::int32_t, 4> nodes{};
			for (int corner = 0; corner < 4; ++corner)
				nodes[corner] = mesh.cells[element][corners[local_face][corner]];
			auto key = nodes;
			std::sort(key.begin(), key.end());
			auto& face = faces[key];
			if (face.uses == 0) {
				face.element = element;
				face.local_face = local_face;
				face.nodes = nodes;
			}
			if (++face.uses > 2) throw std::runtime_error("non-manifold control-mesh face");
		}
	std::vector<BoundaryFace> result;
	for (const auto& item : faces)
		if (item.second.uses == 1) result.push_back(item.second);
	return result;
}

inline bool IsWallFace(const BoundaryFace& face, const std::vector<int>& labels)
{
	return std::all_of(face.nodes.begin(), face.nodes.end(), [&](std::int32_t node) {
		return labels.at(static_cast<std::size_t>(node)) == 0;
	});
}

inline std::array<int, 16> FaceBezierColumns(int local_face)
{
	std::array<int, 16> columns{};
	int count = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const bool on_face = (local_face == 0 && k == 0)
					|| (local_face == 1 && j == 0) || (local_face == 2 && i == 3)
					|| (local_face == 3 && j == 3) || (local_face == 4 && i == 0)
					|| (local_face == 5 && k == 3);
				if (on_face) columns[static_cast<std::size_t>(count++)] = i + 4*j + 16*k;
			}
	if (count != 16) throw std::logic_error("invalid Bezier face-column count");
	return columns;
}

inline std::set<std::int32_t> WallTraceBasis(Database& database, const LabeledHexMesh& mesh)
{
	std::set<std::int32_t> result;
	for (const auto& face : ExternalFaces(mesh)) {
		if (!IsWallFace(face, mesh.labels)) continue;
		const auto element = database.Load(face.element);
		if (element.id != face.element)
			throw std::runtime_error("database element IDs do not match VTK cell order");
		for (std::size_t row = 0; row < element.extraction.size(); ++row)
			for (int column : FaceBezierColumns(face.local_face))
				if (element.extraction[row][static_cast<std::size_t>(column)] != 0.0) {
					result.insert(element.connectivity[row]);
					break;
				}
	}
	return result;
}

} // namespace iga

#endif

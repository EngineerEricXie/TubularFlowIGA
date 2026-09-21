#ifndef NATIVE_TET_MESH_COMPONENTS_HPP
#define NATIVE_TET_MESH_COMPONENTS_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

inline std::size_t CountNativeTetFaceComponents(const NativeTetMesh& mesh)
{
	if (mesh.cells.empty())
		throw std::invalid_argument("native tetra mesh has no cells");
	std::vector<std::size_t> parent(mesh.cells.size());
	std::iota(parent.begin(), parent.end(), 0);
	auto find = [&](std::size_t index) {
		while (parent[index] != index) {
			parent[index] = parent[parent[index]];
			index = parent[index];
		}
		return index;
	};
	std::map<std::array<std::uint32_t, 3>,
		std::pair<std::size_t, unsigned>> face_owners;
	std::set<std::array<std::uint32_t, 4>> unique_tetrahedra;
	std::set<std::uint64_t> unique_cell_ids;
	for (std::size_t index = 0; index < mesh.cells.size(); ++index) {
		const auto& cell = mesh.cells[index];
		const auto& nodes = cell.nodes;
		for (const auto node : nodes)
			if (node >= mesh.points.size())
				throw std::invalid_argument("native tetra node index is out of range");
		auto unique_nodes = nodes;
		std::sort(unique_nodes.begin(), unique_nodes.end());
		if (std::adjacent_find(unique_nodes.begin(), unique_nodes.end()) != unique_nodes.end())
			throw std::invalid_argument("native tetra has repeated node index");
		if (!unique_tetrahedra.insert(unique_nodes).second)
			throw std::invalid_argument("native tetra geometry is duplicated");
		if (!unique_cell_ids.insert(cell.id).second)
			throw std::invalid_argument("native tetra cell ID is duplicated");
		for (std::size_t omitted = 0; omitted < 4; ++omitted) {
			std::array<std::uint32_t, 3> face{};
			std::size_t face_node = 0;
			for (std::size_t local = 0; local < 4; ++local)
				if (local != omitted) face[face_node++] = nodes[local];
			std::sort(face.begin(), face.end());
			const auto inserted = face_owners.emplace(face,
				std::make_pair(index, 1u));
			if (!inserted.second) {
				if (++inserted.first->second.second > 2)
					throw std::invalid_argument("native tetra face has more than two owners");
				const auto first = find(inserted.first->second.first);
				const auto second = find(index);
				parent[second] = first;
			}
		}
	}
	std::size_t components = 0;
	for (std::size_t index = 0; index < parent.size(); ++index)
		if (find(index) == index) ++components;
	return components;
}

inline void RequireNativeTetSingleFaceComponent(const NativeTetMesh& mesh,
		const char* role)
{
	const auto count = CountNativeTetFaceComponents(mesh);
	if (count != 1)
		throw std::invalid_argument(std::string(role)+" tetra mesh has "
			+std::to_string(count)+" disconnected face components");
}

} // namespace iga

#endif

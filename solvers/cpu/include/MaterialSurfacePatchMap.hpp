#ifndef IGA_MATERIAL_SURFACE_PATCH_MAP_HPP
#define IGA_MATERIAL_SURFACE_PATCH_MAP_HPP

// Immutable, explicit P1 subset binding from a structural patch layout to an
// otherwise closed material surface.  It deliberately has no coordinate
// lookup authority: callers supply every node/triangle correspondence and we
// prove it against directed source topology.
#include "MaterialSurfaceKinematics.hpp"
#include "DistributedSurfaceInterface.hpp"

#include <algorithm>
#include <map>
#include <limits>
#include <set>
#include <stdexcept>

namespace iga {

class MaterialSurfacePatchMap {
public:
	using GlobalToSourceVertex = std::pair<std::uint64_t, std::uint32_t>;
	static std::string BuildReferenceIdentitySha256(const MaterialSurfaceKinematics& full_reference,
		std::uint32_t patch_label, const std::vector<GlobalToSourceVertex>& global_to_source_vertex,
		const std::vector<std::array<std::uint64_t, 3>>& layout_reference_triangles,
		const std::vector<std::uint32_t>& layout_triangle_to_source_triangle,
		const std::vector<std::uint64_t>& configured_clamped_global_node_ids)
	{
		ValidateReferenceInputs(full_reference, patch_label, global_to_source_vertex,
			layout_reference_triangles, layout_triangle_to_source_triangle, configured_clamped_global_node_ids);
		return HashPatchReference(full_reference, patch_label, global_to_source_vertex,
			layout_reference_triangles, layout_triangle_to_source_triangle);
	}

	static MaterialSurfacePatchMap Create(DistributedSurfaceInterface patch_interface,
		DistributedSurfaceLayout patch_layout, MaterialSurfaceKinematics full_reference,
		std::uint32_t patch_label, std::vector<GlobalToSourceVertex> global_to_source_vertex,
		std::vector<std::uint32_t> layout_triangle_to_source_triangle,
		std::vector<std::uint64_t> configured_clamped_global_node_ids)
	{
		ValidateDistributedSurfaceInterface(patch_interface);
		ValidateDistributedSurfaceLayout(patch_layout);
		full_reference.Validate();
		ValidateReferenceInputs(full_reference, patch_label, global_to_source_vertex,
			patch_layout.reference_triangles, layout_triangle_to_source_triangle, configured_clamped_global_node_ids);
		if (patch_interface.provides != std::vector<SurfaceFieldQuantity>{
			SurfaceFieldQuantity::Displacement, SurfaceFieldQuantity::Velocity}
			|| patch_interface.requires != std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::TractionOnStructure})
			throw std::runtime_error("material patch map requires the exact structural displacement/velocity and traction role");
		for (std::size_t vertex = 0; vertex < full_reference.ReferenceMaterialVerticesM().size(); ++vertex)
			if (!EqualNormalized(full_reference.ReferenceMaterialVerticesM()[vertex], full_reference.SourceVerticesM()[vertex])
				|| !EqualNormalized(full_reference.SourceVertexVelocitiesMPerS()[vertex], {{0.0, 0.0, 0.0}}))
				throw std::runtime_error("material patch map requires an undeformed, zero-velocity full reference state");
		// The reference map remains complete on each partition. Publication
		// ownership is a separate slice and is certified by the MPI caller.
		if (patch_interface.boundary_labels.size() != 1
			|| patch_interface.boundary_labels[0] < 0
			|| static_cast<std::uint64_t>(patch_interface.boundary_labels[0]) != patch_label)
			throw std::runtime_error("material patch map requires exactly its one selected label");
		if (layout_triangle_to_source_triangle.size() != patch_layout.reference_triangles.size())
			throw std::runtime_error("material patch map triangle mapping does not align with layout");

		std::set<std::uint32_t> selected_source_triangles;
		for (std::size_t index = 0; index < full_reference.SourceTriangles().size(); ++index)
			if (full_reference.SourceTriangles()[index].boundary_id == patch_label)
				selected_source_triangles.insert(static_cast<std::uint32_t>(index));
		if (selected_source_triangles.empty()) throw std::runtime_error("material patch label selects no source triangles");
		std::set<std::uint32_t> selected_source_vertices;
		for (const auto source_triangle : selected_source_triangles)
			for (const auto source_vertex : full_reference.SourceTriangles()[source_triangle].source_vertex_indices)
				selected_source_vertices.insert(source_vertex);
		if (global_to_source_vertex.size() != selected_source_vertices.size()
			|| patch_layout.reference_positions.size() != global_to_source_vertex.size())
			throw std::runtime_error("material patch map has missing or extra patch vertices");
		for (const auto& mapping : global_to_source_vertex) {
			if (mapping.second >= full_reference.ReferenceMaterialVerticesM().size())
				throw std::runtime_error("material patch map vertex mapping is invalid");
			const auto position = FindPosition(patch_layout, mapping.first);
			if (!EqualNormalized(position, full_reference.ReferenceMaterialVerticesM()[mapping.second]))
				throw std::runtime_error("material patch layout reference coordinate does not exactly match material reference");
		}
		for (std::size_t layout_triangle = 0; layout_triangle < patch_layout.reference_triangles.size(); ++layout_triangle) {
			const auto source_triangle = layout_triangle_to_source_triangle[layout_triangle];
			if (source_triangle >= full_reference.SourceTriangles().size())
				throw std::runtime_error("material patch source triangle is out of range");
			std::array<std::uint32_t, 3> mapped{};
			for (std::size_t corner = 0; corner < 3; ++corner) mapped[corner] = FindSource(global_to_source_vertex,
				patch_layout.reference_triangles[layout_triangle][corner]);
			if (!EvenPermutation(mapped, full_reference.SourceTriangles()[source_triangle].source_vertex_indices))
				throw std::runtime_error("material patch triangle winding is reversed or nonconforming");
		}

		const auto inferred_boundary = InferBoundarySourceVertices(full_reference, selected_source_triangles);
		if (inferred_boundary.empty()) throw std::runtime_error("material patch requires a nonempty one-sided boundary");
		std::vector<std::uint64_t> inferred_clamps;
		for (const auto& mapping : global_to_source_vertex)
			if (inferred_boundary.count(mapping.second)) inferred_clamps.push_back(mapping.first);
		if (inferred_clamps != configured_clamped_global_node_ids)
			throw std::runtime_error("material patch inferred and configured clamps differ");
		for (const auto& mapping : global_to_source_vertex)
			if (VertexTouchesNonpatch(full_reference, mapping.second, selected_source_triangles)
				&& !std::binary_search(configured_clamped_global_node_ids.begin(), configured_clamped_global_node_ids.end(), mapping.first))
				throw std::runtime_error("material patch seam vertex is not clamped");
		if (!ConnectedByEdges(patch_layout.reference_triangles))
			throw std::runtime_error("material patch must be connected through edges");

		const auto reference_identity = HashPatchReference(full_reference, patch_label, global_to_source_vertex,
			patch_layout.reference_triangles, layout_triangle_to_source_triangle);
		if (patch_layout.reference_mesh_identity_sha256 != reference_identity
			|| patch_interface.reference_mesh_identity_sha256 != reference_identity)
			throw std::runtime_error("material patch layout/interface reference identity is not the derived patch identity");
		MaterialSurfacePatchMap result(std::move(patch_interface), std::move(patch_layout), std::move(full_reference),
			patch_label, std::move(global_to_source_vertex), std::move(layout_triangle_to_source_triangle),
			std::move(configured_clamped_global_node_ids), reference_identity);
		result.identity_sha256_ = result.HashIdentity();
		return result;
	}

	const DistributedSurfaceInterface& Interface() const noexcept { return interface_; }
	const DistributedSurfaceLayout& Layout() const noexcept { return layout_; }
	const MaterialSurfaceKinematics& FullReference() const noexcept { return full_reference_; }
	std::uint32_t PatchLabel() const noexcept { return patch_label_; }
	const std::vector<GlobalToSourceVertex>& GlobalToSourceVertices() const noexcept { return global_to_source_vertex_; }
	const std::vector<std::uint32_t>& LayoutTriangleToSourceTriangles() const noexcept { return layout_triangle_to_source_triangle_; }
	const std::vector<std::uint64_t>& ConfiguredClampedGlobalNodeIds() const noexcept { return configured_clamped_global_node_ids_; }
	const std::string& ReferenceIdentitySha256() const noexcept { return reference_identity_sha256_; }
	const std::string& IdentitySha256() const noexcept { return identity_sha256_; }
	std::uint32_t SourceVertexForGlobalNode(std::uint64_t id) const { return FindSource(global_to_source_vertex_, id); }
	std::uint32_t SourceTriangleForLayoutTriangle(std::size_t index) const
	{ if (index >= layout_triangle_to_source_triangle_.size()) throw std::out_of_range("patch layout triangle is out of range"); return layout_triangle_to_source_triangle_[index]; }
	std::uint32_t CanonicalVertexForGlobalNode(std::uint64_t id) const { return CanonicalVertexForSource(SourceVertexForGlobalNode(id)); }
	std::uint32_t CanonicalTriangleForLayoutTriangle(std::size_t index) const { return CanonicalTriangleForSource(SourceTriangleForLayoutTriangle(index)); }


	// Canonical sorting depends on current coordinates; only source triangle IDs
	// remain stable across material motion.
	std::vector<std::size_t> LayoutTrianglesByCanonical(const MaterialSurfaceKinematics& current) const
	{
		if (current.MaterialIdentitySha256() != full_reference_.MaterialIdentitySha256()
			|| current.TopologyIdentitySha256() != full_reference_.TopologyIdentitySha256())
			throw std::invalid_argument("current material does not match patch reference topology");
		const auto absent = std::numeric_limits<std::size_t>::max();
		std::vector<std::size_t> layout_for_source(current.SourceTriangles().size(), absent);
		for (std::size_t i = 0; i < layout_triangle_to_source_triangle_.size(); ++i)
			layout_for_source.at(layout_triangle_to_source_triangle_[i]) = i;
		std::vector<std::size_t> result;
		result.reserve(current.CanonicalTriangleProvenance().size());
		for (const auto& triangle : current.CanonicalTriangleProvenance())
			result.push_back(layout_for_source.at(triangle.source_triangle));
		return result;
	}

private:
	MaterialSurfacePatchMap(DistributedSurfaceInterface interface, DistributedSurfaceLayout layout,
		MaterialSurfaceKinematics full_reference, std::uint32_t patch_label,
		std::vector<GlobalToSourceVertex> vertices, std::vector<std::uint32_t> triangles,
		std::vector<std::uint64_t> clamps, std::string reference_identity)
		: interface_(std::move(interface)), layout_(std::move(layout)), full_reference_(std::move(full_reference)),
		  patch_label_(patch_label), global_to_source_vertex_(std::move(vertices)),
		  layout_triangle_to_source_triangle_(std::move(triangles)), configured_clamped_global_node_ids_(std::move(clamps)),
		  reference_identity_sha256_(std::move(reference_identity)) {}
	static bool EqualNormalized(const std::array<double, 3>& a, const std::array<double, 3>& b)
	{ for (std::size_t i = 0; i < 3; ++i) if ((a[i] == 0.0 ? 0.0 : a[i]) != (b[i] == 0.0 ? 0.0 : b[i])) return false; return true; }
	static const std::array<double, 3>& FindPosition(const DistributedSurfaceLayout& layout, std::uint64_t id)
	{ const auto it = std::lower_bound(layout.reference_positions.begin(), layout.reference_positions.end(), id, [](const SurfaceReferencePosition& p, std::uint64_t value) { return p.global_node_id < value; }); if (it == layout.reference_positions.end() || it->global_node_id != id) throw std::runtime_error("patch global node is absent from layout positions"); return it->position_m; }
	static std::uint32_t FindSource(const std::vector<GlobalToSourceVertex>& mapping, std::uint64_t id)
	{ const auto it = std::lower_bound(mapping.begin(), mapping.end(), id, [](const GlobalToSourceVertex& p, std::uint64_t value) { return p.first < value; }); if (it == mapping.end() || it->first != id) throw std::runtime_error("patch global node has no source vertex mapping"); return it->second; }
	static bool EvenPermutation(const std::array<std::uint32_t, 3>& candidate, const std::array<std::uint32_t, 3>& source)
	{ std::array<std::uint32_t, 3> p{}; for (std::size_t i = 0; i < 3; ++i) { const auto it = std::find(source.begin(), source.end(), candidate[i]); if (it == source.end()) return false; p[i] = static_cast<std::uint32_t>(it-source.begin()); } return p[0] != p[1] && p[0] != p[2] && p[1] != p[2] && ((p[0] > p[1]) + (p[0] > p[2]) + (p[1] > p[2]))%2 == 0; }
	static std::set<std::uint32_t> InferBoundarySourceVertices(const MaterialSurfaceKinematics& full, const std::set<std::uint32_t>& selected)
	{ std::map<std::array<std::uint32_t, 2>, unsigned> edges; for (const auto triangle_index : selected) { const auto& t = full.SourceTriangles()[triangle_index].source_vertex_indices; for (std::size_t i = 0; i < 3; ++i) { std::array<std::uint32_t, 2> edge{{t[i], t[(i+1)%3]}}; std::sort(edge.begin(), edge.end()); ++edges[edge]; } } std::set<std::uint32_t> result; for (const auto& e : edges) { if (e.second > 2) throw std::runtime_error("material patch is nonmanifold"); if (e.second == 1) { result.insert(e.first[0]); result.insert(e.first[1]); } } return result; }
	static bool VertexTouchesNonpatch(const MaterialSurfaceKinematics& full, std::uint32_t vertex, const std::set<std::uint32_t>& selected)
	{ for (std::size_t i = 0; i < full.SourceTriangles().size(); ++i) if (!selected.count(static_cast<std::uint32_t>(i))) for (const auto v : full.SourceTriangles()[i].source_vertex_indices) if (v == vertex) return true; return false; }
	static bool ConnectedByEdges(const std::vector<std::array<std::uint64_t, 3>>& triangles)
	{ std::vector<bool> visited(triangles.size()); std::vector<std::size_t> queue{0}; visited[0] = true; for (std::size_t cursor = 0; cursor < queue.size(); ++cursor) for (std::size_t other = 0; other < triangles.size(); ++other) if (!visited[other]) { unsigned shared = 0; for (const auto a : triangles[queue[cursor]]) for (const auto b : triangles[other]) if (a == b) ++shared; if (shared >= 2) { visited[other] = true; queue.push_back(other); } } return std::find(visited.begin(), visited.end(), false) == visited.end(); }
	static void ValidateReferenceInputs(const MaterialSurfaceKinematics& full, std::uint32_t label,
		const std::vector<GlobalToSourceVertex>& vertices,
		const std::vector<std::array<std::uint64_t, 3>>& layout_triangles,
		const std::vector<std::uint32_t>& triangles,
		const std::vector<std::uint64_t>& clamps)
	{
		full.Validate();
		if (!std::is_sorted(vertices.begin(), vertices.end(), [](const GlobalToSourceVertex& left, const GlobalToSourceVertex& right) { return left.first < right.first; })
			|| std::adjacent_find(vertices.begin(), vertices.end(), [](const GlobalToSourceVertex& left, const GlobalToSourceVertex& right) { return left.first == right.first; }) != vertices.end())
			throw std::runtime_error("material patch map global/source vertex IDs must be sorted and unique");
		if (!std::is_sorted(clamps.begin(), clamps.end()) || std::adjacent_find(clamps.begin(), clamps.end()) != clamps.end())
			throw std::runtime_error("material patch map configured clamps must be sorted and unique");
		std::set<std::uint32_t> selected_triangles;
		for (std::size_t index = 0; index < full.SourceTriangles().size(); ++index)
			if (full.SourceTriangles()[index].boundary_id == label) selected_triangles.insert(static_cast<std::uint32_t>(index));
		if (selected_triangles.empty()) throw std::runtime_error("material patch label selects no source triangles");
		if (triangles.size() != selected_triangles.size() || layout_triangles.size() != triangles.size()) throw std::runtime_error("material patch triangles must be the exact selected source-label subset");
		std::set<std::uint32_t> mapped_triangles;
		for (const auto triangle : triangles) {
			if (triangle >= full.SourceTriangles().size()) throw std::runtime_error("material patch source triangle is out of range");
			mapped_triangles.insert(triangle);
		}
		if (mapped_triangles.size() != triangles.size() || mapped_triangles != selected_triangles)
			throw std::runtime_error("material patch triangles must be the exact selected source-label subset");
		std::set<std::uint32_t> selected_vertices;
		for (const auto triangle : selected_triangles)
			for (const auto vertex : full.SourceTriangles()[triangle].source_vertex_indices) selected_vertices.insert(vertex);
		if (vertices.size() != selected_vertices.size()) throw std::runtime_error("material patch map has missing or extra patch vertices");
		std::set<std::uint32_t> mapped_vertices;
		for (const auto& vertex : vertices) {
			if (vertex.second >= full.ReferenceMaterialVerticesM().size() || !mapped_vertices.insert(vertex.second).second)
				throw std::runtime_error("material patch map vertex mapping is invalid");
		}
		if (mapped_vertices != selected_vertices) throw std::runtime_error("material patch map vertices are not the exact selected subset");
		for (std::size_t layout_triangle = 0; layout_triangle < layout_triangles.size(); ++layout_triangle) {
			std::array<std::uint32_t, 3> mapped{};
			for (std::size_t corner = 0; corner < 3; ++corner) mapped[corner] = FindSource(vertices, layout_triangles[layout_triangle][corner]);
			if (!EvenPermutation(mapped, full.SourceTriangles()[triangles[layout_triangle]].source_vertex_indices))
				throw std::runtime_error("material patch triangle winding is reversed or nonconforming");
		}
		for (const auto clamp : clamps) {
			const auto mapped = std::lower_bound(vertices.begin(), vertices.end(), clamp,
				[](const GlobalToSourceVertex& vertex, std::uint64_t id) { return vertex.first < id; });
			if (mapped == vertices.end() || mapped->first != clamp)
				throw std::runtime_error("material patch clamp has no mapped vertex");
		}
	}
	static std::string HashPatchReference(const MaterialSurfaceKinematics& full, std::uint32_t label,
		const std::vector<GlobalToSourceVertex>& vertices,
		const std::vector<std::array<std::uint64_t, 3>>& layout_triangles,
		const std::vector<std::uint32_t>& triangles)
	{
		Sha256 hash; distributed_surface_detail::AppendString(hash, "MaterialSurfacePatchMap/reference/v2"); hash.AppendLittleEndian32(label);
		hash.AppendLittleEndian64(vertices.size());
		for (const auto& vertex : vertices) { hash.AppendLittleEndian64(vertex.first); for (const double x : full.ReferenceMaterialVerticesM()[vertex.second]) hash.AppendNormalizedDouble(x); }
		hash.AppendLittleEndian64(triangles.size());
		for (std::size_t index = 0; index < triangles.size(); ++index) {
			const auto triangle = triangles[index];
			const auto& source = full.SourceTriangles()[triangle]; hash.AppendLittleEndian32(source.boundary_id);
			for (const auto global : layout_triangles[index]) hash.AppendLittleEndian64(global);
		}
		return hash.Hex();
	}
	std::uint32_t CanonicalVertexForSource(std::uint32_t source) const
	{ for (std::size_t c = 0; c < full_reference_.CanonicalTriangleProvenance().size(); ++c) { const auto& p = full_reference_.CanonicalTriangleProvenance()[c]; const auto& t = full_reference_.Surface().Triangles()[c]; for (std::size_t corner = 0; corner < 3; ++corner) if (p.source_vertex_indices[p.canonical_corner_to_source_corner[corner]] == source) return t.indices[corner]; } throw std::runtime_error("full source vertex has no canonical vertex"); }
	std::uint32_t CanonicalTriangleForSource(std::uint32_t source) const
	{ for (std::size_t c = 0; c < full_reference_.CanonicalTriangleProvenance().size(); ++c) if (full_reference_.CanonicalTriangleProvenance()[c].source_triangle == source) return static_cast<std::uint32_t>(c); throw std::runtime_error("full source triangle has no canonical triangle"); }
	std::string HashIdentity() const
	{ Sha256 hash; distributed_surface_detail::AppendString(hash, "MaterialSurfacePatchMap/v2"); distributed_surface_detail::AppendString(hash, reference_identity_sha256_); distributed_surface_detail::AppendString(hash, layout_.layout_identity_sha256); distributed_surface_detail::AppendString(hash, BuildDistributedSurfacePartitionIdentitySha256(layout_)); distributed_surface_detail::AppendString(hash, full_reference_.MaterialIdentitySha256()); distributed_surface_detail::AppendString(hash, full_reference_.TopologyIdentitySha256()); distributed_surface_detail::AppendString(hash, BuildDistributedSurfaceInterfaceIdentitySha256(interface_)); hash.AppendLittleEndian64(global_to_source_vertex_.size()); for (const auto& vertex : global_to_source_vertex_) { hash.AppendLittleEndian64(vertex.first); hash.AppendLittleEndian32(vertex.second); } hash.AppendLittleEndian64(layout_triangle_to_source_triangle_.size()); for (const auto triangle : layout_triangle_to_source_triangle_) hash.AppendLittleEndian32(triangle); hash.AppendLittleEndian64(configured_clamped_global_node_ids_.size()); for (const auto clamp : configured_clamped_global_node_ids_) hash.AppendLittleEndian64(clamp); return hash.Hex(); }
	DistributedSurfaceInterface interface_; DistributedSurfaceLayout layout_; MaterialSurfaceKinematics full_reference_;
	std::uint32_t patch_label_ = 0; std::vector<GlobalToSourceVertex> global_to_source_vertex_; std::vector<std::uint32_t> layout_triangle_to_source_triangle_; std::vector<std::uint64_t> configured_clamped_global_node_ids_; std::string reference_identity_sha256_, identity_sha256_;
};

} // namespace iga

#endif

#ifndef IGA_MATERIAL_SURFACE_KINEMATICS_HPP
#define IGA_MATERIAL_SURFACE_KINEMATICS_HPP

// Producer-neutral immutable material-surface state for one exact geometry
// epoch.  It deliberately owns no distributed field layout or FSI graph
// endpoint: those contracts are separate from the geometry consumed by the
// immersed cut, wall, and transient-runtime paths.
#include "SurfaceGeometry.hpp"
#include "Sha256.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

class PrescribedSurfaceMotion;

struct SourceTriangleProvenance {
	std::uint32_t source_triangle = 0;
	std::array<std::uint32_t, 3> source_vertex_indices{{0, 0, 0}};
	std::uint32_t boundary_id = 0;
	// For canonical-triangle provenance only, maps a canonical triangle corner
	// to the corresponding source/material triangle corner.
	std::array<std::uint32_t, 3> canonical_corner_to_source_corner{{0, 1, 2}};
	bool operator==(const SourceTriangleProvenance& other) const
	{
		return source_triangle == other.source_triangle && source_vertex_indices == other.source_vertex_indices
			&& boundary_id == other.boundary_id
			&& canonical_corner_to_source_corner == other.canonical_corner_to_source_corner;
	}
};

class MaterialSurfaceKinematics {
public:
	MaterialSurfaceKinematics(const MaterialSurfaceKinematics&) = default;
	MaterialSurfaceKinematics(MaterialSurfaceKinematics&&) noexcept = default;
	MaterialSurfaceKinematics& operator=(const MaterialSurfaceKinematics&) = delete;
	MaterialSurfaceKinematics& operator=(MaterialSurfaceKinematics&&) = delete;

	// A structural producer may construct the same neutral payload after it has
	// validated/canonicalized its own material surface.  Every supplied identity
	// is checked against a domain-separated digest of owned fields: generic
	// producers use the content identity as their epoch identity, while the
	// prescribed adapter may retain its separately verified historical epoch
	// digest for compatibility.
	static MaterialSurfaceKinematics Create(ClosedTriangulatedSurface surface,
		std::vector<std::array<double, 3>> reference_material_vertices_m,
		std::vector<std::array<double, 3>> source_vertices_m,
		std::vector<std::array<double, 3>> source_vertex_velocities_m_per_s,
		std::vector<SourceTriangleProvenance> canonical_triangle_provenance,
		std::vector<SourceTriangleProvenance> source_triangles, double evaluated_time_s,
		double step_start_s, double step_end_s, std::string geometry_epoch_identity_sha256,
		std::string material_identity_sha256, std::string topology_identity_sha256)
	{
		MaterialSurfaceKinematics result(std::move(surface));
		result.reference_material_vertices_m_ = std::move(reference_material_vertices_m);
		result.source_vertices_m_ = std::move(source_vertices_m);
		result.source_vertex_velocities_m_per_s_ = std::move(source_vertex_velocities_m_per_s);
		result.canonical_triangle_provenance_ = std::move(canonical_triangle_provenance);
		result.source_triangles_ = std::move(source_triangles);
		result.evaluated_time_s_ = evaluated_time_s;
		result.step_start_s_ = step_start_s;
		result.step_end_s_ = step_end_s;
		result.identity_sha256_ = std::move(geometry_epoch_identity_sha256);
		result.material_identity_sha256_ = std::move(material_identity_sha256);
		result.topology_identity_sha256_ = std::move(topology_identity_sha256);
		result.content_identity_sha256_ = result.HashContentIdentity();
		result.Validate();
		return result;
	}

	// Preferred producer-neutral route.  A future membrane producer need not
	// know any hashing format: all identities are derived after the complete
	// payload has been moved into this immutable owner.
	static MaterialSurfaceKinematics Create(ClosedTriangulatedSurface surface,
		std::vector<std::array<double, 3>> reference_material_vertices_m,
		std::vector<std::array<double, 3>> source_vertices_m,
		std::vector<std::array<double, 3>> source_vertex_velocities_m_per_s,
		std::vector<SourceTriangleProvenance> canonical_triangle_provenance,
		std::vector<SourceTriangleProvenance> source_triangles, double evaluated_time_s,
		double step_start_s, double step_end_s)
	{
		MaterialSurfaceKinematics result(std::move(surface));
		result.reference_material_vertices_m_ = std::move(reference_material_vertices_m);
		result.source_vertices_m_ = std::move(source_vertices_m);
		result.source_vertex_velocities_m_per_s_ = std::move(source_vertex_velocities_m_per_s);
		result.canonical_triangle_provenance_ = std::move(canonical_triangle_provenance);
		result.source_triangles_ = std::move(source_triangles);
		result.evaluated_time_s_ = evaluated_time_s;
		result.step_start_s_ = step_start_s;
		result.step_end_s_ = step_end_s;
		result.material_identity_sha256_ = result.HashMaterialIdentity();
		result.topology_identity_sha256_ = result.HashTopologyIdentity();
		result.content_identity_sha256_ = result.HashContentIdentity();
		result.identity_sha256_ = result.content_identity_sha256_;
		result.Validate();
		return result;
	}

	const ClosedTriangulatedSurface& Surface() const noexcept { return surface_; }
	// Immutable material/reference coordinates, distinct from the evaluated
	// current SourceVerticesM() coordinates.
	const std::vector<std::array<double, 3>>& ReferenceMaterialVerticesM() const noexcept
	{ return reference_material_vertices_m_; }
	const std::vector<std::array<double, 3>>& SourceVerticesM() const noexcept { return source_vertices_m_; }
	const std::vector<std::array<double, 3>>& SourceVertexVelocitiesMPerS() const noexcept
	{ return source_vertex_velocities_m_per_s_; }
	const std::vector<SourceTriangleProvenance>& CanonicalTriangleProvenance() const noexcept
	{ return canonical_triangle_provenance_; }
	const std::vector<SourceTriangleProvenance>& SourceTriangles() const noexcept { return source_triangles_; }
	const std::string& IdentitySha256() const noexcept { return identity_sha256_; }
	// The epoch identity is intentionally the historical evaluated-state hash.
	// This keeps existing geometry identities bit-for-bit stable while making
	// the epoch role explicit to producer-neutral consumers.
	const std::string& GeometryEpochIdentitySha256() const noexcept { return identity_sha256_; }
	// Complete current payload identity.  Unlike the legacy prescribed epoch
	// digest, this includes the step envelope and is always canonical here.
	const std::string& ContentIdentitySha256() const noexcept { return content_identity_sha256_; }
	const std::string& MaterialIdentitySha256() const noexcept { return material_identity_sha256_; }
	const std::string& TopologyIdentitySha256() const noexcept { return topology_identity_sha256_; }
	double EvaluatedTimeS() const noexcept { return evaluated_time_s_; }
	double StepStartS() const noexcept { return step_start_s_; }
	double StepEndS() const noexcept { return step_end_s_; }
	double DtS() const noexcept { return step_end_s_-step_start_s_; }

	void Validate() const
	{
		if (!std::isfinite(evaluated_time_s_) || !std::isfinite(step_start_s_) || !std::isfinite(step_end_s_)
			|| !(step_end_s_ > step_start_s_) || evaluated_time_s_ < step_start_s_
			|| evaluated_time_s_ > step_end_s_ || identity_sha256_.empty()
			|| content_identity_sha256_.empty() || material_identity_sha256_.empty() || topology_identity_sha256_.empty())
			throw std::invalid_argument("material surface kinematics time or identity is invalid");
		if (reference_material_vertices_m_.empty()
			|| reference_material_vertices_m_.size() != source_vertices_m_.size()
			|| source_vertices_m_.size() != surface_.Vertices().size()
			|| source_vertices_m_.size() != source_vertex_velocities_m_per_s_.size()
			|| canonical_triangle_provenance_.size() != surface_.Triangles().size() || source_triangles_.empty())
			throw std::invalid_argument("material surface kinematics geometry or provenance is invalid");
		for (const auto& vertex : reference_material_vertices_m_)
			for (const double value : vertex)
				if (!std::isfinite(value)) throw std::invalid_argument("material surface kinematics reference vertex is nonfinite");
		for (const auto& vertex : source_vertices_m_)
			for (const double value : vertex)
				if (!std::isfinite(value)) throw std::invalid_argument("material surface kinematics vertex is nonfinite");
		for (const auto& velocity : source_vertex_velocities_m_per_s_)
			for (const double value : velocity)
				if (!std::isfinite(value)) throw std::invalid_argument("material surface kinematics velocity is nonfinite");
		std::vector<bool> source_triangle_used(source_triangles_.size(), false);
		std::vector<bool> source_vertex_used(source_vertices_m_.size(), false);
		const auto unmapped_vertex = std::numeric_limits<std::uint32_t>::max();
		std::vector<std::uint32_t> canonical_vertex_for_source(source_vertices_m_.size(), unmapped_vertex);
		std::vector<std::uint32_t> source_vertex_for_canonical(surface_.Vertices().size(), unmapped_vertex);
		for (std::size_t source_index = 0; source_index < source_triangles_.size(); ++source_index) {
			const auto& source = source_triangles_[source_index];
			if (source.source_triangle != source_index
				|| source.canonical_corner_to_source_corner != std::array<std::uint32_t, 3>{{0, 1, 2}})
				throw std::invalid_argument("material surface kinematics source triangle index or orientation is inconsistent");
			for (std::size_t corner = 0; corner < 3; ++corner) {
				const auto source_vertex = source.source_vertex_indices[corner];
				if (source_vertex >= source_vertices_m_.size())
					throw std::invalid_argument("material surface kinematics source vertex is out of range");
				if (std::find(source.source_vertex_indices.begin(), source.source_vertex_indices.begin()+corner,
					source_vertex) != source.source_vertex_indices.begin()+corner)
					throw std::invalid_argument("material surface kinematics source triangle is degenerate");
				source_vertex_used[source_vertex] = true;
			}
		}
		for (std::size_t canonical_triangle = 0; canonical_triangle < canonical_triangle_provenance_.size(); ++canonical_triangle) {
			const auto& provenance = canonical_triangle_provenance_[canonical_triangle];
			if (provenance.source_triangle >= source_triangles_.size())
				throw std::invalid_argument("material surface kinematics triangle provenance is out of range");
			const auto& source = source_triangles_[provenance.source_triangle];
			if (source.source_triangle != provenance.source_triangle
				|| source.source_vertex_indices != provenance.source_vertex_indices
				|| source.boundary_id != provenance.boundary_id)
				throw std::invalid_argument("material surface kinematics triangle provenance is inconsistent");
			if (provenance.boundary_id != surface_.Triangles()[canonical_triangle].boundary_id)
				throw std::invalid_argument("material surface kinematics canonical/source boundary label is inconsistent");
			if (source_triangle_used[provenance.source_triangle])
				throw std::invalid_argument("material surface kinematics source triangle mapping is not bijective");
			source_triangle_used[provenance.source_triangle] = true;
			std::array<bool, 3> source_corner_used{{false, false, false}};
			for (const auto source_vertex : provenance.source_vertex_indices)
				if (source_vertex >= source_vertices_m_.size())
					throw std::invalid_argument("material surface kinematics source vertex is out of range");
			for (const auto source_corner : provenance.canonical_corner_to_source_corner) {
				if (source_corner >= 3) throw std::invalid_argument("material surface kinematics corner provenance is out of range");
				if (source_corner_used[source_corner])
					throw std::invalid_argument("material surface kinematics corner provenance is not a permutation");
				source_corner_used[source_corner] = true;
			}
			const auto& permutation = provenance.canonical_corner_to_source_corner;
			const unsigned inversions = (permutation[0] > permutation[1]) + (permutation[0] > permutation[2])
				+ (permutation[1] > permutation[2]);
			if (inversions%2 != 0)
				throw std::invalid_argument("material surface kinematics triangle orientation is reversed");
		}
		if (std::find(source_triangle_used.begin(), source_triangle_used.end(), false) != source_triangle_used.end()
			|| std::find(source_vertex_used.begin(), source_vertex_used.end(), false) != source_vertex_used.end())
			throw std::invalid_argument("material surface kinematics has unreferenced source topology");
		for (std::size_t canonical_triangle = 0; canonical_triangle < surface_.Triangles().size(); ++canonical_triangle) {
			const auto& triangle = surface_.Triangles()[canonical_triangle];
			const auto& provenance = canonical_triangle_provenance_[canonical_triangle];
			for (std::size_t canonical_corner = 0; canonical_corner < 3; ++canonical_corner) {
				const auto source_corner = provenance.canonical_corner_to_source_corner[canonical_corner];
				const auto source_vertex = provenance.source_vertex_indices[source_corner];
				const auto canonical_vertex = triangle.indices[canonical_corner];
				if (source_vertices_m_[source_vertex] != surface_.Vertices()[canonical_vertex])
					throw std::invalid_argument("material surface kinematics canonical/material geometry is inconsistent");
				if ((canonical_vertex_for_source[source_vertex] != unmapped_vertex
					&& canonical_vertex_for_source[source_vertex] != canonical_vertex)
					|| (source_vertex_for_canonical[canonical_vertex] != unmapped_vertex
						&& source_vertex_for_canonical[canonical_vertex] != source_vertex))
					throw std::invalid_argument("material surface kinematics vertex mapping is not one-to-one");
				canonical_vertex_for_source[source_vertex] = canonical_vertex;
				source_vertex_for_canonical[canonical_vertex] = source_vertex;
			}
		}
		if (std::find(canonical_vertex_for_source.begin(), canonical_vertex_for_source.end(), unmapped_vertex)
			!= canonical_vertex_for_source.end()
			|| std::find(source_vertex_for_canonical.begin(), source_vertex_for_canonical.end(), unmapped_vertex)
				!= source_vertex_for_canonical.end())
			throw std::invalid_argument("material surface kinematics vertex mapping is incomplete");
		if (material_identity_sha256_ != HashMaterialIdentity()
			|| topology_identity_sha256_ != HashTopologyIdentity()
			|| content_identity_sha256_ != HashContentIdentity()
			|| (identity_sha256_ != content_identity_sha256_
				&& identity_sha256_ != HashLegacyPrescribedEpochIdentity()))
			throw std::invalid_argument("material surface kinematics identity does not match its owned payload");
	}

	void ValidateNextEpoch(const MaterialSurfaceKinematics& previous, double expected_time_s,
		double expected_dt_s) const
	{
		Validate(); previous.Validate();
		if (!std::isfinite(expected_time_s) || !std::isfinite(expected_dt_s) || !(expected_dt_s > 0.0)
			|| evaluated_time_s_ != expected_time_s || step_start_s_ != previous.evaluated_time_s_
			|| step_end_s_ != expected_time_s || step_end_s_ != step_start_s_+expected_dt_s
			|| material_identity_sha256_ != previous.material_identity_sha256_
			|| topology_identity_sha256_ != previous.topology_identity_sha256_
			|| ContentIdentitySha256() == previous.ContentIdentitySha256())
			throw std::invalid_argument("material surface kinematics is stale or does not match the exact next epoch");
	}

	std::array<double, 3> WallVelocity(std::uint32_t canonical_triangle,
		const std::array<double, 3>& barycentric) const
	{
		if (canonical_triangle >= surface_.Triangles().size()
			|| canonical_triangle >= canonical_triangle_provenance_.size())
			throw std::out_of_range("canonical triangle is out of range");
		const auto& provenance = canonical_triangle_provenance_[canonical_triangle];
		if (provenance.source_triangle >= source_triangles_.size())
			throw std::out_of_range("canonical triangle source provenance is out of range");
		const auto& source_triangle = source_triangles_[provenance.source_triangle];
		if (source_triangle.source_triangle != provenance.source_triangle
			|| source_triangle.source_vertex_indices != provenance.source_vertex_indices
			|| source_triangle.boundary_id != provenance.boundary_id)
			throw std::invalid_argument("canonical triangle source provenance is inconsistent");
		constexpr double barycentric_tolerance = 64.0*std::numeric_limits<double>::epsilon();
		double weight_sum = 0.0;
		for (double weight : barycentric) {
			if (!std::isfinite(weight)) throw std::invalid_argument("barycentric weight is nonfinite");
			if (weight < -barycentric_tolerance || weight > 1.0+barycentric_tolerance)
				throw std::invalid_argument("barycentric weight is outside the simplex");
			weight_sum += weight;
			if (!std::isfinite(weight_sum)) throw std::invalid_argument("barycentric weight sum is nonfinite");
		}
		if (std::fabs(weight_sum-1.0) > barycentric_tolerance)
			throw std::invalid_argument("barycentric weights do not sum to one");
		std::array<double, 3> result{{0.0, 0.0, 0.0}};
		for (std::size_t canonical_corner = 0; canonical_corner < 3; ++canonical_corner) {
			const std::uint32_t source_corner = provenance.canonical_corner_to_source_corner[canonical_corner];
			if (source_corner >= 3) throw std::out_of_range("canonical triangle corner provenance is out of range");
			const std::uint32_t source_vertex = provenance.source_vertex_indices[source_corner];
			if (source_vertex >= source_vertex_velocities_m_per_s_.size())
				throw std::out_of_range("source vertex velocity is out of range");
			for (std::size_t axis = 0; axis < 3; ++axis)
				result[axis] += barycentric[canonical_corner]*source_vertex_velocities_m_per_s_[source_vertex][axis];
		}
		for (double value : result)
			if (!std::isfinite(value)) throw std::invalid_argument("wall velocity is nonfinite");
		return result;
	}

private:
	friend class PrescribedSurfaceMotion;
	explicit MaterialSurfaceKinematics(ClosedTriangulatedSurface checked) : surface_(std::move(checked)) {}
	double evaluated_time_s_ = 0.0, step_start_s_ = 0.0, step_end_s_ = 0.0;
	ClosedTriangulatedSurface surface_;
	std::vector<std::array<double, 3>> reference_material_vertices_m_;
	std::vector<std::array<double, 3>> source_vertices_m_;
	std::vector<std::array<double, 3>> source_vertex_velocities_m_per_s_;
	std::vector<SourceTriangleProvenance> canonical_triangle_provenance_;
	static void AppendCount(Sha256& hash, std::size_t count)
	{
		if (count > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max()))
			throw std::overflow_error("material surface kinematics digest count is unrepresentable");
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(count));
	}
	static void AppendString(Sha256& hash, const std::string& value)
	{
		AppendCount(hash, value.size()); hash.Append(value.data(), value.size());
	}
	void AppendSourceTopology(Sha256& hash) const
	{
		AppendCount(hash, source_triangles_.size());
		for (const auto& triangle : source_triangles_) {
			hash.AppendLittleEndian32(triangle.source_triangle); hash.AppendLittleEndian32(triangle.boundary_id);
			for (const auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
		}
	}
	void AppendTopology(Sha256& hash) const
	{
		AppendSourceTopology(hash);
		AppendCount(hash, surface_.Triangles().size());
		for (const auto& triangle : surface_.Triangles()) {
			for (const auto index : triangle.indices) hash.AppendLittleEndian32(index);
			hash.AppendLittleEndian32(triangle.boundary_id);
		}
		AppendCount(hash, canonical_triangle_provenance_.size());
		for (const auto& triangle : canonical_triangle_provenance_) {
			hash.AppendLittleEndian32(triangle.source_triangle); hash.AppendLittleEndian32(triangle.boundary_id);
			for (const auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
			for (const auto corner : triangle.canonical_corner_to_source_corner) hash.AppendLittleEndian32(corner);
		}
	}
	std::string HashTopologyIdentity() const
	{
		Sha256 hash; static constexpr char domain[] = "MaterialSurfaceKinematics/topology/v1";
		hash.Append(domain, sizeof(domain)-1); AppendSourceTopology(hash); return hash.Hex();
	}
	std::string HashMaterialIdentity() const
	{
		Sha256 hash; static constexpr char domain[] = "MaterialSurfaceKinematics/material/v1";
		hash.Append(domain, sizeof(domain)-1); AppendCount(hash, reference_material_vertices_m_.size());
		for (const auto& vertex : reference_material_vertices_m_)
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendString(hash, HashTopologyIdentity()); return hash.Hex();
	}
	std::string HashContentIdentity() const
	{
		Sha256 hash; static constexpr char domain[] = "MaterialSurfaceKinematics/content/v1";
		hash.Append(domain, sizeof(domain)-1);
		hash.AppendNormalizedDouble(evaluated_time_s_); hash.AppendNormalizedDouble(step_start_s_); hash.AppendNormalizedDouble(step_end_s_);
		AppendCount(hash, source_vertices_m_.size()); for (const auto& vertex : source_vertices_m_)
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendCount(hash, source_vertex_velocities_m_per_s_.size()); for (const auto& velocity : source_vertex_velocities_m_per_s_)
			for (const auto value : velocity) hash.AppendNormalizedDouble(value);
		AppendCount(hash, surface_.Vertices().size()); for (const auto& vertex : surface_.Vertices())
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendTopology(hash); AppendString(hash, surface_.CanonicalSha256()); return hash.Hex();
	}
	std::string HashLegacyPrescribedEpochIdentity() const
	{
		Sha256 hash; static constexpr char domain[] = "iga-prescribed-surface-motion-v2";
		hash.Append(domain, sizeof(domain)-1); hash.AppendNormalizedDouble(evaluated_time_s_);
		AppendCount(hash, source_vertices_m_.size()); for (const auto& vertex : source_vertices_m_)
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendCount(hash, source_vertex_velocities_m_per_s_.size()); for (const auto& velocity : source_vertex_velocities_m_per_s_)
			for (const auto value : velocity) hash.AppendNormalizedDouble(value);
		AppendCount(hash, source_triangles_.size()); for (const auto& triangle : source_triangles_) {
			hash.AppendLittleEndian32(triangle.source_triangle); hash.AppendLittleEndian32(triangle.boundary_id);
			for (const auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
		}
		AppendCount(hash, surface_.Vertices().size()); for (const auto& vertex : surface_.Vertices())
			for (const auto value : vertex) hash.AppendNormalizedDouble(value);
		AppendCount(hash, surface_.Triangles().size()); for (const auto& triangle : surface_.Triangles()) {
			for (const auto index : triangle.indices) hash.AppendLittleEndian32(index);
			hash.AppendLittleEndian32(triangle.boundary_id);
		}
		AppendString(hash, surface_.CanonicalSha256()); AppendCount(hash, canonical_triangle_provenance_.size());
		for (const auto& triangle : canonical_triangle_provenance_) {
			hash.AppendLittleEndian32(triangle.source_triangle); hash.AppendLittleEndian32(triangle.boundary_id);
			for (const auto index : triangle.source_vertex_indices) hash.AppendLittleEndian32(index);
			for (const auto corner : triangle.canonical_corner_to_source_corner) hash.AppendLittleEndian32(corner);
		}
		return hash.Hex();
	}
	std::string identity_sha256_, content_identity_sha256_, material_identity_sha256_, topology_identity_sha256_;
	std::vector<SourceTriangleProvenance> source_triangles_;
};

} // namespace iga

#endif

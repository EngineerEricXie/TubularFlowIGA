#ifndef IGA_DISTRIBUTED_SURFACE_INTERFACE_HPP
#define IGA_DISTRIBUTED_SURFACE_INTERFACE_HPP

// Dependency-free, ownership-oriented contracts for field-valued FSI surface
// exchange.  These are deliberately separate from CouplingPort, whose P/Q
// values are scalar boundary summaries.
#include "Sha256.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace iga {

enum class SurfaceFieldQuantity : std::uint8_t {
	Displacement,
	Velocity,
	TractionOnStructure
};

inline bool IsKnownSurfaceFieldQuantity(SurfaceFieldQuantity quantity)
{
	return quantity == SurfaceFieldQuantity::Displacement
		|| quantity == SurfaceFieldQuantity::Velocity
		|| quantity == SurfaceFieldQuantity::TractionOnStructure;
}

struct SurfaceInterfaceRef {
	std::string domain_id;
	// A domain may host more than one independently stepped solver.  Surface
	// endpoint identity therefore includes the concrete subsystem, not only a
	// user-facing interface label.
	std::string subsystem_id;
	std::string interface_id;
};

// A replicated reference position is keyed by immutable material node ID.
// Only owned_global_node_ids have publication authority on a rank; all other
// reference data support local geometry checks and may be used as scratch.
struct SurfaceReferencePosition {
	std::uint64_t global_node_id = 0;
	std::array<double, 3> position_m{};
};

struct DistributedSurfaceLayout {
	std::string reference_mesh_identity_sha256;
	std::string layout_identity_sha256;
	std::uint64_t global_node_count = 0;
	std::uint64_t partition_count = 0;
	std::uint64_t partition_rank = 0;
	std::vector<std::uint64_t> owned_global_node_ids;
	std::vector<SurfaceReferencePosition> reference_positions;
	std::vector<std::array<std::uint64_t, 3>> reference_triangles;
	// Aligned with owned_global_node_ids, in m^2.  These are the immutable
	// reference weights used for the first strong-coupling slice.
	std::vector<double> owned_reference_lumped_areas_m2;
};

struct SurfaceFieldStamp {
	double time_s = 0.0;
	std::uint64_t step = 0;
	std::uint64_t coupling_iteration = 0;
	std::string reference_mesh_identity_sha256;
	std::string layout_identity_sha256;
	// Binds the publication to this exact ownership slice, not merely to the
	// replicated global topology.  Equal-sized partitions are not fungible.
	std::string partition_identity_sha256;
	std::string producer_state_identity_sha256;
};

// A consumer can schedule the shape of an output publication before the
// producer has solved, but cannot know its state-derived identity.  This
// envelope deliberately binds every producer-neutral part of a stamp and is
// validated against the exact owned publication slice.
struct SurfaceFieldStampEnvelope {
	double time_s = 0.0;
	std::uint64_t step = 0;
	std::uint64_t coupling_iteration = 0;
	std::string reference_mesh_identity_sha256;
	std::string layout_identity_sha256;
	std::string partition_identity_sha256;
};

struct SurfaceKinematics {
	SurfaceInterfaceRef interface;
	SurfaceFieldStamp stamp;
	// Full Cartesian values, aligned with layout.owned_global_node_ids.
	std::vector<std::array<double, 3>> displacement_m;
	std::vector<std::array<double, 3>> velocity_m_per_s;
};

struct SurfaceTraction {
	SurfaceInterfaceRef interface;
	SurfaceFieldStamp stamp;
	// Fluid-on-structure Cauchy traction [Pa] and consistent nodal force [N],
	// both full Cartesian and aligned with owned_global_node_ids.
	std::vector<std::array<double, 3>> traction_on_structure_pa;
	std::vector<std::array<double, 3>> consistent_nodal_force_n;
	std::string projection_identity_sha256;
};

struct DistributedSurfaceInterface {
	SurfaceInterfaceRef id;
	std::string subsystem_id;
	std::string locator_kind = "material_surface";
	std::vector<std::int64_t> boundary_labels;
	std::string reference_mesh_identity_sha256;
	std::vector<SurfaceFieldQuantity> provides;
	std::vector<SurfaceFieldQuantity> requires;
};

namespace distributed_surface_detail {

inline void AppendString(Sha256& hash, const std::string& value)
{
	if (value.size() > std::numeric_limits<std::uint64_t>::max())
		throw std::runtime_error("surface identity string is too large");
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(value.size()));
	hash.Append(value.data(), value.size());
}

inline void AppendVector(Sha256& hash, const std::array<double, 3>& value)
{
	for (const double component : value) hash.AppendNormalizedDouble(component);
}

inline bool IsLowerHexSha256(const std::string& value)
{
	if (value.size() != 64) return false;
	for (const char character : value)
		if (!((character >= '0' && character <= '9')
			|| (character >= 'a' && character <= 'f'))) return false;
	return true;
}

inline void RequireFiniteVector(const std::array<double, 3>& value, const char* name)
{
	for (const double component : value)
		if (!std::isfinite(component)) throw std::runtime_error(std::string(name)+" must be finite");
}

inline void ValidateLayoutContents(const DistributedSurfaceLayout& layout)
{
	if (!IsLowerHexSha256(layout.reference_mesh_identity_sha256))
		throw std::runtime_error("surface layout reference mesh identity must be a lowercase SHA-256 hash");
	if (layout.global_node_count == 0) throw std::runtime_error("surface layout global node count must be positive");
	if (layout.partition_count == 0 || layout.partition_rank >= layout.partition_count)
		throw std::runtime_error("surface layout partition metadata is invalid");
	// Empty publication slices are valid on a multi-rank communicator. Global
	// coverage is checked collectively; a single rank must still own all nodes.
	if (!std::is_sorted(layout.owned_global_node_ids.begin(), layout.owned_global_node_ids.end())
		|| std::adjacent_find(layout.owned_global_node_ids.begin(), layout.owned_global_node_ids.end()) != layout.owned_global_node_ids.end())
		throw std::runtime_error("surface layout owned node IDs must be sorted and unique");
	if (layout.owned_global_node_ids.size() > layout.global_node_count)
		throw std::runtime_error("surface layout owned node count exceeds global node count");
	if (layout.partition_count == 1 && layout.owned_global_node_ids.size() != layout.global_node_count)
		throw std::runtime_error("single-rank surface layout must own every global node");
	if (layout.reference_positions.size() != layout.global_node_count)
		throw std::runtime_error("surface layout must provide every reference position");
	std::vector<std::uint64_t> position_ids;
	position_ids.reserve(layout.reference_positions.size());
	for (const auto& position : layout.reference_positions) {
		position_ids.push_back(position.global_node_id);
		RequireFiniteVector(position.position_m, "surface reference position");
	}
	if (!std::is_sorted(position_ids.begin(), position_ids.end())
		|| std::adjacent_find(position_ids.begin(), position_ids.end()) != position_ids.end())
		throw std::runtime_error("surface reference position node IDs must be sorted and unique");
	for (const auto id : layout.owned_global_node_ids)
		if (!std::binary_search(position_ids.begin(), position_ids.end(), id))
			throw std::runtime_error("surface layout owns an unknown reference node");
	if (layout.reference_triangles.empty()) throw std::runtime_error("surface layout reference triangles must be nonempty");
	for (const auto& triangle : layout.reference_triangles) {
		if (triangle[0] == triangle[1] || triangle[0] == triangle[2] || triangle[1] == triangle[2])
			throw std::runtime_error("surface layout triangle has duplicate nodes");
		for (const auto id : triangle)
			if (!std::binary_search(position_ids.begin(), position_ids.end(), id))
				throw std::runtime_error("surface layout triangle references an unknown node");
		const auto find_position = [&layout](std::uint64_t id) -> const std::array<double, 3>& {
			const auto iterator = std::lower_bound(layout.reference_positions.begin(),
				layout.reference_positions.end(), id, [](const SurfaceReferencePosition& position, std::uint64_t value) {
					return position.global_node_id < value;
				});
			return iterator->position_m;
		};
		const auto& first = find_position(triangle[0]);
		const auto& second = find_position(triangle[1]);
		const auto& third = find_position(triangle[2]);
		const std::array<double, 3> first_edge{{second[0]-first[0], second[1]-first[1], second[2]-first[2]}};
		const std::array<double, 3> second_edge{{third[0]-first[0], third[1]-first[1], third[2]-first[2]}};
		const std::array<double, 3> cross{{
			first_edge[1]*second_edge[2]-first_edge[2]*second_edge[1],
			first_edge[2]*second_edge[0]-first_edge[0]*second_edge[2],
			first_edge[0]*second_edge[1]-first_edge[1]*second_edge[0]}};
		const double squared_area_four = cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2];
		if (!std::isfinite(squared_area_four) || !(squared_area_four > 0.0))
			throw std::runtime_error("surface layout triangle is degenerate");
	}
	if (layout.owned_reference_lumped_areas_m2.size() != layout.owned_global_node_ids.size())
		throw std::runtime_error("surface layout lumped areas must align with owned node IDs");
	for (const double area : layout.owned_reference_lumped_areas_m2)
		if (!std::isfinite(area) || !(area > 0.0))
			throw std::runtime_error("surface layout reference lumped areas must be finite and positive");
}

} // namespace distributed_surface_detail

inline bool IsLowercaseSha256(const std::string& value)
{
	return distributed_surface_detail::IsLowerHexSha256(value);
}

inline void ValidateSurfaceInterfaceRef(const SurfaceInterfaceRef& reference)
{
	if (reference.domain_id.empty() || reference.subsystem_id.empty()
		|| reference.interface_id.empty())
		throw std::runtime_error("surface interface domain, subsystem, and interface IDs must be nonempty");
}

inline bool operator==(const SurfaceInterfaceRef& first,
	const SurfaceInterfaceRef& second)
{
	return first.domain_id == second.domain_id
		&& first.subsystem_id == second.subsystem_id
		&& first.interface_id == second.interface_id;
}

inline bool operator<(const SurfaceInterfaceRef& first,
	const SurfaceInterfaceRef& second)
{
	return std::tie(first.domain_id, first.subsystem_id, first.interface_id)
		< std::tie(second.domain_id, second.subsystem_id, second.interface_id);
}

inline std::string BuildDistributedSurfaceLayoutIdentitySha256(const DistributedSurfaceLayout& layout)
{
	distributed_surface_detail::ValidateLayoutContents(layout);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "DistributedSurfaceLayout/v1");
	distributed_surface_detail::AppendString(hash, layout.reference_mesh_identity_sha256);
	hash.AppendLittleEndian64(layout.global_node_count);
	hash.AppendLittleEndian64(layout.partition_count);
	// This is deliberately global: every rank can independently derive the
	// same identity from replicated material topology, while ownership remains
	// rank-local and is validated separately.
	for (const auto& position : layout.reference_positions) {
		hash.AppendLittleEndian64(position.global_node_id);
		distributed_surface_detail::AppendVector(hash, position.position_m);
	}
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(layout.reference_triangles.size()));
	for (const auto& triangle : layout.reference_triangles)
		for (const auto id : triangle) hash.AppendLittleEndian64(id);
	return hash.Hex();
}

inline void ValidateDistributedSurfaceLayout(const DistributedSurfaceLayout& layout)
{
	distributed_surface_detail::ValidateLayoutContents(layout);
	if (!IsLowercaseSha256(layout.layout_identity_sha256))
		throw std::runtime_error("surface layout identity must be a lowercase SHA-256 hash");
	if (layout.layout_identity_sha256 != BuildDistributedSurfaceLayoutIdentitySha256(layout))
		throw std::runtime_error("surface layout identity does not match layout contents");
}

// A rank-local identity is available where a transaction must bind its exact
// publication ownership and reference weights in addition to global topology.
inline std::string BuildDistributedSurfacePartitionIdentitySha256(const DistributedSurfaceLayout& layout)
{
	ValidateDistributedSurfaceLayout(layout);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "DistributedSurfacePartition/v1");
	distributed_surface_detail::AppendString(hash, layout.layout_identity_sha256);
	hash.AppendLittleEndian64(layout.partition_rank);
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(layout.owned_global_node_ids.size()));
	for (const auto id : layout.owned_global_node_ids) hash.AppendLittleEndian64(id);
	for (const double area : layout.owned_reference_lumped_areas_m2) hash.AppendNormalizedDouble(area);
	return hash.Hex();
}

inline void ValidateSurfaceFieldStamp(const SurfaceFieldStamp& stamp)
{
	if (!std::isfinite(stamp.time_s)) throw std::runtime_error("surface field stamp time must be finite");
	if (!IsLowercaseSha256(stamp.reference_mesh_identity_sha256)
		|| !IsLowercaseSha256(stamp.layout_identity_sha256)
		|| !IsLowercaseSha256(stamp.partition_identity_sha256)
		|| !IsLowercaseSha256(stamp.producer_state_identity_sha256))
		throw std::runtime_error("surface field stamp identities must be lowercase SHA-256 hashes");
}

inline void ValidateSurfaceFieldStamp(const SurfaceFieldStamp& stamp,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceFieldStamp(stamp);
	ValidateDistributedSurfaceLayout(layout);
	if (stamp.reference_mesh_identity_sha256 != layout.reference_mesh_identity_sha256
		|| stamp.layout_identity_sha256 != layout.layout_identity_sha256
		|| stamp.partition_identity_sha256 != BuildDistributedSurfacePartitionIdentitySha256(layout))
		throw std::runtime_error("surface field stamp does not match its layout");
}

inline void ValidateSurfaceFieldStampEnvelope(const SurfaceFieldStampEnvelope& envelope)
{
	if (!std::isfinite(envelope.time_s))
		throw std::runtime_error("surface field stamp envelope time must be finite");
	if (!IsLowercaseSha256(envelope.reference_mesh_identity_sha256)
		|| !IsLowercaseSha256(envelope.layout_identity_sha256)
		|| !IsLowercaseSha256(envelope.partition_identity_sha256))
		throw std::runtime_error("surface field stamp envelope identities must be lowercase SHA-256 hashes");
}

inline void ValidateSurfaceFieldStampEnvelope(const SurfaceFieldStampEnvelope& envelope,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceFieldStampEnvelope(envelope);
	ValidateDistributedSurfaceLayout(layout);
	if (envelope.reference_mesh_identity_sha256 != layout.reference_mesh_identity_sha256
		|| envelope.layout_identity_sha256 != layout.layout_identity_sha256
		|| envelope.partition_identity_sha256 != BuildDistributedSurfacePartitionIdentitySha256(layout))
		throw std::runtime_error("surface field stamp envelope does not match its layout");
}

inline SurfaceFieldStampEnvelope MakeSurfaceFieldStampEnvelope(const SurfaceFieldStamp& stamp)
{
	ValidateSurfaceFieldStamp(stamp);
	SurfaceFieldStampEnvelope envelope;
	envelope.time_s = stamp.time_s;
	envelope.step = stamp.step;
	envelope.coupling_iteration = stamp.coupling_iteration;
	envelope.reference_mesh_identity_sha256 = stamp.reference_mesh_identity_sha256;
	envelope.layout_identity_sha256 = stamp.layout_identity_sha256;
	envelope.partition_identity_sha256 = stamp.partition_identity_sha256;
	return envelope;
}

inline void ValidateSurfaceFieldStampMatchesEnvelope(const SurfaceFieldStamp& stamp,
	const SurfaceFieldStampEnvelope& envelope, const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceFieldStamp(stamp, layout);
	ValidateSurfaceFieldStampEnvelope(envelope, layout);
	if (stamp.time_s != envelope.time_s || stamp.step != envelope.step
		|| stamp.coupling_iteration != envelope.coupling_iteration
		|| stamp.reference_mesh_identity_sha256 != envelope.reference_mesh_identity_sha256
		|| stamp.layout_identity_sha256 != envelope.layout_identity_sha256
		|| stamp.partition_identity_sha256 != envelope.partition_identity_sha256)
		throw std::runtime_error("surface field stamp does not match its producer-neutral envelope");
}

inline std::string BuildSurfaceFieldStampIdentitySha256(const SurfaceFieldStamp& stamp,
	const DistributedSurfaceLayout& layout)
{
	// A stamp identity is meaningful only for the exact owned publication
	// slice it is being bound to.
	ValidateSurfaceFieldStamp(stamp, layout);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "SurfaceFieldStamp/v1");
	hash.AppendNormalizedDouble(stamp.time_s);
	hash.AppendLittleEndian64(stamp.step);
	hash.AppendLittleEndian64(stamp.coupling_iteration);
	distributed_surface_detail::AppendString(hash, stamp.reference_mesh_identity_sha256);
	distributed_surface_detail::AppendString(hash, stamp.layout_identity_sha256);
	distributed_surface_detail::AppendString(hash, stamp.partition_identity_sha256);
	distributed_surface_detail::AppendString(hash, stamp.producer_state_identity_sha256);
	return hash.Hex();
}

inline void ValidateSurfaceKinematics(const SurfaceKinematics& kinematics,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceInterfaceRef(kinematics.interface);
	ValidateSurfaceFieldStamp(kinematics.stamp, layout);
	const auto count = layout.owned_global_node_ids.size();
	if (kinematics.displacement_m.size() != count || kinematics.velocity_m_per_s.size() != count)
		throw std::runtime_error("surface kinematics vectors must align with owned node IDs");
	for (const auto& value : kinematics.displacement_m)
		distributed_surface_detail::RequireFiniteVector(value, "surface displacement");
	for (const auto& value : kinematics.velocity_m_per_s)
		distributed_surface_detail::RequireFiniteVector(value, "surface velocity");
}

inline void ValidateSurfaceTraction(const SurfaceTraction& traction,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceInterfaceRef(traction.interface);
	ValidateSurfaceFieldStamp(traction.stamp, layout);
	if (!IsLowercaseSha256(traction.projection_identity_sha256))
		throw std::runtime_error("surface traction projection identity must be a lowercase SHA-256 hash");
	const auto count = layout.owned_global_node_ids.size();
	if (traction.traction_on_structure_pa.size() != count || traction.consistent_nodal_force_n.size() != count)
		throw std::runtime_error("surface traction vectors must align with owned node IDs");
	for (const auto& value : traction.traction_on_structure_pa)
		distributed_surface_detail::RequireFiniteVector(value, "surface traction");
	for (const auto& value : traction.consistent_nodal_force_n)
		distributed_surface_detail::RequireFiniteVector(value, "surface nodal force");
}

inline std::string BuildSurfaceKinematicsIdentitySha256(const SurfaceKinematics& value,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceKinematics(value, layout);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "SurfaceKinematics/v2");
	distributed_surface_detail::AppendString(hash, value.interface.domain_id);
	distributed_surface_detail::AppendString(hash, value.interface.subsystem_id);
	distributed_surface_detail::AppendString(hash, value.interface.interface_id);
	distributed_surface_detail::AppendString(hash, BuildSurfaceFieldStampIdentitySha256(value.stamp, layout));
	for (const auto& vector : value.displacement_m) distributed_surface_detail::AppendVector(hash, vector);
	for (const auto& vector : value.velocity_m_per_s) distributed_surface_detail::AppendVector(hash, vector);
	return hash.Hex();
}

inline std::string BuildSurfaceTractionIdentitySha256(const SurfaceTraction& value,
	const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceTraction(value, layout);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "SurfaceTraction/v2");
	distributed_surface_detail::AppendString(hash, value.interface.domain_id);
	distributed_surface_detail::AppendString(hash, value.interface.subsystem_id);
	distributed_surface_detail::AppendString(hash, value.interface.interface_id);
	distributed_surface_detail::AppendString(hash, BuildSurfaceFieldStampIdentitySha256(value.stamp, layout));
	distributed_surface_detail::AppendString(hash, value.projection_identity_sha256);
	for (const auto& vector : value.traction_on_structure_pa) distributed_surface_detail::AppendVector(hash, vector);
	for (const auto& vector : value.consistent_nodal_force_n) distributed_surface_detail::AppendVector(hash, vector);
	return hash.Hex();
}

inline void ValidateDistributedSurfaceInterface(const DistributedSurfaceInterface& surface)
{
	ValidateSurfaceInterfaceRef(surface.id);
	if (surface.subsystem_id.empty()) throw std::runtime_error("surface interface subsystem ID must be nonempty");
	if (surface.id.subsystem_id != surface.subsystem_id)
		throw std::runtime_error("surface interface declaration subsystem does not match its reference");
	if (surface.locator_kind != "material_surface")
		throw std::runtime_error("surface interface locator kind must be material_surface");
	if (surface.boundary_labels.empty())
		throw std::runtime_error("material-surface interface boundary labels must be nonempty");
	if (!std::is_sorted(surface.boundary_labels.begin(), surface.boundary_labels.end())
		|| std::adjacent_find(surface.boundary_labels.begin(), surface.boundary_labels.end()) != surface.boundary_labels.end())
		throw std::runtime_error("surface interface boundary labels must be sorted and unique");
	for (const auto label : surface.boundary_labels)
		if (label < 0) throw std::runtime_error("surface interface boundary labels must be nonnegative");
	if (!IsLowercaseSha256(surface.reference_mesh_identity_sha256))
		throw std::runtime_error("surface interface reference mesh identity must be a lowercase SHA-256 hash");
	// The first fluid and structure surface contracts are both bidirectional.
	// Requiring each declaration prevents a one-way declaration from silently
	// satisfying the aggregate nonempty check below.
	if (surface.provides.empty()) throw std::runtime_error("surface interface provides list must be nonempty");
	if (surface.requires.empty()) throw std::runtime_error("surface interface requires list must be nonempty");
	std::vector<SurfaceFieldQuantity> quantities = surface.provides;
	quantities.insert(quantities.end(), surface.requires.begin(), surface.requires.end());
	for (const auto quantity : quantities)
		if (!IsKnownSurfaceFieldQuantity(quantity)) throw std::runtime_error("surface interface has an unknown quantity");
	for (const auto& declarations : {&surface.provides, &surface.requires}) {
		if (!std::is_sorted(declarations->begin(), declarations->end(), [](SurfaceFieldQuantity first, SurfaceFieldQuantity second) {
			return static_cast<std::uint8_t>(first) < static_cast<std::uint8_t>(second);
		})) throw std::runtime_error("surface interface quantity declarations must be canonical sorted");
		if (std::adjacent_find(declarations->begin(), declarations->end()) != declarations->end())
			throw std::runtime_error("surface interface has duplicate quantity declarations");
	}
	for (const auto provided : surface.provides)
		if (std::find(surface.requires.begin(), surface.requires.end(), provided) != surface.requires.end())
			throw std::runtime_error("surface interface cannot provide and require the same quantity");
}

inline std::string BuildDistributedSurfaceInterfaceIdentitySha256(const DistributedSurfaceInterface& surface)
{
	ValidateDistributedSurfaceInterface(surface);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "DistributedSurfaceInterface/v2");
	distributed_surface_detail::AppendString(hash, surface.id.domain_id);
	distributed_surface_detail::AppendString(hash, surface.id.subsystem_id);
	distributed_surface_detail::AppendString(hash, surface.id.interface_id);
	distributed_surface_detail::AppendString(hash, surface.subsystem_id);
	distributed_surface_detail::AppendString(hash, surface.locator_kind);
	distributed_surface_detail::AppendString(hash, surface.reference_mesh_identity_sha256);
	distributed_surface_detail::AppendString(hash, "boundary_labels");
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(surface.boundary_labels.size()));
	for (const auto label : surface.boundary_labels) hash.AppendLittleEndian64(static_cast<std::uint64_t>(label));
	distributed_surface_detail::AppendString(hash, "provides");
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(surface.provides.size()));
	for (const auto quantity : surface.provides) hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
	distributed_surface_detail::AppendString(hash, "requires");
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(surface.requires.size()));
	for (const auto quantity : surface.requires) hash.AppendLittleEndian32(static_cast<std::uint32_t>(quantity));
	return hash.Hex();
}

inline std::vector<double> MakeCartesianSurfaceLumpedAreaWeights(const DistributedSurfaceLayout& layout)
{
	ValidateDistributedSurfaceLayout(layout);
	std::vector<double> result;
	result.reserve(3*layout.owned_reference_lumped_areas_m2.size());
	for (const double area : layout.owned_reference_lumped_areas_m2)
		for (int component = 0; component < 3; ++component) result.push_back(area);
	return result;
}

} // namespace iga

#endif

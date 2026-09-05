#ifndef IGA_FSI_COUPLING_EDGE_HPP
#define IGA_FSI_COUPLING_EDGE_HPP

// Typed, field-valued FSI edge contracts.  These deliberately do not reuse
// CouplingEdge or PortRef: scalar pressure/flow edges and material-surface
// traction/kinematics edges have different endpoint and law semantics.
#include "DistributedSurfaceInterface.hpp"

#include <algorithm>
#include <cstdint>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace iga {

enum class FsiCouplingLaw : std::uint8_t {
	FluidStructureTractionKinematics
};

inline bool IsKnownFsiCouplingLaw(FsiCouplingLaw law)
{
	return law == FsiCouplingLaw::FluidStructureTractionKinematics;
}

inline const char* FsiCouplingLawName(FsiCouplingLaw law)
{
	if (law == FsiCouplingLaw::FluidStructureTractionKinematics)
		return "fluid_structure_traction_kinematics";
	return "unknown";
}

struct FsiCouplingEdge {
	std::string id;
	// Endpoint roles are explicit and never canonicalized or reordered.
	SurfaceInterfaceRef fluid;
	SurfaceInterfaceRef structure;
	FsiCouplingLaw law = FsiCouplingLaw::FluidStructureTractionKinematics;

	FsiCouplingEdge() = default;
	FsiCouplingEdge(std::string edge_id, SurfaceInterfaceRef fluid_endpoint,
		SurfaceInterfaceRef structure_endpoint, FsiCouplingLaw edge_law)
		: id(std::move(edge_id)), fluid(std::move(fluid_endpoint)),
		  structure(std::move(structure_endpoint)), law(edge_law) {}
};

inline bool operator==(const FsiCouplingEdge& first, const FsiCouplingEdge& second)
{
	return first.id == second.id && first.fluid == second.fluid
		&& first.structure == second.structure && first.law == second.law;
}

inline bool operator<(const FsiCouplingEdge& first, const FsiCouplingEdge& second)
{
	return std::tie(first.id, first.fluid, first.structure, first.law)
		< std::tie(second.id, second.fluid, second.structure, second.law);
}

inline void ValidateFsiCouplingEdge(const FsiCouplingEdge& edge)
{
	if (edge.id.empty()) throw std::runtime_error("FSI coupling edge id must be nonempty");
	ValidateSurfaceInterfaceRef(edge.fluid);
	ValidateSurfaceInterfaceRef(edge.structure);
	if (edge.fluid.domain_id == edge.structure.domain_id)
		throw std::runtime_error("FSI coupling edge fluid and structure endpoints must belong to distinct domains");
	if (!IsKnownFsiCouplingLaw(edge.law))
		throw std::runtime_error("FSI coupling edge has an unsupported law");
}

inline std::string BuildFsiCouplingEdgeIdentitySha256(const FsiCouplingEdge& edge)
{
	ValidateFsiCouplingEdge(edge);
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "FsiCouplingEdge/v2");
	distributed_surface_detail::AppendString(hash, edge.id);
	distributed_surface_detail::AppendString(hash, edge.fluid.domain_id);
	distributed_surface_detail::AppendString(hash, edge.fluid.subsystem_id);
	distributed_surface_detail::AppendString(hash, edge.fluid.interface_id);
	distributed_surface_detail::AppendString(hash, edge.structure.domain_id);
	distributed_surface_detail::AppendString(hash, edge.structure.subsystem_id);
	distributed_surface_detail::AppendString(hash, edge.structure.interface_id);
	hash.AppendLittleEndian32(static_cast<std::uint32_t>(edge.law));
	return hash.Hex();
}

inline void ValidateFsiFluidSurfaceInterface(const DistributedSurfaceInterface& surface)
{
	ValidateDistributedSurfaceInterface(surface);
	const std::vector<SurfaceFieldQuantity> expected_provides{
		SurfaceFieldQuantity::TractionOnStructure};
	const std::vector<SurfaceFieldQuantity> expected_requires{
		SurfaceFieldQuantity::Displacement, SurfaceFieldQuantity::Velocity};
	if (surface.provides != expected_provides || surface.requires != expected_requires)
		throw std::runtime_error("FSI fluid surface must provide traction and require displacement and velocity");
}

inline void ValidateFsiStructureSurfaceInterface(const DistributedSurfaceInterface& surface)
{
	ValidateDistributedSurfaceInterface(surface);
	const std::vector<SurfaceFieldQuantity> expected_provides{
		SurfaceFieldQuantity::Displacement, SurfaceFieldQuantity::Velocity};
	const std::vector<SurfaceFieldQuantity> expected_requires{
		SurfaceFieldQuantity::TractionOnStructure};
	if (surface.provides != expected_provides || surface.requires != expected_requires)
		throw std::runtime_error("FSI structure surface must provide displacement and velocity and require traction");
}

inline const DistributedSurfaceInterface& FindFsiSurfaceInterface(
	const std::vector<DistributedSurfaceInterface>& catalog, const SurfaceInterfaceRef& reference)
{
	ValidateSurfaceInterfaceRef(reference);
	const auto found = std::find_if(catalog.begin(), catalog.end(),
		[&reference](const DistributedSurfaceInterface& surface) { return surface.id == reference; });
	if (found == catalog.end()) throw std::runtime_error("FSI edge endpoint is absent from the surface interface catalog");
	return *found;
}

inline void ValidateFsiSurfaceInterfaceCatalog(
	const std::vector<DistributedSurfaceInterface>& catalog, const std::string& domain_id)
{
	if (domain_id.empty()) throw std::runtime_error("FSI surface catalog domain id must be nonempty");
	std::set<SurfaceInterfaceRef> references;
	for (const auto& surface : catalog) {
		ValidateDistributedSurfaceInterface(surface);
		if (surface.id.domain_id != domain_id)
			throw std::runtime_error("FSI surface catalog contains an interface from another domain");
		if (!references.insert(surface.id).second)
			throw std::runtime_error("FSI surface catalog contains duplicate interface IDs");
	}
}

inline void ValidateFsiCouplingEdgeInterfaces(const FsiCouplingEdge& edge,
	const DistributedSurfaceInterface& fluid_surface,
	const DistributedSurfaceInterface& structure_surface)
{
	ValidateFsiCouplingEdge(edge);
	ValidateFsiFluidSurfaceInterface(fluid_surface);
	ValidateFsiStructureSurfaceInterface(structure_surface);
	if (!(fluid_surface.id == edge.fluid) || !(structure_surface.id == edge.structure))
		throw std::runtime_error("FSI edge endpoint does not match its directional surface interface");
	if (fluid_surface.reference_mesh_identity_sha256
		!= structure_surface.reference_mesh_identity_sha256)
		throw std::runtime_error("FSI edge surface interfaces must have identical reference mesh identities");
	// Boundary labels are the topology selector available in this first slice.
	// Label zero is legal and participates in exact equality like any other label.
	if (fluid_surface.boundary_labels != structure_surface.boundary_labels)
		throw std::runtime_error("FSI edge surface interfaces have incompatible boundary topology");
}

inline void ValidateFsiCouplingEdgeLayouts(const FsiCouplingEdge& edge,
	const DistributedSurfaceInterface& fluid_surface,
	const DistributedSurfaceInterface& structure_surface,
	const DistributedSurfaceLayout& fluid_layout,
	const DistributedSurfaceLayout& structure_layout)
{
	ValidateFsiCouplingEdgeInterfaces(edge, fluid_surface, structure_surface);
	ValidateDistributedSurfaceLayout(fluid_layout);
	ValidateDistributedSurfaceLayout(structure_layout);
	if (fluid_layout.reference_mesh_identity_sha256 != fluid_surface.reference_mesh_identity_sha256
		|| structure_layout.reference_mesh_identity_sha256 != structure_surface.reference_mesh_identity_sha256)
		throw std::runtime_error("FSI surface layout does not match its interface reference mesh");
	if (fluid_layout.layout_identity_sha256 != structure_layout.layout_identity_sha256)
		throw std::runtime_error("FSI edge surface layouts have incompatible reference topology");
	if (BuildDistributedSurfacePartitionIdentitySha256(fluid_layout)
		!= BuildDistributedSurfacePartitionIdentitySha256(structure_layout))
		throw std::runtime_error("FSI edge surface layouts have incompatible rank-local ownership or reference weights");
}

inline void ValidateFsiCouplingEdgeBindings(const std::vector<FsiCouplingEdge>& edges,
	const std::vector<DistributedSurfaceInterface>& fluid_catalog,
	const std::vector<DistributedSurfaceInterface>& structure_catalog)
{
	if (edges.empty()) return;
	const std::string& fluid_domain_id = edges.front().fluid.domain_id;
	const std::string& structure_domain_id = edges.front().structure.domain_id;
	ValidateFsiSurfaceInterfaceCatalog(fluid_catalog, fluid_domain_id);
	ValidateFsiSurfaceInterfaceCatalog(structure_catalog, structure_domain_id);
	std::set<std::string> edge_ids;
	std::set<SurfaceInterfaceRef> bound_endpoints;
	for (const auto& edge : edges) {
		ValidateFsiCouplingEdge(edge);
		if (edge.fluid.domain_id != fluid_domain_id
			|| edge.structure.domain_id != structure_domain_id)
			throw std::runtime_error("FSI binding catalog cannot span multiple fluid or structure domains");
		if (!edge_ids.insert(edge.id).second)
			throw std::runtime_error("FSI coupling edge IDs must be unique");
		const auto& fluid_surface = FindFsiSurfaceInterface(fluid_catalog, edge.fluid);
		const auto& structure_surface = FindFsiSurfaceInterface(structure_catalog, edge.structure);
		ValidateFsiCouplingEdgeInterfaces(edge, fluid_surface, structure_surface);
		if (!bound_endpoints.insert(edge.fluid).second
			|| !bound_endpoints.insert(edge.structure).second)
			throw std::runtime_error("FSI surface endpoints must have one-to-one edge bindings");
	}
}

inline void ValidateFsiSurfaceFieldStampMatches(const SurfaceFieldStamp& value,
	const SurfaceFieldStamp& expected, const DistributedSurfaceLayout& layout)
{
	ValidateSurfaceFieldStamp(value, layout);
	ValidateSurfaceFieldStamp(expected, layout);
	if (BuildSurfaceFieldStampIdentitySha256(value, layout)
		!= BuildSurfaceFieldStampIdentitySha256(expected, layout))
		throw std::runtime_error("FSI surface field stamp is stale or belongs to another producer state");
}

inline void ValidateFsiFluidKinematicsInput(const FsiCouplingEdge& edge,
	const SurfaceKinematics& kinematics, const DistributedSurfaceLayout& fluid_layout,
	const SurfaceFieldStamp& expected_stamp)
{
	ValidateFsiCouplingEdge(edge);
	ValidateSurfaceKinematics(kinematics, fluid_layout);
	// Kinematics are published by the structure and consumed at the local
	// fluid endpoint selected by SetSurfaceKinematics(interface_id, ...).
	if (!(kinematics.interface == edge.structure))
		throw std::runtime_error("FSI fluid kinematics input names the wrong producer surface interface");
	ValidateFsiSurfaceFieldStampMatches(kinematics.stamp, expected_stamp, fluid_layout);
}

inline void ValidateFsiStructureTractionInput(const FsiCouplingEdge& edge,
	const SurfaceTraction& traction, const DistributedSurfaceLayout& structure_layout,
	const SurfaceFieldStamp& expected_stamp)
{
	ValidateFsiCouplingEdge(edge);
	ValidateSurfaceTraction(traction, structure_layout);
	// Traction is published by the fluid and consumed at the local structure
	// endpoint selected by SetSurfaceTraction(interface_id, ...).
	if (!(traction.interface == edge.fluid))
		throw std::runtime_error("FSI structure traction input names the wrong producer surface interface");
	ValidateFsiSurfaceFieldStampMatches(traction.stamp, expected_stamp, structure_layout);
}

} // namespace iga

#endif

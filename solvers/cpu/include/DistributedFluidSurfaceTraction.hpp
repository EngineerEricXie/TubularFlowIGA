#ifndef IGA_DISTRIBUTED_FLUID_SURFACE_TRACTION_HPP
#define IGA_DISTRIBUTED_FLUID_SURFACE_TRACTION_HPP

#include "OwnedFluidSurfaceTractionPoints.hpp"

namespace iga {
// Complete traction publication boundary. The caller owns the transaction and
// supplies independently captured expected state identities. No state commits
// occur here. Global conservation diagnostics remain a separate acceptance gate.
inline SurfaceTraction BuildDistributedFluidSurfaceTraction(MPI_Comm comm,
	const CartesianDomainClassification& domain,const ImmersedSurfaceQuadratureCatalog& catalog,
	const MaterialSurfaceKinematics& material,const MaterialSurfacePatchMap& map,
	const DistributedSurfaceInterface& surface,double viscosity,
	const std::vector<std::uint64_t>& owned_cells,const std::vector<FluidSurfaceElementState>& state,
	const FsiTrialContext& expected_context,const std::string& expected_material_identity,
	const std::string& expected_local_fluid_identity,int projection_owner,
	std::size_t maximum_nodes=4096,PointIdentityLimits limits={})
{
	std::string interface_identity;
	CollectiveLocalStage(comm,"distributed traction publication interface",[&] {
		ValidateDistributedSurfaceInterface(surface);
		if(surface.reference_mesh_identity_sha256!=map.ReferenceIdentitySha256()
			||surface.boundary_labels!=map.Interface().boundary_labels
			||surface.provides!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::TractionOnStructure}
			||surface.requires!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity})
			throw std::invalid_argument("distributed traction interface does not match mapped fluid role");
		interface_identity=BuildDistributedSurfaceInterfaceIdentitySha256(surface);
	});
	RequireCollectiveSameText(comm,"distributed traction publication interface",interface_identity);
	const auto points=BuildTrialFluidSurfaceTractionPoints(comm,domain,catalog,material,map,viscosity,owned_cells,state,
		expected_context,expected_material_identity,expected_local_fluid_identity,limits);
	auto values=AssembleDistributedSurfaceTraction(comm,surface.id,map.Layout(),domain.Cells().size(),points,projection_owner,maximum_nodes,limits);
	int ranks=0;MPI_Comm_size(comm,&ranks);
	std::string local_identity,all_identities;
	CollectiveLocalStage(comm,"distributed traction contributing identities",[&] {
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedFluidSurfaceTraction/contributor/v1");
		distributed_surface_detail::AppendString(hash,expected_local_fluid_identity);
		distributed_surface_detail::AppendString(hash,map.IdentitySha256());
		hash.AppendLittleEndian64(owned_cells.size());for(auto id:owned_cells)hash.AppendLittleEndian64(id);
		local_identity=hash.Hex();
		if(static_cast<std::size_t>(ranks)>all_identities.max_size()/64)throw std::overflow_error("traction contributor metadata overflow");
		all_identities.resize(64*static_cast<std::size_t>(ranks));
	});
	MPI_Allgather(local_identity.data(),64,MPI_CHAR,all_identities.data(),64,MPI_CHAR,comm);
	SurfaceTraction result;
	CollectiveLocalStage(comm,"distributed traction stamped result",[&] {
		result.interface=surface.id;
		result.stamp.time_s=expected_context.EndTime();result.stamp.step=expected_context.step;
		result.stamp.coupling_iteration=expected_context.coupling_iteration;
		result.stamp.reference_mesh_identity_sha256=map.Layout().reference_mesh_identity_sha256;
		result.stamp.layout_identity_sha256=map.Layout().layout_identity_sha256;
		result.stamp.partition_identity_sha256=BuildDistributedSurfacePartitionIdentitySha256(map.Layout());
		result.stamp.producer_state_identity_sha256=expected_local_fluid_identity;
		result.consistent_nodal_force_n=std::move(values.owned_force_n);
		result.traction_on_structure_pa=std::move(values.owned_traction_pa);
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedFluidSurfaceTraction/projection/v1");
		distributed_surface_detail::AppendString(hash,interface_identity);
		distributed_surface_detail::AppendString(hash,BuildSurfaceFieldStampIdentitySha256(result.stamp,map.Layout()));
		distributed_surface_detail::AppendString(hash,all_identities);
		hash.AppendLittleEndian64(projection_owner);hash.AppendLittleEndian64(maximum_nodes);
		for(const auto* field:{&result.consistent_nodal_force_n,&result.traction_on_structure_pa})
			for(const auto& tuple:*field)for(double value:tuple)hash.AppendNormalizedDouble(value);
		result.projection_identity_sha256=hash.Hex();
		ValidateSurfaceTraction(result,map.Layout());
	});
	return result;
}
} // namespace iga
#endif

#ifndef IGA_SINGLE_OWNER_SURFACE_TRACTION_HPP
#define IGA_SINGLE_OWNER_SURFACE_TRACTION_HPP

#include "SingleOwnerSurfaceLayout.hpp"
#include "SurfaceGhostTraction.hpp"

namespace iga {
struct SingleOwnerSurfaceTraction {
	DistributedSurfaceLayout layout;
	SurfaceTraction traction;
};

// Validate distributed publications FIRST, then derive a distinct internal
// single-owner input. Its producer/projection identities bind every original
// publication and the newly constructed numerical layout; no stamp is reused.
inline std::optional<SingleOwnerSurfaceTraction> GatherSurfaceTractionAtOwner(MPI_Comm comm,
	const DistributedSurfaceLayout& partition,const SurfaceTraction& publication,
	const SurfaceInterfaceRef& expected_interface,const SurfaceFieldStamp& expected_stamp,
	const std::string& expected_projection_identity,int owner,std::size_t maximum_nodes=4096,
	PointIdentityLimits limits={})
{
	auto layout=GatherSurfaceLayoutAtOwner(comm,expected_interface,partition,owner,maximum_nodes,limits);
	std::vector<std::uint64_t> requested;
	CollectiveLocalStage(comm,"single owner traction requests",[&] {
		if(layout)requested=layout->owned_global_node_ids;
	});
	auto values=FetchSurfaceGhostTraction(comm,partition,publication,expected_interface,expected_stamp,expected_projection_identity,requested,limits);
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	std::string local_identity,identities;
	CollectiveLocalStage(comm,"single owner traction provenance",[&] {
		local_identity=BuildSurfaceTractionIdentitySha256(publication,partition);
		if(rank==owner) {
			if(static_cast<std::size_t>(ranks)>identities.max_size()/64)throw std::overflow_error("single owner traction provenance overflow");
			identities.resize(64*static_cast<std::size_t>(ranks));
		}
	});
	MPI_Gather(local_identity.data(),64,MPI_CHAR,rank==owner?identities.data():nullptr,64,MPI_CHAR,owner,comm);
	std::optional<SingleOwnerSurfaceTraction> result;
	CollectiveLocalStage(comm,"single owner traction input",[&] {
		if(rank!=owner)return;
		result.emplace();result->layout=std::move(*layout);
		auto& traction=result->traction;
		traction.interface=expected_interface;
		traction.stamp.time_s=expected_stamp.time_s;traction.stamp.step=expected_stamp.step;
		traction.stamp.coupling_iteration=expected_stamp.coupling_iteration;
		traction.stamp.reference_mesh_identity_sha256=result->layout.reference_mesh_identity_sha256;
		traction.stamp.layout_identity_sha256=result->layout.layout_identity_sha256;
		traction.stamp.partition_identity_sha256=BuildDistributedSurfacePartitionIdentitySha256(result->layout);
		traction.traction_on_structure_pa=std::move(values.traction_on_structure_pa);
		traction.consistent_nodal_force_n=std::move(values.consistent_nodal_force_n);
		Sha256 producer;distributed_surface_detail::AppendString(producer,"SingleOwnerSurfaceTraction/producer/v1");
		distributed_surface_detail::AppendString(producer,identities);
		distributed_surface_detail::AppendString(producer,traction.stamp.partition_identity_sha256);
		producer.AppendLittleEndian64(owner);traction.stamp.producer_state_identity_sha256=producer.Hex();
		Sha256 projection;distributed_surface_detail::AppendString(projection,"SingleOwnerSurfaceTraction/projection/v1");
		distributed_surface_detail::AppendString(projection,BuildSurfaceFieldStampIdentitySha256(traction.stamp,result->layout));
		distributed_surface_detail::AppendString(projection,identities);
		for(const auto* field:{&traction.traction_on_structure_pa,&traction.consistent_nodal_force_n})
			for(const auto& tuple:*field)for(double value:tuple)projection.AppendNormalizedDouble(value);
		traction.projection_identity_sha256=projection.Hex();ValidateSurfaceTraction(traction,result->layout);
	});
	return result;
}
} // namespace iga
#endif

#ifndef IGA_SINGLE_OWNER_SURFACE_LAYOUT_HPP
#define IGA_SINGLE_OWNER_SURFACE_LAYOUT_HPP

#include "SurfaceOwnershipValidation.hpp"
#include "OwnedPointValues.hpp"
#include <optional>

namespace iga {
// Explicitly create a NEW single-partition solve layout at one chosen owner.
// Preserve physical IDs, reference geometry and the authoritative owned area
// values. The changed layout/partition hashes cannot be used for incoming MPI
// publications; their identities must be verified before a numerical bridge.
inline std::optional<DistributedSurfaceLayout> GatherSurfaceLayoutAtOwner(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& partition,int owner,
	std::size_t maximum_nodes=4096,PointIdentityLimits limits={})
{
	int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
	RequireCollectiveSameInt(comm,"single surface layout owner",owner);
	std::vector<std::uint64_t> requested;
	CollectiveLocalStage(comm,"single surface layout request",[&] {
		ValidateDistributedSurfaceLayout(partition);
		if(owner<0||owner>=ranks||!maximum_nodes||partition.global_node_count>maximum_nodes)
			throw std::invalid_argument("invalid single surface layout owner or node limit");
		if(rank==owner)for(const auto& node:partition.reference_positions)requested.push_back(node.global_node_id);
	});
	ValidateSurfacePublicationOwnership(reference,partition,comm);
	const auto areas=FetchOwnedPointValues(comm,partition.owned_global_node_ids,
		partition.owned_reference_lumped_areas_m2,requested,1,limits);
	std::optional<DistributedSurfaceLayout> result;
	CollectiveLocalStage(comm,"single surface solve layout",[&] {
		if(rank!=owner)return;
		result=partition;result->partition_count=1;result->partition_rank=0;
		result->owned_global_node_ids=std::move(requested);result->owned_reference_lumped_areas_m2=areas;
		result->layout_identity_sha256=BuildDistributedSurfaceLayoutIdentitySha256(*result);
		ValidateDistributedSurfaceLayout(*result);
	});
	return result;
}
} // namespace iga
#endif

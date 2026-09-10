#ifndef IGA_DISTRIBUTED_SURFACE_FORCES_HPP
#define IGA_DISTRIBUTED_SURFACE_FORCES_HPP

#include "SurfaceOwnershipValidation.hpp"
#include "OwnedVectorContributions.hpp"

namespace iga {
// Each uniquely owned reference triangle supplies three already-integrated P1
// corner forces, in reference triangle order. This does not approximate a
// traction projection or substitute lumped areas for the consistent mass solve.
inline std::vector<std::array<double,3>> AssembleOwnedSurfaceForces(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& layout,
	const std::vector<std::uint64_t>& owned_triangles,
	const std::vector<std::array<std::array<double,3>,3>>& corner_forces_n,
	PointIdentityLimits limits={})
{
	std::vector<std::pair<std::uint64_t,std::array<double,3>>> contributions;
	CollectiveLocalStage(comm,"surface corner force input",[&] {
		ValidateDistributedSurfaceLayout(layout);
		if(corner_forces_n.size()!=owned_triangles.size())throw std::invalid_argument("surface corner force size mismatch");
		if(!limits.max_local_occurrences||layout.owned_global_node_ids.size()>limits.max_local_occurrences
			||owned_triangles.size()>(limits.max_local_occurrences-layout.owned_global_node_ids.size())/3)
			throw std::invalid_argument("surface corner forces exceed record limit");
		contributions.reserve(3*owned_triangles.size());
		for(std::size_t row=0;row<owned_triangles.size();++row) {
			if(owned_triangles[row]>=layout.reference_triangles.size())throw std::invalid_argument("unknown force triangle");
			for(int corner=0;corner<3;++corner)
				contributions.emplace_back(layout.reference_triangles[owned_triangles[row]][corner],corner_forces_n[row][corner]);
		}
	});
	ValidateSurfacePublicationOwnership(reference,layout,comm);
	ValidateUniqueEntityCoverage(owned_triangles,layout.reference_triangles.size(),comm,"surface force triangles");
	return SumOwnedVectorContributions(comm,layout.owned_global_node_ids,contributions,limits);
}
} // namespace iga
#endif

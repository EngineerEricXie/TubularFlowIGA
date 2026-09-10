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

struct SurfaceCellForces {
	std::uint64_t cell_id=0;
	std::vector<std::pair<std::uint64_t,std::array<double,3>>> nodal_contributions_n;
};

// Every fluid cell in the dense background catalog has one owner, including
// cells with no retained interface quadrature (which supply an empty vector).
// A material node/triangle may receive contributions from many fluid cells.
// The caller integrates and stamps each cell's quadrature before this routing
// layer; this function certifies cell coverage, not quadrature provenance.
inline std::vector<std::array<double,3>> AssembleOwnedSurfaceCellForces(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& layout,
	std::uint64_t global_cell_count,const std::vector<SurfaceCellForces>& cells,
	PointIdentityLimits limits={})
{
	std::vector<std::uint64_t> owned_cells;
	std::vector<std::pair<std::uint64_t,std::array<double,3>>> contributions;
	CollectiveLocalStage(comm,"surface cell force input",[&] {
		ValidateDistributedSurfaceLayout(layout);
		if(!limits.max_local_occurrences||cells.size()>limits.max_local_occurrences
			||layout.owned_global_node_ids.size()>limits.max_local_occurrences)
			throw std::invalid_argument("surface cell forces exceed record limit");
		std::size_t count=layout.owned_global_node_ids.size();
		for(const auto& cell:cells) {
			if(cell.nodal_contributions_n.size()>limits.max_local_occurrences-count)
				throw std::invalid_argument("surface cell contributions exceed record limit");
			count+=cell.nodal_contributions_n.size();
		}
		owned_cells.reserve(cells.size());contributions.reserve(count-layout.owned_global_node_ids.size());
		for(const auto& cell:cells) {
			owned_cells.push_back(cell.cell_id);
			contributions.insert(contributions.end(),cell.nodal_contributions_n.begin(),cell.nodal_contributions_n.end());
		}
	});
	ValidateSurfacePublicationOwnership(reference,layout,comm);
	ValidateUniqueEntityCoverage(owned_cells,global_cell_count,comm,"surface force fluid cells");
	return SumOwnedVectorContributions(comm,layout.owned_global_node_ids,contributions,limits);
}
} // namespace iga
#endif

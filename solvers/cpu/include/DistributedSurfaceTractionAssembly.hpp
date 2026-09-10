#ifndef IGA_DISTRIBUTED_SURFACE_TRACTION_ASSEMBLY_HPP
#define IGA_DISTRIBUTED_SURFACE_TRACTION_ASSEMBLY_HPP

#include "DistributedSurfaceForces.hpp"
#include "DistributedSurfaceProjection.hpp"
#include "SurfaceP1TractionContribution.hpp"

namespace iga {
struct SurfaceP1TractionPoint {
	// Node IDs and barycentric values use the same canonical corner order.
	std::array<std::uint64_t,3> node_ids{};
	std::array<double,3> barycentric{},traction_pa{},physical_m{};
	double weight_m2=0.;
};
struct SurfaceCellTractionPoints {
	std::uint64_t cell_id=0;
	std::vector<SurfaceP1TractionPoint> points;
};
struct DistributedSurfaceTractionValues {
	std::vector<std::array<double,3>> owned_force_n,owned_traction_pa;
};

// Integrate force and consistent mass from the SAME local point records,
// certify unique fluid-cell coverage, then assemble/project to node owners.
// The fluid caller must certify point provenance, Cauchy stress and trial
// identity. This algebraic result deliberately is not a stamped publication.
inline DistributedSurfaceTractionValues AssembleDistributedSurfaceTraction(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& layout,
	std::uint64_t global_cell_count,const std::vector<SurfaceCellTractionPoints>& cells,
	int projection_owner,std::size_t maximum_nodes=4096,PointIdentityLimits limits={})
{
	std::vector<SurfaceCellForces> forces;
	std::vector<SurfaceMassEntry> mass;
	CollectiveLocalStage(comm,"surface P1 cell integration",[&] {
		ValidateDistributedSurfaceLayout(layout);
		if(!limits.max_local_occurrences||cells.size()>limits.max_local_occurrences)
			throw std::invalid_argument("surface traction cell record limit");
		std::size_t count=0;
		for(const auto& cell:cells) {
			if(cell.points.size()>(limits.max_local_occurrences-count)/9)
				throw std::invalid_argument("surface traction mass record limit");
			count+=9*cell.points.size();
		}
		forces.reserve(cells.size());mass.reserve(count);
		for(const auto& cell:cells) {
			SurfaceCellForces integrated;integrated.cell_id=cell.cell_id;
			integrated.nodal_contributions_n.reserve(3*cell.points.size());
			for(const auto& point:cell.points) {
				const auto contribution=BuildSurfaceP1TractionContribution(point.barycentric,point.traction_pa,point.weight_m2);
				for(int row=0;row<3;++row) {
					integrated.nodal_contributions_n.emplace_back(point.node_ids[row],contribution.corner_force_n[row]);
					for(int col=0;col<3;++col)mass.push_back({point.node_ids[row],point.node_ids[col],contribution.consistent_mass_m2[row][col]});
				}
			}
			forces.push_back(std::move(integrated));
		}
	});
	DistributedSurfaceTractionValues result;
	result.owned_force_n=AssembleOwnedSurfaceCellForces(comm,reference,layout,global_cell_count,forces,limits);
	result.owned_traction_pa=ProjectDistributedSurfaceTraction(comm,reference,layout,result.owned_force_n,mass,projection_owner,maximum_nodes,limits);
	return result;
}
} // namespace iga
#endif

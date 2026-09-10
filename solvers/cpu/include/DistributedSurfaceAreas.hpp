#ifndef IGA_DISTRIBUTED_SURFACE_AREAS_HPP
#define IGA_DISTRIBUTED_SURFACE_AREAS_HPP

#include "SurfaceOwnershipValidation.hpp"
#include "OwnedScalarContributions.hpp"

namespace iga {
struct DistributedSurfaceAreas {
	std::vector<double> owned_lumped_areas_m2;
	double global_area_m2=0.;
};

// Replicated reference geometry is retained, but each triangle is integrated
// by exactly one rank. Triangle indices refer to the agreed reference array;
// node ownership is independent of triangle ownership. Existing layout area
// values are not used in the calculation. The input/stamps remain unchanged.
inline DistributedSurfaceAreas ComputeDistributedSurfaceAreas(MPI_Comm comm,
	const SurfaceInterfaceRef& reference,const DistributedSurfaceLayout& layout,
	const std::vector<std::uint64_t>& owned_triangles,PointIdentityLimits limits={})
{
	CollectiveLocalStage(comm,"surface area record limits",[&] {
		if(!limits.max_local_occurrences||layout.owned_global_node_ids.size()>limits.max_local_occurrences
			||owned_triangles.size()>(limits.max_local_occurrences-layout.owned_global_node_ids.size())/3)
			throw std::invalid_argument("surface area contributions exceed record limit");
	});
	ValidateSurfacePublicationOwnership(reference,layout,comm);
	ValidateUniqueEntityCoverage(owned_triangles,layout.reference_triangles.size(),comm,"reference surface triangles");
	std::vector<std::pair<std::uint64_t,double>> contributions;
	long double local_area=0.;
	CollectiveLocalStage(comm,"surface reference triangle areas",[&] {
		contributions.reserve(3*owned_triangles.size());
		const auto position=[&](std::uint64_t id)->const std::array<double,3>& {
			const auto found=std::lower_bound(layout.reference_positions.begin(),layout.reference_positions.end(),id,
				[](const SurfaceReferencePosition& item,std::uint64_t value){return item.global_node_id<value;});
			if(found==layout.reference_positions.end()||found->global_node_id!=id)throw std::runtime_error("unknown surface area node");
			return found->position_m;
		};
		for(auto index:owned_triangles) {
			const auto& triangle=layout.reference_triangles[index];
			const auto& a=position(triangle[0]);const auto& b=position(triangle[1]);const auto& c=position(triangle[2]);
			std::array<double,3> u{},v{},cross{};
			for(int i=0;i<3;++i){u[i]=b[i]-a[i];v[i]=c[i]-a[i];}
			cross={{u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]}};
			const double area=.5*std::hypot(cross[0],cross[1],cross[2]);
			if(!std::isfinite(area)||!(area>0.))throw std::runtime_error("invalid reference triangle area");
			local_area+=static_cast<long double>(area);
			for(auto id:triangle)contributions.emplace_back(id,area/3.);
		}
	});
	DistributedSurfaceAreas result;
	result.owned_lumped_areas_m2=SumOwnedScalarContributions(comm,layout.owned_global_node_ids,contributions,limits);
	long double totals[2]{local_area,0.},global[2]{};
	CollectiveLocalStage(comm,"surface nodal area validation",[&] {
		for(double area:result.owned_lumped_areas_m2) {
			if(!std::isfinite(area)||!(area>0.))throw std::runtime_error("surface node has no positive reference area");
			totals[1]+=static_cast<long double>(area);
		}
		if(!std::isfinite(totals[0])||!std::isfinite(totals[1]))throw std::overflow_error("surface area sum overflow");
	});
	MPI_Allreduce(totals,global,2,MPI_LONG_DOUBLE,MPI_SUM,comm);
	CollectiveLocalStage(comm,"surface global area conservation",[&] {
		if(!std::isfinite(global[0])||!std::isfinite(global[1])||!(global[0]>0.)
			||global[0]>static_cast<long double>(std::numeric_limits<double>::max())
			||std::abs(global[0]-global[1])>128.L*std::numeric_limits<double>::epsilon()*global[0])
			throw std::runtime_error("reference triangle and nodal areas disagree");
		result.global_area_m2=static_cast<double>(global[0]);
	});
	return result;
}
} // namespace iga
#endif

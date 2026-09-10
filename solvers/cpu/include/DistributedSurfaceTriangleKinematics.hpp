#ifndef IGA_DISTRIBUTED_SURFACE_TRIANGLE_KINEMATICS_HPP
#define IGA_DISTRIBUTED_SURFACE_TRIANGLE_KINEMATICS_HPP

#include "SurfaceGhostKinematics.hpp"

namespace iga {
struct SurfaceTriangleKinematics {
	std::uint64_t reference_triangle_index=0;
	std::array<std::uint64_t,3> node_ids{};
	std::array<std::array<double,3>,3> positions_m{},velocities_m_per_s{};
};

// Reference topology remains replicated. Only nodes used by locally owned
// triangles are requested; triangle ownership is independent of node ownership.
// This evaluates P1 geometry, not closed-surface validity or time integration.
inline std::vector<SurfaceTriangleKinematics> BuildDistributedSurfaceTriangleKinematics(
	MPI_Comm comm,const DistributedSurfaceLayout& layout,const SurfaceKinematics& publication,
	const SurfaceInterfaceRef& expected_interface,const SurfaceFieldStamp& expected_stamp,
	const std::vector<std::uint64_t>& owned_triangles,PointIdentityLimits limits={})
{
	std::vector<std::uint64_t> requested;
	CollectiveLocalStage(comm,"surface triangle kinematics requests",[&] {
		ValidateDistributedSurfaceLayout(layout);
		if(!limits.max_local_occurrences||owned_triangles.size()>limits.max_local_occurrences/3)
			throw std::invalid_argument("surface triangle requests exceed record limit");
		requested.reserve(3*owned_triangles.size());
		for(auto index:owned_triangles) {
			if(index>=layout.reference_triangles.size())throw std::invalid_argument("unknown reference triangle");
			for(auto id:layout.reference_triangles[index])requested.push_back(id);
		}
		std::sort(requested.begin(),requested.end());
		requested.erase(std::unique(requested.begin(),requested.end()),requested.end());
	});
	ValidateUniqueEntityCoverage(owned_triangles,layout.reference_triangles.size(),comm,"kinematic surface triangles");
	const auto ghosts=FetchSurfaceGhostKinematics(comm,layout,publication,expected_interface,expected_stamp,requested,limits);
	std::vector<SurfaceTriangleKinematics> result;
	CollectiveLocalStage(comm,"surface triangle kinematics evaluation",[&] {
		result.resize(owned_triangles.size());
		for(std::size_t item=0;item<owned_triangles.size();++item) {
			auto& triangle=result[item];triangle.reference_triangle_index=owned_triangles[item];
			triangle.node_ids=layout.reference_triangles[owned_triangles[item]];
			for(int corner=0;corner<3;++corner) {
				const auto id=triangle.node_ids[corner];
				const auto row=std::lower_bound(requested.begin(),requested.end(),id)-requested.begin();
				const auto reference=std::lower_bound(layout.reference_positions.begin(),layout.reference_positions.end(),id,
					[](const SurfaceReferencePosition& value,std::uint64_t key){return value.global_node_id<key;});
				if(reference==layout.reference_positions.end()||reference->global_node_id!=id)
					throw std::runtime_error("surface triangle reference node missing");
				for(int axis=0;axis<3;++axis) {
					triangle.positions_m[corner][axis]=reference->position_m[axis]+ghosts.displacement_m[row][axis];
					if(!std::isfinite(triangle.positions_m[corner][axis]))throw std::runtime_error("nonfinite deformed surface position");
				}
				triangle.velocities_m_per_s[corner]=ghosts.velocity_m_per_s[row];
			}
		}
	});
	return result;
}
} // namespace iga
#endif

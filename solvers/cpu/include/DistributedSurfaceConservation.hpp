#ifndef IGA_DISTRIBUTED_SURFACE_CONSERVATION_HPP
#define IGA_DISTRIBUTED_SURFACE_CONSERVATION_HPP

#include "DistributedSurfaceTractionAssembly.hpp"

namespace iga {
struct DistributedSurfaceConservation {
	// Resultant xyz [N], moment xyz about origin [N m], power [W].
	std::array<long double,7> quadrature{},nodal{};
};

// Called after unique cell/node ownership has been certified. Physical point
// positions must come from the integration catalog. Power uses the material
// P1 nodal velocity field, i.e. discrete interface work, not fluid dissipation.
inline DistributedSurfaceConservation CheckDistributedSurfaceConservation(MPI_Comm comm,
	const std::vector<SurfaceCellTractionPoints>& cells,const MaterialSurfaceKinematics& material,
	const MaterialSurfacePatchMap& map,const std::vector<std::array<double,3>>& owned_force)
{
	long double local[28]{},global[28]{};
	CollectiveLocalStage(comm,"surface conservation local contributions",[&] {
		material.Validate();
		if(owned_force.size()!=map.Layout().owned_global_node_ids.size())throw std::invalid_argument("surface conservation force size mismatch");
		const auto accumulate=[&](int offset,const std::array<double,3>& x,const std::array<double,3>& v,
			const std::array<double,3>& f,double weight) {
			std::array<long double,3> force{};
			for(int axis=0;axis<3;++axis) {
				if(!std::isfinite(x[axis])||!std::isfinite(v[axis])||!std::isfinite(f[axis])||!std::isfinite(weight))
					throw std::invalid_argument("nonfinite surface conservation input");
				force[axis]=static_cast<long double>(f[axis])*weight;
			}
			const long double values[7]{force[0],force[1],force[2],
				x[1]*force[2]-x[2]*force[1],x[2]*force[0]-x[0]*force[2],x[0]*force[1]-x[1]*force[0],
				v[0]*force[0]+v[1]*force[1]+v[2]*force[2]};
			for(int i=0;i<7;++i) {
				local[offset+i]+=values[i];local[14+offset+i]+=std::abs(values[i]);
			}
		};
		for(const auto& cell:cells)for(const auto& point:cell.points) {
			std::array<double,3> velocity{};
			for(int corner=0;corner<3;++corner) {
				if(!std::isfinite(point.barycentric[corner]))throw std::invalid_argument("nonfinite work interpolation shape");
				const auto source=map.SourceVertexForGlobalNode(point.node_ids[corner]);
				for(int axis=0;axis<3;++axis)velocity[axis]+=point.barycentric[corner]*material.SourceVertexVelocitiesMPerS()[source][axis];
			}
			accumulate(0,point.physical_m,velocity,point.traction_pa,point.weight_m2);
		}
		for(std::size_t row=0;row<owned_force.size();++row) {
			const auto source=map.SourceVertexForGlobalNode(map.Layout().owned_global_node_ids[row]);
			accumulate(7,material.SourceVerticesM()[source],material.SourceVertexVelocitiesMPerS()[source],owned_force[row],1.);
		}
		for(auto value:local)if(!std::isfinite(value))throw std::overflow_error("surface conservation accumulation overflow");
	});
	MPI_Allreduce(local,global,28,MPI_LONG_DOUBLE,MPI_SUM,comm);
	DistributedSurfaceConservation result;
	CollectiveLocalStage(comm,"surface global force moment and work conservation",[&] {
		for(auto value:global)if(!std::isfinite(value))throw std::overflow_error("surface conservation global overflow");
		for(int i=0;i<7;++i) {
			result.quadrature[i]=global[i];result.nodal[i]=global[7+i];
			const long double tolerance=1.e-12L+2.e-10L*(global[14+i]+global[21+i]);
			if(std::abs(global[i]-global[7+i])>tolerance)throw std::runtime_error("distributed surface force, moment or power is not conserved");
		}
	});
	return result;
}
} // namespace iga
#endif

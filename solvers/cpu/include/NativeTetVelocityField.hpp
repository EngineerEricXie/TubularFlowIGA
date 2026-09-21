#ifndef IGA_NATIVE_TET_VELOCITY_FIELD_HPP
#define IGA_NATIVE_TET_VELOCITY_FIELD_HPP

#include "NativeTetFem.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

namespace iga {

// Evaluate the native P2 velocity at a physical point without extrapolating
// beyond the fitted tetrahedral domain. The first containing cell wins on faces.
class NativeTetVelocityField
{
public:
	NativeTetVelocityField(const NativeTetMesh& mesh,
		const NativeTaylorHoodTopology& topology,const std::vector<double>& state)
		: mesh_(mesh),topology_(topology),state_(state)
	{
		if(topology_.cell_velocity_nodes.size()!=mesh_.cells.size()
			||state_.size()!=3*(mesh_.points.size()+topology_.edges.size())+mesh_.points.size())
			throw std::invalid_argument("native tetrahedral velocity field dimensions are inconsistent");
		geometry_.reserve(mesh_.cells.size());
		for(const auto& cell:mesh_.cells)geometry_.push_back(EvaluateNativeTetGeometry(mesh_,cell));
	}

	std::optional<std::array<double,3>> At(const std::array<double,3>& point) const
	{
		for(double coordinate:point)if(!std::isfinite(coordinate))
			throw std::invalid_argument("native tetrahedral velocity query is nonfinite");
		constexpr double tolerance=2e-10;
		for(std::size_t cell=0;cell<mesh_.cells.size();++cell){
			const auto& tetra=mesh_.cells[cell];
			const auto& origin=mesh_.points[tetra.nodes[0]];
			std::array<double,4> barycentric{};
			for(std::size_t node=0;node<4;++node){
				barycentric[node]=node==0?1.0:0.0;
				for(int component=0;component<3;++component)
					barycentric[node]+=geometry_[cell].barycentric_gradients[node][component]
						*(point[component]-origin[component]);
			}
			bool inside=true;
			for(double value:barycentric)inside=inside&&value>=-tolerance&&value<=1.0+tolerance;
			if(!inside)continue;
			const auto basis=EvaluateNativeTaylorHoodBasis(
				barycentric[1],barycentric[2],barycentric[3]);
			std::array<double,3> velocity{};
			for(std::size_t local=0;local<10;++local){
				const auto global=topology_.cell_velocity_nodes[cell][local];
				for(int component=0;component<3;++component)
					velocity[component]+=basis.velocity[local]*state_[3*global+component];
			}
			return velocity;
		}
		return std::nullopt;
	}

private:
	const NativeTetMesh& mesh_;
	const NativeTaylorHoodTopology& topology_;
	const std::vector<double>& state_;
	std::vector<NativeTetGeometry> geometry_;
};

} // namespace iga

#endif

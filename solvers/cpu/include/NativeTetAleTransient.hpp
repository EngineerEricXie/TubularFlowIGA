#ifndef IGA_NATIVE_TET_ALE_TRANSIENT_HPP
#define IGA_NATIVE_TET_ALE_TRANSIENT_HPP

#include "NativeTetFem.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

// Backward-Euler material time derivative at fixed reference coordinate,
// integrated on the current geometry. The steady ALE operator supplies
// relative convection, viscosity, pressure, and incompressibility.
inline NativeTaylorHoodElementSystem BuildNativeTaylorHoodAleTransientElement(
	const NativeTetMesh& current_mesh,const NativeTetCell& cell,
	const std::array<double,NativeTaylorHoodElementSystem::dofs>& current_state,
	const std::array<double,NativeTaylorHoodElementSystem::dofs>& committed_state,
	const std::array<std::array<double,3>,4>& mesh_velocity_m_s,
	const NativeNavierStokesParameters& parameters,double dt_s)
{
	if(!(dt_s>0.0)||!std::isfinite(dt_s))
		throw std::invalid_argument("native ALE transient dt must be finite and positive");
	auto result=BuildNativeTaylorHoodAleNavierStokesElement(current_mesh,cell,
		current_state,mesh_velocity_m_s,parameters);
	const auto geometry=EvaluateNativeTetGeometry(current_mesh,cell);
	for(const auto& quadrature:NativeTetDegreeFiveQuadrature()) {
		const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(geometry,quadrature.reference);
		const double weighted_volume=quadrature.weight*geometry.determinant;
		for(std::size_t row_node=0;row_node<10;++row_node)
			for(int component=0;component<3;++component) {
				const std::size_t row=3*row_node+component;
				double acceleration=0.0;
				for(std::size_t column_node=0;column_node<10;++column_node) {
					const std::size_t column=3*column_node+component;
					acceleration+=basis.velocity[column_node]
						*(current_state[column]-committed_state[column])/dt_s;
					result.jacobian[row*NativeTaylorHoodElementSystem::dofs+column]
						+=weighted_volume*parameters.density*basis.velocity[row_node]
							*basis.velocity[column_node]/dt_s;
				}
				result.residual[row]+=weighted_volume*parameters.density
					*basis.velocity[row_node]*acceleration;
			}
	}
	return result;
}

} // namespace iga

#endif

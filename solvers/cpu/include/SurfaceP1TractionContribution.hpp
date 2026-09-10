#ifndef IGA_SURFACE_P1_TRACTION_CONTRIBUTION_HPP
#define IGA_SURFACE_P1_TRACTION_CONTRIBUTION_HPP

#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {
struct SurfaceP1TractionContribution {
	std::array<std::array<double,3>,3> corner_force_n{};
	std::array<std::array<double,3>,3> consistent_mass_m2{};
};

// One quadrature point in canonical triangle corner order. Quadrature/catalog
// authority validates the shape functions and signed integration rule; this
// arithmetic kernel checks finite inputs and products and preserves their order.
inline SurfaceP1TractionContribution BuildSurfaceP1TractionContribution(
	const std::array<double,3>& barycentric,const std::array<double,3>& traction_pa,
	double weight_m2)
{
	if(!std::isfinite(weight_m2))throw std::invalid_argument("nonfinite P1 traction weight");
	for(double value:barycentric)if(!std::isfinite(value))throw std::invalid_argument("nonfinite P1 traction shape");
	for(double value:traction_pa)if(!std::isfinite(value))throw std::invalid_argument("nonfinite P1 traction");
	SurfaceP1TractionContribution result;
	for(int left=0;left<3;++left) {
		for(int component=0;component<3;++component) {
			const double value=traction_pa[component]*barycentric[left]*weight_m2;
			if(!std::isfinite(value))throw std::overflow_error("nonfinite P1 corner force");
			result.corner_force_n[left][component]=value;
		}
		for(int right=0;right<3;++right) {
			const double value=barycentric[left]*barycentric[right]*weight_m2;
			if(!std::isfinite(value))throw std::overflow_error("nonfinite P1 consistent mass");
			result.consistent_mass_m2[left][right]=value;
		}
	}
	return result;
}
} // namespace iga
#endif

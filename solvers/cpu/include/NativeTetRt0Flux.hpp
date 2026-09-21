#ifndef IGA_NATIVE_TET_RT0_FLUX_HPP
#define IGA_NATIVE_TET_RT0_FLUX_HPP

#include "NativeTetFem.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace iga {

// The face opposite local vertex i has RT0 basis (x-x_i)/(3V).
// face_flow[cell][i] is its integrated outward flow in m^3/s.
inline std::array<double,3> EvaluateNativeTetRt0Flux(
	const NativeTetMesh& mesh,
	const std::vector<std::array<double,4>>& face_flow_m3_s,
	std::size_t cell_index,const std::array<double,3>& position_m)
{
	if(face_flow_m3_s.size()!=mesh.cells.size())
		throw std::invalid_argument("native tetra RT0 face-flow size is invalid");
	for(const double component:position_m)
		if(!std::isfinite(component))
			throw std::invalid_argument("native tetra RT0 position is nonfinite");
	const auto& cell=mesh.cells.at(cell_index);
	const auto geometry=EvaluateNativeTetGeometry(mesh,cell);
	std::array<double,3> flux{};
	for(std::size_t local=0;local<4;++local){
		const auto& vertex=mesh.points.at(cell.nodes[local]);
		const double flow=face_flow_m3_s[cell_index][local];
		if(!std::isfinite(flow))
			throw std::invalid_argument("native tetra RT0 face flow is nonfinite");
		for(int axis=0;axis<3;++axis)
			flux[axis]+=flow*(position_m[axis]-vertex[axis])
				*2./geometry.determinant;
	}
	return flux;
}

} // namespace iga

#endif

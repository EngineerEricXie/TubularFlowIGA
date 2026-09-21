#include "NativeTetAleMeshMotion.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {

std::array<double,3> AffineDisplacement(const std::array<double,3>& point)
{
	return {{0.1+0.2*point[0]-0.1*point[1],
		-0.05+0.3*point[1]+0.04*point[2],0.02-0.2*point[0]+0.1*point[2]}};
}

} // namespace

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.25,0.25,0.25}}}};
	mesh.cells = {iga::NativeTetCell{1,{{4,1,2,3}}},
		iga::NativeTetCell{2,{{0,4,2,3}}},iga::NativeTetCell{3,{{0,1,4,3}}},
		iga::NativeTetCell{4,{{0,1,2,4}}}};
	std::map<std::uint32_t,std::array<double,3>> boundary;
	for (std::uint32_t node = 0; node < 4; ++node)
		boundary[node] = AffineDisplacement(mesh.points[node]);
	const auto motion = iga::SolveNativeTetAleHarmonicMotion(mesh,boundary);
	assert(motion.constrained_nodes == 4);
	assert(motion.free_nodes == 1);
	assert(motion.maximum_free_residual < 1.0e-14);
	const auto expected = AffineDisplacement(mesh.points[4]);
	for (int component = 0; component < 3; ++component)
		assert(std::abs(motion.displacement_m[4][component]-expected[component]) < 1.0e-14);

	bool rejected = false;
	try {
		(void)iga::SolveNativeTetAleHarmonicMotion(mesh,
			{{9,{{0,0,0}}}});
	} catch (const std::invalid_argument&) { rejected = true; }
	assert(rejected);
	std::cout << "native tetrahedral ALE harmonic mesh-motion tests passed\n";
	return 0;
}

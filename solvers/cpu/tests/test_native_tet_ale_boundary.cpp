#include "NativeTetAleBoundary.hpp"

#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {

template<class Function> void Reject(Function&& function)
{
	bool rejected=false;
	try { function(); } catch (const std::runtime_error&) { rejected=true; }
	assert(rejected);
}

} // namespace

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.boundary_triangles={iga::NativeTetTriangle{1,{{0,1,2}},7},
		iga::NativeTetTriangle{2,{{0,1,3}},8}};
	const std::vector<iga::NativeTetAleBoundaryRule> rules{
		{7,iga::NativeTetAleBoundaryMode::PrescribedDisplacement},
		{8,iga::NativeTetAleBoundaryMode::Fixed}};
	const std::map<std::uint32_t,std::array<double,3>> compatible{
		{0,{{0,0,0}}},{1,{{0,0,0}}},{2,{{0.1,0,0}}}};
	const auto constraints=iga::ResolveNativeTetAleBoundaryConstraints(mesh,rules,compatible);
	assert(constraints.applied_labels.size()==2);
	assert(constraints.full_displacement_m.size()==4);
	assert(constraints.full_displacement_m.at(2)[0]==0.1);
	assert(constraints.full_displacement_m.at(3)[0]==0.0);

	auto conflicting=compatible;conflicting[0]={{0.1,0,0}};
	Reject([&] { (void)iga::ResolveNativeTetAleBoundaryConstraints(mesh,rules,conflicting); });
	Reject([&] { (void)iga::ResolveNativeTetAleBoundaryConstraints(mesh,
		{{7,iga::NativeTetAleBoundaryMode::SlidingNormalConstraint},
		 {8,iga::NativeTetAleBoundaryMode::Fixed}},compatible); });
	Reject([&] { (void)iga::ResolveNativeTetAleBoundaryConstraints(mesh,
		{{7,iga::NativeTetAleBoundaryMode::PrescribedDisplacement}},compatible); });
	std::cout<<"native tetrahedral ALE explicit boundary-rule tests passed\n";
	return 0;
}

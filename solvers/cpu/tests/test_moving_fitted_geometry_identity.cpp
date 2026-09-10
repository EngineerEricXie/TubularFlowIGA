#include "MovingCutGeometry.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include <iostream>
#include <cstring>
#include <stdexcept>
iga::RawSurfaceSoup Cube()
{
	iga::RawSurfaceSoup s;const double lo=.15,hi=.75;
	s.vertices={{{lo,lo,lo}},{{hi,lo,lo}},{{hi,hi,lo}},{{lo,hi,lo}},{{lo,lo,hi}},{{hi,lo,hi}},{{hi,hi,hi}},{{lo,hi,hi}}};
	for(const auto& t:std::vector<std::array<std::int64_t,3>>{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}})s.triangles.push_back({t,7});
	return s;
}
void Require(bool value,const char* message)
{
	if(!value)throw std::runtime_error(message);
}
int main()
{
	try {
		const auto cube=Cube();iga::PrescribedSurfaceMotion motion({{0.,cube},{1.,cube}});
		const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{1,1,1}}};
		iga::MovingCutGeometryOptions options;options.volume.max_depth=2;
		const auto legacy=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0.,0.,1.),options);
		// Captured with the identical fixture at pre-catalog revision 53b3238.
		Require(legacy->GeometryIdentitySha256()=="f2b2333e955112adb2663fb95214d01981b9fa2f195bd382b8809d44ee8140b1","default geometry identity changed");
		options.volume_fitting.emplace();
		const auto fitted=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0.,0.,1.),options);
		const auto repeat=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0.,0.,1.),options);
		Require(fitted->GeometryIdentitySha256()!=legacy->GeometryIdentitySha256(),"fitting policy was not bound by geometry identity");
		Require(fitted->GeometryIdentitySha256()==repeat->GeometryIdentitySha256(),"fitted geometry identity is not reproducible");
		auto changed_options=options;++changed_options.volume_fitting->fit.max_iterations;
		const auto changed=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0.,0.,1.),changed_options);
		Require(changed->GeometryIdentitySha256()!=fitted->GeometryIdentitySha256(),"fitting limit was not bound by geometry identity");
		const auto& original=fitted->Volume().Cell(0).rule.Points();const auto& same_rule=changed->Volume().Cell(0).rule.Points();
		Require(original.size()==same_rule.size(),"unused iteration budget changed the rule size");
		for(std::size_t i=0;i<original.size();++i)
			Require(std::memcmp(original[i].parametric.data(),same_rule[i].parametric.data(),3*sizeof(double))==0
				&&std::memcmp(&original[i].weight,&same_rule[i].weight,sizeof(double))==0,"unused iteration budget changed rule values");
		const auto identity=fitted->GeometryIdentitySha256();
		const auto saved_points=original;
		changed_options.volume_fitting->max_point_queries=0;bool rejected=false;
		try { iga::MovingCutGeometry::Build(grid,motion.Evaluate(1.,0.,1.),changed_options,fitted.get()); }
		catch(const std::exception&) { rejected=true; }
		Require(rejected&&identity==fitted->GeometryIdentitySha256(),"failed fitted build changed the prior geometry");
		const auto& preserved=fitted->Volume().Cell(0).rule.Points();
		Require(preserved.size()==saved_points.size(),"failed fitted build changed prior rule size");
		for(std::size_t i=0;i<saved_points.size();++i)
			Require(std::memcmp(saved_points[i].parametric.data(),preserved[i].parametric.data(),3*sizeof(double))==0
				&&std::memcmp(&saved_points[i].weight,&preserved[i].weight,sizeof(double))==0,"failed fitted build changed prior rule values");
		std::cout<<"moving_fitted_geometry_identity_test: PASS legacy golden identity, fitted determinism, policy binding and failed-build isolation\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}

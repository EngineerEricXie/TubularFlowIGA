#include "NativeIgaShellSurface.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace {

bool Near(double a,double b,double tolerance=2.0e-11)
{
	return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

iga::NativeNurbsShellPatch Plane()
{
	iga::NativeNurbsShellPatch patch;
	patch.degree_u=patch.degree_v=2;
	patch.knots_u=patch.knots_v={0,0,0,1,1,1};
	patch.controls_u=patch.controls_v=3;
	for(std::size_t j=0;j<3;++j)
		for(std::size_t i=0;i<3;++i) {
			patch.control_points.push_back({{0.5*i,0.5*j,0.0}});
			patch.weights.push_back(1.0);
		}
	return patch;
}

}

int main()
{
	const auto plane=Plane();
	for(const auto parameter:std::array<std::array<double,2>,3>{{{{0.17,0.23}},{{0.5,0.5}},{{0.83,0.71}}}}) {
		const auto point=iga::EvaluateNativeNurbsShellSurface(plane,parameter[0],parameter[1]);
		assert(Near(std::accumulate(point.basis.begin(),point.basis.end(),0.0),1.0));
		assert(Near(std::accumulate(point.basis_u.begin(),point.basis_u.end(),0.0),0.0));
		assert(Near(std::accumulate(point.basis_v.begin(),point.basis_v.end(),0.0),0.0));
		assert(Near(std::accumulate(point.basis_uu.begin(),point.basis_uu.end(),0.0),0.0));
		assert(Near(std::accumulate(point.basis_uv.begin(),point.basis_uv.end(),0.0),0.0));
		assert(Near(std::accumulate(point.basis_vv.begin(),point.basis_vv.end(),0.0),0.0));
		assert(Near(point.position[0],parameter[0])&&Near(point.position[1],parameter[1]));
		assert(Near(point.tangent_u[0],1.0)&&Near(point.tangent_v[1],1.0));
		assert(Near(point.surface_jacobian,1.0)&&Near(point.unit_normal[2],1.0));
		for(double value:point.second_uu) assert(Near(value,0.0));
		for(double value:point.second_uv) assert(Near(value,0.0));
		for(double value:point.second_vv) assert(Near(value,0.0));
		const double h=1.0e-5;
		const auto plus=iga::EvaluateNativeNurbsShellSurface(plane,parameter[0]+h,parameter[1]);
		const auto minus=iga::EvaluateNativeNurbsShellSurface(plane,parameter[0]-h,parameter[1]);
		for(std::size_t i=0;i<point.basis.size();++i) {
			assert(Near((plus.basis[i]-minus.basis[i])/(2*h),point.basis_u[i],2.0e-9));
			assert(Near((plus.basis_u[i]-minus.basis_u[i])/(2*h),point.basis_uu[i],2.0e-8));
		}
	}
	auto quarter=Plane();
	const double q=std::sqrt(0.5);
	quarter.control_points.clear();quarter.weights.clear();
	for(std::size_t j=0;j<3;++j) {
		quarter.control_points.push_back({{1,0,0.5*j}});quarter.weights.push_back(1);
		quarter.control_points.push_back({{1,1,0.5*j}});quarter.weights.push_back(q);
		quarter.control_points.push_back({{0,1,0.5*j}});quarter.weights.push_back(1);
	}
	const auto cylinder=iga::EvaluateNativeNurbsShellSurface(quarter,0.5,0.37);
	assert(Near(cylinder.position[0],q)&&Near(cylinder.position[1],q));
	assert(Near(cylinder.position[2],0.37));
	assert(cylinder.surface_jacobian>0.0);
	auto invalid=plane;invalid.weights[0]=0.0;
	try { iga::ValidateNativeNurbsShellPatch(invalid);assert(false); }
	catch(const std::invalid_argument&) {}
	invalid=plane;invalid.knots_u={0,0,0,0.5,0.5,1,1,1};invalid.controls_u=5;
	invalid.control_points.resize(15);invalid.weights.assign(15,1.0);
	try { iga::ValidateNativeNurbsShellPatch(invalid);assert(false); }
	catch(const std::invalid_argument&) {}
	std::cout<<"native IGA shell surface basis tests passed\n";
	return 0;
}

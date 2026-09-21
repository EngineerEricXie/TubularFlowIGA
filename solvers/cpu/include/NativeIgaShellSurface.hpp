#ifndef IGA_NATIVE_IGA_SHELL_SURFACE_HPP
#define IGA_NATIVE_IGA_SHELL_SURFACE_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeSplineBasisDerivatives
{
	std::vector<double> value;
	std::vector<double> first;
	std::vector<double> second;
};

inline void ValidateNativeShellKnotVector(const std::vector<double>& knots,
	int degree, bool require_c1=true)
{
	if(degree<2) throw std::invalid_argument("native shell spline degree must be at least two");
	if(knots.size()<static_cast<std::size_t>(2*degree+2))
		throw std::invalid_argument("native shell knot vector is too short");
	for(std::size_t i=0;i<knots.size();++i) {
		if(!std::isfinite(knots[i])) throw std::invalid_argument("native shell knot is nonfinite");
		if(i&&knots[i]<knots[i-1])
			throw std::invalid_argument("native shell knot vector is not nondecreasing");
	}
	for(int i=1;i<=degree;++i)
		if(knots[static_cast<std::size_t>(i)]!=knots.front()
			||knots[knots.size()-1-static_cast<std::size_t>(i)]!=knots.back())
			throw std::invalid_argument("native shell knot vector must be open");
	if(!(knots[static_cast<std::size_t>(degree)]
		<knots[knots.size()-1-static_cast<std::size_t>(degree)]))
		throw std::invalid_argument("native shell knot domain is empty");
	if(require_c1) {
		std::size_t i=static_cast<std::size_t>(degree+1);
		const std::size_t end=knots.size()-static_cast<std::size_t>(degree+1);
		while(i<end) {
			std::size_t next=i+1;
			while(next<end&&knots[next]==knots[i]) ++next;
			if(next-i>static_cast<std::size_t>(degree-1))
				throw std::invalid_argument("native Kirchhoff-Love shell has a non-C1 internal knot");
			i=next;
		}
	}
}

inline NativeSplineBasisDerivatives EvaluateNativeSplineBasis(
	const std::vector<double>& knots,int degree,double parameter)
{
	ValidateNativeShellKnotVector(knots,degree,false);
	const std::size_t count=knots.size()-static_cast<std::size_t>(degree)-1;
	const double lower=knots[static_cast<std::size_t>(degree)];
	const double upper=knots[count];
	if(!std::isfinite(parameter)||parameter<lower||parameter>upper)
		throw std::invalid_argument("native shell parameter is outside the knot domain");
	if(parameter==upper) parameter=std::nextafter(upper,lower);
	std::vector<double> previous_value(knots.size(),0.0), previous_first(knots.size(),0.0);
	for(std::size_t i=0;i+1<knots.size();++i)
		if(knots[i]<=parameter&&parameter<knots[i+1]) previous_value[i]=1.0;
	for(int p=1;p<=degree;++p) {
		std::vector<double> value(knots.size(),0.0),first(knots.size(),0.0),second(knots.size(),0.0);
		for(std::size_t i=0;i+static_cast<std::size_t>(p)+1<knots.size();++i) {
			const double left=knots[i+static_cast<std::size_t>(p)]-knots[i];
			const double right=knots[i+static_cast<std::size_t>(p)+1]-knots[i+1];
			if(left>0.0) {
				value[i]+=(parameter-knots[i])/left*previous_value[i];
				first[i]+=p/left*previous_value[i];
				second[i]+=p/left*previous_first[i];
			}
			if(right>0.0) {
				value[i]+=(knots[i+static_cast<std::size_t>(p)+1]-parameter)/right
					*previous_value[i+1];
				first[i]-=p/right*previous_value[i+1];
				second[i]-=p/right*previous_first[i+1];
			}
		}
		previous_value.swap(value);
		previous_first.swap(first);
		if(p==degree) {
			NativeSplineBasisDerivatives result;
			result.value.assign(previous_value.begin(),previous_value.begin()+count);
			result.first.assign(previous_first.begin(),previous_first.begin()+count);
			result.second.assign(second.begin(),second.begin()+count);
			return result;
		}
	}
	throw std::logic_error("native shell basis evaluation did not reach its degree");
}

struct NativeNurbsShellPatch
{
	int degree_u=0,degree_v=0;
	std::vector<double> knots_u,knots_v;
	std::size_t controls_u=0,controls_v=0;
	std::vector<std::array<double,3>> control_points;
	std::vector<double> weights;
};

struct NativeNurbsShellSurfacePoint
{
	std::vector<double> basis,basis_u,basis_v,basis_uu,basis_uv,basis_vv;
	std::array<double,3> position{},tangent_u{},tangent_v{};
	std::array<double,3> second_uu{},second_uv{},second_vv{},unit_normal{};
	double surface_jacobian=0.0;
};

inline void ValidateNativeNurbsShellPatch(const NativeNurbsShellPatch& patch)
{
	ValidateNativeShellKnotVector(patch.knots_u,patch.degree_u);
	ValidateNativeShellKnotVector(patch.knots_v,patch.degree_v);
	const std::size_t expected_u=patch.knots_u.size()-static_cast<std::size_t>(patch.degree_u)-1;
	const std::size_t expected_v=patch.knots_v.size()-static_cast<std::size_t>(patch.degree_v)-1;
	if(patch.controls_u!=expected_u||patch.controls_v!=expected_v)
		throw std::invalid_argument("native shell control-grid dimensions disagree with knots");
	if(patch.controls_u&&patch.controls_v>
		std::numeric_limits<std::size_t>::max()/patch.controls_u)
		throw std::overflow_error("native shell control-grid size overflows");
	const std::size_t count=patch.controls_u*patch.controls_v;
	if(patch.control_points.size()!=count||patch.weights.size()!=count)
		throw std::invalid_argument("native shell control points or weights have the wrong size");
	for(const auto& point:patch.control_points)
		for(const double value:point)
			if(!std::isfinite(value)) throw std::invalid_argument("native shell control point is nonfinite");
	for(const double weight:patch.weights)
		if(!(weight>0.0)||!std::isfinite(weight))
			throw std::invalid_argument("native shell NURBS weight must be finite and positive");
}

inline NativeNurbsShellSurfacePoint EvaluateNativeNurbsShellSurface(
	const NativeNurbsShellPatch& patch,double u,double v)
{
	ValidateNativeNurbsShellPatch(patch);
	const auto bu=EvaluateNativeSplineBasis(patch.knots_u,patch.degree_u,u);
	const auto bv=EvaluateNativeSplineBasis(patch.knots_v,patch.degree_v,v);
	const std::size_t count=patch.controls_u*patch.controls_v;
	std::vector<double> a(count),au(count),av(count),auu(count),auv(count),avv(count);
	double w=0.0,wu=0.0,wv=0.0,wuu=0.0,wuv=0.0,wvv=0.0;
	for(std::size_t j=0;j<patch.controls_v;++j)
		for(std::size_t i=0;i<patch.controls_u;++i) {
			const std::size_t index=j*patch.controls_u+i;
			const double weight=patch.weights[index];
			a[index]=weight*bu.value[i]*bv.value[j];
			au[index]=weight*bu.first[i]*bv.value[j];
			av[index]=weight*bu.value[i]*bv.first[j];
			auu[index]=weight*bu.second[i]*bv.value[j];
			auv[index]=weight*bu.first[i]*bv.first[j];
			avv[index]=weight*bu.value[i]*bv.second[j];
			w+=a[index];wu+=au[index];wv+=av[index];
			wuu+=auu[index];wuv+=auv[index];wvv+=avv[index];
		}
	if(!(w>0.0)||!std::isfinite(w)) throw std::runtime_error("native shell rational weight sum is invalid");
	NativeNurbsShellSurfacePoint result;
	result.basis.resize(count);result.basis_u.resize(count);result.basis_v.resize(count);
	result.basis_uu.resize(count);result.basis_uv.resize(count);result.basis_vv.resize(count);
	for(std::size_t index=0;index<count;++index) {
		result.basis[index]=a[index]/w;
		result.basis_u[index]=(au[index]-result.basis[index]*wu)/w;
		result.basis_v[index]=(av[index]-result.basis[index]*wv)/w;
		result.basis_uu[index]=(auu[index]-2.0*result.basis_u[index]*wu
			-result.basis[index]*wuu)/w;
		result.basis_uv[index]=(auv[index]-result.basis_u[index]*wv
			-result.basis_v[index]*wu-result.basis[index]*wuv)/w;
		result.basis_vv[index]=(avv[index]-2.0*result.basis_v[index]*wv
			-result.basis[index]*wvv)/w;
		for(std::size_t component=0;component<3;++component) {
			const double coordinate=patch.control_points[index][component];
			result.position[component]+=result.basis[index]*coordinate;
			result.tangent_u[component]+=result.basis_u[index]*coordinate;
			result.tangent_v[component]+=result.basis_v[index]*coordinate;
			result.second_uu[component]+=result.basis_uu[index]*coordinate;
			result.second_uv[component]+=result.basis_uv[index]*coordinate;
			result.second_vv[component]+=result.basis_vv[index]*coordinate;
		}
	}
	const auto& x=result.tangent_u;const auto& y=result.tangent_v;
	std::array<double,3> normal{{x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],
		x[0]*y[1]-x[1]*y[0]}};
	result.surface_jacobian=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	if(!(result.surface_jacobian>0.0)||!std::isfinite(result.surface_jacobian))
		throw std::runtime_error("native shell midsurface has a singular metric");
	for(std::size_t component=0;component<3;++component)
		result.unit_normal[component]=normal[component]/result.surface_jacobian;
	return result;
}

} // namespace iga

#endif

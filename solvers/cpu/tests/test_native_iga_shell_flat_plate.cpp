#include "NativeIgaShellFlatPlate.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

std::vector<double> Knots(int elements,int degree)
{
	std::vector<double> result(static_cast<std::size_t>(degree+1),0.0);
	for(int i=1;i<elements;++i) result.push_back(static_cast<double>(i)/elements);
	result.insert(result.end(),static_cast<std::size_t>(degree+1),1.0);return result;
}

iga::NativeNurbsShellPatch Plane(int elements)
{
	iga::NativeNurbsShellPatch patch;patch.degree_u=patch.degree_v=3;
	patch.knots_u=patch.knots_v=Knots(elements,3);
	patch.controls_u=patch.controls_v=static_cast<std::size_t>(elements+3);
	for(std::size_t j=0;j<patch.controls_v;++j) for(std::size_t i=0;i<patch.controls_u;++i) {
		double x=0.0,y=0.0;
		for(int k=1;k<=3;++k) { x+=patch.knots_u[i+static_cast<std::size_t>(k)]/3.0;
			y+=patch.knots_v[j+static_cast<std::size_t>(k)]/3.0; }
		patch.control_points.push_back({{x,y,0}});patch.weights.push_back(1);
	}
	return patch;
}

struct Error { double l2=0.0,h2=0.0; };

Error Run(int elements)
{
	const auto patch=Plane(elements);
	const iga::NativeKirchhoffLoveMaterial material{2e6,0.3,0.01,1000};
	const double pi=std::acos(-1.0);
	const double rigidity=material.young_modulus_pa*std::pow(material.thickness_m,3)/
		(12.0*(1.0-material.poisson_ratio*material.poisson_ratio));
	const auto solution=iga::SolveNativeIgaSimplySupportedFlatPlate(patch,material,
		[&](double x,double y) { return 4.0*std::pow(pi,4)*rigidity*std::sin(pi*x)*std::sin(pi*y); });
	const std::array<double,4> gx{{-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526}};
	const std::array<double,4> gw{{0.3478548451374538,0.6521451548625461,0.6521451548625461,0.3478548451374538}};
	double e0=0.0,r0=0.0,e2=0.0,r2=0.0;
	for(int ey=0;ey<elements;++ey) for(int ex=0;ex<elements;++ex)
		for(std::size_t j=0;j<4;++j) for(std::size_t i=0;i<4;++i) {
			const double x=(ex+0.5*(1+gx[i]))/elements,y=(ey+0.5*(1+gx[j]))/elements;
			const auto point=iga::EvaluateNativeNurbsShellSurface(patch,x,y);
			double value=0.0,xx=0.0,xy=0.0,yy=0.0;
			for(std::size_t a=0;a<point.basis.size();++a) {
				const double coefficient=solution.transverse_control_displacement_m[a];
				value+=coefficient*point.basis[a];xx+=coefficient*point.basis_uu[a];
				xy+=coefficient*point.basis_uv[a];yy+=coefficient*point.basis_vv[a];
			}
			const double exact=std::sin(pi*x)*std::sin(pi*y);
			const double exact_xx=-pi*pi*exact,exact_yy=exact_xx;
			const double exact_xy=pi*pi*std::cos(pi*x)*std::cos(pi*y);
			const double weight=gw[i]*gw[j]/(4.0*elements*elements);
			e0+=weight*(value-exact)*(value-exact);r0+=weight*exact*exact;
			e2+=weight*((xx-exact_xx)*(xx-exact_xx)+2*(xy-exact_xy)*(xy-exact_xy)
				+(yy-exact_yy)*(yy-exact_yy));
			r2+=weight*(exact_xx*exact_xx+2*exact_xy*exact_xy+exact_yy*exact_yy);
		}
	return {std::sqrt(e0/r0),std::sqrt(e2/r2)};
}

}

int main()
{
	const auto coarse=Run(2),medium=Run(4),fine=Run(8);
	assert(coarse.l2>medium.l2&&medium.l2>fine.l2);
	assert(coarse.h2>medium.h2&&medium.h2>fine.h2);
	const double l2a=std::log(coarse.l2/medium.l2)/std::log(2.0);
	const double l2b=std::log(medium.l2/fine.l2)/std::log(2.0);
	const double h2a=std::log(coarse.h2/medium.h2)/std::log(2.0);
	const double h2b=std::log(medium.h2/fine.h2)/std::log(2.0);
	assert(l2a>3.0&&l2b>3.0&&h2a>1.5&&h2b>1.5);
	std::cout<<"native IGA loaded flat-plate convergence passed L2="<<coarse.l2<<','
		<<medium.l2<<','<<fine.l2<<" H2="<<coarse.h2<<','<<medium.h2<<','<<fine.h2
		<<" orders="<<l2a<<','<<l2b<<','<<h2a<<','<<h2b<<'\n';
	return 0;
}

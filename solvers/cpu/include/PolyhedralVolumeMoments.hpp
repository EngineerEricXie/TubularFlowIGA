#ifndef IGA_POLYHEDRAL_VOLUME_MOMENTS_HPP
#define IGA_POLYHEDRAL_VOLUME_MOMENTS_HPP

#include "SurfaceGeometry.hpp"

namespace iga {

namespace polyhedral_moment_detail {
inline long double Power(long double x, unsigned exponent)
{
	long double value=1;
	for(unsigned i=0;i<exponent;++i) value*=x;
	return value;
}

// Gauss-Legendre on [0,1]. The bounded order covers tensor-cubic products
// and their antiderivatives without an external quadrature dependency.
inline std::vector<std::pair<long double,long double>> Gauss(unsigned count)
{
	std::vector<std::pair<long double,long double>> result(count);
	for(unsigned i=0;i<(count+1)/2;++i) {
		long double x=std::cos(std::acos(-1.L)*(i+.75L)/(count+.5L));
		long double derivative=0;
		bool converged=false;
		for(unsigned iteration=0;iteration<64;++iteration) {
			long double previous=1, value=x;
			for(unsigned k=2;k<=count;++k) {
				const long double next=((2*k-1)*x*value-(k-1)*previous)/k;
				previous=value; value=next;
			}
			derivative=count*(x*value-previous)/(x*x-1);
			const long double correction=value/derivative;
			if(std::abs(correction)<=8*std::numeric_limits<long double>::epsilon()) { converged=true; break; }
			x-=correction;
		}
		if(!converged) throw std::runtime_error("polyhedral moment Gauss root did not converge");
		const long double weight=1/((1-x*x)*derivative*derivative);
		result[i]={(1-x)/2,weight}; result[count-1-i]={(1+x)/2,weight};
	}
	return result;
}
}

// Physical integral of product_q ((x_q-origin_q)/scale_q)^exponent_q.
// The validated closed surface supplies orientation and topology. No convexity
// assumption is made. This computes target moments, not positive volume rules.
inline long double PolyhedralVolumeMoment(const ClosedTriangulatedSurface& surface,
	const std::array<unsigned,3>& exponent,
	const std::array<double,3>& origin={{0,0,0}},
	const std::array<double,3>& scale={{1,1,1}})
{
	unsigned degree=0;
	for(unsigned q=0;q<3;++q) {
		if(exponent[q]>6) throw std::invalid_argument("polyhedral moment exponent exceeds supported tensor degree six");
		if(!std::isfinite(origin[q])||!std::isfinite(scale[q])||!(scale[q]>0))
			throw std::invalid_argument("polyhedral moment coordinate frame is invalid");
		degree+=exponent[q];
	}
	// x antiderivative raises total degree by one; the Duffy Jacobian
	// contributes another degree in its first coordinate.
	const auto rule=polyhedral_moment_detail::Gauss((degree+4)/2);
	long double sum=0, compensation=0;
	for(const auto& triangle:surface.Triangles()) {
		const auto& a=surface.Vertices()[triangle.indices[0]];
		const auto& b=surface.Vertices()[triangle.indices[1]];
		const auto& c=surface.Vertices()[triangle.indices[2]];
		const long double cross_x=(static_cast<long double>(b[1])-a[1])*(static_cast<long double>(c[2])-a[2])
			-(static_cast<long double>(b[2])-a[2])*(static_cast<long double>(c[1])-a[1]);
		for(const auto& u:rule) for(const auto& v:rule) {
			long double value=cross_x*scale[0]/(exponent[0]+1)*u.second*v.second*(1-u.first);
			for(unsigned q=0;q<3;++q) {
				const long double coordinate=((static_cast<long double>(a[q])-origin[q])
					+u.first*(static_cast<long double>(b[q])-a[q])
					+(1-u.first)*v.first*(static_cast<long double>(c[q])-a[q]))/scale[q];
				value*=polyhedral_moment_detail::Power(coordinate,exponent[q]+(q==0?1:0));
			}
			if(!std::isfinite(value)) throw std::overflow_error("polyhedral moment contribution is not finite");
			const long double next=sum+value;
			compensation+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;
			sum=next;
		}
	}
	const long double result=sum+compensation;
	if(!std::isfinite(result)) throw std::overflow_error("polyhedral moment is not finite");
	return result;
}

}

#endif

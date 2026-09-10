#ifndef IGA_POLYHEDRAL_BOX_VOLUME_MOMENTS_HPP
#define IGA_POLYHEDRAL_BOX_VOLUME_MOMENTS_HPP

#include "PolyhedralVolumeMoments.hpp"

namespace iga {

struct PolyhedralBoxMomentOptions {
	std::size_t max_triangles=1000000;
	std::size_t max_quadrature_points=10000000;
	// Frame in normalized box coordinates. Centering/scaling the polynomial
	// on the retained support avoids ill-conditioned monomials on small cuts.
	std::array<double,3> coordinate_origin{{0,0,0}};
	std::array<double,3> coordinate_scale{{1,1,1}};
};

namespace polyhedral_box_moment_detail {
using Point=std::array<long double,3>;
using Polygon=std::vector<Point>;

inline Polygon Clip(const Polygon& input,unsigned axis,long double plane,bool greater)
{
	Polygon output;
	if(input.empty()) return output;
	auto previous=input.back();
	bool previous_inside=greater?previous[axis]>=plane:previous[axis]<=plane;
	for(const auto& current:input) {
		const bool current_inside=greater?current[axis]>=plane:current[axis]<=plane;
		if(current_inside!=previous_inside) {
			const long double t=(plane-previous[axis])/(current[axis]-previous[axis]);
			Point intersection;
			for(unsigned q=0;q<3;++q) intersection[q]=previous[q]+t*(current[q]-previous[q]);
			intersection[axis]=plane;
			output.push_back(intersection);
		}
		if(current_inside) output.push_back(current);
		previous=current;previous_inside=current_inside;
	}
	return output;
}
}

// Integral on surface-interior intersected with [lower,upper], in physical
// volume units, of product ((xi-coordinate_origin)/coordinate_scale)^exponent,
// where xi=(x-lower)/(upper-lower). The default frame is the unit box. A clipped
// antiderivative along x avoids constructing and triangulating section caps,
// including sections with multiple components or holes. Original facets to
// the RIGHT of the cell must be retained: their saturated antiderivative
// supplies the cap contribution. Cell-AABB triangle queries are insufficient.
// Clipping uses long double arithmetic; this is not an exact-predicate API.
inline long double PolyhedralBoxVolumeMoment(const ClosedTriangulatedSurface& surface,
	const std::array<double,3>& lower,const std::array<double,3>& upper,
	const std::array<unsigned,3>& exponent,PolyhedralBoxMomentOptions options={})
{
	using namespace polyhedral_box_moment_detail;
	if(!options.max_triangles||!options.max_quadrature_points)
		throw std::invalid_argument("polyhedral box moment caps must be positive");
	if(surface.Triangles().size()>options.max_triangles)
		throw std::runtime_error("polyhedral box moment triangle cap reached");
	Point extent;
	unsigned degree=0;
	long double determinant=1;
	for(unsigned q=0;q<3;++q) {
		if(!std::isfinite(lower[q])||!std::isfinite(upper[q])||!(upper[q]>lower[q])||exponent[q]>6
			||!std::isfinite(options.coordinate_origin[q])||!std::isfinite(options.coordinate_scale[q])||!(options.coordinate_scale[q]>0))
			throw std::invalid_argument("invalid polyhedral box moment bounds or exponent");
		extent[q]=static_cast<long double>(upper[q])-lower[q];
		determinant*=extent[q];degree+=exponent[q];
	}
	if(!std::isfinite(determinant)) throw std::overflow_error("polyhedral box moment determinant is not finite");
	const auto rule=polyhedral_moment_detail::Gauss((degree+4)/2);
	std::size_t points=0;
	long double sum=0,compensation=0;
	const auto integrate=[&](const Polygon& polygon,bool saturated) {
		for(std::size_t i=1;i+1<polygon.size();++i) {
			const auto& a=polygon[0];const auto& b=polygon[i];const auto& c=polygon[i+1];
			const long double cross_x=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]);
			if(cross_x==0) continue;
			for(const auto& u:rule) for(const auto& v:rule) {
				if(points==options.max_quadrature_points) throw std::runtime_error("polyhedral box moment quadrature cap reached");
				++points;
				long double value=cross_x*determinant*u.second*v.second*(1-u.first);
				for(unsigned q=0;q<3;++q) {
					const long double coordinate=q==0&&saturated?1:a[q]+u.first*(b[q]-a[q])+(1-u.first)*v.first*(c[q]-a[q]);
					const long double normalized=(coordinate-options.coordinate_origin[q])/options.coordinate_scale[q];
					if(q==0) value*=options.coordinate_scale[0]/(exponent[0]+1.L)
						*(polyhedral_moment_detail::Power(normalized,exponent[0]+1)
						-polyhedral_moment_detail::Power(-static_cast<long double>(options.coordinate_origin[0])/options.coordinate_scale[0],exponent[0]+1));
					else value*=polyhedral_moment_detail::Power(normalized,exponent[q]);
				}
				if(!std::isfinite(value)) throw std::overflow_error("polyhedral box moment contribution is not finite");
				const long double next=sum+value;
				compensation+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;
				sum=next;
			}
		}
	};
	for(const auto& triangle:surface.Triangles()) {
		Polygon polygon;
		for(auto index:triangle.indices) {
			Point point;
			for(unsigned q=0;q<3;++q) {
				point[q]=(static_cast<long double>(surface.Vertices()[index][q])-lower[q])/extent[q];
				if(!std::isfinite(point[q])) throw std::overflow_error("polyhedral box moment vertex is not finite");
			}
			polygon.push_back(point);
		}
		for(unsigned q=1;q<3;++q) { polygon=Clip(polygon,q,0,true);polygon=Clip(polygon,q,1,false); }
		polygon=Clip(polygon,0,0,true);
		if(polygon.size()<3) continue;
		long double minimum=polygon.front()[0],maximum=minimum;
		for(const auto& point:polygon) { minimum=std::min(minimum,point[0]);maximum=std::max(maximum,point[0]); }
		// Assign facets exactly on x=upper to one branch, never both.
		if(minimum>=1) integrate(polygon,true);
		else {
			integrate(Clip(polygon,0,1,false),false);
			if(maximum>1) integrate(Clip(polygon,0,1,true),true);
		}
	}
	const long double result=sum+compensation;
	if(!std::isfinite(result)) throw std::overflow_error("polyhedral box moment is not finite");
	return result;
}

}

#endif

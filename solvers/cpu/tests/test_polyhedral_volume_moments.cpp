#include "PolyhedralVolumeMoments.hpp"
#include "PolyhedralBoxVolumeMoments.hpp"
#include <iostream>

namespace {
void Check(bool value,const char* message)
{
	if(!value) throw std::runtime_error(message);
}

iga::RawSurfaceSoup ConcavePrism()
{
	iga::RawSurfaceSoup result;
	for(double z:{0.,1.}) for(const std::array<double,2>& point:
		{std::array<double,2>{{0,0}},{{2,0}},{{2,1}},{{1,1}},{{1,2}},{{0,2}}})
		result.vertices.push_back({{point[0],point[1],z}});
	for(const std::array<std::int64_t,3>& face:
		{std::array<std::int64_t,3>{{0,1,3}},{{1,2,3}},{{0,3,5}},{{3,4,5}}}) {
		result.triangles.push_back({{{face[2],face[1],face[0]}},0});
		result.triangles.push_back({{{face[0]+6,face[1]+6,face[2]+6}},0});
	}
	for(std::int64_t i=0;i<6;++i) {
		const auto next=(i+1)%6;
		result.triangles.push_back({{{i,next,next+6}},0});
		result.triangles.push_back({{{i,next+6,i+6}},0});
	}
	return result;
}

long double Factorial(unsigned n)
{
	long double value=1;for(unsigned i=2;i<=n;++i)value*=i;return value;
}

long double IntervalMoment(long double lower,long double upper,long double origin,long double width,unsigned degree)
{
	if(upper<=lower) return 0;
	return (iga::polyhedral_moment_detail::Power(upper-origin,degree+1)
		-iga::polyhedral_moment_detail::Power(lower-origin,degree+1))/((degree+1)*iga::polyhedral_moment_detail::Power(width,degree));
}

template<class Function> void Reject(Function function)
{
	bool rejected=false;try{function();}catch(const std::exception&){rejected=true;}
	Check(rejected,"invalid moment input was accepted");
}
}

int main()
{
	try {
		auto soup=ConcavePrism();
		const auto concave=iga::ClosedTriangulatedSurface::Build(soup);
		const std::array<double,3> origin{{1048576,-2097152,4194304}},scale{{.125,2,.5}};
		for(auto& point:soup.vertices)for(unsigned q=0;q<3;++q)point[q]=origin[q]+scale[q]*point[q];
		// Input winding and triangle order must not change the validated integral.
		for(auto& triangle:soup.triangles)std::swap(triangle.indices[0],triangle.indices[1]);
		std::reverse(soup.triangles.begin(),soup.triangles.end());
		const auto shifted=iga::ClosedTriangulatedSurface::Build(soup);
		iga::RawSurfaceSoup tetra;
		tetra.vertices={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}};
		tetra.triangles={{{{0,2,1}},0},{{{0,1,3}},0},{{{0,3,2}},0},{{{1,2,3}},0}};
		const auto simplex=iga::ClosedTriangulatedSurface::Build(tetra);
		long double maximum_relative=0;
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			const std::array<unsigned,3> exponent{{a,b,c}};
			const long double x=std::ldexp(1.L,a+1),y=std::ldexp(1.L,b+1);
			const long double expected=(x*y-(x-1)*(y-1))/((a+1)*(b+1)*(c+1));
			const long double actual=iga::PolyhedralVolumeMoment(concave,exponent);
			const long double translated=iga::PolyhedralVolumeMoment(shifted,exponent,origin,scale);
			const long double simplex_expected=Factorial(a)*Factorial(b)*Factorial(c)/Factorial(a+b+c+3);
			const long double simplex_actual=iga::PolyhedralVolumeMoment(simplex,exponent);
			for(long double error:{std::abs(actual/expected-1),std::abs(translated/(expected*.125L)-1),std::abs(simplex_actual/simplex_expected-1)}) {
				maximum_relative=std::max(maximum_relative,error);
				Check(error<2e-13L,"analytic polyhedral moment mismatch");
			}
		}
		Reject([&]{iga::PolyhedralVolumeMoment(concave,{{7,0,0}});});
		Reject([&]{iga::PolyhedralVolumeMoment(concave,{{0,0,0}},{{0,0,0}},{{1,0,1}});});
		Reject([&]{iga::PolyhedralVolumeMoment(concave,{{0,0,0}},{{std::numeric_limits<double>::infinity(),0,0}});});
		long double maximum_box_relative=0;
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			const std::array<unsigned,3> exponent{{a,b,c}};
			const long double expected=(IntervalMoment(.25L,1.75L,.25L,1.5L,a)*IntervalMoment(.5L,1,.5L,1,b)
				+IntervalMoment(.25L,1,.25L,1.5L,a)*IntervalMoment(1,1.5L,.5L,1,b))*.75L/(c+1);
			const long double actual=iga::PolyhedralBoxVolumeMoment(concave,{{.25,.5,.125}},{{1.75,1.5,.875}},exponent);
			// Analytic simplex x-slab: integrate x^a (1-x)^(b+c+2)
			// using its binomial expansion; y and z integrals are beta moments.
			long double slab=0;
			const unsigned n=b+c+2;
			for(unsigned k=0;k<=n;++k) slab+=(k%2?-1.L:1.L)*Factorial(n)/(Factorial(k)*Factorial(n-k))
				*iga::polyhedral_moment_detail::Power(.25L,a+k+1)/(a+k+1);
			slab*=Factorial(b)*Factorial(c)/Factorial(n)/iga::polyhedral_moment_detail::Power(.25L,a);
			const long double clipped=iga::PolyhedralBoxVolumeMoment(simplex,{{0,0,0}},{{.25,1,1}},exponent);
			for(long double error:{std::abs(actual/expected-1),std::abs(clipped/slab-1)}) {
				maximum_box_relative=std::max(maximum_box_relative,error);
				Check(error<2e-13L,"analytic clipped polyhedral moment mismatch");
			}
		}
		long double partition_volume=0,partition_x=0;
		for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)for(unsigned k=0;k<2;++k) {
			const std::array<double,3> lo{{.5*i,.5*j,.5*k}},hi{{.5*(i+1),.5*(j+1),.5*(k+1)}};
			const auto volume=iga::PolyhedralBoxVolumeMoment(concave,lo,hi,{{0,0,0}});
			partition_volume+=volume;
			partition_x+=lo[0]*volume+.5L*iga::PolyhedralBoxVolumeMoment(concave,lo,hi,{{1,0,0}});
		}
		Check(std::abs(partition_volume-3)<2e-14L&&std::abs(partition_x-2.5L)<2e-14L,"box partition moments do not reconcile");
		Check(std::abs(iga::PolyhedralBoxVolumeMoment(concave,{{3,0,0}},{{4,1,1}},{{0,0,0}}))<2e-14L,"outside box is not empty");
		Check(std::abs(iga::PolyhedralBoxVolumeMoment(concave,{{.25,.25,.25}},{{.5,.5,.5}},{{0,0,0}})-.015625L)<2e-14L,"interior box lost right-hand cap contribution");
		Reject([&]{iga::PolyhedralBoxVolumeMoment(concave,{{0,0,0}},{{0,1,1}},{{0,0,0}});});
		Reject([&]{iga::PolyhedralBoxVolumeMoment(concave,{{0,0,0}},{{1,1,1}},{{0,0,0}},{1,10000});});
		Reject([&]{iga::PolyhedralBoxVolumeMoment(concave,{{0,0,0}},{{1,1,1}},{{0,0,0}},{10000,1});});
		std::cout<<"polyhedral_volume_moments_test: PASS 1029 analytic moments; max_relative="<<maximum_relative<<'\n';
		std::cout<<"polyhedral_box_volume_moments: PASS 686 analytic moments and partition/cap gates; max_relative="<<maximum_box_relative<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}

#include "PolyhedralVolumeMoments.hpp"
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
		std::cout<<"polyhedral_volume_moments_test: PASS 1029 analytic moments; max_relative="<<maximum_relative<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}

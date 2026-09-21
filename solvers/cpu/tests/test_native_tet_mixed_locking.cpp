#include "NativeTetMixedSolidStaticSolver.hpp"
#include "NativeTetAleKinematics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

namespace {
constexpr double length=4.,height=1.,width=1.;
std::uint32_t Node(int nx,int ny,int nz,int i,int j,int k)
{(void)nx;return static_cast<std::uint32_t>((i*(ny+1)+j)*(nz+1)+k);}
iga::NativeTetMesh Mesh(int nx,int ny,int nz)
{
	iga::NativeTetMesh mesh;
	for(int i=0;i<=nx;++i)for(int j=0;j<=ny;++j)for(int k=0;k<=nz;++k)
		mesh.points.push_back({{length*i/nx,height*j/ny,width*k/nz}});
	std::uint64_t id=1;
	for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k){
		const auto a=Node(nx,ny,nz,i,j,k),b=Node(nx,ny,nz,i+1,j,k);
		const auto c=Node(nx,ny,nz,i,j+1,k),d=Node(nx,ny,nz,i+1,j+1,k);
		const auto e=Node(nx,ny,nz,i,j,k+1),f=Node(nx,ny,nz,i+1,j,k+1);
		const auto g=Node(nx,ny,nz,i,j+1,k+1),h=Node(nx,ny,nz,i+1,j+1,k+1);
		const std::array<std::array<std::uint32_t,4>,6> cells{{{{a,b,d,h}},{{a,d,c,h}},{{a,c,g,h}},
			{{a,g,e,h}},{{a,e,f,h}},{{a,f,b,h}}}};
		for(auto nodes:cells){if(iga::NativeAleDeterminant(mesh.points[nodes[0]],mesh.points[nodes[1]],
			mesh.points[nodes[2]],mesh.points[nodes[3]])<0.)std::swap(nodes[1],nodes[2]);
			mesh.cells.push_back({id++,nodes});}
	}
	return mesh;
}
struct Observation{double displacement=0.,volume_change=0.;};
Observation Run(int transverse,double poisson=.4999)
{
	const int nx=2*transverse,ny=transverse,nz=transverse;const auto mesh=Mesh(nx,ny,nz);
	std::map<std::size_t,double> constraints;std::vector<std::size_t> tip;
	for(std::size_t node=0;node<mesh.points.size();++node){
		if(mesh.points[node][0]==0.)for(int component=0;component<3;++component)constraints[3*node+component]=0.;
		if(mesh.points[node][0]==length)tip.push_back(node);
	}
	std::vector<double> load(3*mesh.points.size(),0.);for(auto node:tip)load[3*node+1]=-.01/tip.size();
	const auto material=iga::NativeTetMixedSolidMaterialFromYoungPoisson(1000.,poisson,1000.,.1);
	iga::NativeTetSolidStaticOptions options;options.maximum_dofs=3000;options.maximum_iterations=20;
	options.relative_tolerance=1e-8;options.absolute_tolerance_n=1e-10;
	const auto result=iga::SolveNativeTetMixedSolidStatic(mesh,material,load,constraints,{}, {},options);
	double mean=0.;for(auto node:tip)mean+=result.displacement_m[3*node+1]/tip.size();
	assert(result.converged&&result.minimum_deformation_jacobian>.9&&mean<0.);
	return {mean,result.maximum_absolute_volume_change};
}
}
int main()
{
	const std::array<int,6> refinement{{1,2,3,4,5,6}};std::array<Observation,6> observation{};
	for(std::size_t level=0;level<refinement.size();++level){std::cout<<"mixed locking level="<<refinement[level]<<std::endl;
		observation[level]=Run(refinement[level]);}
	for(std::size_t level=1;level<observation.size();++level)
		assert(std::abs(observation[level].displacement)>std::abs(observation[level-1].displacement));
	assert(std::abs(observation.back().displacement)>.002
		&&std::abs(observation[5].displacement-observation[4].displacement)
			/std::abs(observation[5].displacement)<.06);
	const std::array<double,3> poisson{{.49,.499,.4999}};std::array<Observation,3> poisson_observation{};
	for(std::size_t i=0;i<poisson.size();++i)poisson_observation[i]=Run(3,poisson[i]);
	assert(std::abs(poisson_observation[2].displacement-poisson_observation[0].displacement)
		/std::abs(poisson_observation[2].displacement)<.005);
	assert(poisson_observation[2].volume_change<2e-4);
	std::cout<<"native mixed near-incompressible bending displacements="<<observation[0].displacement<<','
		<<observation[1].displacement<<','<<observation[2].displacement<<','<<observation[3].displacement<<','
		<<observation[4].displacement<<','<<observation[5].displacement<<" poisson_sweep="
		<<poisson_observation[0].displacement<<','<<poisson_observation[1].displacement<<','
		<<poisson_observation[2].displacement<<" max_volume_change="
		<<poisson_observation[2].volume_change<<'\n';return 0;
}

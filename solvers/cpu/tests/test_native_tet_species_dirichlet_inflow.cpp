#include "NativeTetMovingSpeciesPetscRuntime.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

// Channel [0, L] x [0, 0.2] x [0, 0.2] of `segments` cubes, six tetrahedra each.
// Labels: x = 0 inlet (1), x = L outlet (2), lateral walls (0).
iga::NativeTetMesh Channel(int segments,double length)
{
	iga::NativeTetMesh mesh;
	const double h=length/segments,width=0.2;
	const auto id=[](int i,int j,int k){return static_cast<std::uint32_t>(4*i+2*j+k);};
	for(int i=0;i<=segments;++i)
		for(int j=0;j<2;++j)
			for(int k=0;k<2;++k)mesh.points.push_back({{i*h,j*width,k*width}});
	std::uint64_t cell_id=1;
	for(int i=0;i<segments;++i){
		const std::uint32_t v0=id(i,0,0),v1=id(i+1,0,0),v2=id(i,1,0),v3=id(i+1,1,0),
			v4=id(i,0,1),v5=id(i+1,0,1),v6=id(i,1,1),v7=id(i+1,1,1);
		const std::array<std::array<std::uint32_t,4>,6> tets{{
			{{v0,v1,v3,v7}},{{v0,v1,v5,v7}},{{v0,v2,v3,v7}},
			{{v0,v2,v6,v7}},{{v0,v4,v5,v7}},{{v0,v4,v6,v7}}}};
		for(auto tet:tets){
			const auto& a=mesh.points[tet[0]];
			std::array<Vector,3> e{};
			for(int r=0;r<3;++r)for(int c=0;c<3;++c)e[r][c]=mesh.points[tet[r+1]][c]-a[c];
			const double det=e[0][0]*(e[1][1]*e[2][2]-e[1][2]*e[2][1])
				-e[0][1]*(e[1][0]*e[2][2]-e[1][2]*e[2][0])+e[0][2]*(e[1][0]*e[2][1]-e[1][1]*e[2][0]);
			if(det<0.)std::swap(tet[2],tet[3]);
			mesh.cells.push_back({cell_id++,tet});
		}
	}
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,unsigned> uses;
	for(const auto& cell:mesh.cells)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[entry++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++uses[face];
		}
	std::uint64_t face_id=1;
	for(const auto& item:uses){
		if(item.second!=1)continue;
		int label=0;
		if(std::all_of(item.first.begin(),item.first.end(),
			[&](std::uint32_t n){return mesh.points[n][0]<1e-12;}))label=1;
		else if(std::all_of(item.first.begin(),item.first.end(),
			[&](std::uint32_t n){return mesh.points[n][0]>length-1e-12;}))label=2;
		mesh.boundary_triangles.push_back({face_id++,item.first,label});
	}
	return mesh;
}

std::vector<Vector> Uniform(const iga::NativeTetMesh& mesh,double speed)
{
	return std::vector<Vector>(mesh.points.size()
		+iga::BuildNativeTaylorHoodTopology(mesh).edges.size(),Vector{{speed,0,0}});
}

bool Rejected(const std::function<void()>& call)
{
	try{call();}catch(const std::invalid_argument&){return true;}
	return false;
}

}

int main(int argc,char** argv)
{
	if(PetscInitialize(&argc,&argv,nullptr,nullptr))return 1;
	int rank=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	{
		const double length=1.,dt=0.1,diffusivity=0.05,area=0.04;
		const auto mesh=Channel(10,length);
		const std::vector<Vector> grid(mesh.points.size(),Vector{{0,0,0}});
		const auto velocity=Uniform(mesh,1.);
		const std::map<int,double> inflow{{1,1.}};
		iga::NativeTetMovingSpeciesStepOptions strong;
		strong.dirichlet_labels={1};
		std::vector<double> older,previous(mesh.points.size(),0.);
		iga::NativeTetMovingSpeciesStep last;
		for(int step=0;step<80;++step){
			auto options=strong;
			if(step>0){
				options.time_integration=iga::NativeTetSpeciesTimeIntegration::Bdf2;
				options.older_concentration_mol_m3=older;
				options.allow_negative_concentration=true;
			}
			const auto dense=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,velocity,grid,
				previous,inflow,diffusivity,0.,dt,false,0.,{},{},{},options);
			const auto parallel=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,mesh,velocity,grid,
				previous,inflow,diffusivity,0.,dt,false,0.,{},{},{},{},options);
			assert(std::abs(dense.balance_defect_mol_s)<=1e-11);
			assert(std::abs(parallel.step.balance_defect_mol_s)<=1e-9);
			for(std::size_t n=0;n<mesh.points.size();++n){
				assert(std::abs(dense.concentration_mol_m3[n]-parallel.step.concentration_mol_m3[n])<=1e-9);
				if(mesh.points[n][0]<1e-12)assert(dense.concentration_mol_m3[n]==1.);
			}
			assert(std::abs(dense.dirichlet_reaction_by_label_mol_s.at(1)
				-parallel.step.dirichlet_reaction_by_label_mol_s.at(1))<=1e-9);
			older=previous;
			previous=dense.concentration_mol_m3;
			last=dense;
		}
		// Steady state c = 1: inlet transfer is -u A c_in (advection only), outlet +u A.
		for(const double value:previous)assert(std::abs(value-1.)<1e-6);
		assert(std::abs(last.outward_advective_flux_by_label_mol_s.at(1)+area)<1e-6);
		assert(std::abs(last.outward_advective_flux_by_label_mol_s.at(2)-area)<1e-6);
		if(rank==0)
			std::cout<<"inlet transfer="<<last.outward_advective_flux_by_label_mol_s.at(1)
				<<" outlet transfer="<<last.outward_advective_flux_by_label_mol_s.at(2)
				<<" reaction="<<last.dirichlet_reaction_by_label_mol_s.at(1)<<'\n';

		// Weak inflow stays the default and differs from the strong condition while the front passes.
		const std::vector<double> empty(mesh.points.size(),0.);
		const auto weak=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,velocity,grid,
			empty,inflow,diffusivity,0.,dt);
		assert(weak.dirichlet_reaction_by_label_mol_s.empty());
		assert(std::abs(weak.balance_defect_mol_s)<=1e-12);
		const auto strong_first=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,velocity,grid,
			empty,inflow,diffusivity,0.,dt,false,0.,{},{},{},strong);
		assert(std::abs(weak.concentration_mol_m3.front()-strong_first.concentration_mol_m3.front())>1e-3);

		// Invalid uses are rejected.
		assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,velocity,grid,
			empty,{},diffusivity,0.,dt,false,0.,{},{},{},strong);}));
		iga::NativeTetMovingSpeciesStepOptions missing;
		missing.dirichlet_labels={7};
		assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,velocity,grid,
			empty,{{1,1.},{7,1.}},diffusivity,0.,dt,false,0.,{},{},{},missing);}));
		const std::vector<Vector> still(velocity.size(),Vector{{0,0,0}});
		assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,still,grid,
			empty,{{0,1.}},diffusivity,0.,dt,false,0.,{{0,{1e-3,0.}}},{},{},
			[]{iga::NativeTetMovingSpeciesStepOptions o;o.dirichlet_labels={0};return o;}());}));
		if(rank==0)std::cout<<"native Tet species Dirichlet inflow: PASS\n";
	}
	PetscFinalize();
	return 0;
}

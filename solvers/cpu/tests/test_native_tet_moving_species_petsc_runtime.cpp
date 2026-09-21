#include "NativeTetMovingSpeciesPetscRuntime.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh SmallMesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0,0,-1}}};
	mesh.cells={{1,{{0,1,2,3}}},{2,{{0,2,1,4}}}};
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,unsigned> uses;
	for(const auto& cell:mesh.cells)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[entry++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++uses[face];
		}
	std::uint64_t id=1;
	for(const auto& item:uses)
		if(item.second==1)mesh.boundary_triangles.push_back({id++,item.first,1});
	return mesh;
}

iga::NativeTetMesh ManyCells()
{
	iga::NativeTetMesh mesh;
	for(std::uint32_t cell=0;cell<40;++cell){
		const auto base=static_cast<std::uint32_t>(mesh.points.size());
		const double x=3.*cell;
		mesh.points.push_back({{x,0,0}});
		mesh.points.push_back({{x+1,0,0}});
		mesh.points.push_back({{x,1,0}});
		mesh.points.push_back({{x,0,1}});
		mesh.cells.push_back({cell+1,{{base,base+1,base+2,base+3}}});
		mesh.boundary_triangles.push_back({4*cell+1,{{base+1,base+2,base+3}},1});
		mesh.boundary_triangles.push_back({4*cell+2,{{base,base+2,base+3}},1});
		mesh.boundary_triangles.push_back({4*cell+3,{{base,base+1,base+3}},1});
		mesh.boundary_triangles.push_back({4*cell+4,{{base,base+1,base+2}},1});
	}
	return mesh;
}

std::vector<Vector> P2Velocity(const iga::NativeTetMesh& mesh,const Vector& velocity)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	return std::vector<Vector>(mesh.points.size()+topology.edges.size(),velocity);
}

std::vector<Vector> P2VertexField(const iga::NativeTetMesh& mesh,
	const std::vector<Vector>& vertices)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	std::vector<Vector> result=vertices;
	for(const auto& edge:topology.edges){
		Vector value{};
		for(int component=0;component<3;++component)
			value[component]=0.5*(vertices[edge[0]][component]
				+vertices[edge[1]][component]);
		result.push_back(value);
	}
	return result;
}

void Close(double actual,double expected,double tolerance=1e-10)
{
	assert(std::isfinite(actual));
	assert(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)));
}

}

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const auto previous=SmallMesh();
	auto translated=previous;
	for(auto& point:translated.points)point[0]+=0.1;
	const std::vector<Vector> grid(previous.points.size(),Vector{{1,0,0}});
	const std::vector<double> empty(previous.points.size(),0.);
	const auto high_peclet_dense=iga::SolveNativeTetMovingSpeciesDenseStep(
		previous,previous,P2Velocity(previous,Vector{{1,0,0}}),
		std::vector<Vector>(previous.points.size(),Vector{{0,0,0}}),
		empty,{{1,1.}},0.,0.,0.01,true);
	const auto high_peclet_parallel=iga::SolveNativeTetMovingSpeciesPetscStep(
		previous,previous,P2Velocity(previous,Vector{{1,0,0}}),
		std::vector<Vector>(previous.points.size(),Vector{{0,0,0}}),
		empty,{{1,1.}},0.,0.,0.01,true);
	for(std::size_t node=0;node<empty.size();++node){
		assert(high_peclet_parallel.step.concentration_mol_m3[node]>=-1e-12);
		Close(high_peclet_parallel.step.concentration_mol_m3[node],
			high_peclet_dense.concentration_mol_m3[node]);
	}
	Close(high_peclet_parallel.step.balance_defect_mol_s,0.);
	const auto velocity=P2Velocity(translated,Vector{{2,0,0}});
	const auto dense=iga::SolveNativeTetMovingSpeciesDenseStep(previous,translated,
		velocity,grid,empty,{{1,2.}},0.01,0.,0.1);
	const auto parallel=iga::SolveNativeTetMovingSpeciesPetscStep(previous,translated,
		velocity,grid,empty,{{1,2.}},0.01,0.,0.1);
	assert(parallel.converged_reason>0);
	for(std::size_t node=0;node<empty.size();++node)
		Close(parallel.step.concentration_mol_m3[node],
			dense.concentration_mol_m3[node]);
	Close(parallel.step.current_inventory_mol,dense.current_inventory_mol);
	Close(parallel.step.outward_advective_flux_mol_s,
		dense.outward_advective_flux_mol_s);
	Close(parallel.step.balance_defect_mol_s,0.);
	const std::vector<Vector> fixed_grid(previous.points.size(),Vector{{0,0,0}});
	const std::vector<double> wall_uniform(previous.points.size(),2.);
	const std::map<int,iga::NativeTetWallExchange> wall{{1,{0.5,0.}}};
	const auto dense_wall=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		P2Velocity(previous,Vector{{0,0,0}}),fixed_grid,wall_uniform,{},
		0.,0.,0.1,false,0.,wall);
	const auto parallel_wall=iga::SolveNativeTetMovingSpeciesPetscStep(previous,
		previous,P2Velocity(previous,Vector{{0,0,0}}),fixed_grid,wall_uniform,{},
		0.,0.,0.1,false,0.,wall);
	for(std::size_t node=0;node<wall_uniform.size();++node)
		Close(parallel_wall.step.concentration_mol_m3[node],
			dense_wall.concentration_mol_m3[node]);
	Close(parallel_wall.step.outward_wall_exchange_mol_s,
		dense_wall.outward_wall_exchange_mol_s);
	Close(parallel_wall.step.balance_defect_mol_s,0.);
	auto labelled=previous;
	for(auto& face:labelled.boundary_triangles){
		const bool inlet=std::all_of(face.nodes.begin(),face.nodes.end(),
			[&](std::uint32_t node){return labelled.points[node][0]==0.;});
		face.boundary_label=inlet?2:3;
	}
	const std::vector<Vector> stationary_grid(labelled.points.size(),Vector{{0,0,0}});
	const std::vector<double> uniform_ports(labelled.points.size(),2.);
	const auto dense_ports=iga::SolveNativeTetMovingSpeciesDenseStep(labelled,labelled,
		P2Velocity(labelled,Vector{{1,0,0}}),stationary_grid,uniform_ports,
		{{2,2.}},0.,0.,0.1);
	const auto parallel_ports=iga::SolveNativeTetMovingSpeciesPetscStep(labelled,labelled,
		P2Velocity(labelled,Vector{{1,0,0}}),stationary_grid,uniform_ports,
		{{2,2.}},0.,0.,0.1);
	Close(parallel_ports.step.outward_advective_flux_by_label_mol_s.at(2),-2.);
	Close(parallel_ports.step.outward_advective_flux_by_label_mol_s.at(3),2.);
	Close(parallel_ports.step.outward_relative_flow_by_label_m3_s.at(2),-1.);
	Close(parallel_ports.step.outward_relative_flow_by_label_m3_s.at(3),1.);
	for(const auto& item:dense_ports.outward_advective_flux_by_label_mol_s)
		Close(parallel_ports.step.outward_advective_flux_by_label_mol_s.at(item.first),
			item.second);
	for(const auto& item:dense_ports.outward_relative_flow_by_label_m3_s)
		Close(parallel_ports.step.outward_relative_flow_by_label_m3_s.at(item.first),
			item.second);
	Close(parallel_ports.step.outward_advective_flux_by_label_mol_s.at(2)
		+parallel_ports.step.outward_advective_flux_by_label_mol_s.at(3),
		parallel_ports.step.outward_advective_flux_mol_s);
	bool rejected=false;
	try{iga::SolveNativeTetMovingSpeciesPetscStep(previous,translated,
		velocity,grid,empty,{},0.01,0.,0.1);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	const std::vector<double> uniform(previous.points.size(),2.);
	const auto first=iga::SolveNativeTetMovingSpeciesPetscStep(previous,translated,
		P2Velocity(translated,Vector{{1,0,0}}),grid,uniform,{},0.,0.,0.1);
	auto expanded=translated;
	for(auto& point:expanded.points)point[2]*=1.2;
	std::vector<Vector> expansion(expanded.points.size());
	for(std::size_t node=0;node<expansion.size();++node)
		expansion[node][2]=(expanded.points[node][2]
			-translated.points[node][2])/0.1;
	const auto second=iga::SolveNativeTetMovingSpeciesPetscStep(translated,expanded,
		P2VertexField(expanded,expansion),expansion,
		first.step.concentration_mol_m3,{},0.,0.,0.1);
	for(const double value:second.step.concentration_mol_m3)Close(value,2./1.2);
	Close(first.step.current_inventory_mol,second.step.current_inventory_mol);
	Close(second.step.balance_defect_mol_s,0.);
	const auto large=ManyCells();
	assert(large.points.size()>128);
	const std::vector<Vector> no_motion(large.points.size(),Vector{{0,0,0}});
	const std::vector<double> constant(large.points.size(),2.);
	const auto scaled=iga::SolveNativeTetMovingSpeciesPetscStep(large,large,
		P2Velocity(large,Vector{{0,0,0}}),no_motion,constant,{},0.,0.,0.1);
	assert(scaled.converged_reason>0);
	for(const double value:scaled.step.concentration_mol_m3)Close(value,2.);
	Close(scaled.step.balance_defect_mol_s,0.);
	const auto scaled_source=iga::SolveNativeTetMovingSpeciesPetscStep(large,large,
		P2Velocity(large,Vector{{0,0,0}}),no_motion,constant,{},0.,3.,0.1);
	for(const double value:scaled_source.step.concentration_mol_m3)Close(value,2.3);
	Close(scaled_source.step.source_mol_s,20.);
	Close(scaled_source.step.balance_defect_mol_s,0.);
	const auto scaled_decay=iga::SolveNativeTetMovingSpeciesPetscStep(large,large,
		P2Velocity(large,Vector{{0,0,0}}),no_motion,constant,{},0.,3.,0.1,false,2.);
	for(const double value:scaled_decay.step.concentration_mol_m3)
		Close(value,(2.+0.1*3.)/(1.+0.1*2.));
	Close(scaled_decay.step.reaction_sink_mol_s,
		2.*scaled_decay.step.current_inventory_mol);
	Close(scaled_decay.step.balance_defect_mol_s,0.);
	if(rank==0)std::cout<<"native moving species PETSc passed ranks="<<ranks
		<<" large_vertices="<<large.points.size()
		<<" small_iterations="<<parallel.linear_iterations
		<<" large_iterations="<<scaled.linear_iterations
		<<" balance_mol_s="<<parallel.step.balance_defect_mol_s<<'\n';
	PetscFinalize();
}

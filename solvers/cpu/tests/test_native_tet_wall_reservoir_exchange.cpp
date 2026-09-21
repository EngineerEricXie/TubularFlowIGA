#include "NativeTetWallReservoirExchange.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh TwoCellMesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0,0,-1}}};
	mesh.cells={{1,{{0,1,2,3}}},{2,{{0,2,1,4}}}};
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,unsigned> uses;
	for(const auto& cell:mesh.cells)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t index=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[index++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++uses[face];
		}
	std::uint64_t id=1;
	for(const auto& item:uses)
		if(item.second==1)mesh.boundary_triangles.push_back({id++,item.first,1});
	return mesh;
}

std::vector<Vector> Velocity(const iga::NativeTetMesh& mesh,const Vector& value)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	return std::vector<Vector>(mesh.points.size()+topology.edges.size(),value);
}

void Close(double actual,double expected,double tolerance=1e-10)
{
	assert(std::isfinite(actual));
	assert(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)));
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const auto mesh=TwoCellMesh();
	const std::vector<double> blood(mesh.points.size(),2.);
	const std::vector<Vector> fixed(mesh.points.size(),Vector{{0,0,0}});
	const iga::NativeTetWallReservoirModel model{1,0.5,1.};
	const iga::NativeTetWallReservoirState empty{0.};
	const auto uptake=iga::SolveNativeTetWallReservoirStep(mesh,mesh,
		Velocity(mesh,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
		model,empty,true);
	assert(uptake.outward_wall_flux_mol_s>0.);
	assert(uptake.reservoir.amount_mol>0.);
	assert(uptake.vessel.step.current_inventory_mol
		<uptake.vessel.step.previous_inventory_mol);
	Close(uptake.reservoir.amount_mol-empty.amount_mol,
		0.1*uptake.outward_wall_flux_mol_s);
	Close(uptake.vessel.step.current_inventory_mol
		+uptake.reservoir.amount_mol,
		uptake.vessel.step.previous_inventory_mol+empty.amount_mol);
	Close(uptake.vessel.step.balance_defect_mol_s,0.);
	Close(uptake.reservoir_balance_defect_mol_s,0.);
	Close(uptake.combined_balance_defect_mol_s,0.);
	for(const double value:uptake.vessel.step.concentration_mol_m3)
		assert(value>=0.&&value<=2.+1e-10);
	const iga::NativeTetWallReservoirState donor{3.};
	const auto release=iga::SolveNativeTetWallReservoirStep(mesh,mesh,
		Velocity(mesh,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
		model,donor,true);
	assert(release.outward_wall_flux_mol_s<0.);
	assert(release.reservoir.amount_mol<donor.amount_mol);
	assert(release.vessel.step.current_inventory_mol
		>release.vessel.step.previous_inventory_mol);
	Close(release.vessel.step.current_inventory_mol
		+release.reservoir.amount_mol,
		release.vessel.step.previous_inventory_mol+donor.amount_mol);
	Close(release.reservoir_balance_defect_mol_s,0.);
	Close(release.combined_balance_defect_mol_s,0.);
	auto translated=mesh;
	for(auto& point:translated.points)point[0]+=0.1;
	const std::vector<Vector> moving(mesh.points.size(),Vector{{1,0,0}});
	const auto moved=iga::SolveNativeTetWallReservoirStep(mesh,translated,
		Velocity(translated,Vector{{1,0,0}}),moving,blood,{},0.,0.,0.1,
		model,empty,true);
	for(std::size_t node=0;node<blood.size();++node)
		Close(moved.vessel.step.concentration_mol_m3[node],
			uptake.vessel.step.concentration_mol_m3[node]);
	Close(moved.reservoir.amount_mol,uptake.reservoir.amount_mol);
	Close(moved.outward_wall_flux_mol_s,uptake.outward_wall_flux_mol_s);
	Close(moved.combined_balance_defect_mol_s,0.);
	const auto material=iga::SolveNativeTetWallReservoirStep(mesh,mesh,
		Velocity(mesh,Vector{{0,0,0}}),fixed,blood,{},0.,0.5,0.1,
		model,empty,true,0.2);
	Close(material.combined_balance_defect_mol_s,0.);
	Close(material.vessel.step.balance_defect_mol_s,0.);
	auto two_regions=mesh;
	for(auto& face:two_regions.boundary_triangles)
		face.boundary_label=std::find(face.nodes.begin(),face.nodes.end(),4)
			!=face.nodes.end()?2:1;
	const std::map<int,iga::NativeTetWallReservoirRegion> regions{
		{1,{{1,0.5,1.},{0.}}},
		{2,{{2,0.5,1.},{3.}}}};
	const auto paired=iga::SolveNativeTetWallReservoirsStep(two_regions,two_regions,
		Velocity(two_regions,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
		regions,true);
	assert(paired.outward_wall_flux_mol_s.at(1)>0.);
	assert(paired.outward_wall_flux_mol_s.at(2)<0.);
	assert(paired.reservoirs.at(1).amount_mol>0.);
	assert(paired.reservoirs.at(2).amount_mol<3.);
	Close(paired.reservoirs.at(1).amount_mol,
		0.1*paired.outward_wall_flux_mol_s.at(1));
	Close(paired.reservoirs.at(2).amount_mol-3.,
		0.1*paired.outward_wall_flux_mol_s.at(2));
	Close(paired.vessel.step.current_inventory_mol
		+paired.reservoirs.at(1).amount_mol+paired.reservoirs.at(2).amount_mol,
		paired.vessel.step.previous_inventory_mol+3.);
	Close(paired.combined_balance_defect_mol_s,0.);
	const auto both_empty=iga::SolveNativeTetWallReservoirsStep(
		two_regions,two_regions,Velocity(two_regions,Vector{{0,0,0}}),fixed,
		blood,{},0.,0.,0.1,
		{{1,{{1,0.5,1.},{0.}}},{2,{{2,0.5,1.},{0.}}}},true);
	assert(paired.outward_wall_flux_mol_s.at(1)
		>both_empty.outward_wall_flux_mol_s.at(1)+1e-8);
	const auto paired_again=iga::SolveNativeTetWallReservoirsStep(
		two_regions,two_regions,Velocity(two_regions,Vector{{0,0,0}}),fixed,
		paired.vessel.step.concentration_mol_m3,{},0.,0.,0.1,
		{{1,{{1,0.5,1.},paired.reservoirs.at(1)}},
			{2,{{2,0.5,1.},paired.reservoirs.at(2)}}},true);
	Close(paired_again.combined_balance_defect_mol_s,0.);
	Close(paired.reservoirs.at(1).amount_mol,
		0.1*paired.outward_wall_flux_mol_s.at(1));
	Close(regions.at(1).committed.amount_mol,0.);
	Close(regions.at(2).committed.amount_mol,3.);
	bool region_rejected=false;
	try{iga::SolveNativeTetWallReservoirsStep(two_regions,two_regions,
			Velocity(two_regions,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
			{{1,{{2,0.5,1.},{0.}}}},true);}
	catch(const std::invalid_argument&){region_rejected=true;}
	assert(region_rejected);
	bool rejected=false;
	try{iga::SolveNativeTetWallReservoirStep(mesh,mesh,
			Velocity(mesh,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
			{99,0.5,1.},empty,true);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetWallReservoirStep(mesh,mesh,
			Velocity(mesh,Vector{{0,0,0}}),fixed,blood,{},0.,0.,0.1,
			model,{-1.},true);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetWallReservoirStep(mesh,mesh,
			Velocity(mesh,Vector{{1,0,0}}),fixed,blood,{},0.,0.,0.1,
			model,empty,true);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	Close(empty.amount_mol,0.);
	if(argc>1){
		std::ifstream input(argv[1]);
		if(!input)throw std::runtime_error("native wall reservoir mesh is missing");
		const auto pipe=iga::ReadNativeTetMeshGmsh41(input);
		const std::vector<double> pipe_blood(pipe.points.size(),2.);
		const std::vector<Vector> pipe_grid(pipe.points.size(),Vector{{0,0,0}});
		const iga::NativeTetWallReservoirModel pipe_model{0,1e-4,1e-6};
		const auto first=iga::SolveNativeTetWallReservoirStep(pipe,pipe,
			Velocity(pipe,Vector{{0,0,0}}),pipe_grid,pipe_blood,{},
			1e-5,0.,0.01,pipe_model,empty,true);
		assert(first.outward_wall_flux_mol_s>0.);
		Close(first.combined_balance_defect_mol_s,0.);
		const auto second=iga::SolveNativeTetWallReservoirStep(pipe,pipe,
			Velocity(pipe,Vector{{0,0,0}}),pipe_grid,
			first.vessel.step.concentration_mol_m3,{},1e-5,0.,0.01,
			pipe_model,first.reservoir,true);
		assert(second.reservoir.amount_mol>first.reservoir.amount_mol);
		Close(second.combined_balance_defect_mol_s,0.);
		const auto stronger=iga::SolveNativeTetWallReservoirStep(pipe,pipe,
			Velocity(pipe,Vector{{0,0,0}}),pipe_grid,pipe_blood,{},
			1e-5,0.,0.01,{0,2e-4,1e-6},empty,true);
		assert(stronger.reservoir.amount_mol>first.reservoir.amount_mol+1e-13);
		Close(stronger.combined_balance_defect_mol_s,0.);
		const auto reverse=iga::SolveNativeTetWallReservoirStep(pipe,pipe,
			Velocity(pipe,Vector{{0,0,0}}),pipe_grid,pipe_blood,{},
			1e-5,0.,0.01,pipe_model,{3e-6},true);
		assert(reverse.outward_wall_flux_mol_s<0.);
		Close(reverse.combined_balance_defect_mol_s,0.);
		if(rank==0)std::cout<<"native wall reservoir fTetWild pipe: PASS tetrahedra="
			<<pipe.cells.size()<<" uptake_mol="<<first.reservoir.amount_mol
			<<" stronger_uptake_mol="<<stronger.reservoir.amount_mol
			<<" second_uptake_mol="<<second.reservoir.amount_mol<<'\n';
	}
	if(rank==0)std::cout<<"native wall reservoir exchange: PASS ranks="<<ranks
		<<" uptake_mol="<<uptake.reservoir.amount_mol
		<<" release_mol="<<release.reservoir.amount_mol
		<<" combined_residual_mol_s="<<uptake.combined_balance_defect_mol_s<<'\n';
	PetscFinalize();
}

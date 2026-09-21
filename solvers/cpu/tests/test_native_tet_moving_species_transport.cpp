#include "NativeTetMovingSpeciesTransport.hpp"
#include "NativeTetSpeciesVisualization.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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

std::vector<Vector> FluidVelocity(const iga::NativeTetMesh& mesh,
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

void Close(double actual,double expected,double tolerance=1e-11)
{
	assert(std::isfinite(actual));
	assert(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)));
}

}

int main()
{
	const auto previous=TwoCellMesh();
	const double dt=0.1;
	std::vector<Vector> zero(previous.points.size(),Vector{{0,0,0}});
	const std::vector<double> uniform(previous.points.size(),2.0);
	const auto stationary=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt);
	for(const double value:stationary.concentration_mol_m3)Close(value,2.);
	Close(stationary.current_inventory_mol,stationary.previous_inventory_mol);
	Close(stationary.balance_defect_mol_s,0.);
	const auto localized=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,
		std::vector<double>(previous.points.size(),0.),{},0.,0.,dt,
		false,0.,{},{},{1.,0.});
	Close(localized.source_mol_s,1./6.);
	Close(localized.current_inventory_mol,dt/6.);
	Close(localized.balance_defect_mol_s,0.);
	bool invalid_cell_source=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
			FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,
			false,0.,{},{},{1.});}
	catch(const std::invalid_argument&){invalid_cell_source=true;}
	assert(invalid_cell_source);

	auto translated=previous;
	for(auto& point:translated.points)point[0]+=0.1;
	const std::vector<Vector> translation(previous.points.size(),Vector{{1,0,0}});
	const std::vector<double> varying{1.,2.,3.,4.,5.};
	const auto rigid=iga::SolveNativeTetMovingSpeciesDenseStep(previous,translated,
		FluidVelocity(translated,translation),translation,varying,{},0.,0.,dt);
	for(std::size_t node=0;node<varying.size();++node)
		Close(rigid.concentration_mol_m3[node],varying[node]);
	Close(rigid.balance_defect_mol_s,0.);
	Close(rigid.outward_advective_flux_mol_s,0.);

	auto expanded=previous;
	for(auto& point:expanded.points)point[2]*=1.2;
	std::vector<Vector> expansion(previous.points.size());
	for(std::size_t node=0;node<expansion.size();++node)
		expansion[node][2]=(expanded.points[node][2]-previous.points[node][2])/dt;
	const auto moving=iga::SolveNativeTetMovingSpeciesDenseStep(previous,expanded,
		FluidVelocity(expanded,expansion),expansion,uniform,{},0.,0.,dt);
	for(const double value:moving.concentration_mol_m3)Close(value,2./1.2);
	Close(moving.current_inventory_mol,moving.previous_inventory_mol);
	Close(moving.outward_advective_flux_mol_s,0.);
	Close(moving.balance_defect_mol_s,0.);
	for(int rank=0;rank<2;++rank){
		const auto piece=iga::BuildNativeTetSpeciesVtkPartition(previous,expanded,
			moving.concentration_mol_m3,rank,2);
		assert(piece.cell_ids.size()==1&&piece.grid.types==std::vector<unsigned int>{10});
		assert(piece.point_arrays.size()==3);
		assert(piece.point_arrays[0].name=="concentration_mol_m3");
		assert(piece.point_arrays[1].name=="reference_position_m");
		assert(piece.point_arrays[2].name=="displacement_m");
		for(std::size_t local=0;local<piece.point_ids.size();++local){
			const auto global=static_cast<std::size_t>(piece.point_ids[local]);
			Close(piece.point_arrays[0].values[local],2./1.2);
			for(int axis=0;axis<3;++axis){
				Close(piece.grid.points[3*local+axis],expanded.points[global][axis]);
				Close(piece.point_arrays[1].values[3*local+axis],previous.points[global][axis]);
				Close(piece.point_arrays[2].values[3*local+axis],
					expanded.points[global][axis]-previous.points[global][axis]);
			}
		}
	}

	const std::vector<Vector> flow(previous.points.size(),Vector{{1,0,0}});
	const auto through=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,flow),zero,uniform,{{1,2.}},0.,0.,dt);
	for(const double value:through.concentration_mol_m3)Close(value,2.);
	Close(through.outward_advective_flux_mol_s,0.);
	Close(through.balance_defect_mol_s,0.);
	auto labelled=previous;
	for(auto& face:labelled.boundary_triangles){
		const bool inlet=std::all_of(face.nodes.begin(),face.nodes.end(),
			[&](std::uint32_t node){return labelled.points[node][0]==0.;});
		face.boundary_label=inlet?2:3;
	}
	const auto ports=iga::SolveNativeTetMovingSpeciesDenseStep(labelled,labelled,
		FluidVelocity(labelled,flow),zero,uniform,{{2,2.}},0.,0.,dt);
	Close(ports.outward_advective_flux_by_label_mol_s.at(2),-2.);
	Close(ports.outward_advective_flux_by_label_mol_s.at(3),2.);
	Close(ports.outward_relative_flow_by_label_m3_s.at(2),-1.);
	Close(ports.outward_relative_flow_by_label_m3_s.at(3),1.);
	Close(ports.outward_advective_flux_by_label_mol_s.at(2)
		+ports.outward_advective_flux_by_label_mol_s.at(3),
		ports.outward_advective_flux_mol_s);
	auto moved_labelled=labelled;
	for(auto& point:moved_labelled.points)point[0]+=0.1;
	const std::vector<Vector> moving_flow(labelled.points.size(),Vector{{2,0,0}});
	const std::vector<Vector> moving_grid(labelled.points.size(),Vector{{1,0,0}});
	const auto moving_ports=iga::SolveNativeTetMovingSpeciesDenseStep(
		labelled,moved_labelled,FluidVelocity(moved_labelled,moving_flow),
		moving_grid,uniform,{{2,2.}},0.,0.,dt);
	for(const auto& item:ports.outward_advective_flux_by_label_mol_s)
		Close(moving_ports.outward_advective_flux_by_label_mol_s.at(item.first),
			item.second);
	for(const auto& item:ports.outward_relative_flow_by_label_m3_s)
		Close(moving_ports.outward_relative_flow_by_label_m3_s.at(item.first),
			item.second);
	const auto topology=iga::BuildNativeTaylorHoodTopology(previous);
	std::vector<double> native_flow_state(3*(previous.points.size()
		+topology.edges.size())+previous.points.size(),0.);
	for(std::size_t node=0;node<previous.points.size()+topology.edges.size();++node)
		native_flow_state[3*node]=1.;
	const auto adapted=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		native_flow_state,zero,uniform,{{1,2.}},0.,0.,dt);
	for(std::size_t node=0;node<uniform.size();++node)
		Close(adapted.concentration_mol_m3[node],through.concentration_mol_m3[node]);
	Close(adapted.balance_defect_mol_s,through.balance_defect_mol_s);
	const std::vector<double> empty(previous.points.size(),0.);
	const auto high_peclet_galerkin=iga::SolveNativeTetMovingSpeciesDenseStep(
		previous,previous,FluidVelocity(previous,flow),zero,empty,{{1,1.}},
		0.,0.,0.01);
	assert(*std::min_element(high_peclet_galerkin.concentration_mol_m3.begin(),
		high_peclet_galerkin.concentration_mol_m3.end())<-0.1);
	const auto high_peclet_monotone=iga::SolveNativeTetMovingSpeciesDenseStep(
		previous,previous,FluidVelocity(previous,flow),zero,empty,{{1,1.}},
		0.,0.,0.01,true);
	assert(*std::min_element(high_peclet_monotone.concentration_mol_m3.begin(),
		high_peclet_monotone.concentration_mol_m3.end())>=-1e-12);
	Close(high_peclet_monotone.balance_defect_mol_s,0.);
	auto high_peclet_shifted=previous;
	for(auto& point:high_peclet_shifted.points)point[0]+=0.01;
	const std::vector<Vector> high_peclet_grid(previous.points.size(),Vector{{1,0,0}});
	const std::vector<Vector> high_peclet_flow(previous.points.size(),Vector{{2,0,0}});
	const auto high_peclet_moving=iga::SolveNativeTetMovingSpeciesDenseStep(
		previous,high_peclet_shifted,
		FluidVelocity(high_peclet_shifted,high_peclet_flow),high_peclet_grid,
		empty,{{1,1.}},0.,0.,0.01,true);
	for(std::size_t node=0;node<empty.size();++node)
		Close(high_peclet_moving.concentration_mol_m3[node],
			high_peclet_monotone.concentration_mol_m3[node]);
	Close(high_peclet_moving.balance_defect_mol_s,0.);
	const auto galerkin_matrix=iga::AssembleNativeTetMovingSpeciesStep(
		previous,previous,FluidVelocity(previous,flow),zero,empty,{{1,1.}},
		0.,0.,0.01);
	const auto monotone_matrix=iga::AssembleNativeTetMovingSpeciesStep(
		previous,previous,FluidVelocity(previous,flow),zero,empty,{{1,1.}},
		0.,0.,0.01,0,1,true);
	for(std::size_t column=0;column<previous.points.size();++column){
		double galerkin_sum=0.,monotone_sum=0.;
		for(std::size_t row=0;row<previous.points.size();++row){
			const auto key=std::make_pair(static_cast<std::uint32_t>(row),
				static_cast<std::uint32_t>(column));
			const auto galerkin=galerkin_matrix.matrix.find(key);
			const auto stabilized=monotone_matrix.matrix.find(key);
			if(galerkin!=galerkin_matrix.matrix.end())galerkin_sum+=galerkin->second;
			if(stabilized!=monotone_matrix.matrix.end()){
				monotone_sum+=stabilized->second;
				if(row!=column)assert(stabilized->second<=1e-12);
			}
		}
		Close(monotone_sum,galerkin_sum);
	}
	const auto incoming=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,flow),zero,empty,{{1,2.}},0.,0.,dt);
	assert(incoming.current_inventory_mol>0.);
	assert(incoming.outward_advective_flux_mol_s<0.);
	Close(incoming.balance_defect_mol_s,0.);
	const std::vector<Vector> doubled_flow(previous.points.size(),Vector{{2,0,0}});
	const auto moved_incoming=iga::SolveNativeTetMovingSpeciesDenseStep(previous,translated,
		FluidVelocity(translated,doubled_flow),translation,empty,{{1,2.}},0.,0.,dt);
	for(std::size_t node=0;node<empty.size();++node)
		Close(moved_incoming.concentration_mol_m3[node],
			incoming.concentration_mol_m3[node]);
	Close(moved_incoming.outward_advective_flux_mol_s,
		incoming.outward_advective_flux_mol_s);
	Close(moved_incoming.balance_defect_mol_s,0.);

	const auto diffused=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,varying,{},0.1,0.,dt);
	Close(diffused.current_inventory_mol,diffused.previous_inventory_mol);
	Close(diffused.balance_defect_mol_s,0.);
	assert(*std::max_element(diffused.concentration_mol_m3.begin(),
		diffused.concentration_mol_m3.end())<5.);
	assert(*std::min_element(diffused.concentration_mol_m3.begin(),
		diffused.concentration_mol_m3.end())>1.);
	const auto sourced=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,3.,dt);
	for(const double value:sourced.concentration_mol_m3)Close(value,2.3);
	Close(sourced.source_mol_s,1.);
	Close(sourced.balance_defect_mol_s,0.);
	const auto decayed=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.5,dt,false,2.);
	for(const double value:decayed.concentration_mol_m3)
		Close(value,(2.+dt*0.5)/(1.+dt*2.));
	Close(decayed.reaction_sink_mol_s,2.*decayed.current_inventory_mol);
	Close(decayed.balance_defect_mol_s,0.);
	const auto monotone_decay=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.5,dt,true,2.);
	for(const double value:monotone_decay.concentration_mol_m3)
		Close(value,(2.+dt*0.5)/(1.+dt*2.));
	Close(monotone_decay.reaction_sink_mol_s,2.*monotone_decay.current_inventory_mol);
	Close(monotone_decay.balance_defect_mol_s,0.);
	const auto no_decay_matrix=iga::AssembleNativeTetMovingSpeciesStep(previous,
		previous,FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt);
	const auto decay_matrix=iga::AssembleNativeTetMovingSpeciesStep(previous,
		previous,FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,
		0,1,false,2.);
	std::map<std::pair<std::uint32_t,std::uint32_t>,double> reaction_mass;
	for(const auto& cell:previous.cells){
		const double volume=iga::EvaluateNativeTetGeometry(previous,cell).determinant/6.;
		for(std::size_t row=0;row<4;++row)
			for(std::size_t column=0;column<4;++column)
				reaction_mass[{cell.nodes[row],cell.nodes[column]}]
					+=2.*volume*(row==column?2.:1.)/20.;
	}
	for(const auto& entry:reaction_mass)
		Close(decay_matrix.matrix.at(entry.first)
			-no_decay_matrix.matrix.at(entry.first),entry.second);
	const std::map<int,iga::NativeTetWallExchange> wall_sink{{1,{0.5,0.}}};
	const auto exchanged=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,false,0.,wall_sink);
	assert(exchanged.current_inventory_mol<exchanged.previous_inventory_mol);
	assert(exchanged.outward_wall_exchange_mol_s>0.);
	Close(exchanged.outward_wall_exchange_by_label_mol_s.at(1),
		exchanged.outward_wall_exchange_mol_s);
	Close(exchanged.balance_defect_mol_s,0.);
	const std::map<int,iga::NativeTetWallExchange> wall_donor{{1,{0.5,3.}}};
	const auto supplied=iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,false,0.,wall_donor);
	assert(supplied.current_inventory_mol>supplied.previous_inventory_mol);
	assert(supplied.outward_wall_exchange_mol_s<0.);
	Close(supplied.balance_defect_mol_s,0.);
	const auto wall_matrix=iga::AssembleNativeTetMovingSpeciesStep(previous,
		previous,FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,
		0,1,false,0.,wall_sink);
	const auto donor_matrix=iga::AssembleNativeTetMovingSpeciesStep(previous,
		previous,FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,
		0,1,false,0.,wall_donor);
	std::map<std::pair<std::uint32_t,std::uint32_t>,double> expected_wall_matrix;
	std::vector<double> expected_wall_rhs(previous.points.size(),0.);
	for(const auto& face:previous.boundary_triangles){
		const auto& x=previous.points[face.nodes[0]];
		const auto& y=previous.points[face.nodes[1]];
		const auto& z=previous.points[face.nodes[2]];
		const Vector a{{y[0]-x[0],y[1]-x[1],y[2]-x[2]}};
		const Vector b{{z[0]-x[0],z[1]-x[1],z[2]-x[2]}};
		const Vector cross{{a[1]*b[2]-a[2]*b[1],
			a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
		const double area=0.5*std::hypot(cross[0],cross[1],cross[2]);
		for(std::size_t row=0;row<3;++row){
			expected_wall_rhs[face.nodes[row]]+=0.5*3.*area/3.;
			for(std::size_t column=0;column<3;++column)
				expected_wall_matrix[{face.nodes[row],face.nodes[column]}]
					+=0.5*area*(row==column?2.:1.)/12.;
		}
	}
	for(const auto& entry:expected_wall_matrix)
		Close(wall_matrix.matrix.at(entry.first)
			-no_decay_matrix.matrix.at(entry.first),entry.second);
	for(std::size_t node=0;node<expected_wall_rhs.size();++node)
		Close(donor_matrix.rhs[node]-no_decay_matrix.rhs[node],expected_wall_rhs[node]);
	const auto moved_exchange=iga::SolveNativeTetMovingSpeciesDenseStep(previous,
		translated,FluidVelocity(translated,translation),translation,uniform,{},
		0.,0.,dt,false,0.,wall_sink);
	for(std::size_t node=0;node<uniform.size();++node)
		Close(moved_exchange.concentration_mol_m3[node],
			exchanged.concentration_mol_m3[node]);
	Close(moved_exchange.outward_wall_exchange_mol_s,
		exchanged.outward_wall_exchange_mol_s);
	Close(moved_exchange.balance_defect_mol_s,0.);
	bool invalid_wall=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
			FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,false,0.,
			{{99,{0.5,0.}}});}
	catch(const std::invalid_argument&){invalid_wall=true;}
	assert(invalid_wall);
	invalid_wall=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
			FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,false,0.,
			{{1,{-0.5,0.}}});}
	catch(const std::invalid_argument&){invalid_wall=true;}
	assert(invalid_wall);
	invalid_wall=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
			FluidVelocity(previous,translation),zero,uniform,{},0.,0.,dt,false,0.,
			wall_sink);}
	catch(const std::invalid_argument&){invalid_wall=true;}
	assert(invalid_wall);
	bool invalid_decay=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
			FluidVelocity(previous,zero),zero,uniform,{},0.,0.,dt,false,-1.);}
	catch(const std::invalid_argument&){invalid_decay=true;}
	assert(invalid_decay);

	bool rejected=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,flow),zero,uniform,{},0.,0.,dt);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,translated,
		FluidVelocity(translated,translation),zero,uniform,{},0.,0.,dt);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		std::vector<double>{1.,2.},zero,uniform,{{1,2.}},0.,0.,dt);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,previous,
		FluidVelocity(previous,zero),zero,uniform,{},0.,0.,0.);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	auto changed=previous;std::swap(changed.cells[0].nodes[1],changed.cells[0].nodes[2]);
	rejected=false;
	try{iga::SolveNativeTetMovingSpeciesDenseStep(previous,changed,
		FluidVelocity(changed,zero),zero,uniform,{},0.,0.,dt);}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	std::cout<<"native moving-tetra species transport passed volume_ratio=1.2"
		<<" inventory_mol="<<moving.current_inventory_mol
		<<" open_inflow_inventory_mol="<<incoming.current_inventory_mol<<'\n';
}

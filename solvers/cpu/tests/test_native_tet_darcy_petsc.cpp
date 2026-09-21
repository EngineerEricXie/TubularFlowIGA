#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetVesselTissueSourceMap.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

iga::NativeTetMesh Cube(unsigned subdivisions)
{
	iga::NativeTetMesh mesh;
	const auto node=[&](unsigned x,unsigned y,unsigned z){
		return static_cast<std::uint32_t>((z*(subdivisions+1)+y)
			*(subdivisions+1)+x);
	};
	for(unsigned z=0;z<=subdivisions;++z)
		for(unsigned y=0;y<=subdivisions;++y)
			for(unsigned x=0;x<=subdivisions;++x)
				mesh.points.push_back({{double(x)/subdivisions,
					double(y)/subdivisions,double(z)/subdivisions}});
	const std::array<std::array<unsigned,4>,6> pattern{{
		{{0,1,3,7}},{{0,3,2,7}},{{0,2,6,7}},
		{{0,6,4,7}},{{0,4,5,7}},{{0,5,1,7}}}};
	std::map<std::array<std::uint32_t,3>,unsigned> uses;
	std::uint64_t id=1;
	for(unsigned z=0;z<subdivisions;++z)
		for(unsigned y=0;y<subdivisions;++y)
			for(unsigned x=0;x<subdivisions;++x){
				const std::array<std::uint32_t,8> corners{{
					node(x,y,z),node(x+1,y,z),node(x,y+1,z),node(x+1,y+1,z),
					node(x,y,z+1),node(x+1,y,z+1),node(x,y+1,z+1),
					node(x+1,y+1,z+1)}};
				for(const auto& local:pattern){
					iga::NativeTetCell cell;
					cell.id=id++;
					for(std::size_t i=0;i<4;++i)cell.nodes[i]=corners[local[i]];
					mesh.cells.push_back(cell);
					for(std::size_t opposite=0;opposite<4;++opposite){
						std::array<std::uint32_t,3> face{};
						std::size_t entry=0;
						for(std::size_t local_node=0;local_node<4;++local_node)
							if(local_node!=opposite)face[entry++]=cell.nodes[local_node];
						std::sort(face.begin(),face.end());++uses[face];
					}
				}
			}
	for(const auto& item:uses){
		if(item.second!=1)continue;
		const auto& face=item.first;
		const bool plane0=std::all_of(face.begin(),face.end(),[&](auto vertex){
			return mesh.points[vertex][0]==0.;});
		const bool plane1=std::all_of(face.begin(),face.end(),[&](auto vertex){
			return mesh.points[vertex][0]==1.;});
		mesh.boundary_triangles.push_back({id++,face,plane0?1:(plane1?2:3)});
	}
	return mesh;
}

void Close(double actual,double expected,double tolerance=1e-10)
{
	assert(std::isfinite(actual));
	if(std::abs(actual-expected)>tolerance*std::max(1.,std::abs(expected)))
		std::cerr<<"Darcy mismatch actual="<<actual<<" expected="<<expected<<'\n';
	assert(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)));
}

double QuadraticError(const iga::NativeTetMesh& mesh,
	const iga::NativeTetDarcyResult& solution)
{
	double squared=0.,total_volume=0.;
	for(const auto& cell:mesh.cells){
		const double volume=iga::EvaluateNativeTetGeometry(mesh,cell).determinant/6.;
		double x=0.,computed=0.;
		for(const auto node:cell.nodes){
			x+=mesh.points[node][0]/4.;
			computed+=solution.pressure_pa[node]/4.;
		}
		const double exact=0.5*x*(1.-x);
		squared+=volume*std::pow(computed-exact,2);
		total_volume+=volume;
	}
	return std::sqrt(squared/total_volume);
}

void CheckConservativeFlux(const iga::NativeTetMesh& mesh,
	const iga::NativeTetDarcyResult& solution,const std::vector<double>& source)
{
	assert(solution.conservative_face_flow_m3_s.size()==mesh.cells.size());
	assert(solution.maximum_cell_balance_defect_m3_s<1e-10);
	std::map<std::array<std::uint32_t,3>,double> first_flux;
	std::map<std::array<std::uint32_t,3>,unsigned> uses;
	for(std::size_t cell_index=0;cell_index<mesh.cells.size();++cell_index){
		const auto& cell=mesh.cells[cell_index];
		double balance=0.;
		for(std::size_t opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> face{};
			std::size_t next=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[next++]=cell.nodes[local];
			std::sort(face.begin(),face.end());
			const double flow=solution.conservative_face_flow_m3_s[cell_index][opposite];
			std::array<double,3> face_centroid{},area_normal{};
			for(const auto vertex:face)
				for(int axis=0;axis<3;++axis)
					face_centroid[axis]+=mesh.points[vertex][axis]/3.;
			const auto& a=mesh.points[face[0]];
			const auto& b=mesh.points[face[1]];
			const auto& c=mesh.points[face[2]];
			const auto& interior=mesh.points[cell.nodes[opposite]];
			area_normal={{
				0.5*((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])),
				0.5*((b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2])),
				0.5*((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))}};
		if(area_normal[0]*(interior[0]-a[0])
			+area_normal[1]*(interior[1]-a[1])
			+area_normal[2]*(interior[2]-a[2])>0.)
			for(auto& component:area_normal)component=-component;
		const auto rt0=iga::EvaluateNativeTetDarcyRt0Flux(mesh,solution,
			cell_index,face_centroid);
		const double integrated=rt0[0]*area_normal[0]
			+rt0[1]*area_normal[1]+rt0[2]*area_normal[2];
		Close(integrated,flow,1e-12);
			balance+=flow;
			if(uses[face]++)Close(flow,-first_flux.at(face));
			else first_flux[face]=flow;
		}
		const double volume=iga::EvaluateNativeTetGeometry(mesh,cell).determinant/6.;
		Close(balance,source[cell_index]*volume,1e-10);
		std::array<double,3> centroid{};
		for(const auto node:cell.nodes)
			for(int axis=0;axis<3;++axis)
				centroid[axis]+=mesh.points[node][axis]/4.;
		const auto central=iga::EvaluateNativeTetDarcyRt0Flux(mesh,solution,
			cell_index,centroid);
		const double step=0.125*std::cbrt(volume);
		double divergence=0.;
		for(int axis=0;axis<3;++axis){
			auto shifted=centroid;
			shifted[axis]+=step;
			const auto next_flux=iga::EvaluateNativeTetDarcyRt0Flux(mesh,
				solution,cell_index,shifted);
			divergence+=(next_flux[axis]-central[axis])/step;
		}
		Close(divergence,source[cell_index],1e-9);
	}
	for(const auto& item:uses)assert(item.second<=2);
	double boundary_total=0.;
	for(const auto& item:solution.conservative_outward_boundary_flow_m3_s)
		boundary_total+=item.second;
	Close(boundary_total,solution.volume_source_m3_s,1e-10);
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const auto mesh=Cube(2);
	const std::vector<double> mobility(mesh.cells.size(),2.);
	const std::vector<double> zero(mesh.cells.size(),0.);
	const auto linear=iga::SolveNativeTetDarcyPetsc(mesh,mobility,zero,
		{{1,1.},{2,0.}});
	CheckConservativeFlux(mesh,linear,zero);
	for(std::size_t node=0;node<mesh.points.size();++node)
		Close(linear.pressure_pa[node],1.-mesh.points[node][0]);
	for(const auto& flux:linear.cell_flux_m_s){
		Close(flux[0],2.);
		Close(flux[1],0.);
		Close(flux[2],0.);
	}
	for(std::size_t index=0;index<mesh.cells.size();++index){
		const auto& cell=mesh.cells[index];
		std::array<double,3> centroid{};
		for(const auto node:cell.nodes)
			for(int axis=0;axis<3;++axis)
				centroid[axis]+=mesh.points[node][axis]/4.;
		const auto recovered=iga::EvaluateNativeTetDarcyRt0Flux(mesh,linear,
			index,centroid);
		Close(recovered[0],2.);
		Close(recovered[1],0.);
		Close(recovered[2],0.);
	}
	Close(linear.outward_boundary_flow_m3_s.at(1),-2.);
	Close(linear.outward_boundary_flow_m3_s.at(2),2.);
	Close(linear.outward_boundary_flow_m3_s.at(3),0.);
	Close(linear.volume_source_m3_s,0.);
	Close(linear.conservative_outward_boundary_flow_m3_s.at(1),-2.);
	Close(linear.conservative_outward_boundary_flow_m3_s.at(2),2.);
	Close(linear.conservative_outward_boundary_flow_m3_s.at(3),0.);
	// The same native FEM transport assembly consumes RT0 directly at cell
	// and boundary quadrature points, without converting it to shared P2 nodes.
	const std::vector<std::array<double,3>> stationary_mesh_velocity(
		mesh.points.size());
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::vector<std::array<double,3>> p2_velocity(
		mesh.points.size()+topology.edges.size(),{{2.,0.,0.}});
	const std::vector<double> previous_tracer(mesh.points.size(),1.);
	const auto rt0_tracer=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,mesh,
		{},stationary_mesh_velocity,previous_tracer,{{1,2.}},0.1,0.,0.01,
		false,0.,{},linear.conservative_face_flow_m3_s);
	const auto p2_tracer=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,mesh,
		p2_velocity,stationary_mesh_velocity,previous_tracer,{{1,2.}},
		0.1,0.,0.01);
	for(std::size_t node=0;node<mesh.points.size();++node)
		Close(rt0_tracer.step.concentration_mol_m3[node],
			p2_tracer.step.concentration_mol_m3[node],1e-10);
	for(const auto& item:linear.conservative_outward_boundary_flow_m3_s)
		Close(rt0_tracer.step.outward_relative_flow_by_label_m3_s.at(item.first),
			item.second,1e-10);
	Close(rt0_tracer.step.balance_defect_mol_s,0.,1e-10);
	const auto flux_boundary=iga::SolveNativeTetDarcyPetsc(mesh,mobility,zero,
		{{1,1.}},{{2,2.}});
	CheckConservativeFlux(mesh,flux_boundary,zero);
	Close(flux_boundary.conservative_outward_boundary_flow_m3_s.at(2),2.);
	for(std::size_t node=0;node<mesh.points.size();++node)
		Close(flux_boundary.pressure_pa[node],linear.pressure_pa[node]);
	Close(flux_boundary.outward_boundary_flow_m3_s.at(2),2.);
	const auto stronger=iga::SolveNativeTetDarcyPetsc(mesh,
		std::vector<double>(mesh.cells.size(),4.),zero,{{1,1.},{2,0.}});
	for(std::size_t node=0;node<mesh.points.size();++node)
		Close(stronger.pressure_pa[node],linear.pressure_pa[node]);
	Close(stronger.outward_boundary_flow_m3_s.at(2),4.);
	const auto coarse=iga::SolveNativeTetDarcyPetsc(mesh,
		std::vector<double>(mesh.cells.size(),1.),
		std::vector<double>(mesh.cells.size(),1.),{{1,0.},{2,0.}});
	CheckConservativeFlux(mesh,coarse,
		std::vector<double>(mesh.cells.size(),1.));
	Close(coarse.conservative_outward_boundary_flow_m3_s.at(3),0.);
	const auto refined_mesh=Cube(4);
	const auto fine=iga::SolveNativeTetDarcyPetsc(refined_mesh,
		std::vector<double>(refined_mesh.cells.size(),1.),
		std::vector<double>(refined_mesh.cells.size(),1.),{{1,0.},{2,0.}});
	CheckConservativeFlux(refined_mesh,fine,
		std::vector<double>(refined_mesh.cells.size(),1.));
	const double coarse_error=QuadraticError(mesh,coarse);
	const double fine_error=QuadraticError(refined_mesh,fine);
	assert(coarse_error<0.04&&fine_error<0.01&&fine_error<coarse_error);
	Close(coarse.volume_source_m3_s,1.);
	Close(fine.volume_source_m3_s,1.);
	const auto mapped=iga::MapNativeTetVesselFlowToTissueSource(mesh,{
		{"vessel_a",0.6,{{1,0.25},{2,0.75}}},
		{"vessel_b",-0.2,{{3,1.}}}});
	Close(mapped.tissue_source_s_inv[0],7.2);
	Close(mapped.tissue_source_s_inv[1],21.6);
	Close(mapped.tissue_source_s_inv[2],-9.6);
	Close(mapped.total_vessel_outward_flow_m3_s,0.4);
	Close(mapped.total_tissue_inward_flow_m3_s,0.4);
	Close(mapped.maximum_port_balance_defect_m3_s,0.);
	const auto mapped_darcy=iga::SolveNativeTetDarcyPetsc(mesh,
		std::vector<double>(mesh.cells.size(),1.),mapped.tissue_source_s_inv,
		{{1,0.},{2,0.}});
	CheckConservativeFlux(mesh,mapped_darcy,mapped.tissue_source_s_inv);
	Close(mapped_darcy.volume_source_m3_s,0.4);
	// A manufactured constant project-owned P2/P1 vessel field supplies the
	// measured port Q directly; no externally entered flow number is used.
	const std::size_t vessel_velocity_nodes=mesh.points.size()+topology.edges.size();
	std::vector<double> vessel_state(3*vessel_velocity_nodes+mesh.points.size(),0.);
	for(std::size_t node=0;node<vessel_velocity_nodes;++node)
		vessel_state[3*node]=0.6;
	const auto measured_ports=iga::MapNativeTetVesselStateToTissueSource(
		mesh,topology,vessel_state,mesh,
		{{"outward_x1",2,{{1,0.25},{2,0.75}}},
			{"inward_x0",1,{{3,1.}}}});
	Close(measured_ports.vessel_outward_flow_m3_s.at("outward_x1"),0.6);
	Close(measured_ports.vessel_outward_flow_m3_s.at("inward_x0"),-0.6);
	Close(measured_ports.total_tissue_inward_flow_m3_s,0.);
	Close(measured_ports.maximum_port_balance_defect_m3_s,0.);
	const auto measured_darcy=iga::SolveNativeTetDarcyPetsc(mesh,
		std::vector<double>(mesh.cells.size(),1.),measured_ports.tissue_source_s_inv,
		{{1,0.},{2,0.}});
	CheckConservativeFlux(mesh,measured_darcy,measured_ports.tissue_source_s_inv);
	Close(measured_darcy.volume_source_m3_s,0.);
	// Two independently numbered tetrahedra share one exact triangular face.
	// The vessel lies at x<0, tissue at x>0, and their outward normals oppose.
	iga::NativeTetMesh adjacent_vessel,adjacent_tissue;
	adjacent_vessel.points={{{0.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}},{{-1.,0.,0.}}};
	adjacent_vessel.cells.push_back({10,{{0,2,1,3}}});
	adjacent_vessel.boundary_triangles.push_back({1,{{0,1,2}},10});
	adjacent_tissue.points={{{0.,0.,0.}},{{0.,1.,0.}},{{0.,0.,1.}},{{1.,0.,0.}}};
	adjacent_tissue.cells.push_back({11,{{0,1,2,3}}});
	adjacent_tissue.boundary_triangles={
		{1,{{0,1,2}},20},{2,{{1,2,3}},2},
		{3,{{0,2,3}},3},{4,{{0,1,3}},3}};
	const auto adjacent_topology=iga::BuildNativeTaylorHoodTopology(adjacent_vessel);
	const std::size_t adjacent_velocity_nodes=adjacent_vessel.points.size()
		+adjacent_topology.edges.size();
	std::vector<double> adjacent_state(3*adjacent_velocity_nodes
		+adjacent_vessel.points.size(),0.);
	for(std::size_t node=0;node<adjacent_velocity_nodes;++node)
		adjacent_state[3*node]=0.6;
	const auto adjacent_map=iga::MapNativeTetMatchingInterfaceToTissueSource(
		adjacent_vessel,adjacent_topology,adjacent_state,adjacent_tissue,
		{{"matched_face",10,20}});
	Close(iga::EvaluateNativeTetBoundaryFlows(adjacent_vessel,adjacent_topology,
		adjacent_state).at(10).outward_flow_m3_s,0.3);
	Close(adjacent_map.source.vessel_outward_flow_m3_s.at("matched_face"),0.3);
	Close(adjacent_map.source.tissue_inward_flow_m3_s.at("matched_face"),0.3);
	Close(adjacent_map.source.maximum_port_balance_defect_m3_s,0.);
	Close(adjacent_map.source.tissue_source_s_inv.at(0),1.8);
	assert(adjacent_map.matched_facets.at("matched_face")==1);
	assert(adjacent_map.facet_transfers.size()==1);
	assert(adjacent_map.facet_transfers[0].port_name=="matched_face");
	assert(adjacent_map.facet_transfers[0].vessel_cell_id==adjacent_vessel.cells[0].id);
	assert(adjacent_map.facet_transfers[0].tissue_cell_id==adjacent_tissue.cells[0].id);
	Close(adjacent_map.facet_transfers[0].vessel_outward_flow_m3_s,0.3);
	const auto adjacent_darcy=iga::SolveNativeTetDarcyPetsc(adjacent_tissue,
		{1.},adjacent_map.source.tissue_source_s_inv,{{2,0.}});
	CheckConservativeFlux(adjacent_tissue,adjacent_darcy,
		adjacent_map.source.tissue_source_s_inv);
	Close(adjacent_darcy.volume_source_m3_s,0.3);
	Close(adjacent_darcy.conservative_outward_boundary_flow_m3_s.at(2),0.3);
	for(std::size_t node=0;node<adjacent_velocity_nodes;++node) {
		std::array<double,3> point{};
		if(node<adjacent_vessel.points.size()) point=adjacent_vessel.points[node];
		else {
			const auto& edge=adjacent_topology.edges[node-adjacent_vessel.points.size()];
			for(int axis=0;axis<3;++axis)
				point[axis]=0.5*(adjacent_vessel.points[edge[0]][axis]
					+adjacent_vessel.points[edge[1]][axis]);
		}
		adjacent_state[3*node]=0.2+0.4*point[1]*point[1]
			+0.3*point[1]*point[2];
	}
	const auto quadratic_map=iga::MapNativeTetMatchingInterfaceToTissueSource(
		adjacent_vessel,adjacent_topology,adjacent_state,adjacent_tissue,
		{{"matched_face",10,20}});
	Close(quadratic_map.source.vessel_outward_flow_m3_s.at("matched_face"),
		iga::EvaluateNativeTetBoundaryFlows(adjacent_vessel,adjacent_topology,
			adjacent_state).at(10).outward_flow_m3_s,1e-12);
	// Two nonuniform interface facets must feed their own adjacent tissue cells.
	auto two_face_vessel=adjacent_vessel;
	auto two_face_tissue=adjacent_tissue;
	const auto append_shifted_pair=[](iga::NativeTetMesh& mesh){
		const auto original=mesh;
		const auto offset=static_cast<std::uint32_t>(original.points.size());
		for(auto point:original.points){point[1]+=2.;mesh.points.push_back(point);}
		for(auto cell:original.cells){
			cell.id+=100;
			for(auto& node:cell.nodes)node+=offset;
			mesh.cells.push_back(cell);
		}
		for(auto face:original.boundary_triangles){
			face.id+=100;
			for(auto& node:face.nodes)node+=offset;
			mesh.boundary_triangles.push_back(face);
		}
	};
	append_shifted_pair(two_face_vessel);
	append_shifted_pair(two_face_tissue);
	const auto two_face_topology=iga::BuildNativeTaylorHoodTopology(two_face_vessel);
	const std::size_t two_face_velocity_nodes=two_face_vessel.points.size()
		+two_face_topology.edges.size();
	std::vector<double> two_face_state(3*two_face_velocity_nodes
		+two_face_vessel.points.size(),0.);
	for(std::size_t node=0;node<two_face_velocity_nodes;++node){
		double y=0.;
		if(node<two_face_vessel.points.size())y=two_face_vessel.points[node][1];
		else {
			const auto& edge=two_face_topology.edges[node-two_face_vessel.points.size()];
			y=0.5*(two_face_vessel.points[edge[0]][1]
				+two_face_vessel.points[edge[1]][1]);
		}
		two_face_state[3*node]=0.6+0.2*y;
	}
	const auto two_face_map=iga::MapNativeTetMatchingInterfaceToTissueSource(
		two_face_vessel,two_face_topology,two_face_state,two_face_tissue,
		{{"two_faces",10,20}});
	assert(two_face_map.matched_facets.at("two_faces")==2);
	assert(two_face_map.facet_transfers.size()==2);
	Close(two_face_map.facet_transfers[0].vessel_outward_flow_m3_s
		+two_face_map.facet_transfers[1].vessel_outward_flow_m3_s,13./15.,1e-12);
	Close(two_face_map.source.tissue_source_s_inv[0]/6.,1./3.,1e-12);
	Close(two_face_map.source.tissue_source_s_inv[1]/6.,8./15.,1e-12);
	Close(two_face_map.source.vessel_outward_flow_m3_s.at("two_faces"),13./15.,1e-12);
	for(std::size_t node=0;node<adjacent_velocity_nodes;++node)
		adjacent_state[3*node]=-0.6;
	const auto reversed_map=iga::MapNativeTetMatchingInterfaceToTissueSource(
		adjacent_vessel,adjacent_topology,adjacent_state,adjacent_tissue,
		{{"matched_face",10,20}});
	Close(reversed_map.source.vessel_outward_flow_m3_s.at("matched_face"),-0.3);
	Close(reversed_map.source.tissue_source_s_inv.at(0),-1.8);
	const auto reversed_darcy=iga::SolveNativeTetDarcyPetsc(adjacent_tissue,
		{1.},reversed_map.source.tissue_source_s_inv,{{2,0.}});
	CheckConservativeFlux(adjacent_tissue,reversed_darcy,
		reversed_map.source.tissue_source_s_inv);
	Close(reversed_darcy.conservative_outward_boundary_flow_m3_s.at(2),-0.3);
	bool invalid_interface=false;
	auto same_side_tissue=adjacent_tissue;
	same_side_tissue.points[3]={{-1.,0.,0.}};
	same_side_tissue.cells[0].nodes={{0,2,1,3}};
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,adjacent_state,same_side_tissue,{{"same_side",10,20}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	auto displaced_tissue=adjacent_tissue;
	displaced_tissue.points[0][0]=0.1;
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,adjacent_state,displaced_tissue,{{"gap",10,20}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,adjacent_state,adjacent_tissue,
		{{"one",10,20},{"two",10,20}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	auto duplicated_port_vessel=adjacent_vessel;
	auto duplicated_port_tissue=adjacent_tissue;
	duplicated_port_vessel.boundary_triangles.push_back({5,{{0,1,2}},11});
	duplicated_port_tissue.boundary_triangles.push_back({5,{{0,1,2}},21});
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(duplicated_port_vessel,
		adjacent_topology,adjacent_state,duplicated_port_tissue,
		{{"first",10,20},{"second",11,21}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,adjacent_state,adjacent_tissue,{{"missing",10,99}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	auto duplicate_tissue=adjacent_tissue;
	duplicate_tissue.boundary_triangles.push_back({5,{{0,1,2}},20});
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,adjacent_state,duplicate_tissue,{{"duplicate",10,20}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	invalid_interface=false;
	auto invalid_state=adjacent_state;
	invalid_state[0]=std::numeric_limits<double>::quiet_NaN();
	try{iga::MapNativeTetMatchingInterfaceToTissueSource(adjacent_vessel,
		adjacent_topology,invalid_state,adjacent_tissue,{{"nan",10,20}});}
	catch(const std::invalid_argument&){invalid_interface=true;}
	assert(invalid_interface);
	bool invalid_map=false;
	try{iga::MapNativeTetVesselFlowToTissueSource(mesh,
		{{"vessel_a",1.,{{1,0.4},{2,0.4}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	invalid_map=false;
	try{iga::MapNativeTetVesselStateToTissueSource(mesh,topology,vessel_state,mesh,
		{{"missing",99,{{1,1.}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	invalid_map=false;
	try{iga::MapNativeTetVesselStateToTissueSource(mesh,topology,vessel_state,mesh,
		{{"a",2,{{1,1.}}},{"b",2,{{2,1.}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	invalid_map=false;
	auto nonfinite_state=vessel_state;nonfinite_state[0]=std::numeric_limits<double>::quiet_NaN();
	try{iga::MapNativeTetVesselStateToTissueSource(mesh,topology,nonfinite_state,mesh,
		{{"a",2,{{1,1.}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	invalid_map=false;
	try{iga::MapNativeTetVesselFlowToTissueSource(mesh,
		{{"vessel_a",1.,{{999,1.}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	invalid_map=false;
	try{iga::MapNativeTetVesselFlowToTissueSource(mesh,
		{{"vessel_a",1.,{{1,1.}}},{"vessel_a",1.,{{2,1.}}}});}
	catch(const std::invalid_argument&){invalid_map=true;}
	assert(invalid_map);
	bool rejected=false;
	try{iga::SolveNativeTetDarcyPetsc(mesh,mobility,zero,{{99,0.}});}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetDarcyPetsc(mesh,mobility,zero,{{1,0.},{2,0.}},
			{{2,1.}});}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetDarcyPetsc(mesh,mobility,zero,{{1,0.},{3,1.}});}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	rejected=false;
	try{iga::SolveNativeTetDarcyPetsc(mesh,
			std::vector<double>(mesh.cells.size(),0.),zero,{{1,1.},{2,0.}});}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	if(argc>1){
		std::ifstream input(argv[1]);
		if(!input)throw std::runtime_error("native Darcy external mesh is missing");
		const auto pipe=iga::ReadNativeTetMeshGmsh41(input);
		const std::vector<double> pipe_mobility(pipe.cells.size(),1e-4);
		const auto pipe_solution=iga::SolveNativeTetDarcyPetsc(pipe,pipe_mobility,
			std::vector<double>(pipe.cells.size(),0.),{{1,1.},{2,0.}});
		CheckConservativeFlux(pipe,pipe_solution,
			std::vector<double>(pipe.cells.size(),0.));
		Close(pipe_solution.conservative_outward_boundary_flow_m3_s.at(0),0.);
		const auto bounds=std::minmax_element(pipe.points.begin(),pipe.points.end(),
			[](const auto& left,const auto& right){return left[0]<right[0];});
		const double first=bounds.first->at(0),last=bounds.second->at(0);
		double max_pressure_error=0.,max_flux_error=0.;
		for(std::size_t node=0;node<pipe.points.size();++node)
			max_pressure_error=std::max(max_pressure_error,std::abs(
				pipe_solution.pressure_pa[node]
				-(last-pipe.points[node][0])/(last-first)));
		for(const auto& flux:pipe_solution.cell_flux_m_s){
			max_flux_error=std::max({max_flux_error,
				std::abs(flux[0]-1e-4/(last-first)),std::abs(flux[1]),
				std::abs(flux[2])});
		}
		if(rank==0)std::cerr<<"pipe diagnostic pressure_error="
			<<max_pressure_error<<" flux_error="<<max_flux_error
			<<" inlet="<<pipe_solution.outward_boundary_flow_m3_s.at(1)
			<<" outlet="<<pipe_solution.outward_boundary_flow_m3_s.at(2)
			<<" wall="<<pipe_solution.outward_boundary_flow_m3_s.at(0)<<'\n';
		if(max_pressure_error>5e-3||max_flux_error>2e-4)
			throw std::runtime_error("native Darcy generated-pipe field gate failed");
		assert(pipe_solution.outward_boundary_flow_m3_s.at(1)<0.);
		assert(pipe_solution.outward_boundary_flow_m3_s.at(2)>0.);
		assert(std::abs(pipe_solution.outward_boundary_flow_m3_s.at(1)
			+pipe_solution.outward_boundary_flow_m3_s.at(2))<1e-9);
		assert(std::abs(pipe_solution.outward_boundary_flow_m3_s.at(0))<1e-9);
	}
	if(rank==0)std::cout<<"native_tet_darcy: PASS ranks="<<ranks
		<<" coarse_error="<<coarse_error<<" fine_error="<<fine_error
		<<" outlet_m3_s="<<linear.outward_boundary_flow_m3_s.at(2)<<'\n';
	PetscFinalize();
}

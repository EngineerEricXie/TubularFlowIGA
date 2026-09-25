#include "NativeTetAleFsiRuntime.hpp"
#include "NativeTetFsiCheckpoint.hpp"
#include "NativeTetSolidFsiRuntime.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
constexpr double length=.02,height=.004,width=.004;
std::string Hash(const char* text){iga::Sha256 hash;hash.Append(text,std::char_traits<char>::length(text));return hash.Hex();}
std::uint32_t Node(int i,int j,int k){return static_cast<std::uint32_t>((i*2+j)*3+k);}
double Determinant(const iga::NativeTetMesh& mesh,const std::array<std::uint32_t,4>& n)
{return iga::NativeAleDeterminant(mesh.points[n[0]],mesh.points[n[1]],mesh.points[n[2]],mesh.points[n[3]]);}
iga::NativeTetMesh FluidMesh()
{
	iga::NativeTetMesh mesh;
	for(int i=0;i<3;++i)for(int j=0;j<2;++j)for(int k=0;k<3;++k)
		mesh.points.push_back({{length*i/2.,height*j,width*k/2.}});
	std::uint64_t id=1;
	for(int i=0;i<2;++i)for(int k=0;k<2;++k){
		const auto a=Node(i,0,k),b=Node(i+1,0,k),c=Node(i,1,k),d=Node(i+1,1,k);
		const auto e=Node(i,0,k+1),f=Node(i+1,0,k+1),g=Node(i,1,k+1),h=Node(i+1,1,k+1);
		const std::array<std::array<std::uint32_t,4>,6> cells{{{{a,b,d,h}},{{a,d,c,h}},{{a,c,g,h}},
			{{a,g,e,h}},{{a,e,f,h}},{{a,f,b,h}}}};
		for(auto nodes:cells){if(Determinant(mesh,nodes)<0.)std::swap(nodes[1],nodes[2]);
			mesh.cells.push_back({id++,nodes});}
	}
	using Face=std::array<std::uint32_t,3>;std::map<Face,std::pair<int,Face>> faces;
	for(const auto& cell:mesh.cells)for(std::size_t opposite=0;opposite<4;++opposite){
		Face oriented{};std::size_t p=0;for(std::size_t local=0;local<4;++local)if(local!=opposite)oriented[p++]=cell.nodes[local];
		auto key=oriented;std::sort(key.begin(),key.end());auto& record=faces[key];++record.first;record.second=oriented;
	}
	id=1;
	for(const auto& item:faces)if(item.second.first==1){const auto face=item.second.second;
		auto all=[&](int component,double value){for(auto node:face)if(std::abs(mesh.points[node][component]-value)>1e-14)return false;return true;};
		int label=0;if(all(0,0.))label=1;else if(all(0,length))label=2;else if(all(1,height))label=7;
		mesh.boundary_triangles.push_back({id++,face,label});
	}
	return mesh;
}
iga::DistributedSurfaceLayout Layout(const iga::NativeTetMesh& mesh)
{
	iga::DistributedSurfaceLayout value;value.reference_mesh_identity_sha256=Hash("native-compliant-channel");
	std::set<std::uint64_t> ids;
	for(const auto& triangle:mesh.boundary_triangles)if(triangle.boundary_label==7){
		value.reference_triangles.push_back({{triangle.nodes[0],triangle.nodes[1],triangle.nodes[2]}});
		for(auto node:triangle.nodes)ids.insert(node);
	}
	value.global_node_count=ids.size();value.partition_count=1;value.partition_rank=0;
	value.owned_global_node_ids.assign(ids.begin(),ids.end());
	for(auto node:ids)value.reference_positions.push_back({node,mesh.points[node]});
	value.owned_reference_lumped_areas_m2.assign(ids.size(),0.0);
	for(const auto& triangle:value.reference_triangles){
		const auto&a=mesh.points[triangle[0]],&b=mesh.points[triangle[1]],&c=mesh.points[triangle[2]];
		std::array<double,3> ab{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}},ac{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
		const double area=.5*iga::NativeTetInterfaceNorm(iga::NativeTetInterfaceCross(ab,ac));
		for(auto node:triangle){const auto row=std::lower_bound(value.owned_global_node_ids.begin(),value.owned_global_node_ids.end(),node)-value.owned_global_node_ids.begin();
			value.owned_reference_lumped_areas_m2[row]+=area/3.;}
	}
	value.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(value);return value;
}
iga::NativeTetMesh SolidMesh(const iga::NativeTetMesh& fluid)
{
	iga::NativeTetMesh solid;solid.points=fluid.points;std::uint64_t cell_id=1,triangle_id=1;
	for(const auto& triangle:fluid.boundary_triangles)if(triangle.boundary_label==7){
		std::array<double,3> apex{};for(auto node:triangle.nodes)for(int c=0;c<3;++c)apex[c]+=fluid.points[node][c]/3.;apex[1]+=.001;
		const auto apex_id=static_cast<std::uint32_t>(solid.points.size());solid.points.push_back(apex);
		auto nodes=std::array<std::uint32_t,4>{{triangle.nodes[0],triangle.nodes[1],triangle.nodes[2],apex_id}};
		if(Determinant(solid,nodes)<0.)std::swap(nodes[1],nodes[2]);
		solid.cells.push_back({cell_id++,nodes});
		solid.boundary_triangles.push_back({triangle_id++,triangle.nodes,7});
		for(const auto pair:{std::array<int,2>{{0,1}},{{1,2}},{{2,0}}})
			solid.boundary_triangles.push_back({triangle_id++,{{triangle.nodes[pair[0]],triangle.nodes[pair[1]],apex_id}},0});
	}
	return solid;
}
iga::DistributedSurfaceInterface Surface(const std::string& domain,const std::string& subsystem,
	const std::string& mesh,bool fluid)
{
	iga::DistributedSurfaceInterface value;value.id={domain,subsystem,"top-wall"};value.subsystem_id=subsystem;
	value.boundary_labels={7};value.reference_mesh_identity_sha256=mesh;
	if(fluid){value.provides={iga::SurfaceFieldQuantity::TractionOnStructure};
		value.requires={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};}
	else{value.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
		value.requires={iga::SurfaceFieldQuantity::TractionOnStructure};}return value;
}
std::array<double,3> VelocityAt(const std::array<double,3>& point)
{
	const double shape_y=4.*point[1]*(height-point[1])/(height*height);
	const double shape_z=4.*point[2]*(width-point[2])/(width*width);
	return {{.002*std::max(0.,shape_y)*std::max(0.,shape_z),0,0}};
}
std::map<std::uint32_t,std::array<double,3>> EssentialVelocity(const iga::NativeTetMesh& mesh)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);std::set<std::uint32_t> nodes;
	for(int label:{0,1})nodes.insert(topology.boundary_velocity_nodes.at(label).begin(),topology.boundary_velocity_nodes.at(label).end());
	std::map<std::uint32_t,std::array<double,3>> result;
	for(auto node:nodes){std::array<double,3> point{};
		if(node<mesh.points.size())point=mesh.points[node];else{const auto& edge=topology.edges[node-mesh.points.size()];
			for(int c=0;c<3;++c)point[c]=.5*(mesh.points[edge[0]][c]+mesh.points[edge[1]][c]);}
		result[node]=VelocityAt(point);
	}return result;
}
std::map<std::size_t,double> SolidConstraints(const iga::NativeTetMesh& solid,const iga::DistributedSurfaceLayout& layout)
{
	std::set<std::uint64_t> interface(layout.owned_global_node_ids.begin(),layout.owned_global_node_ids.end());
	std::map<std::size_t,double> result;
	for(std::size_t node=0;node<solid.points.size();++node){
		const bool unused_or_apex=!interface.count(node);
		const bool edge=interface.count(node)&&(std::abs(solid.points[node][0])<1e-14
			||std::abs(solid.points[node][0]-length)<1e-14||std::abs(solid.points[node][2])<1e-14
			||std::abs(solid.points[node][2]-width)<1e-14);
		if(unused_or_apex||edge)for(int c=0;c<3;++c)result.emplace(3*node+c,0.0);
	}return result;
}
}

int main()
{
	const auto fluid_mesh=FluidMesh();const auto layout=Layout(fluid_mesh);const auto solid_mesh=SolidMesh(fluid_mesh);
	const auto fluid_surface=Surface("fluid","native-channel-ale",layout.reference_mesh_identity_sha256,true);
	const auto solid_surface=Surface("solid","native-channel-solid",layout.reference_mesh_identity_sha256,false);
	const iga::FsiCouplingEdge edge("native-compliant-top",fluid_surface.id,solid_surface.id,
		iga::FsiCouplingLaw::FluidStructureTractionKinematics);
	const auto essential_velocity=EssentialVelocity(fluid_mesh);
	iga::NativeTetAleFsiRuntime fluid("fluid","native-channel-ale",edge,fluid_surface,solid_surface,
		layout,layout,fluid_mesh,7,{1000.,.004},Node(2,0,1),{},0.,1e-7,15,
		essential_velocity,{2});
	iga::NativeTetSolidStaticOptions solid_options;solid_options.relative_tolerance=1e-5;
	solid_options.absolute_tolerance_n=1e-9;solid_options.maximum_iterations=40;
	iga::NativeTetSolidFsiRuntime solid("solid","native-channel-solid",edge,solid_surface,fluid_surface,
		layout,layout,solid_mesh,7,{2.e6,.3,1000.},SolidConstraints(solid_mesh,layout),solid_options);
	iga::StrongFluidStructureCouplingOptions options;options.maximum_iterations=12;
	options.absolute_displacement_tolerance_m=2e-6;options.relative_displacement_tolerance=2e-2;
	options.absolute_traction_residual_tolerance_pa=.5;options.absolute_nodal_force_residual_tolerance_n=1e-4;
	options.fluid_residual_tolerance=1e-7;options.structure_residual_tolerance_n=1e-8;
	options.interface_power_defect_tolerance_w=1e-10;options.aitken_controls.initial_relaxation=.8;
	iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,iga::NativeTetSolidFsiRuntime>
		coordinator(fluid,solid,edge,layout,options);
	const auto result=coordinator.Execute({1,0.,.01});assert(result.converged&&result.iterations>=2);
	const auto checkpoint=iga::CaptureNativeTetFsiCheckpoint(1,fluid,solid);
	double maximum_displacement=0.;for(double value:checkpoint.solid_displacement_m)
		maximum_displacement=std::max(maximum_displacement,std::abs(value));
	const auto bytes=iga::SerializeNativeTetFsiCheckpoint(checkpoint);
	const auto restored_state=iga::ParseNativeTetFsiCheckpoint(bytes);
	auto corrupted=bytes;corrupted[corrupted.size()/2]^=1;bool rejected_corruption=false;
	try{(void)iga::ParseNativeTetFsiCheckpoint(corrupted);}catch(const std::runtime_error&){rejected_corruption=true;}
	assert(rejected_corruption);
	iga::NativeTetAleFsiRuntime restored_fluid("fluid","native-channel-ale",edge,fluid_surface,solid_surface,
		layout,layout,fluid_mesh,7,{1000.,.004},Node(2,0,1),restored_state.fluid_state,
		restored_state.time_s,1e-7,15,essential_velocity,{2},restored_state.ale_displacement_m,
		restored_state.interface_displacement_m);
	iga::NativeTetSolidFsiRuntime restored_solid("solid","native-channel-solid",edge,solid_surface,fluid_surface,
		layout,layout,solid_mesh,7,{2.e6,.3,1000.},SolidConstraints(solid_mesh,layout),solid_options,
		restored_state.solid_displacement_m,restored_state.solid_velocity_m_per_s);
	iga::ValidateRestoredNativeTetFsiCheckpoint(restored_state,restored_fluid,restored_solid);
	restored_fluid.BeginMacroStep({2,.01,.01});bool rejected_active_capture=false;
	try{(void)iga::CaptureNativeTetFsiCheckpoint(1,restored_fluid,restored_solid);}
	catch(const std::runtime_error&){rejected_active_capture=true;}
	assert(rejected_active_capture);restored_fluid.AbortStep();
	iga::StrongFluidStructureCoupling<iga::NativeTetAleFsiRuntime,iga::NativeTetSolidFsiRuntime>
		restored_coordinator(restored_fluid,restored_solid,edge,layout,options);
	const auto continued=coordinator.Execute({2,.01,.01});
	const auto restarted=restored_coordinator.Execute({2,.01,.01});
	assert(continued.converged&&restarted.converged);
	assert(fluid.CommittedCompositionIdentitySha256()==restored_fluid.CommittedCompositionIdentitySha256());
	assert(solid.CommittedStateIdentitySha256()==restored_solid.CommittedStateIdentitySha256());
	for(int step=3;step<=5;++step){
		const double start=fluid.CommittedAleTimeS();
		assert(start==restored_fluid.CommittedAleTimeS());
		const auto uninterrupted=coordinator.Execute({step,start,.01});
		const auto replayed=restored_coordinator.Execute({step,start,.01});
		assert(uninterrupted.converged&&replayed.converged);
		assert(std::abs(fluid.CommittedAleTimeS()-static_cast<double>(step)*.01)<1e-14);
		assert(fluid.CommittedCompositionIdentitySha256()==restored_fluid.CommittedCompositionIdentitySha256());
		assert(solid.CommittedStateIdentitySha256()==restored_solid.CommittedStateIdentitySha256());
		const auto& gate=uninterrupted.history.back();
		assert(gate.displacement_gate_passed&&gate.traction_gate_passed
			&&gate.nodal_force_gate_passed&&gate.fluid_gate_passed
			&&gate.structure_gate_passed&&gate.interface_power_gate_passed);
	}
	assert(maximum_displacement>1e-12);const auto& accepted=result.history.back();
	assert(accepted.displacement_gate_passed&&accepted.traction_gate_passed&&accepted.nodal_force_gate_passed
		&&accepted.fluid_gate_passed&&accepted.structure_gate_passed&&accepted.interface_power_gate_passed);
	std::cout<<"native compliant-channel FSI functional test passed iterations="<<result.iterations
		<<" max_displacement_m="<<maximum_displacement<<" displacement_residual_m="
		<<accepted.area_weighted_rms_residual_m<<" traction_residual_pa="<<*accepted.traction_residual_pa
		<<" nodal_force_residual_n="<<*accepted.nodal_force_residual_n
		<<" fluid_residual="<<*accepted.fluid_acceptance_residual<<" solid_residual_n="
		<<*accepted.structure_acceptance_residual_n<<" power_defect_w="
		<<*accepted.interface_power_defect_w<<" checkpoint_restart_parity=1"
		<<" accepted_steps=5 final_time_s="<<fluid.CommittedAleTimeS()<<"\n";return 0;
}

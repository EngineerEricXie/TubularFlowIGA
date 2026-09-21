#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetVesselTissueSourceMap.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::array<double,3> Coordinate(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,std::uint32_t node)
{
	if(node<mesh.points.size()) return mesh.points[node];
	const auto& edge=topology.edges[node-mesh.points.size()];
	return {{0.5*(mesh.points[edge[0]][0]+mesh.points[edge[1]][0]),
		0.5*(mesh.points[edge[0]][1]+mesh.points[edge[1]][1]),
		0.5*(mesh.points[edge[0]][2]+mesh.points[edge[1]][2])}};
}

iga::NativeTetMesh BuildMatchedOutletTissue(const iga::NativeTetMesh& vessel)
{
	using Face=std::array<std::uint32_t,3>;
	const auto owners=iga::native_tet_matching_detail::BoundaryOwners(vessel);
	iga::NativeTetMesh tissue;
	std::map<Face,int> interface_labels;
	std::uint64_t cell_id=1;
	for(const auto& labels:std::array<std::pair<int,int>,2>{{{2,20},{3,21}}}) {
		std::vector<iga::NativeTetTriangle> cap;
		std::array<double,3> normal_sum{{0.,0.,0.}},centroid_sum{{0.,0.,0.}};
		double total_area=0.;
		for(const auto& triangle:vessel.boundary_triangles) {
			if(triangle.boundary_label!=labels.first) continue;
			cap.push_back(triangle);
			auto face=triangle.nodes;std::sort(face.begin(),face.end());
			const auto owner=owners.find(face);
			if(owner==owners.end())
				throw std::runtime_error("Y outlet facet lacks one vessel owner");
			const auto normal=iga::native_tet_matching_detail::OutwardAreaVector(
				vessel,face,owner->second);
			const double area=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]
				+normal[2]*normal[2]);
			total_area+=area;
			for(int axis=0;axis<3;++axis) {
				normal_sum[axis]+=normal[axis];
				centroid_sum[axis]+=area*(vessel.points[face[0]][axis]
					+vessel.points[face[1]][axis]
					+vessel.points[face[2]][axis])/3.;
			}
		}
		const double normal_length=std::sqrt(normal_sum[0]*normal_sum[0]
			+normal_sum[1]*normal_sum[1]+normal_sum[2]*normal_sum[2]);
		if(cap.empty()||!(total_area>0.)||!(normal_length>0.))
			throw std::runtime_error("Y outlet cap is absent or degenerate");
		std::map<std::uint32_t,std::uint32_t> node_map;
		for(const auto& triangle:cap) for(const auto node:triangle.nodes)
			if(node_map.count(node)==0) {
				const auto next=static_cast<std::uint32_t>(tissue.points.size());
				node_map.emplace(node,next);
				tissue.points.push_back(vessel.points[node]);
			}
		std::array<double,3> apex{};
		for(int axis=0;axis<3;++axis)
			apex[axis]=centroid_sum[axis]/total_area
				+0.004*normal_sum[axis]/normal_length;
		const auto apex_node=static_cast<std::uint32_t>(tissue.points.size());
		tissue.points.push_back(apex);
		for(const auto& triangle:cap) {
			auto vessel_face=triangle.nodes;
			std::sort(vessel_face.begin(),vessel_face.end());
			const auto outward=iga::native_tet_matching_detail::OutwardAreaVector(
				vessel,vessel_face,owners.at(vessel_face));
			double signed_apex_distance=0.;
			for(int axis=0;axis<3;++axis)
				signed_apex_distance+=outward[axis]
					*(apex[axis]-vessel.points[vessel_face[0]][axis]);
			if(!(signed_apex_distance>0.)||!std::isfinite(signed_apex_distance))
				throw std::runtime_error("Y manufactured tissue apex is not outside outlet face");
			iga::NativeTetCell cell;
			cell.id=cell_id++;
			for(std::size_t local=0;local<3;++local)
				cell.nodes[local]=node_map.at(triangle.nodes[local]);
			cell.nodes[3]=apex_node;
			const auto& a=tissue.points[cell.nodes[0]];
			const auto& b=tissue.points[cell.nodes[1]];
			const auto& c=tissue.points[cell.nodes[2]];
			const auto& d=tissue.points[cell.nodes[3]];
			const double det=(b[0]-a[0])*((c[1]-a[1])*(d[2]-a[2])
				-(c[2]-a[2])*(d[1]-a[1]))
				-(b[1]-a[1])*((c[0]-a[0])*(d[2]-a[2])
				-(c[2]-a[2])*(d[0]-a[0]))
				+(b[2]-a[2])*((c[0]-a[0])*(d[1]-a[1])
				-(c[1]-a[1])*(d[0]-a[0]));
			if(!std::isfinite(det)||det==0.)
				throw std::runtime_error("Y matched tissue tetrahedron is degenerate");
			if(det<0.) std::swap(cell.nodes[0],cell.nodes[1]);
			iga::EvaluateNativeTetGeometry(tissue,cell);
			tissue.cells.push_back(cell);
			Face base{{cell.nodes[0],cell.nodes[1],cell.nodes[2]}};
			std::sort(base.begin(),base.end());
			if(!interface_labels.emplace(base,labels.second).second)
				throw std::runtime_error("Y matched tissue has duplicate cap triangle");
		}
	}
	std::map<Face,unsigned> face_uses;
	for(const auto& cell:tissue.cells)
		for(std::size_t opposite=0;opposite<4;++opposite) {
			Face face{};std::size_t next=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[next++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++face_uses[face];
		}
	std::uint64_t face_id=1;
	for(const auto& item:face_uses) {
		if(item.second>2)
			throw std::runtime_error("Y matched tissue is nonmanifold");
		if(item.second!=1)continue;
		const auto label=interface_labels.find(item.first);
		tissue.boundary_triangles.push_back({face_id++,item.first,
			label==interface_labels.end()?30:label->second});
	}
	for(const auto& item:interface_labels)
		if(face_uses.at(item.first)!=1)
			throw std::runtime_error("Y matched tissue cap is not exterior");
	return tissue;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int exit_code=0;
	try {
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		PetscInt requested_steps=5;PetscReal requested_dt=0.01;
		PetscBool steady_mode=PETSC_FALSE;
		PetscOptionsGetInt(nullptr,nullptr,"-bifurcation_steps",&requested_steps,nullptr);
		PetscOptionsGetReal(nullptr,nullptr,"-bifurcation_dt",&requested_dt,nullptr);
		PetscOptionsGetBool(nullptr,nullptr,"-bifurcation_steady",&steady_mode,nullptr);
		if(requested_steps<1||!(requested_dt>0.0))
			throw std::invalid_argument("bifurcation steps and dt must be positive");
		if(argc<2)
			throw std::invalid_argument("usage: native_tet_bifurcation_runtime_test y.msh [result.json]");
		std::ifstream input(argv[1]);if(!input) throw std::runtime_error("cannot open Y mesh");
		const auto mesh=iga::ReadNativeTetMeshGmsh41(input);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		for(const int label:{0,1,2,3})
			if(topology.boundary_velocity_nodes.count(label)==0)
				throw std::runtime_error("Y mesh requires wall/inlet/two-outlet labels 0/1/2/3");
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		const std::size_t dofs=3*velocity_nodes+mesh.points.size();
		const double radius=0.005,maximum_inlet_velocity=0.01,dt=requested_dt;
		std::map<std::uint32_t,std::array<double,3>> boundary;
		for(const auto node:topology.boundary_velocity_nodes.at(1)) {
			const auto point=Coordinate(mesh,topology,node);
			const double profile=std::max(0.0,maximum_inlet_velocity
				*(1.0-(point[1]*point[1]+point[2]*point[2])/(radius*radius)));
			boundary[node]={{profile,0,0}};
		}
		for(const auto node:topology.boundary_velocity_nodes.at(0)) boundary[node]={{0,0,0}};
		std::vector<std::array<double,3>> grid(mesh.points.size(),{{0,0,0}});
		std::vector<double> committed(dofs,0.0);
		iga::NativeTetAleBoundaryConditions outlet_conditions;
		outlet_conditions.prescribed_pressure_pa={{2,0.},{3,0.}};
		iga::NativeTetAlePetscSolveResult result;
		const std::size_t steps=steady_mode?1:static_cast<std::size_t>(requested_steps);
		if(steady_mode)
			result=iga::SolveNativeTetAlePetscSteady(mesh,grid,committed,boundary,
				std::numeric_limits<std::uint32_t>::max(),{1000.0,0.004},
				1e-12,10,outlet_conditions),
				committed=result.replicated_state;
		else for(std::size_t step=0;step<steps;++step)
			result=iga::SolveNativeTetAlePetscTransient(mesh,grid,committed,committed,
				boundary,std::numeric_limits<std::uint32_t>::max(),
				{1000.0,0.004},dt,1e-12,10,outlet_conditions),
				committed=result.replicated_state;
		const auto flows=iga::EvaluateNativeTetBoundaryFlows(mesh,topology,committed);
		const auto tractions=iga::EvaluateNativeTetBoundaryTractions(mesh,topology,committed,0.004);
		const auto recovered_gradients=iga::RecoverNativeTetVelocityGradients(
			mesh,topology,committed);
		const auto recovered_tractions=iga::EvaluateNativeTetBoundaryTractions(
			mesh,topology,committed,0.004,{},&recovered_gradients);
		const auto regular_wall_tractions=iga::EvaluateNativeTetBoundaryTractions(
			mesh,topology,committed,0.004,[](int label,const std::array<double,3>& centroid) {
				return label==0&&centroid[0]>=0.01&&centroid[0]<=0.025;
			});
		const auto recovered_regular_wall_tractions=iga::EvaluateNativeTetBoundaryTractions(
			mesh,topology,committed,0.004,[](int label,const std::array<double,3>& centroid) {
				return label==0&&centroid[0]>=0.01&&centroid[0]<=0.025;
			},&recovered_gradients);
		const auto volume_velocity=iga::EvaluateNativeTetVolumeVelocity(mesh,topology,committed,1000.0);
		const double inlet=flows.at(1).outward_flow_m3_s;
		const double upper=flows.at(2).outward_flow_m3_s;
		const double lower=flows.at(3).outward_flow_m3_s;
		const double wall_flow=flows.at(0).outward_flow_m3_s;
		iga::NativeTetMesh manufactured_tissue;
		manufactured_tissue.points={{{0.,0.,0.}},{{1.,0.,0.}},
			{{0.,1.,0.}},{{0.,0.,1.}}};
		manufactured_tissue.cells.push_back({1,{{0,1,2,3}}});
		manufactured_tissue.boundary_triangles={
			{1,{{1,2,3}},1},{2,{{0,2,3}},2},
			{3,{{0,1,3}},2},{4,{{0,1,2}},2}};
		const auto transferred=iga::MapNativeTetVesselStateToTissueSource(
			mesh,topology,committed,manufactured_tissue,
			{{"upper",2,{{1,1.}}},{"lower",3,{{1,1.}}}});
		const double exchange_flow=upper+lower;
		const double mapped_flow=transferred.tissue_source_s_inv[0]/6.;
		const double exchange_tolerance=1e-12*std::max(1e-12,std::abs(exchange_flow));
		if(std::abs(mapped_flow-exchange_flow)>exchange_tolerance
			||transferred.maximum_port_balance_defect_m3_s>exchange_tolerance)
			throw std::runtime_error("solved vessel outlet-to-tissue source balance failed");
		const auto tissue=iga::SolveNativeTetDarcyPetsc(manufactured_tissue,
			{1.},transferred.tissue_source_s_inv,{{1,0.}});
		if(tissue.converged_reason<=0
			||std::abs(tissue.volume_source_m3_s-exchange_flow)>exchange_tolerance
			||std::abs(tissue.conservative_outward_boundary_flow_m3_s.at(1)
				-exchange_flow)>exchange_tolerance
			||tissue.maximum_cell_balance_defect_m3_s>exchange_tolerance)
			throw std::runtime_error("solved vessel-to-Darcy conservative balance failed");
		const auto matched_tissue=BuildMatchedOutletTissue(mesh);
		const auto matched=iga::MapNativeTetMatchingInterfaceToTissueSource(
			mesh,topology,committed,matched_tissue,
			{{"upper_matching",2,20},{"lower_matching",3,21}});
		const double matched_tolerance=1e-8*std::max(1e-12,std::abs(exchange_flow));
		if(matched.matched_facets.at("upper_matching")==0
			||matched.matched_facets.at("lower_matching")==0
			||std::abs(matched.source.vessel_outward_flow_m3_s.at("upper_matching")
				-upper)>matched_tolerance
			||std::abs(matched.source.vessel_outward_flow_m3_s.at("lower_matching")
				-lower)>matched_tolerance
			||matched.source.maximum_port_balance_defect_m3_s>matched_tolerance)
			throw std::runtime_error("solved Y matching-face map failed its flux gates");
		const auto matched_darcy=iga::SolveNativeTetDarcyPetsc(matched_tissue,
			std::vector<double>(matched_tissue.cells.size(),1.),
			matched.source.tissue_source_s_inv,{{30,0.}});
		if(matched_darcy.converged_reason<=0
			||std::abs(matched_darcy.volume_source_m3_s-exchange_flow)>matched_tolerance
			||std::abs(matched_darcy.conservative_outward_boundary_flow_m3_s.at(30)
				-exchange_flow)>matched_tolerance
			||matched_darcy.maximum_cell_balance_defect_m3_s>matched_tolerance)
			throw std::runtime_error("solved Y matching-face Darcy balance failed");
		const double imbalance=std::abs(inlet+upper+lower)/std::abs(inlet);
		const double split=upper/(upper+lower);
		const auto& wall_traction=tractions.at(0);
		const auto& recovered_wall_traction=recovered_tractions.at(0);
		const auto& regular_wall_traction=regular_wall_tractions.at(0);
		const auto& recovered_regular_wall_traction=recovered_regular_wall_tractions.at(0);
		if(rank==0) std::cout<<"native bifurcation diagnostics residual="
			<<result.final_residual_l2<<" mode="<<(steady_mode?"steady":"backward_euler")
			<<" inlet="<<inlet<<" upper="<<upper
			<<" lower="<<lower<<" wall="<<wall_flow
			<<" split="<<split<<" imbalance="<<imbalance
			<<" wall_mean_wss="<<wall_traction.mean_tangential_traction_pa
			<<" rms_speed="<<volume_velocity.rms_speed_m_s<<'\n';
		if(!(inlet<0.0&&upper>0.0&&lower>0.0&&imbalance<1e-6
			&&split>0.45&&split<0.55&&result.final_residual_l2<1e-8))
			throw std::runtime_error("native Y-bifurcation functional gates failed");
		if(rank==0) {
			std::cout<<"native tetrahedral bifurcation runtime passed ranks="<<ranks
				<<" cells="<<mesh.cells.size()<<" inlet="<<inlet<<" outlets="<<upper<<','<<lower
				<<" split="<<split<<" imbalance="<<imbalance<<'\n';
			if(argc>=3) {
				std::ofstream output(argv[2]);if(!output) throw std::runtime_error("cannot open Y result");
				output<<std::setprecision(17)<<"{\n"
					<<"  \"schema_version\": 1,\n"
					<<"  \"classification\": \"functional_smoke_not_mesh_convergence\",\n"
					<<"  \"backend\": \"native_cpp_petsc_tetrahedral_fem\",\n"
					<<"  \"mpi_ranks\": "<<ranks<<",\n"
					<<"  \"time_integration\": \""
					<<(steady_mode?"steady":"backward_euler")<<"\",\n"
					<<"  \"time_steps\": "<<steps<<",\n"
					<<"  \"dt_s\": "<<(steady_mode?0.0:dt)<<",\n"
					<<"  \"tetrahedra\": "<<mesh.cells.size()<<",\n"
					<<"  \"mixed_dofs\": "<<dofs<<",\n"
					<<"  \"newton_iterations\": "<<result.newton_iterations<<",\n"
					<<"  \"total_linear_iterations\": "<<result.total_linear_iterations<<",\n"
					<<"  \"final_linear_converged_reason\": "
					<<result.final_linear_converged_reason<<",\n"
					<<"  \"linear_solver_type\": \""<<result.linear_solver_type<<"\",\n"
					<<"  \"preconditioner_type\": \""<<result.preconditioner_type<<"\",\n"
					<<"  \"system_scaling_enabled\": "
					<<(result.system_scaling_enabled?"true":"false")<<",\n"
					<<"  \"user_schur_enabled\": "
					<<(result.user_schur_enabled?"true":"false")<<",\n"
					<<"  \"final_residual_l2\": "<<result.final_residual_l2<<",\n"
					<<"  \"inlet_outward_flow_m3_s\": "<<inlet<<",\n"
					<<"  \"upper_outlet_flow_m3_s\": "<<upper<<",\n"
					<<"  \"lower_outlet_flow_m3_s\": "<<lower<<",\n"
					<<"  \"upper_flow_fraction\": "<<split<<",\n"
					<<"  \"relative_mass_imbalance\": "<<imbalance<<",\n"
					<<"  \"prescribed_outlet_pressure_pa\": 0.0,\n"
					<<"  \"mapped_outlet_exchange_flow_m3_s\": "<<exchange_flow<<",\n"
					<<"  \"mapped_tissue_source_flow_m3_s\": "<<mapped_flow<<",\n"
					<<"  \"maximum_port_balance_defect_m3_s\": "
					<<transferred.maximum_port_balance_defect_m3_s<<",\n"
					<<"  \"tissue_darcy_conservative_outward_flow_m3_s\": "
					<<tissue.conservative_outward_boundary_flow_m3_s.at(1)<<",\n"
					<<"  \"tissue_darcy_maximum_cell_balance_defect_m3_s\": "
					<<tissue.maximum_cell_balance_defect_m3_s<<",\n"
					<<"  \"tissue_darcy_linear_converged_reason\": "
					<<tissue.converged_reason<<",\n"
					<<"  \"matched_tissue_tetrahedra\": "<<matched_tissue.cells.size()<<",\n"
					<<"  \"matched_upper_interface_facets\": "
					<<matched.matched_facets.at("upper_matching")<<",\n"
					<<"  \"matched_lower_interface_facets\": "
					<<matched.matched_facets.at("lower_matching")<<",\n"
					<<"  \"matched_tissue_source_flow_m3_s\": "
					<<matched_darcy.volume_source_m3_s<<",\n"
					<<"  \"matched_tissue_conservative_outward_flow_m3_s\": "
					<<matched_darcy.conservative_outward_boundary_flow_m3_s.at(30)<<",\n"
					<<"  \"matched_tissue_maximum_cell_balance_defect_m3_s\": "
					<<matched_darcy.maximum_cell_balance_defect_m3_s<<",\n"
					<<"  \"matched_tissue_linear_converged_reason\": "
					<<matched_darcy.converged_reason<<",\n"
					<<"  \"volume_mean_speed_m_s\": "<<volume_velocity.mean_speed_m_s<<",\n"
					<<"  \"volume_rms_speed_m_s\": "<<volume_velocity.rms_speed_m_s<<",\n"
					<<"  \"kinetic_energy_j\": "<<volume_velocity.kinetic_energy_j<<",\n"
					<<"  \"inlet_average_pressure_pa\": "<<flows.at(1).average_pressure_pa<<",\n"
					<<"  \"upper_average_pressure_pa\": "<<flows.at(2).average_pressure_pa<<",\n"
					<<"  \"lower_average_pressure_pa\": "<<flows.at(3).average_pressure_pa<<",\n"
					<<"  \"wall_mean_wss_pa\": "
					<<wall_traction.mean_tangential_traction_pa<<",\n"
					<<"  \"wall_rms_wss_pa\": "
					<<wall_traction.rms_tangential_traction_pa<<",\n"
					<<"  \"wall_mean_normal_traction_pa\": "
					<<wall_traction.mean_normal_traction_pa<<",\n"
					<<"  \"wall_rms_normal_traction_pa\": "
					<<wall_traction.rms_normal_traction_pa<<",\n"
					<<"  \"wall_mean_total_traction_pa\": "
					<<wall_traction.mean_total_traction_pa<<",\n"
					<<"  \"wall_rms_total_traction_pa\": "
					<<wall_traction.rms_total_traction_pa<<",\n"
					<<"  \"regular_wall_area_m2\": "<<regular_wall_traction.area_m2<<",\n"
					<<"  \"regular_wall_mean_wss_pa\": "
					<<regular_wall_traction.mean_tangential_traction_pa<<",\n"
					<<"  \"regular_wall_rms_wss_pa\": "
					<<regular_wall_traction.rms_tangential_traction_pa<<",\n"
					<<"  \"regular_wall_mean_total_traction_pa\": "
					<<regular_wall_traction.mean_total_traction_pa<<",\n"
					<<"  \"regular_wall_rms_total_traction_pa\": "
					<<regular_wall_traction.rms_total_traction_pa<<",\n"
					<<"  \"recovered_wall_mean_wss_pa\": "
					<<recovered_wall_traction.mean_tangential_traction_pa<<",\n"
					<<"  \"recovered_wall_rms_wss_pa\": "
					<<recovered_wall_traction.rms_tangential_traction_pa<<",\n"
					<<"  \"recovered_wall_mean_total_traction_pa\": "
					<<recovered_wall_traction.mean_total_traction_pa<<",\n"
					<<"  \"recovered_wall_rms_total_traction_pa\": "
					<<recovered_wall_traction.rms_total_traction_pa<<",\n"
					<<"  \"recovered_regular_wall_mean_wss_pa\": "
					<<recovered_regular_wall_traction.mean_tangential_traction_pa<<",\n"
					<<"  \"recovered_regular_wall_rms_wss_pa\": "
					<<recovered_regular_wall_traction.rms_tangential_traction_pa<<",\n"
					<<"  \"recovered_regular_wall_mean_total_traction_pa\": "
					<<recovered_regular_wall_traction.mean_total_traction_pa<<",\n"
					<<"  \"recovered_regular_wall_rms_total_traction_pa\": "
					<<recovered_regular_wall_traction.rms_total_traction_pa<<",\n"
					<<"  \"wall_integrated_traction_n\": ["
					<<wall_traction.integrated_traction_n[0]<<','
					<<wall_traction.integrated_traction_n[1]<<','
					<<wall_traction.integrated_traction_n[2]<<"],\n"
					<<"  \"wall_integrated_tangential_traction_n\": ["
					<<wall_traction.integrated_tangential_traction_n[0]<<','
					<<wall_traction.integrated_tangential_traction_n[1]<<','
					<<wall_traction.integrated_tangential_traction_n[2]<<"],\n"
					<<"  \"upper_integrated_traction_n\": ["
					<<tractions.at(2).integrated_traction_n[0]<<','
					<<tractions.at(2).integrated_traction_n[1]<<','
					<<tractions.at(2).integrated_traction_n[2]<<"],\n"
					<<"  \"lower_integrated_traction_n\": ["
					<<tractions.at(3).integrated_traction_n[0]<<','
					<<tractions.at(3).integrated_traction_n[1]<<','
					<<tractions.at(3).integrated_traction_n[2]<<"]\n}\n";
			}
		}
	} catch(const std::exception& error) {
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0) std::cerr<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}

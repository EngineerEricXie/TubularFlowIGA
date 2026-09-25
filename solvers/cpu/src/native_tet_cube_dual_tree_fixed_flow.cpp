#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetDarcyVisualization.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetLinearElasticPetsc.hpp"
#include "NativeTetLinearElasticVisualization.hpp"
#include "NativeTetMatchingFsiInterface.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetMeshComponents.hpp"
#include "ParallelVtkOutput.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

iga::NativeTetMesh Read(const char* path)
{
	std::ifstream input(path);
	if(!input)throw std::runtime_error(std::string("cannot open mesh: ")+path);
	return iga::ReadNativeTetMeshGmsh41(input);
}

int Label(const char* value)
{
	std::size_t used=0;const int label=std::stoi(value,&used);
	if(used!=std::string(value).size()||label<0)
		throw std::invalid_argument("boundary label is invalid");
	return label;
}

double Positive(const char* value,const char* name)
{
	std::size_t used=0;const double result=std::stod(value,&used);
	if(used!=std::string(value).size()||!std::isfinite(result)||!(result>0.))
		throw std::invalid_argument(std::string(name)+" must be finite and positive");
	return result;
}

std::vector<double> ZeroState(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology)
{
	return std::vector<double>(3*(mesh.points.size()+topology.edges.size())
		+mesh.points.size(),0.);
}

std::map<std::uint32_t,std::array<double,3>> WallVelocity(
	const iga::NativeTaylorHoodTopology& topology,int label)
{
	const auto found=topology.boundary_velocity_nodes.find(label);
	if(found==topology.boundary_velocity_nodes.end())
		throw std::invalid_argument("fluid wall label is absent");
	std::map<std::uint32_t,std::array<double,3>> result;
	for(auto node:found->second)result.emplace(node,std::array<double,3>{{0,0,0}});
	return result;
}

std::map<std::uint32_t,std::array<double,3>> ArterialVelocity(
	const iga::NativeTaylorHoodTopology& topology,int wall,int inlet,double speed)
{
	auto result=WallVelocity(topology,wall);
	const auto found=topology.boundary_velocity_nodes.find(inlet);
	if(found==topology.boundary_velocity_nodes.end())
		throw std::invalid_argument("arterial inlet label is absent");
	std::size_t nonzero=0;
	for(auto node:found->second)
		if(result.emplace(node,std::array<double,3>{{speed,0,0}}).second)++nonzero;
	if(nonzero==0)throw std::invalid_argument("arterial inlet has no free P2 nodes");
	return result;
}

double Flow(const std::map<int,iga::NativeTetBoundaryFlowValue>& values,int label)
{
	const auto found=values.find(label);
	if(found==values.end()||!(found->second.area_m2>0.))
		throw std::runtime_error("fluid boundary label has no positive area");
	return found->second.outward_flow_m3_s;
}

using GeometricFace=std::array<std::array<double,3>,3>;

GeometricFace FaceCoordinates(const iga::NativeTetMesh& mesh,
	const std::array<std::uint32_t,3>& nodes)
{
	GeometricFace result{{mesh.points.at(nodes[0]),mesh.points.at(nodes[1]),
		mesh.points.at(nodes[2])}};
	std::sort(result.begin(),result.end());return result;
}

iga::NativeTetLinearElasticResult WallResponse(
	const iga::NativeTetMesh& reference_fluid,const iga::NativeTetMesh& current_fluid,
	const iga::NativeTaylorHoodTopology& topology,
	const std::vector<double>& state,const iga::NativeTetMesh& wall,
	int fluid_wall_label,int tissue_wall_label,int root_wall_label,
	double viscosity,double young,double poisson)
{
	std::map<GeometricFace,std::array<std::uint32_t,3>> solid_faces;
	std::map<std::uint64_t,std::array<std::uint32_t,3>> fluid_faces;
	std::map<std::uint64_t,double> fluid_areas;
	for(const auto& triangle:wall.boundary_triangles)
		if(triangle.boundary_label==fluid_wall_label){
			if(!solid_faces.emplace(FaceCoordinates(wall,triangle.nodes),
				triangle.nodes).second)
				throw std::runtime_error("duplicate solid interface facet");
		}
	for(const auto& triangle:reference_fluid.boundary_triangles)
		if(triangle.boundary_label==fluid_wall_label){
			fluid_faces.emplace(triangle.id,triangle.nodes);
			const auto&a=current_fluid.points.at(triangle.nodes[0]);
			const auto&b=current_fluid.points.at(triangle.nodes[1]);
			const auto&c=current_fluid.points.at(triangle.nodes[2]);
			std::array<double,3> ab{},ac{};
			for(int axis=0;axis<3;++axis){ab[axis]=b[axis]-a[axis];
				ac[axis]=c[axis]-a[axis];}
			fluid_areas.emplace(triangle.id,.5*iga::NativeTetInterfaceNorm(
				iga::NativeTetInterfaceCross(ab,ac)));
		}
	if(solid_faces.size()!=fluid_faces.size()||solid_faces.empty())
		throw std::runtime_error("solid-fluid interface facet count differs");
	std::vector<double> load(3*wall.points.size(),0.);
	const auto tractions=iga::EvaluateNativeTetFluidTractionOnMatchingStructure(
		current_fluid,topology,state,viscosity,fluid_wall_label);
	for(const auto& item:tractions){
		const auto& fluid_face=fluid_faces.at(item.triangle_id);
		const auto match=solid_faces.find(FaceCoordinates(reference_fluid,fluid_face));
		if(match==solid_faces.end())
			throw std::runtime_error("solid-fluid interface coordinates differ");
		const double area=fluid_areas.at(item.triangle_id);
		for(auto node:match->second)for(int axis=0;axis<3;++axis)
			load[3*node+axis]+=item.traction_on_structure_pa[axis]*area/3.;
	}
	std::map<std::uint32_t,std::array<double,3>> supports;
	for(const auto& triangle:wall.boundary_triangles)
		if(triangle.boundary_label==tissue_wall_label||
			triangle.boundary_label==root_wall_label)
			for(auto node:triangle.nodes)
				supports.emplace(node,std::array<double,3>{{0,0,0}});
	if(supports.empty()||supports.size()==wall.points.size())
		throw std::runtime_error("wall has no valid free/support split");
	return iga::SolveNativeTetLinearElasticPetsc(wall,young,poisson,load,supports);
}

iga::NativeTetLinearElasticResult MeshMotion(
	const iga::NativeTetMesh& fluid,const iga::NativeTetMesh& wall,
	const iga::NativeTetLinearElasticResult& wall_state,int wall_label,
	int root_label,const std::array<int,4>& terminal_labels)
{
	std::map<std::array<double,3>,std::array<double,3>> wall_displacement;
	for(const auto& face:wall.boundary_triangles)
		if(face.boundary_label==wall_label)
			for(auto node:face.nodes)
				wall_displacement.emplace(wall.points.at(node),
					wall_state.displacement_m.at(node));
	std::set<int> fixed_labels{root_label};
	fixed_labels.insert(terminal_labels.begin(),terminal_labels.end());
	std::map<std::uint32_t,std::array<double,3>> prescribed;
	for(const auto& face:fluid.boundary_triangles)
		if(fixed_labels.count(face.boundary_label))
			for(auto node:face.nodes)
				prescribed.emplace(node,std::array<double,3>{{0,0,0}});
	std::size_t moved=0;
	for(const auto& face:fluid.boundary_triangles)
		if(face.boundary_label==wall_label)
			for(auto node:face.nodes){
				const auto match=wall_displacement.find(fluid.points.at(node));
				if(match==wall_displacement.end())
					throw std::runtime_error("fluid-solid wall vertex lacks exact match");
				const auto inserted=prescribed.emplace(node,match->second);
				if(!inserted.second){
					for(int axis=0;axis<3;++axis)
						if(std::abs(inserted.first->second[axis]-match->second[axis])>1e-14)
							throw std::runtime_error("wall movement conflicts with fixed cap");
				}else ++moved;
			}
	if(moved==0)throw std::runtime_error("fluid wall has no moving vertices");
	return iga::SolveNativeTetLinearElasticPetsc(fluid,1.,.3,
		std::vector<double>(3*fluid.points.size(),0.),prescribed);
}

iga::NativeTetMesh Moved(const iga::NativeTetMesh& reference,
	const iga::NativeTetLinearElasticResult& motion)
{
	if(motion.displacement_m.size()!=reference.points.size())
		throw std::invalid_argument("ALE movement shape mismatch");
	auto current=reference;
	for(std::size_t node=0;node<current.points.size();++node)
		for(int axis=0;axis<3;++axis)
			current.points[node][axis]+=motion.displacement_m[node][axis];
	return current;
}

double RelativeDisplacementDifference(const iga::NativeTetLinearElasticResult& before,
	const iga::NativeTetLinearElasticResult& after)
{
	if(before.displacement_m.size()!=after.displacement_m.size())
		throw std::invalid_argument("FSI wall displacement size changed");
	double difference=0.;
	for(std::size_t node=0;node<before.displacement_m.size();++node)
		for(int axis=0;axis<3;++axis)
			difference=std::max(difference,std::abs(
				after.displacement_m[node][axis]-before.displacement_m[node][axis]));
	return difference/std::max(before.maximum_displacement_m,1e-18);
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1,status=0;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	try {
		if(argc!=20&&argc!=21&&argc!=28&&argc!=29)throw std::invalid_argument(
			"usage: native_tet_cube_dual_tree_fixed_flow artery.msh tissue.msh vein.msh "
			"artery_inlet artery_wall artery_tip0 artery_tip1 artery_tip2 artery_tip3 "
			"vein_wall vein_tip0 vein_tip1 vein_tip2 vein_tip3 vein_outlet "
			"inlet_speed_m_s density_kg_m3 viscosity_pa_s mobility_m2_pa_s "
			"[artery_wall.msh vein_wall.msh artery_wall_tissue vein_wall_tissue "
			"artery_root_wall vein_root_wall wall_young_pa wall_poisson_ratio] "
			"[field_output_dir]");
		const auto artery=Read(argv[1]),tissue=Read(argv[2]),vein=Read(argv[3]);
		iga::RequireNativeTetSingleFaceComponent(artery,"arterial fluid");
		iga::RequireNativeTetSingleFaceComponent(vein,"venous fluid");
		iga::RequireNativeTetSingleFaceComponent(tissue,"Darcy tissue");
		const int artery_inlet=Label(argv[4]),artery_wall=Label(argv[5]);
		std::array<int,4> artery_tips{},vein_tips{};
		for(int i=0;i<4;++i)artery_tips[i]=Label(argv[6+i]);
		const int vein_wall=Label(argv[10]);
		for(int i=0;i<4;++i)vein_tips[i]=Label(argv[11+i]);
		const int vein_outlet=Label(argv[15]);
		const double speed=Positive(argv[16],"inlet speed");
		const double density=Positive(argv[17],"fluid density");
		const double viscosity=Positive(argv[18],"fluid viscosity");
		const double mobility=Positive(argv[19],"Darcy mobility");
		const std::array<int,12> labels{{artery_inlet,artery_wall,
			artery_tips[0],artery_tips[1],artery_tips[2],artery_tips[3],
			vein_wall,vein_tips[0],vein_tips[1],vein_tips[2],vein_tips[3],vein_outlet}};
		if(std::set<int>(labels.begin(),labels.end()).size()!=labels.size())
			throw std::invalid_argument("fluid boundary labels must be distinct");
		const auto artery_topology=iga::BuildNativeTaylorHoodTopology(artery);
		const auto artery_bc=ArterialVelocity(artery_topology,artery_wall,
			artery_inlet,speed);
		iga::NativeTetAleBoundaryConditions artery_ports;
		for(int label:artery_tips)artery_ports.prescribed_pressure_pa.emplace(label,0.);
		const auto artery_flow=iga::SolveNativeTetAlePetscSteady(artery,
			std::vector<std::array<double,3>>(artery.points.size(),{{0,0,0}}),
			ZeroState(artery,artery_topology),artery_bc,
			std::numeric_limits<std::uint32_t>::max(),{density,viscosity},
			1e-12,15,artery_ports);
		if(artery_flow.final_linear_converged_reason<=0||
			artery_flow.final_residual_l2>1e-12)
			throw std::runtime_error("arterial native FEM solve did not converge");
		const auto artery_boundary=iga::EvaluateNativeTetBoundaryFlows(artery,
			artery_topology,artery_flow.replicated_state);
		const double inlet_flow=Flow(artery_boundary,artery_inlet);
		std::vector<iga::NativeTetMatchingInterfacePort> ports;
		for(int i=0;i<4;++i)ports.push_back({"arterial_tip_"+std::to_string(i),
			artery_tips[i],artery_tips[i]});
		const auto mapped=iga::MapNativeTetMatchingInterfaceToTissueSource(
			artery,artery_topology,artery_flow.replicated_state,tissue,ports);
		double arterial_exchange=0.;
		for(int i=0;i<4;++i){
			const double port_flow=Flow(artery_boundary,artery_tips[i]);
			const double mapped_flow=mapped.source.vessel_outward_flow_m3_s.at(ports[i].name);
			if(!(port_flow>0.)||std::abs(port_flow-mapped_flow)>
				1e-10*std::max(std::abs(port_flow),1e-18))
				throw std::runtime_error("arterial terminal flux mapping failed");
			arterial_exchange+=port_flow;
		}
		const double scale=std::max(std::abs(inlet_flow),1e-18);
		if(!(inlet_flow<0.)||std::abs(inlet_flow+arterial_exchange)>1e-5*scale)
			throw std::runtime_error("arterial global mass gate failed");
		std::map<int,double> tissue_outlet_pressure;
		for(int label:vein_tips)tissue_outlet_pressure.emplace(label,0.);
		const std::vector<double> mobility_cells(tissue.cells.size(),mobility);
		const auto darcy=iga::SolveNativeTetDarcyPetsc(tissue,mobility_cells,
			mapped.source.tissue_source_s_inv,tissue_outlet_pressure);
		if(darcy.converged_reason<=0||darcy.maximum_cell_balance_defect_m3_s>
			1e-6*scale||std::abs(darcy.volume_source_m3_s-arterial_exchange)>1e-8*scale)
			throw std::runtime_error("Darcy convergence or source gate failed");
		double tissue_to_vein=0.;
		iga::NativeTetAleBoundaryConditions vein_ports;
		vein_ports.prescribed_pressure_pa.emplace(vein_outlet,0.);
		for(int label:vein_tips){
			const double flow=darcy.conservative_outward_boundary_flow_m3_s.at(label);
			if(!(flow>0.))throw std::runtime_error("Darcy venous terminal has reverse flow");
			vein_ports.flow_rate_controls.push_back({label,-flow});
			tissue_to_vein+=flow;
		}
		if(std::abs(tissue_to_vein-arterial_exchange)>1e-6*scale)
			throw std::runtime_error("Darcy global mass gate failed");
		const auto vein_topology=iga::BuildNativeTaylorHoodTopology(vein);
		const auto vein_flow=iga::SolveNativeTetAlePetscSteady(vein,
			std::vector<std::array<double,3>>(vein.points.size(),{{0,0,0}}),
			ZeroState(vein,vein_topology),WallVelocity(vein_topology,vein_wall),
			std::numeric_limits<std::uint32_t>::max(),{density,viscosity},
			1e-12,15,vein_ports);
		if(vein_flow.final_linear_converged_reason<=0||
			vein_flow.final_residual_l2>1e-12)
			throw std::runtime_error("venous native FEM solve did not converge: reason="+
				std::to_string(vein_flow.final_linear_converged_reason)+" residual="+
				std::to_string(vein_flow.final_residual_l2));
		const auto vein_boundary=iga::EvaluateNativeTetBoundaryFlows(vein,
			vein_topology,vein_flow.replicated_state);
		const double vein_exit_flow=Flow(vein_boundary,vein_outlet);
		for(std::size_t i=0;i<vein_tips.size();++i)
			if(std::abs(Flow(vein_boundary,vein_tips[i])-
				vein_ports.flow_rate_controls[i].target_outward_flow_m3_s)>1e-6*scale)
				throw std::runtime_error("venous per-terminal mass gate failed");
		if(!(vein_exit_flow>0.)||std::abs(vein_exit_flow-tissue_to_vein)>1e-5*scale)
			throw std::runtime_error("venous/global mass gate failed");
		std::optional<iga::NativeTetLinearElasticResult> artery_wall_response;
		std::optional<iga::NativeTetLinearElasticResult> vein_wall_response;
		std::optional<iga::NativeTetLinearElasticResult> artery_motion,vein_motion;
		std::optional<iga::NativeTetMesh> moved_artery,moved_vein;
		std::optional<iga::NativeTetAlePetscSolveResult> artery_feedback,vein_feedback;
		std::optional<iga::NativeTetDarcyResult> darcy_feedback;
		std::optional<iga::NativeTetMatchingInterfaceSourceResult> mapped_feedback;
		double fsi_displacement_relative_residual=0.;
		double feedback_arterial_exchange=0.,feedback_tissue_to_vein=0.;
		double feedback_vein_outlet=0.;
		std::array<double,4> feedback_arterial_tips{},feedback_tissue_tips{},
			feedback_venous_tips{};
		std::optional<iga::NativeTetMesh> artery_wall_mesh,vein_wall_mesh;
		if(argc>=28){
			artery_wall_mesh=Read(argv[20]);vein_wall_mesh=Read(argv[21]);
			iga::RequireNativeTetSingleFaceComponent(*artery_wall_mesh,"arterial wall");
			iga::RequireNativeTetSingleFaceComponent(*vein_wall_mesh,"venous wall");
			const int artery_tissue_wall=Label(argv[22]);
			const int vein_tissue_wall=Label(argv[23]);
			const int artery_root_wall=Label(argv[24]);
			const int vein_root_wall=Label(argv[25]);
			const double young=Positive(argv[26],"wall Young modulus");
			const double poisson=Positive(argv[27],"wall Poisson ratio");
			if(!(poisson<.45))throw std::invalid_argument("wall Poisson ratio must be <0.45");
			artery_wall_response=WallResponse(artery,artery,artery_topology,
				artery_flow.replicated_state,*artery_wall_mesh,artery_wall,
				artery_tissue_wall,artery_root_wall,viscosity,young,poisson);
			vein_wall_response=WallResponse(vein,vein,vein_topology,
				vein_flow.replicated_state,*vein_wall_mesh,vein_wall,
				vein_tissue_wall,vein_root_wall,viscosity,young,poisson);
			if(!(artery_wall_response->maximum_displacement_m>0.)||
				!(vein_wall_response->maximum_displacement_m>0.)||
				artery_wall_response->converged_reason<=0||
				vein_wall_response->converged_reason<=0)
				throw std::runtime_error("one-way wall response failed displacement/convergence gate");
			artery_motion=MeshMotion(artery,*artery_wall_mesh,*artery_wall_response,
				artery_wall,artery_inlet,artery_tips);
			vein_motion=MeshMotion(vein,*vein_wall_mesh,*vein_wall_response,
				vein_wall,vein_outlet,vein_tips);
			if(artery_motion->maximum_displacement_m<=0.||
				vein_motion->maximum_displacement_m<=0.||
				artery_motion->minimum_deformation_jacobian<=.5||
				vein_motion->minimum_deformation_jacobian<=.5)
				throw std::runtime_error("FSI ALE movement or Jacobian gate failed");
			moved_artery=Moved(artery,*artery_motion);
			moved_vein=Moved(vein,*vein_motion);
			artery_feedback=iga::SolveNativeTetAlePetscSteady(*moved_artery,
				std::vector<std::array<double,3>>(artery.points.size(),{{0,0,0}}),
				ZeroState(artery,artery_topology),artery_bc,
				std::numeric_limits<std::uint32_t>::max(),{density,viscosity},
				1e-12,15,artery_ports);
			if(artery_feedback->final_linear_converged_reason<=0||
				artery_feedback->final_residual_l2>1e-12)
				throw std::runtime_error("FSI arterial feedback fluid failed to converge: reason="+
					std::to_string(artery_feedback->final_linear_converged_reason)+
					" residual="+std::to_string(artery_feedback->final_residual_l2));
			const auto artery_feedback_boundary=iga::EvaluateNativeTetBoundaryFlows(
				*moved_artery,artery_topology,artery_feedback->replicated_state);
			mapped_feedback=iga::MapNativeTetMatchingInterfaceToTissueSource(
				*moved_artery,artery_topology,artery_feedback->replicated_state,
				tissue,ports);
			for(int i=0;i<4;++i){
				const double flow=Flow(artery_feedback_boundary,artery_tips[i]);
				feedback_arterial_tips[i]=flow;
				if(!(flow>0.)||std::abs(flow-
					mapped_feedback->source.vessel_outward_flow_m3_s.at(ports[i].name))>
					1e-10*scale)
					throw std::runtime_error("FSI arterial terminal mapping failed");
				feedback_arterial_exchange+=flow;
			}
			if(std::abs(Flow(artery_feedback_boundary,artery_inlet)+
				feedback_arterial_exchange)>1e-5*scale)
				throw std::runtime_error("FSI arterial mass gate failed");
			darcy_feedback=iga::SolveNativeTetDarcyPetsc(tissue,mobility_cells,
				mapped_feedback->source.tissue_source_s_inv,tissue_outlet_pressure);
			if(darcy_feedback->converged_reason<=0||
				darcy_feedback->maximum_cell_balance_defect_m3_s>1e-6*scale||
				std::abs(darcy_feedback->volume_source_m3_s-
					feedback_arterial_exchange)>1e-8*scale)
				throw std::runtime_error("FSI Darcy feedback source gate failed");
			iga::NativeTetAleBoundaryConditions vein_feedback_ports;
			vein_feedback_ports.prescribed_pressure_pa.emplace(vein_outlet,0.);
			for(std::size_t i=0;i<vein_tips.size();++i){
				const int label=vein_tips[i];
				const double flow=darcy_feedback->conservative_outward_boundary_flow_m3_s.at(label);
				feedback_tissue_tips[i]=flow;
				if(!(flow>0.))throw std::runtime_error("FSI Darcy venous tip reversed");
				vein_feedback_ports.flow_rate_controls.push_back({label,-flow});
				feedback_tissue_to_vein+=flow;
			}
			if(std::abs(feedback_tissue_to_vein-feedback_arterial_exchange)>1e-6*scale)
				throw std::runtime_error("FSI Darcy global mass gate failed");
			vein_feedback=iga::SolveNativeTetAlePetscSteady(*moved_vein,
				std::vector<std::array<double,3>>(vein.points.size(),{{0,0,0}}),
				ZeroState(vein,vein_topology),WallVelocity(vein_topology,vein_wall),
				std::numeric_limits<std::uint32_t>::max(),{density,viscosity},
				1e-12,15,vein_feedback_ports);
			if(vein_feedback->final_linear_converged_reason<=0||
				vein_feedback->final_residual_l2>1e-12)
				throw std::runtime_error("FSI venous feedback fluid failed to converge");
			const auto vein_feedback_boundary=iga::EvaluateNativeTetBoundaryFlows(
				*moved_vein,vein_topology,vein_feedback->replicated_state);
			for(std::size_t i=0;i<vein_tips.size();++i){
				feedback_venous_tips[i]=-Flow(vein_feedback_boundary,vein_tips[i]);
				if(std::abs(Flow(vein_feedback_boundary,vein_tips[i])-
					vein_feedback_ports.flow_rate_controls[i].target_outward_flow_m3_s)>1e-6*scale)
					throw std::runtime_error("FSI venous per-tip mass gate failed");
			}
			feedback_vein_outlet=Flow(vein_feedback_boundary,vein_outlet);
			if(!(feedback_vein_outlet>0.)||
				std::abs(feedback_vein_outlet-feedback_tissue_to_vein)>1e-5*scale)
				throw std::runtime_error("FSI venous/global mass gate failed");
			const auto artery_second_wall=WallResponse(artery,*moved_artery,
				artery_topology,artery_feedback->replicated_state,*artery_wall_mesh,
				artery_wall,artery_tissue_wall,artery_root_wall,
				viscosity,young,poisson);
			const auto vein_second_wall=WallResponse(vein,*moved_vein,
				vein_topology,vein_feedback->replicated_state,*vein_wall_mesh,
				vein_wall,vein_tissue_wall,vein_root_wall,
				viscosity,young,poisson);
			fsi_displacement_relative_residual=std::max(
				RelativeDisplacementDifference(*artery_wall_response,artery_second_wall),
				RelativeDisplacementDifference(*vein_wall_response,vein_second_wall));
			if(fsi_displacement_relative_residual>1e-3)
				throw std::runtime_error("FSI displacement fixed-point gate failed");
		}
		if(argc==21||argc==29){
			const std::filesystem::path output(argv[argc-1]);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"dual-tree output directory",[&]{
				if(rank!=0)return;
				if(!output.parent_path().empty())
					std::filesystem::create_directories(output.parent_path());
				if(!std::filesystem::create_directory(output))
					throw std::runtime_error("dual-tree output already exists");
			});
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"arterial",
				iga::BuildNativeTetHydraulicVtkPartition(artery,artery,
					artery_flow.replicated_state,viscosity,rank,ranks),0.);
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"tissue",
				iga::BuildNativeTetDarcyVtkPartition(tissue,darcy,mobility_cells,
					mapped.source.tissue_source_s_inv,rank,ranks),0.);
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"venous",
				iga::BuildNativeTetHydraulicVtkPartition(vein,vein,
					vein_flow.replicated_state,viscosity,rank,ranks),0.);
			if(artery_wall_response&&vein_wall_response){
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"arterial_wall",
					iga::BuildNativeTetLinearElasticVtkPartition(*artery_wall_mesh,
						*artery_wall_response,rank,ranks),0.);
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"venous_wall",
					iga::BuildNativeTetLinearElasticVtkPartition(*vein_wall_mesh,
						*vein_wall_response,rank,ranks),0.);
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"arterial_fsi",
					iga::BuildNativeTetHydraulicVtkPartition(artery,*moved_artery,
						artery_feedback->replicated_state,viscosity,rank,ranks),0.);
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"tissue_fsi",
					iga::BuildNativeTetDarcyVtkPartition(tissue,*darcy_feedback,
						mobility_cells,mapped_feedback->source.tissue_source_s_inv,
						rank,ranks),0.);
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"venous_fsi",
					iga::BuildNativeTetHydraulicVtkPartition(vein,*moved_vein,
						vein_feedback->replicated_state,viscosity,rank,ranks),0.);
			}
		}
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native_cube_dual_tree_fixed_flow: PASS artery_inlet_m3_s="<<inlet_flow
			<<" artery_to_tissue_m3_s="<<arterial_exchange
			<<" tissue_to_vein_m3_s="<<tissue_to_vein
			<<" vein_outlet_m3_s="<<vein_exit_flow
			<<" artery_linear_reason="<<artery_flow.final_linear_converged_reason
			<<" darcy_linear_reason="<<darcy.converged_reason
			<<" vein_linear_reason="<<vein_flow.final_linear_converged_reason
			<<" darcy_max_cell_defect_m3_s="<<darcy.maximum_cell_balance_defect_m3_s
			<<" coupled_pressure_continuity=0 vessel_wall_fsi="
				<<(artery_feedback?1:0)
			<<" fsi_displacement_relative_residual="
				<<fsi_displacement_relative_residual
			<<" fsi_artery_to_tissue_m3_s="<<feedback_arterial_exchange
			<<" fsi_tissue_to_vein_m3_s="<<feedback_tissue_to_vein
			<<" fsi_vein_outlet_m3_s="<<feedback_vein_outlet
			<<" fsi_artery_linear_reason="
				<<(artery_feedback?artery_feedback->final_linear_converged_reason:0)
			<<" fsi_darcy_linear_reason="
				<<(darcy_feedback?darcy_feedback->converged_reason:0)
			<<" fsi_vein_linear_reason="
				<<(vein_feedback?vein_feedback->final_linear_converged_reason:0)
			<<" fsi_darcy_max_cell_defect_m3_s="
				<<(darcy_feedback?darcy_feedback->maximum_cell_balance_defect_m3_s:0.)
			<<" fsi_minimum_fluid_jacobian="
				<<(artery_motion&&vein_motion?std::min(
					artery_motion->minimum_deformation_jacobian,
					vein_motion->minimum_deformation_jacobian):0.)
			<<" arterial_wall_max_displacement_m="
				<<(artery_wall_response?artery_wall_response->maximum_displacement_m:0.)
			<<" venous_wall_max_displacement_m="
				<<(vein_wall_response?vein_wall_response->maximum_displacement_m:0.)
			<<" arterial_wall_linear_reason="
				<<(artery_wall_response?artery_wall_response->converged_reason:0)
			<<" venous_wall_linear_reason="
				<<(vein_wall_response?vein_wall_response->converged_reason:0)<<'\n';
		if(rank==0)for(std::size_t i=0;i<4;++i)
			std::cout<<std::setprecision(17)<<"terminal_"<<i
				<<" artery_to_tissue_m3_s="<<Flow(artery_boundary,artery_tips[i])
				<<" tissue_to_vein_m3_s="
				<<darcy.conservative_outward_boundary_flow_m3_s.at(vein_tips[i])
				<<" vein_inward_m3_s="<<-Flow(vein_boundary,vein_tips[i])<<'\n';
		if(rank==0&&artery_feedback)for(std::size_t i=0;i<4;++i)
			std::cout<<std::setprecision(17)<<"terminal_fsi_"<<i
				<<" artery_to_tissue_m3_s="<<feedback_arterial_tips[i]
				<<" tissue_to_vein_m3_s="<<feedback_tissue_tips[i]
				<<" vein_inward_m3_s="<<feedback_venous_tips[i]<<'\n';
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_cube_dual_tree_fixed_flow: ERROR: "
			<<error.what()<<'\n';status=2;
	}
	PetscFinalize();return status;
}

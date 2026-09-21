#ifndef IGA_NATIVE_TET_ALE_GRAPH_PORTS_HPP
#define IGA_NATIVE_TET_ALE_GRAPH_PORTS_HPP

#include "NativeTetAleBackflow.hpp"
#include "NativeTetAleBoundaryControl.hpp"
#include "CoupledDomainRuntime.hpp"
#include "FlowDomainPortMetadata.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Hydraulic part of a native tetra pressure/flow graph adapter. The native
// flow solver uses absolute u at its boundary controller, while graph ports
// exchange the material-domain relative flow (u-w). No FEM assembly is done
// by the graph layer.
class NativeTetAleGraphPorts
{
public:
	NativeTetAleGraphPorts(std::string domain_id,const NativeTetMesh& reference_mesh,
		std::vector<CouplingPort> ports)
		:domain_id_(std::move(domain_id)),ports_(std::move(ports))
	{
		ValidateThreeDBodyFittedFlowDomainMetadata(domain_id_,ports_);
		std::set<int> labels;
		for(const auto& face:reference_mesh.boundary_triangles)
			labels.insert(face.boundary_label);
		for(const auto& port:ports_){
			if(port.orientation.native_to_outward_sign!=1
				||!port.provides.count(PortQuantity::Area)
				||!port.provides.count(PortQuantity::FlowRate)
				||!port.provides.count(PortQuantity::MeanPressure)
				||port.requires.count(PortQuantity::FlowRate)
					+port.requires.count(PortQuantity::MeanPressure)!=1
				||port.requires.count(PortQuantity::MeanNormalTraction))
				throw std::invalid_argument("native tetra graph hydraulic port metadata is unsupported");
			const int label=ParseThreeDFlowBoundaryLabel(port);
			if(!labels.count(label)||!label_to_port_.emplace(label,port.id).second)
				throw std::invalid_argument("native tetra graph hydraulic label is absent or duplicate");
		}
	}

	NativeTetAleBoundaryConditions BoundaryConditions(const NativeTetMesh& current_mesh,
		const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
		const std::map<std::string,PortBoundaryData>& inputs,double time_s) const
	{
		if(!std::isfinite(time_s))
			throw std::invalid_argument("native tetra graph hydraulic time is invalid");
		const auto mesh_flow=MeshFlowByLabel(current_mesh,mesh_velocity_nodes_m_s);
		const auto operators=BuildNativeTetAleBoundaryFluxOperators(current_mesh,
			BuildNativeTaylorHoodTopology(current_mesh));
		NativeTetAleBoundaryConditions result;
		for(const auto& input:inputs)
			if(!HasPort(input.first))
				throw std::invalid_argument("native tetra graph hydraulic input port is unknown");
		for(const auto& port:ports_){
			const auto input=inputs.find(port.id);
			if(input==inputs.end())
				throw std::invalid_argument("native tetra graph hydraulic input is missing");
			ValidatePortBoundaryData(input->second);
			if(std::abs(input->second.time_s-time_s)>
				1e-12*std::max({1.,std::abs(input->second.time_s),std::abs(time_s)})
				||input->second.mean_normal_traction_pa
				||input->second.total_pressure_pa
				||!input->second.concentration.empty()
				||!input->second.outward_species_flux.empty())
				throw std::invalid_argument("native tetra graph hydraulic input has wrong time or quantity");
			const int label=ParseThreeDFlowBoundaryLabel(port);
			if(port.requires.count(PortQuantity::MeanPressure)){
				if(!input->second.mean_pressure_pa||input->second.outward_flow_m3_s)
					throw std::invalid_argument("native tetra graph pressure input is invalid");
				result.prescribed_pressure_pa.emplace(label,*input->second.mean_pressure_pa);
			}else{
				if(!input->second.outward_flow_m3_s||input->second.mean_pressure_pa)
					throw std::invalid_argument("native tetra graph flow input is invalid");
				result.flow_rate_controls.push_back({label,
					*input->second.outward_flow_m3_s+mesh_flow.at(label)});
			}
		}
		ValidateNativeTetAleBoundaryConditions(result,operators);
		return result;
	}

	std::map<std::string,PortState> Observe(const NativeTetMesh& current_mesh,
		const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
		const std::vector<double>& fluid_state,double time_s) const
	{
		if(!std::isfinite(time_s))
			throw std::invalid_argument("native tetra graph hydraulic observation time is invalid");
		const auto topology=BuildNativeTaylorHoodTopology(current_mesh);
		const std::size_t velocity_nodes=current_mesh.points.size()+topology.edges.size();
		if(fluid_state.size()!=3*velocity_nodes+current_mesh.points.size())
			throw std::invalid_argument("native tetra graph hydraulic state size is invalid");
		for(const double value:fluid_state)
			if(!std::isfinite(value))
				throw std::invalid_argument("native tetra graph hydraulic state is nonfinite");
		const auto operators=BuildNativeTetAleBoundaryFluxOperators(current_mesh,topology);
		const auto mesh_flow=MeshFlowByLabel(current_mesh,mesh_velocity_nodes_m_s);
		const auto pressure=MeanPressureByLabel(current_mesh,fluid_state,3*velocity_nodes);
		std::map<std::string,PortState> result;
		for(const auto& item:label_to_port_){
			const auto& operator_value=operators.at(item.first);
			PortState state;state.time_s=time_s;
			state.area_m2=operator_value.area_m2;
			state.outward_flow_m3_s=operator_value.OutwardFlowM3S(fluid_state)
				-mesh_flow.at(item.first);
			state.mean_pressure_pa=pressure.at(item.first);
			ValidatePortState(state);
			result.emplace(item.second,std::move(state));
		}
		return result;
	}

private:
	static std::map<int,double> MeshFlowByLabel(const NativeTetMesh& mesh,
		const std::vector<std::array<double,3>>& mesh_velocity)
	{
		if(mesh_velocity.size()!=mesh.points.size())
			throw std::invalid_argument("native tetra graph mesh velocity size is invalid");
		for(const auto& velocity:mesh_velocity)
			for(const double value:velocity)
				if(!std::isfinite(value))
					throw std::invalid_argument("native tetra graph mesh velocity is nonfinite");
		const auto owners=NativeTetAleBoundaryFaceOwners(mesh);
		std::map<int,double> result;
		for(std::size_t index=0;index<mesh.boundary_triangles.size();++index){
			const auto& face=mesh.boundary_triangles[index];
			const auto& cell=mesh.cells.at(owners[index]);
			const auto& a=mesh.points.at(face.nodes[0]);
			const auto& b=mesh.points.at(face.nodes[1]);
			const auto& c=mesh.points.at(face.nodes[2]);
			const std::array<double,3> u{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> v{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			std::array<double,3> area{{0.5*(u[1]*v[2]-u[2]*v[1]),
				0.5*(u[2]*v[0]-u[0]*v[2]),0.5*(u[0]*v[1]-u[1]*v[0])}};
			std::uint32_t opposite=cell.nodes[0];
			for(const auto node:cell.nodes)
				if(std::find(face.nodes.begin(),face.nodes.end(),node)==face.nodes.end())
					opposite=node;
			double inward=0.;
			for(int axis=0;axis<3;++axis)
				inward+=area[axis]*(mesh.points.at(opposite)[axis]-a[axis]);
			if(inward>0.)for(double& value:area)value=-value;
			std::array<double,3> mean{};
			for(const auto node:face.nodes)
				for(int axis=0;axis<3;++axis)
					mean[axis]+=mesh_velocity.at(node)[axis]/3.;
			for(int axis=0;axis<3;++axis)
				result[face.boundary_label]+=area[axis]*mean[axis];
		}
		return result;
	}

	static std::map<int,double> MeanPressureByLabel(const NativeTetMesh& mesh,
		const std::vector<double>& fluid_state,std::size_t pressure_offset)
	{
		std::map<int,double> area,integral;
		for(const auto& face:mesh.boundary_triangles){
			const auto& a=mesh.points.at(face.nodes[0]);
			const auto& b=mesh.points.at(face.nodes[1]);
			const auto& c=mesh.points.at(face.nodes[2]);
			const std::array<double,3> u{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> v{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			const std::array<double,3> cross{{u[1]*v[2]-u[2]*v[1],
				u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]}};
			const double patch=0.5*std::sqrt(cross[0]*cross[0]
				+cross[1]*cross[1]+cross[2]*cross[2]);
			if(!(patch>0.)||!std::isfinite(patch))
				throw std::invalid_argument("native tetra graph pressure face is invalid");
			area[face.boundary_label]+=patch;
			for(const auto node:face.nodes)
				integral[face.boundary_label]+=patch*fluid_state.at(
					pressure_offset+node)/3.;
		}
		for(auto& item:integral)item.second/=area.at(item.first);
		return integral;
	}

	bool HasPort(const std::string& id) const
	{
		for(const auto& port:ports_)if(port.id==id)return true;
		return false;
	}
	std::string domain_id_;
	std::vector<CouplingPort> ports_;
	std::map<int,std::string> label_to_port_;
};

} // namespace iga

#endif

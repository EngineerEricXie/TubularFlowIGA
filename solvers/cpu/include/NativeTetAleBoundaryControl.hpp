#ifndef IGA_NATIVE_TET_ALE_BOUNDARY_CONTROL_HPP
#define IGA_NATIVE_TET_ALE_BOUNDARY_CONTROL_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetAleFlowRateControl
{
	int boundary_label = 0;
	double target_outward_flow_m3_s = 0.0;
	// Optional affine dependence of the target on the label's area-mean P1
	// pressure: target + slope*mean_pressure. Zero keeps a fixed target.
	double target_pressure_slope_m3_s_pa = 0.0;
};

// Implicit resistive pressure port: the boundary pressure is an extra unknown
// pi with (absolute outward flow - mesh_flow) - (pi - offset)/resistance = 0,
// i.e. pi = offset + resistance*relative outward flow. An implicit 0D RCR step
// reduces to this affine relation, so a monolithic solve needs no outer loop.
struct NativeTetAleResistanceControl
{
	int boundary_label = 0;
	double resistance_pa_s_m3 = 0.0;
	double offset_pa = 0.0;
	double mesh_flow_m3_s = 0.0;
};

struct NativeTetAleBoundaryConditions
{
	std::map<int,double> prescribed_pressure_pa;
	std::vector<NativeTetAleFlowRateControl> flow_rate_controls;
	std::vector<NativeTetAleResistanceControl> resistance_controls;
	// Optional kinetic-energy stabilization on reverse relative flow. This
	// changes the pressure-port model and is never enabled implicitly.
	std::map<int,double> backflow_stabilization_beta;
};

// Sparse, ascending-DOF boundary functionals of one label: the outward flux of
// the P2 velocity and the area-mean of the P1 pressure (pressure node ids).
struct NativeTetAleBoundaryFluxOperator
{
	double area_m2 = 0.0;
	std::vector<std::size_t> velocity_dofs;
	std::vector<double> velocity_coefficients;
	std::vector<std::size_t> pressure_nodes;
	std::vector<double> pressure_weights;

	double OutwardFlowM3S(const std::vector<double>& physical_state) const
	{
		double result=0.0;
		for(std::size_t entry=0;entry<velocity_dofs.size();++entry){
			if(velocity_dofs[entry]>=physical_state.size())
				throw std::invalid_argument("native ALE boundary state is too short");
			result+=velocity_coefficients[entry]*physical_state[velocity_dofs[entry]];
		}
		return result;
	}
};

inline std::map<int,NativeTetAleBoundaryFluxOperator>
BuildNativeTetAleBoundaryFluxOperators(const NativeTetMesh& mesh,
	const NativeTaylorHoodTopology& topology)
{
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	using Face=std::array<std::uint32_t,3>;
	struct Owner { std::size_t cell=0;std::size_t opposite=0; };
	std::map<Face,Owner> owners;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[entry++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			const auto inserted=owners.emplace(face,Owner{cell,opposite});
			if(!inserted.second)owners.erase(inserted.first);
		}
	const std::array<std::array<double,3>,3> quadrature{{
		{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},
		{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
	std::map<int,NativeTetAleBoundaryFluxOperator> result;
	std::map<int,std::map<std::size_t,double>> velocity,pressure;
	for(const auto& triangle:mesh.boundary_triangles){
		auto face=triangle.nodes;std::sort(face.begin(),face.end());
		const auto owner=owners.find(face);
		if(owner==owners.end())
			throw std::runtime_error("native ALE boundary-control face has no unique owner");
		const auto& cell=mesh.cells[owner->second.cell];
		const auto& p0=mesh.points[triangle.nodes[0]];
		const auto& p1=mesh.points[triangle.nodes[1]];
		const auto& p2=mesh.points[triangle.nodes[2]];
		std::array<double,3> first{},second{},area_vector{};
		for(int component=0;component<3;++component){
			first[component]=p1[component]-p0[component];
			second[component]=p2[component]-p0[component];
		}
		area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		const auto& interior=mesh.points[cell.nodes[owner->second.opposite]];
		double inward=0.0;
		for(int component=0;component<3;++component)
			inward+=area_vector[component]*(interior[component]-p0[component]);
		if(inward>0.0)for(double& component:area_vector)component=-component;
		const double area=std::sqrt(area_vector[0]*area_vector[0]
			+area_vector[1]*area_vector[1]+area_vector[2]*area_vector[2]);
		if(!(area>0.0)||!std::isfinite(area))
			throw std::runtime_error("native ALE boundary-control face is degenerate");
		std::array<std::size_t,3> local_vertex{};
		for(std::size_t vertex=0;vertex<3;++vertex){
			local_vertex[vertex]=4;
			for(std::size_t local=0;local<4;++local)
				if(cell.nodes[local]==triangle.nodes[vertex])local_vertex[vertex]=local;
			if(local_vertex[vertex]>=4)
				throw std::runtime_error("native ALE boundary-control vertex is absent from owner");
		}
		auto& output=result[triangle.boundary_label];
		auto& coefficients=velocity[triangle.boundary_label];
		output.area_m2+=area;
		for(const auto node:triangle.nodes)
			pressure[triangle.boundary_label][node]+=area/3.0;
		for(const auto& face_lambda:quadrature){
			std::array<double,4> lambda{{0,0,0,0}};
			for(std::size_t vertex=0;vertex<3;++vertex)
				lambda[local_vertex[vertex]]=face_lambda[vertex];
			const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
			for(std::size_t local=0;local<10;++local){
				const auto node=topology.cell_velocity_nodes[owner->second.cell][local];
				for(int component=0;component<3;++component)
					coefficients[3*node+component]+=basis.velocity[local]*area_vector[component]/3.0;
			}
		}
	}
	for(auto& item:result){
		for(const auto& entry:velocity.at(item.first)){
			if(entry.first>=3*velocity_nodes)
				throw std::runtime_error("native ALE boundary-control DOF is out of range");
			if(entry.second==0.0)continue;
			item.second.velocity_dofs.push_back(entry.first);
			item.second.velocity_coefficients.push_back(entry.second);
		}
		for(const auto& entry:pressure.at(item.first)){
			item.second.pressure_nodes.push_back(entry.first);
			item.second.pressure_weights.push_back(entry.second/item.second.area_m2);
		}
	}
	return result;
}

inline void ValidateNativeTetAleBoundaryConditions(
	const NativeTetAleBoundaryConditions& conditions,
	const std::map<int,NativeTetAleBoundaryFluxOperator>& operators)
{
	std::map<int,bool> controlled;
	for(const auto& pressure:conditions.prescribed_pressure_pa)
		if(!std::isfinite(pressure.second)||operators.count(pressure.first)==0)
			throw std::invalid_argument("native ALE prescribed-pressure boundary is invalid");
	for(const auto& control:conditions.flow_rate_controls){
		if(!std::isfinite(control.target_outward_flow_m3_s)
			||!std::isfinite(control.target_pressure_slope_m3_s_pa)
			||operators.count(control.boundary_label)==0
			||!controlled.emplace(control.boundary_label,true).second)
			throw std::invalid_argument("native ALE flow-rate boundary is invalid");
		if(conditions.prescribed_pressure_pa.count(control.boundary_label))
			throw std::invalid_argument("native ALE boundary cannot prescribe pressure and flow");
	}
	for(const auto& control:conditions.resistance_controls){
		if(!(control.resistance_pa_s_m3>0.0)||!std::isfinite(control.resistance_pa_s_m3)
			||!std::isfinite(control.offset_pa)||!std::isfinite(control.mesh_flow_m3_s)
			||operators.count(control.boundary_label)==0
			||!controlled.emplace(control.boundary_label,true).second
			||conditions.prescribed_pressure_pa.count(control.boundary_label))
			throw std::invalid_argument("native ALE resistance boundary is invalid");
	}
	for(const auto& stabilization:conditions.backflow_stabilization_beta){
		bool resistance=false;
		for(const auto& control:conditions.resistance_controls)
			resistance=resistance||control.boundary_label==stabilization.first;
		if(!(stabilization.second>=0.0)||!std::isfinite(stabilization.second)
			||(conditions.prescribed_pressure_pa.count(stabilization.first)==0&&!resistance))
			throw std::invalid_argument("native ALE backflow stabilization needs a pressure port");
	}
}

} // namespace iga

#endif

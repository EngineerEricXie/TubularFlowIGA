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
};

struct NativeTetAleBoundaryConditions
{
	std::map<int,double> prescribed_pressure_pa;
	std::vector<NativeTetAleFlowRateControl> flow_rate_controls;
	// Optional kinetic-energy stabilization on reverse relative flow. This
	// changes the pressure-port model and is never enabled implicitly.
	std::map<int,double> backflow_stabilization_beta;
};

struct NativeTetAleBoundaryFluxOperator
{
	double area_m2 = 0.0;
	std::vector<double> velocity_coefficients;

	double OutwardFlowM3S(const std::vector<double>& physical_state) const
	{
		if(physical_state.size()<velocity_coefficients.size())
			throw std::invalid_argument("native ALE boundary state is too short");
		double result=0.0;
		for(std::size_t dof=0;dof<velocity_coefficients.size();++dof)
			result+=velocity_coefficients[dof]*physical_state[dof];
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
		if(output.velocity_coefficients.empty())
			output.velocity_coefficients.assign(3*velocity_nodes,0.0);
		output.area_m2+=area;
		for(const auto& face_lambda:quadrature){
			std::array<double,4> lambda{{0,0,0,0}};
			for(std::size_t vertex=0;vertex<3;++vertex)
				lambda[local_vertex[vertex]]=face_lambda[vertex];
			const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
			for(std::size_t local=0;local<10;++local){
				const auto node=topology.cell_velocity_nodes[owner->second.cell][local];
				for(int component=0;component<3;++component)
					output.velocity_coefficients[3*node+component]
						+=basis.velocity[local]*area_vector[component]/3.0;
			}
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
			||operators.count(control.boundary_label)==0
			||!controlled.emplace(control.boundary_label,true).second)
			throw std::invalid_argument("native ALE flow-rate boundary is invalid");
		if(conditions.prescribed_pressure_pa.count(control.boundary_label))
			throw std::invalid_argument("native ALE boundary cannot prescribe pressure and flow");
	}
	for(const auto& stabilization:conditions.backflow_stabilization_beta)
		if(!(stabilization.second>=0.0)||!std::isfinite(stabilization.second)
			||conditions.prescribed_pressure_pa.count(stabilization.first)==0)
			throw std::invalid_argument("native ALE backflow stabilization needs a pressure port");
}

} // namespace iga

#endif

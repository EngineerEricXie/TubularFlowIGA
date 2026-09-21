#ifndef IGA_NATIVE_TET_BOUNDARY_FLOW_HPP
#define IGA_NATIVE_TET_BOUNDARY_FLOW_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetBoundaryFlowValue
{
	double outward_flow_m3_s = 0.0;
	double area_m2 = 0.0;
	double average_pressure_pa = 0.0;
};

struct NativeTetBoundaryTractionValue
{
	double area_m2 = 0.0;
	std::array<double,3> integrated_traction_n{{0,0,0}};
	std::array<double,3> integrated_tangential_traction_n{{0,0,0}};
	double mean_normal_traction_pa = 0.0;
	double rms_normal_traction_pa = 0.0;
	double mean_total_traction_pa = 0.0;
	double rms_total_traction_pa = 0.0;
	double mean_tangential_traction_pa = 0.0;
	double rms_tangential_traction_pa = 0.0;
};

using NativeTetBoundaryTriangleFilter = std::function<bool(int,
	const std::array<double,3>&)>;
using NativeTetVelocityGradient = std::array<std::array<double,3>,3>;

inline std::vector<NativeTetVelocityGradient> RecoverNativeTetVelocityGradients(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& state)
{
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	if(state.size()!=3*velocity_nodes+mesh.points.size())
		throw std::invalid_argument("native velocity-gradient recovery state size is invalid");
	const std::array<std::array<double,3>,10> local_coordinates{{
		{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0.5,0,0}},
		{{0,0.5,0}},{{0,0,0.5}},{{0.5,0.5,0}},{{0.5,0,0.5}},{{0,0.5,0.5}}}};
	std::vector<NativeTetVelocityGradient> recovered(velocity_nodes);
	std::vector<double> weights(velocity_nodes,0.0);
	for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
		const auto geometry=EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		for(std::size_t sample=0;sample<10;++sample) {
			const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(
				geometry,local_coordinates[sample]);
			NativeTetVelocityGradient gradient{};
			for(std::size_t local=0;local<10;++local) {
				const auto global=topology.cell_velocity_nodes[cell][local];
				for(int component=0;component<3;++component)
					for(int derivative=0;derivative<3;++derivative)
						gradient[component][derivative]+=
							basis.velocity_gradients[local][derivative]
							*state[3*global+component];
			}
			const auto node=topology.cell_velocity_nodes[cell][sample];
			weights[node]+=geometry.determinant;
			for(int component=0;component<3;++component)
				for(int derivative=0;derivative<3;++derivative)
					recovered[node][component][derivative]+=
						geometry.determinant*gradient[component][derivative];
		}
	}
	for(std::size_t node=0;node<velocity_nodes;++node) {
		if(!(weights[node]>0.0))
			throw std::runtime_error("native velocity-gradient recovery found an orphan node");
		for(auto& row:recovered[node]) for(auto& value:row) value/=weights[node];
	}
	return recovered;
}

struct NativeTetVolumeVelocityValue
{
	double volume_m3 = 0.0;
	double mean_speed_m_s = 0.0;
	double rms_speed_m_s = 0.0;
	double kinetic_energy_j = 0.0;
};

inline NativeTetVolumeVelocityValue EvaluateNativeTetVolumeVelocity(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& state,double density_kg_m3)
{
	if(!(density_kg_m3>0.0)||!std::isfinite(density_kg_m3))
		throw std::invalid_argument("native volume-velocity density is invalid");
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	if(state.size()!=3*velocity_nodes+mesh.points.size())
		throw std::invalid_argument("native volume-velocity state size is invalid");
	NativeTetVolumeVelocityValue result;double speed_integral=0.0,speed_squared_integral=0.0;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell) {
		const auto geometry=EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		for(const auto& quadrature:NativeTetDegreeFiveQuadrature()) {
			const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(geometry,quadrature.reference);
			std::array<double,3> velocity{{0,0,0}};
			for(std::size_t local=0;local<10;++local) {
				const auto global=topology.cell_velocity_nodes[cell][local];
				for(int component=0;component<3;++component)
					velocity[component]+=basis.velocity[local]*state[3*global+component];
			}
			const double speed_squared=velocity[0]*velocity[0]+velocity[1]*velocity[1]
				+velocity[2]*velocity[2];
			const double weighted_volume=quadrature.weight*geometry.determinant;
			result.volume_m3+=weighted_volume;
			speed_integral+=std::sqrt(speed_squared)*weighted_volume;
			speed_squared_integral+=speed_squared*weighted_volume;
		}
	}
	result.mean_speed_m_s=speed_integral/result.volume_m3;
	result.rms_speed_m_s=std::sqrt(speed_squared_integral/result.volume_m3);
	result.kinetic_energy_j=0.5*density_kg_m3*speed_squared_integral;
	return result;
}

inline std::map<int,NativeTetBoundaryFlowValue> EvaluateNativeTetBoundaryFlows(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& state)
{
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	if(state.size()!=3*velocity_nodes+mesh.points.size())
		throw std::invalid_argument("native boundary-flow state size is invalid");
	using Face=std::array<std::uint32_t,3>;
	struct Owner { std::size_t cell=0;std::size_t opposite=0; };
	std::map<Face,Owner> owners;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		for(std::size_t opposite=0;opposite<4;++opposite) {
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite) face[entry++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			const auto found=owners.find(face);
			if(found==owners.end()) owners.emplace(face,Owner{cell,opposite});
			else owners.erase(found);
		}
	const std::array<std::array<double,3>,3> triangle_points{{
		{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},
		{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
	std::map<int,NativeTetBoundaryFlowValue> result;
	std::map<int,double> pressure_integral;
	for(const auto& triangle:mesh.boundary_triangles) {
		auto face=triangle.nodes;std::sort(face.begin(),face.end());
		const auto owner=owners.find(face);
		if(owner==owners.end())
			throw std::runtime_error("native boundary-flow triangle has no unique owner");
		const auto cell=owner->second.cell;
		const auto& p0=mesh.points[triangle.nodes[0]];
		const auto& p1=mesh.points[triangle.nodes[1]];
		const auto& p2=mesh.points[triangle.nodes[2]];
		std::array<double,3> first{},second{},area_vector{};
		for(int component=0;component<3;++component) {
			first[component]=p1[component]-p0[component];
			second[component]=p2[component]-p0[component];
		}
		area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		const auto& interior=mesh.points[mesh.cells[cell].nodes[owner->second.opposite]];
		double inward=0.0;
		for(int component=0;component<3;++component)
			inward+=area_vector[component]*(interior[component]-p0[component]);
		if(inward>0.0) for(auto& component:area_vector) component=-component;
		const double area=std::sqrt(area_vector[0]*area_vector[0]
			+area_vector[1]*area_vector[1]+area_vector[2]*area_vector[2]);
		std::array<int,3> local_vertex{};
		for(std::size_t vertex=0;vertex<3;++vertex) {
			local_vertex[vertex]=-1;
			for(std::size_t local=0;local<4;++local)
				if(mesh.cells[cell].nodes[local]==triangle.nodes[vertex])
					local_vertex[vertex]=static_cast<int>(local);
			if(local_vertex[vertex]<0)
				throw std::runtime_error("native boundary-flow face is absent from owner");
		}
		auto& value=result[triangle.boundary_label];value.area_m2+=area;
		for(const auto& face_lambda:triangle_points) {
			std::array<double,4> lambda{{0,0,0,0}};
			for(std::size_t vertex=0;vertex<3;++vertex)
				lambda[static_cast<std::size_t>(local_vertex[vertex])]=face_lambda[vertex];
			const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
			std::array<double,3> velocity{{0,0,0}};double pressure=0.0;
			for(std::size_t local=0;local<10;++local) {
				const auto global=topology.cell_velocity_nodes[cell][local];
				for(int component=0;component<3;++component)
					velocity[component]+=basis.velocity[local]*state[3*global+component];
			}
			for(std::size_t local=0;local<4;++local)
				pressure+=basis.pressure[local]
					*state[3*velocity_nodes+mesh.cells[cell].nodes[local]];
			for(int component=0;component<3;++component)
				value.outward_flow_m3_s+=velocity[component]*area_vector[component]/3.0;
			pressure_integral[triangle.boundary_label]+=pressure*area/3.0;
		}
	}
	for(auto& labelled:result)
		labelled.second.average_pressure_pa=pressure_integral[labelled.first]
			/labelled.second.area_m2;
	return result;
}

inline std::map<int,NativeTetBoundaryTractionValue> EvaluateNativeTetBoundaryTractions(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,
	const std::vector<double>& state,double dynamic_viscosity_pa_s,
	const NativeTetBoundaryTriangleFilter& triangle_filter={},
	const std::vector<NativeTetVelocityGradient>* recovered_gradients=nullptr)
{
	if(!(dynamic_viscosity_pa_s>0.0)||!std::isfinite(dynamic_viscosity_pa_s))
		throw std::invalid_argument("native boundary traction viscosity is invalid");
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	if(state.size()!=3*velocity_nodes+mesh.points.size())
		throw std::invalid_argument("native boundary-traction state size is invalid");
	if(recovered_gradients&&recovered_gradients->size()!=velocity_nodes)
		throw std::invalid_argument("native recovered velocity-gradient size is invalid");
	using Face=std::array<std::uint32_t,3>;
	struct Owner { std::size_t cell=0;std::size_t opposite=0; };
	std::map<Face,Owner> owners;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell)
		for(std::size_t opposite=0;opposite<4;++opposite) {
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite) face[entry++]=mesh.cells[cell].nodes[local];
			std::sort(face.begin(),face.end());
			const auto found=owners.find(face);
			if(found==owners.end()) owners.emplace(face,Owner{cell,opposite});
			else owners.erase(found);
		}
	const std::array<std::array<double,3>,3> triangle_points{{
		{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},
		{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
	std::map<int,NativeTetBoundaryTractionValue> result;
	std::map<int,double> normal_integral,normal_squared_integral,total_integral,
		total_squared_integral,tangential_integral,tangential_squared_integral;
	for(const auto& triangle:mesh.boundary_triangles) {
		std::array<double,3> centroid{};
		for(const auto node:triangle.nodes)
			for(int component=0;component<3;++component)
				centroid[component]+=mesh.points[node][component]/3.0;
		if(triangle_filter&&!triangle_filter(triangle.boundary_label,centroid)) continue;
		auto face=triangle.nodes;std::sort(face.begin(),face.end());
		const auto owner=owners.find(face);
		if(owner==owners.end())
			throw std::runtime_error("native boundary-traction triangle has no unique owner");
		const auto cell=owner->second.cell;
		const auto geometry=EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		const auto& p0=mesh.points[triangle.nodes[0]];
		const auto& p1=mesh.points[triangle.nodes[1]];
		const auto& p2=mesh.points[triangle.nodes[2]];
		std::array<double,3> first{},second{},area_vector{},normal{};
		for(int component=0;component<3;++component) {
			first[component]=p1[component]-p0[component];
			second[component]=p2[component]-p0[component];
		}
		area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		const auto& interior=mesh.points[mesh.cells[cell].nodes[owner->second.opposite]];
		double inward=0.0;
		for(int component=0;component<3;++component)
			inward+=area_vector[component]*(interior[component]-p0[component]);
		if(inward>0.0) for(auto& component:area_vector) component=-component;
		const double area=std::sqrt(area_vector[0]*area_vector[0]
			+area_vector[1]*area_vector[1]+area_vector[2]*area_vector[2]);
		for(int component=0;component<3;++component) normal[component]=area_vector[component]/area;
		std::array<int,3> local_vertex{};
		for(std::size_t vertex=0;vertex<3;++vertex) {
			local_vertex[vertex]=-1;
			for(std::size_t local=0;local<4;++local)
				if(mesh.cells[cell].nodes[local]==triangle.nodes[vertex])
					local_vertex[vertex]=static_cast<int>(local);
			if(local_vertex[vertex]<0)
				throw std::runtime_error("native boundary-traction face is absent from owner");
		}
		auto& aggregate=result[triangle.boundary_label];aggregate.area_m2+=area;
		for(const auto& face_lambda:triangle_points) {
			std::array<double,4> lambda{{0,0,0,0}};
			for(std::size_t vertex=0;vertex<3;++vertex)
				lambda[static_cast<std::size_t>(local_vertex[vertex])]=face_lambda[vertex];
			const auto basis=EvaluateNativeTaylorHoodPhysicalBasis(geometry,
				{{lambda[1],lambda[2],lambda[3]}});
			std::array<std::array<double,3>,3> gradient{};double pressure=0.0;
			for(std::size_t local=0;local<10;++local) {
				const auto global=topology.cell_velocity_nodes[cell][local];
				for(int component=0;component<3;++component)
					for(int derivative=0;derivative<3;++derivative)
						gradient[component][derivative]+=
							recovered_gradients?basis.velocity[local]
								*(*recovered_gradients)[global][component][derivative]
							:basis.velocity_gradients[local][derivative]
								*state[3*global+component];
			}
			for(std::size_t local=0;local<4;++local)
				pressure+=basis.pressure[local]
					*state[3*velocity_nodes+mesh.cells[cell].nodes[local]];
			std::array<double,3> traction{};
			for(int row=0;row<3;++row)
				for(int column=0;column<3;++column)
					traction[row]+=((row==column?-pressure:0.0)
						+dynamic_viscosity_pa_s*(gradient[row][column]
							+gradient[column][row]))*normal[column];
			double normal_traction=0.0;
			for(int component=0;component<3;++component)
				normal_traction+=traction[component]*normal[component];
			double tangential_squared=0.0,total_squared=0.0;
			for(int component=0;component<3;++component) {
				const double tangential=traction[component]-normal_traction*normal[component];
				tangential_squared+=tangential*tangential;
				total_squared+=traction[component]*traction[component];
				aggregate.integrated_traction_n[component]+=traction[component]*area/3.0;
				aggregate.integrated_tangential_traction_n[component]+=tangential*area/3.0;
			}
			const double magnitude=std::sqrt(tangential_squared);
			normal_integral[triangle.boundary_label]+=normal_traction*area/3.0;
			normal_squared_integral[triangle.boundary_label]+=
				normal_traction*normal_traction*area/3.0;
			total_integral[triangle.boundary_label]+=std::sqrt(total_squared)*area/3.0;
			total_squared_integral[triangle.boundary_label]+=total_squared*area/3.0;
			tangential_integral[triangle.boundary_label]+=magnitude*area/3.0;
			tangential_squared_integral[triangle.boundary_label]+=tangential_squared*area/3.0;
		}
	}
	for(auto& labelled:result) {
		labelled.second.mean_normal_traction_pa=
			normal_integral[labelled.first]/labelled.second.area_m2;
		labelled.second.rms_normal_traction_pa=std::sqrt(
			normal_squared_integral[labelled.first]/labelled.second.area_m2);
		labelled.second.mean_total_traction_pa=
			total_integral[labelled.first]/labelled.second.area_m2;
		labelled.second.rms_total_traction_pa=std::sqrt(
			total_squared_integral[labelled.first]/labelled.second.area_m2);
		labelled.second.mean_tangential_traction_pa=
			tangential_integral[labelled.first]/labelled.second.area_m2;
		labelled.second.rms_tangential_traction_pa=std::sqrt(
			tangential_squared_integral[labelled.first]/labelled.second.area_m2);
	}
	return result;
}

} // namespace iga

#endif

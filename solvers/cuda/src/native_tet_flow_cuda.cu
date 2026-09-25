#include "CudaRuntime.hpp"
#include "NativeTetFlowKernels.cuh"
#include "NativeTetAleTransient.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "PartitionedVtkOutput.hpp"
#include "VtkOutput.hpp"
#include "CaseConfig.hpp"
#include "ZeroDFlowDomain.hpp"

#include <cusolverDn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>

#include <sys/resource.h>

namespace {

std::vector<iga::cuda::NativeTetFlowQuadrature> Quadrature(
	const iga::NativeTetMesh& mesh)
{
	std::vector<iga::cuda::NativeTetFlowQuadrature> result;
	result.reserve(mesh.cells.size()*iga::cuda::kTetFlowQuadrature);
	for(const auto& cell:mesh.cells){
		const auto geometry=iga::EvaluateNativeTetGeometry(mesh,cell);
		for(const auto& point:iga::NativeTetDegreeFiveQuadrature()){
			const auto basis=iga::EvaluateNativeTaylorHoodPhysicalBasis(
				geometry,point.reference);
			iga::cuda::NativeTetFlowQuadrature item{};
			item.weight=point.weight*geometry.determinant;
			for(int node=0;node<10;++node){
				item.velocity[node]=basis.velocity[node];
				for(int axis=0;axis<3;++axis)
					item.gradient[node][axis]=basis.velocity_gradients[node][axis];
			}
			for(int node=0;node<4;++node){
				item.pressure[node]=basis.pressure[node];
				for(int axis=0;axis<3;++axis)
					item.position[axis]+=basis.pressure[node]
						*mesh.points[cell.nodes[node]][axis];
			}
			result.push_back(item);
		}
	}
	return result;
}

std::vector<int> ElementRows(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology)
{
	std::vector<int> rows;
	rows.reserve(mesh.cells.size()*iga::cuda::kTetFlowDofs);
	const int pressure_start=3*static_cast<int>(mesh.points.size()+topology.edges.size());
	for(std::size_t cell=0;cell<mesh.cells.size();++cell){
		for(const auto node:topology.cell_velocity_nodes[cell])
			for(int axis=0;axis<3;++axis)rows.push_back(3*static_cast<int>(node)+axis);
		for(const auto node:mesh.cells[cell].nodes)
			rows.push_back(pressure_start+static_cast<int>(node));
	}
	return rows;
}

struct FlowCase
{
	struct T7Outlet
	{
		int boundary_label=-1;
		iga::ZeroDFlowModel model;
		iga::ZeroDFlowState state;
	};
	std::filesystem::path mesh_file;
	double density=0.,viscosity=0.,dt_s=0.,tolerance=0.;
	int steps=0,maximum_iterations=0;
	std::array<double,3> initial_velocity{};
	std::map<int,std::array<double,3>> velocity_by_label;
	std::map<int,double> pressure_by_label;
	bool t7_enabled=false;
	int t7_inlet_label=-1;
	iga::ZeroDFlowModel t7_source;
	iga::ZeroDFlowState t7_source_state;
	std::vector<T7Outlet> t7_outlets;
};

std::string TextFile(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if(!input)throw std::runtime_error("cannot open FEM CUDA input: "+path.string());
	return iga::ReadCheckedText(input);
}

std::array<double,3> Vector3(const iga::config_detail::JsonValue& value,
	const std::string& context)
{
	const auto& array=iga::config_detail::RequireArray(value,context);
	if(array.size()!=3)throw std::invalid_argument(context+" requires three components");
	std::array<double,3> result{};
	for(int axis=0;axis<3;++axis)
		result[axis]=iga::config_detail::RequireNumber(array[axis],context);
	return result;
}

FlowCase ParseCase(const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	const auto parsed=JsonParser(TextFile(path)).Parse();
	const auto& root=RequireObject(parsed,"FEM CUDA case");
	auto required=[](const auto& object,const char* key)->const JsonValue& {
		const auto* item=Find(object,key);
		if(!item)throw std::invalid_argument(std::string("FEM CUDA case requires ")+key);
		return *item;
	};
	if(RequireInteger(required(root,"schema_version"),"schema_version")!=1)
		throw std::invalid_argument("unsupported FEM CUDA case schema");
	FlowCase result;
	result.mesh_file=path.parent_path()/RequireString(
		required(root,"mesh_file"),"mesh_file");
	const auto& fluid=RequireObject(required(root,"fluid"),"fluid");
	result.density=RequireNumber(required(fluid,"density_kg_m3"),"density_kg_m3");
	result.viscosity=RequireNumber(required(fluid,"dynamic_viscosity_pa_s"),
		"dynamic_viscosity_pa_s");
	result.tolerance=RequireNumber(required(fluid,"nonlinear_tolerance"),
		"nonlinear_tolerance");
	result.maximum_iterations=RequireInteger(required(fluid,"maximum_iterations"),
		"maximum_iterations");
	if(const auto* initial=Find(fluid,"initial_velocity_m_s"))
		result.initial_velocity=Vector3(*initial,"initial_velocity_m_s");
	const auto& boundaries=RequireObject(required(root,"boundaries"),"boundaries");
	const auto& velocities=RequireObject(required(boundaries,"velocity_by_label_m_s"),
		"velocity_by_label_m_s");
	for(const auto& item:velocities)
		result.velocity_by_label.emplace(std::stoi(item.first),
			Vector3(item.second,"velocity_by_label_m_s"));
	const auto& pressures=RequireObject(required(boundaries,"pressure_by_label_pa"),
		"pressure_by_label_pa");
	for(const auto& item:pressures)
		result.pressure_by_label.emplace(std::stoi(item.first),
			RequireNumber(item.second,"pressure_by_label_pa"));
	const auto& time=RequireObject(required(root,"time"),"time");
	result.dt_s=RequireNumber(required(time,"dt_s"),"dt_s");
	result.steps=RequireInteger(required(time,"steps"),"steps");
	if(const auto* ports=Find(root,"t7_ports")){
		const auto& object=RequireObject(*ports,"t7_ports");
		result.t7_enabled=true;
		result.t7_inlet_label=RequireInteger(required(object,"inlet_boundary_label"),
			"t7_ports.inlet_boundary_label");
		const auto& source=RequireObject(required(object,"source"),"t7_ports.source");
		result.t7_source.role=iga::ZeroDFlowRole::SourceReservoir;
		result.t7_source.source={
			RequireNumber(required(source,"capacitance_m3_pa"),"source capacitance"),
			RequireNumber(required(source,"resistance_pa_s_m3"),"source resistance"),
			RequireNumber(required(source,"prescribed_flow_m3_s"),"source prescribed flow")};
		result.t7_source_state.stored_pressure_pa=
			RequireNumber(required(source,"initial_pressure_pa"),"source initial pressure");
		for(const auto& item:RequireArray(required(object,"outlets"),"t7_ports.outlets")){
			const auto& terminal=RequireObject(item,"t7_ports.outlets[]");
			FlowCase::T7Outlet outlet;
			outlet.boundary_label=RequireInteger(required(terminal,"boundary_label"),
				"terminal boundary label");
			outlet.model.role=iga::ZeroDFlowRole::TerminalRcr;
			outlet.model.terminal={
				RequireNumber(required(terminal,"proximal_resistance_pa_s_m3"),
					"terminal proximal resistance"),
				RequireNumber(required(terminal,"distal_resistance_pa_s_m3"),
					"terminal distal resistance"),
				RequireNumber(required(terminal,"capacitance_m3_pa"),"terminal capacitance"),
				RequireNumber(required(terminal,"distal_pressure_pa"),"terminal distal pressure")};
			outlet.state.stored_pressure_pa=RequireNumber(
				required(terminal,"initial_pressure_pa"),"terminal initial pressure");
			result.t7_outlets.push_back(outlet);
		}
		iga::ValidateZeroDFlowModel(result.t7_source);
		iga::ValidateZeroDFlowState(result.t7_source_state);
		if(result.velocity_by_label.count(result.t7_inlet_label)!=1
			||result.t7_outlets.empty())
			throw std::invalid_argument("FEM CUDA T7 inlet or outlet set is invalid");
		for(const auto& outlet:result.t7_outlets){
			iga::ValidateZeroDFlowModel(outlet.model);
			iga::ValidateZeroDFlowState(outlet.state);
			if(result.pressure_by_label.count(outlet.boundary_label)!=1)
				throw std::invalid_argument("FEM CUDA T7 outlet lacks a pressure boundary");
		}
	}
	if(!(result.density>0.)||!(result.viscosity>0.)||!(result.dt_s>0.)
		||!(result.tolerance>0.)||result.maximum_iterations<1||result.steps<1
		||result.velocity_by_label.empty()||result.pressure_by_label.empty())
		throw std::invalid_argument("invalid FEM CUDA flow case parameters");
	return result;
}

std::map<int,std::array<double,3>> BoundaryAreaVectors(const iga::NativeTetMesh& mesh)
{
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,std::uint32_t> opposite;
	for(const auto& cell:mesh.cells)
		for(int excluded=0;excluded<4;++excluded){
			Face face{};int count=0;
			for(int local=0;local<4;++local)
				if(local!=excluded)face[count++]=cell.nodes[local];
			std::sort(face.begin(),face.end());
			if(!opposite.emplace(face,cell.nodes[excluded]).second)opposite.erase(face);
		}
	std::map<int,std::array<double,3>> result;
	for(const auto& face:mesh.boundary_triangles){
		auto key=face.nodes;std::sort(key.begin(),key.end());
		const auto& inside=mesh.points[opposite.at(key)];
		const auto& a=mesh.points[face.nodes[0]];
		const auto& b=mesh.points[face.nodes[1]];
		const auto& c=mesh.points[face.nodes[2]];
		std::array<double,3> first{},second{},area{};
		for(int axis=0;axis<3;++axis){
			first[axis]=b[axis]-a[axis];second[axis]=c[axis]-a[axis];
		}
		area={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		double inward=0.;
		for(int axis=0;axis<3;++axis)inward+=area[axis]*(inside[axis]-a[axis]);
		if(inward>0.)for(auto& value:area)value=-value;
		auto& total=result[face.boundary_label];
		for(int axis=0;axis<3;++axis)total[axis]+=area[axis];
	}
	return result;
}

std::array<double,3> UniformVelocityForFlow(const std::array<double,3>& area,
	double outward_flow_m3_s)
{
	const double squared=area[0]*area[0]+area[1]*area[1]+area[2]*area[2];
	std::array<double,3> result{};
	for(int axis=0;axis<3;++axis)result[axis]=outward_flow_m3_s*area[axis]/squared;
	return result;
}

std::vector<double> BoundaryLoad(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,
	const std::map<int,double>& pressure,int dofs)
{
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,std::uint32_t> opposite;
	for(const auto& cell:mesh.cells)
		for(int excluded=0;excluded<4;++excluded){
			Face face{};int count=0;
			for(int local=0;local<4;++local)
				if(local!=excluded)face[count++]=cell.nodes[local];
			std::sort(face.begin(),face.end());
			if(!opposite.emplace(face,cell.nodes[excluded]).second)
				opposite.erase(face);
		}
	std::map<std::array<std::uint32_t,2>,std::uint32_t> edge_ids;
	for(std::size_t index=0;index<topology.edges.size();++index)
		edge_ids.emplace(topology.edges[index],
			static_cast<std::uint32_t>(mesh.points.size()+index));
	std::vector<double> load(dofs,0.);
	for(const auto& face:mesh.boundary_triangles){
		const auto condition=pressure.find(face.boundary_label);
		if(condition==pressure.end())continue;
		auto key=face.nodes;std::sort(key.begin(),key.end());
		const auto& inside=mesh.points[opposite.at(key)];
		const auto& a=mesh.points[face.nodes[0]];
		const auto& b=mesh.points[face.nodes[1]];
		const auto& c=mesh.points[face.nodes[2]];
		std::array<double,3> first{},second{},area_vector{};
		for(int axis=0;axis<3;++axis){
			first[axis]=b[axis]-a[axis];second[axis]=c[axis]-a[axis];
		}
		area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),
			0.5*(first[2]*second[0]-first[0]*second[2]),
			0.5*(first[0]*second[1]-first[1]*second[0])}};
		double inward=0.;
		for(int axis=0;axis<3;++axis)
			inward+=area_vector[axis]*(inside[axis]-a[axis]);
		if(inward>0.)for(auto& value:area_vector)value=-value;
		const std::array<std::array<int,2>,3> pairs{{{{0,1}},{{0,2}},{{1,2}}}};
		for(const auto pair:pairs){
			std::array<std::uint32_t,2> edge{{face.nodes[pair[0]],face.nodes[pair[1]]}};
			std::sort(edge.begin(),edge.end());
			const auto node=edge_ids.at(edge);
			for(int axis=0;axis<3;++axis)
				load[3*node+axis]+=condition->second*area_vector[axis]/3.;
		}
	}
	return load;
}

void VerifyElement(const std::filesystem::path& mesh_path)
{
	std::ifstream input(mesh_path);
	if(!input)throw std::runtime_error("cannot open FEM mesh");
	const auto mesh=iga::ReadNativeTetMeshGmsh41(input);
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const auto quadrature=Quadrature(mesh);
	const auto rows=ElementRows(mesh,topology);
	const int dofs=3*static_cast<int>(mesh.points.size()+topology.edges.size())
		+static_cast<int>(mesh.points.size());
	std::vector<double> state(dofs),previous(dofs);
	for(int index=0;index<dofs;++index){
		state[index]=0.03*std::sin(0.3*index);
		previous[index]=0.02*std::cos(0.2*index);
	}
	std::vector<double> velocity(mesh.cells.size()*12);
	for(std::size_t index=0;index<velocity.size();++index)
		velocity[index]=0.001*std::sin(0.5*index);
	const std::array<double,12> acceleration{{0.02,0.1,-0.2,0.3,
		-0.01,0.04,0.05,-0.06,0.03,-0.02,0.01,0.07}};
	iga::NativeNavierStokesParameters parameters;
	parameters.density=1000.;parameters.dynamic_viscosity=0.004;
	for(int component=0;component<3;++component)
		for(int coefficient=0;coefficient<4;++coefficient)
			parameters.body_acceleration_m_s2[component][coefficient]=
				acceleration[4*component+coefficient];
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowQuadrature> device_quadrature(
		quadrature.size());
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowField> fields(quadrature.size());
	iga::cuda::DeviceBuffer<int> device_rows(rows.size());
	iga::cuda::DeviceBuffer<double> device_state(state.size()),device_previous(previous.size());
	iga::cuda::DeviceBuffer<double> device_velocity(velocity.size());
	iga::cuda::DeviceBuffer<double> device_acceleration(acceleration.size());
	iga::cuda::DeviceBuffer<double> device_residual(mesh.cells.size()*34);
	iga::cuda::DeviceBuffer<double> device_jacobian(mesh.cells.size()*34*34);
	device_quadrature.CopyFromHost(quadrature.data(),quadrature.size());
	device_rows.CopyFromHost(rows.data(),rows.size());
	device_state.CopyFromHost(state.data(),state.size());
	device_previous.CopyFromHost(previous.data(),previous.size());
	device_velocity.CopyFromHost(velocity.data(),velocity.size());
	device_acceleration.CopyFromHost(acceleration.data(),acceleration.size());
	const int field_count=static_cast<int>(quadrature.size());
	const int entry_count=static_cast<int>(mesh.cells.size())
		*iga::cuda::kTetFlowEntries;
	iga::cuda::EvaluateNativeTetFlowFields<<<(field_count+127)/128,128>>>(
		device_quadrature.data(),device_rows.data(),device_state.data(),
		device_previous.data(),device_velocity.data(),
		static_cast<int>(mesh.cells.size()),fields.data());
	iga::cuda::CheckKernel("FEM CUDA field evaluation");
	iga::cuda::AssembleNativeTetFlowElements<<<(entry_count+127)/128,128>>>(
		device_quadrature.data(),fields.data(),device_acceleration.data(),
		parameters.density,parameters.dynamic_viscosity,0.05,
		static_cast<int>(mesh.cells.size()),device_residual.data(),device_jacobian.data());
	iga::cuda::CheckKernel("FEM CUDA element assembly");
	std::vector<double> residual(device_residual.size()),jacobian(device_jacobian.size());
	device_residual.CopyToHost(residual.data(),residual.size());
	device_jacobian.CopyToHost(jacobian.data(),jacobian.size());
	double max_residual_error=0.,max_jacobian_error=0.;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell){
		std::array<double,34> current{},committed{};
		std::array<std::array<double,3>,4> grid{};
		for(int local=0;local<34;++local){
			current[local]=state[rows[cell*34+local]];
			committed[local]=previous[rows[cell*34+local]];
		}
		for(int node=0;node<4;++node)
			for(int axis=0;axis<3;++axis)
				grid[node][axis]=velocity[cell*12+3*node+axis];
		const auto expected=iga::BuildNativeTaylorHoodAleTransientElement(
			mesh,mesh.cells[cell],current,committed,grid,parameters,0.05);
		for(int local=0;local<34;++local)
			max_residual_error=std::max(max_residual_error,
				std::abs(residual[cell*34+local]-expected.residual[local]));
		for(int entry=0;entry<34*34;++entry)
			max_jacobian_error=std::max(max_jacobian_error,
				std::abs(jacobian[cell*34*34+entry]-expected.jacobian[entry]));
	}
	std::cout<<"native_tet_flow_cuda_element cells="<<mesh.cells.size()
		<<" max_residual_error="<<max_residual_error
		<<" max_jacobian_error="<<max_jacobian_error<<'\n';
	if(max_residual_error>1e-9||max_jacobian_error>1e-9)
		throw std::runtime_error("FEM CUDA element differs from CPU P2/P1 element");
}

void CheckSolver(cusolverStatus_t status,const char* operation)
{
	if(status!=CUSOLVER_STATUS_SUCCESS)
		throw std::runtime_error(std::string(operation)+" failed with status "
			+std::to_string(static_cast<int>(status)));
}

class DenseSolver
{
public:
	DenseSolver(int dofs,double* matrix)
	{
		CheckSolver(cusolverDnCreate(&handle_),"cusolverDnCreate");
		int workspace=0;
		CheckSolver(cusolverDnDgetrf_bufferSize(handle_,dofs,dofs,matrix,dofs,
			&workspace),"cusolverDnDgetrf_bufferSize");
		work_.Allocate(workspace);
		pivots_.Allocate(dofs);
		info_.Allocate(1);
	}
	~DenseSolver(){if(handle_)cusolverDnDestroy(handle_);}
	void Solve(int dofs,double* matrix,double* right_hand_side)
	{
		CheckSolver(cusolverDnDgetrf(handle_,dofs,dofs,matrix,dofs,work_.data(),
			pivots_.data(),info_.data()),"cusolverDnDgetrf");
		int info=0;info_.CopyToHost(&info,1);
		if(info!=0)throw std::runtime_error("FEM CUDA Jacobian factorization failed");
		CheckSolver(cusolverDnDgetrs(handle_,CUBLAS_OP_N,dofs,1,matrix,dofs,
			pivots_.data(),right_hand_side,dofs,info_.data()),
			"cusolverDnDgetrs");
		info_.CopyToHost(&info,1);
		if(info!=0)throw std::runtime_error("FEM CUDA Newton solve failed");
	}
private:
	cusolverDnHandle_t handle_=nullptr;
	iga::cuda::DeviceBuffer<double> work_;
	iga::cuda::DeviceBuffer<int> pivots_,info_;
};

void RunCase(const std::filesystem::path& case_path,
	const std::filesystem::path& output)
{
	auto specification=ParseCase(case_path);
	std::ifstream input(specification.mesh_file);
	if(!input)throw std::runtime_error("cannot open FEM CUDA tetra mesh");
	const auto mesh=iga::ReadNativeTetMeshGmsh41(input);
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const auto quadrature=Quadrature(mesh);
	const auto rows=ElementRows(mesh,topology);
	const int cells=static_cast<int>(mesh.cells.size());
	const int velocity_nodes=static_cast<int>(mesh.points.size()+topology.edges.size());
	const int dofs=3*velocity_nodes+static_cast<int>(mesh.points.size());
	std::vector<double> state(dofs,0.);
	for(int node=0;node<velocity_nodes;++node)
		for(int axis=0;axis<3;++axis)
			state[3*node+axis]=specification.initial_velocity[axis];
	const auto area_vectors=BoundaryAreaVectors(mesh);
	std::vector<int> constrained(dofs,0);
	std::vector<double> prescribed(dofs,0.);
	for(const auto& boundary:specification.velocity_by_label){
		const auto& nodes=topology.boundary_velocity_nodes.at(boundary.first);
		for(const auto node:nodes)
			for(int axis=0;axis<3;++axis){
				const int row=3*static_cast<int>(node)+axis;
				constrained[row]=1;
				prescribed[row]=boundary.second[axis];
				state[row]=boundary.second[axis];
			}
	}
	auto boundary_load=BoundaryLoad(mesh,topology,
		specification.pressure_by_label,dofs);
	std::vector<double> mesh_velocity(cells*12,0.);
	std::array<double,12> acceleration{};
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowQuadrature> device_quadrature(
		quadrature.size());
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetFlowField> fields(quadrature.size());
	iga::cuda::DeviceBuffer<int> device_rows(rows.size()),device_constrained(dofs);
	iga::cuda::DeviceBuffer<double> device_state(dofs),device_previous(dofs);
	iga::cuda::DeviceBuffer<double> device_mesh_velocity(mesh_velocity.size());
	iga::cuda::DeviceBuffer<double> device_acceleration(acceleration.size());
	iga::cuda::DeviceBuffer<double> device_element_residual(cells*34);
	iga::cuda::DeviceBuffer<double> device_element_jacobian(cells*34*34);
	iga::cuda::DeviceBuffer<double> device_residual(dofs),device_jacobian(
		static_cast<std::size_t>(dofs)*dofs),device_rhs(dofs);
	iga::cuda::DeviceBuffer<double> device_prescribed(dofs),device_boundary_load(dofs);
	device_quadrature.CopyFromHost(quadrature.data(),quadrature.size());
	device_rows.CopyFromHost(rows.data(),rows.size());
	device_constrained.CopyFromHost(constrained.data(),constrained.size());
	device_state.CopyFromHost(state.data(),state.size());
	device_mesh_velocity.CopyFromHost(mesh_velocity.data(),mesh_velocity.size());
	device_acceleration.CopyFromHost(acceleration.data(),acceleration.size());
	device_prescribed.CopyFromHost(prescribed.data(),prescribed.size());
	device_boundary_load.CopyFromHost(boundary_load.data(),boundary_load.size());
	DenseSolver solver(dofs,device_jacobian.data());
	std::vector<double> increment(dofs);
	std::vector<std::pair<double,std::filesystem::path>> snapshots;
	std::filesystem::create_directories(output);
	double assembly_seconds=0.,solve_seconds=0.;
	auto boundary_state=iga::EvaluateNativeTetBoundaryFlows(mesh,topology,state);
	for(int step=1;step<=specification.steps;++step){
		std::optional<iga::ZeroDFlowTrial> source_trial;
		std::vector<iga::ZeroDFlowTrial> terminal_trials;
		if(specification.t7_enabled){
			source_trial=iga::EvaluateZeroDFlowTrial(specification.t7_source,
				specification.t7_source_state,
				boundary_state.at(specification.t7_inlet_label).average_pressure_pa,
				specification.dt_s,(step-1)*specification.dt_s);
			specification.velocity_by_label[specification.t7_inlet_label]=
				UniformVelocityForFlow(area_vectors.at(specification.t7_inlet_label),
					-*source_trial->port.outward_flow_m3_s);
			for(const auto& outlet:specification.t7_outlets){
				terminal_trials.push_back(iga::EvaluateZeroDFlowTrial(outlet.model,
					outlet.state,-boundary_state.at(outlet.boundary_label).outward_flow_m3_s,
					specification.dt_s,(step-1)*specification.dt_s));
				specification.pressure_by_label[outlet.boundary_label]=
					*terminal_trials.back().port.mean_pressure_pa;
			}
			for(const auto& boundary:specification.velocity_by_label)
				for(const auto node:topology.boundary_velocity_nodes.at(boundary.first))
					for(int axis=0;axis<3;++axis){
						const int row=3*static_cast<int>(node)+axis;
						prescribed[row]=boundary.second[axis];
					}
			boundary_load=BoundaryLoad(mesh,topology,
				specification.pressure_by_label,dofs);
			device_prescribed.CopyFromHost(prescribed.data(),prescribed.size());
			device_boundary_load.CopyFromHost(boundary_load.data(),boundary_load.size());
		}
		iga::cuda::Check(cudaMemcpy(device_previous.data(),device_state.data(),
			static_cast<std::size_t>(dofs)*sizeof(double),cudaMemcpyDeviceToDevice),
			"copy previous FEM CUDA state");
		bool converged=false;
		for(int iteration=1;iteration<=specification.maximum_iterations;++iteration){
			const auto assembly_start=std::chrono::steady_clock::now();
			device_residual.Clear();device_jacobian.Clear();
			const int field_count=static_cast<int>(quadrature.size());
			const int entry_count=cells*iga::cuda::kTetFlowEntries;
			iga::cuda::EvaluateNativeTetFlowFields<<<(field_count+127)/128,128>>>(
				device_quadrature.data(),device_rows.data(),device_state.data(),
				device_previous.data(),device_mesh_velocity.data(),cells,fields.data());
			iga::cuda::AssembleNativeTetFlowElements<<<(entry_count+127)/128,128>>>(
				device_quadrature.data(),fields.data(),device_acceleration.data(),
				specification.density,specification.viscosity,specification.dt_s,
				cells,device_element_residual.data(),device_element_jacobian.data());
			iga::cuda::ScatterNativeTetFlowElements<<<(entry_count+127)/128,128>>>(
				device_rows.data(),device_element_residual.data(),
				device_element_jacobian.data(),cells,dofs,
				device_residual.data(),device_jacobian.data());
			iga::cuda::ApplyNativeTetFlowBoundaryRows<<<
				(static_cast<std::size_t>(dofs)*dofs+127)/128,128>>>(
				device_constrained.data(),dofs,device_jacobian.data());
			iga::cuda::NativeTetFlowNewtonRightHandSide<<<(dofs+127)/128,128>>>(
				device_residual.data(),device_boundary_load.data(),
				device_state.data(),device_constrained.data(),
				device_prescribed.data(),dofs,device_rhs.data());
			iga::cuda::CheckKernel("FEM CUDA global assembly");
			iga::cuda::Check(cudaDeviceSynchronize(),"FEM CUDA assembly sync");
			assembly_seconds+=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-assembly_start).count();
			const auto solve_start=std::chrono::steady_clock::now();
			solver.Solve(dofs,device_jacobian.data(),device_rhs.data());
			device_rhs.CopyToHost(increment.data(),increment.size());
			solve_seconds+=std::chrono::duration<double>(
				std::chrono::steady_clock::now()-solve_start).count();
			iga::cuda::UpdateNativeTetFlowState<<<(dofs+127)/128,128>>>(
				device_state.data(),device_rhs.data(),dofs);
			iga::cuda::CheckKernel("FEM CUDA state update");
			const double maximum_increment=*std::max_element(increment.begin(),
				increment.end(),[](double a,double b){return std::abs(a)<std::abs(b);});
			std::cout<<std::setprecision(17)<<"native_tet_flow_cuda_iteration step="
				<<step<<" iteration="<<iteration<<" max_increment="
				<<std::abs(maximum_increment)<<'\n';
			if(std::abs(maximum_increment)<specification.tolerance){
				converged=true;break;
			}
		}
		if(!converged)throw std::runtime_error("FEM CUDA Newton iteration did not converge");
		device_state.CopyToHost(state.data(),state.size());
		boundary_state=iga::EvaluateNativeTetBoundaryFlows(mesh,topology,state);
		if(specification.t7_enabled){
			const auto source_initial=specification.t7_source_state;
			specification.t7_source_state=source_trial->state;
			const auto source_accounting=iga::MakeZeroDFlowStepAccounting(
				specification.t7_source,source_initial,*source_trial);
			std::cout<<std::setprecision(17)<<"native_tet_flow_cuda_t7 step="<<step
				<<" port=inlet boundary_label="<<specification.t7_inlet_label
				<<" pressure_pa="<<boundary_state.at(specification.t7_inlet_label).average_pressure_pa
				<<" flow_m3_s="<<boundary_state.at(specification.t7_inlet_label).outward_flow_m3_s
				<<" peer_flow_m3_s="<<*source_trial->port.outward_flow_m3_s
				<<" storage_residual_m3="<<source_accounting.residual_m3<<'\n';
			for(std::size_t index=0;index<specification.t7_outlets.size();++index){
				auto& outlet=specification.t7_outlets[index];
				const auto initial=outlet.state;
				const auto accepted=iga::EvaluateZeroDFlowTrial(outlet.model,initial,
					-boundary_state.at(outlet.boundary_label).outward_flow_m3_s,
					specification.dt_s,(step-1)*specification.dt_s);
				const auto accounting=iga::MakeZeroDFlowStepAccounting(
					outlet.model,initial,accepted);
				outlet.state=accepted.state;
				std::cout<<std::setprecision(17)<<"native_tet_flow_cuda_t7 step="<<step
					<<" port=outlet boundary_label="<<outlet.boundary_label
					<<" pressure_pa="<<boundary_state.at(outlet.boundary_label).average_pressure_pa
					<<" flow_m3_s="<<boundary_state.at(outlet.boundary_label).outward_flow_m3_s
					<<" next_rcr_pressure_pa="<<*accepted.port.mean_pressure_pa
					<<" storage_residual_m3="<<accounting.residual_m3<<'\n';
			}
		}
		const auto velocity_file=output/("flow_velocity_step_"
			+std::to_string(step)+".bin");
		std::ofstream velocity_stream(velocity_file,std::ios::binary);
		velocity_stream.write(reinterpret_cast<const char*>(state.data()),
			static_cast<std::streamsize>(3*velocity_nodes*sizeof(double)));
		if(!velocity_stream)throw std::runtime_error("cannot write FEM CUDA flow velocity");
		const auto snapshot=output/("flow_step_"+std::to_string(step)+".vtu");
		const auto piece=iga::BuildNativeTetHydraulicVtkPartition(mesh,mesh,state,
			specification.viscosity,0,1);
		iga::WriteVtuPartition(snapshot,piece,step*specification.dt_s);
		snapshots.push_back({step*specification.dt_s,snapshot});
		iga::WritePvd(output/"flow.pvd",snapshots);
		std::cout<<"native_tet_flow_cuda_step step="<<step
			<<" time_s="<<step*specification.dt_s<<'\n';
	}
	const auto usage=[](){struct rusage value{};getrusage(RUSAGE_SELF,&value);
		return value.ru_maxrss;}();
	std::cout<<"native_tet_flow_cuda_complete assembly_s="<<assembly_seconds
		<<" solve_s="<<solve_seconds<<" host_peak_rss_kib="<<usage
		<<" cuda_peak_bytes="<<iga::cuda::DeviceAllocationCounter::Peak()<<'\n';
}

} // namespace

int main(int argc,char** argv)
{
	try{
		if(argc==3&&std::string(argv[1])=="--verify-element")
			VerifyElement(argv[2]);
		else if(argc==3)RunCase(argv[1],argv[2]);
		else throw std::runtime_error(
			"usage: native_tet_flow_cuda CASE.json OUTPUT_DIR | --verify-element MESH.msh");
		return 0;
	}catch(const std::exception& error){
		std::cerr<<"native_tet_flow_cuda: "<<error.what()<<'\n';
		return 1;
	}
}

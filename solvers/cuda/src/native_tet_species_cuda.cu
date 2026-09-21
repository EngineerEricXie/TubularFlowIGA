#include "CudaRuntime.hpp"
#include "NativeTetSpeciesKernels.cuh"
#include "NativeTetSpeciesVisualization.hpp"
#include "VtkOutput.hpp"
#include "CaseConfig.hpp"

#include <cusolverDn.h>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/resource.h>

namespace {

struct Case
{
	std::filesystem::path mesh_file;
	std::string species_id;
	std::array<double,3> fluid_velocity{},mesh_velocity{};
	std::map<int,double> inflow;
	struct Exchange {double transfer=0.,external=0.;};
	std::map<int,Exchange> wall_exchange;
	double initial=0.,diffusivity=0.,source=0.,decay=0.,dt=0.;
	int steps=0;
};

Case ReadCase(const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error("cannot open FEM CUDA species case");
	const auto parsed=JsonParser(iga::ReadCheckedText(stream)).Parse();
	const auto& root=RequireObject(parsed,"FEM CUDA species case");
	auto field=[](const auto& object,const char* key)->const JsonValue& {
		const auto* value=Find(object,key);
		if(!value)throw std::invalid_argument(std::string("species case requires ")+key);
		return *value;
	};
	auto number=[&](const auto& object,const char* key){
		return RequireNumber(field(object,key),key);};
	if(RequireInteger(field(root,"schema_version"),"schema_version")!=1)
		throw std::invalid_argument("unsupported FEM CUDA species schema");
	Case result;
	result.mesh_file=path.parent_path()/RequireString(field(root,"mesh_file"),"mesh_file");
	result.species_id=RequireString(field(root,"species_id"),"species_id");
	result.initial=number(root,"initial_concentration_mol_m3");
	result.diffusivity=number(root,"diffusivity_m2_s");
	result.source=number(root,"source_mol_m3_s");
	if(const auto* value=Find(root,"first_order_decay_rate_s_inv"))
		result.decay=RequireNumber(*value,"first_order_decay_rate_s_inv");
	for(const auto& specification:std::array<std::pair<const char*,std::array<double,3>*>,2>{{
		{"fluid_velocity_m_s",&result.fluid_velocity},
		{"mesh_velocity_m_s",&result.mesh_velocity}}}){
		const auto& values=RequireArray(field(root,specification.first),specification.first);
		if(values.size()!=3)throw std::invalid_argument("species velocity requires three components");
		for(int axis=0;axis<3;++axis)
			(*specification.second)[axis]=RequireNumber(values[axis],specification.first);
	}
	const auto& inflow=RequireObject(field(root,"inflow_concentration_by_label_mol_m3"),
		"inflow_concentration_by_label_mol_m3");
	for(const auto& item:inflow)
		result.inflow.emplace(std::stoi(item.first),RequireNumber(item.second,
			"inflow concentration"));
	if(const auto* value=Find(root,"wall_exchange_by_label"))
		for(const auto& item:RequireObject(*value,"wall_exchange_by_label")){
			const auto& exchange=RequireObject(item.second,"wall_exchange");
			result.wall_exchange.emplace(std::stoi(item.first),Case::Exchange{
				number(exchange,"transfer_coefficient_m_s"),
				number(exchange,"external_concentration_mol_m3")});
		}
	if(const auto* value=Find(root,"monotone"))
		if(RequireBoolean(*value,"monotone"))
			throw std::invalid_argument("FEM CUDA species monotone mode is not implemented");
	if(Find(root,"finite_wall_reservoir")||Find(root,"finite_wall_reservoirs_by_label"))
		throw std::invalid_argument("FEM CUDA finite wall reservoir is not implemented");
	for(double component:result.mesh_velocity)
		if(component!=0.)throw std::invalid_argument("FEM CUDA species requires a fixed mesh");
	const auto& time=RequireObject(field(root,"time"),"time");
	result.dt=number(time,"dt_s");
	result.steps=RequireInteger(field(time,"steps"),"steps");
	return result;
}

void CheckSolver(cusolverStatus_t status,const char* operation)
{
	if(status!=CUSOLVER_STATUS_SUCCESS)
		throw std::runtime_error(std::string(operation)+" failed");
}

void RunCase(const std::filesystem::path& case_path,const std::filesystem::path& output,
	const std::filesystem::path& velocity_series,
	const std::filesystem::path& darcy_flux_directory)
{
	const auto specification=ReadCase(case_path);
	std::ifstream input(specification.mesh_file);
	if(!input)throw std::runtime_error("cannot open FEM CUDA species mesh");
	const auto mesh=iga::ReadNativeTetMeshGmsh41(input);
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const int dofs=static_cast<int>(mesh.points.size());
	const int cell_count=static_cast<int>(mesh.cells.size());
	const int velocity_nodes=static_cast<int>(mesh.points.size()+topology.edges.size());
	std::array<double,40> moments{};
	for(const auto& point:iga::NativeTetDegreeFiveQuadrature()){
		const auto basis=iga::EvaluateNativeTaylorHoodBasis(point.reference[0],
			point.reference[1],point.reference[2]);
		for(int column=0;column<4;++column)
			for(int local=0;local<10;++local)
				moments[10*column+local]+=6.*point.weight
					*basis.pressure[column]*basis.velocity[local];
	}
	std::vector<iga::cuda::NativeTetSpeciesCell> cells(cell_count);
	for(int index=0;index<cell_count;++index){
		const auto& cell=mesh.cells[index];
		const auto geometry=iga::EvaluateNativeTetGeometry(mesh,cell);
		auto& device=cells[index];
		device.volume=geometry.determinant/6.;
		for(int node=0;node<4;++node){
			device.nodes[node]=static_cast<int>(cell.nodes[node]);
			for(int axis=0;axis<3;++axis)
				device.gradient[node][axis]=geometry.barycentric_gradients[node][axis];
			for(int axis=0;axis<3;++axis)
				device.position[node][axis]=mesh.points[cell.nodes[node]][axis];
		}
		for(int node=0;node<10;++node)
			device.velocity_nodes[node]=static_cast<int>(
				topology.cell_velocity_nodes[index][node]);
	}
	std::vector<iga::cuda::NativeTetSpeciesFace> faces;
	for(const auto& triangle:mesh.boundary_triangles){
		const auto& x0=mesh.points[triangle.nodes[0]];
		const auto& x1=mesh.points[triangle.nodes[1]];
		const auto& x2=mesh.points[triangle.nodes[2]];
		const std::array<double,3> a{{x1[0]-x0[0],x1[1]-x0[1],x1[2]-x0[2]}};
		const std::array<double,3> b{{x2[0]-x0[0],x2[1]-x0[1],x2[2]-x0[2]}};
		std::array<double,3> area_vector{{
			0.5*(a[1]*b[2]-a[2]*b[1]),
			0.5*(a[2]*b[0]-a[0]*b[2]),
			0.5*(a[0]*b[1]-a[1]*b[0])}};
		int owner_index=0,owner_opposite=0;
		for(const auto& cell:mesh.cells){
			int opposite=-1;
			for(int local=0;local<4;++local)
				if(cell.nodes[local]!=triangle.nodes[0]
					&&cell.nodes[local]!=triangle.nodes[1]
					&&cell.nodes[local]!=triangle.nodes[2])opposite=local;
			if(opposite<0)continue;
			int shared=0;
			for(const auto node:cell.nodes)
				for(const auto face_node:triangle.nodes)
					if(node==face_node)++shared;
			if(shared!=3){++owner_index;continue;}
			const auto& interior=mesh.points[cell.nodes[opposite]];
			double inward=0.;
			for(int axis=0;axis<3;++axis)
				inward+=area_vector[axis]*(interior[axis]-x0[axis]);
			if(inward>0.)for(auto& component:area_vector)component=-component;
			owner_opposite=opposite;
			break;
		}
		iga::cuda::NativeTetSpeciesFace face{};
		const auto& owner=mesh.cells[owner_index];
		face.cell_index=owner_index;
		face.opposite=owner_opposite;
		for(int node=0;node<3;++node){
			face.nodes[node]=static_cast<int>(triangle.nodes[node]);
			for(int local=0;local<4;++local)
				if(owner.nodes[local]==triangle.nodes[node])face.local_nodes[node]=local;
		}
		for(int node=0;node<10;++node)
			face.velocity_nodes[node]=static_cast<int>(
				topology.cell_velocity_nodes[owner_index][node]);
		face.area=std::hypot(area_vector[0],area_vector[1],area_vector[2]);
		for(int axis=0;axis<3;++axis)
			face.area_vector[axis]=area_vector[axis];
		double normal_flow=0.;
		for(int axis=0;axis<3;++axis)
			normal_flow+=specification.fluid_velocity[axis]*area_vector[axis];
		const auto inflow=specification.inflow.find(triangle.boundary_label);
		if(inflow!=specification.inflow.end())
			face.inlet_concentration=inflow->second;
		if(normal_flow<0.){
			if(inflow==specification.inflow.end())
				throw std::invalid_argument("species inflow boundary has no concentration");
		}
		const auto wall=specification.wall_exchange.find(triangle.boundary_label);
		if(wall!=specification.wall_exchange.end()){
			if(normal_flow!=0.)
				throw std::invalid_argument("species wall exchange has relative flow");
			face.transfer=wall->second.transfer;
			face.external_concentration=wall->second.external;
		}
		faces.push_back(face);
	}
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetSpeciesCell> device_cells(cells.size());
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetSpeciesFace> device_faces(faces.size());
	iga::cuda::DeviceBuffer<double> matrix(static_cast<std::size_t>(dofs)*dofs);
	iga::cuda::DeviceBuffer<double> rhs(dofs),previous(dofs);
	iga::cuda::DeviceBuffer<double> device_velocity(static_cast<std::size_t>(velocity_nodes)*3);
	iga::cuda::DeviceBuffer<double> device_moments(moments.size());
	iga::cuda::DeviceBuffer<double> device_rt0(static_cast<std::size_t>(cell_count)*4);
	iga::cuda::DeviceBuffer<int> pivots(dofs),info(1);
	device_cells.CopyFromHost(cells.data(),cells.size());
	device_faces.CopyFromHost(faces.data(),faces.size());
	device_moments.CopyFromHost(moments.data(),moments.size());
	if(!darcy_flux_directory.empty()){
		std::vector<double> face_flow(static_cast<std::size_t>(cell_count)*4);
		std::ifstream stream(darcy_flux_directory/"darcy_rt0_face_flow.bin",
			std::ios::binary);
		if(!stream.read(reinterpret_cast<char*>(face_flow.data()),
			static_cast<std::streamsize>(face_flow.size()*sizeof(double))))
			throw std::runtime_error("cannot read FEM CUDA Darcy RT0 face flow");
		device_rt0.CopyFromHost(face_flow.data(),face_flow.size());
	}
	cusolverDnHandle_t handle=nullptr;
	CheckSolver(cusolverDnCreate(&handle),"species cusolverDnCreate");
	int workspace=0;
	CheckSolver(cusolverDnDgetrf_bufferSize(handle,dofs,dofs,matrix.data(),dofs,
		&workspace),"species getrf buffer size");
	iga::cuda::DeviceBuffer<double> work(workspace);
	std::filesystem::create_directories(output);
	std::vector<std::pair<double,std::filesystem::path>> snapshots;
	std::vector<double> concentration(dofs,specification.initial);
	std::vector<double> velocity(static_cast<std::size_t>(velocity_nodes)*3);
	double assembly_seconds=0.,solve_seconds=0.;
	for(int step=1;step<=specification.steps;++step){
		if(velocity_series.empty()){
			for(int node=0;node<velocity_nodes;++node)
				for(int axis=0;axis<3;++axis)
					velocity[3*node+axis]=specification.fluid_velocity[axis];
		}else{
			std::ifstream stream(velocity_series/(
				"flow_velocity_step_"+std::to_string(step)+".bin"),std::ios::binary);
			if(!stream.read(reinterpret_cast<char*>(velocity.data()),
				static_cast<std::streamsize>(velocity.size()*sizeof(double))))
				throw std::runtime_error("cannot read FEM CUDA flow velocity step");
		}
		const auto assembly_start=std::chrono::steady_clock::now();
		matrix.Clear();rhs.Clear();
		previous.CopyFromHost(concentration.data(),concentration.size());
		device_velocity.CopyFromHost(velocity.data(),velocity.size());
		iga::cuda::AssembleNativeTetSpeciesCells<<<(cell_count*16+127)/128,128>>>(
			device_cells.data(),cell_count,dofs,specification.dt,
			specification.diffusivity,specification.source,specification.decay,
			device_velocity.data(),device_moments.data(),
			darcy_flux_directory.empty()?nullptr:device_rt0.data(),previous.data(),
			matrix.data(),rhs.data());
		iga::cuda::AssembleNativeTetSpeciesFaces<<<(faces.size()*9+127)/128,128>>>(
			device_faces.data(),static_cast<int>(faces.size()),dofs,
			device_velocity.data(),
			darcy_flux_directory.empty()?nullptr:device_rt0.data(),
			matrix.data(),rhs.data());
		iga::cuda::CheckKernel("FEM CUDA species assembly");
		iga::cuda::Check(cudaDeviceSynchronize(),"FEM CUDA species assembly sync");
		assembly_seconds+=std::chrono::duration<double>(
			std::chrono::steady_clock::now()-assembly_start).count();
		const auto solve_start=std::chrono::steady_clock::now();
		CheckSolver(cusolverDnDgetrf(handle,dofs,dofs,matrix.data(),dofs,work.data(),
			pivots.data(),info.data()),"species getrf");
		CheckSolver(cusolverDnDgetrs(handle,CUBLAS_OP_N,dofs,1,matrix.data(),dofs,
			pivots.data(),rhs.data(),dofs,info.data()),"species getrs");
		int result_info=0;info.CopyToHost(&result_info,1);
		if(result_info!=0)throw std::runtime_error("FEM CUDA species matrix is singular");
		solve_seconds+=std::chrono::duration<double>(
			std::chrono::steady_clock::now()-solve_start).count();
		rhs.CopyToHost(concentration.data(),concentration.size());
		const auto snapshot=output/("species_step_"+std::to_string(step)+".vtu");
		const auto piece=iga::BuildNativeTetSpeciesVtkPartition(mesh,mesh,
			concentration,0,1);
		iga::WriteVtuPartition(snapshot,piece,step*specification.dt);
		snapshots.push_back({step*specification.dt,snapshot});
		iga::WritePvd(output/"species.pvd",snapshots);
	}
	const auto usage=[](){struct rusage value{};getrusage(RUSAGE_SELF,&value);
		return value.ru_maxrss;}();
	std::cout<<"native_tet_species_cuda_complete species="<<specification.species_id
		<<" cells="<<cell_count<<" steps="<<specification.steps
		<<" assembly_s="<<assembly_seconds<<" solve_s="<<solve_seconds
		<<" host_peak_rss_kib="<<usage
		<<" cuda_peak_bytes="<<iga::cuda::DeviceAllocationCounter::Peak()<<'\n';
	CheckSolver(cusolverDnDestroy(handle),"species cusolverDnDestroy");
}

} // namespace

int main(int argc,char** argv)
{
	try{
		if(argc!=3&&argc!=5)
			throw std::invalid_argument("usage: native_tet_species_cuda CASE.json OUTPUT_DIR [--velocity-series DIR | --darcy-flux DIR]");
		if(argc==5&&std::string(argv[3])!="--velocity-series"
				&&std::string(argv[3])!="--darcy-flux")
			throw std::invalid_argument("species field option is invalid");
		RunCase(argv[1],argv[2],
			argc==5&&std::string(argv[3])=="--velocity-series"?argv[4]:"",
			argc==5&&std::string(argv[3])=="--darcy-flux"?argv[4]:"");
		return 0;
	}catch(const std::exception& error){
		std::cerr<<"native_tet_species_cuda: "<<error.what()<<'\n';
		return 1;
	}
}

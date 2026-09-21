#include "CudaRuntime.hpp"
#include "NativeTetDarcyKernels.cuh"
#include "NativeTetDarcyVisualization.hpp"
#include "NativeTetVesselTissueSourceMap.hpp"
#include "VtkOutput.hpp"
#include "CaseConfig.hpp"

#include <cusolverDn.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/resource.h>

namespace {

struct Case
{
	std::filesystem::path mesh_file;
	double mobility=0.,source=0.;
	std::map<std::uint64_t,double> mobility_by_cell,source_by_cell;
	std::map<int,double> pressure_by_label,flux_by_label;
};

struct FlowSourceMap
{
	std::filesystem::path vessel_mesh_file,flow_velocity_file;
	std::vector<iga::NativeTetVesselTissueBoundaryPort> ports;
};

FlowSourceMap ReadFlowSourceMap(const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error("cannot open FEM CUDA flow-Darcy source map");
	const auto parsed=JsonParser(iga::ReadCheckedText(stream)).Parse();
	const auto& root=RequireObject(parsed,"FEM CUDA flow-Darcy source map");
	auto required=[](const auto& object,const char* key)->const JsonValue&
	{
		const auto* value=Find(object,key);
		if(!value)throw std::invalid_argument(std::string("flow-Darcy source map requires ")+key);
		return *value;
	};
	if(RequireInteger(required(root,"schema_version"),"schema_version")!=1)
		throw std::invalid_argument("unsupported FEM CUDA flow-Darcy source map schema");
	FlowSourceMap result;
	result.vessel_mesh_file=path.parent_path()/RequireString(
		required(root,"vessel_mesh_file"),"vessel_mesh_file");
	result.flow_velocity_file=path.parent_path()/RequireString(
		required(root,"flow_velocity_file"),"flow_velocity_file");
	for(const auto& item:RequireArray(required(root,"ports"),"ports"))
	{
		const auto& object=RequireObject(item,"flow-Darcy port");
		iga::NativeTetVesselTissueBoundaryPort port;
		port.name=RequireString(required(object,"name"),"port name");
		port.vessel_boundary_label=RequireInteger(
			required(object,"vessel_boundary_label"),"vessel_boundary_label");
		for(const auto& target:RequireObject(required(object,"tissue_cell_id_weights"),
			"tissue_cell_id_weights"))
			port.tissue_cell_id_weights.push_back({std::stoull(target.first),
				RequireNumber(target.second,"tissue cell weight")});
		result.ports.push_back(std::move(port));
	}
	return result;
}

iga::NativeTetVesselTissueSourceResult MapFlowSource(
	const FlowSourceMap& specification,const iga::NativeTetMesh& tissue_mesh)
{
	std::ifstream mesh_stream(specification.vessel_mesh_file);
	if(!mesh_stream)throw std::runtime_error("cannot open FEM CUDA source vessel mesh");
	const auto vessel_mesh=iga::ReadNativeTetMeshGmsh41(mesh_stream);
	const auto topology=iga::BuildNativeTaylorHoodTopology(vessel_mesh);
	const std::size_t velocity_nodes=vessel_mesh.points.size()+topology.edges.size();
	std::vector<double> state(3*velocity_nodes+vessel_mesh.points.size(),0.);
	std::ifstream velocity_stream(specification.flow_velocity_file,std::ios::binary);
	velocity_stream.read(reinterpret_cast<char*>(state.data()),
		static_cast<std::streamsize>(3*velocity_nodes*sizeof(double)));
	if(!velocity_stream)throw std::runtime_error("cannot read FEM CUDA source velocity");
	return iga::MapNativeTetVesselStateToTissueSource(vessel_mesh,topology,state,
		tissue_mesh,specification.ports);
}

Case ReadCase(const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	std::ifstream stream(path);
	if(!stream)throw std::runtime_error("cannot open FEM CUDA Darcy case");
	const auto parsed=JsonParser(iga::ReadCheckedText(stream)).Parse();
	const auto& root=RequireObject(parsed,"FEM CUDA Darcy case");
	auto required=[](const auto& object,const char* key)->const JsonValue& {
		const auto* value=Find(object,key);
		if(!value)throw std::invalid_argument(std::string("Darcy case requires ")+key);
		return *value;
	};
	if(RequireInteger(required(root,"schema_version"),"schema_version")!=1)
		throw std::invalid_argument("unsupported FEM CUDA Darcy schema");
	Case result;
	result.mesh_file=path.parent_path()/RequireString(required(root,"mesh_file"),"mesh_file");
	result.mobility=RequireNumber(required(root,"mobility_m2_pa_s"),"mobility_m2_pa_s");
	result.source=RequireNumber(required(root,"source_s_inv"),"source_s_inv");
	const auto& pressures=RequireObject(required(root,"pressure_by_boundary_label_pa"),
		"pressure_by_boundary_label_pa");
	for(const auto& item:pressures)
		result.pressure_by_label.emplace(std::stoi(item.first),
			RequireNumber(item.second,"Darcy boundary pressure"));
	if(const auto* flux=Find(root,"outward_flux_by_boundary_label_m_s"))
		for(const auto& item:RequireObject(*flux,"outward_flux_by_boundary_label_m_s"))
			result.flux_by_label.emplace(std::stoi(item.first),
				RequireNumber(item.second,"Darcy boundary flux"));
	if(const auto* mobility=Find(root,"mobility_by_cell_id_m2_pa_s"))
		for(const auto& item:RequireObject(*mobility,"mobility_by_cell_id_m2_pa_s"))
			result.mobility_by_cell.emplace(std::stoull(item.first),
				RequireNumber(item.second,"Darcy cell mobility"));
	if(const auto* source=Find(root,"source_by_cell_id_s_inv"))
		for(const auto& item:RequireObject(*source,"source_by_cell_id_s_inv"))
			result.source_by_cell.emplace(std::stoull(item.first),
				RequireNumber(item.second,"Darcy cell source"));
	if(!(result.mobility>0.)||result.pressure_by_label.empty())
		throw std::invalid_argument("invalid FEM CUDA Darcy material or pressure boundary");
	return result;
}

void CheckSolver(cusolverStatus_t status,const char* operation)
{
	if(status!=CUSOLVER_STATUS_SUCCESS)
		throw std::runtime_error(std::string(operation)+" failed with status "
			+std::to_string(static_cast<int>(status)));
}

void RunCase(const std::filesystem::path& case_path,
	const std::filesystem::path& output,
	const std::filesystem::path& source_map_path)
{
	const auto specification=ReadCase(case_path);
	std::ifstream input(specification.mesh_file);
	if(!input)throw std::runtime_error("cannot open FEM CUDA Darcy mesh");
	const auto mesh=iga::ReadNativeTetMeshGmsh41(input);
	const int dofs=static_cast<int>(mesh.points.size());
	const int cell_count=static_cast<int>(mesh.cells.size());
	std::vector<iga::cuda::NativeTetDarcyCell> cells(cell_count);
	std::vector<double> mobility(cell_count),source(cell_count);
	for(int index=0;index<cell_count;++index){
		const auto& cell=mesh.cells[index];
		const auto geometry=iga::EvaluateNativeTetGeometry(mesh,cell);
		const auto material=specification.mobility_by_cell.find(cell.id);
		const auto forcing=specification.source_by_cell.find(cell.id);
		mobility[index]=material==specification.mobility_by_cell.end()
			?specification.mobility:material->second;
		source[index]=forcing==specification.source_by_cell.end()
			?specification.source:forcing->second;
		auto& device=cells[index];
		device.volume=geometry.determinant/6.;
		device.mobility=mobility[index];device.source=source[index];
		for(int node=0;node<4;++node){
			device.nodes[node]=static_cast<int>(cell.nodes[node]);
			for(int axis=0;axis<3;++axis)
				device.gradient[node][axis]=geometry.barycentric_gradients[node][axis];
		}
	}
	std::optional<iga::NativeTetVesselTissueSourceResult> mapped_source;
	if(!source_map_path.empty())
	{
		mapped_source=MapFlowSource(ReadFlowSourceMap(source_map_path),mesh);
		for(int index=0;index<cell_count;++index)
		{
			source[index]+=mapped_source->tissue_source_s_inv[index];
			cells[index].source=source[index];
		}
	}
	std::vector<int> fixed(dofs,0);
	std::vector<double> prescribed(dofs,0.),boundary_source(dofs,0.);
	for(const auto& face:mesh.boundary_triangles){
		const auto pressure=specification.pressure_by_label.find(face.boundary_label);
		if(pressure!=specification.pressure_by_label.end())
			for(const auto node:face.nodes){
				fixed[node]=1;prescribed[node]=pressure->second;
			}
		const auto flux=specification.flux_by_label.find(face.boundary_label);
		if(flux!=specification.flux_by_label.end()){
			const auto& a=mesh.points[face.nodes[0]];
			const auto& b=mesh.points[face.nodes[1]];
			const auto& c=mesh.points[face.nodes[2]];
			const std::array<double,3> first{{b[0]-a[0],b[1]-a[1],b[2]-a[2]}};
			const std::array<double,3> second{{c[0]-a[0],c[1]-a[1],c[2]-a[2]}};
			const std::array<double,3> cross{{
				first[1]*second[2]-first[2]*second[1],
				first[2]*second[0]-first[0]*second[2],
				first[0]*second[1]-first[1]*second[0]}};
			const double area=0.5*std::sqrt(cross[0]*cross[0]
				+cross[1]*cross[1]+cross[2]*cross[2]);
			for(const auto node:face.nodes)
				boundary_source[node]-=flux->second*area/3.;
		}
	}
	iga::cuda::DeviceBuffer<iga::cuda::NativeTetDarcyCell> device_cells(cells.size());
	iga::cuda::DeviceBuffer<double> matrix(static_cast<std::size_t>(dofs)*dofs);
	iga::cuda::DeviceBuffer<double> right_hand_side(dofs),device_pressure(dofs);
	iga::cuda::DeviceBuffer<double> device_flux(static_cast<std::size_t>(cell_count)*3);
	iga::cuda::DeviceBuffer<int> device_fixed(dofs),pivots(dofs),info(1);
	device_cells.CopyFromHost(cells.data(),cells.size());
	device_fixed.CopyFromHost(fixed.data(),fixed.size());
	device_pressure.CopyFromHost(prescribed.data(),prescribed.size());
	matrix.Clear();
	right_hand_side.CopyFromHost(boundary_source.data(),boundary_source.size());
	const auto assembly_start=std::chrono::steady_clock::now();
	iga::cuda::AssembleNativeTetDarcy<<<(cell_count*20+127)/128,128>>>(
		device_cells.data(),cell_count,dofs,matrix.data(),right_hand_side.data());
	iga::cuda::ApplyNativeTetDarcyBoundaryRows<<<
		(static_cast<std::size_t>(dofs)*dofs+127)/128,128>>>(
		device_fixed.data(),dofs,matrix.data());
	iga::cuda::SetNativeTetDarcyBoundaryValues<<<(dofs+127)/128,128>>>(
		device_fixed.data(),device_pressure.data(),dofs,right_hand_side.data());
	iga::cuda::CheckKernel("FEM CUDA Darcy assembly");
	iga::cuda::Check(cudaDeviceSynchronize(),"FEM CUDA Darcy assembly sync");
	const double assembly_seconds=std::chrono::duration<double>(
		std::chrono::steady_clock::now()-assembly_start).count();
	cusolverDnHandle_t handle=nullptr;
	CheckSolver(cusolverDnCreate(&handle),"Darcy cusolverDnCreate");
	int workspace=0;
	CheckSolver(cusolverDnDgetrf_bufferSize(handle,dofs,dofs,matrix.data(),dofs,
		&workspace),"Darcy getrf buffer size");
	iga::cuda::DeviceBuffer<double> work(workspace);
	const auto solve_start=std::chrono::steady_clock::now();
	CheckSolver(cusolverDnDgetrf(handle,dofs,dofs,matrix.data(),dofs,work.data(),
		pivots.data(),info.data()),"Darcy getrf");
	CheckSolver(cusolverDnDgetrs(handle,CUBLAS_OP_N,dofs,1,matrix.data(),dofs,
		pivots.data(),right_hand_side.data(),dofs,info.data()),"Darcy getrs");
	int result_info=0;info.CopyToHost(&result_info,1);
	if(result_info!=0)throw std::runtime_error("FEM CUDA Darcy matrix is singular");
	const double solve_seconds=std::chrono::duration<double>(
		std::chrono::steady_clock::now()-solve_start).count();
	iga::cuda::EvaluateNativeTetDarcyCellFlux<<<(cell_count+127)/128,128>>>(
		device_cells.data(),cell_count,right_hand_side.data(),device_flux.data());
	iga::cuda::CheckKernel("FEM CUDA Darcy flux");
	iga::NativeTetDarcyResult result;
	result.pressure_pa.resize(dofs);
	right_hand_side.CopyToHost(result.pressure_pa.data(),result.pressure_pa.size());
	std::vector<double> flat_flux(static_cast<std::size_t>(cell_count)*3);
	device_flux.CopyToHost(flat_flux.data(),flat_flux.size());
	result.cell_flux_m_s.resize(cell_count);
	for(int cell=0;cell<cell_count;++cell)
		for(int axis=0;axis<3;++axis)
			result.cell_flux_m_s[cell][axis]=flat_flux[3*cell+axis];
	iga::CompleteNativeTetDarcyFlux(mesh,mobility,source,
		specification.pressure_by_label,specification.flux_by_label,result);
	std::filesystem::create_directories(output);
	std::ofstream flow_stream(output/"darcy_rt0_face_flow.bin",std::ios::binary);
	flow_stream.write(reinterpret_cast<const char*>(
		result.conservative_face_flow_m3_s.data()),
		static_cast<std::streamsize>(result.conservative_face_flow_m3_s.size()
			*4*sizeof(double)));
	if(!flow_stream)throw std::runtime_error("cannot write FEM CUDA Darcy RT0 flow");
	const auto piece=iga::BuildNativeTetDarcyVtkPartition(mesh,result,
		mobility,source,0,1);
	iga::WriteVtuPartition(output/"tissue.vtu",piece,0.);
	iga::WritePvd(output/"tissue.pvd",{{0.,output/"tissue.vtu"}});
	const auto usage=[](){struct rusage value{};getrusage(RUSAGE_SELF,&value);
		return value.ru_maxrss;}();
	std::cout<<"native_tet_darcy_cuda_complete cells="<<cell_count
		<<" max_cell_balance_defect_m3_s="<<result.maximum_cell_balance_defect_m3_s
		<<" flow_source_m3_s="<<(mapped_source
			?mapped_source->total_vessel_outward_flow_m3_s:0.)
		<<" source_map_defect_m3_s="<<(mapped_source
			?mapped_source->maximum_port_balance_defect_m3_s:0.)
		<<" darcy_volume_source_m3_s="<<result.volume_source_m3_s
		<<" assembly_s="<<assembly_seconds<<" solve_s="<<solve_seconds
		<<" host_peak_rss_kib="<<usage
		<<" cuda_peak_bytes="<<iga::cuda::DeviceAllocationCounter::Peak()<<'\n';
	CheckSolver(cusolverDnDestroy(handle),"Darcy cusolverDnDestroy");
}

} // namespace

int main(int argc,char** argv)
{
	try{
		if(argc!=3&&(argc!=5||std::string(argv[3])!="--source-map"))
			throw std::invalid_argument(
				"usage: native_tet_darcy_cuda CASE.json OUTPUT_DIR [--source-map MAP.json]");
		RunCase(argv[1],argv[2],argc==5?argv[4]:"");
		return 0;
	}catch(const std::exception& error){
		std::cerr<<"native_tet_darcy_cuda: "<<error.what()<<'\n';
		return 1;
	}
}

#include "CaseConfig.hpp"
#include "CollectiveFailure.hpp"
#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetDarcyVisualization.hpp"
#include "ParallelVtkOutput.hpp"
#include "Sha256.hpp"

#include <petscksp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Object=std::map<std::string,iga::config_detail::JsonValue>;

const iga::config_detail::JsonValue& Field(const Object& object,
	const std::string& key,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	if(!value)throw std::invalid_argument(context+" requires "+key);
	return *value;
}

int BoundaryLabel(const std::string& text)
{
	std::size_t used=0;
	int value=-1;
	try{value=std::stoi(text,&used);}
	catch(const std::exception&){throw std::invalid_argument("native Darcy boundary label is invalid");}
	if(value<0||used!=text.size()||std::to_string(value)!=text)
		throw std::invalid_argument("native Darcy boundary label is invalid");
	return value;
}

std::uint64_t CellId(const std::string& text)
{
	std::size_t used=0;
	std::uint64_t value=0;
	try{value=std::stoull(text,&used);}
	catch(const std::exception&){throw std::invalid_argument("native Darcy cell id is invalid");}
	if(value==0||used!=text.size()||std::to_string(value)!=text)
		throw std::invalid_argument("native Darcy cell id is invalid");
	return value;
}

struct Case
{
	std::filesystem::path mesh_file;
	double mobility=0.,source=0.;
	std::map<std::uint64_t,double> mobility_by_cell,source_by_cell;
	std::map<int,double> pressure_by_label,flux_by_label;
};

Case ParseCase(const std::string& text,const std::filesystem::path& case_path)
{
	using namespace iga::config_detail;
	const auto parsed=JsonParser(text).Parse();
	const auto& root=RequireObject(parsed,"native Darcy case");
	RequireKnownKeys(root,{"schema_version","mesh_file","mobility_m2_pa_s",
		"source_s_inv","mobility_by_cell_id_m2_pa_s","source_by_cell_id_s_inv",
		"pressure_by_boundary_label_pa","outward_flux_by_boundary_label_m_s"},
		"native Darcy case");
	if(RequireInteger(Field(root,"schema_version","native Darcy case"),
		"schema_version")!=1)
		throw std::invalid_argument("unsupported native Darcy case schema");
	Case result;
	result.mesh_file=case_path.parent_path()/RequireString(
		Field(root,"mesh_file","native Darcy case"),"mesh_file");
	result.mobility=RequireNumber(Field(root,"mobility_m2_pa_s","native Darcy case"),
		"mobility_m2_pa_s");
	result.source=RequireNumber(Field(root,"source_s_inv","native Darcy case"),
		"source_s_inv");
	if(result.mesh_file.empty()||!(result.mobility>0.))
		throw std::invalid_argument("native Darcy mesh or mobility is invalid");
	const auto& pressure=RequireObject(Field(root,"pressure_by_boundary_label_pa",
		"native Darcy case"),"pressure_by_boundary_label_pa");
	for(const auto& item:pressure)
		result.pressure_by_label.emplace(BoundaryLabel(item.first),
			RequireNumber(item.second,"Darcy boundary pressure"));
	if(result.pressure_by_label.empty())
		throw std::invalid_argument("native Darcy requires a pressure boundary");
	if(const auto* value=Find(root,"outward_flux_by_boundary_label_m_s")){
		const auto& flux=RequireObject(*value,"outward_flux_by_boundary_label_m_s");
		for(const auto& item:flux)
			result.flux_by_label.emplace(BoundaryLabel(item.first),
				RequireNumber(item.second,"Darcy boundary flux"));
	}
	for(const auto& item:result.flux_by_label)
		if(result.pressure_by_label.count(item.first))
			throw std::invalid_argument("native Darcy boundary label has both pressure and flux");
	if(const auto* value=Find(root,"mobility_by_cell_id_m2_pa_s")){
		const auto& cells=RequireObject(*value,"mobility_by_cell_id_m2_pa_s");
		for(const auto& item:cells){
			const double mobility=RequireNumber(item.second,"Darcy cell mobility");
			if(!(mobility>0.))
				throw std::invalid_argument("native Darcy cell mobility is invalid");
			result.mobility_by_cell.emplace(CellId(item.first),mobility);
		}
	}
	if(const auto* value=Find(root,"source_by_cell_id_s_inv")){
		const auto& cells=RequireObject(*value,"source_by_cell_id_s_inv");
		for(const auto& item:cells)
			result.source_by_cell.emplace(CellId(item.first),
				RequireNumber(item.second,"Darcy cell source"));
	}
	return result;
}

std::string Hash(const std::string& value)
{
	iga::Sha256 digest;
	digest.Append(value.data(),value.size());
	return digest.Hex();
}

struct Input
{
	Case specification;
	iga::NativeTetMesh mesh;
	std::vector<double> mobility,source;
	std::string case_hash,mesh_hash;
};

Input ReadInput(const std::filesystem::path& case_path)
{
	std::ifstream case_stream(case_path);
	if(!case_stream)throw std::runtime_error("native Darcy case is missing");
	const std::string case_text=iga::ReadCheckedText(case_stream);
	Input result;
	result.specification=ParseCase(case_text,case_path);
	std::ifstream mesh_stream(result.specification.mesh_file);
	if(!mesh_stream)throw std::runtime_error("native Darcy mesh is missing");
	const std::string mesh_text=iga::ReadCheckedText(mesh_stream);
	std::istringstream input(mesh_text);
	result.mesh=iga::ReadNativeTetMeshGmsh41(input);
	result.case_hash=Hash(case_text);
	result.mesh_hash=Hash(mesh_text);
	std::set<std::uint64_t> ids;
	result.mobility.reserve(result.mesh.cells.size());
	result.source.reserve(result.mesh.cells.size());
	for(const auto& cell:result.mesh.cells){
		ids.insert(cell.id);
		const auto mobility=result.specification.mobility_by_cell.find(cell.id);
		const auto source=result.specification.source_by_cell.find(cell.id);
		result.mobility.push_back(mobility==result.specification.mobility_by_cell.end()
			?result.specification.mobility:mobility->second);
		result.source.push_back(source==result.specification.source_by_cell.end()
			?result.specification.source:source->second);
	}
	for(const auto& item:result.specification.mobility_by_cell)
		if(!ids.count(item.first))
			throw std::invalid_argument("native Darcy mobility cell id is absent");
	for(const auto& item:result.specification.source_by_cell)
		if(!ids.count(item.first))
			throw std::invalid_argument("native Darcy source cell id is absent");
	return result;
}

void WriteSummary(const std::filesystem::path& output,const Input& input,
	const iga::NativeTetDarcyResult& result,int ranks)
{
	const auto temporary=output/"run_summary.json.pending";
	std::ofstream stream(temporary);
	if(!stream)throw std::runtime_error("cannot create native Darcy summary");
	stream<<std::setprecision(17)
		<<"{\n  \"schema_version\": 1,\n"
		<<"  \"kind\": \"native_tet_p1_darcy_steady\",\n"
		<<"  \"case_sha256\": \""<<input.case_hash<<"\",\n"
		<<"  \"mesh_sha256\": \""<<input.mesh_hash<<"\",\n"
		<<"  \"mpi_ranks\": "<<ranks<<",\n"
		<<"  \"vertices\": "<<input.mesh.points.size()<<",\n"
		<<"  \"tetrahedra\": "<<input.mesh.cells.size()<<",\n"
		<<"  \"volume_source_m3_s\": "<<result.volume_source_m3_s<<",\n"
		<<"  \"linear_iterations\": "<<result.linear_iterations<<",\n"
		<<"  \"linear_residual_norm\": "<<result.residual_norm<<",\n"
		<<"  \"boundary_flow_kind\": \"P1_cell_gradient_diagnostic\",\n"
		<<"  \"conservative_face_flow_kind\": \"least_squares_cell_balance_recovery\",\n"
		<<"  \"recovered_vector_field_kind\": \"rt0_hdiv_from_conservative_faces\",\n"
		<<"  \"flux_recovery_iterations\": "<<result.flux_recovery_iterations<<",\n"
		<<"  \"maximum_cell_balance_defect_m3_s\": "
		<<result.maximum_cell_balance_defect_m3_s<<",\n"
		<<"  \"outward_boundary_flow_m3_s\": {";
	bool first=true;
	for(const auto& item:result.outward_boundary_flow_m3_s){
		if(!first)stream<<',';
		first=false;
		stream<<"\n    \""<<item.first<<"\": "<<item.second;
	}
	stream<<"\n  },\n  \"conservative_outward_boundary_flow_m3_s\": {";
	first=true;
	for(const auto& item:result.conservative_outward_boundary_flow_m3_s){
		if(!first)stream<<',';
		first=false;
		stream<<"\n    \""<<item.first<<"\": "<<item.second;
	}
	stream<<"\n  }\n}\n";
	stream.close();
	if(!stream)throw std::runtime_error("cannot write native Darcy summary");
	std::filesystem::rename(temporary,output/"run_summary.json");
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	int status=0;
	try{
		if(argc<3||argc>4)
			throw std::invalid_argument("usage: native_tet_darcy CASE.json (--check-input | --output-dir DIR)");
		const std::filesystem::path case_path=argv[1];
		const std::string option=argv[2];
		const bool check_only=option=="--check-input";
		std::filesystem::path output;
		if(check_only){
			if(argc!=3)throw std::invalid_argument("--check-input cannot select output");
		}else if(option=="--output-dir"&&argc==4){
			output=argv[3];
			if(output.empty())throw std::invalid_argument("native Darcy output path is empty");
		}else throw std::invalid_argument("native Darcy command option is invalid");
		Input input;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native Darcy input",[&]{
			input=ReadInput(case_path);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native Darcy case identity",
			input.case_hash);
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native Darcy mesh identity",
			input.mesh_hash);
		const auto solution=iga::SolveNativeTetDarcyPetsc(input.mesh,input.mobility,
			input.source,input.specification.pressure_by_label,
			input.specification.flux_by_label);
		if(check_only){
			if(rank==0)std::cout<<"native_darcy_input: PASS tetrahedra="
				<<input.mesh.cells.size()<<" vertices="<<input.mesh.points.size()<<'\n';
		}else{
			iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native Darcy output path",
				std::filesystem::absolute(output).lexically_normal().string());
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native Darcy output directory",[&]{
				if(rank!=0)return;
				if(!output.parent_path().empty())
					std::filesystem::create_directories(output.parent_path());
				if(!std::filesystem::create_directory(output))
					throw std::runtime_error("native Darcy output directory already exists");
			});
			iga::VtkPartition piece;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native Darcy VTU build",[&]{
				piece=iga::BuildNativeTetDarcyVtkPartition(input.mesh,solution,
					input.mobility,input.source,rank,ranks);
			});
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"fields",piece,0.);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native Darcy summary",[&]{
				if(rank==0)WriteSummary(output,input,solution,ranks);
			});
			if(rank==0)std::cout<<std::setprecision(17)
				<<"native_darcy: PASS tetrahedra="<<input.mesh.cells.size()
				<<" vertices="<<input.mesh.points.size()
				<<" volume_source_m3_s="<<solution.volume_source_m3_s
				<<" linear_iterations="<<solution.linear_iterations<<'\n';
		}
	}catch(const std::exception& error){
		if(rank==0)std::cerr<<"native_tet_darcy: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();
	return status;
}

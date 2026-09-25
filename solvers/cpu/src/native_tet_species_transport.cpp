#include "CaseConfig.hpp"
#include "CollectiveFailure.hpp"
#include "NativeTetMovingSpeciesCheckpoint.hpp"
#include "NativeTetMovingSpeciesPetscRuntime.hpp"
#include "NativeTetSpeciesVisualization.hpp"
#include "NativeTetWallReservoirExchange.hpp"
#include "ParallelVtkOutput.hpp"
#include "Sha256.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include <utility>

#include <cerrno>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace {

using Object=std::map<std::string,iga::config_detail::JsonValue>;

const iga::config_detail::JsonValue& Field(const Object& object,
	const std::string& key,const std::string& context)
{
	const auto* value=iga::config_detail::Find(object,key);
	if(!value)throw std::invalid_argument(context+" requires "+key);
	return *value;
}

double Number(const Object& object,const std::string& key,const std::string& context)
{
	return iga::config_detail::RequireNumber(Field(object,key,context),context+"."+key);
}

std::array<double,3> Vector3(const Object& object,const std::string& key)
{
	const auto& values=iga::config_detail::RequireArray(Field(object,key,"species case"),key);
	if(values.size()!=3)throw std::invalid_argument(key+" requires three SI components");
	std::array<double,3> result{};
	for(std::size_t component=0;component<3;++component)
		result[component]=iga::config_detail::RequireNumber(values[component],key);
	return result;
}

int BoundaryLabel(const std::string& text)
{
	std::size_t consumed=0;
	int label=-1;
	try{label=std::stoi(text,&consumed);}
	catch(const std::exception&){throw std::invalid_argument("native species boundary label is invalid");}
	if(label<0||consumed!=text.size()||std::to_string(label)!=text)
		throw std::invalid_argument("native species boundary label is invalid");
	return label;
}

struct Case
{
	std::filesystem::path mesh_file;
	std::string species_id;
	std::array<double,3> fluid_velocity{},mesh_velocity{};
	std::map<int,double> inflow;
	std::map<int,iga::NativeTetWallExchange> wall_exchange;
	std::optional<iga::NativeTetWallReservoirModel> finite_wall_reservoir;
	std::map<int,iga::NativeTetWallReservoirRegion> finite_wall_reservoirs;
	double initial_reservoir_concentration=0.;
	double initial_concentration=0.,diffusivity=0.,source=0.,decay_rate=0.,dt=0.;
	int steps=0;
	bool monotone=false;
	bool HasTissue() const
	{
		return finite_wall_reservoir.has_value()||!finite_wall_reservoirs.empty();
	}
};

Case ParseCase(const std::string& text,const std::filesystem::path& path)
{
	using namespace iga::config_detail;
	const auto parsed=JsonParser(text).Parse();
	const auto& root=RequireObject(parsed,"species case");
	RequireKnownKeys(root,{"schema_version","mesh_file","species_id",
		"initial_concentration_mol_m3","diffusivity_m2_s","source_mol_m3_s",
		"first_order_decay_rate_s_inv",
		"fluid_velocity_m_s","mesh_velocity_m_s",
		"inflow_concentration_by_label_mol_m3","wall_exchange_by_label",
		"finite_wall_reservoir","finite_wall_reservoirs_by_label",
		"monotone","time"},"species case");
	if(RequireInteger(Field(root,"schema_version","species case"),"schema_version")!=1)
		throw std::invalid_argument("unsupported native species case schema");
	Case result;
	result.mesh_file=path.parent_path()/RequireString(Field(root,"mesh_file","species case"),
		"mesh_file");
	result.species_id=RequireString(Field(root,"species_id","species case"),"species_id");
	result.initial_concentration=Number(root,"initial_concentration_mol_m3","species case");
	result.diffusivity=Number(root,"diffusivity_m2_s","species case");
	result.source=Number(root,"source_mol_m3_s","species case");
	if(Find(root,"first_order_decay_rate_s_inv"))
		result.decay_rate=Number(root,"first_order_decay_rate_s_inv","species case");
	result.fluid_velocity=Vector3(root,"fluid_velocity_m_s");
	result.mesh_velocity=Vector3(root,"mesh_velocity_m_s");
	result.monotone=RequireBoolean(Field(root,"monotone","species case"),"monotone");
	const auto& inflow=RequireObject(Field(root,"inflow_concentration_by_label_mol_m3",
		"species case"),"inflow_concentration_by_label_mol_m3");
	for(const auto& item:inflow){
		const int label=BoundaryLabel(item.first);
		const double concentration=RequireNumber(item.second,"inflow concentration");
		if(concentration<0.)throw std::invalid_argument("native species inflow is negative");
		result.inflow.emplace(label,concentration);
	}
	if(const auto* value=Find(root,"wall_exchange_by_label")){
		const auto& boundaries=RequireObject(*value,"wall_exchange_by_label");
		for(const auto& item:boundaries){
			const int label=BoundaryLabel(item.first);
			const auto& exchange=RequireObject(item.second,"wall exchange");
			RequireKnownKeys(exchange,{"transfer_coefficient_m_s",
				"external_concentration_mol_m3"},"wall exchange");
			iga::NativeTetWallExchange condition;
			condition.transfer_coefficient_m_s=Number(exchange,
				"transfer_coefficient_m_s","wall exchange");
			condition.external_concentration_mol_m3=Number(exchange,
				"external_concentration_mol_m3","wall exchange");
			if(condition.transfer_coefficient_m_s<0.
				||condition.external_concentration_mol_m3<0.)
				throw std::invalid_argument("native species wall exchange is negative");
			result.wall_exchange.emplace(label,condition);
		}
	}
	if(const auto* value=Find(root,"finite_wall_reservoir")){
		if(Find(root,"wall_exchange_by_label")||Find(root,"finite_wall_reservoirs_by_label"))
			throw std::invalid_argument("finite wall reservoir cannot combine with prescribed wall exchange");
		const auto& reservoir=RequireObject(*value,"finite_wall_reservoir");
		RequireKnownKeys(reservoir,{"boundary_label","transfer_coefficient_m_s",
			"volume_m3","initial_concentration_mol_m3"},"finite_wall_reservoir");
		iga::NativeTetWallReservoirModel model;
		model.boundary_label=RequireInteger(Field(reservoir,"boundary_label",
			"finite_wall_reservoir"),"finite_wall_reservoir.boundary_label");
		model.transfer_coefficient_m_s=Number(reservoir,
			"transfer_coefficient_m_s","finite_wall_reservoir");
		model.volume_m3=Number(reservoir,"volume_m3","finite_wall_reservoir");
		result.initial_reservoir_concentration=Number(reservoir,
			"initial_concentration_mol_m3","finite_wall_reservoir");
		if(model.boundary_label<0||!(model.transfer_coefficient_m_s>0.)
			||!(model.volume_m3>0.)||result.initial_reservoir_concentration<0.
			||!std::isfinite(model.volume_m3*result.initial_reservoir_concentration))
			throw std::invalid_argument("native species finite wall reservoir is invalid");
		result.finite_wall_reservoir=model;
	}
	if(const auto* value=Find(root,"finite_wall_reservoirs_by_label")){
		if(Find(root,"wall_exchange_by_label"))
			throw std::invalid_argument("finite wall reservoirs cannot combine with prescribed wall exchange");
		const auto& regions=RequireObject(*value,"finite_wall_reservoirs_by_label");
		if(regions.empty()||regions.size()>32)
			throw std::invalid_argument("native species finite wall reservoir count is invalid");
		for(const auto& item:regions){
			const int label=BoundaryLabel(item.first);
			const auto& input=RequireObject(item.second,"finite wall reservoir region");
			RequireKnownKeys(input,{"transfer_coefficient_m_s","volume_m3",
				"initial_concentration_mol_m3"},"finite wall reservoir region");
			iga::NativeTetWallReservoirModel model;
			model.boundary_label=label;
			model.transfer_coefficient_m_s=Number(input,"transfer_coefficient_m_s",
				"finite wall reservoir region");
			model.volume_m3=Number(input,"volume_m3","finite wall reservoir region");
			const double concentration=Number(input,"initial_concentration_mol_m3",
				"finite wall reservoir region");
			const double amount=model.volume_m3*concentration;
			if(!(model.transfer_coefficient_m_s>0.)||!(model.volume_m3>0.)
				||concentration<0.||!std::isfinite(amount))
				throw std::invalid_argument("native species finite wall reservoir region is invalid");
			result.finite_wall_reservoirs.emplace(label,
				iga::NativeTetWallReservoirRegion{model,{amount}});
		}
	}
	const auto& time=RequireObject(Field(root,"time","species case"),"time");
	RequireKnownKeys(time,{"dt_s","steps"},"time");
	result.dt=Number(time,"dt_s","time");
	result.steps=RequireInteger(Field(time,"steps","time"),"time.steps");
	if(result.mesh_file.empty()||result.species_id.empty()
		||result.initial_concentration<0.||result.diffusivity<0.
		||result.decay_rate<0.
		||!(result.dt>0.)||result.steps<1||result.steps>100000)
		throw std::invalid_argument("native species case has invalid material, time or input");
	for(const unsigned char character:result.species_id)
		if(!((character>='a'&&character<='z')
			||(character>='A'&&character<='Z')
			||(character>='0'&&character<='9')
			||character=='_'||character=='-'||character=='.'))
			throw std::invalid_argument("native species id must be an ASCII identifier");
	return result;
}

std::string Hash(const std::string& contents)
{
	iga::Sha256 hash;hash.Append(contents.data(),contents.size());return hash.Hex();
}

struct PublishedStep
{
	int step=0,ranks=0;
	double time_s=0.,inventory_mol=0.,balance_defect_mol_s=0.;
	std::string case_hash,mesh_hash,checkpoint_hash,tissue_checkpoint_hash;
};

std::filesystem::path CheckpointPath(const std::filesystem::path& output,int step)
{
	return output/("checkpoint_step_"+std::to_string(step)+".bin");
}

std::filesystem::path TissueCheckpointPath(const std::filesystem::path& output,int step)
{
	return output/("tissue_step_"+std::to_string(step)+".bin");
}

std::string SerializeTissueCheckpoint(int step,double time_s,double amount_mol,
	const std::string& case_hash,const std::string& mesh_hash)
{
	if(step<1||!(time_s>0.)||!std::isfinite(time_s)
		||!(amount_mol>=0.)||!std::isfinite(amount_mol))
		throw std::invalid_argument("native species tissue checkpoint state is invalid");
	iga::checkpoint_metadata::Writer payload;
	payload.Text("NativeTetWallReservoirState/v1");
	payload.Unsigned(static_cast<std::uint64_t>(step));
	payload.Real(time_s);payload.Real(amount_mol);
	payload.Text(case_hash);payload.Text(mesh_hash);
	iga::checkpoint_metadata::Writer checksum;
	checksum.Text(Hash(payload.Bytes()));
	return payload.Bytes()+checksum.Bytes();
}

double ParseTissueCheckpoint(const std::string& bytes,int step,double time_s,
	const std::string& case_hash,const std::string& mesh_hash)
{
	constexpr std::size_t checksum_bytes=8+64;
	if(bytes.size()<checksum_bytes||bytes.size()>iga::checkpoint_metadata::maximum_bytes)
		throw std::runtime_error("native species tissue checkpoint size is invalid");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);
	iga::checkpoint_metadata::Reader checksum(std::string_view(bytes).substr(
		bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();
	if(expected!=Hash(payload))
		throw std::runtime_error("native species tissue checkpoint checksum differs");
	iga::checkpoint_metadata::Reader input(payload);
	if(input.Text()!="NativeTetWallReservoirState/v1"
		||input.Unsigned()!=static_cast<std::uint64_t>(step))
		throw std::runtime_error("native species tissue checkpoint schema or step differs");
	const double saved_time=input.Real(),amount=input.Real();
	if(input.Text()!=case_hash||input.Text()!=mesh_hash)
		throw std::runtime_error("native species tissue checkpoint model differs");
	input.Finish();
	if(std::abs(saved_time-time_s)>1e-12*std::max(1.,std::abs(time_s))
		||!(amount>=0.))
		throw std::runtime_error("native species tissue checkpoint clock or amount differs");
	return amount;
}

std::string SerializeTissueRegionsCheckpoint(int step,double time_s,
	const std::map<int,iga::NativeTetWallReservoirRegion>& regions,
	const std::string& case_hash,const std::string& mesh_hash)
{
	if(step<1||!(time_s>0.)||!std::isfinite(time_s)
		||regions.empty()||regions.size()>32)
		throw std::invalid_argument("native species tissue regions checkpoint state is invalid");
	iga::checkpoint_metadata::Writer payload;
	payload.Text("NativeTetWallReservoirRegions/v2");
	payload.Unsigned(static_cast<std::uint64_t>(step));
	payload.Real(time_s);
	payload.Text(case_hash);payload.Text(mesh_hash);
	payload.Unsigned(static_cast<std::uint64_t>(regions.size()));
	for(const auto& item:regions){
		if(item.first<0||item.first!=item.second.model.boundary_label
			||!(item.second.committed.amount_mol>=0.)
			||!std::isfinite(item.second.committed.amount_mol))
			throw std::invalid_argument("native species tissue region checkpoint is invalid");
		payload.Unsigned(static_cast<std::uint64_t>(item.first));
		payload.Real(item.second.committed.amount_mol);
	}
	iga::checkpoint_metadata::Writer checksum;
	checksum.Text(Hash(payload.Bytes()));
	return payload.Bytes()+checksum.Bytes();
}

std::map<int,double> ParseTissueRegionsCheckpoint(const std::string& bytes,
	int step,double time_s,const std::string& case_hash,const std::string& mesh_hash,
	const std::map<int,iga::NativeTetWallReservoirRegion>& regions)
{
	constexpr std::size_t checksum_bytes=8+64;
	if(bytes.size()<checksum_bytes||bytes.size()>iga::checkpoint_metadata::maximum_bytes)
		throw std::runtime_error("native species tissue regions checkpoint size is invalid");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);
	iga::checkpoint_metadata::Reader checksum(std::string_view(bytes).substr(
		bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();
	if(expected!=Hash(payload))
		throw std::runtime_error("native species tissue regions checkpoint checksum differs");
	iga::checkpoint_metadata::Reader input(payload);
	if(input.Text()!="NativeTetWallReservoirRegions/v2"
		||input.Unsigned()!=static_cast<std::uint64_t>(step))
		throw std::runtime_error("native species tissue regions checkpoint schema or step differs");
	const double saved_time=input.Real();
	if(input.Text()!=case_hash||input.Text()!=mesh_hash
		||std::abs(saved_time-time_s)>1e-12*std::max(1.,std::abs(time_s))
		||input.Unsigned()!=regions.size())
		throw std::runtime_error("native species tissue regions checkpoint identity differs");
	std::map<int,double> amounts;
	for(const auto& item:regions){
		if(input.Unsigned()!=static_cast<std::uint64_t>(item.first))
			throw std::runtime_error("native species tissue regions checkpoint labels differ");
		const double amount=input.Real();
		if(!(amount>=0.))
			throw std::runtime_error("native species tissue region amount is invalid");
		amounts.emplace(item.first,amount);
	}
	input.Finish();
	return amounts;
}

std::string ReadBinary(const std::filesystem::path& path)
{
	if(!std::filesystem::is_regular_file(path)
		||std::filesystem::file_size(path)>iga::checkpoint_metadata::maximum_bytes)
		throw std::runtime_error("native species checkpoint file is missing or too large");
	std::ifstream input(path,std::ios::binary);
	if(!input)throw std::runtime_error("cannot open native species checkpoint");
	const std::string bytes(std::istreambuf_iterator<char>(input),{});
	if(!input.eof()&&input.fail())throw std::runtime_error("cannot read native species checkpoint");
	return bytes;
}

PublishedStep ReadPublishedStep(const std::filesystem::path& output)
{
	std::ifstream input(output/"run_state.json");
	if(!input)throw std::runtime_error("native species resume state is missing");
	using namespace iga::config_detail;
	const auto parsed=JsonParser(iga::ReadCheckedText(input)).Parse();
	const auto& root=RequireObject(parsed,"native species run state");
	RequireKnownKeys(root,{"schema_version","accepted_steps","mpi_ranks","time_s",
		"inventory_mol","balance_defect_mol_s","case_sha256","mesh_sha256",
		"checkpoint_sha256","tissue_checkpoint_sha256"},"native species run state");
	if(RequireInteger(Field(root,"schema_version","run state"),"schema_version")!=1)
		throw std::runtime_error("unsupported native species run state schema");
	PublishedStep state;
	state.step=RequireInteger(Field(root,"accepted_steps","run state"),"accepted_steps");
	state.ranks=RequireInteger(Field(root,"mpi_ranks","run state"),"mpi_ranks");
	state.time_s=Number(root,"time_s","run state");
	state.inventory_mol=Number(root,"inventory_mol","run state");
	state.balance_defect_mol_s=Number(root,"balance_defect_mol_s","run state");
	state.case_hash=RequireString(Field(root,"case_sha256","run state"),"case_sha256");
	state.mesh_hash=RequireString(Field(root,"mesh_sha256","run state"),"mesh_sha256");
	state.checkpoint_hash=RequireString(Field(root,"checkpoint_sha256","run state"),
		"checkpoint_sha256");
	if(const auto* tissue=Find(root,"tissue_checkpoint_sha256"))
		state.tissue_checkpoint_hash=RequireString(*tissue,"tissue_checkpoint_sha256");
	return state;
}

void PublishStep(const std::filesystem::path& output,const PublishedStep& state,
	const std::string& checkpoint,const std::string& tissue_checkpoint={})
{
	if((state.tissue_checkpoint_hash.empty())!=tissue_checkpoint.empty()
		||(!tissue_checkpoint.empty()&&Hash(tissue_checkpoint)!=state.tissue_checkpoint_hash))
		throw std::runtime_error("native species tissue publication hash differs");
	const auto path=CheckpointPath(output,state.step);
	const auto pending=std::filesystem::path(path.string()+".pending");
	if(std::filesystem::exists(path)||std::filesystem::exists(pending))
		throw std::runtime_error("native species checkpoint would overwrite an existing file");
	{
		std::ofstream stream(pending,std::ios::binary);
		if(!stream)throw std::runtime_error("cannot create native species checkpoint");
		stream.write(checkpoint.data(),static_cast<std::streamsize>(checkpoint.size()));
		stream.close();
		if(!stream)throw std::runtime_error("cannot write native species checkpoint");
	}
	std::filesystem::rename(pending,path);
	if(!tissue_checkpoint.empty()){
		const auto tissue_path=TissueCheckpointPath(output,state.step);
		const auto tissue_pending=std::filesystem::path(tissue_path.string()+".pending");
		if(std::filesystem::exists(tissue_path)||std::filesystem::exists(tissue_pending))
			throw std::runtime_error("native species tissue checkpoint would overwrite");
		std::ofstream stream(tissue_pending,std::ios::binary);
		if(!stream)throw std::runtime_error("cannot create native species tissue checkpoint");
		stream.write(tissue_checkpoint.data(),
			static_cast<std::streamsize>(tissue_checkpoint.size()));
		stream.close();
		if(!stream)throw std::runtime_error("cannot write native species tissue checkpoint");
		std::filesystem::rename(tissue_pending,tissue_path);
	}
	const auto manifest=output/"run_state.json";
	const auto temporary=output/"run_state.json.pending";
	if(std::filesystem::exists(temporary))
		throw std::runtime_error("native species pending run state already exists");
	std::ofstream report(temporary);
	if(!report)throw std::runtime_error("cannot create native species run state");
	report<<std::setprecision(17)
		<<"{\n  \"schema_version\": 1,\n"
		<<"  \"accepted_steps\": "<<state.step<<",\n"
		<<"  \"mpi_ranks\": "<<state.ranks<<",\n"
		<<"  \"time_s\": "<<state.time_s<<",\n"
		<<"  \"inventory_mol\": "<<state.inventory_mol<<",\n"
		<<"  \"balance_defect_mol_s\": "<<state.balance_defect_mol_s<<",\n"
		<<"  \"case_sha256\": \""<<state.case_hash<<"\",\n"
		<<"  \"mesh_sha256\": \""<<state.mesh_hash<<"\",\n"
		<<"  \"checkpoint_sha256\": \""<<state.checkpoint_hash<<"\"";
	if(!state.tissue_checkpoint_hash.empty())
		report<<",\n  \"tissue_checkpoint_sha256\": \""
			<<state.tissue_checkpoint_hash<<"\"";
	report<<"\n}\n";
	report.close();
	if(!report)throw std::runtime_error("cannot write native species run state");
	std::filesystem::rename(temporary,manifest);
}

struct OutputLock
{
	int descriptor=-1;
	~OutputLock(){if(descriptor>=0){flock(descriptor,LOCK_UN);close(descriptor);}}
	void Acquire(const std::filesystem::path& output)
	{
		const auto path=output/"run.lock";
		descriptor=open(path.c_str(),O_CREAT|O_RDWR,0600);
		if(descriptor<0||flock(descriptor,LOCK_EX|LOCK_NB)!=0)
			throw std::runtime_error("native species output is already locked or inaccessible");
	}
};

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int exit_code=0;
	try{
		if(argc<3)
			throw std::invalid_argument("usage: native_tet_species_transport case.json --check-input | --output-dir DIR [--stop-after-step N] | --resume DIR");
		const std::filesystem::path case_path=argv[1];
		bool check_only=false,output_mode=false,resume=false;
		std::filesystem::path output;
		int stop_after=0;
		bool stop_requested=false;
		for(int index=2;index<argc;++index){
			const std::string option=argv[index];
			if(option=="--check-input"&&!check_only&&!output_mode&&!resume){
				check_only=true;
			}else if((option=="--output-dir"||option=="--resume")
				&&index+1<argc&&!check_only&&!output_mode&&!resume){
				output_mode=option=="--output-dir";
				resume=option=="--resume";
				output=argv[++index];
			}else if(option=="--stop-after-step"&&index+1<argc
				&&output_mode&&!stop_requested){
				stop_requested=true;
				try{
					std::size_t consumed=0;
					const std::string value=argv[++index];
					stop_after=std::stoi(value,&consumed);
					if(consumed!=value.size())
						throw std::invalid_argument("native species stop step has trailing text");
				}
				catch(const std::exception&){
					throw std::invalid_argument("native species stop step is invalid");
				}
			}else throw std::invalid_argument("native species command option is invalid");
		}
		if((!check_only&&!output_mode&&!resume)||((output_mode||resume)&&output.empty()))
			throw std::invalid_argument("native species requires input check, new output, or resume");
		if(check_only&&argc!=3)
			throw std::invalid_argument("native species input check has extra options");
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		std::string case_text,mesh_text;
		Case specification;
		iga::NativeTetMesh reference;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species input",[&]{
			std::ifstream source(case_path);
			if(!source)throw std::runtime_error("cannot open native species case");
			case_text=iga::ReadCheckedText(source);
			specification=ParseCase(case_text,case_path);
			std::ifstream mesh_source(specification.mesh_file);
			if(!mesh_source)throw std::runtime_error("cannot open native species tetra mesh");
			mesh_text=iga::ReadCheckedText(mesh_source);
			std::istringstream stream(mesh_text);
			reference=iga::ReadNativeTetMeshGmsh41(stream);
			std::set<int> labels;
			for(const auto& face:reference.boundary_triangles)
				labels.insert(face.boundary_label);
			for(const auto& item:specification.inflow)
				if(!labels.count(item.first))
					throw std::invalid_argument("native species inflow label is absent");
			const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
			auto next=reference;
			for(auto& point:next.points)
				for(int axis=0;axis<3;++axis)
					point[axis]+=specification.mesh_velocity[axis]*specification.dt;
			const std::vector<std::array<double,3>> fluid(
				reference.points.size()+topology.edges.size(),specification.fluid_velocity);
			const std::vector<std::array<double,3>> grid(reference.points.size(),
				specification.mesh_velocity);
			const std::vector<double> initial(reference.points.size(),
				specification.initial_concentration);
			auto preflight_walls=specification.wall_exchange;
			if(specification.finite_wall_reservoir){
				const auto& model=*specification.finite_wall_reservoir;
				preflight_walls.emplace(model.boundary_label,iga::NativeTetWallExchange{
					model.transfer_coefficient_m_s,
					specification.initial_reservoir_concentration});
			}
			for(const auto& item:specification.finite_wall_reservoirs)
				preflight_walls.emplace(item.first,iga::NativeTetWallExchange{
					item.second.model.transfer_coefficient_m_s,
					item.second.committed.amount_mol/item.second.model.volume_m3});
			(void)iga::AssembleNativeTetMovingSpeciesStep(reference,next,fluid,grid,initial,
				specification.inflow,specification.diffusivity,specification.source,
				specification.dt,0,1,specification.monotone,specification.decay_rate,
				preflight_walls);
		});
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native species case identity",
			Hash(case_text));
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native species mesh identity",
			Hash(mesh_text));
		const auto case_hash=Hash(case_text),mesh_hash=Hash(mesh_text);
		if(check_only){
			if(rank==0)std::cout<<"native_species_input: PASS tetrahedra="
				<<reference.cells.size()<<" vertices="<<reference.points.size()<<'\n';
		}else{
			if(stop_requested&&(stop_after<1||stop_after>specification.steps))
				throw std::invalid_argument("native species stop step is outside the case time range");
			iga::RequireCollectiveSameText(PETSC_COMM_WORLD,"native species output path",
				std::filesystem::absolute(output).lexically_normal().string());
			OutputLock lock;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species output directory",[&]{
				if(rank!=0)return;
				if(resume){
					if(!std::filesystem::is_directory(output))
						throw std::runtime_error("native species resume directory is missing");
				}else{
					if(!output.parent_path().empty())
						std::filesystem::create_directories(output.parent_path());
					if(!std::filesystem::create_directory(output))
						throw std::runtime_error("native species output directory already exists");
				}
				lock.Acquire(output);
			});
			const auto topology=iga::BuildNativeTaylorHoodTopology(reference);
			const std::vector<std::array<double,3>> fluid(
				reference.points.size()+topology.edges.size(),specification.fluid_velocity);
			const std::vector<std::array<double,3>> grid(reference.points.size(),
				specification.mesh_velocity);
			auto committed=iga::InitializeNativeTetMovingSpeciesState(reference,
				std::vector<double>(reference.points.size(),specification.initial_concentration),
				specification.diffusivity);
			iga::NativeTetWallReservoirState tissue;
			auto tissue_regions=specification.finite_wall_reservoirs;
			if(specification.finite_wall_reservoir)
				tissue.amount_mol=specification.finite_wall_reservoir->volume_m3
					*specification.initial_reservoir_concentration;
			double last_inventory=0.,last_balance_defect=0.,last_reaction_sink=0.,
				last_wall_exchange=0.,last_combined_defect=0.;
			int first_step=1;
			if(resume){
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species resume validation",[&]{
					const auto published=ReadPublishedStep(output);
					const bool has_tissue=specification.HasTissue();
					if(published.step<1||published.step>=specification.steps
						||published.ranks!=ranks||published.case_hash!=case_hash
						||published.mesh_hash!=mesh_hash
						||published.tissue_checkpoint_hash.empty()==has_tissue
						||std::abs(published.time_s-published.step*specification.dt)
							>1e-12*std::max(1.,std::abs(published.time_s))
						||std::filesystem::exists(output/"run_summary.json")
						||std::filesystem::exists(output/"run_state.json.pending")
						||std::filesystem::exists(output/("step_"+std::to_string(published.step+1)))
						||std::filesystem::exists(CheckpointPath(output,published.step+1))
						||std::filesystem::exists(TissueCheckpointPath(output,published.step+1))
						||std::filesystem::exists(std::filesystem::path(
							TissueCheckpointPath(output,published.step+1).string()+".pending"))
						||std::filesystem::exists(std::filesystem::path(
							CheckpointPath(output,published.step+1).string()+".pending")))
						throw std::runtime_error("native species resume identity, clock or publication differs");
					for(int step=1;step<=published.step;++step){
						const auto snapshot=output/("step_"+std::to_string(step));
						if(!std::filesystem::is_regular_file(snapshot/"snapshot.pvtu")
							||!std::filesystem::is_regular_file(CheckpointPath(output,step)))
							throw std::runtime_error("native species resume is missing a published step");
						if(std::filesystem::is_regular_file(TissueCheckpointPath(output,step))
							!=has_tissue
							||std::filesystem::exists(std::filesystem::path(
								TissueCheckpointPath(output,step).string()+".pending")))
							throw std::runtime_error("native species resume tissue publication differs");
						for(int peer=0;peer<ranks;++peer)
							if(!std::filesystem::is_regular_file(snapshot/(
								"rank"+std::to_string(peer)+".vtu")))
								throw std::runtime_error("native species resume is missing a VTU rank piece");
					}
					const auto bytes=ReadBinary(CheckpointPath(output,published.step));
					if(Hash(bytes)!=published.checkpoint_hash)
						throw std::runtime_error("native species resume checkpoint hash differs");
					const auto checkpoint=iga::ParseNativeTetMovingSpeciesCheckpoint(bytes);
					if(checkpoint.accepted_steps!=static_cast<std::uint64_t>(published.step)
						||std::abs(checkpoint.time_s-published.time_s)>1e-12)
						throw std::runtime_error("native species resume checkpoint clock differs");
					committed=iga::RestoreNativeTetMovingSpeciesCheckpoint(checkpoint,
						reference,specification.diffusivity);
					if(has_tissue){
						const auto tissue_bytes=ReadBinary(TissueCheckpointPath(output,published.step));
						if(Hash(tissue_bytes)!=published.tissue_checkpoint_hash)
							throw std::runtime_error("native species tissue checkpoint hash differs");
						if(specification.finite_wall_reservoir)
							tissue.amount_mol=ParseTissueCheckpoint(tissue_bytes,published.step,
								published.time_s,case_hash,mesh_hash);
						else{
							const auto amounts=ParseTissueRegionsCheckpoint(tissue_bytes,
								published.step,published.time_s,case_hash,mesh_hash,tissue_regions);
							for(const auto& item:amounts)
								tissue_regions.at(item.first).committed.amount_mol=item.second;
						}
					}
					for(std::size_t node=0;node<reference.points.size();++node)
						for(int axis=0;axis<3;++axis)
							if(std::abs(committed.current_mesh.points[node][axis]
								-reference.points[node][axis]
								-published.step*specification.dt*specification.mesh_velocity[axis])
								>1e-12*std::max(1.,std::abs(reference.points[node][axis])))
								throw std::runtime_error("native species resume ALE motion differs");
					first_step=published.step+1;
					last_inventory=published.inventory_mol;
					last_balance_defect=published.balance_defect_mol_s;
				});
			}
			const int final_step=stop_requested?stop_after:specification.steps;
			if(final_step<first_step)
				throw std::invalid_argument("native species requested stop precedes resumed state");
			for(int step=first_step;step<=final_step;++step){
				auto current=committed.current_mesh;
				for(auto& point:current.points)
					for(int axis=0;axis<3;++axis)
						point[axis]+=specification.mesh_velocity[axis]*specification.dt;
				iga::NativeTetMovingSpeciesPetscResult solved;
				std::optional<iga::NativeTetWallReservoirResult> tissue_result;
				std::optional<iga::NativeTetWallReservoirsResult> regions_result;
				if(specification.finite_wall_reservoir){
					tissue_result=iga::SolveNativeTetWallReservoirStep(
						committed.current_mesh,current,fluid,grid,
						committed.concentration_mol_m3,specification.inflow,
						specification.diffusivity,specification.source,specification.dt,
						*specification.finite_wall_reservoir,tissue,
						specification.monotone,specification.decay_rate);
					solved=tissue_result->vessel;
				}else if(!tissue_regions.empty()){
					regions_result=iga::SolveNativeTetWallReservoirsStep(
						committed.current_mesh,current,fluid,grid,
						committed.concentration_mol_m3,specification.inflow,
						specification.diffusivity,specification.source,specification.dt,
						tissue_regions,specification.monotone,specification.decay_rate);
					solved=regions_result->vessel;
				}else solved=iga::SolveNativeTetMovingSpeciesPetscStep(
					committed.current_mesh,current,fluid,grid,
					committed.concentration_mol_m3,specification.inflow,
					specification.diffusivity,specification.source,specification.dt,
					specification.monotone,specification.decay_rate,
					specification.wall_exchange);
				const double scale=std::max({std::abs(solved.step.previous_inventory_mol)
					/specification.dt,std::abs(solved.step.current_inventory_mol)
					/specification.dt,std::abs(solved.step.outward_advective_flux_mol_s),
					std::abs(solved.step.source_mol_s),
					std::abs(solved.step.reaction_sink_mol_s),
					std::abs(solved.step.outward_wall_exchange_mol_s)});
				if(std::abs(solved.step.balance_defect_mol_s)>1e-12+1e-8*scale)
					throw std::runtime_error("native species accepted-step balance failed");
				if(tissue_result&&std::abs(tissue_result->combined_balance_defect_mol_s)
					>1e-12+1e-8*scale)
					throw std::runtime_error("native species paired tissue balance failed");
				if(regions_result&&std::abs(regions_result->combined_balance_defect_mol_s)
					>1e-12+1e-8*scale)
					throw std::runtime_error("native species paired tissue regions balance failed");
				for(const double value:solved.step.concentration_mol_m3)
					if(!std::isfinite(value)||value<0.)
						throw std::runtime_error("native species concentration is invalid");
				auto candidate=committed;
				candidate.current_mesh=current;
				candidate.concentration_mol_m3=solved.step.concentration_mol_m3;
				candidate.accepted_steps=static_cast<std::uint64_t>(step);
				candidate.time_s=step*specification.dt;
				std::string checkpoint_bytes;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species checkpoint build",[&]{
					checkpoint_bytes=iga::SerializeNativeTetMovingSpeciesCheckpoint(
						iga::CaptureNativeTetMovingSpeciesCheckpoint(candidate));
				});
				iga::RequireCollectiveSameText(PETSC_COMM_WORLD,
					"native species checkpoint agreement",Hash(checkpoint_bytes));
				std::string tissue_bytes;
				if(tissue_result){
					tissue_bytes=SerializeTissueCheckpoint(step,candidate.time_s,
						tissue_result->reservoir.amount_mol,case_hash,mesh_hash);
				}else if(regions_result){
					auto candidate_regions=tissue_regions;
					for(const auto& item:regions_result->reservoirs)
						candidate_regions.at(item.first).committed=item.second;
					tissue_bytes=SerializeTissueRegionsCheckpoint(step,candidate.time_s,
						candidate_regions,case_hash,mesh_hash);
				}
				if(!tissue_bytes.empty()){
					iga::RequireCollectiveSameText(PETSC_COMM_WORLD,
						"native species tissue checkpoint agreement",Hash(tissue_bytes));
				}
				iga::VtkPartition piece;
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species field build",[&]{
					piece=iga::BuildNativeTetSpeciesVtkPartition(reference,current,
						solved.step.concentration_mol_m3,rank,ranks);
				});
				iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,
					output/("step_"+std::to_string(step)),piece,step*specification.dt);
				const PublishedStep published{step,ranks,candidate.time_s,
					solved.step.current_inventory_mol,solved.step.balance_defect_mol_s,
					case_hash,mesh_hash,Hash(checkpoint_bytes),
					tissue_bytes.empty()?std::string{}:Hash(tissue_bytes)};
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species step publication",[&]{
					if(rank==0)PublishStep(output,published,checkpoint_bytes,tissue_bytes);
				});
				if(rank==0)std::cout<<std::setprecision(17)
					<<"native_species_step step="<<step
					<<" time_s="<<step*specification.dt
					<<" inventory_mol="<<solved.step.current_inventory_mol
					<<" outward_flux_mol_s="<<solved.step.outward_advective_flux_mol_s
					<<" source_mol_s="<<solved.step.source_mol_s
					<<" reaction_sink_mol_s="<<solved.step.reaction_sink_mol_s
					<<" outward_wall_exchange_mol_s="
					<<solved.step.outward_wall_exchange_mol_s
					<<" balance_defect_mol_s="<<solved.step.balance_defect_mol_s
					<<" linear_iterations="<<solved.linear_iterations;
				if(rank==0&&tissue_result)std::cout<<std::setprecision(17)
					<<" tissue_amount_mol="<<tissue_result->reservoir.amount_mol
					<<" combined_balance_defect_mol_s="
					<<tissue_result->combined_balance_defect_mol_s;
				if(rank==0&&regions_result)std::cout<<std::setprecision(17)
					<<" tissue_regions="<<regions_result->reservoirs.size()
					<<" combined_balance_defect_mol_s="
					<<regions_result->combined_balance_defect_mol_s;
				if(rank==0)std::cout<<'\n';
				last_inventory=solved.step.current_inventory_mol;
				last_balance_defect=solved.step.balance_defect_mol_s;
				last_reaction_sink=solved.step.reaction_sink_mol_s;
				last_wall_exchange=solved.step.outward_wall_exchange_mol_s;
				if(tissue_result){
					last_combined_defect=tissue_result->combined_balance_defect_mol_s;
					tissue=tissue_result->reservoir;
				}
				if(regions_result){
					last_combined_defect=regions_result->combined_balance_defect_mol_s;
					for(const auto& item:regions_result->reservoirs)
						tissue_regions.at(item.first).committed=item.second;
				}
				committed=std::move(candidate);
			}
			if(final_step==specification.steps)
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"native species summary publication",[&]{
				if(rank!=0)return;
				const auto temporary=output/"run_summary.json.pending";
				std::ofstream report(temporary);
				if(!report)throw std::runtime_error("cannot create native species summary");
				report<<std::setprecision(17)
					<<"{\n  \"schema_version\": 1,\n"
					<<"  \"kind\": \"native_tet_species_prescribed_velocity\",\n"
					<<"  \"species_id\": \""<<specification.species_id<<"\",\n"
					<<"  \"case_sha256\": \""<<case_hash<<"\",\n"
					<<"  \"mesh_sha256\": \""<<mesh_hash<<"\",\n"
					<<"  \"mpi_ranks\": "<<ranks<<",\n"
					<<"  \"accepted_steps\": "<<specification.steps<<",\n"
					<<"  \"last_time_s\": "<<specification.steps*specification.dt<<",\n"
					<<"  \"last_inventory_mol\": "<<last_inventory<<",\n"
					<<"  \"last_reaction_sink_mol_s\": "<<last_reaction_sink<<",\n"
					<<"  \"last_outward_wall_exchange_mol_s\": "
					<<last_wall_exchange<<",\n"
					<<"  \"last_balance_defect_mol_s\": "<<last_balance_defect;
				if(specification.finite_wall_reservoir)
					report<<",\n  \"tissue_volume_m3\": "
						<<specification.finite_wall_reservoir->volume_m3
						<<",\n  \"last_tissue_amount_mol\": "<<tissue.amount_mol
						<<",\n  \"last_tissue_concentration_mol_m3\": "
						<<tissue.amount_mol/specification.finite_wall_reservoir->volume_m3
						<<",\n  \"last_combined_balance_defect_mol_s\": "
						<<last_combined_defect;
				if(!tissue_regions.empty()){
					report<<",\n  \"tissue_regions_by_label\": {";
					bool first=true;
					for(const auto& item:tissue_regions){
						if(!first)report<<",";
						first=false;
						report<<"\n    \""<<item.first<<"\": {\"volume_m3\": "
							<<item.second.model.volume_m3
							<<", \"amount_mol\": "<<item.second.committed.amount_mol
							<<", \"concentration_mol_m3\": "
							<<item.second.committed.amount_mol/item.second.model.volume_m3
							<<"}";
					}
					report<<"\n  },\n  \"last_combined_balance_defect_mol_s\": "
						<<last_combined_defect;
				}
				report<<"\n}\n";
				report.close();
				if(!report)throw std::runtime_error("cannot write native species summary");
				std::filesystem::rename(temporary,output/"run_summary.json");
			});
		}
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<"native_tet_species_transport: "<<error.what()<<'\n';
		exit_code=1;
	}
	PetscFinalize();return exit_code;
}

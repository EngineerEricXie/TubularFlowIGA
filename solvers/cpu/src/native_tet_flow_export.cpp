#include "CaseConfig.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "TemporalUnstructuredVtkHdf.hpp"
#include <cstring>
#include <iostream>
#include <memory>

// Offline conversion of native P2 velocity / P1 pressure states. No MPI or
// solver is initialized. The manifest declares every time and pressure offset.
int main(int argc,char** argv)
{
	try{
		if(argc!=3)throw std::invalid_argument("usage: native_tet_flow_export series.json output.vtkhdf");
		using namespace iga::config_detail;
		const std::filesystem::path manifest=std::filesystem::absolute(argv[1]);
		std::ifstream input(manifest);
		if(!input)throw std::runtime_error("cannot open native field series manifest");
		const std::string text((std::istreambuf_iterator<char>(input)),{});
		const auto parsed=JsonParser(text).Parse();
		const auto& root=RequireObject(parsed,"field series");
		RequireKnownKeys(root,{"mesh_file","dynamic_viscosity_pa_s","pressure_offset_pa","snapshots"},"field series");
		const auto member=[](const auto& object,const std::string& name)->const JsonValue&{
			const auto* value=Find(object,name);
			if(!value)throw std::invalid_argument("field series requires "+name);
			return *value;
		};
		const auto mesh_path=manifest.parent_path()/RequireString(member(root,"mesh_file"),"mesh_file");
		std::ifstream source(mesh_path);
		if(!source)throw std::runtime_error("cannot open native mesh");
		const auto mesh=iga::ReadNativeTetMeshGmsh41(source);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const auto velocity_dofs=3*(mesh.points.size()+topology.edges.size());
		const auto count=velocity_dofs+mesh.points.size();
		const double viscosity=RequireNumber(member(root,"dynamic_viscosity_pa_s"),"viscosity");
		const double pressure_offset=RequireNumber(member(root,"pressure_offset_pa"),"pressure offset");
		const auto& snapshots=RequireArray(member(root,"snapshots"),"snapshots");
		if(snapshots.empty())throw std::invalid_argument("native field series has no snapshots");
		std::unique_ptr<iga::TemporalUnstructuredVtkHdfWriter> writer;
		for(const auto& item:snapshots){
			const auto& snapshot=RequireObject(item,"snapshot");
			RequireKnownKeys(snapshot,{"time_s","state_file"},"snapshot");
			const auto path=manifest.parent_path()/RequireString(member(snapshot,"state_file"),"state_file");
			const double time=RequireNumber(member(snapshot,"time_s"),"time_s");
			std::ifstream state_file(path,std::ios::binary);
			if(!state_file||std::filesystem::file_size(path)!=count*8)
				throw std::runtime_error("native FP64 state size differs: "+path.string());
			std::vector<double> state(count);
			static_assert(sizeof(double)==8&&std::numeric_limits<double>::is_iec559,"IEEE FP64 required");
			std::vector<unsigned char> bytes(count*8);
			state_file.read(reinterpret_cast<char*>(bytes.data()),bytes.size());
			if(!state_file)throw std::runtime_error("cannot read complete native state");
			for(std::size_t index=0;index<count;++index){
				std::uint64_t bits=0;
				for(int byte=0;byte<8;++byte)bits|=std::uint64_t(bytes[8*index+byte])<<(8*byte);
				std::memcpy(&state[index],&bits,8);
				if(index>=velocity_dofs)state[index]+=pressure_offset;
			}
			const auto piece=iga::BuildNativeTetHydraulicVtkPartition(mesh,mesh,state,viscosity,0,1);
			if(!writer)writer=std::make_unique<iga::TemporalUnstructuredVtkHdfWriter>(argv[2],piece);
			writer->Append(time,piece);
			std::cout<<"exported_time_s="<<std::setprecision(17)<<time<<'\n';
		}
		writer->Close();
		std::cout<<"native_tet_flow_export: PASS snapshots="<<snapshots.size()
			<<" points="<<velocity_dofs/3<<" cells="<<mesh.cells.size()<<'\n';
		return 0;
	}catch(const std::exception& error){
		std::cerr<<"native_tet_flow_export: "<<error.what()<<'\n';return 1;
	}
}

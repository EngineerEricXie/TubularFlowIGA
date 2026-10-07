#ifndef IGA_PARALLEL_VTK_OUTPUT_HPP
#define IGA_PARALLEL_VTK_OUTPUT_HPP

#include "CollectiveFailure.hpp"
#include "PartitionedVtkOutput.hpp"
#include "Sha256.hpp"
#include <regex>

namespace iga {

// The caller supplies an owned-cell partition and globally assigned point ids.
// Only schema/path/time metadata are exchanged here, never mesh or field arrays.
// A snapshot directory is immutable and must not already exist. Failed attempts
// leave diagnostic pieces behind; retry with a fresh directory. All members of
// communicator enter together. MPI process loss and durable fsync are excluded.
inline void WriteParallelVtkSnapshot(MPI_Comm communicator,const std::filesystem::path& directory,
	const VtkPartition& piece,double physical_time)
{
	int rank=0,ranks=0;MPI_Comm_rank(communicator,&rank);MPI_Comm_size(communicator,&ranks);
	std::filesystem::path root,local_path;
	std::vector<VtkArraySchema> point_schema,cell_schema;
	std::string identity;
	CollectiveLocalStage(communicator,"parallel VTK preparation",[&] {
		if(directory.empty())throw std::invalid_argument("parallel VTK snapshot directory is empty");
		root=std::filesystem::absolute(directory).lexically_normal();
		ValidateVtkPartition(piece,physical_time);
		point_schema=partitioned_vtk_detail::Schema(piece.point_arrays,piece.point_ids.size());
		cell_schema=partitioned_vtk_detail::Schema(piece.cell_arrays,piece.cell_ids.size());
		local_path=root/("rank"+std::to_string(rank)+".vtu");
		Sha256 hash;
		const auto append=[&](const std::string& text) {
			hash.AppendLittleEndian64(text.size());hash.Append(text.data(),text.size());
		};
		append(root.string());hash.AppendNormalizedDouble(physical_time);
		for(const auto* schema:{&point_schema,&cell_schema}) {
			hash.AppendLittleEndian64(schema->size());
			for(const auto& array:*schema) { append(array.name);hash.AppendLittleEndian32(static_cast<std::uint32_t>(array.components)); }
		}
		identity=hash.Hex();
	});
	RequireCollectiveSameText(communicator,"parallel VTK schema/path/time agreement",identity);
	CollectiveLocalStage(communicator,"parallel VTK directory publication",[&] {
		if(rank!=0)return;
		std::filesystem::create_directories(root.parent_path());
		if(!std::filesystem::create_directory(root))throw std::runtime_error("parallel VTK snapshot directory already exists");
	});
	CollectiveLocalStage(communicator,"parallel VTK piece write",[&] {
		WriteVtuPartition(local_path,piece,physical_time);
	});
	// The successful piece-write agreement is the publication barrier. Only
	// root builds the O(ranks) manifest, after every stream has closed cleanly.
	CollectiveLocalStage(communicator,"parallel VTK index publication",[&] {
		if(rank!=0)return;
		std::vector<std::filesystem::path> names;names.reserve(static_cast<std::size_t>(ranks));
		for(int peer=0;peer<ranks;++peer)names.emplace_back("rank"+std::to_string(peer)+".vtu");
		const auto temporary=root/"snapshot.pvtu.pending";
		WritePvtu(temporary,names,point_schema,cell_schema);
		std::filesystem::rename(temporary,root/"snapshot.pvtu");
	});
}

// Read the collection format emitted below, preserving its exact times and
// references. Reject incomplete indexes before a resumed run can replace them.
inline std::vector<std::pair<double,std::filesystem::path>> ReadParallelVtkSeries(
	const std::filesystem::path& path)
{
	std::ifstream input(path);
	if(!input)throw std::runtime_error("cannot read parallel VTK series");
	std::string line;
	if(!std::getline(input,line)||line!="<?xml version=\"1.0\"?>"
		||!std::getline(input,line)||line!="<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>")
		throw std::runtime_error("invalid parallel VTK series header");
	const std::regex entry(R"vtk(<DataSet timestep="([^"]+)" group="" part="0" file="([^"]+)"/>)vtk");
	std::vector<std::pair<double,std::filesystem::path>> snapshots;
	bool complete=false;
	double previous=-std::numeric_limits<double>::infinity();
	while(std::getline(input,line)){
		if(line=="</Collection></VTKFile>"){complete=true;break;}
		std::smatch match;
		if(!std::regex_match(line,match,entry))throw std::runtime_error("invalid parallel VTK series entry");
		const auto clock=match[1].str();std::size_t consumed=0;
		const double time=std::stod(clock,&consumed);
		if(consumed!=clock.size()||!std::isfinite(time)||time<=previous)
			throw std::runtime_error("invalid parallel VTK series time");
		previous=time;
		auto name=match[2].str();
		// Decode ampersands last so an escaped literal entity stays literal.
		for(const auto& entity:std::vector<std::pair<std::string,std::string>>{
			{"&quot;","\""},{"&apos;","'"},{"&lt;","<"},{"&gt;",">"},{"&amp;","&"}}){
			std::size_t position=0;
			while((position=name.find(entity.first,position))!=std::string::npos){
				name.replace(position,entity.first.size(),entity.second);position+=entity.second.size();
			}
		}
		const std::filesystem::path relative=name;
		if(relative.empty()||relative.is_absolute())throw std::runtime_error("invalid parallel VTK series reference");
		for(const auto& part:relative)if(part=="..")throw std::runtime_error("parallel VTK series reference escapes directory");
		const auto snapshot=path.parent_path()/relative;
		if(!std::filesystem::is_regular_file(snapshot))throw std::runtime_error("parallel VTK series snapshot is missing");
		snapshots.push_back({time,snapshot});
	}
	if(!complete||snapshots.empty()||input.bad())throw std::runtime_error("incomplete parallel VTK series");
	while(std::getline(input,line))if(line.find_first_not_of(" \t\r")!=std::string::npos)
		throw std::runtime_error("trailing parallel VTK series content");
	if(input.bad())throw std::runtime_error("cannot read complete parallel VTK series");
	return snapshots;
}

// Publish a time index after its immutable snapshots are complete. A single
// caller owns the series path. A failed write leaves the previous index intact.
inline void WriteParallelVtkSeries(MPI_Comm comm,const std::filesystem::path& path,
	const std::vector<std::pair<double,std::filesystem::path>>& snapshots)
{
	int rank=0;MPI_Comm_rank(comm,&rank);std::string text;
	CollectiveLocalStage(comm,"parallel VTK series preparation",[&] {
		if(path.empty()||snapshots.empty())throw std::invalid_argument("parallel VTK series is empty");
		std::ostringstream output;output.exceptions(std::ios::badbit|std::ios::failbit);
		output<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"><Collection>\n";
		const auto parent=std::filesystem::absolute(path).parent_path().lexically_normal();
		double previous=-std::numeric_limits<double>::infinity();
		for(const auto& snapshot:snapshots) {
			if(!std::isfinite(snapshot.first)||snapshot.first<=previous)throw std::invalid_argument("parallel VTK times must increase");
			previous=snapshot.first;
			const auto relative=std::filesystem::absolute(snapshot.second).lexically_normal().lexically_relative(parent);
			if(relative.empty()||relative.is_absolute())throw std::invalid_argument("invalid parallel VTK series reference");
			for(const auto& component:relative)if(component=="..")throw std::invalid_argument("parallel VTK snapshot is outside series directory");
			for(unsigned char c:relative.string())if(c<32||c==127)throw std::invalid_argument("control character in parallel VTK series reference");
			output<<"<DataSet timestep=\""<<snapshot.first<<"\" group=\"\" part=\"0\" file=\""<<EscapeVtkXml(relative.generic_string())<<"\"/>\n";
		}
		output<<"</Collection></VTKFile>\n";text=output.str();
	});
	RequireCollectiveSameText(comm,"parallel VTK series path agreement",path.string());
	RequireCollectiveSameText(comm,"parallel VTK series content agreement",text);
	CollectiveLocalStage(comm,"parallel VTK series publication",[&] {
		if(rank!=0)return;
		const auto temporary=std::filesystem::path(path.string()+".pending");
		for(const auto& target:{path,temporary})if(std::filesystem::exists(target)&&!std::filesystem::is_regular_file(target))
			throw std::runtime_error("parallel VTK series target is not a regular file");
		if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
		std::ofstream output(temporary);if(!output)throw std::runtime_error("cannot create parallel VTK series");
		output<<text;output.close();if(!output)throw std::runtime_error("cannot write parallel VTK series");
		std::filesystem::rename(temporary,path);
	});
}

} // namespace iga
#endif

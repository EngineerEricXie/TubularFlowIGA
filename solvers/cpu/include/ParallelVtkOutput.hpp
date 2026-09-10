#ifndef IGA_PARALLEL_VTK_OUTPUT_HPP
#define IGA_PARALLEL_VTK_OUTPUT_HPP

#include "CollectiveFailure.hpp"
#include "PartitionedVtkOutput.hpp"
#include "Sha256.hpp"

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

} // namespace iga
#endif

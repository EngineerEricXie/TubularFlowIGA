#ifndef IGA_NATIVE_TET_MESH_PARTITION_HPP
#define IGA_NATIVE_TET_MESH_PARTITION_HPP

// Partition a tetrahedral mesh over an MPI communicator. Cells are assigned
// to ranks with PETSc's MatPartitioning (ParMETIS when available, otherwise
// contiguous chunks); a vertex belongs to the rank owning most of its cells.
// Each rank keeps a local submesh of its owned cells plus the halo cells that
// touch its owned vertices, numbered with owned vertices first, so that
// assembly over owned cells and evaluation on owned entities need no further
// communication. The global cell owner and vertex owner arrays are replicated
// (one int per entity); everything else is local.

#include "NativeTetFem.hpp"

#include <petscmat.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct NativeTetMeshPartition
{
	int rank=0,ranks=1;
	std::vector<int> cell_owner;                  // global cell -> rank
	std::vector<int> node_owner;                  // global vertex -> rank
	NativeTetMesh local;                          // owned + halo cells, local numbering
	std::vector<std::uint32_t> local_to_global_node;   // local vertex -> global vertex
	std::vector<std::uint32_t> local_to_global_cell;   // local cell -> global cell
	std::vector<std::uint8_t> local_cell_owned;        // 1 for owned cells
	std::vector<std::uint8_t> local_triangle_owned;    // 1 when the adjacent cell is owned
	std::uint32_t owned_node_count=0;                  // local vertices [0, count) are owned
	std::vector<std::uint32_t> owned_node_begin;       // per rank: first owned vertex in the
	                                                   // owner-major global vertex numbering

	// Vertex index in the owner-major global numbering (owned blocks per rank).
	std::uint32_t OwnerMajorNode(std::uint32_t global_node)const
	{
		return owner_major_node_[global_node];
	}
	bool OwnsNode(std::uint32_t local_node)const{return local_node<owned_node_count;}

	std::vector<std::uint32_t> owner_major_node_;      // global vertex -> owner-major index
};

namespace native_tet_partition_detail {

inline void Check(PetscErrorCode code,const char* context)
{
	if(code)throw std::runtime_error(std::string("native tetra partition: ")+context
		+" failed with PETSc error "+std::to_string(static_cast<int>(code)));
}

// Owner of every cell, replicated on all ranks.
inline std::vector<int> PartitionCells(const NativeTetMesh& mesh,MPI_Comm communicator)
{
	int rank=0,ranks=1;MPI_Comm_rank(communicator,&rank);MPI_Comm_size(communicator,&ranks);
	const std::size_t cells=mesh.cells.size();
	std::vector<int> owner(cells,0);
	if(ranks==1||cells<static_cast<std::size_t>(ranks))return owner;
	// Cells adjacent through a shared face, from a global face map.
	std::map<std::array<std::uint32_t,3>,std::array<std::int64_t,2>> faces;
	for(std::size_t cell=0;cell<cells;++cell){
		const auto& nodes=mesh.cells[cell].nodes;
		for(int opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> face{};
			int entry=0;
			for(int vertex=0;vertex<4;++vertex)if(vertex!=opposite)face[entry++]=nodes[vertex];
			std::sort(face.begin(),face.end());
			auto& adjacent=faces.try_emplace(face,std::array<std::int64_t,2>{{-1,-1}}).first->second;
			if(adjacent[0]<0)adjacent[0]=static_cast<std::int64_t>(cell);
			else if(adjacent[1]<0)adjacent[1]=static_cast<std::int64_t>(cell);
			else throw std::invalid_argument("native tetra partition: a face belongs to three cells");
		}
	}
	std::vector<std::vector<PetscInt>> adjacency(cells);
	for(const auto& face:faces)
		if(face.second[1]>=0){
			adjacency[static_cast<std::size_t>(face.second[0])].push_back(static_cast<PetscInt>(face.second[1]));
			adjacency[static_cast<std::size_t>(face.second[1])].push_back(static_cast<PetscInt>(face.second[0]));
		}
	faces.clear();
	// This rank's contiguous cell range of the distributed adjacency matrix.
	const std::size_t begin=cells*static_cast<std::size_t>(rank)/ranks;
	const std::size_t end=cells*static_cast<std::size_t>(rank+1)/ranks;
	std::vector<PetscInt> ia(end-begin+1,0),ja;
	for(std::size_t cell=begin;cell<end;++cell){
		auto& row=adjacency[cell];
		std::sort(row.begin(),row.end());
		ja.insert(ja.end(),row.begin(),row.end());
		ia[cell-begin+1]=static_cast<PetscInt>(ja.size());
	}
	adjacency.clear();
	// MatCreateMPIAdj takes ownership of copies; keep ours alive until applied.
	PetscInt* ia_copy=nullptr;PetscInt* ja_copy=nullptr;
	Check(PetscMalloc1(ia.size(),&ia_copy),"PetscMalloc1 ia");
	Check(PetscMalloc1(ja.size()+1,&ja_copy),"PetscMalloc1 ja");
	std::copy(ia.begin(),ia.end(),ia_copy);std::copy(ja.begin(),ja.end(),ja_copy);
	Mat adjacency_matrix=nullptr;
	Check(MatCreateMPIAdj(communicator,static_cast<PetscInt>(end-begin),static_cast<PetscInt>(cells),
		ia_copy,ja_copy,nullptr,&adjacency_matrix),"MatCreateMPIAdj");
	MatPartitioning partitioning=nullptr;
	Check(MatPartitioningCreate(communicator,&partitioning),"MatPartitioningCreate");
	Check(MatPartitioningSetAdjacency(partitioning,adjacency_matrix),"MatPartitioningSetAdjacency");
#ifdef PETSC_HAVE_PARMETIS
	Check(MatPartitioningSetType(partitioning,MATPARTITIONINGPARMETIS),"MatPartitioningSetType");
#endif
	Check(MatPartitioningSetFromOptions(partitioning),"MatPartitioningSetFromOptions");
	IS assignment=nullptr;
	Check(MatPartitioningApply(partitioning,&assignment),"MatPartitioningApply");
	const PetscInt* local_owner=nullptr;
	Check(ISGetIndices(assignment,&local_owner),"ISGetIndices");
	std::vector<int> local(end-begin);
	for(std::size_t cell=0;cell<local.size();++cell)local[cell]=static_cast<int>(local_owner[cell]);
	Check(ISRestoreIndices(assignment,&local_owner),"ISRestoreIndices");
	ISDestroy(&assignment);MatPartitioningDestroy(&partitioning);MatDestroy(&adjacency_matrix);
	std::vector<int> counts(ranks),displacements(ranks);
	for(int r=0;r<ranks;++r){
		const std::size_t first=cells*static_cast<std::size_t>(r)/ranks,last=cells*static_cast<std::size_t>(r+1)/ranks;
		counts[r]=static_cast<int>(last-first);displacements[r]=static_cast<int>(first);
	}
	if(MPI_Allgatherv(local.data(),static_cast<int>(local.size()),MPI_INT,owner.data(),counts.data(),
		displacements.data(),MPI_INT,communicator)!=MPI_SUCCESS)
		throw std::runtime_error("native tetra partition: owner gather failed");
	for(const int value:owner)if(value<0||value>=ranks)
		throw std::runtime_error("native tetra partition: ParMETIS returned an invalid rank");
	return owner;
}

} // namespace native_tet_partition_detail

inline NativeTetMeshPartition PartitionNativeTetMesh(const NativeTetMesh& mesh,MPI_Comm communicator)
{
	using namespace native_tet_partition_detail;
	NativeTetMeshPartition result;
	MPI_Comm_rank(communicator,&result.rank);MPI_Comm_size(communicator,&result.ranks);
	if(mesh.points.empty()||mesh.cells.empty())
		throw std::invalid_argument("native tetra partition needs a nonempty mesh");
	const std::size_t nodes=mesh.points.size(),cells=mesh.cells.size();
	result.cell_owner=PartitionCells(mesh,communicator);
	// A vertex goes to the rank with most adjacent cells; ties take the lower rank.
	result.node_owner.assign(nodes,-1);
	{
		std::vector<std::map<int,int>> votes(nodes);
		for(std::size_t cell=0;cell<cells;++cell)
			for(const auto node:mesh.cells[cell].nodes)++votes[node][result.cell_owner[cell]];
		for(std::size_t node=0;node<nodes;++node){
			if(votes[node].empty())throw std::invalid_argument("native tetra partition: isolated vertex");
			int best=-1,best_count=0;
			for(const auto& vote:votes[node])if(vote.second>best_count){best=vote.first;best_count=vote.second;}
			result.node_owner[node]=best;
		}
	}
	// Owner-major global vertex numbering: rank blocks in rank order, global id order inside.
	std::vector<std::uint32_t> counts(static_cast<std::size_t>(result.ranks),0);
	for(const int owner:result.node_owner)++counts[static_cast<std::size_t>(owner)];
	result.owned_node_begin.assign(static_cast<std::size_t>(result.ranks),0);
	for(int r=1;r<result.ranks;++r)result.owned_node_begin[r]=result.owned_node_begin[r-1]+counts[r-1];
	result.owner_major_node_.resize(nodes);
	{
		std::vector<std::uint32_t> next=result.owned_node_begin;
		for(std::size_t node=0;node<nodes;++node)
			result.owner_major_node_[node]=next[static_cast<std::size_t>(result.node_owner[node])]++;
	}
	// Local cells: owned, plus halo cells touching an owned vertex.
	std::vector<std::uint8_t> keep(cells,0);
	for(std::size_t cell=0;cell<cells;++cell){
		if(result.cell_owner[cell]==result.rank)keep[cell]=1;
		else for(const auto node:mesh.cells[cell].nodes)
			if(result.node_owner[node]==result.rank){keep[cell]=1;break;}
	}
	std::vector<std::uint32_t> global_to_local(nodes,std::numeric_limits<std::uint32_t>::max());
	std::vector<std::uint32_t> owned,ghost;
	for(std::size_t node=0;node<nodes;++node)
		if(result.node_owner[node]==result.rank)owned.push_back(static_cast<std::uint32_t>(node));
	std::vector<std::uint8_t> referenced(nodes,0);
	for(std::size_t cell=0;cell<cells;++cell)if(keep[cell])
		for(const auto node:mesh.cells[cell].nodes)referenced[node]=1;
	for(std::size_t node=0;node<nodes;++node)
		if(referenced[node]&&result.node_owner[node]!=result.rank)ghost.push_back(static_cast<std::uint32_t>(node));
	result.owned_node_count=static_cast<std::uint32_t>(owned.size());
	result.local_to_global_node=owned;
	result.local_to_global_node.insert(result.local_to_global_node.end(),ghost.begin(),ghost.end());
	for(std::size_t local=0;local<result.local_to_global_node.size();++local)
		global_to_local[result.local_to_global_node[local]]=static_cast<std::uint32_t>(local);
	result.local.points.reserve(result.local_to_global_node.size());
	for(const auto node:result.local_to_global_node)result.local.points.push_back(mesh.points[node]);
	for(std::size_t cell=0;cell<cells;++cell)if(keep[cell]){
		NativeTetCell local_cell=mesh.cells[cell];
		for(auto& node:local_cell.nodes)node=global_to_local[node];
		result.local.cells.push_back(local_cell);
		result.local_to_global_cell.push_back(static_cast<std::uint32_t>(cell));
		result.local_cell_owned.push_back(result.cell_owner[cell]==result.rank?1:0);
	}
	// Boundary triangles whose three vertices are local; owned when their cell is.
	std::map<std::array<std::uint32_t,3>,std::size_t> face_cell;
	for(std::size_t local_cell=0;local_cell<result.local.cells.size();++local_cell){
		const auto& cell_nodes=result.local.cells[local_cell].nodes;
		for(int opposite=0;opposite<4;++opposite){
			std::array<std::uint32_t,3> face{};int entry=0;
			for(int vertex=0;vertex<4;++vertex)if(vertex!=opposite)face[entry++]=cell_nodes[vertex];
			std::sort(face.begin(),face.end());
			face_cell.emplace(face,local_cell);
		}
	}
	for(const auto& triangle:mesh.boundary_triangles){
		std::array<std::uint32_t,3> face{};bool local=true;
		for(int vertex=0;vertex<3;++vertex){
			const auto mapped=global_to_local[triangle.nodes[vertex]];
			if(mapped==std::numeric_limits<std::uint32_t>::max()){local=false;break;}
			face[vertex]=mapped;
		}
		if(!local)continue;
		auto sorted=face;std::sort(sorted.begin(),sorted.end());
		const auto found=face_cell.find(sorted);
		if(found==face_cell.end())continue;
		auto local_triangle=triangle;
		local_triangle.nodes=face;
		result.local.boundary_triangles.push_back(local_triangle);
		result.local_triangle_owned.push_back(result.local_cell_owned[found->second]);
	}
	return result;
}

} // namespace iga

#endif

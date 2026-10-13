#ifndef IGA_NATIVE_TET_DISTRIBUTED_DOFS_HPP
#define IGA_NATIVE_TET_DISTRIBUTED_DOFS_HPP

// Global PETSc row layout of the Taylor-Hood unknowns on a partitioned mesh.
// Rank r owns one contiguous block: velocity rows of its owned vertices
// (three per vertex, global-id order), then of its owned edges (three per
// edge, (min, max) global-endpoint order), then pressure rows of its owned
// vertices. An edge belongs to the rank that owns its lower-global-id
// endpoint, so every cell that touches an owned row is local. Ghost rows of
// the other local vertices and edges are resolved from the replicated owner
// arrays (vertices) or by asking the owning rank (edges), and listed in the
// order a ghosted PETSc vector uses.

#include "NativeTetFem.hpp"
#include "NativeTetMeshPartition.hpp"

#include <petscsys.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct NativeTetDistributedDofs
{
	NativeTaylorHoodTopology topology;        // of partition.local (local vertex ids)
	std::uint32_t owned_vertices=0,owned_edges=0;
	std::vector<std::uint8_t> edge_owned;     // local edge -> owned by this rank
	PetscInt row_begin=0,row_end=0,global_rows=0;
	PetscInt velocity_rows_end=0;             // owned velocity rows are [row_begin, velocity_rows_end)
	std::vector<PetscInt> vertex_velocity_row;  // local vertex -> global row of component 0
	std::vector<PetscInt> edge_velocity_row;    // local edge -> global row of component 0
	std::vector<PetscInt> vertex_pressure_row;  // local vertex -> global pressure row
	std::vector<PetscInt> ghost_rows;           // non-owned local rows, ghosted-vector order
	std::vector<PetscInt> vertex_velocity_local,edge_velocity_local,vertex_pressure_local;  // index into the ghosted local array (component 0)
	std::vector<std::array<std::uint32_t,2>> edge_global;  // local edge -> (min, max) global vertex ids

	PetscInt OwnedRows()const{return row_end-row_begin;}
	PetscInt LocalSize()const{return OwnedRows()+static_cast<PetscInt>(ghost_rows.size());}
	// Global row of the 34 unknowns of one local cell in element order.
	std::array<PetscInt,34> CellRows(std::size_t cell)const
	{
		std::array<PetscInt,34> rows{};
		std::size_t entry=0;
		const std::size_t vertices=vertex_velocity_row.size();
		for(const auto node:topology.cell_velocity_nodes.at(cell)){
			const PetscInt base=node<vertices?vertex_velocity_row[node]:edge_velocity_row[node-vertices];
			for(int component=0;component<3;++component)rows[entry++]=base+component;
		}
		for(const auto node:local_cells_->at(cell).nodes)rows[entry++]=vertex_pressure_row[node];
		return rows;
	}
	// Ghosted-local index of the same 34 unknowns.
	std::array<PetscInt,34> CellLocal(std::size_t cell)const
	{
		std::array<PetscInt,34> rows{};
		std::size_t entry=0;
		const std::size_t vertices=vertex_velocity_local.size();
		for(const auto node:topology.cell_velocity_nodes.at(cell)){
			const PetscInt base=node<vertices?vertex_velocity_local[node]:edge_velocity_local[node-vertices];
			for(int component=0;component<3;++component)rows[entry++]=base+component;
		}
		for(const auto node:local_cells_->at(cell).nodes)rows[entry++]=vertex_pressure_local[node];
		return rows;
	}

	const std::vector<NativeTetCell>* local_cells_=nullptr;
};

namespace native_tet_distributed_detail {

inline void CheckMpi(int code,const char* context)
{
	if(code!=MPI_SUCCESS)throw std::runtime_error(std::string("native tetra distributed dofs: ")+context+" failed");
}

} // namespace native_tet_distributed_detail

inline NativeTetDistributedDofs BuildNativeTetDistributedDofs(const NativeTetMeshPartition& partition,
	MPI_Comm communicator)
{
	using namespace native_tet_distributed_detail;
	NativeTetDistributedDofs result;
	result.local_cells_=&partition.local.cells;
	result.topology=BuildNativeTaylorHoodTopology(partition.local);
	const int rank=partition.rank,ranks=partition.ranks;
	const std::size_t vertices=partition.local.points.size(),edges=result.topology.edges.size();
	result.owned_vertices=partition.owned_node_count;
	// Edge ownership from the lower-global-id endpoint.
	result.edge_global.resize(edges);
	result.edge_owned.assign(edges,0);
	std::vector<int> edge_owner(edges);
	for(std::size_t edge=0;edge<edges;++edge){
		const auto& ends=result.topology.edges[edge];
		const std::uint32_t a=partition.local_to_global_node[ends[0]],b=partition.local_to_global_node[ends[1]];
		result.edge_global[edge]={{std::min(a,b),std::max(a,b)}};
		edge_owner[edge]=partition.node_owner[result.edge_global[edge][0]];
		result.edge_owned[edge]=edge_owner[edge]==rank?1:0;
	}
	std::vector<std::uint32_t> owned_edge_list;
	for(std::size_t edge=0;edge<edges;++edge)if(result.edge_owned[edge])owned_edge_list.push_back(static_cast<std::uint32_t>(edge));
	std::sort(owned_edge_list.begin(),owned_edge_list.end(),[&](std::uint32_t x,std::uint32_t y){
		return result.edge_global[x]<result.edge_global[y];});
	result.owned_edges=static_cast<std::uint32_t>(owned_edge_list.size());
	// Row blocks per rank.
	const PetscInt local_rows=static_cast<PetscInt>(4u*result.owned_vertices+3u*result.owned_edges);
	std::vector<PetscInt> rows_of(ranks),vertices_of(ranks),edges_of(ranks);
	const PetscInt owned_vertices_int=static_cast<PetscInt>(result.owned_vertices),owned_edges_int=static_cast<PetscInt>(result.owned_edges);
	CheckMpi(MPI_Allgather(&local_rows,1,MPIU_INT,rows_of.data(),1,MPIU_INT,communicator),"row count gather");
	CheckMpi(MPI_Allgather(&owned_vertices_int,1,MPIU_INT,vertices_of.data(),1,MPIU_INT,communicator),"vertex count gather");
	CheckMpi(MPI_Allgather(&owned_edges_int,1,MPIU_INT,edges_of.data(),1,MPIU_INT,communicator),"edge count gather");
	std::vector<PetscInt> row_begin_of(ranks,0);
	for(int r=1;r<ranks;++r)row_begin_of[r]=row_begin_of[r-1]+rows_of[r-1];
	result.row_begin=row_begin_of[rank];result.row_end=result.row_begin+local_rows;
	result.global_rows=row_begin_of[ranks-1]+rows_of[ranks-1];
	result.velocity_rows_end=result.row_begin+3*(owned_vertices_int+owned_edges_int);
	// Vertex rows: owned from the local order, ghosts from the owner's block.
	result.vertex_velocity_row.assign(vertices,-1);
	result.vertex_pressure_row.assign(vertices,-1);
	for(std::size_t vertex=0;vertex<vertices;++vertex){
		const std::uint32_t global=partition.local_to_global_node[vertex];
		const int owner=partition.node_owner[global];
		const PetscInt index=static_cast<PetscInt>(partition.OwnerMajorNode(global)-partition.owned_node_begin[owner]);
		result.vertex_velocity_row[vertex]=row_begin_of[owner]+3*index;
		result.vertex_pressure_row[vertex]=row_begin_of[owner]+3*(vertices_of[owner]+edges_of[owner])+index;
	}
	// Edge rows: owned by position; others by asking their owner.
	result.edge_velocity_row.assign(edges,-1);
	std::map<std::array<std::uint32_t,2>,PetscInt> owned_edge_row;
	for(std::size_t j=0;j<owned_edge_list.size();++j){
		const auto edge=owned_edge_list[j];
		result.edge_velocity_row[edge]=result.row_begin+3*(owned_vertices_int+static_cast<PetscInt>(j));
		owned_edge_row.emplace(result.edge_global[edge],result.edge_velocity_row[edge]);
	}
	std::vector<std::vector<std::uint32_t>> queries(ranks);
	std::vector<std::vector<std::uint32_t>> queried_edges(ranks);
	for(std::size_t edge=0;edge<edges;++edge)if(!result.edge_owned[edge]){
		queries[edge_owner[edge]].push_back(result.edge_global[edge][0]);
		queries[edge_owner[edge]].push_back(result.edge_global[edge][1]);
		queried_edges[edge_owner[edge]].push_back(static_cast<std::uint32_t>(edge));
	}
	std::vector<int> send_counts(ranks),receive_counts(ranks),send_offsets(ranks),receive_offsets(ranks);
	for(int r=0;r<ranks;++r)send_counts[r]=static_cast<int>(queries[r].size());
	CheckMpi(MPI_Alltoall(send_counts.data(),1,MPI_INT,receive_counts.data(),1,MPI_INT,communicator),"query counts");
	for(int r=1;r<ranks;++r){send_offsets[r]=send_offsets[r-1]+send_counts[r-1];receive_offsets[r]=receive_offsets[r-1]+receive_counts[r-1];}
	std::vector<std::uint32_t> send_buffer(static_cast<std::size_t>(send_offsets[ranks-1]+send_counts[ranks-1]));
	for(int r=0;r<ranks;++r)std::copy(queries[r].begin(),queries[r].end(),send_buffer.begin()+send_offsets[r]);
	std::vector<std::uint32_t> receive_buffer(static_cast<std::size_t>(receive_offsets[ranks-1]+receive_counts[ranks-1]));
	CheckMpi(MPI_Alltoallv(send_buffer.data(),send_counts.data(),send_offsets.data(),MPI_UINT32_T,
		receive_buffer.data(),receive_counts.data(),receive_offsets.data(),MPI_UINT32_T,communicator),"edge queries");
	// Answer: one row per queried pair, in the order received.
	std::vector<PetscInt> answers(receive_buffer.size()/2);
	for(std::size_t pair=0;pair<answers.size();++pair){
		const std::array<std::uint32_t,2> key{{receive_buffer[2*pair],receive_buffer[2*pair+1]}};
		const auto found=owned_edge_row.find(key);
		if(found==owned_edge_row.end())
			throw std::runtime_error("native tetra distributed dofs: a queried edge is not owned by the asked rank");
		answers[pair]=found->second;
	}
	std::vector<int> answer_send(ranks),answer_receive(ranks),answer_send_offsets(ranks),answer_receive_offsets(ranks);
	for(int r=0;r<ranks;++r){answer_send[r]=receive_counts[r]/2;answer_receive[r]=send_counts[r]/2;}
	for(int r=1;r<ranks;++r){answer_send_offsets[r]=answer_send_offsets[r-1]+answer_send[r-1];
		answer_receive_offsets[r]=answer_receive_offsets[r-1]+answer_receive[r-1];}
	std::vector<PetscInt> received_rows(send_buffer.size()/2);
	CheckMpi(MPI_Alltoallv(answers.data(),answer_send.data(),answer_send_offsets.data(),MPIU_INT,
		received_rows.data(),answer_receive.data(),answer_receive_offsets.data(),MPIU_INT,communicator),"edge answers");
	for(int r=0;r<ranks;++r)
		for(std::size_t k=0;k<queried_edges[r].size();++k)
			result.edge_velocity_row[queried_edges[r][k]]=received_rows[static_cast<std::size_t>(answer_receive_offsets[r])+k];
	// Ghosted-local indices: owned rows map to [0, local_rows); ghosts follow.
	result.vertex_velocity_local.assign(vertices,-1);
	result.vertex_pressure_local.assign(vertices,-1);
	result.edge_velocity_local.assign(edges,-1);
	PetscInt next=local_rows;
	const auto place=[&](PetscInt global_row,int count,PetscInt& local_index){
		if(global_row>=result.row_begin&&global_row<result.row_end){local_index=global_row-result.row_begin;return;}
		local_index=next;
		for(int component=0;component<count;++component){result.ghost_rows.push_back(global_row+component);++next;}
	};
	for(std::size_t vertex=0;vertex<vertices;++vertex){
		place(result.vertex_velocity_row[vertex],3,result.vertex_velocity_local[vertex]);
	}
	for(std::size_t edge=0;edge<edges;++edge)place(result.edge_velocity_row[edge],3,result.edge_velocity_local[edge]);
	for(std::size_t vertex=0;vertex<vertices;++vertex)place(result.vertex_pressure_row[vertex],1,result.vertex_pressure_local[vertex]);
	return result;
}

} // namespace iga

#endif

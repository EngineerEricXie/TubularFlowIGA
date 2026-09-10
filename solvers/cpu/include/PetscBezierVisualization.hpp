#ifndef IGA_PETSC_BEZIER_VISUALIZATION_HPP
#define IGA_PETSC_BEZIER_VISUALIZATION_HPP

#include "ParallelBezierVisualization.hpp"
#include "PetscGather.hpp"

namespace iga {
struct PetscBezierPartition {
	VtkPartition piece;
	std::size_t requested_nodes=0,requested_rows=0;
	std::uint64_t global_rows=0;
};
namespace petsc_bezier_detail {
struct Indices {
	IS source=nullptr,target=nullptr;
	~Indices() { if(source)ISDestroy(&source);if(target)ISDestroy(&target); }
	void Close(MPI_Comm comm)
	{
		RequireCollectivePetscSuccess(comm,"Bezier source index cleanup",ISDestroy(&source));
		RequireCollectivePetscSuccess(comm,"Bezier target index cleanup",ISDestroy(&target));
	}
};
}

// Reads only rows referenced by this rank's owned elements. The source vector
// is borrowed, interleaved by node, and is never modified. All source ranks
// enter together, including ranks with no owned visualization elements.
inline PetscBezierPartition BuildPetscBezierPartition(Vec state,const std::vector<Element>& owned_elements,
	std::uint64_t global_elements,std::uint64_t global_nodes,const std::vector<VtkArraySchema>& fields,
	PointIdentityLimits limits={})
{
	const auto comm=PetscObjectComm(reinterpret_cast<PetscObject>(state));
	std::vector<std::int32_t> nodes;std::vector<PetscInt> rows;
	int components=0;PetscBezierPartition result;
	CollectiveLocalStage(comm,"PETSc Bezier extraction preparation",[&] {
		partitioned_vtk_detail::ValidateSchema(fields);
		for(const auto& field:fields) {
			if(field.components>std::numeric_limits<int>::max()-components)throw std::overflow_error("Bezier state component count overflows");
			components+=field.components;
		}
		if(!components||global_nodes>static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max())/static_cast<unsigned>(components))
			throw std::invalid_argument("invalid Bezier state dimensions");
		PetscInt size=0;if(VecGetSize(state,&size))throw std::runtime_error("cannot query Bezier source vector");
		result.global_rows=global_nodes*static_cast<unsigned>(components);
		if(size<0||static_cast<std::uint64_t>(size)!=result.global_rows)throw std::invalid_argument("Bezier state size differs from control mesh");
		std::set<std::int32_t> required;
		for(const auto& element:owned_elements)for(auto node:element.connectivity) {
			if(node<0||static_cast<std::uint64_t>(node)>=global_nodes)throw std::invalid_argument("Bezier control node is out of range");
			required.insert(node);
		}
		nodes.assign(required.begin(),required.end());
		if(nodes.size()>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())/static_cast<unsigned>(components))
			throw std::overflow_error("Bezier selected row count exceeds PETSc range");
		result.requested_nodes=nodes.size();result.requested_rows=nodes.size()*static_cast<unsigned>(components);
		rows.reserve(result.requested_rows);
		for(auto node:nodes)for(int component=0;component<components;++component)rows.push_back(static_cast<PetscInt>(node)*components+component);
	});
	RequireCollectiveSameInt(comm,"PETSc Bezier component agreement",components);
	petsc_bezier_detail::Indices indices;
	PetscGatherObjects objects;
	RequireCollectivePetscSuccess(comm,"Bezier selected vector create",VecCreateSeq(PETSC_COMM_SELF,static_cast<PetscInt>(rows.size()),&objects.all));
	RequireCollectivePetscSuccess(comm,"Bezier source index create",ISCreateGeneral(PETSC_COMM_SELF,static_cast<PetscInt>(rows.size()),rows.data(),PETSC_COPY_VALUES,&indices.source));
	RequireCollectivePetscSuccess(comm,"Bezier target index create",ISCreateStride(PETSC_COMM_SELF,static_cast<PetscInt>(rows.size()),0,1,&indices.target));
	RequireCollectivePetscSuccess(comm,"Bezier selected scatter create",VecScatterCreate(state,indices.source,objects.all,indices.target,&objects.scatter));
	RequireCollectivePetscSuccess(comm,"Bezier selected scatter begin",VecScatterBegin(objects.scatter,state,objects.all,INSERT_VALUES,SCATTER_FORWARD));
	RequireCollectivePetscSuccess(comm,"Bezier selected scatter end",VecScatterEnd(objects.scatter,state,objects.all,INSERT_VALUES,SCATTER_FORWARD));
	PetscReadArray view;
	CollectiveLocalStage(comm,"Bezier selected state read",[&] { view.Acquire(objects.all); });
	result.piece=BuildParallelBezierPartition(comm,owned_elements,global_elements,fields,[&](std::int32_t node,int component) {
		const auto found=std::lower_bound(nodes.begin(),nodes.end(),node);
		if(found==nodes.end()||*found!=node)throw std::runtime_error("missing selected Bezier control node");
		return PetscRealPart(view.Data()[static_cast<std::size_t>(found-nodes.begin())*static_cast<unsigned>(components)+component]);
	},limits);
	CollectiveLocalStage(comm,"Bezier selected state restore",[&] { view.Restore(); });
	objects.Close(comm);indices.Close(comm);
	return result;
}
} // namespace iga
#endif

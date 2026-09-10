#ifndef IGA_DISTRIBUTED_MATERIAL_SURFACE_PATCH_KINEMATICS_HPP
#define IGA_DISTRIBUTED_MATERIAL_SURFACE_PATCH_KINEMATICS_HPP

#include "MaterialSurfacePatchKinematics.hpp"
#include "SurfaceGhostKinematics.hpp"
#include <optional>

namespace iga {

// Owned publications stay partitioned. The closed geometry owner currently
// needs the complete patch coordinates on every rank; this bounded fetch does
// not replicate any fluid solution or perform a numerical solve.
inline MaterialSurfacePatchKinematics::Result ComposeDistributedMaterialSurfacePatch(
	MPI_Comm comm,const MaterialSurfacePatchMap& map,const MaterialSurfaceKinematics& committed,
	const SurfaceKinematics& publication,const SurfaceFieldStamp& expected_stamp,
	const std::string& expected_committed_identity,const FsiTrialContext& context,
	std::size_t maximum_patch_nodes=4096,std::size_t maximum_material_vertices=1000000,
	PointIdentityLimits limits={})
{
	std::string common,local_publication,contributors;
	std::vector<std::uint64_t> requested;
	int ranks=0;MPI_Comm_size(comm,&ranks);
	CollectiveLocalStage(comm,"distributed material patch input",[&] {
		ValidateFsiTrialContext(context);committed.Validate();
		if(committed.ContentIdentitySha256()!=expected_committed_identity)
			throw std::invalid_argument("material patch committed epoch differs from authority");
		if(expected_stamp.time_s!=context.EndTime() || expected_stamp.step!=context.step
			|| expected_stamp.coupling_iteration!=context.coupling_iteration)
			throw std::invalid_argument("material patch expected stamp differs from context");
		if(!maximum_patch_nodes || !maximum_material_vertices
			|| map.GlobalToSourceVertices().size()>maximum_patch_nodes
			|| committed.SourceVerticesM().size()>maximum_material_vertices)
			throw std::invalid_argument("material patch geometry exceeds configured limits");
		for(const auto& vertex:map.GlobalToSourceVertices())requested.push_back(vertex.first);
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedMaterialSurfacePatch/input/v1");
		for(const auto* value:{&expected_committed_identity,&map.ReferenceIdentitySha256(),&map.FullReference().ContentIdentitySha256()})
			distributed_surface_detail::AppendString(hash,*value);
		hash.AppendLittleEndian64(maximum_patch_nodes);hash.AppendLittleEndian64(maximum_material_vertices);
		hash.AppendLittleEndian64(map.ConfiguredClampedGlobalNodeIds().size());
		for(auto id:map.ConfiguredClampedGlobalNodeIds())hash.AppendLittleEndian64(id);
		hash.AppendLittleEndian64(context.step);hash.AppendNormalizedDouble(context.start_time_s);
		hash.AppendNormalizedDouble(context.dt_s);hash.AppendLittleEndian64(context.coupling_iteration);
		common=hash.Hex();
		local_publication=BuildSurfaceKinematicsIdentitySha256(publication,map.Layout());
		if(static_cast<std::size_t>(ranks)>contributors.max_size()/64)
			throw std::overflow_error("material patch contributor metadata overflows");
		contributors.resize(64*static_cast<std::size_t>(ranks));
	});
	RequireCollectiveSameText(comm,"distributed material patch reference agreement",common);
	const auto fields=FetchSurfaceGhostKinematics(comm,map.Layout(),publication,map.Interface().id,expected_stamp,requested,limits);
	MPI_Allgather(local_publication.data(),64,MPI_CHAR,contributors.data(),64,MPI_CHAR,comm);
	std::optional<MaterialSurfacePatchKinematics::Result> result;
	CollectiveLocalStage(comm,"distributed material patch geometry",[&] {
		auto target=MaterialSurfacePatchKinematics::ComposeCompletePatch(map,committed,context,
			fields.requested_node_ids,fields.displacement_m,fields.velocity_m_per_s);
		Sha256 hash;distributed_surface_detail::AppendString(hash,"DistributedMaterialSurfacePatch/composition/v1");
		for(const auto* value:std::array<const std::string*,3>{{&common,&contributors,&target.ContentIdentitySha256()}})
			distributed_surface_detail::AppendString(hash,*value);
		result.emplace(MaterialSurfacePatchKinematics::Result{std::move(target),hash.Hex()});
	});
	RequireCollectiveSameText(comm,"distributed material patch target agreement",result->target.ContentIdentitySha256());
	return std::move(*result);
}

} // namespace iga
#endif

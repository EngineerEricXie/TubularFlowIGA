#ifndef IGA_MOVING_FSI_CHECKPOINT_BUNDLE_HPP
#define IGA_MOVING_FSI_CHECKPOINT_BUNDLE_HPP

#include "MovingCheckpointRestore.hpp"
#include "ImmersedMovingDistributedFsiRuntime.hpp"
#include "SingleOwnerMembraneRuntime.hpp"

namespace iga {

inline std::vector<CoupledCheckpointShardSpec> MovingFsiCheckpointCatalog(
	const MovingCheckpointLayout& layout,std::uint32_t world_size,std::uint32_t membrane_owner)
{
	auto catalog=MovingCheckpointCatalog(layout,world_size);
	coupled_checkpoint_detail::Require(std::find(layout.world_ranks.begin(),layout.world_ranks.end(),membrane_owner)!=layout.world_ranks.end(),
		"membrane checkpoint owner outside fluid communicator");
	const auto prefix=layout.domain_id+".";
	catalog.push_back({prefix+"membrane.metadata","single-owner-membrane-metadata-v1",membrane_owner});
	catalog.push_back({prefix+"membrane.numerical","membrane-state-v1",membrane_owner});
	for(auto rank:layout.world_ranks)catalog.push_back({prefix+"fsi.rank-"+std::to_string(rank),"moving-fsi-publication-v1",rank});
	std::sort(catalog.begin(),catalog.end(),[](const auto& a,const auto& b) { return a.id<b.id; });
	coupled_checkpoint_detail::ValidateCatalog(catalog,world_size);return catalog;
}

namespace moving_fsi_checkpoint_detail {
inline void RequireContext(const FsiTrialContext& fluid,const FsiTrialContext& membrane,const CoupledCheckpointEpoch& epoch)
{
	ValidateFsiTrialContext(fluid);ValidateFsiTrialContext(membrane);
	checkpoint_metadata::Require(fluid.step==membrane.step&&fluid.start_time_s==membrane.start_time_s
		&&fluid.dt_s==membrane.dt_s&&fluid.coupling_iteration==membrane.coupling_iteration
		&&fluid.step==epoch.accepted_steps&&fluid.EndTime()==epoch.time_s&&fluid.dt_s==epoch.dt_s,
		"paired FSI checkpoint context differs");
}
inline std::string Hash(std::string_view bytes)
{
	Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex();
}
inline std::map<std::string,std::size_t> Indices(const CoupledCheckpointManifest& manifest,
	const MovingCheckpointLayout& layout,std::uint32_t owner)
{
	std::map<std::string,std::size_t> result;
	for(const auto& spec:MovingFsiCheckpointCatalog(layout,manifest.epoch.compatibility.ranks,owner)) {
		const auto found=std::find_if(manifest.shards.begin(),manifest.shards.end(),[&](const auto& value) { return value.spec.id==spec.id; });
		checkpoint_metadata::Require(found!=manifest.shards.end()&&coupled_checkpoint_detail::Same(found->spec,spec),"incomplete paired FSI catalog");
		result.emplace(spec.id.substr(layout.domain_id.size()+1),static_cast<std::size_t>(found-manifest.shards.begin()));
	}
	return result;
}
inline MovingFsiPublicationCheckpoint RepartitionPublication(const std::filesystem::path& root,
	const CoupledCheckpointManifest& manifest,const MovingCheckpointLayout& source_layout,
	const std::map<std::string,std::size_t>& indices,const DistributedSurfaceLayout& target,
	const SurfaceInterfaceRef& interface,const std::string& configuration,const std::string& material_identity,
	const std::string& geometry_identity,std::size_t maximum_nodes)
{
	using checkpoint_metadata::Require;
	ValidateDistributedSurfaceLayout(target);Require(target.global_node_count<=maximum_nodes,"FSI target surface exceeds bound");
	std::map<std::uint64_t,std::pair<std::array<double,3>,std::array<double,3>>> records;
	std::vector<std::string> payload_identities;std::set<std::string> partitions;
	std::optional<MovingFsiPublicationRecords> authority;
	for(auto source_rank:source_layout.world_ranks) {
		const auto index=indices.at("fsi.rank-"+std::to_string(source_rank));
		const auto bytes=moving_checkpoint_bundle_detail::Metadata(root,manifest,index);
		const auto identity=Hash(bytes);payload_identities.push_back(identity);
		auto shard=DecodeMovingFsiPublicationRecords(bytes,identity,interface,configuration,material_identity,
			geometry_identity,manifest.epoch.accepted_steps,manifest.epoch.time_s,maximum_nodes);
		Require(shard.reference_identity==target.reference_mesh_identity_sha256,"FSI source reference mesh differs");
		Require(partitions.insert(shard.partition_identity).second,"FSI source partition identity repeats");
		if(!authority)authority.emplace(shard);
		else Require(shard.layout_identity==authority->layout_identity
			&&shard.state.context.step==authority->state.context.step
			&&shard.state.context.start_time_s==authority->state.context.start_time_s
			&&shard.state.context.dt_s==authority->state.context.dt_s
			&&shard.state.context.coupling_iteration==authority->state.context.coupling_iteration
			&&shard.state.composition_identity==authority->state.composition_identity,"FSI source slice epoch differs");
		for(std::size_t row=0;row<shard.node_ids.size();++row) {
			Require(records.emplace(shard.node_ids[row],std::make_pair(shard.state.traction.traction_on_structure_pa[row],
				shard.state.traction.consistent_nodal_force_n[row])).second,"FSI source surface ownership overlaps");
			Require(records.size()<=target.global_node_count,"FSI source surface coverage exceeds target");
		}
	}
	Require(authority.has_value()&&records.size()==target.global_node_count,"FSI source surface coverage differs");
	for(const auto& reference:target.reference_positions)Require(records.count(reference.global_node_id)==1,"FSI source surface node missing");
	MovingFsiPublicationCheckpoint result;result.context=authority->state.context;
	result.composition_identity=authority->state.composition_identity;result.material_identity=material_identity;result.geometry_identity=geometry_identity;
	auto& traction=result.traction;traction.interface=interface;traction.stamp.time_s=result.context.EndTime();
	traction.stamp.step=result.context.step;traction.stamp.coupling_iteration=result.context.coupling_iteration;
	traction.stamp.reference_mesh_identity_sha256=target.reference_mesh_identity_sha256;
	traction.stamp.layout_identity_sha256=target.layout_identity_sha256;
	traction.stamp.partition_identity_sha256=BuildDistributedSurfacePartitionIdentitySha256(target);
	Sha256 producer;distributed_surface_detail::AppendString(producer,"MovingFsiCheckpoint/repartitioned-producer/v1");
	for(const auto& identity:payload_identities)distributed_surface_detail::AppendString(producer,identity);
	distributed_surface_detail::AppendString(producer,traction.stamp.partition_identity_sha256);
	traction.stamp.producer_state_identity_sha256=producer.Hex();
	for(auto id:target.owned_global_node_ids) {
		traction.traction_on_structure_pa.push_back(records.at(id).first);
		traction.consistent_nodal_force_n.push_back(records.at(id).second);
	}
	Sha256 projection;distributed_surface_detail::AppendString(projection,"MovingFsiCheckpoint/repartitioned-projection/v1");
	distributed_surface_detail::AppendString(projection,traction.stamp.producer_state_identity_sha256);
	for(const auto* values:{&traction.traction_on_structure_pa,&traction.consistent_nodal_force_n})
		for(const auto& tuple:*values)for(double value:tuple)projection.AppendNormalizedDouble(value);
	traction.projection_identity_sha256=projection.Hex();ValidateSurfaceTraction(traction,target);return result;
}
}

// Collective capture and local shard writes; the caller must gather receipts
// and publish the complete paired catalog, never the fluid-only subset.
inline std::vector<CoupledCheckpointShard> WriteMovingFsiCheckpointShards(MPI_Comm comm,
	const std::filesystem::path& root,const CoupledCheckpointEpoch& epoch,const MovingCheckpointLayout& layout,
	std::uint32_t membrane_owner,const std::string& configuration,
	ImmersedMovingTransientDistributedRuntime& flow,ImmersedMovingDistributedFsiRuntime& fluid,
	SingleOwnerMembraneRuntime& membrane)
{
	int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
	std::string agreement;
	CollectiveLocalStage(comm,"paired checkpoint capture preflight",[&] {
		for(auto runtime_comm:{flow.Communicator(),fluid.Communicator(),membrane.Communicator()}) {
			int comparison=MPI_UNEQUAL;MPI_Comm_compare(comm,runtime_comm,&comparison);
			checkpoint_metadata::Require(comparison==MPI_IDENT||comparison==MPI_CONGRUENT,"paired checkpoint communicator differs");
		}
		checkpoint_metadata::Require(layout.world_ranks.size()==static_cast<std::size_t>(size),"paired checkpoint mapping differs");
		checkpoint_metadata::Require(fluid.UsesFlowRuntime(flow),"paired checkpoint adapter belongs to another fluid owner");
		coupled_checkpoint_detail::Validate(epoch);
		(void)MovingFsiCheckpointCatalog(layout,epoch.compatibility.ranks,membrane_owner);
		moving_fsi_checkpoint_detail::RequireContext(fluid.CommittedContext(),membrane.CommittedContext(),epoch);
		Sha256 hash;
		for(const auto& text:{epoch.id,epoch.previous_id,epoch.compatibility.case_sha256,epoch.compatibility.configuration_sha256,
			epoch.compatibility.execution_sha256,configuration,root.generic_string(),layout.domain_id})distributed_surface_detail::AppendString(hash,text);
		hash.AppendLittleEndian64(epoch.compatibility.ranks);hash.AppendLittleEndian64(epoch.accepted_steps);
		hash.AppendNormalizedDouble(epoch.time_s);hash.AppendNormalizedDouble(epoch.dt_s);
		hash.AppendLittleEndian64(membrane_owner);hash.AppendLittleEndian64(fluid.CommittedContext().coupling_iteration);
		for(auto member:layout.world_ranks)hash.AppendLittleEndian64(member);
		agreement=hash.Hex();
	});
	RequireCollectiveSameText(comm,"paired checkpoint capture agreement",agreement);
	const auto state=flow.CaptureAcceptedCheckpoint();const auto publication=fluid.CaptureCheckpoint();
	const auto structure=membrane.CaptureCheckpoint();std::vector<CoupledCheckpointShard> receipts;
	CollectiveLocalStage(comm,"paired checkpoint shard writes",[&] {
		const auto owner=layout.world_ranks.at(rank);
		checkpoint_metadata::Require(structure.has_value()==(owner==membrane_owner),"paired checkpoint membrane owner differs");
		receipts=WriteMovingCheckpointShards(root,epoch,layout,rank,state,configuration);
		for(const auto& spec:MovingFsiCheckpointCatalog(layout,epoch.compatibility.ranks,membrane_owner)) {
			if(spec.owner_rank!=owner)continue;
			const std::string* bytes=nullptr;
			if(spec.id==layout.domain_id+".fsi.rank-"+std::to_string(owner))bytes=&publication;
			else if(spec.id==layout.domain_id+".membrane.metadata")bytes=&structure->metadata;
			else if(spec.id==layout.domain_id+".membrane.numerical")bytes=&structure->numerical;
			if(bytes)receipts.push_back(WriteCoupledCheckpointShard(root,epoch,spec,bytes->size(),[&](auto& out) { out.Write(bytes->data(),bytes->size()); }));
		}
		std::sort(receipts.begin(),receipts.end(),[](const auto& a,const auto& b) { return a.spec.id<b.spec.id; });
	});
	return receipts;
}

// Member order preserves borrowed lifetimes. The map and communicator remain
// caller-owned. Close is collective; destruction alone does not close PETSc.
struct RestoredMovingFsiPair {
	std::unique_ptr<ImmersedMovingTransientDistributedRuntime> flow;
	std::unique_ptr<ImmersedMovingDistributedFsiRuntime> fluid;
	std::unique_ptr<SingleOwnerMembraneRuntime> membrane;
	void Close() { if(flow)flow->Close(); }
};

inline std::unique_ptr<RestoredMovingFsiPair> RestoreMovingFsiCheckpoint(MPI_Comm comm,
	const std::filesystem::path& root,const CoupledCheckpointManifest& manifest,
	const MovingCheckpointLayout& source_layout,std::uint32_t source_membrane_owner,const std::string& configuration,
	const CubicCartesianGridSpec& grid,const MovingCutGeometryOptions& geometry_options,
	const ImmersedTransientFlowOptions& flow_options,const MaterialSurfacePatchMap& map,
	const DistributedSurfaceInterface& fluid_interface,const ImmersedMovingDistributedFsiRuntime::GeometryProvider& provider,
	int target_owner,ImmersedMovingDistributedFsiRuntime::Policy policy,
	PretensionedMembraneMaterial material,PretensionedMembraneOptions membrane_options={},PointIdentityLimits limits={})
{
	int rank=0,size=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&size);
	std::map<std::string,std::size_t> indices;std::unique_ptr<RestoredMovingFsiPair> result;std::string agreement;
	CollectiveLocalStage(comm,"paired checkpoint restore catalog",[&] {
		// Fluid fields and surface publications are restored independently into
		// target ownership; saved surface stamps are never silently relabeled.
		checkpoint_metadata::Require(target_owner>=0&&target_owner<size,"paired surface restore requires valid target owner");
		indices=moving_fsi_checkpoint_detail::Indices(manifest,source_layout,source_membrane_owner);
		Sha256 hash;hash.AppendLittleEndian64(source_membrane_owner);hash.AppendLittleEndian64(target_owner);
		distributed_surface_detail::AppendString(hash,SerializeCoupledCheckpointManifest(manifest));
		distributed_surface_detail::AppendString(hash,configuration);
		distributed_surface_detail::AppendString(hash,root.generic_string());
		distributed_surface_detail::AppendString(hash,source_layout.domain_id);
		for(auto member:source_layout.world_ranks)hash.AppendLittleEndian64(member);
		agreement=hash.Hex();result=std::make_unique<RestoredMovingFsiPair>();
	});
	RequireCollectiveSameText(comm,"paired checkpoint restore agreement",agreement);
	result->flow=RestoreMovingCheckpointRuntime(comm,root,manifest,source_layout,configuration,grid,geometry_options,flow_options);
	try {
		result->fluid=AllocateCollectiveRuntime<ImmersedMovingDistributedFsiRuntime>(comm,*result->flow,map,fluid_interface,provider,target_owner,policy,limits);
		std::optional<DistributedSurfaceLayout> partition;
		std::optional<DistributedSurfaceInterface> structure_interface,membrane_fluid_interface;
		std::vector<std::uint64_t> clamps;std::string metadata_identity;
		MovingFsiPublicationCheckpoint publication;
		std::optional<SingleOwnerMembraneCheckpoint> structure;
		CollectiveLocalStage(comm,"paired checkpoint restore payloads",[&] {
			partition.emplace(map.Layout());structure_interface.emplace(map.Interface());membrane_fluid_interface.emplace(fluid_interface);
			clamps=map.ConfiguredClampedGlobalNodeIds();
			const auto& geometry=result->flow->CommittedGeometry();
			bool direct=false;std::string bytes;
			if(source_layout.world_ranks.size()==static_cast<std::size_t>(size)) {
				bytes=moving_checkpoint_bundle_detail::Metadata(root,manifest,
					indices.at("fsi.rank-"+std::to_string(source_layout.world_ranks.at(rank))));
				const auto records=DecodeMovingFsiPublicationRecords(bytes,moving_fsi_checkpoint_detail::Hash(bytes),fluid_interface.id,
					result->fluid->CheckpointConfigurationIdentity(),geometry.Evaluation().ContentIdentitySha256(),geometry.GeometryIdentitySha256(),
					manifest.epoch.accepted_steps,manifest.epoch.time_s,policy.maximum_patch_nodes);
				direct=records.layout_identity==map.Layout().layout_identity_sha256
					&&records.partition_identity==BuildDistributedSurfacePartitionIdentitySha256(map.Layout())
					&&records.node_ids==map.Layout().owned_global_node_ids;
			}
			if(direct) {
				publication=ParseMovingFsiPublicationCheckpoint(bytes,moving_fsi_checkpoint_detail::Hash(bytes),map.Layout(),fluid_interface.id,
					result->fluid->CheckpointConfigurationIdentity(),geometry.Evaluation().ContentIdentitySha256(),geometry.GeometryIdentitySha256(),
					manifest.epoch.accepted_steps,manifest.epoch.time_s,policy.maximum_patch_nodes);
			} else publication=moving_fsi_checkpoint_detail::RepartitionPublication(root,manifest,source_layout,indices,map.Layout(),fluid_interface.id,
				result->fluid->CheckpointConfigurationIdentity(),geometry.Evaluation().ContentIdentitySha256(),
				geometry.GeometryIdentitySha256(),policy.maximum_patch_nodes);
			const auto metadata=moving_checkpoint_bundle_detail::Metadata(root,manifest,indices.at("membrane.metadata"));
			metadata_identity=moving_fsi_checkpoint_detail::Hash(metadata);
			if(rank==target_owner)structure.emplace(SingleOwnerMembraneCheckpoint{metadata,
				moving_checkpoint_bundle_detail::Metadata(root,manifest,indices.at("membrane.numerical"))});
		});
		result->membrane=AllocateCollectiveRuntime<SingleOwnerMembraneRuntime>(comm,comm,std::move(*partition),target_owner,
			std::move(*structure_interface),std::move(*membrane_fluid_interface),material,std::move(clamps),membrane_options,limits);
		result->fluid->RestoreCheckpointState(std::move(publication));
		result->membrane->RestoreCheckpoint(structure?&*structure:nullptr,metadata_identity);
		CollectiveLocalStage(comm,"paired checkpoint final context gate",[&] {
			moving_fsi_checkpoint_detail::RequireContext(result->fluid->CommittedContext(),result->membrane->CommittedContext(),manifest.epoch);
		});
		return result;
	} catch(...) {
		const auto failure=std::current_exception();
		try { result->Close(); }catch(...) {}
		std::rethrow_exception(failure);
	}
}

} // namespace iga
#endif

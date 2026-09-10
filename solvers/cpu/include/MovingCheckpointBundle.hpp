#ifndef IGA_MOVING_CHECKPOINT_BUNDLE_HPP
#define IGA_MOVING_CHECKPOINT_BUNDLE_HPP

#include "MovingAcceptedCheckpointMetadata.hpp"
#include "MaterialSurfaceCheckpoint.hpp"
#include "CoupledCheckpointBundle.hpp"

namespace iga {
struct MovingCheckpointLayout {
	std::string domain_id;
	std::vector<std::uint32_t> world_ranks;
};
inline std::vector<CoupledCheckpointShardSpec> MovingCheckpointCatalog(const MovingCheckpointLayout& layout,std::uint32_t world_size)
{
	using namespace coupled_checkpoint_detail;
	Require(Identifier(layout.domain_id)&&!layout.world_ranks.empty(),"invalid moving checkpoint layout");
	std::set<std::uint32_t> unique;
	for(auto rank:layout.world_ranks)Require(rank<world_size&&unique.insert(rank).second,"invalid moving rank mapping");
	const auto owner=layout.world_ranks.front();const auto prefix=layout.domain_id+".";
	std::vector<CoupledCheckpointShardSpec> catalog{{prefix+"metadata","moving-flow-v1",owner},
		{prefix+"conservation","moving-conservation-v1",owner},{prefix+"material.binding","moving-material-binding-v1",owner},
		{prefix+"material.previous","material-surface-v1",owner},{prefix+"material.current","material-surface-v1",owner}};
	for(auto rank:layout.world_ranks)catalog.push_back({prefix+"field.rank-"+std::to_string(rank),"owned-real-field-v1",rank});
	std::sort(catalog.begin(),catalog.end(),[](const auto& a,const auto& b) { return a.id<b.id; });
	ValidateCatalog(catalog,world_size);return catalog;
}
namespace moving_checkpoint_bundle_detail {
inline std::map<std::string,std::size_t> Indices(const CoupledCheckpointManifest& manifest,const MovingCheckpointLayout& layout)
{
	using namespace coupled_checkpoint_detail;Validate(manifest.epoch);
	std::map<std::string,std::size_t> indices;
	for(const auto& spec:MovingCheckpointCatalog(layout,manifest.epoch.compatibility.ranks)) {
		const auto found=std::find_if(manifest.shards.begin(),manifest.shards.end(),[&](const auto& shard) { return shard.spec.id==spec.id; });
		Require(found!=manifest.shards.end()&&Same(found->spec,spec),"moving manifest catalog differs");
		indices.emplace(spec.id.substr(layout.domain_id.size()+1),static_cast<std::size_t>(found-manifest.shards.begin()));
	}
	return indices;
}
inline std::string Metadata(const std::filesystem::path& root,const CoupledCheckpointManifest& manifest,std::size_t index)
{
	using checkpoint_metadata::Require;
	Require(manifest.shards.at(index).payload_bytes<=checkpoint_metadata::maximum_bytes,"moving metadata exceeds limit");
	std::string bytes;
	ReadCoupledCheckpointShard(root,manifest,index,[&](const void* data,std::size_t size) {
		Require(size<=checkpoint_metadata::maximum_bytes-bytes.size(),"moving metadata exceeds limit");bytes.append(static_cast<const char*>(data),size);
	});return bytes;
}
}
// LOCAL I/O only. The coordinator agrees on all receipts before publishing the
// complete bundle. No full field is gathered by these producers or consumers.
inline std::vector<CoupledCheckpointShard> WriteMovingCheckpointShards(const std::filesystem::path& root,
	const CoupledCheckpointEpoch& epoch,const MovingCheckpointLayout& layout,std::size_t local_rank,
	const ImmersedMovingAcceptedCheckpoint& state,const std::string& configuration_identity)
{
	using namespace coupled_checkpoint_detail;Validate(epoch);
	const auto catalog=MovingCheckpointCatalog(layout,epoch.compatibility.ranks);
	Require(IsLowercaseSha256(configuration_identity)&&local_rank<layout.world_ranks.size()&&state.previous_material&&state.current_material,"missing moving producer ownership or materials");
	Require(state.step==epoch.accepted_steps&&state.time_s==epoch.time_s
		&&state.conservation.source_time_s+epoch.dt_s==state.time_s,"moving producer epoch differs");
	Require(state.current_material->ContentIdentitySha256()==state.material_identity_sha256
		&&state.current_material->EvaluatedTimeS()==state.time_s
		&&state.previous_material->EvaluatedTimeS()==state.conservation.source_time_s
		&&state.current_material->StepStartS()==state.previous_material->EvaluatedTimeS()
		&&state.current_material->StepEndS()==state.time_s
		&&state.current_material->MaterialIdentitySha256()==state.previous_material->MaterialIdentitySha256()
		&&state.current_material->TopologyIdentitySha256()==state.previous_material->TopologyIdentitySha256(),"moving producer material binding differs");
	Require(state.field.row_begin<=state.field.row_end&&state.field.row_end<=state.field.global_rows
		&&state.field.values.size()==state.field.row_end-state.field.row_begin,"invalid moving producer field ownership");
	const auto owner=layout.world_ranks[local_rank];const auto prefix=layout.domain_id+".";
	std::vector<CoupledCheckpointShard> receipts;
	for(const auto& spec:catalog) {
		if(spec.owner_rank!=owner)continue;
		if(spec.id==prefix+"field.rank-"+std::to_string(owner)) {
			receipts.push_back(WriteCoupledCheckpointShard(root,epoch,spec,OwnedCheckpointFieldBytes(state.field),[&](auto& output) {
				WriteOwnedCheckpointField(state.field,[&](const void* bytes,std::size_t size) { output.Write(bytes,size); });
			}));continue;
		}
		std::string bytes;
		if(spec.id==prefix+"metadata")bytes=SerializeMovingAcceptedCheckpointMetadata(state,configuration_identity);
		else if(spec.id==prefix+"conservation")bytes=SerializeMovingConservationCheckpoint(state.conservation);
		else if(spec.id==prefix+"material.previous")bytes=SerializeMaterialSurfaceCheckpoint(*state.previous_material);
		else if(spec.id==prefix+"material.current")bytes=SerializeMaterialSurfaceCheckpoint(*state.current_material);
		else {
			checkpoint_metadata::Writer output;output.Text("IGA_MOVING_MATERIAL_BINDING/1");output.Text(configuration_identity);
			output.Unsigned(state.step);output.Real(state.time_s);output.Text(state.previous_material->ContentIdentitySha256());
			output.Text(state.current_material->ContentIdentitySha256());bytes=output.Bytes();
		}
		receipts.push_back(WriteCoupledCheckpointShard(root,epoch,spec,bytes.size(),[&](auto& output) { output.Write(bytes.data(),bytes.size()); }));
	}
	return receipts;
}
inline std::map<std::string,double> LoadMovingCheckpointPortControls(const std::filesystem::path& root,
	const CoupledCheckpointManifest& manifest,const MovingCheckpointLayout& source_layout,const std::string& expected_configuration)
{
	using namespace moving_checkpoint_bundle_detail;
	const auto indices=Indices(manifest,source_layout);
	const auto saved=DecodeMovingCheckpointMetadata(Metadata(root,manifest,indices.at("metadata")),expected_configuration);
	checkpoint_metadata::Require(saved.step==manifest.epoch.accepted_steps&&saved.time_s==manifest.epoch.time_s,
		"moving controls epoch differs from bundle");
	return saved.port_control_values;
}
struct MovingCheckpointMaterials {
	MaterialSurfaceKinematics previous,current;
};
inline MovingCheckpointMaterials LoadMovingCheckpointMaterials(const std::filesystem::path& root,
	const CoupledCheckpointManifest& manifest,const MovingCheckpointLayout& source_layout,const std::string& expected_configuration)
{
	using namespace moving_checkpoint_bundle_detail;using checkpoint_metadata::Require;
	const auto indices=Indices(manifest,source_layout);
	const auto bytes=Metadata(root,manifest,indices.at("material.binding"));checkpoint_metadata::Reader input(bytes);
	Require(IsLowercaseSha256(expected_configuration)&&input.Text()=="IGA_MOVING_MATERIAL_BINDING/1"
		&&input.Text()==expected_configuration,"moving material configuration differs");
	Require(input.Unsigned()==manifest.epoch.accepted_steps&&input.Real()==manifest.epoch.time_s,"moving material epoch differs");
	const auto previous_identity=input.Text(),current_identity=input.Text();input.Finish();
	auto previous=ParseMaterialSurfaceCheckpoint(Metadata(root,manifest,indices.at("material.previous")),previous_identity);
	auto current=ParseMaterialSurfaceCheckpoint(Metadata(root,manifest,indices.at("material.current")),current_identity);
	Require(current.EvaluatedTimeS()==manifest.epoch.time_s&&current.StepEndS()==manifest.epoch.time_s
		&&previous.EvaluatedTimeS()+manifest.epoch.dt_s==manifest.epoch.time_s
		&&current.StepStartS()==previous.EvaluatedTimeS(),"moving material interval differs");
	return {std::move(previous),std::move(current)};
}
// Target geometry/material identities and clock are supplied by the fresh
// runtime construction. Source layout is the saved mapping, not the new ranks.
// On failure discard the unpublished target; no runtime is modified here.
inline void LoadMovingCheckpointFields(const std::filesystem::path& root,const CoupledCheckpointManifest& manifest,
	const MovingCheckpointLayout& source_layout,const std::string& expected_configuration,ImmersedMovingAcceptedCheckpoint& target)
{
	using namespace moving_checkpoint_bundle_detail;using checkpoint_metadata::Require;
	const auto indices=Indices(manifest,source_layout);
	Require(target.step==manifest.epoch.accepted_steps&&target.time_s==manifest.epoch.time_s,"moving target epoch differs from bundle");
	const auto authority=ParseMovingAcceptedCheckpointMetadata(Metadata(root,manifest,indices.at("metadata")),target,expected_configuration);
	target.conservation=ParseMovingConservationCheckpoint(Metadata(root,manifest,indices.at("conservation")),authority,
		target.geometry_identity_sha256,target.publication_identity_sha256,target.step,target.time_s);
	Require(target.conservation.source_time_s+manifest.epoch.dt_s==target.time_s,"moving conservation interval differs from bundle");
	OwnedCheckpointRepartitionReader reader(target.field);
	for(auto rank:source_layout.world_ranks) {
		const auto index=indices.at("field.rank-"+std::to_string(rank));
		reader.ReadShard(manifest.shards[index].payload_bytes,[&](auto consume) { ReadCoupledCheckpointShard(root,manifest,index,consume); });
	}
	reader.Finish();
}
} // namespace iga
#endif

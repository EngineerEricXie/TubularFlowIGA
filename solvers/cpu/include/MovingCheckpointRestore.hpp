#ifndef IGA_MOVING_CHECKPOINT_RESTORE_HPP
#define IGA_MOVING_CHECKPOINT_RESTORE_HPP

#include "MovingCheckpointBundle.hpp"

namespace iga {
// Validate the complete saved catalog before changing any candidate option.
inline void ApplyMovingCheckpointPortControls(ImmersedTransientFlowOptions& options,const std::map<std::string,double>& controls)
{
	if(controls.size()!=options.ports.size())throw std::runtime_error("moving saved port catalog differs from target");
	std::set<std::string> ids;
	for(const auto& port:options.ports) {
		const auto found=controls.find(port.id);
		if(found==controls.end()||!std::isfinite(found->second)||!ids.insert(port.id).second)
			throw std::runtime_error("invalid moving saved port or target catalog");
	}
	for(auto& port:options.ports)port.value=controls.at(port.id);
}
// Collective fresh-runtime factory. The caller supplies a verified immutable
// manifest and case configuration; numerical fields never pass through root.
// Publish this returned owner only after other coupled candidates also succeed.
inline std::unique_ptr<ImmersedMovingTransientDistributedRuntime> RestoreMovingCheckpointRuntime(
	MPI_Comm comm,const std::filesystem::path& root,const CoupledCheckpointManifest& manifest,
	const MovingCheckpointLayout& source_layout,const std::string& configuration_identity,
	const CubicCartesianGridSpec& grid,const MovingCutGeometryOptions& geometry_options,
	const ImmersedTransientFlowOptions& flow_options)
{
	std::string agreement;
	CollectiveLocalStage(comm,"moving checkpoint factory input",[&] {
		if(!IsLowercaseSha256(configuration_identity))throw std::invalid_argument("invalid moving restore configuration identity");
		(void)MovingCheckpointCatalog(source_layout,manifest.epoch.compatibility.ranks);
		Sha256 hash;const auto bytes=SerializeCoupledCheckpointManifest(manifest);
		distributed_surface_detail::AppendString(hash,bytes);
		distributed_surface_detail::AppendString(hash,root.generic_string());
		distributed_surface_detail::AppendString(hash,configuration_identity);
		distributed_surface_detail::AppendString(hash,source_layout.domain_id);
		hash.AppendLittleEndian64(source_layout.world_ranks.size());
		for(auto rank:source_layout.world_ranks)hash.AppendLittleEndian32(rank);
		agreement=hash.Hex();
	});
	RequireCollectiveSameText(comm,"moving checkpoint factory agreement",agreement);
	std::optional<MovingCheckpointMaterials> materials;
	std::optional<ImmersedTransientFlowOptions> restored_options;
	std::unique_ptr<MovingCutGeometry> geometry;
	std::string source_geometry,target_geometry;
	CollectiveLocalStage(comm,"moving checkpoint factory geometry",[&] {
		materials.emplace(LoadMovingCheckpointMaterials(root,manifest,source_layout,configuration_identity));
		const auto controls=LoadMovingCheckpointPortControls(root,manifest,source_layout,configuration_identity);
		restored_options.emplace(flow_options);
		ApplyMovingCheckpointPortControls(*restored_options,controls);
		auto previous=MovingCutGeometry::Build(grid,materials->previous,geometry_options);
		source_geometry=previous->GeometryIdentitySha256();
		geometry=MovingCutGeometry::Build(grid,materials->current,geometry_options,previous.get());
		const auto indices=moving_checkpoint_bundle_detail::Indices(manifest,source_layout);
		const auto saved=DecodeMovingCheckpointMetadata(moving_checkpoint_bundle_detail::Metadata(root,manifest,indices.at("metadata")),configuration_identity);
		// Providers may publish an independently built target without predecessor
		// transition metadata. Accept only the exact authenticated saved identity.
		if(geometry->PublicationIdentitySha256()!=saved.identities[1])
			geometry=MovingCutGeometry::Build(grid,materials->current,geometry_options);
		if(geometry->GeometryIdentitySha256()!=saved.identities[0]||geometry->PublicationIdentitySha256()!=saved.identities[1])
			throw std::runtime_error("moving checkpoint geometry reconstruction differs from saved identity");
		target_geometry=geometry->GeometryIdentitySha256();
	});
	RequireCollectiveSameText(comm,"moving checkpoint source geometry agreement",source_geometry);
	RequireCollectiveSameText(comm,"moving checkpoint target geometry agreement",target_geometry);
	auto runtime=AllocateCollectiveRuntime<ImmersedMovingTransientDistributedRuntime>(comm,comm,
		std::move(geometry),*restored_options,manifest.epoch.accepted_steps);
	try {
		ImmersedMovingAcceptedCheckpoint state;
		CollectiveLocalStage(comm,"moving checkpoint factory fields",[&] {
			state.field=CaptureOwnedCheckpointVector(runtime->CommittedState());
			for(const auto& port:runtime->PortDefinitions())state.port_control_values.emplace(port.id,port.value);
			state.geometry_identity_sha256=runtime->CommittedGeometry().GeometryIdentitySha256();
			state.publication_identity_sha256=runtime->CommittedGeometry().PublicationIdentitySha256();
			state.layout_identity_sha256=runtime->CommittedLayout().HashSha256();
			state.material_identity_sha256=materials->current.ContentIdentitySha256();
			state.step=runtime->Clock().index;state.time_s=runtime->Clock().time_s;
			state.previous_material.emplace(std::move(materials->previous));state.current_material.emplace(std::move(materials->current));
			LoadMovingCheckpointFields(root,manifest,source_layout,configuration_identity,state);
			if(state.conservation.source_geometry_identity_sha256!=source_geometry)
				throw std::runtime_error("moving checkpoint predecessor geometry differs from conservation");
		});
		runtime->RestoreAcceptedCheckpoint(state);return runtime;
	} catch(...) {
		const auto failure=std::current_exception();
		try { runtime->Close(); } catch(...) {}
		std::rethrow_exception(failure);
	}
}
} // namespace iga
#endif

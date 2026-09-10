#ifndef IGA_MOVING_CONSERVATION_CHECKPOINT_HPP
#define IGA_MOVING_CONSERVATION_CHECKPOINT_HPP

#include "ImmersedMovingDistributedConservation.hpp"
#include "CheckpointMetadataCodec.hpp"

namespace iga {
// Replicated transition diagnostics only; numerical fields use owned shards.
// The external bundle hash authenticates every saved aggregate and identity.
inline std::string SerializeMovingConservationCheckpoint(const ImmersedMovingDistributedConservation& state)
{
	using checkpoint_metadata::Require;
	Require(state.source_index<std::numeric_limits<std::uint64_t>::max()
		&&state.target_index==state.source_index+1&&state.source_time_s>=0
		&&state.dt_s>0&&state.target_time_s==state.source_time_s+state.dt_s,
		"invalid moving conservation epoch");
	Require(state.source_audited_volume_m3>=0&&state.target_audited_volume_m3>=0
		&&state.normalization_scale_m3_s>0,"invalid moving conservation scale or volume");
	for(double value:{state.normalized_divergence_theorem_defect,state.normalized_reynolds_defect,
		state.normalized_moving_mass_defect,state.normalized_wall_relative_leakage,
		state.endpoint.absolute_surface_flux_sum_m3_s,state.endpoint.normalized_open_balance,
		state.endpoint.normalized_wall_leakage,state.endpoint.normalized_discrete_moving_wall_continuity_defect})
		Require(value>=0,"negative moving conservation diagnostic");
	Require(state.endpoint.discrete_moving_wall_continuity_normalization_scale_m3_s>0,
		"invalid endpoint continuity normalization scale");
	checkpoint_metadata::Writer output;output.Text("IGA_MOVING_CONSERVATION/1");
	Require(IsLowercaseSha256(state.source_geometry_identity_sha256),"invalid moving conservation identity");output.Text(state.source_geometry_identity_sha256);
	Require(IsLowercaseSha256(state.target_geometry_identity_sha256),"invalid moving conservation identity");output.Text(state.target_geometry_identity_sha256);
	Require(IsLowercaseSha256(state.source_publication_identity_sha256),"invalid moving conservation identity");output.Text(state.source_publication_identity_sha256);
	Require(IsLowercaseSha256(state.target_publication_identity_sha256),"invalid moving conservation identity");output.Text(state.target_publication_identity_sha256);
	output.Unsigned(state.source_index);output.Unsigned(state.target_index);
	output.Real(state.source_time_s);
	output.Real(state.target_time_s);
	output.Real(state.dt_s);
	output.Real(state.source_audited_volume_m3);
	output.Real(state.target_audited_volume_m3);
	output.Real(state.backward_euler_volume_rate_m3_s);
	output.Real(state.reynolds_defect_m3_s);
	output.Real(state.moving_mass_defect_m3_s);
	output.Real(state.normalization_scale_m3_s);
	output.Real(state.normalized_divergence_theorem_defect);
	output.Real(state.normalized_reynolds_defect);
	output.Real(state.normalized_moving_mass_defect);
	output.Real(state.normalized_wall_relative_leakage);
	output.Unsigned(state.endpoint.surface_flux_term_count);
	output.Real(state.endpoint.absolute_surface_flux_sum_m3_s);
	output.Real(state.endpoint.endpoint_volume_divergence_m3_s);
	output.Real(state.endpoint.total_surface_outward_flow_m3_s);
	output.Real(state.endpoint.open_port_outward_flow_m3_s);
	output.Real(state.endpoint.wall_outward_flow_m3_s);
	output.Real(state.endpoint.total_material_surface_outward_flow_m3_s);
	output.Real(state.endpoint.total_material_wall_outward_flow_m3_s);
	output.Real(state.endpoint.wall_relative_leakage_m3_s);
	output.Real(state.endpoint.divergence_theorem_defect_m3_s);
	output.Real(state.endpoint.discrete_moving_wall_continuity_defect_m3_s);
	output.Real(state.endpoint.discrete_moving_wall_continuity_normalization_scale_m3_s);
	output.Real(state.endpoint.normalized_open_balance);
	output.Real(state.endpoint.normalized_wall_leakage);
	output.Real(state.endpoint.normalized_discrete_moving_wall_continuity_defect);
	Require(state.endpoint.surface_flow_by_boundary_label_m3_s.size()<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
	output.Unsigned(state.endpoint.surface_flow_by_boundary_label_m3_s.size());
	for(const auto& entry:state.endpoint.surface_flow_by_boundary_label_m3_s) {
		Require(entry.first>=0,"negative conservation label");output.Unsigned(entry.first);output.Real(entry.second);
	}
	Require(state.endpoint.material_surface_outward_flow_by_boundary_label_m3_s.size()<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
	output.Unsigned(state.endpoint.material_surface_outward_flow_by_boundary_label_m3_s.size());
	for(const auto& entry:state.endpoint.material_surface_outward_flow_by_boundary_label_m3_s) {
		Require(entry.first>=0,"negative conservation label");output.Unsigned(entry.first);output.Real(entry.second);
	}
	Require(state.endpoint.material_wall_outward_flow_by_boundary_label_m3_s.size()<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
	output.Unsigned(state.endpoint.material_wall_outward_flow_by_boundary_label_m3_s.size());
	for(const auto& entry:state.endpoint.material_wall_outward_flow_by_boundary_label_m3_s) {
		Require(entry.first>=0,"negative conservation label");output.Unsigned(entry.first);output.Real(entry.second);
	}
	return output.Bytes();
}

inline ImmersedMovingDistributedConservation ParseMovingConservationCheckpoint(std::string_view bytes,
	const std::string& expected_payload_identity,const std::string& expected_geometry,
	const std::string& expected_publication,std::uint64_t expected_step,double expected_time)
{
	using checkpoint_metadata::Require;
	Require(bytes.size()<=checkpoint_metadata::maximum_bytes,"moving conservation exceeds byte limit");
	Sha256 hash;hash.Append(bytes.data(),bytes.size());
	Require(IsLowercaseSha256(expected_payload_identity)&&hash.Hex()==expected_payload_identity,"moving conservation payload authority differs");
	checkpoint_metadata::Reader input(bytes);Require(input.Text()=="IGA_MOVING_CONSERVATION/1","unsupported moving conservation checkpoint");
	ImmersedMovingDistributedConservation state;
	state.source_geometry_identity_sha256=input.Text();
	state.target_geometry_identity_sha256=input.Text();
	state.source_publication_identity_sha256=input.Text();
	state.target_publication_identity_sha256=input.Text();
	state.source_index=input.Unsigned();state.target_index=input.Unsigned();
	state.source_time_s=input.Real();
	state.target_time_s=input.Real();
	state.dt_s=input.Real();
	state.source_audited_volume_m3=input.Real();
	state.target_audited_volume_m3=input.Real();
	state.backward_euler_volume_rate_m3_s=input.Real();
	state.reynolds_defect_m3_s=input.Real();
	state.moving_mass_defect_m3_s=input.Real();
	state.normalization_scale_m3_s=input.Real();
	state.normalized_divergence_theorem_defect=input.Real();
	state.normalized_reynolds_defect=input.Real();
	state.normalized_moving_mass_defect=input.Real();
	state.normalized_wall_relative_leakage=input.Real();
	state.endpoint.surface_flux_term_count=input.Unsigned();
	state.endpoint.absolute_surface_flux_sum_m3_s=input.Real();
	state.endpoint.endpoint_volume_divergence_m3_s=input.Real();
	state.endpoint.total_surface_outward_flow_m3_s=input.Real();
	state.endpoint.open_port_outward_flow_m3_s=input.Real();
	state.endpoint.wall_outward_flow_m3_s=input.Real();
	state.endpoint.total_material_surface_outward_flow_m3_s=input.Real();
	state.endpoint.total_material_wall_outward_flow_m3_s=input.Real();
	state.endpoint.wall_relative_leakage_m3_s=input.Real();
	state.endpoint.divergence_theorem_defect_m3_s=input.Real();
	state.endpoint.discrete_moving_wall_continuity_defect_m3_s=input.Real();
	state.endpoint.discrete_moving_wall_continuity_normalization_scale_m3_s=input.Real();
	state.endpoint.normalized_open_balance=input.Real();
	state.endpoint.normalized_wall_leakage=input.Real();
	state.endpoint.normalized_discrete_moving_wall_continuity_defect=input.Real();
	{
		const auto count=input.Unsigned();Require(count<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
		int previous=-1;
		for(std::uint64_t i=0;i<count;++i) {
			const auto label=input.Unsigned();Require(label<=static_cast<std::uint64_t>(std::numeric_limits<int>::max())
				&&static_cast<int>(label)>previous,"invalid conservation label order");
			previous=static_cast<int>(label);state.endpoint.surface_flow_by_boundary_label_m3_s.emplace(previous,input.Real());
		}
	}
	{
		const auto count=input.Unsigned();Require(count<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
		int previous=-1;
		for(std::uint64_t i=0;i<count;++i) {
			const auto label=input.Unsigned();Require(label<=static_cast<std::uint64_t>(std::numeric_limits<int>::max())
				&&static_cast<int>(label)>previous,"invalid conservation label order");
			previous=static_cast<int>(label);state.endpoint.material_surface_outward_flow_by_boundary_label_m3_s.emplace(previous,input.Real());
		}
	}
	{
		const auto count=input.Unsigned();Require(count<=checkpoint_metadata::maximum_map_entries,"too many conservation labels");
		int previous=-1;
		for(std::uint64_t i=0;i<count;++i) {
			const auto label=input.Unsigned();Require(label<=static_cast<std::uint64_t>(std::numeric_limits<int>::max())
				&&static_cast<int>(label)>previous,"invalid conservation label order");
			previous=static_cast<int>(label);state.endpoint.material_wall_outward_flow_by_boundary_label_m3_s.emplace(previous,input.Real());
		}
	}
	input.Finish();
	Require(state.target_geometry_identity_sha256==expected_geometry
		&&state.target_publication_identity_sha256==expected_publication
		&&state.target_index==expected_step&&state.target_time_s==expected_time,"moving conservation target differs");
	Require(SerializeMovingConservationCheckpoint(state)==bytes,"noncanonical moving conservation payload");
	return state;
}
} // namespace iga
#endif

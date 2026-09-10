#ifndef IGA_MOVING_ACCEPTED_CHECKPOINT_METADATA_HPP
#define IGA_MOVING_ACCEPTED_CHECKPOINT_METADATA_HPP

#include "ImmersedMovingTransientDistributedRuntime.hpp"
#include "OwnedCheckpointFieldStream.hpp"

namespace iga {
// Shared metadata only. Configuration comes from the outer case/bundle
// authority; field shards and this metadata must also pass bundle SHA checks.
inline std::string SerializeMovingAcceptedCheckpointMetadata(const ImmersedMovingAcceptedCheckpoint& state,
	const std::string& configuration_identity)
{
	using checkpoint_metadata::Require;
	Require(IsLowercaseSha256(configuration_identity),"invalid moving configuration identity");
	Require(state.step>0&&state.time_s>=0&&state.field.global_rows>0,"invalid moving accepted epoch or shape");
	Require(state.conservation.target_index==state.step&&state.conservation.target_time_s==state.time_s
		&&state.conservation.target_geometry_identity_sha256==state.geometry_identity_sha256
		&&state.conservation.target_publication_identity_sha256==state.publication_identity_sha256,
		"moving accepted conservation binding differs");
	const auto conservation=SerializeMovingConservationCheckpoint(state.conservation);
	Sha256 hash;hash.Append(conservation.data(),conservation.size());
	checkpoint_metadata::Writer output;output.Text("IGA_MOVING_ACCEPTED/2");output.Text(configuration_identity);
	for(const auto* identity:{&state.geometry_identity_sha256,&state.publication_identity_sha256,
		&state.layout_identity_sha256,&state.material_identity_sha256}) {
		Require(IsLowercaseSha256(*identity),"invalid moving accepted identity");output.Text(*identity);
	}
	output.Unsigned(state.step);output.Real(state.time_s);output.Unsigned(state.field.global_rows);
	output.Text(hash.Hex());output.Reals(state.port_control_values);return output.Bytes();
}

// Descriptor can be read before constructing the fresh numerical owner.
struct MovingCheckpointMetadata {
	std::array<std::string,4> identities;
	std::uint64_t step=0,global_rows=0;
	double time_s=0;
	std::string conservation_identity;
	std::map<std::string,double> port_control_values;
};
inline MovingCheckpointMetadata DecodeMovingCheckpointMetadata(std::string_view bytes,const std::string& expected_configuration_identity)
{
	using checkpoint_metadata::Require;
	Require(IsLowercaseSha256(expected_configuration_identity),"invalid expected moving configuration identity");
	checkpoint_metadata::Reader input(bytes);const auto version=input.Text();
	Require(version=="IGA_MOVING_ACCEPTED/1"||version=="IGA_MOVING_ACCEPTED/2","unsupported moving accepted checkpoint");
	Require(input.Text()==expected_configuration_identity,"moving checkpoint configuration differs");
	MovingCheckpointMetadata result;
	for(auto& identity:result.identities) { identity=input.Text();Require(IsLowercaseSha256(identity),"invalid moving target identity"); }
	result.step=input.Unsigned();result.time_s=input.Real();result.global_rows=input.Unsigned();
	Require(result.step>0&&result.time_s>=0&&result.global_rows>0,"invalid moving accepted epoch or shape");
	result.conservation_identity=input.Text();Require(IsLowercaseSha256(result.conservation_identity),"invalid moving conservation authority");
	// Legacy /1 did not preserve controls, so it can restore only port-free cases.
	if(version=="IGA_MOVING_ACCEPTED/2")result.port_control_values=input.Reals();
	input.Finish();return result;
}
// Target comes from verified geometry/configuration at the saved epoch, with
// saved controls applied during fresh construction. No state is published here.
inline std::string ParseMovingAcceptedCheckpointMetadata(std::string_view bytes,
	const ImmersedMovingAcceptedCheckpoint& target,const std::string& expected_configuration_identity)
{
	using checkpoint_metadata::Require;
	const auto saved=DecodeMovingCheckpointMetadata(bytes,expected_configuration_identity);
	const std::array<std::string,4> identities{{target.geometry_identity_sha256,target.publication_identity_sha256,
		target.layout_identity_sha256,target.material_identity_sha256}};
	Require(saved.identities==identities,"moving checkpoint target identity differs");
	Require(saved.step==target.step&&saved.time_s==target.time_s&&saved.global_rows==target.field.global_rows,
		"moving checkpoint target epoch or global shape differs");
	Require(saved.port_control_values==target.port_control_values,"moving checkpoint target port controls differ");
	return saved.conservation_identity;
}
} // namespace iga
#endif

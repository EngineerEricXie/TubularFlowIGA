#ifndef IGA_ZERO_D_FLOW_SPECIES_CHECKPOINT_HPP
#define IGA_ZERO_D_FLOW_SPECIES_CHECKPOINT_HPP

#include "CheckpointMetadataCodec.hpp"
#include "CoupledCheckpointManifest.hpp"
#include "ZeroDFlowDomain.hpp"
#include "ZeroDSpeciesReservoir.hpp"

#include <set>
#include <string>
#include <string_view>

namespace iga {

struct ZeroDFlowSpeciesCheckpointState
{
	std::string domain_id;
	std::string model_identity_sha256;
	DomainStepContext accepted_step;
	ZeroDFlowState hydraulic;
	ZeroDSpeciesReservoirState species;
};

inline std::string BuildZeroDFlowSpeciesModelIdentitySha256(
	const ZeroDFlowModel& hydraulic,
	const ZeroDSpeciesReservoirModel& species,
	const std::map<std::string,double>& external_donor)
{
	ValidateZeroDFlowModel(hydraulic);
	ValidateZeroDSpeciesReservoirModel(species);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash,"ZeroDFlowSpeciesModel/v1");
	zero_d_flow_detail::AppendString(hash,BuildZeroDFlowModelIdentitySha256(hydraulic));
	hash.AppendLittleEndian64(species.port_ids.size());
	for(const auto& id:species.port_ids)
		zero_d_flow_detail::AppendString(hash,id);
	hash.AppendLittleEndian64(species.species_ids.size());
	for(const auto& id:species.species_ids)
		zero_d_flow_detail::AppendString(hash,id);
	hash.AppendLittleEndian64(species.source_rate_mol_s.size());
	for(const auto& item:species.source_rate_mol_s){
		zero_d_flow_detail::AppendString(hash,item.first);
		hash.AppendNormalizedDouble(item.second);
	}
	hash.AppendLittleEndian64(external_donor.size());
	for(const auto& item:external_donor){
		zero_d_flow_detail::AppendString(hash,item.first);
		hash.AppendNormalizedDouble(item.second);
	}
	return hash.Hex();
}

inline void ValidateZeroDFlowSpeciesCheckpointState(
	const ZeroDFlowSpeciesCheckpointState& value)
{
	checkpoint_metadata::Require(!value.domain_id.empty()
		&&value.model_identity_sha256.size()==64
		&&value.model_identity_sha256.find_first_not_of("0123456789abcdef")
			==std::string::npos,"invalid 0D species checkpoint identity");
	value.accepted_step.Validate();
	ValidateZeroDFlowState(value.hydraulic);
	checkpoint_metadata::Require(value.species.volume_m3>0.0
		&&std::isfinite(value.species.volume_m3)
		&&!value.species.amount_mol.empty(),"invalid 0D species checkpoint volume or species set");
	for(const auto& item:value.species.amount_mol)
		checkpoint_metadata::Require(!item.first.empty()&&item.second>=0.0
			&&std::isfinite(item.second),"invalid 0D species checkpoint amount");
}

inline std::string SerializeZeroDFlowSpeciesCheckpoint(
	const ZeroDFlowSpeciesCheckpointState& value)
{
	ValidateZeroDFlowSpeciesCheckpointState(value);
	checkpoint_metadata::Writer output;
	output.Text("IGA_ZERO_D_FLOW_SPECIES_ACCEPTED");output.Unsigned(1);
	output.Text(value.domain_id);output.Text(value.model_identity_sha256);
	checkpoint_metadata::WriteStep(output,value.accepted_step);
	output.Real(value.hydraulic.stored_pressure_pa);
	output.Real(value.species.volume_m3);
	output.Reals(value.species.amount_mol);
	Sha256 hash;hash.Append(output.Bytes().data(),output.Bytes().size());
	checkpoint_metadata::Writer checksum;checksum.Text(hash.Hex());
	checkpoint_metadata::Require(output.Bytes().size()+checksum.Bytes().size()
		<=checkpoint_metadata::maximum_bytes,"0D species checkpoint exceeds size limit");
	return output.Bytes()+checksum.Bytes();
}

inline ZeroDFlowSpeciesCheckpointState ParseZeroDFlowSpeciesCheckpoint(
	std::string_view bytes)
{
	using namespace checkpoint_metadata;
	constexpr std::size_t checksum_bytes=8+64;
	Require(bytes.size()>=checksum_bytes&&bytes.size()<=maximum_bytes,
		"0D species checkpoint size is invalid");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);
	Reader checksum(bytes.substr(bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();
	Sha256 hash;hash.Append(payload.data(),payload.size());
	Require(expected==hash.Hex(),"0D species checkpoint checksum mismatch");
	Reader input(payload);
	Require(input.Text()=="IGA_ZERO_D_FLOW_SPECIES_ACCEPTED"&&input.Unsigned()==1,
		"unsupported 0D species checkpoint schema");
	ZeroDFlowSpeciesCheckpointState result;
	result.domain_id=input.Text();result.model_identity_sha256=input.Text();
	result.accepted_step=ReadStep(input);
	result.hydraulic.stored_pressure_pa=input.Real();
	result.species.volume_m3=input.Real();
	result.species.amount_mol=input.Reals();
	input.Finish();ValidateZeroDFlowSpeciesCheckpointState(result);
	return result;
}

inline ZeroDFlowSpeciesCheckpointState ParseZeroDFlowSpeciesCheckpoint(
	std::string_view bytes,const CoupledCheckpointEpoch& epoch)
{
	coupled_checkpoint_detail::Validate(epoch);
	auto result=ParseZeroDFlowSpeciesCheckpoint(bytes);
	checkpoint_metadata::RequireAcceptedStep(result.accepted_step,
		epoch.accepted_steps,epoch.time_s,epoch.dt_s);
	return result;
}

} // namespace iga

#endif

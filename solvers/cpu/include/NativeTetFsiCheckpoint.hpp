#ifndef IGA_NATIVE_TET_FSI_CHECKPOINT_HPP
#define IGA_NATIVE_TET_FSI_CHECKPOINT_HPP

// Single-payload checkpoint for the currently single-partition native
// tetrahedral ALE-fluid/solid vertical slice. The payload is captured only
// between macro steps and binds both model and committed-state identities.
#include "CheckpointMetadataCodec.hpp"
#include "NativeTetAleFsiRuntime.hpp"
#include "NativeTetSolidFsiRuntime.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace iga {

struct NativeTetFsiCheckpoint
{
	std::uint64_t accepted_steps=0;
	double time_s=0.0;
	std::string fluid_model_identity_sha256,solid_model_identity_sha256;
	std::string fluid_state_identity_sha256,solid_state_identity_sha256;
	std::vector<double> fluid_state,solid_displacement_m,solid_velocity_m_per_s;
	std::vector<std::array<double,3>> ale_displacement_m,interface_displacement_m;
};

namespace native_tet_fsi_checkpoint_detail {

constexpr std::uint64_t maximum_scalars=1000000;

inline bool IsSha256(const std::string& value)
{
	if(value.size()!=64)return false;
	for(char c:value)if(!((c>='0'&&c<='9')||(c>='a'&&c<='f')))return false;
	return true;
}

inline void Validate(const NativeTetFsiCheckpoint& value)
{
	using checkpoint_metadata::Require;
	Require(std::isfinite(value.time_s)&&value.time_s>=0.0,"native FSI checkpoint time is invalid");
	Require(IsSha256(value.fluid_model_identity_sha256)&&IsSha256(value.solid_model_identity_sha256)
		&&IsSha256(value.fluid_state_identity_sha256)&&IsSha256(value.solid_state_identity_sha256),
		"native FSI checkpoint identity is invalid");
	const std::uint64_t scalar_count=value.fluid_state.size()+value.solid_displacement_m.size()
		+value.solid_velocity_m_per_s.size()+3*value.ale_displacement_m.size()
		+3*value.interface_displacement_m.size();
	Require(scalar_count>0&&scalar_count<=maximum_scalars,"native FSI checkpoint state exceeds limit");
	Require(!value.fluid_state.empty()&&!value.solid_displacement_m.empty()
		&&value.solid_displacement_m.size()==value.solid_velocity_m_per_s.size()
		&&!value.ale_displacement_m.empty()&&!value.interface_displacement_m.empty(),
		"native FSI checkpoint state shape is invalid");
	for(double scalar:value.fluid_state)Require(std::isfinite(scalar),"native FSI checkpoint contains nonfinite fluid state");
	for(double scalar:value.solid_displacement_m)Require(std::isfinite(scalar),"native FSI checkpoint contains nonfinite solid state");
	for(double scalar:value.solid_velocity_m_per_s)Require(std::isfinite(scalar),"native FSI checkpoint contains nonfinite solid velocity");
	for(const auto& vector:value.ale_displacement_m)for(double scalar:vector)
		Require(std::isfinite(scalar),"native FSI checkpoint contains nonfinite ALE state");
	for(const auto& vector:value.interface_displacement_m)for(double scalar:vector)
		Require(std::isfinite(scalar),"native FSI checkpoint contains nonfinite interface state");
}

inline void WriteScalars(checkpoint_metadata::Writer& output,const std::vector<double>& values)
{
	output.Unsigned(values.size());for(double value:values)output.Real(value);
}
inline std::vector<double> ReadScalars(checkpoint_metadata::Reader& input)
{
	const auto count=input.Unsigned();checkpoint_metadata::Require(count<=maximum_scalars,"native FSI checkpoint vector exceeds limit");
	std::vector<double> result(static_cast<std::size_t>(count));for(double& value:result)value=input.Real();return result;
}
inline void WriteVectors(checkpoint_metadata::Writer& output,const std::vector<std::array<double,3>>& values)
{
	output.Unsigned(values.size());for(const auto& value:values)for(double scalar:value)output.Real(scalar);
}
inline std::vector<std::array<double,3>> ReadVectors(checkpoint_metadata::Reader& input)
{
	const auto count=input.Unsigned();checkpoint_metadata::Require(count<=maximum_scalars/3,"native FSI checkpoint vector field exceeds limit");
	std::vector<std::array<double,3>> result(static_cast<std::size_t>(count));
	for(auto& value:result)for(double& scalar:value)scalar=input.Real();
	return result;
}
inline std::string Payload(const NativeTetFsiCheckpoint& value)
{
	Validate(value);checkpoint_metadata::Writer output;output.Text("IGA_NATIVE_TET_FSI");output.Unsigned(1);
	output.Unsigned(value.accepted_steps);output.Real(value.time_s);
	output.Text(value.fluid_model_identity_sha256);output.Text(value.solid_model_identity_sha256);
	output.Text(value.fluid_state_identity_sha256);output.Text(value.solid_state_identity_sha256);
	WriteScalars(output,value.fluid_state);WriteVectors(output,value.ale_displacement_m);
	WriteVectors(output,value.interface_displacement_m);WriteScalars(output,value.solid_displacement_m);
	WriteScalars(output,value.solid_velocity_m_per_s);return output.Bytes();
}
inline std::string Digest(std::string_view bytes)
{
	Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex();
}

} // namespace native_tet_fsi_checkpoint_detail

inline NativeTetFsiCheckpoint CaptureNativeTetFsiCheckpoint(std::uint64_t accepted_steps,
	const NativeTetAleFsiRuntime& fluid,const NativeTetSolidFsiRuntime& solid)
{
	if(fluid.Lifecycle().HasActiveStep()||solid.Lifecycle().HasActiveStep())
		throw std::runtime_error("native tetrahedral FSI checkpoint requires an idle accepted-step boundary");
	NativeTetFsiCheckpoint value;value.accepted_steps=accepted_steps;value.time_s=fluid.CommittedAleTimeS();
	value.fluid_model_identity_sha256=fluid.ModelIdentitySha256();
	value.solid_model_identity_sha256=solid.ModelIdentitySha256();
	value.fluid_state_identity_sha256=fluid.CommittedCompositionIdentitySha256();
	value.solid_state_identity_sha256=solid.CommittedStateIdentitySha256();
	value.fluid_state=fluid.CommittedFluidState();value.ale_displacement_m=fluid.CommittedAleDisplacementM();
	value.interface_displacement_m=fluid.CommittedInterfaceDisplacementM();
	value.solid_displacement_m=solid.CommittedDisplacementM();value.solid_velocity_m_per_s=solid.CommittedVelocityMPerS();
	native_tet_fsi_checkpoint_detail::Validate(value);return value;
}

inline std::string SerializeNativeTetFsiCheckpoint(const NativeTetFsiCheckpoint& value)
{
	auto payload=native_tet_fsi_checkpoint_detail::Payload(value);checkpoint_metadata::Writer checksum;
	checksum.Text(native_tet_fsi_checkpoint_detail::Digest(payload));payload+=checksum.Bytes();return payload;
}

inline NativeTetFsiCheckpoint ParseNativeTetFsiCheckpoint(std::string_view bytes)
{
	using checkpoint_metadata::Require;constexpr std::size_t checksum_bytes=72;
	Require(bytes.size()>=checksum_bytes,"native FSI checkpoint is truncated");
	const auto payload=bytes.substr(0,bytes.size()-checksum_bytes);checkpoint_metadata::Reader checksum(bytes.substr(bytes.size()-checksum_bytes));
	const auto expected=checksum.Text();checksum.Finish();Require(expected==native_tet_fsi_checkpoint_detail::Digest(payload),
		"native FSI checkpoint checksum mismatch");
	checkpoint_metadata::Reader input(payload);Require(input.Text()=="IGA_NATIVE_TET_FSI"&&input.Unsigned()==1,
		"unsupported native FSI checkpoint schema");
	NativeTetFsiCheckpoint value;value.accepted_steps=input.Unsigned();value.time_s=input.Real();
	value.fluid_model_identity_sha256=input.Text();value.solid_model_identity_sha256=input.Text();
	value.fluid_state_identity_sha256=input.Text();value.solid_state_identity_sha256=input.Text();
	value.fluid_state=native_tet_fsi_checkpoint_detail::ReadScalars(input);
	value.ale_displacement_m=native_tet_fsi_checkpoint_detail::ReadVectors(input);
	value.interface_displacement_m=native_tet_fsi_checkpoint_detail::ReadVectors(input);
	value.solid_displacement_m=native_tet_fsi_checkpoint_detail::ReadScalars(input);
	value.solid_velocity_m_per_s=native_tet_fsi_checkpoint_detail::ReadScalars(input);
	input.Finish();native_tet_fsi_checkpoint_detail::Validate(value);return value;
}

inline void ValidateRestoredNativeTetFsiCheckpoint(const NativeTetFsiCheckpoint& value,
	const NativeTetAleFsiRuntime& fluid,const NativeTetSolidFsiRuntime& solid)
{
	native_tet_fsi_checkpoint_detail::Validate(value);
	if(fluid.Lifecycle().HasActiveStep()||solid.Lifecycle().HasActiveStep()
		||value.time_s!=fluid.CommittedAleTimeS()
		||value.fluid_model_identity_sha256!=fluid.ModelIdentitySha256()
		||value.solid_model_identity_sha256!=solid.ModelIdentitySha256()
		||value.fluid_state_identity_sha256!=fluid.CommittedCompositionIdentitySha256()
		||value.solid_state_identity_sha256!=solid.CommittedStateIdentitySha256())
		throw std::runtime_error("restored native tetrahedral FSI pair differs from its checkpoint");
}

} // namespace iga
#endif

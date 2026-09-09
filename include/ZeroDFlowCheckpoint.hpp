#ifndef IGA_ZERO_D_FLOW_CHECKPOINT_HPP
#define IGA_ZERO_D_FLOW_CHECKPOINT_HPP

#include "CheckpointMetadataCodec.hpp"
#include "CoupledCheckpointManifest.hpp"
#include "ZeroDFlowDomain.hpp"

namespace iga {

inline void ValidateZeroDCheckpointIdentity(const ZeroDFlowCheckpointState& value)
{
	checkpoint_metadata::Require(!value.domain_id.empty() && value.model_identity_sha256.size() == 64
		&& value.model_identity_sha256.find_first_not_of("0123456789abcdef") == std::string::npos,
		"invalid 0D domain or model identity");
}

inline std::string SerializeZeroDFlowCheckpoint(const ZeroDFlowCheckpointState& value)
{
	using namespace checkpoint_metadata; ValidateZeroDCheckpointIdentity(value);
	ValidateZeroDFlowState(value.state); ValidateZeroDFlowStepAccounting(value.accounting);
	Writer output; output.Text("IGA_ZERO_D_ACCEPTED"); output.Unsigned(1);
	output.Text(value.domain_id); output.Text(value.model_identity_sha256); WriteStep(output, value.accepted_step);
	output.Real(value.state.stored_pressure_pa); WritePort(output, value.port);
	const auto& accounting = value.accounting;
	for (double amount : {accounting.initial_stored_volume_m3, accounting.final_stored_volume_m3,
		accounting.prescribed_source_amount_m3, accounting.distal_sink_amount_m3,
		accounting.outward_graph_port_amount_m3, accounting.residual_m3}) output.Real(amount);
	return output.Bytes();
}

inline ZeroDFlowCheckpointState ParseZeroDFlowCheckpoint(std::string_view bytes)
{
	using namespace checkpoint_metadata; Reader input(bytes);
	Require(input.Text() == "IGA_ZERO_D_ACCEPTED" && input.Unsigned() == 1, "unsupported 0D checkpoint schema");
	ZeroDFlowCheckpointState result; result.domain_id = input.Text(); result.model_identity_sha256 = input.Text();
	result.accepted_step = ReadStep(input); result.state.stored_pressure_pa = input.Real(); result.port = ReadPort(input);
	auto& accounting = result.accounting;
	accounting.initial_stored_volume_m3 = input.Real(); accounting.final_stored_volume_m3 = input.Real();
	accounting.prescribed_source_amount_m3 = input.Real(); accounting.distal_sink_amount_m3 = input.Real();
	accounting.outward_graph_port_amount_m3 = input.Real(); accounting.residual_m3 = input.Real();
	input.Finish(); ValidateZeroDCheckpointIdentity(result); return result;
}

inline ZeroDFlowCheckpointState ParseZeroDFlowCheckpoint(std::string_view bytes, const CoupledCheckpointEpoch& epoch)
{
	coupled_checkpoint_detail::Validate(epoch);
	auto result = ParseZeroDFlowCheckpoint(bytes);
	checkpoint_metadata::RequireAcceptedStep(result.accepted_step, epoch.accepted_steps, epoch.time_s, epoch.dt_s);
	return result;
}

} // namespace iga
#endif

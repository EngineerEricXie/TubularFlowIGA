#ifndef IGA_BODY_FITTED_ACCEPTED_CHECKPOINT_HPP
#define IGA_BODY_FITTED_ACCEPTED_CHECKPOINT_HPP

#include "TransientFlowRuntime.hpp"
#include "TransientTransportRuntime.hpp"
#include "CheckpointMetadataCodec.hpp"
#include "OwnedCheckpointFieldStream.hpp"
#include "CoupledCheckpointManifest.hpp"

namespace iga {
namespace body_fitted_checkpoint_detail {
using checkpoint_metadata::Require;
inline void Count(checkpoint_metadata::Reader& input, std::uint64_t expected)
{
	Require(input.Unsigned() == expected, "body-fitted checkpoint shape or catalog differs");
}
inline int PositiveCount(checkpoint_metadata::Reader& input)
{
	const auto count = input.Unsigned();
	Require(count > 0 && count <= static_cast<std::uint64_t>(INT_MAX), "invalid body-fitted accepted count");
	return static_cast<int>(count);
}
inline void Epoch(int count, double time, double dt, const CoupledCheckpointEpoch& epoch)
{
	Require(static_cast<std::uint64_t>(count) == epoch.accepted_steps && time == epoch.time_s && dt == epoch.dt_s,
		"body-fitted checkpoint clock differs from bundle");
}
} // namespace body_fitted_checkpoint_detail

// Shared metadata and replicated boundary arrays are stored ONCE per domain;
// each rank writes its own field shard. No full distributed field is gathered.
inline std::string SerializeBodyFittedFlowMetadata(const FlowAcceptedCheckpointState& state)
{
	using namespace body_fitted_checkpoint_detail;
	ValidateCheckpointConfigurationIdentity(state.configuration_identity_sha256);
	Require(state.accepted_steps > 0 && state.total_linear_iterations >= 0, "invalid flow counters");
	checkpoint_metadata::Writer output; output.Text("IGA_BODY_FITTED_FLOW/1"); output.Text(state.configuration_identity_sha256);
	output.Unsigned(state.accepted_steps); output.Real(state.accepted_time_s); output.Real(state.macro_dt_s);
	output.Unsigned(static_cast<std::uint64_t>(state.total_linear_iterations)); output.Unsigned(state.field.global_rows);
	output.Unsigned(state.boundaries.velocity.size()); output.Unsigned(state.boundaries.pressure.size());
	Require(state.pressure_tractions.size() <= checkpoint_metadata::maximum_map_entries, "too many pressure tractions");
	output.Unsigned(state.pressure_tractions.size());
	for (const auto& traction : state.pressure_tractions) {
		Require(traction.first >= 0, "negative traction label"); output.Unsigned(traction.first); output.Real(traction.second);
	}
	Require(state.outlets.size() <= checkpoint_metadata::maximum_map_entries, "too many outlets");
	output.Unsigned(state.outlets.size());
	for (const auto& outlet : state.outlets) {
		Require(outlet.label >= 0, "negative outlet label"); output.Unsigned(outlet.label); output.Unsigned(static_cast<std::uint64_t>(outlet.kind));
		for (double scalar : {outlet.resistance, outlet.proximal_resistance, outlet.distal_resistance, outlet.capacitance,
			outlet.reference_pressure, outlet.capacitor_pressure, outlet.flow, outlet.pressure}) output.Real(scalar);
	}
	return output.Bytes();
}

// Candidate comes from CreateCheckpointRestoreCandidate; shapes and immutable
// constraint masks, scalar boundary data and outlet models are already bound to
// the verified target configuration. Do not publish until every stream succeeds.
inline void ParseBodyFittedFlowMetadata(std::string_view bytes, FlowAcceptedCheckpointState& candidate,
	const CoupledCheckpointEpoch& epoch)
{
	using namespace body_fitted_checkpoint_detail;
	checkpoint_metadata::Reader input(bytes);
	Require(candidate.accepted_steps == 0 && input.Text() == "IGA_BODY_FITTED_FLOW/1", "invalid flow metadata target or version");
	ValidateCheckpointConfigurationIdentity(candidate.configuration_identity_sha256);
	Require(input.Text() == candidate.configuration_identity_sha256, "flow configuration identity differs");
	const int steps = PositiveCount(input); const double time = input.Real(), dt = input.Real();
	Epoch(steps, time, dt, epoch); Require(dt == candidate.macro_dt_s, "flow macro dt differs from target");
	const auto iterations = input.Unsigned(); Require(iterations <= static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max()), "flow iteration count overflows");
	Count(input, candidate.field.global_rows); Count(input, candidate.boundaries.velocity.size()); Count(input, candidate.boundaries.pressure.size());
	const auto count = input.Unsigned(); Require(count <= checkpoint_metadata::maximum_map_entries, "traction count exceeds limit");
	std::map<int, double> tractions; int previous = -1;
	for (std::uint64_t i = 0; i < count; ++i) {
		const auto label = input.Unsigned(); Require(label <= INT_MAX && static_cast<int>(label) > previous, "traction labels must be sorted and unique");
		previous = static_cast<int>(label); tractions.emplace(previous, input.Real());
	}
	Count(input, candidate.outlets.size());
	for (auto& outlet : candidate.outlets) {
		Count(input, static_cast<std::uint64_t>(outlet.label)); Count(input, static_cast<std::uint64_t>(outlet.kind));
		for (double scalar : {outlet.resistance, outlet.proximal_resistance, outlet.distal_resistance, outlet.capacitance, outlet.reference_pressure})
			Require(input.Real() == scalar, "outlet model differs from target");
		outlet.capacitor_pressure = input.Real(); outlet.flow = input.Real(); outlet.pressure = input.Real();
	}
	input.Finish(); candidate.pressure_tractions.swap(tractions); candidate.accepted_steps = steps; candidate.accepted_time_s = time;
	candidate.total_linear_iterations = static_cast<PetscInt>(iterations);
}

inline std::string SerializeBodyFittedTransportMetadata(const TransportAcceptedCheckpointState& state, double dt)
{
	body_fitted_checkpoint_detail::Require(state.accepted_steps > 0 && std::isfinite(dt) && dt > 0.0, "invalid transport clock");
	ValidateCheckpointConfigurationIdentity(state.configuration_identity_sha256);
	checkpoint_metadata::Writer output; output.Text("IGA_BODY_FITTED_TRANSPORT/1"); output.Text(state.configuration_identity_sha256);
	output.Unsigned(state.accepted_steps); output.Real(dt); output.Unsigned(state.field.global_rows); return output.Bytes();
}
inline void ParseBodyFittedTransportMetadata(std::string_view bytes, TransportAcceptedCheckpointState& candidate,
	double target_dt, const CoupledCheckpointEpoch& epoch)
{
	using namespace body_fitted_checkpoint_detail;
	checkpoint_metadata::Reader input(bytes);
	Require(candidate.accepted_steps == 0 && input.Text() == "IGA_BODY_FITTED_TRANSPORT/1", "invalid transport target or version");
	ValidateCheckpointConfigurationIdentity(candidate.configuration_identity_sha256);
	Require(input.Text() == candidate.configuration_identity_sha256, "transport configuration identity differs");
	const int count = PositiveCount(input); const double dt = input.Real();
	Require(dt == target_dt && dt == epoch.dt_s && static_cast<std::uint64_t>(count) == epoch.accepted_steps, "transport clock differs from bundle or target");
	const double expected_time = count*dt;
	Require(dt > 0.0 && std::isfinite(expected_time) && std::isfinite(epoch.time_s)
		&& std::abs(expected_time-epoch.time_s) <= 1e-12*std::max({1.0, expected_time, epoch.time_s}), "transport bundle time differs from count");
	Count(input, candidate.field.global_rows); input.Finish(); candidate.accepted_steps = count;
}

inline std::uint64_t BodyFittedBoundaryWords(const FlowAcceptedCheckpointState& state)
{
	checkpoint_stream::Require(state.boundaries.velocity.size() <= UINT64_MAX/3, "boundary shape overflows");
	return checkpoint_stream::Add(2, checkpoint_stream::Add(state.boundaries.velocity.size()*3, state.boundaries.pressure.size()));
}
inline std::uint64_t BodyFittedBoundaryBytes(const FlowAcceptedCheckpointState& state)
{
	return checkpoint_stream::Bytes(BodyFittedBoundaryWords(state));
}
template<class Sink> void WriteBodyFittedBoundaries(const FlowAcceptedCheckpointState& state, Sink sink)
{
	checkpoint_stream::Writer<Sink> output(std::move(sink));
	output.Word(state.boundaries.velocity.size()); output.Word(state.boundaries.pressure.size());
	for (const auto& value : state.boundaries.velocity) for (double scalar : value) output.Real(scalar);
	for (double value : state.boundaries.pressure) output.Real(value);
	output.Finish();
}
inline auto ReadBodyFittedBoundaries(FlowAcceptedCheckpointState& state)
{
	return checkpoint_stream::Reader(BodyFittedBoundaryWords(state), [&state](std::uint64_t index, std::uint64_t word) {
		if (index == 0) checkpoint_stream::Require(word == state.boundaries.velocity.size(), "velocity boundary shape differs");
		else if (index == 1) checkpoint_stream::Require(word == state.boundaries.pressure.size(), "pressure boundary shape differs");
		else {
			const auto offset = index-2, velocities = state.boundaries.velocity.size()*3;
			const auto value = checkpoint_stream::Real(word);
			if (offset < velocities) state.boundaries.velocity.at(static_cast<std::size_t>(offset/3))[offset%3] = value;
			else state.boundaries.pressure.at(static_cast<std::size_t>(offset-velocities)) = value;
		}
	});
}

} // namespace iga
#endif

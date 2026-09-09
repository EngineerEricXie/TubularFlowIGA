#ifndef IGA_PRESSURE_FLOW_CHECKPOINT_CONTROLS_HPP
#define IGA_PRESSURE_FLOW_CHECKPOINT_CONTROLS_HPP

#include "CheckpointMetadataCodec.hpp"
#include "CoupledCheckpointManifest.hpp"
#include "PressureFlowComponentExecutor.hpp"

namespace iga {

struct PressureFlowCheckpointControls {
	DomainStepContext accepted_step;
	std::map<std::string, double> next_pressure_pa;
	std::map<std::pair<std::string, std::string>, SpeciesDonor> donors;

	DomainStepContext NextStep() const
	{
		accepted_step.Validate();
		checkpoint_metadata::Require(accepted_step.step_index < std::numeric_limits<int>::max(), "next graph step overflows");
		DomainStepContext result{accepted_step.step_index+1, accepted_step.EndTime(), accepted_step.dt_s}; result.Validate(); return result;
	}
};

inline void ValidatePressureFlowCheckpointControls(const PressureFlowCheckpointControls& value,
	const SimulationGraph& graph, const CoupledCheckpointEpoch& epoch)
{
	using checkpoint_metadata::Require; coupled_checkpoint_detail::Validate(epoch);
	checkpoint_metadata::RequireAcceptedStep(value.accepted_step, epoch.accepted_steps, epoch.time_s, epoch.dt_s);
	Require(value.next_pressure_pa.size() == graph.Edges().size(), "checkpoint pressure edge catalog differs");
	std::size_t donor_count = 0;
	for (const auto& edge : graph.Edges()) {
		const auto pressure = value.next_pressure_pa.find(edge.id);
		Require(pressure != value.next_pressure_pa.end() && std::isfinite(pressure->second), "missing or nonfinite next pressure");
		for (const auto& species : edge.species) {
			++donor_count; const auto donor = value.donors.find({edge.id, species});
			Require(donor != value.donors.end() && (donor->second == SpeciesDonor::First || donor->second == SpeciesDonor::Second),
				"missing or invalid checkpoint donor");
		}
	}
	Require(value.donors.size() == donor_count, "unknown checkpoint donor keys");
}

// Only call after the executor returns an accepted result. The next guess is
// the measured pressure, exactly as in the native runner's accepted loop.
inline PressureFlowCheckpointControls MakePressureFlowCheckpointControls(const DomainStepContext& accepted_step,
	const PressureFlowIterationState& last_iteration,
	std::map<std::pair<std::string, std::string>, SpeciesDonor> donors = {})
{
	accepted_step.Validate(); PressureFlowCheckpointControls result; result.accepted_step = accepted_step;
	for (const auto& edge : last_iteration.edges)
		checkpoint_metadata::Require(!edge.edge_id.empty() && std::isfinite(edge.measured_pressure_pa)
			&& result.next_pressure_pa.emplace(edge.edge_id, edge.measured_pressure_pa).second, "invalid accepted pressure edge");
	result.donors = std::move(donors); return result;
}

inline std::string SerializePressureFlowCheckpointControls(const PressureFlowCheckpointControls& value)
{
	using namespace checkpoint_metadata; Writer output;
	output.Text("IGA_PRESSURE_FLOW_CONTROLS"); output.Unsigned(1); WriteStep(output, value.accepted_step);
	output.Reals(value.next_pressure_pa);
	Require(value.donors.size() <= maximum_map_entries, "donor count exceeds limit"); output.Unsigned(value.donors.size());
	for (const auto& entry : value.donors) {
		Require(!entry.first.first.empty() && !entry.first.second.empty(), "empty donor key");
		Require(entry.second == SpeciesDonor::First || entry.second == SpeciesDonor::Second, "unknown donor enum");
		output.Text(entry.first.first); output.Text(entry.first.second); output.Unsigned(entry.second == SpeciesDonor::First ? 0 : 1);
	}
	return output.Bytes();
}

inline PressureFlowCheckpointControls ParsePressureFlowCheckpointControls(std::string_view bytes,
	const SimulationGraph& graph, const CoupledCheckpointEpoch& epoch)
{
	using namespace checkpoint_metadata; Reader input(bytes);
	Require(input.Text() == "IGA_PRESSURE_FLOW_CONTROLS" && input.Unsigned() == 1, "unsupported graph control schema");
	PressureFlowCheckpointControls result; result.accepted_step = ReadStep(input); result.next_pressure_pa = input.Reals();
	const auto count = input.Unsigned(); Require(count <= maximum_map_entries, "donor count exceeds limit");
	std::pair<std::string, std::string> previous;
	for (std::uint64_t i = 0; i < count; ++i) {
		const auto edge = input.Text(); const auto species = input.Text(); const auto donor = input.Unsigned();
		const auto key = std::make_pair(edge, species);
		Require(!edge.empty() && !species.empty() && previous < key && donor <= 1, "invalid or duplicate donor entry");
		result.donors.emplace(key, donor == 0 ? SpeciesDonor::First : SpeciesDonor::Second); previous = key;
	}
	input.Finish(); ValidatePressureFlowCheckpointControls(result, graph, epoch); return result;
}

} // namespace iga
#endif

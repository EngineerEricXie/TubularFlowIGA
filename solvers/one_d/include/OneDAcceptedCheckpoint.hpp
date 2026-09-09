#ifndef IGA_ONE_D_ACCEPTED_CHECKPOINT_HPP
#define IGA_ONE_D_ACCEPTED_CHECKPOINT_HPP

#include "OneDRuntime.hpp"
#include "CheckpointMetadataCodec.hpp"
#include "CoupledCheckpointManifest.hpp"

namespace iga {
namespace one_d_checkpoint_detail {
using checkpoint_metadata::Require;

// The traversal is the field-shard schema. Metadata carries every shape, and
// the decoder obtains storage only from the verified, initialized target model.
template<class State, class Visit>
void Fields(State& state, Visit&& visit)
{
	visit(state.flow.area); visit(state.flow.flow); visit(state.flow.pressure);
	visit(state.flow.node_pressure); visit(state.flow.segment_flow); visit(state.segment_radii_m);
	for (auto& transport : state.transports)
		for (auto& species : transport.species) visit(species.concentration);
}

inline int Integer(checkpoint_metadata::Reader& input)
{
	const auto value = input.Unsigned();
	Require(value <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()), "1d integer overflows");
	return static_cast<int>(value);
}
inline bool Boolean(checkpoint_metadata::Reader& input)
{
	const auto value = input.Unsigned(); Require(value <= 1, "invalid 1d boolean"); return value != 0;
}
// Presence determines whether the scalar participates in validation. Preserve
// inactive storage bits (including signed zero or NaN) without interpreting it.
inline void WriteInletScalar(checkpoint_metadata::Writer& output, bool present, double value)
{
	output.Unsigned(present);
	if (present) output.Real(value);
	else { std::uint64_t bits; std::memcpy(&bits, &value, 8); output.Unsigned(bits); }
}
inline double ReadInletScalar(checkpoint_metadata::Reader& input, bool& present)
{
	present = Boolean(input);
	if (present) return input.Real();
	const auto bits = input.Unsigned(); double value; std::memcpy(&value, &bits, 8); return value;
}
inline void WriteInlet(checkpoint_metadata::Writer& output, const VascularInletState& inlet)
{
	ValidateVascularInletState(inlet); output.Real(inlet.time_s);
	WriteInletScalar(output, inlet.has_flow, inlet.flow_m3_s);
	WriteInletScalar(output, inlet.has_pressure, inlet.pressure_pa);
	WriteInletScalar(output, inlet.has_temperature, inlet.temperature_c);
	WriteInletScalar(output, inlet.has_hematocrit, inlet.hematocrit_percent);
	output.Reals(inlet.species);
	Require(inlet.metadata.size() <= checkpoint_metadata::maximum_map_entries, "1d inlet metadata exceeds limit");
	output.Unsigned(inlet.metadata.size());
	for (const auto& entry : inlet.metadata) {
		output.Text(entry.first); output.Text(entry.second);
	}
}
inline VascularInletState ReadInlet(checkpoint_metadata::Reader& input)
{
	VascularInletState inlet; inlet.time_s = input.Real();
	inlet.flow_m3_s = ReadInletScalar(input, inlet.has_flow);
	inlet.pressure_pa = ReadInletScalar(input, inlet.has_pressure);
	inlet.temperature_c = ReadInletScalar(input, inlet.has_temperature);
	inlet.hematocrit_percent = ReadInletScalar(input, inlet.has_hematocrit);
	inlet.species = input.Reals(); const auto count = input.Unsigned();
	Require(count <= checkpoint_metadata::maximum_map_entries, "1d inlet metadata exceeds limit");
	std::string previous;
	for (std::uint64_t i = 0; i < count; ++i) {
		auto key = input.Text(); Require(i == 0 || previous < key, "1d metadata keys must be sorted and unique");
		auto value = input.Text(); previous = key; inlet.metadata.emplace(std::move(key), std::move(value));
	}
	ValidateVascularInletState(inlet); return inlet;
}
inline void WriteOutletMap(checkpoint_metadata::Writer& output, const std::map<int, double>& values)
{
	Require(values.size() <= checkpoint_metadata::maximum_map_entries, "1d outlet map exceeds limit");
	output.Unsigned(values.size());
	for (const auto& entry : values) {
		Require(entry.first >= 0, "negative 1d outlet index"); output.Unsigned(entry.first); output.Real(entry.second);
	}
}
inline std::map<int, double> ReadOutletMap(checkpoint_metadata::Reader& input, const OneDNetwork& network)
{
	const auto count = input.Unsigned();
	Require(count == network.outlet_nodes.size() && count <= checkpoint_metadata::maximum_map_entries, "1d outlet coverage differs");
	std::map<int, double> result; int previous = -1;
	for (std::uint64_t i = 0; i < count; ++i) {
		const int node = Integer(input);
		Require(node > previous && std::find(network.outlet_nodes.begin(), network.outlet_nodes.end(), node) != network.outlet_nodes.end(),
			"1d outlet indices must be sorted, unique and known");
		previous = node; result.emplace(node, input.Real());
	}
	return result;
}
} // namespace one_d_checkpoint_detail

// This small shard and the streamed fields shard must share one bundle epoch.
// Static constitutive data comes from the target whose SHA identity is checked
// here and by RestoreCheckpointState. This does not load legacy 1d checkpoints.
inline std::string SerializeOneDAcceptedCheckpointMetadata(const OneDFlowCheckpointState& value)
{
	using namespace one_d_checkpoint_detail;
	checkpoint_metadata::Writer output; output.Text("IGA_ONE_D_ACCEPTED/1");
	output.Text(value.configuration_identity_sha256);
	Require(value.accepted_macro_steps > 0 && value.flow.completed_step > 0 && value.flow.internal_substeps >= 0, "invalid 1d accepted counters");
	output.Unsigned(value.accepted_macro_steps); output.Real(value.last_macro_start_s); output.Real(value.last_macro_dt_s);
	output.Unsigned(value.flow.completed_step); output.Real(value.flow.physical_time);
	output.Unsigned(static_cast<std::uint64_t>(value.flow.internal_substeps)); output.Real(value.flow.inlet_flow);
	WriteInlet(output, value.last_inlet); for (double entry : value.blood_state) output.Real(entry);
	output.Unsigned(value.flow.outlets.size());
	for (const auto& outlet : value.flow.outlets) {
		Require(outlet.node >= 0, "negative 1d outlet index"); output.Unsigned(outlet.node);
		output.Real(outlet.pressure); output.Real(outlet.reference_pressure);
		output.Real(outlet.capacitor_pressure); output.Real(outlet.flow);
	}
	output.Unsigned(value.transports.size());
	for (const auto& transport : value.transports) {
		output.Text(transport.name); output.Unsigned(transport.species.size());
		for (const auto& species : transport.species) {
			output.Text(species.definition.field); output.Real(species.inlet_value); output.Text(species.inlet_waveform);
			output.Unsigned(species.boundary_flux_valid); output.Real(species.root_native_flux); WriteOutletMap(output, species.outlet_native_flux);
			output.Unsigned(species.step_accounting_valid); output.Real(species.step_initial_mass);
			output.Real(species.step_root_native_amount); output.Real(species.step_source_amount); WriteOutletMap(output, species.step_outlet_native_amount);
		}
	}
	Fields(value, [&](const auto& field) { output.Unsigned(field.size()); });
	return output.Bytes();
}

inline OneDFlowCheckpointState ParseOneDAcceptedCheckpointMetadata(std::string_view bytes,
	const OneDFlowRuntime& target, const CoupledCheckpointEpoch& epoch)
{
	using namespace one_d_checkpoint_detail;
	Require(target.CurrentPhase() == OneDFlowRuntime::Phase::Ready && target.FlowState().completed_step == 0,
		"1d metadata requires a fresh initialized target");
	checkpoint_metadata::Reader input(bytes);
	Require(input.Text() == "IGA_ONE_D_ACCEPTED/1", "unsupported 1d checkpoint version");
	OneDFlowCheckpointState value; value.configuration_identity_sha256 = input.Text();
	Require(!target.CheckpointConfigurationIdentity().empty() && value.configuration_identity_sha256 == target.CheckpointConfigurationIdentity(),
		"1d checkpoint configuration identity differs or is unbound");
	value.accepted_macro_steps = Integer(input); value.last_macro_start_s = input.Real(); value.last_macro_dt_s = input.Real();
	Require(value.accepted_macro_steps > 0, "1d checkpoint is not an accepted step");
	checkpoint_metadata::RequireAcceptedStep({value.accepted_macro_steps-1, value.last_macro_start_s, value.last_macro_dt_s},
		epoch.accepted_steps, epoch.time_s, epoch.dt_s);
	value.flow = target.FlowState(); value.transports = target.Transports();
	for (const auto& segment : target.Network().segments) value.segment_radii_m.push_back(segment.radius0);
	value.flow.completed_step = Integer(input); value.flow.physical_time = input.Real();
	const auto internal = input.Unsigned(); Require(internal <= static_cast<std::uint64_t>(std::numeric_limits<long long>::max()), "1d internal substeps overflow");
	value.flow.internal_substeps = static_cast<long long>(internal); value.flow.inlet_flow = input.Real();
	value.last_inlet = ReadInlet(input); for (auto& entry : value.blood_state) entry = input.Real();
	Require(input.Unsigned() == value.flow.outlets.size(), "1d outlet catalog differs");
	for (auto& outlet : value.flow.outlets) {
		Require(Integer(input) == outlet.node, "1d outlet order differs");
		outlet.pressure = input.Real(); outlet.reference_pressure = input.Real();
		outlet.capacitor_pressure = input.Real(); outlet.flow = input.Real();
	}
	Require(input.Unsigned() == value.transports.size(), "1d transport catalog differs");
	for (auto& transport : value.transports) {
		Require(input.Text() == transport.name && input.Unsigned() == transport.species.size(), "1d species catalog differs");
		for (auto& species : transport.species) {
			Require(input.Text() == species.definition.field, "1d species order differs");
			species.inlet_value = input.Real(); species.inlet_waveform = input.Text();
			species.boundary_flux_valid = Boolean(input); species.root_native_flux = input.Real(); species.outlet_native_flux = ReadOutletMap(input, target.Network());
			species.step_accounting_valid = Boolean(input); species.step_initial_mass = input.Real();
			species.step_root_native_amount = input.Real(); species.step_source_amount = input.Real(); species.step_outlet_native_amount = ReadOutletMap(input, target.Network());
		}
	}
	Fields(value, [&](const auto& field) { Require(input.Unsigned() == field.size(), "1d checkpoint field shape differs"); });
	input.Finish(); return value;
}

inline std::uint64_t OneDAcceptedCheckpointFieldBytes(const OneDFlowCheckpointState& value)
{
	std::uint64_t bytes = 0;
	one_d_checkpoint_detail::Fields(value, [&](const auto& field) {
		checkpoint_metadata::Require(field.size() <= (UINT64_MAX-bytes)/8, "1d field payload overflows"); bytes += field.size()*8;
	});
	return bytes;
}

// A producer-compatible sink receives at most 64 KiB per call. No full encoded
// field image is allocated; each binary64 value uses canonical little endian.
template<class Sink>
void WriteOneDAcceptedCheckpointFields(const OneDFlowCheckpointState& value, Sink&& sink)
{
	std::array<char, 65536> buffer; std::size_t used = 0;
	one_d_checkpoint_detail::Fields(value, [&](const auto& field) {
		for (double entry : field) {
			checkpoint_metadata::Require(std::isfinite(entry), "nonfinite 1d field");
			std::uint64_t bits; std::memcpy(&bits, &entry, 8);
			for (unsigned byte = 0; byte < 8; ++byte) buffer[used++] = static_cast<char>((bits>>(8*byte)) & 255);
			if (used == buffer.size()) { sink(buffer.data(), used); used = 0; }
		}
	});
	if (used) sink(buffer.data(), used);
}

// Own the candidate outside this consumer, keep it unpublished, and do not
// resize its vectors until Finish AND bundle SHA verification have succeeded.
// Failure leaves only the candidate dirty; no live runtime is modified.
class OneDAcceptedCheckpointFieldsReader {
public:
	explicit OneDAcceptedCheckpointFieldsReader(OneDFlowCheckpointState& candidate)
		: expected_(OneDAcceptedCheckpointFieldBytes(candidate))
	{
		one_d_checkpoint_detail::Fields(candidate, [&](auto& field) { if (!field.empty()) fields_.push_back(&field); });
	}
	OneDAcceptedCheckpointFieldsReader(const OneDAcceptedCheckpointFieldsReader&) = delete;
	OneDAcceptedCheckpointFieldsReader& operator=(const OneDAcceptedCheckpointFieldsReader&) = delete;
	void Consume(const void* data, std::size_t size)
	{
		checkpoint_metadata::Require(!failed_, "1d field reader already failed");
		try {
			checkpoint_metadata::Require(size <= expected_-consumed_, "trailing 1d field bytes");
			const auto* bytes = static_cast<const unsigned char*>(data);
			for (std::size_t i = 0; i < size; ++i) {
				bits_ |= std::uint64_t(bytes[i])<<(8*partial_); ++partial_; ++consumed_;
				if (partial_ == 8) {
					double value; std::memcpy(&value, &bits_, 8);
					checkpoint_metadata::Require(std::isfinite(value), "nonfinite 1d field payload");
					auto& field = *fields_.at(field_); field[element_++] = value;
					if (element_ == field.size()) { ++field_; element_ = 0; }
					partial_ = 0; bits_ = 0;
				}
			}
		} catch (...) { failed_ = true; throw; }
	}
	void Finish() const
	{
		checkpoint_metadata::Require(!failed_ && consumed_ == expected_ && partial_ == 0 && field_ == fields_.size(), "incomplete 1d field payload");
	}
private:
	std::vector<std::vector<double>*> fields_;
	std::uint64_t expected_ = 0, consumed_ = 0, bits_ = 0;
	std::size_t field_ = 0, element_ = 0;
	unsigned partial_ = 0;
	bool failed_ = false;
};

} // namespace iga
#endif

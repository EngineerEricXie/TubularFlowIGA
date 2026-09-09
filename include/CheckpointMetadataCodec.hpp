#ifndef IGA_CHECKPOINT_METADATA_CODEC_HPP
#define IGA_CHECKPOINT_METADATA_CODEC_HPP

// Bounded metadata codec, separate from streamed field payloads. Integers and
// IEEE binary64 bits use little-endian encoding; strings are length-prefixed.
#include "CoupledDomainRuntime.hpp"
#include <cstdint>
#include <cstring>
#include <limits>
#include <string_view>

namespace iga {
namespace checkpoint_metadata {

constexpr std::size_t maximum_bytes = 16*1024*1024;
constexpr std::size_t maximum_string_bytes = 4096;
constexpr std::size_t maximum_map_entries = 16384;

inline void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(std::string("checkpoint metadata: ")+message);
}

class Writer {
public:
	void Unsigned(std::uint64_t value)
	{
		char bytes[8]; for (unsigned i = 0; i < 8; ++i) bytes[i] = static_cast<char>((value>>(8*i)) & 255);
		Append(bytes, sizeof(bytes));
	}
	void Real(double value)
	{
		static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "checkpoint requires IEEE binary64");
		Require(std::isfinite(value), "nonfinite scalar");
		std::uint64_t bits; std::memcpy(&bits, &value, sizeof(bits)); Unsigned(bits);
	}
	void Text(const std::string& value)
	{
		Require(value.size() <= maximum_string_bytes, "string exceeds limit");
		Unsigned(value.size()); Append(value.data(), value.size());
	}
	void Optional(const std::optional<double>& value)
	{
		Unsigned(value ? 1 : 0); if (value) Real(*value);
	}
	void Reals(const std::map<std::string, double>& values)
	{
		Require(values.size() <= maximum_map_entries, "map exceeds limit"); Unsigned(values.size());
		for (const auto& value : values) { Require(!value.first.empty(), "empty map key"); Text(value.first); Real(value.second); }
	}
	const std::string& Bytes() const noexcept { return bytes_; }
private:
	void Append(const char* data, std::size_t count)
	{
		Require(count <= maximum_bytes-bytes_.size(), "payload exceeds metadata limit"); bytes_.append(data, count);
	}
	std::string bytes_;
};

class Reader {
public:
	explicit Reader(std::string_view bytes) : bytes_(bytes)
	{ Require(bytes.size() <= maximum_bytes, "payload exceeds metadata limit"); }
	std::uint64_t Unsigned()
	{
		Require(bytes_.size()-position_ >= 8, "truncated integer"); std::uint64_t value = 0;
		for (unsigned i = 0; i < 8; ++i) value |= std::uint64_t(static_cast<unsigned char>(bytes_[position_++]))<<(8*i);
		return value;
	}
	double Real()
	{
		const auto bits = Unsigned(); double value; std::memcpy(&value, &bits, sizeof(value));
		Require(std::isfinite(value), "nonfinite scalar"); return value;
	}
	std::string Text()
	{
		const auto size = Unsigned(); Require(size <= maximum_string_bytes && size <= bytes_.size()-position_, "invalid string length");
		std::string value(bytes_.substr(position_, static_cast<std::size_t>(size))); position_ += static_cast<std::size_t>(size); return value;
	}
	std::optional<double> Optional()
	{
		const auto present = Unsigned(); Require(present <= 1, "invalid optional presence");
		if (!present) return std::nullopt;
		return Real();
	}
	std::map<std::string, double> Reals()
	{
		const auto count = Unsigned(); Require(count <= maximum_map_entries, "map exceeds limit");
		std::map<std::string, double> result; std::string previous;
		for (std::uint64_t i = 0; i < count; ++i) {
			auto key = Text(); Require(previous < key, "map keys must be nonempty, sorted and unique");
			const auto value = Real(); previous = key; result.emplace(std::move(key), value);
		}
		return result;
	}
	void Finish() const { Require(position_ == bytes_.size(), "trailing payload bytes"); }
private:
	std::string_view bytes_;
	std::size_t position_ = 0;
};

inline void WritePort(Writer& output, const PortState& port)
{
	ValidatePortState(port); output.Real(port.time_s); output.Optional(port.area_m2);
	output.Optional(port.outward_flow_m3_s); output.Optional(port.mean_pressure_pa);
	output.Optional(port.mean_normal_traction_pa); output.Optional(port.total_pressure_pa);
	output.Reals(port.concentration); output.Reals(port.outward_species_flux);
}

inline PortState ReadPort(Reader& input)
{
	PortState port; port.time_s = input.Real(); port.area_m2 = input.Optional();
	port.outward_flow_m3_s = input.Optional(); port.mean_pressure_pa = input.Optional();
	port.mean_normal_traction_pa = input.Optional(); port.total_pressure_pa = input.Optional();
	port.concentration = input.Reals(); port.outward_species_flux = input.Reals(); ValidatePortState(port); return port;
}

inline void WriteStep(Writer& output, const DomainStepContext& step)
{
	step.Validate(); output.Unsigned(static_cast<std::uint64_t>(step.step_index)); output.Real(step.start_time_s); output.Real(step.dt_s);
}

inline DomainStepContext ReadStep(Reader& input)
{
	const auto index = input.Unsigned(); Require(index <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()), "step index overflows");
	DomainStepContext step; step.step_index = static_cast<int>(index); step.start_time_s = input.Real(); step.dt_s = input.Real(); step.Validate(); return step;
}

inline void RequireAcceptedStep(const DomainStepContext& step, std::uint64_t count, double time, double dt)
{
	step.Validate();
	Require(count == static_cast<std::uint64_t>(step.step_index)+1 && time == step.EndTime() && dt == step.dt_s,
		"accepted step differs from bundle clock");
}

} // namespace checkpoint_metadata
} // namespace iga
#endif

#ifndef IGA_COUPLED_CHECKPOINT_MANIFEST_HPP
#define IGA_COUPLED_CHECKPOINT_MANIFEST_HPP

#include "Sha256.hpp"
#include <charconv>
#include <climits>
#include <cstring>
#include <locale>
#include <sstream>
#include <vector>

namespace iga {

struct CoupledCheckpointCompatibility {
	std::string case_sha256, configuration_sha256, execution_sha256;
	std::uint32_t ranks = 0;
};

struct CoupledCheckpointEpoch {
	std::string id, previous_id;
	CoupledCheckpointCompatibility compatibility;
	std::uint64_t accepted_steps = 0;
	double time_s = 0, dt_s = 0;
};

struct CoupledCheckpointShardSpec {
	// Portable logical IDs, not caller-provided paths. Catalog sorted by id.
	std::string id, format;
	std::uint32_t owner_rank = 0;
};

struct CoupledCheckpointShard {
	CoupledCheckpointShardSpec spec;
	std::uint64_t payload_bytes = 0;
	std::string sha256; // Entire envelope plus payload, binding the epoch.
};

struct CoupledCheckpointManifest {
	CoupledCheckpointEpoch epoch;
	std::vector<CoupledCheckpointShard> shards;
};

namespace coupled_checkpoint_detail {

constexpr std::size_t maximum_manifest_bytes = 4*1024*1024;
constexpr std::size_t maximum_shards = 4096;
constexpr std::uint64_t maximum_payload_bytes = std::uint64_t{1}<<40;

inline void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(std::string("coupled checkpoint: ")+message);
}

inline bool Digest(const std::string& value)
{
	return value.size() == 64 && value.find_first_not_of("0123456789abcdef") == std::string::npos;
}

inline bool Identifier(const std::string& value)
{
	return !value.empty() && value.size() <= 128 && value != "." && value != ".."
		&& value.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-") == std::string::npos;
}

inline void Validate(const CoupledCheckpointCompatibility& value)
{
	Require(Digest(value.case_sha256) && Digest(value.configuration_sha256)
		&& Digest(value.execution_sha256), "invalid compatibility digest");
	Require(value.ranks > 0 && value.ranks <= INT_MAX, "invalid rank count");
}

inline void Validate(const CoupledCheckpointEpoch& value)
{
	Validate(value.compatibility);
	Require(Digest(value.id) && (value.previous_id.empty() || Digest(value.previous_id))
		&& value.previous_id != value.id, "invalid epoch identity");
	Require(value.accepted_steps > 0 && std::isfinite(value.time_s) && value.time_s > 0
		&& std::isfinite(value.dt_s) && value.dt_s > 0, "invalid accepted clock");
}

inline void Validate(const CoupledCheckpointShardSpec& value, std::uint32_t ranks)
{
	Require(Identifier(value.id) && Identifier(value.format), "invalid shard identifier or format");
	Require(value.owner_rank < ranks, "shard owner outside communicator");
}

inline void ValidateCatalog(const std::vector<CoupledCheckpointShardSpec>& specs, std::uint32_t ranks)
{
	Require(!specs.empty() && specs.size() <= maximum_shards, "invalid shard count");
	std::string previous;
	for (const auto& spec : specs) {
		Validate(spec, ranks);
		Require(previous < spec.id, "shard catalog must have sorted unique IDs");
		previous = spec.id;
	}
}

inline bool Same(const CoupledCheckpointShardSpec& a, const CoupledCheckpointShardSpec& b)
{
	return a.id == b.id && a.format == b.format && a.owner_rank == b.owner_rank;
}

inline void RequireCompatible(const CoupledCheckpointCompatibility& actual,
	const CoupledCheckpointCompatibility& expected)
{
	Validate(actual); Validate(expected);
	Require(actual.case_sha256 == expected.case_sha256
		&& actual.configuration_sha256 == expected.configuration_sha256
		&& actual.execution_sha256 == expected.execution_sha256 && actual.ranks == expected.ranks,
		"incompatible case, configuration, execution, or ranks");
}

inline std::string DoubleBits(double value)
{
	static_assert(sizeof(double) == sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559,
		"checkpoint requires IEEE binary64");
	std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits));
	std::string result(16, '0');
	for (int i = 15; i >= 0; --i) { result[i] = "0123456789abcdef"[bits & 15]; bits >>= 4; }
	return result;
}

inline std::uint64_t Unsigned(const std::string& text, int base = 10)
{
	std::uint64_t result = 0;
	const auto parsed = std::from_chars(text.data(), text.data()+text.size(), result, base);
	Require(!text.empty() && parsed.ec == std::errc{} && parsed.ptr == text.data()+text.size(),
		"invalid or overflowing integer");
	return result;
}

inline double ReadDoubleBits(const std::string& text)
{
	Require(text.size() == 16 && text.find_first_not_of("0123456789abcdef") == std::string::npos,
		"invalid binary64 encoding");
	const auto bits = Unsigned(text, 16); double result = 0;
	std::memcpy(&result, &bits, sizeof(bits)); return result;
}

inline std::string Hash(const std::string& text)
{
	Sha256 hash; hash.Append(text.data(), text.size()); return hash.Hex();
}

inline std::string ManifestBody(const CoupledCheckpointManifest& manifest)
{
	const auto& epoch = manifest.epoch; Validate(epoch);
	Require(!manifest.shards.empty() && manifest.shards.size() <= maximum_shards, "invalid shard count");
	std::vector<CoupledCheckpointShardSpec> specs;
	for (const auto& shard : manifest.shards) specs.push_back(shard.spec);
	ValidateCatalog(specs, epoch.compatibility.ranks);
	std::ostringstream out; out.imbue(std::locale::classic());
	out.exceptions(std::ios::badbit | std::ios::failbit);
	out << "IGA_COUPLED_CHECKPOINT 1\nepoch " << epoch.id << "\nprevious "
		<< (epoch.previous_id.empty() ? "-" : epoch.previous_id)
		<< "\ncase " << epoch.compatibility.case_sha256
		<< "\nconfiguration " << epoch.compatibility.configuration_sha256
		<< "\nexecution " << epoch.compatibility.execution_sha256
		<< "\nranks " << epoch.compatibility.ranks << "\naccepted_steps " << epoch.accepted_steps
		<< "\ntime_bits " << DoubleBits(epoch.time_s) << "\ndt_bits " << DoubleBits(epoch.dt_s)
		<< "\nshards " << manifest.shards.size() << '\n';
	for (const auto& shard : manifest.shards) {
		Require(shard.payload_bytes <= maximum_payload_bytes && Digest(shard.sha256), "invalid shard size or checksum");
		out << "shard " << shard.spec.id << ' ' << shard.spec.owner_rank << ' ' << shard.spec.format
			<< ' ' << shard.payload_bytes << ' ' << shard.sha256 << '\n';
	}
	return out.str();
}

} // namespace coupled_checkpoint_detail

inline std::string SerializeCoupledCheckpointManifest(const CoupledCheckpointManifest& manifest)
{
	using namespace coupled_checkpoint_detail;
	const auto body = ManifestBody(manifest);
	const auto text = body+"sha256 "+Hash(body)+"\n";
	Require(text.size() <= maximum_manifest_bytes, "manifest exceeds size limit"); return text;
}

inline CoupledCheckpointManifest ParseCoupledCheckpointManifest(const std::string& text)
{
	using namespace coupled_checkpoint_detail;
	Require(text.size() <= maximum_manifest_bytes, "manifest exceeds size limit");
	std::istringstream input(text); input.imbue(std::locale::classic());
	auto token = [&]() {
		std::string value; Require(static_cast<bool>(input >> value) && value.size() <= 128, "missing or oversized manifest token");
		return value;
	};
	auto expect = [&](const char* key) { Require(token() == key, "unexpected manifest key or version"); };
	auto rank = [&]() { const auto value = Unsigned(token()); Require(value <= INT_MAX, "rank overflows"); return static_cast<std::uint32_t>(value); };
	CoupledCheckpointManifest result; auto& epoch = result.epoch;
	expect("IGA_COUPLED_CHECKPOINT"); expect("1");
	expect("epoch"); epoch.id = token(); expect("previous"); epoch.previous_id = token();
	if (epoch.previous_id == "-") epoch.previous_id.clear();
	expect("case"); epoch.compatibility.case_sha256 = token();
	expect("configuration"); epoch.compatibility.configuration_sha256 = token();
	expect("execution"); epoch.compatibility.execution_sha256 = token();
	expect("ranks"); epoch.compatibility.ranks = rank();
	expect("accepted_steps"); epoch.accepted_steps = Unsigned(token());
	expect("time_bits"); epoch.time_s = ReadDoubleBits(token());
	expect("dt_bits"); epoch.dt_s = ReadDoubleBits(token());
	expect("shards"); const auto count = Unsigned(token());
	Require(count > 0 && count <= maximum_shards, "invalid shard count");
	for (std::uint64_t i = 0; i < count; ++i) {
		CoupledCheckpointShard shard; expect("shard"); shard.spec.id = token(); shard.spec.owner_rank = rank();
		shard.spec.format = token(); shard.payload_bytes = Unsigned(token()); shard.sha256 = token();
		result.shards.push_back(std::move(shard));
	}
	expect("sha256"); const auto checksum = token();
	Require(checksum == Hash(ManifestBody(result)), "manifest checksum mismatch");
	Require(SerializeCoupledCheckpointManifest(result) == text, "noncanonical or trailing manifest data");
	return result;
}

} // namespace iga
#endif

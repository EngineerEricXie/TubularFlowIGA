#ifndef IGA_COUPLED_CHECKPOINT_BUNDLE_HPP
#define IGA_COUPLED_CHECKPOINT_BUNDLE_HPP

// POSIX shared-storage protocol. No MPI calls: the coordinator must supply the
// agreed complete catalog, and invoke Publish after all writers have returned.
// Root must already exist durably. Epochs and published files are immutable.
#include "CoupledCheckpointManifest.hpp"
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <optional>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>

namespace iga {
namespace coupled_checkpoint_detail {

#ifdef IGA_COUPLED_CHECKPOINT_TESTING
inline bool fail_next_sync = false;
inline bool fail_next_write = false;
inline std::function<void(const CoupledCheckpointEpoch&, const CoupledCheckpointShardSpec&, int)> before_shard_payload;
inline std::function<void(const CoupledCheckpointEpoch&)> before_manifest_publication;
#endif

class Descriptor {
public:
	explicit Descriptor(int value) : value_(value)
	{
		if (value < 0) throw std::runtime_error("coupled checkpoint: open failed: "+std::string(std::strerror(errno)));
	}
	Descriptor(const Descriptor&) = delete;
	Descriptor& operator=(const Descriptor&) = delete;
	Descriptor(Descriptor&& other) noexcept : value_(other.value_) { other.value_ = -1; }
	~Descriptor() { if (value_ >= 0) ::close(value_); }
	int Get() const noexcept { return value_; }
	void Sync() const
	{
#ifdef IGA_COUPLED_CHECKPOINT_TESTING
		if (fail_next_sync) { fail_next_sync = false; throw std::runtime_error("injected checkpoint sync failure"); }
#endif
		int status; do { status = ::fsync(value_); } while (status && errno == EINTR);
		Require(status == 0, "fsync failed");
	}
	void Close()
	{
		const int value = value_; value_ = -1;
		Require(::close(value) == 0, "close failed");
	}
private:
	int value_;
};

inline Descriptor Root(const std::filesystem::path& root)
{
	const auto name = root.string();
	Require(!name.empty() && name.find('\0') == std::string::npos, "invalid root path");
	return Descriptor(::open(name.c_str(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC));
}

inline Descriptor Directory(int root, const std::string& id)
{
	Require(Digest(id), "invalid epoch directory");
	return Descriptor(::openat(root, id.c_str(), O_RDONLY | O_DIRECTORY | O_NOFOLLOW | O_CLOEXEC));
}

inline Descriptor Input(int directory, const std::string& name)
{
	Descriptor file(::openat(directory, name.c_str(), O_RDONLY | O_NONBLOCK | O_NOFOLLOW | O_CLOEXEC));
	struct stat info{};
	Require(::fstat(file.Get(), &info) == 0 && S_ISREG(info.st_mode), "input is not a regular file");
	return file;
}

inline Descriptor Output(int directory, const std::string& name)
{
	return Descriptor(::openat(directory, name.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC, 0600));
}

inline void Write(int fd, const void* data, std::size_t bytes)
{
#ifdef IGA_COUPLED_CHECKPOINT_TESTING
	if (fail_next_write) { fail_next_write = false; throw std::runtime_error("injected checkpoint write failure"); }
#endif
	Require(data || !bytes, "null payload");
	const auto* pointer = static_cast<const char*>(data);
	while (bytes) {
		const auto count = ::write(fd, pointer, std::min<std::size_t>(bytes, 65536));
		if (count < 0 && errno == EINTR) continue;
		Require(count > 0, "write failed"); pointer += count; bytes -= static_cast<std::size_t>(count);
	}
}

inline std::size_t Read(int fd, void* data, std::size_t bytes)
{
	ssize_t count; do { count = ::read(fd, data, bytes); } while (count < 0 && errno == EINTR);
	Require(count >= 0, "read failed"); return static_cast<std::size_t>(count);
}

inline void RequireUnpublished(int directory)
{
	struct stat info{};
	const int status = ::fstatat(directory, "manifest", &info, AT_SYMLINK_NOFOLLOW);
	Require(status != 0 && errno == ENOENT, "epoch is already published or inaccessible");
}

inline void PublishFile(const Descriptor& directory, const std::string& temporary, const std::string& final)
{
	// linkat gives atomic no-replace publication, unlike overwriting rename.
	Require(::linkat(directory.Get(), temporary.c_str(), directory.Get(), final.c_str(), 0) == 0,
		"cannot publish immutable file");
	Require(::unlinkat(directory.Get(), temporary.c_str(), 0) == 0, "cannot remove published temporary link");
	directory.Sync();
}

inline std::string ShardHeader(const std::string& epoch, const CoupledCheckpointShardSpec& spec, std::uint64_t bytes)
{
	return "IGA_CP_SHARD 1\nepoch "+epoch+"\nid "+spec.id+"\nowner "+std::to_string(spec.owner_rank)
		+"\nformat "+spec.format+"\nbytes "+std::to_string(bytes)+"\n";
}

using Consumer = std::function<void(const void*, std::size_t)>;

inline CoupledCheckpointShard Inspect(const Descriptor& directory, const CoupledCheckpointEpoch& epoch,
	const CoupledCheckpointShardSpec& spec, bool sync, const Consumer& consumer = {})
{
	Validate(spec, epoch.compatibility.ranks);
	auto file = Input(directory.Get(), spec.id+".shard");
	std::string header;
	for (int lines = 0; lines < 6;) {
		char value;
		Require(header.size() < 512 && Read(file.Get(), &value, 1) == 1, "truncated or oversized shard header");
		header += value; if (value == '\n') ++lines;
	}
	const auto position = header.rfind("\nbytes ");
	Require(position != std::string::npos, "invalid shard size header");
	const auto bytes = Unsigned(header.substr(position+7, header.size()-position-8));
	Require(bytes <= maximum_payload_bytes, "shard exceeds size limit");
	Require(header == ShardHeader(epoch.id, spec, bytes), "shard epoch, owner, ID, format, or encoding differs");
	struct stat info{};
	Require(::fstat(file.Get(), &info) == 0 && info.st_size >= 0
		&& static_cast<std::uint64_t>(info.st_size) == header.size()+bytes, "truncated or oversized shard payload");
	Sha256 hash; hash.Append(header.data(), header.size());
	std::array<char, 65536> buffer{}; std::uint64_t remaining = bytes;
	while (remaining) {
		const auto count = Read(file.Get(), buffer.data(), static_cast<std::size_t>(std::min<std::uint64_t>(remaining, buffer.size())));
		Require(count > 0, "truncated shard while reading");
		hash.Append(buffer.data(), count); if (consumer) consumer(buffer.data(), count); remaining -= count;
	}
	char extra; Require(Read(file.Get(), &extra, 1) == 0, "shard grew while reading");
	if (sync) file.Sync();
	file.Close(); return {spec, bytes, hash.Hex()};
}

inline std::string ReadManifest(const Descriptor& directory)
{
	auto file = Input(directory.Get(), "manifest"); struct stat info{};
	Require(::fstat(file.Get(), &info) == 0 && info.st_size >= 0
		&& static_cast<std::uint64_t>(info.st_size) <= maximum_manifest_bytes, "manifest exceeds size limit");
	std::string result(static_cast<std::size_t>(info.st_size), '\0'); std::size_t done = 0;
	while (done < result.size()) {
		const auto count = Read(file.Get(), result.data()+done, result.size()-done);
		Require(count > 0, "truncated manifest"); done += count;
	}
	char extra; Require(Read(file.Get(), &extra, 1) == 0, "manifest grew while reading"); file.Close(); return result;
}

} // namespace coupled_checkpoint_detail

class CoupledCheckpointShardOutput {
public:
	CoupledCheckpointShardOutput(int fd, std::uint64_t bytes, const std::string& header) : fd_(fd), remaining_(bytes)
	{ hash_.Append(header.data(), header.size()); }
	void Write(const void* data, std::size_t bytes)
	{
		coupled_checkpoint_detail::Require(bytes <= remaining_, "writer exceeded declared payload size");
		coupled_checkpoint_detail::Write(fd_, data, bytes); hash_.Append(data, bytes); remaining_ -= bytes;
	}
	std::uint64_t Remaining() const noexcept { return remaining_; }
	std::string Sha256Hex() const { return hash_.Hex(); }
private:
	int fd_;
	std::uint64_t remaining_;
	Sha256 hash_;
};

inline void CreateCoupledCheckpointEpoch(const std::filesystem::path& root, const CoupledCheckpointEpoch& epoch)
{
	using namespace coupled_checkpoint_detail; Validate(epoch); auto parent = Root(root);
	Require(::mkdirat(parent.Get(), epoch.id.c_str(), 0700) == 0, "epoch already exists or cannot be created");
	parent.Sync();
}

template<class Writer>
CoupledCheckpointShard WriteCoupledCheckpointShard(const std::filesystem::path& root, const CoupledCheckpointEpoch& epoch,
	const CoupledCheckpointShardSpec& spec, std::uint64_t bytes, Writer&& writer)
{
	using namespace coupled_checkpoint_detail; Validate(epoch); Validate(spec, epoch.compatibility.ranks);
	Require(bytes <= maximum_payload_bytes, "shard exceeds size limit");
	auto parent = Root(root); auto directory = Directory(parent.Get(), epoch.id); RequireUnpublished(directory.Get());
	const auto temporary = spec.id+".tmp"; auto file = Output(directory.Get(), temporary);
	const auto header = ShardHeader(epoch.id, spec, bytes); Write(file.Get(), header.data(), header.size());
	CoupledCheckpointShardOutput output(file.Get(), bytes, header);
#ifdef IGA_COUPLED_CHECKPOINT_TESTING
	if (before_shard_payload) before_shard_payload(epoch, spec, file.Get());
#endif
	std::forward<Writer>(writer)(output);
	Require(output.Remaining() == 0, "writer omitted declared payload bytes"); file.Sync(); file.Close();
	PublishFile(directory, temporary, spec.id+".shard"); return {spec, bytes, output.Sha256Hex()};
}

inline CoupledCheckpointManifest PublishCoupledCheckpoint(const std::filesystem::path& root,
	const CoupledCheckpointEpoch& epoch, const std::vector<CoupledCheckpointShardSpec>& complete_catalog,
	const std::vector<CoupledCheckpointShard>& writer_receipts)
{
	using namespace coupled_checkpoint_detail; Validate(epoch); ValidateCatalog(complete_catalog, epoch.compatibility.ranks);
	auto parent = Root(root); auto directory = Directory(parent.Get(), epoch.id); RequireUnpublished(directory.Get());
	CoupledCheckpointManifest result; result.epoch = epoch;
	Require(writer_receipts.size() == complete_catalog.size(), "missing writer receipts");
	for (std::size_t i = 0; i < complete_catalog.size(); ++i) {
		const auto& expected = writer_receipts[i];
		Require(Same(expected.spec, complete_catalog[i]) && Digest(expected.sha256), "writer receipt catalog differs");
		const auto actual = Inspect(directory, epoch, complete_catalog[i], true);
		Require(actual.payload_bytes == expected.payload_bytes && actual.sha256 == expected.sha256,
			"shard differs from writer receipt");
		result.shards.push_back(actual);
	}
	const auto text = SerializeCoupledCheckpointManifest(result);
	auto file = Output(directory.Get(), "manifest.tmp"); Write(file.Get(), text.data(), text.size()); file.Sync(); file.Close();
#ifdef IGA_COUPLED_CHECKPOINT_TESTING
	// Test binaries can terminate the whole writer after every file is durable
	// but before the completion marker exists. Never enabled in production.
	if (before_manifest_publication) before_manifest_publication(epoch);
	if (std::getenv("IGA_CHECKPOINT_EXIT_BEFORE_MANIFEST")) ::_exit(86);
#endif
	PublishFile(directory, "manifest.tmp", "manifest"); return result;
}

inline CoupledCheckpointManifest LoadCoupledCheckpoint(const std::filesystem::path& root, const std::string& epoch_id,
	const CoupledCheckpointCompatibility& expected, const std::vector<CoupledCheckpointShardSpec>& complete_catalog)
{
	using namespace coupled_checkpoint_detail; Validate(expected); ValidateCatalog(complete_catalog, expected.ranks);
	auto parent = Root(root); auto directory = Directory(parent.Get(), epoch_id);
	auto result = ParseCoupledCheckpointManifest(ReadManifest(directory));
	Require(result.epoch.id == epoch_id, "manifest belongs to a different epoch directory");
	RequireCompatible(result.epoch.compatibility, expected);
	Require(result.shards.size() == complete_catalog.size(), "manifest shard catalog size differs");
	for (std::size_t i = 0; i < complete_catalog.size(); ++i) {
		const auto& stored = result.shards[i]; Require(Same(stored.spec, complete_catalog[i]), "manifest shard catalog differs");
		const auto actual = Inspect(directory, result.epoch, stored.spec, false);
		Require(actual.payload_bytes == stored.payload_bytes && actual.sha256 == stored.sha256, "shard checksum or size mismatch");
	}
	return result;
}

// Consumer must populate an unpublished restore candidate: checksum validation
// completes after streaming, so it must not mutate accepted live state.
inline void ReadCoupledCheckpointShard(const std::filesystem::path& root, const CoupledCheckpointManifest& manifest,
	std::size_t shard_index, const coupled_checkpoint_detail::Consumer& consumer)
{
	using namespace coupled_checkpoint_detail; Validate(manifest.epoch);
	Require(shard_index < manifest.shards.size() && static_cast<bool>(consumer), "invalid shard consumer");
	const auto& expected = manifest.shards[shard_index];
	auto parent = Root(root); auto directory = Directory(parent.Get(), manifest.epoch.id);
	const auto actual = Inspect(directory, manifest.epoch, expected.spec, false, consumer);
	Require(actual.payload_bytes == expected.payload_bytes && actual.sha256 == expected.sha256, "streamed shard checksum or size mismatch");
}

struct CoupledCheckpointDiscovery {
	std::optional<CoupledCheckpointManifest> latest;
	// Caller must report rejected generations when recovering an older one.
	std::vector<std::pair<std::string, std::string>> rejected;
};

inline CoupledCheckpointDiscovery FindLatestCoupledCheckpoint(const std::filesystem::path& root,
	const CoupledCheckpointCompatibility& expected, const std::vector<CoupledCheckpointShardSpec>& complete_catalog)
{
	using namespace coupled_checkpoint_detail; Validate(expected); ValidateCatalog(complete_catalog, expected.ranks);
	auto parent = Root(root); CoupledCheckpointDiscovery result;
	std::vector<CoupledCheckpointEpoch> candidates;
	std::size_t count = 0;
	for (const auto& entry : std::filesystem::directory_iterator(root)) {
		Require(++count <= 65536, "generation scan exceeds limit; select an explicit epoch");
		const auto id = entry.path().filename().string(); if (!Digest(id)) continue;
		try {
			auto directory = Directory(parent.Get(), id);
			auto manifest = ParseCoupledCheckpointManifest(ReadManifest(directory));
			Require(manifest.epoch.id == id, "manifest epoch directory differs");
			RequireCompatible(manifest.epoch.compatibility, expected); candidates.push_back(std::move(manifest.epoch));
		} catch (const std::exception& error) { result.rejected.emplace_back(id, error.what()); }
	}
	std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
		return a.accepted_steps != b.accepted_steps ? a.accepted_steps > b.accepted_steps : a.id < b.id;
	});
	for (const auto& epoch : candidates) {
		if (result.latest && epoch.accepted_steps < result.latest->epoch.accepted_steps) break;
		std::optional<CoupledCheckpointManifest> loaded;
		try { loaded = LoadCoupledCheckpoint(root, epoch.id, expected, complete_catalog); }
		catch (const std::exception& error) { result.rejected.emplace_back(epoch.id, error.what()); }
		if (loaded) {
			Require(!result.latest, "ambiguous completed epochs at the same step; select an explicit epoch");
			result.latest = std::move(loaded);
		}
	}
	return result;
}

} // namespace iga
#endif

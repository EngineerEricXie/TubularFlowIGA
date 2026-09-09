#ifndef IGA_COLLECTIVE_ASSET_INPUT_HPP
#define IGA_COLLECTIVE_ASSET_INPUT_HPP

#include "CollectiveFailure.hpp"
#include "Sha256.hpp"
#include <array>
#include <cerrno>
#include <filesystem>
#include <map>
#include <string>
#include <system_error>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef IGA_ASSET_INPUT_TESTING
#include <functional>
#endif

namespace iga {

using AssetFileCatalog = std::map<std::string, std::filesystem::path>;

#ifdef IGA_ASSET_INPUT_TESTING
inline std::function<void()>& AssetReadProbeForTesting()
{
	static thread_local std::function<void()> probe;
	return probe;
}
#endif

// Local operation. The caller coordinates exceptions before entering MPI.
// Inputs must remain immutable throughout the run. The descriptor checks
// reject changes observed during this read; they are not a file-system lock.
inline std::string ReadAssetFingerprint(const std::filesystem::path& path)
{
	struct Descriptor {
		int value = -1;
		~Descriptor() { if (value >= 0) ::close(value); }
	} descriptor;
	descriptor.value = ::open(path.c_str(), O_RDONLY | O_NONBLOCK | O_CLOEXEC);
	if (descriptor.value < 0)
		throw std::system_error(errno, std::generic_category(), "cannot open asset " + path.string());
	struct stat before{}, after{}, visible{};
	if (::fstat(descriptor.value, &before) != 0)
		throw std::system_error(errno, std::generic_category(), "cannot inspect asset " + path.string());
	if (!S_ISREG(before.st_mode) || before.st_size < 0)
		throw std::runtime_error("asset is not a regular file: " + path.string());
	Sha256 hash;
	std::array<char, 65536> buffer{};
	std::uint64_t bytes = 0;
	const auto expected = static_cast<std::uint64_t>(before.st_size);
	for (;;) {
		const auto count = ::read(descriptor.value, buffer.data(), buffer.size());
		if (count < 0 && errno == EINTR) continue;
		if (count < 0)
			throw std::system_error(errno, std::generic_category(), "cannot read asset " + path.string());
		if (!count) break;
		if (static_cast<std::uint64_t>(count) > expected - bytes)
			throw std::runtime_error("asset grew while being read: " + path.string());
		hash.Append(buffer.data(), static_cast<std::size_t>(count));
		bytes += static_cast<std::uint64_t>(count);
#ifdef IGA_ASSET_INPUT_TESTING
		if (AssetReadProbeForTesting()) AssetReadProbeForTesting()();
#endif
	}
	if (::fstat(descriptor.value, &after) != 0 || ::stat(path.c_str(), &visible) != 0)
		throw std::system_error(errno, std::generic_category(), "cannot recheck asset " + path.string());
	if (bytes != expected || before.st_size != after.st_size
		|| before.st_mtim.tv_sec != after.st_mtim.tv_sec || before.st_mtim.tv_nsec != after.st_mtim.tv_nsec
		|| before.st_ctim.tv_sec != after.st_ctim.tv_sec || before.st_ctim.tv_nsec != after.st_ctim.tv_nsec
		|| before.st_dev != visible.st_dev || before.st_ino != visible.st_ino)
		throw std::runtime_error("asset changed while being read: " + path.string());
	return std::to_string(bytes) + ":" + hash.Hex();
}

// Keys describe logical inputs, not their rank-local paths. Catalog agreement
// precedes the per-file collectives, so a missing key cannot change their order.
inline void RequireCollectiveAssetFiles(MPI_Comm communicator, const AssetFileCatalog& files)
{
	std::string catalog;
	CollectiveLocalStage(communicator, "asset catalog preparation", [&] {
		for (const auto& file : files)
			catalog += std::to_string(file.first.size()) + ":" + file.first;
	});
	RequireCollectiveSameText(communicator, "asset catalog agreement", catalog);
	std::map<std::string, std::string> fingerprints;
	CollectiveLocalStage(communicator, "asset content read", [&] {
		for (const auto& file : files)
			fingerprints.emplace(file.first, ReadAssetFingerprint(file.second));
	});
	for (const auto& fingerprint : fingerprints)
		RequireCollectiveSameText(communicator, fingerprint.first.c_str(), fingerprint.second);
}

} // namespace iga

#endif

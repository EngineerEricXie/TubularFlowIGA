#include "CollectiveAssetInput.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

namespace fs = std::filesystem;

namespace {
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void Write(const fs::path& path, const std::string& bytes)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
	output.close();
	Require(static_cast<bool>(output), "asset fixture write failed");
}

std::size_t DescriptorCount()
{
	return static_cast<std::size_t>(std::distance(fs::directory_iterator("/proc/self/fd"), fs::directory_iterator{}));
}

void RunGroup(MPI_Comm comm, const fs::path& root)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	std::string data;
	const auto directory = root / std::to_string(rank);
	const auto file = directory / "data";
	iga::CollectiveLocalStage(comm, "asset fixture", [&] {
		Require(fs::create_directories(directory), "asset fixture already exists");
		data.resize(196625);
		for (std::size_t i = 0; i < data.size(); ++i) data[i] = static_cast<char>(i % 251);
		Write(file, data); Write(directory / "empty", ""); Write(directory / "abc", "abc");
		Require(iga::ReadAssetFingerprint(directory / "empty") ==
			"0:e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855", "empty file SHA differs");
		Require(iga::ReadAssetFingerprint(directory / "abc") ==
			"3:ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", "abc file SHA differs");
		Require(iga::ReadAssetFingerprint(file) ==
			"196625:4bb7439bc39bc2d0e3d6a915d7c81e38250a9c5bb320a94a19a95bba0d5fe40a", "multi-buffer SHA differs");
		int metadata_changes = 0;
		const auto original_time = fs::last_write_time(file);
		iga::AssetReadProbeForTesting() = [&] {
			++metadata_changes;
			fs::last_write_time(file, original_time + std::chrono::seconds(metadata_changes));
		};
		const auto retried = iga::ReadAssetFingerprint(file);
		iga::AssetReadProbeForTesting() = {};
		Require(metadata_changes >= 2, "metadata retry fixture did not continuously change the timestamp");
		Require(retried == "196625:4bb7439bc39bc2d0e3d6a915d7c81e38250a9c5bb320a94a19a95bba0d5fe40a",
			"timestamp-only metadata retry changed the fingerprint");
	});
	iga::AssetFileCatalog healthy;
	iga::CollectiveLocalStage(comm, "asset healthy catalog", [&] {
		healthy = {{"asset unit data", file}, {"asset unit empty", directory / "empty"}};
	});
	iga::RequireCollectiveAssetFiles(comm, healthy);
	// Empty catalogs are valid and use the same collective protocol.
	iga::RequireCollectiveAssetFiles(comm, {});
	int cases = 0;
	const std::vector<std::string> modes{"different", "shorter", "missing-key", "extra-key", "renamed-key",
		"missing", "directory", "fifo", "device", "grow", "truncate", "overwrite", "replace", "unlink", "allocation"};
	for (const auto& mode : modes) {
		const bool catalog_difference = mode == "missing-key" || mode == "extra-key" || mode == "renamed-key";
		const bool contents_difference = mode == "different" || mode == "shorter";
		if (ranks == 1 && (catalog_difference || contents_difference)) continue;
		iga::AssetFileCatalog files;
		bool changed = false;
		iga::CollectiveLocalStage(comm, "asset failure fixture", [&] {
			files = healthy;
			fs::remove(file); Write(file, data);
			if (rank != ranks - 1) return;
			if (mode == "different") { auto other = data; other.back() ^= 1; Write(file, other); }
			if (mode == "shorter") fs::resize_file(file, data.size() - 1);
			if (mode == "missing-key") files.erase("asset unit data");
			if (mode == "extra-key") files.emplace("another asset", file);
			if (mode == "renamed-key") { files.erase("asset unit data"); files.emplace("renamed asset", file); }
			if (mode == "missing") fs::remove(file);
			if (mode == "directory") { fs::remove(file); fs::create_directory(file); }
			if (mode == "fifo") { fs::remove(file); Require(::mkfifo(file.c_str(), 0600) == 0, "cannot create FIFO"); }
			if (mode == "device") files.at("asset unit data") = "/dev/zero";
			iga::AssetReadProbeForTesting() = [&] {
				if (changed) return;
				changed = true;
				if (mode == "grow") { std::ofstream output(file, std::ios::app); output << 'x'; }
				if (mode == "truncate") fs::resize_file(file, 8);
				if (mode == "overwrite") { auto other = data; other.front() ^= 1; Write(file, other); }
				if (mode == "replace") { Write(directory / "replacement", data); fs::rename(directory / "replacement", file); }
				if (mode == "unlink") fs::remove(file);
				if (mode == "allocation") throw std::bad_alloc();
			};
		});
		const auto descriptors = DescriptorCount();
		std::string diagnostic;
		try { iga::RequireCollectiveAssetFiles(comm, files); }
		catch (const std::exception& error) { diagnostic = error.what(); }
		iga::AssetReadProbeForTesting() = {};
		iga::CollectiveLocalStage(comm, "asset rejection verification", [&] {
			const char* expected = catalog_difference ? "asset catalog agreement"
				: contents_difference ? "asset unit data" : "asset content read";
			if (diagnostic.find(expected) == std::string::npos)
				throw std::runtime_error("asset failure mode " + mode
					+ " was not rejected at expected stage; diagnostic: " + diagnostic);
			Require(DescriptorCount() == descriptors, "asset failure leaked a file descriptor");
			fs::remove(file); Write(file, data);
		});
		iga::RequireCollectiveSameText(comm, "asset common diagnostic", diagnostic);
		iga::RequireCollectiveAssetFiles(comm, healthy);
		++cases;
	}
	if (rank == 0) std::cout << "asset_input ranks=" << ranks << " failures=" << cases
		<< " known_digests=3 replica_paths=passed descriptor_cleanup=passed retry=passed\n";
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1, status = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(argc == 2 && ranks == 3, "usage: mpiexec -np 3 collective_asset_input_test FRESH_OUTPUT");
		RunGroup(PETSC_COMM_WORLD, fs::path(argv[1]) / "world");
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		RunGroup(group, fs::path(argv[1]) / (rank == 0 ? "single" : "pair"));
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "asset test rank " << rank << ": " << error.what() << '\n'; status = 1;
	}
	PetscFinalize();
	return status;
}

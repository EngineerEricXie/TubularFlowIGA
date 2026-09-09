#include "CoupledCheckpointBundle.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <sys/resource.h>
#include <sys/wait.h>

namespace fs = std::filesystem;
using namespace iga;
using coupled_checkpoint_detail::Require;

namespace {
int rejected = 0;
template<class Work> void Reject(Work&& work)
{
	try { work(); } catch (const std::exception&) { ++rejected; return; }
	throw std::runtime_error("expected checkpoint rejection");
}

std::string Digest(const std::string& name) { return coupled_checkpoint_detail::Hash(name); }

CoupledCheckpointEpoch Epoch(const std::string& name, std::uint64_t step)
{
	return {Digest(name), {}, {Digest("case"), Digest("configuration"), Digest("execution"), 2}, step, step*0.1, 0.1};
}

std::vector<CoupledCheckpointShardSpec> Catalog()
{
	return {{"domain-a-rank-0", "owned-f64-v1", 0}, {"graph-state", "graph-v1", 1}};
}

std::string Read(const fs::path& path)
{
	std::ifstream input(path, std::ios::binary); input.exceptions(std::ios::badbit);
	Require(input.is_open(), "test input unavailable");
	return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void Replace(const fs::path& path, const std::string& data)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	output.exceptions(std::ios::failbit | std::ios::badbit);
	output.write(data.data(), static_cast<std::streamsize>(data.size())); output.close();
}

void Child(const std::function<void()>& work, int expected)
{
	const auto pid = ::fork(); Require(pid >= 0, "fork failed");
	if (pid == 0) {
		try { work(); ::_exit(0); }
		catch (const std::exception& error) { std::cerr << error.what() << '\n'; ::_exit(99); }
	}
	int status = 0; Require(::waitpid(pid, &status, 0) == pid && WIFEXITED(status)
		&& WEXITSTATUS(status) == expected, "child exit differed");
}

std::vector<CoupledCheckpointShard> WriteEpoch(const fs::path& root, const CoupledCheckpointEpoch& epoch)
{
	CreateCoupledCheckpointEpoch(root, epoch); std::vector<CoupledCheckpointShard> receipts;
	struct Producer { pid_t pid; int receipt_fd; CoupledCheckpointShardSpec spec; };
	std::vector<Producer> producers;
	for (const auto& spec : Catalog()) {
		// Independent producer process. Receipt travels through a pipe, not a
		// later rehash which could inadvertently bless damaged storage.
		int pipefd[2]; Require(::pipe(pipefd) == 0, "pipe failed");
		const auto pid = ::fork(); Require(pid >= 0, "producer fork failed");
		if (pid == 0) {
		try {
			::close(pipefd[0]); const std::string data = spec.id+std::string("\0payload\n", 9);
			const auto receipt = WriteCoupledCheckpointShard(root, epoch, spec, data.size(), [&](auto& out) {
				out.Write(data.data(), 3); out.Write(data.data()+3, data.size()-3);
			});
			coupled_checkpoint_detail::Write(pipefd[1], receipt.sha256.data(), 64); ::close(pipefd[1]); ::_exit(0);
		} catch (...) { ::_exit(99); }
		}
		::close(pipefd[1]); producers.push_back({pid, pipefd[0], spec});
	}
	for (const auto& producer : producers) {
		int status = 0; Require(::waitpid(producer.pid, &status, 0) == producer.pid && WIFEXITED(status)
			&& WEXITSTATUS(status) == 0, "producer failed");
		std::string checksum(64, '\0');
		Require(coupled_checkpoint_detail::Read(producer.receipt_fd, checksum.data(), 64) == 64, "missing writer receipt");
		::close(producer.receipt_fd); receipts.push_back({producer.spec, producer.spec.id.size()+9, checksum});
	}
	return receipts;
}

void ManifestTests(const CoupledCheckpointManifest& baseline)
{
	const auto text = SerializeCoupledCheckpointManifest(baseline);
	Require(SerializeCoupledCheckpointManifest(ParseCoupledCheckpointManifest(text)) == text, "manifest roundtrip failed");
	auto exact = baseline; exact.epoch.accepted_steps = UINT64_MAX;
	exact.epoch.time_s = 0; for (int i = 0; i < 10; ++i) exact.epoch.time_s += 0.1;
	const auto recovered = ParseCoupledCheckpointManifest(SerializeCoupledCheckpointManifest(exact));
	Require(recovered.epoch.accepted_steps == UINT64_MAX
		&& coupled_checkpoint_detail::DoubleBits(recovered.epoch.time_s) == coupled_checkpoint_detail::DoubleBits(exact.epoch.time_s), "clock precision lost");
	for (std::size_t i = 0; i < text.size(); ++i) Reject([&] { ParseCoupledCheckpointManifest(text.substr(0, i)); });
	for (const auto suffix : {"x", "\n", "shards 1\n"}) Reject([&] { ParseCoupledCheckpointManifest(text+suffix); });
	for (const auto& change : std::vector<std::pair<std::string, std::string>>{
		{"IGA_COUPLED_CHECKPOINT 1", "IGA_COUPLED_CHECKPOINT 2"}, {"accepted_steps 3", "accepted_steps 18446744073709551616"},
		{"ranks 2", "ranks -1"}, {"shards 2", "shards 18446744073709551615"}, {"ranks 2", "ranks 02"}}) {
		auto modified = text; const auto at = modified.find(change.first); Require(at != std::string::npos, "test mutation target absent");
		modified.replace(at, change.first.size(), change.second); Reject([&] { ParseCoupledCheckpointManifest(modified); });
	}
	Reject([&] { ParseCoupledCheckpointManifest(std::string(coupled_checkpoint_detail::maximum_manifest_bytes+1, ' ')); });
	for (const double invalid : {0., -1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
		auto candidate = baseline; candidate.epoch.time_s = invalid;
		Reject([&] { SerializeCoupledCheckpointManifest(candidate); });
	}
	auto candidate = baseline; candidate.shards[1].spec.id = candidate.shards[0].spec.id;
	Reject([&] { SerializeCoupledCheckpointManifest(candidate); });
	candidate = baseline; candidate.shards[0].spec.id = "../escape";
	Reject([&] { SerializeCoupledCheckpointManifest(candidate); });
	candidate = baseline; candidate.shards[0].spec.owner_rank = 2;
	Reject([&] { SerializeCoupledCheckpointManifest(candidate); });
}
} // namespace

int main(int argc, char** argv)
{
	try {
		Require(argc == 2, "test requires an unused output directory"); const fs::path root(argv[1]);
		Require(fs::create_directory(root), "test directory must be new");
		const auto begin = std::chrono::steady_clock::now();
		auto first = Epoch("first", 3); const auto specs = Catalog();
		const auto receipts = WriteEpoch(root, first);
		Reject([&] { LoadCoupledCheckpoint(root, first.id, first.compatibility, specs); });
		const auto baseline = PublishCoupledCheckpoint(root, first, specs, receipts);
		const auto baseline_bytes = Read(root/first.id/"manifest"); ManifestTests(baseline);
		Child([&] {
			const auto loaded = LoadCoupledCheckpoint(root, first.id, first.compatibility, specs);
			for (std::size_t i = 0; i < specs.size(); ++i) {
				std::string data; ReadCoupledCheckpointShard(root, loaded, i, [&](const void* p, std::size_t n) { data.append(static_cast<const char*>(p), n); });
				Require(data == specs[i].id+std::string("\0payload\n", 9), "payload roundtrip failed");
			}
		}, 0);
		Reject([&] { CreateCoupledCheckpointEpoch(root, first); });
		Reject([&] { PublishCoupledCheckpoint(root, first, specs, receipts); });
		Reject([&] { WriteCoupledCheckpointShard(root, first, specs[0], 0, [](auto&) {}); });
		for (int item = 0; item < 4; ++item) {
			auto incompatible = first.compatibility;
			if (item == 0) incompatible.case_sha256 = Digest("wrong");
			if (item == 1) incompatible.configuration_sha256 = Digest("wrong");
			if (item == 2) incompatible.execution_sha256 = Digest("wrong");
			if (item == 3) incompatible.ranks = 3;
			Reject([&] { LoadCoupledCheckpoint(root, first.id, incompatible, specs); });
		}
		auto second = Epoch("second", 4); second.previous_id = first.id;
		const auto second_receipts = WriteEpoch(root, second);
		Child([&] { ::setenv("IGA_CHECKPOINT_EXIT_BEFORE_MANIFEST", "1", 1);
			PublishCoupledCheckpoint(root, second, specs, second_receipts); }, 86);
		Require(fs::exists(root/second.id/"manifest.tmp") && !fs::exists(root/second.id/"manifest"), "early completion marker");
		auto broken = Epoch("partial", 5); CreateCoupledCheckpointEpoch(root, broken);
		Child([&] { WriteCoupledCheckpointShard(root, broken, specs[0], 100, [&](auto& out) {
			out.Write("partial", 7); ::_exit(87); }); }, 87);
		Reject([&] { PublishCoupledCheckpoint(root, broken, specs, receipts); });
		Child([&] { const auto found = FindLatestCoupledCheckpoint(root, first.compatibility, specs);
			Require(found.latest && found.latest->epoch.id == first.id && found.rejected.size() == 2, "restart failed to recover prior complete epoch"); }, 0);
		Require(Read(root/first.id/"manifest") == baseline_bytes, "previous manifest changed");

		auto third = Epoch("third", 6); third.previous_id = first.id;
		const auto third_receipts = WriteEpoch(root, third);
		const auto path = root/third.id/(specs[0].id+".shard"); const auto original = Read(path);
		auto corrupt = original; corrupt.back() ^= 1; Replace(path, corrupt);
		Reject([&] { PublishCoupledCheckpoint(root, third, specs, third_receipts); }); // receipt catches prepublication damage
		Replace(path, original); PublishCoupledCheckpoint(root, third, specs, third_receipts);
		for (const auto& data : {corrupt, original.substr(0, original.size()-1), original+"x", Read(root/first.id/(specs[0].id+".shard"))}) {
			Replace(path, data); Reject([&] { LoadCoupledCheckpoint(root, third.id, third.compatibility, specs); });
			Require(FindLatestCoupledCheckpoint(root, first.compatibility, specs).latest->epoch.id == first.id, "corruption fallback failed");
		}
		fs::remove(path); Reject([&] { LoadCoupledCheckpoint(root, third.id, third.compatibility, specs); });
		fs::create_symlink(root/first.id/(specs[0].id+".shard"), path);
		Reject([&] { LoadCoupledCheckpoint(root, third.id, third.compatibility, specs); });
		fs::remove(path); Replace(path, original);
		Require(FindLatestCoupledCheckpoint(root, first.compatibility, specs).latest->epoch.id == third.id, "latest complete selection failed");
		auto wrong_catalog = specs; wrong_catalog[0].format = "other-v1";
		Reject([&] { LoadCoupledCheckpoint(root, third.id, third.compatibility, wrong_catalog); });
		for (int fault = 0; fault < 5; ++fault) {
			auto failed = Epoch("io-failure-"+std::to_string(fault), 7); CreateCoupledCheckpointEpoch(root, failed);
			if (fault == 0) coupled_checkpoint_detail::fail_next_write = true;
			if (fault == 1) coupled_checkpoint_detail::fail_next_sync = true;
			Reject([&] { WriteCoupledCheckpointShard(root, failed, specs[0], 3, [&](auto& out) {
				if (fault == 2) out.Write("1234", 4);
				else if (fault == 3) out.Write("1", 1);
				else if (fault == 4) throw std::runtime_error("producer failed");
				else out.Write("123", 3);
			}); });
			Require(!fs::exists(root/failed.id/"manifest"), "failed writer published epoch");
		}
		auto missing = Epoch("missing-receipt", 7); const auto missing_receipts = WriteEpoch(root, missing);
		Reject([&] { PublishCoupledCheckpoint(root, missing, specs, {missing_receipts[0]}); });
		coupled_checkpoint_detail::fail_next_sync = true;
		Reject([&] { PublishCoupledCheckpoint(root, missing, specs, missing_receipts); });
		Require(!fs::exists(root/missing.id/"manifest"), "failed sync published epoch");
		coupled_checkpoint_detail::fail_next_write = true;
		Reject([&] { PublishCoupledCheckpoint(root, missing, specs, missing_receipts); });
		Require(!fs::exists(root/missing.id/"manifest"), "failed manifest write published epoch");
		Child([&] {
			const auto found = FindLatestCoupledCheckpoint(root, first.compatibility, specs);
			Require(found.latest && found.latest->epoch.id == third.id, "fresh recovery after I/O failures differs");
		}, 0);
		// Completion marker corruption is rejected, even when all shards are intact.
		const auto manifest_path = root/third.id/"manifest"; const auto good_manifest = Read(manifest_path);
		Replace(manifest_path, good_manifest+"\n");
		Reject([&] { LoadCoupledCheckpoint(root, third.id, third.compatibility, specs); });
		Replace(manifest_path, good_manifest);
		auto duplicate = Epoch("duplicate-step", 6); const auto duplicate_receipts = WriteEpoch(root, duplicate);
		PublishCoupledCheckpoint(root, duplicate, specs, duplicate_receipts);
		Reject([&] { FindLatestCoupledCheckpoint(root, first.compatibility, specs); });
		// Explicit selection remains valid when discovery reports an ambiguous fork.
		LoadCoupledCheckpoint(root, third.id, third.compatibility, specs);

		const auto large_root = root/"large"; fs::create_directory(large_root);
		auto large = Epoch("streaming", 100);
		const std::vector<CoupledCheckpointShardSpec> large_specs{{"empty-rank", "bytes-v1", 1}, {"field", "bytes-v1", 0}};
		CreateCoupledCheckpointEpoch(large_root, large); const auto io_begin = std::chrono::steady_clock::now();
		const auto empty_receipt = WriteCoupledCheckpointShard(large_root, large, large_specs[0], 0, [](auto&) {});
		constexpr std::uint64_t payload_bytes = 32*1024*1024;
		const auto large_receipt = WriteCoupledCheckpointShard(large_root, large, large_specs[1], payload_bytes, [&](auto& out) {
			std::array<char, 65536> chunk{}; for (std::size_t i = 0; i < chunk.size(); ++i) chunk[i] = static_cast<char>(i%127);
			for (std::uint64_t i = 0; i < payload_bytes/chunk.size(); ++i) out.Write(chunk.data(), chunk.size());
		});
		const auto write_end = std::chrono::steady_clock::now();
		const auto large_manifest = PublishCoupledCheckpoint(large_root, large, large_specs, {empty_receipt, large_receipt});
		const auto publish_end = std::chrono::steady_clock::now();
		LoadCoupledCheckpoint(large_root, large.id, large.compatibility, large_specs);
		std::uint64_t consumed = 0;
		ReadCoupledCheckpointShard(large_root, large_manifest, 0, [](const void*, std::size_t) { throw std::runtime_error("empty shard supplied data"); });
		ReadCoupledCheckpointShard(large_root, large_manifest, 1, [&](const void* pointer, std::size_t bytes) {
			const auto* data = static_cast<const char*>(pointer);
			for (std::size_t i = 0; i < bytes; ++i) Require(data[i] == static_cast<char>(((consumed+i)%65536)%127), "large streamed payload differs");
			consumed += bytes;
		});
		Require(consumed == payload_bytes, "streamed payload size differs");
		const auto end = std::chrono::steady_clock::now(); struct rusage usage{}, children{};
		Require(::getrusage(RUSAGE_SELF, &usage) == 0 && ::getrusage(RUSAGE_CHILDREN, &children) == 0, "RSS probe failed");
		std::cout << "checkpoint_bundle status=passed rejections=" << rejected
			<< " payload_bytes=" << payload_bytes << " write_sync_s=" << std::chrono::duration<double>(write_end-io_begin).count()
			<< " verify_publish_sync_s=" << std::chrono::duration<double>(publish_end-write_end).count()
			<< " load_stream_s=" << std::chrono::duration<double>(end-publish_end).count()
			<< " total_s=" << std::chrono::duration<double>(end-begin).count() << " peak_rss_kib=" << usage.ru_maxrss
			<< " max_child_peak_rss_kib=" << children.ru_maxrss << "\n";
		return 0;
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

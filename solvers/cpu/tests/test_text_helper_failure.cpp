#include "OneDCheckpoint.hpp"
#include "TemporalVtkHdf.hpp"
#include "VtkOutput.hpp"
#include "ExecutionResources.hpp"

#include <cstdlib>
#include <iostream>
#include <new>

namespace allocation_test {
thread_local bool active = false;
thread_local std::size_t count = 0, fail_at = 0;
thread_local bool injected = false;
}

[[gnu::noinline]] void* operator new(std::size_t size)
{
	if (allocation_test::active && ++allocation_test::count == allocation_test::fail_at) {
		allocation_test::active = false;
		allocation_test::injected = true;
		throw std::bad_alloc();
	}
	if (auto* value = std::malloc(size ? size : 1)) return value;
	throw std::bad_alloc();
}

[[gnu::noinline]] void operator delete(void* value) noexcept { std::free(value); }
[[gnu::noinline]] void operator delete(void* value, std::size_t) noexcept { std::free(value); }

namespace {
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

// Count the allocations made by the real helper, then fail each one in turn.
// Inputs exist before injection. Agreement and error diagnostics run unarmed.
template<class Work>
void Sweep(MPI_Comm comm, const char* name, Work&& work)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	std::string baseline;
	std::size_t allocations = 0;
	iga::CollectiveLocalStage(comm, "helper baseline", [&] {
		allocation_test::count = allocation_test::fail_at = 0;
		allocation_test::active = true;
		try { baseline = work(); }
		catch (...) { allocation_test::active = false; throw; }
		allocation_test::active = false;
		allocations = allocation_test::count;
		Require(allocations > 0, "fixture did not allocate");
	});
	iga::RequireCollectiveSameInt(comm, "helper allocation count", static_cast<int>(allocations));
	for (std::size_t ordinal = 1; ordinal <= allocations; ++ordinal) {
		bool failed = false;
		allocation_test::injected = false;
		try {
			iga::CollectiveLocalStage(comm, name, [&] {
				allocation_test::count = 0;
				allocation_test::fail_at = ordinal;
				allocation_test::active = rank == ranks-1;
				try { (void)work(); }
				catch (...) { allocation_test::active = false; throw; }
				allocation_test::active = false;
			});
		} catch (const std::exception&) { failed = true; }
		allocation_test::active = false;
		iga::CollectiveLocalStage(comm, "helper fault checks", [&] {
			Require(allocation_test::injected == (rank == ranks-1), "allocation injection missed target");
			Require(failed, "helper swallowed allocation failure");
			Require(work() == baseline, "healthy retry changed helper result");
		});
	}
	iga::RequireCollectiveSameText(comm, "helper baseline agreement", baseline);
	if (rank == 0) std::cout << name << " ranks=" << ranks << " faults=" << allocations
		<< " retries=" << allocations << " passed\n";
}

void Run(MPI_Comm comm)
{
	iga::OneDNetwork network;
	iga::OneDCheckpointMetadata metadata;
	std::vector<iga::VtkPointArray> arrays;
	std::filesystem::path path;
	iga::CollectiveLocalStage(comm, "helper fixtures", [&] {
		network.nodes.resize(128);
		for (std::size_t i = 0; i < network.nodes.size(); ++i) {
			auto& node = network.nodes[i];
			node.id = static_cast<int>(i)+1; node.parent_id = static_cast<int>(i);
			node.position = {{0.125*static_cast<double>(i), 0.0, -0.5}};
			node.radius = 0.0125;
		}
		metadata.state_file = std::string(2048, 'x')+"\"\\\n.state";
		metadata.species = {"oxygen", "carbon_dioxide"};
		arrays.resize(2);
		arrays[0].name = std::string(2048, 'v'); arrays[0].components = 3;
		arrays[1].name = "pressure"; arrays[1].components = 1;
		path = std::string(2048, 'p')+".txt";
	});
	Sweep(comm, "1d network fingerprint", [&] { return std::to_string(iga::OneDNetworkFingerprint(network)); });
	Sweep(comm, "1d checkpoint metadata", [&] { return iga::SerializeOneDCheckpointMetadata(metadata); });
	Sweep(comm, "VTKHDF array schema", [&] { return iga::hdf_detail::ArraySchema(arrays); });
	Sweep(comm, "VTU step filename", [&] { return iga::VtuStepPath(path, 42).string(); });
	const std::array<int, 6> minimum{{1, 0, 1, 1, -1, -1}}, maximum{{1, 0, 1, 2, -1, -1}};
	Sweep(comm, "execution resource text", [&] {
		return iga::execution_resources_detail::FormatReport(3, minimum, maximum);
	});
}
} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Require(ranks == 3, "test requires three ranks");
		Run(PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		Run(group); MPI_Comm_free(&group);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		allocation_test::active = false;
		std::cerr << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1); return 1;
	}
}

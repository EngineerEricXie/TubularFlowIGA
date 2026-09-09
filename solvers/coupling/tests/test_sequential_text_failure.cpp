#define IGA_SEQUENTIAL_NO_MAIN
#include "../src/iga_1d_3d_explicit.cpp"
#include "FlowCheckpoint.hpp"
#include "TransportCheckpoint.hpp"
#include "VcaCheckpoint.hpp"

#include <cstdlib>
#include <new>

namespace allocation_test {
thread_local bool fail_next = false;
thread_local bool injected = false;
}

// Fail the first allocation only inside the selected local operation.
// MPI/PETSc initialization, agreement and diagnostics run with injection off.
[[gnu::noinline]] void* operator new(std::size_t size)
{
	if (allocation_test::fail_next) {
		allocation_test::fail_next = false;
		allocation_test::injected = true;
		throw std::bad_alloc();
	}
	if (auto* value = std::malloc(size ? size : 1)) return value;
	throw std::bad_alloc();
}

[[gnu::noinline]] void operator delete(void* value) noexcept { std::free(value); }
[[gnu::noinline]] void operator delete(void* value, std::size_t) noexcept { std::free(value); }

namespace {

void RequireText(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

void TestText(MPI_Comm comm, const fs::path& root)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	std::string input, expected;
	const auto file = root/"input.json";
	iga::CollectiveLocalStage(comm, "text fixture", [&] {
		input = std::string(9000, 'x')+"\"\\\n\t\x01";
		expected = std::string(9000, 'x')+"\\\"\\\\\\n\\t\\u0001";
		if (rank != 0) return;
		RequireText(fs::create_directories(root), "test output must be new");
		std::ofstream output(file, std::ios::binary);
		output << input; output.close();
		RequireText(bool(output), "cannot write text fixture");
	});
	const std::array<const char*, 5> stages{{"native text read", "native JSON escape",
		"flow checkpoint text", "transport checkpoint text", "VCA checkpoint text"}};
	for (std::size_t operation = 0; operation < stages.size(); ++operation) {
		const char* stage = stages[operation];
		const auto convert = [&]() -> std::string {
			switch (operation) {
			case 0: return ReadText(file);
			case 1: return JsonEscape(input);
			case 2: return iga::SerializeFlowCheckpointMetadata(iga::FlowCheckpointMetadata{});
			case 3: return iga::SerializeTransportCheckpointMetadata(iga::TransportCheckpointMetadata{});
			default: return iga::SerializeVcaCheckpointMetadata(iga::VcaCheckpointMetadata{});
			}
		};
		std::string baseline;
		iga::CollectiveLocalStage(comm, "native text baseline", [&] {
			baseline = convert();
			if (operation < 2) RequireText(baseline == (operation == 1 ? expected : input), "incorrect baseline text");
			else RequireText(baseline.find("schema_version") != std::string::npos, "missing checkpoint text");
		});
		bool failed = false;
		allocation_test::injected = false;
		try {
			iga::CollectiveLocalStage(comm, stage, [&] {
				allocation_test::fail_next = rank == ranks-1;
				try {
					(void)convert();
				} catch (...) {
					allocation_test::fail_next = false;
					throw;
				}
				allocation_test::fail_next = false;
			});
		} catch (const std::exception&) { failed = true; }
		allocation_test::fail_next = false;
		iga::CollectiveLocalStage(comm, "native text failure checks", [&] {
			RequireText(failed, "allocation failure was swallowed by text conversion");
			RequireText(allocation_test::injected == (rank == ranks-1), "allocation injection missed its rank");
		});
		std::string recovered;
		iga::CollectiveLocalStage(comm, "native text retry", [&] {
			recovered = convert();
			RequireText(recovered == baseline, "text retry changed content");
		});
		iga::RequireCollectiveSameText(comm, "native text retry agreement", recovered);
	}
	if (rank == 0) std::cout << "native text ranks=" << ranks << " faults=5 retries=5 passed\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 1;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		RequireText(argc == 2 && ranks == 3, "require three ranks and a new output directory");
		TestText(PETSC_COMM_WORLD, fs::path(argv[1])/"world");
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		TestText(group, fs::path(argv[1])/(rank == 0 ? "single" : "pair"));
		MPI_Comm_free(&group);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		allocation_test::fail_next = false;
		std::cerr << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); return 2;
	}
}

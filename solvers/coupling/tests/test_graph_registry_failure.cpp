#include "MultidomainRunner.hpp"
#include "RuntimeConstruction.hpp"

#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

namespace {

void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

struct Retained {
	std::array<PetscObject, 128> objects{};
	std::size_t count = 0;
	void Add(PetscObject object)
	{
		Check(count < objects.size(), "too many retained runtime objects");
		Check(PetscObjectReference(object) == 0, "cannot retain runtime object");
		objects[count++] = object;
	}
	void Release(MPI_Comm comm)
	{
		iga::CollectiveLocalStage(comm, "graph retained coverage", [&] {
			Check(count > 0, "no runtime objects were observed");
		});
		iga::RequireCollectiveSameInt(comm, "graph retained object count", static_cast<int>(count));
		while (count) {
			const auto object = objects[--count];
			iga::CollectiveLocalStage(comm, "graph released ownership", [&] {
				PetscInt references = 0;
				Check(PetscObjectGetReference(object, &references) == 0 && references == 1,
					"graph runner leaked a runtime ownership reference");
			});
			iga::RequireCollectivePetscSuccess(comm, "graph retained release", PetscObjectDereference(object));
		}
	}
};

int Run(MPI_Comm comm, const fs::path& fixture, const fs::path& output)
{
	std::vector<std::string> args{"graph_registry_failure_test", "--graph-case",
		fixture.string(), "--output-dir", output.string()};
	std::vector<char*> raw;
	for (auto& arg : args) raw.push_back(arg.data());
	return iga::RunMultidomainFlow(static_cast<int>(raw.size()), raw.data(), comm);
}

void RunCases(MPI_Comm comm, const fs::path& parent, bool species, bool zero_d = false)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const auto fixture = parent/(zero_d ? "zero-d" : (species ? "species" : "flow"));
	iga::CollectiveLocalStage(comm, "graph registry fixture", [&] {
		if (rank != 0) return;
		Check(fs::create_directories(fixture), "fixture must be new");
		if (zero_d) {
			for (const auto& id : {"zero_source", "zero_terminal", "zero_junction"}) fs::create_directory(fixture/id);
			WriteZeroDModel(fixture/"zero_source", true);
			WriteZeroDModel(fixture/"zero_terminal", false);
			WriteZeroDThreeDCase(fixture/"zero_junction");
			WriteDatabase(fixture/"zero_junction.ntiga", ranks);
			WriteZeroDClockGraph(fixture);
		} else if (species) {
			for (const auto& id : {"species_source", "species_leaf_a", "species_leaf_b", "junction"})
				fs::create_directory(fixture/id);
			WriteOneDTransportCase(fixture/"species_source", 1e-3, "source_red", "source_blue");
			WriteOneDTransportCase(fixture/"species_leaf_a", 6e-4, "leaf_a_red", "leaf_a_blue");
			WriteOneDTransportCase(fixture/"species_leaf_b", 4e-4, "leaf_b_red", "leaf_b_blue");
			WriteThreeDTransportCase(fixture/"junction");
			WriteDatabase(fixture/"graph.ntiga", ranks);
			WriteSpeciesGraph(fixture, "graph.ntiga");
		} else {
			for (const auto& id : {"three_d", "source", "branch_a", "branch_b"}) fs::create_directory(fixture/id);
			WriteThreeDCase(fixture/"three_d");
			WriteOneDCase(fixture/"source", 1e-3);
			WriteOneDCase(fixture/"branch_a", 6e-4);
			WriteOneDCase(fixture/"branch_b", 4e-4);
			WriteDatabase(fixture/"graph.ntiga", ranks);
			WriteGraph(fixture, "graph.ntiga", "explicit");
		}
	});
	std::vector<std::string> stages{"graph boundary label", "graph boundary coverage",
		"graph initial port catalog", "graph adapter catalog", "registry synchronization ready",
		"registry index ready", "registry allocation"};
	if (species) {
		stages.push_back("graph transport preflight");
		stages.push_back("graph transport node agreement");
	}
	std::size_t observed_objects = 0;
	for (std::size_t i = 0; i < stages.size(); ++i) {
		const auto failed_output = fixture/("failed-"+std::to_string(i));
		const auto retry_output = fixture/("retry-"+std::to_string(i));
		bool injected = false, registry_ready = false;
		Retained retained;
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		iga::RuntimeConstructionHooksForTesting().after_local_stage = [&](const char* stage) {
			if (std::string_view(stage) == "registry synchronization ready") registry_ready = true;
			if (rank != ranks-1 || injected) return;
			if (stages[i] == stage || (stages[i] == "registry allocation" && registry_ready
				&& std::string_view(stage) == "runtime object allocation")) {
				injected = true; throw std::bad_alloc();
			}
		};
		const int failed_status = Run(comm, fixture, failed_output);
		iga::RuntimeConstructionHooksForTesting() = {};
		observed_objects += retained.count;
		retained.Release(comm);
		iga::CollectiveLocalStage(comm, "graph registry expected failure", [&] {
			Check(failed_status == 1, "runner did not reject registry/construction fault");
			if (rank == ranks-1) Check(injected, "target stage was not reached");
			if (rank == 0) Check(!fs::exists(failed_output), "failed construction published output");
		});
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		const int retry_status = Run(comm, fixture, retry_output);
		iga::RuntimeConstructionHooksForTesting() = {};
		observed_objects += retained.count;
		retained.Release(comm);
		iga::CollectiveLocalStage(comm, "graph registry retry validation", [&] {
			Check(retry_status == 0, "healthy graph retry failed");
			if (rank != 0) return;
			if (zero_d) ValidateZeroDClockRun(retry_output);
			else if (species) ValidateSpeciesMultidomain(retry_output);
			else Validate(retry_output, true);
		});
	}
	if (rank == 0) std::cout << "graph registry ranks=" << ranks << " species=" << species << " zero_d=" << zero_d
		<< " faults=" << stages.size() << " retries=" << stages.size()
		<< " observed_objects=" << observed_objects << " retained_objects_released\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Check(argc == 2 && ranks == 3, "usage: mpiexec -np 3 graph_registry_failure_test NEW_OUTPUT_PARENT");
		const fs::path root(argv[1]);
		RunCases(PETSC_COMM_WORLD, root/"world", false);
		RunCases(PETSC_COMM_WORLD, root/"world", true);
		RunCases(PETSC_COMM_WORLD, root/"world", false, true);
		MPI_Comm subgroup = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &subgroup);
		const auto group_root = root/(rank == 0 ? "single" : "pair");
		RunCases(subgroup, group_root, false);
		RunCases(subgroup, group_root, true);
		RunCases(subgroup, group_root, false, true);
		MPI_Comm_free(&subgroup);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); return 2;
	}
}

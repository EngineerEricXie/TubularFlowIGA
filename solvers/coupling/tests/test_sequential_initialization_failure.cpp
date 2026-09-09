#include "MultidomainRunner.hpp"
#include "RuntimeConstruction.hpp"

#define main SequentialSmokeFixtureMain
#include "test_explicit_coupling_smoke.cpp"
#undef main

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
		Check(count < objects.size(), "too many runtime objects");
		Check(PetscObjectReference(object) == 0, "cannot retain runtime object");
		objects[count++] = object;
	}
	void Release(MPI_Comm comm, bool require_objects)
	{
		iga::CollectiveLocalStage(comm, "sequential retained coverage", [&] {
			if (require_objects) Check(count > 0, "no runtime objects observed");
		});
		iga::RequireCollectiveSameInt(comm, "sequential retained count", static_cast<int>(count));
		while (count) {
			const auto object = objects[--count];
			iga::CollectiveLocalStage(comm, "sequential released ownership", [&] {
				PetscInt references = 0;
				Check(PetscObjectGetReference(object, &references) == 0 && references == 1,
					"runner retained a runtime ownership reference");
			});
			iga::RequireCollectivePetscSuccess(comm, "sequential retained release", PetscObjectDereference(object));
		}
	}
};

std::vector<std::string> Arguments(const fs::path& fixture, const fs::path& output, bool graph)
{
	std::vector<std::string> args{"sequential_initialization_failure_test"};
	if (graph) args.insert(args.end(), {"--graph-case", fixture.string()});
	else args.insert(args.end(), {(fixture/"one.ntiga").string(), (fixture/"three_d").string(),
		(fixture/"upstream").string(), (fixture/"downstream").string(), "--upstream-terminal-node", "2"});
	args.insert(args.end(), {"--output-dir", output.string()});
	return args;
}

int Invoke(MPI_Comm comm, std::vector<std::string> args)
{
	std::vector<char*> raw;
	for (auto& arg : args) raw.push_back(arg.data());
	return iga::RunSequentialFlow(static_cast<int>(raw.size()), raw.data(), comm);
}

void ValidateOutput(const fs::path& output, bool graph, const std::string& scheme = "explicit")
{
	if (scheme == "explicit") {
		ValidateRows(ReadHistory(output/"explicit_coupling_history.csv"));
		RequireManifest(output/"explicit_coupling_manifest.json");
	} else {
		ValidateStrongRun(ReadHistory(output/"strong_coupling_history.csv"),
			ReadHistory(output/"strong_coupling_iterations.csv"), output/"strong_coupling_manifest.json", false);
	}
	if (graph) RequireGraphBindingManifest(output/"graph_binding_manifest.json", scheme);
}

void RunCases(MPI_Comm comm, const fs::path& parent, bool graph)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const auto fixture = parent/(graph ? "graph" : "positional");
	iga::CollectiveLocalStage(comm, "sequential test fixture", [&] {
		if (rank != 0) return;
		Check(fs::create_directories(fixture), "fixture must be new");
		for (const auto& name : {"three_d", "upstream", "downstream"}) fs::create_directory(fixture/name);
		WriteThreeDCase(fixture/"three_d");
		WriteOneDCase(fixture/"upstream"); WriteOneDCase(fixture/"downstream");
		WriteUnitDatabase(fixture/"one.ntiga", ranks);
		WriteGraphCase(fixture);
	});
	const std::vector<std::pair<std::string, bool>> targets{
		{"sequential arguments", false}, {"sequential failure control", false},
		{"sequential controls", false}, {"sequential option names", false},
		{"sequential case catalog", false}, {"sequential configurations", false},
		{"sequential selected assets", false}, {"sequential native inputs", false},
		{"runtime object allocation", false}, {"sequential scalar preflight", true},
		{"sequential port identifiers", true}, {"sequential graph binding", true},
		{"sequential port validation", true}, {"sequential initial state input", true},
		{"sequential initial measurements", true}, {"sequential initial lag", true},
		{"sequential history preparation", true}, {"1d trial preparation", true},
		{"1d substep preparation", true}};
	std::vector<std::pair<std::string, bool>> stages = targets;
	if (graph) {
		stages.emplace_back("sequential graph paths", false);
		stages.emplace_back("sequential graph definition", false);
	}
	std::size_t observed = 0;
	for (std::size_t i = 0; i < stages.size(); ++i) {
		bool injected = false;
		Retained retained;
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		iga::RuntimeConstructionHooksForTesting().after_local_stage = [&](const char* stage) {
			if (rank == ranks-1 && !injected && stages[i].first == stage) {
				injected = true;
				if (i == 3) throw 7; // The common stage also handles non-standard exceptions.
				throw std::bad_alloc();
			}
		};
		const auto failed_output = fixture/("failed-"+std::to_string(i));
		const auto retry_output = fixture/("retry-"+std::to_string(i));
		const int failed = Invoke(comm, Arguments(fixture, failed_output, graph));
		iga::RuntimeConstructionHooksForTesting() = {};
		observed += retained.count;
		retained.Release(comm, stages[i].second);
		iga::CollectiveLocalStage(comm, "sequential expected failure", [&] {
			Check(failed == 1, "runner accepted local initialization fault");
			if (rank == ranks-1) Check(injected, "target stage not reached");
			if (rank == 0) Check(!fs::exists(failed_output), "failure published output");
		});
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		const int retry = Invoke(comm, Arguments(fixture, retry_output, graph));
		iga::RuntimeConstructionHooksForTesting() = {};
		observed += retained.count;
		retained.Release(comm, true);
		iga::CollectiveLocalStage(comm, "sequential retry validation", [&] {
			Check(retry == 0, "healthy retry failed");
			if (rank == 0) ValidateOutput(retry_output, graph);
		});
	}
	if (graph) for (const auto& scheme : {"fixed", "aitken"}) {
		iga::CollectiveLocalStage(comm, "sequential strong fixture", [&] {
			if (rank == 0) WriteGraphCase(fixture, scheme);
		});
		const auto output = fixture/scheme;
		const int status = Invoke(comm, Arguments(fixture, output, true));
		iga::CollectiveLocalStage(comm, "sequential strong validation", [&] {
			Check(status == 0, "healthy strong run failed");
			if (rank == 0) ValidateOutput(output, true, scheme);
		});
		iga::CollectiveLocalStage(comm, "restore sequential fixture", [&] {
			if (rank == 0) WriteGraphCase(fixture);
		});
	}
	if (ranks > 1) {
		int index = 0;
		for (const auto& option : {"-ksp_type", "-fieldsplit_0_pc_type", "-sequential_unused_probe"}) {
			std::array<char, 128> previous{};
			PetscBool present = PETSC_FALSE;
			iga::CollectiveLocalStage(comm, "sequential option fault setup", [&] {
				Check(PetscOptionsGetString(nullptr, nullptr, option, previous.data(), previous.size(), &present) == 0,
					"cannot save test option");
				if (rank == ranks-1)
					Check(PetscOptionsSetValue(nullptr, option, index == 0 ? "cg" : "jacobi") == 0,
						"cannot inject effective option difference");
			});
			bool built_native = false;
			iga::RuntimeConstructionHooksForTesting().after_local_stage = [&](const char* stage) {
				if (std::string_view(stage) == "sequential native inputs") built_native = true;
			};
			const auto failed_output = fixture/("option-failed-"+std::to_string(index));
			const int failed = Invoke(comm, Arguments(fixture, failed_output, graph));
			iga::RuntimeConstructionHooksForTesting() = {};
			iga::CollectiveLocalStage(comm, "sequential option fault restore", [&] {
				Check(PetscOptionsClearValue(nullptr, option) == 0, "cannot clear test option");
				if (present) Check(PetscOptionsSetValue(nullptr, option, previous.data()) == 0, "cannot restore test option");
				Check(failed == 1 && !built_native, "effective option difference entered native construction");
				if (rank == 0) Check(!fs::exists(failed_output), "option failure published output");
			});
			const auto retry_output = fixture/("option-retry-"+std::to_string(index));
			const int retry = Invoke(comm, Arguments(fixture, retry_output, graph));
			iga::CollectiveLocalStage(comm, "sequential option retry", [&] {
				Check(retry == 0, "restored option retry failed");
				if (rank == 0) ValidateOutput(retry_output, graph);
			});
			++index;
		}
		if (rank == 0) std::cout << "sequential options ranks=" << ranks << " graph=" << graph
			<< " faults=" << index << " retries=" << index << '\n';
	}
	if (rank == 0) std::cout << "sequential initialization ranks=" << ranks << " graph=" << graph
		<< " faults=" << stages.size() << " retries=" << stages.size()
		<< " observed_objects=" << observed << " retained_objects_released\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Check(argc == 2 && ranks == 3, "usage: mpiexec -np 3 sequential_initialization_failure_test NEW_OUTPUT_PARENT");
		const fs::path root(argv[1]);
		RunCases(PETSC_COMM_WORLD, root/"world", false);
		RunCases(PETSC_COMM_WORLD, root/"world", true);
		MPI_Comm subgroup = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &subgroup);
		const auto group_root = root/(rank == 0 ? "single" : "pair");
		RunCases(subgroup, group_root, false);
		RunCases(subgroup, group_root, true);
		MPI_Comm_free(&subgroup);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); return 2;
	}
}

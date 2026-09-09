#include "MultidomainRunner.hpp"
#include "RuntimeConstruction.hpp"
#include "SequentialExecution.hpp"
#include "OneDRuntime.hpp"
#include "TransientFlowRuntime.hpp"

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
	std::vector<std::string> args{"sequential_strong_failure_test"};
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

void ValidateOutput(const fs::path& output, bool graph, const std::string& scheme = "explicit",
	int first_step_veto_iteration = 0)
{
	if (scheme == "explicit") {
		ValidateRows(ReadHistory(output/"explicit_coupling_history.csv"));
		RequireManifest(output/"explicit_coupling_manifest.json");
	} else {
		ValidateStrongRun(ReadHistory(output/"strong_coupling_history.csv"),
			ReadHistory(output/"strong_coupling_iterations.csv"), output/"strong_coupling_manifest.json",
			scheme == "fixed", first_step_veto_iteration);
	}
	if (graph) RequireGraphBindingManifest(output/"graph_binding_manifest.json", scheme);
}

bool SameFlow(const iga::OneDFlowState& a, const iga::OneDFlowState& b)
{
	if (std::tie(a.area, a.flow, a.pressure, a.node_pressure, a.segment_flow, a.completed_step,
		a.physical_time, a.internal_substeps, a.inlet_flow)
		!= std::tie(b.area, b.flow, b.pressure, b.node_pressure, b.segment_flow, b.completed_step,
			b.physical_time, b.internal_substeps, b.inlet_flow) || a.outlets.size() != b.outlets.size()) return false;
	for (std::size_t i = 0; i < a.outlets.size(); ++i) {
		const auto& x = a.outlets[i]; const auto& y = b.outlets[i];
		if (std::tie(x.node, x.kind, x.pressure, x.resistance, x.proximal_resistance, x.distal_resistance,
			x.capacitance, x.reference_pressure, x.capacitor_pressure, x.flow)
			!= std::tie(y.node, y.kind, y.pressure, y.resistance, y.proximal_resistance, y.distal_resistance,
				y.capacitance, y.reference_pressure, y.capacitor_pressure, y.flow)) return false;
	}
	return true;
}

struct Snapshot {
	iga::OneDFlowState upstream, downstream;
	std::vector<PetscScalar> fluid;
	void Capture(const iga::OneDFlowRuntime& up, const iga::TransientFlowRuntime& flow, const iga::OneDFlowRuntime& down)
	{
		upstream = up.FlowState(); downstream = down.FlowState();
		PetscInt count = 0;
		Check(VecGetLocalSize(flow.State(), &count) == 0, "cannot inspect local state length");
		iga::PetscReadArray array; array.Acquire(flow.State());
		fluid.clear();
		if (count) fluid.assign(array.Data(), array.Data()+count);
		array.Restore();
	}
	void Verify(const iga::OneDFlowRuntime& up, const iga::TransientFlowRuntime& flow, const iga::OneDFlowRuntime& down) const
	{
		Check(up.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready
			&& down.CurrentPhase() == iga::OneDFlowRuntime::Phase::Ready
			&& flow.Phase() == iga::FlowStepPhase::Committed, "abort left a runtime trial open");
		Snapshot current; current.Capture(up, flow, down);
		Check(SameFlow(upstream, current.upstream) && SameFlow(downstream, current.downstream)
			&& fluid == current.fluid, "abort changed the committed physical state");
	}
};

void RunCases(MPI_Comm comm, const fs::path& parent, const std::string& scheme)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const auto fixture = parent/scheme;
	iga::CollectiveLocalStage(comm, "strong test fixture", [&] {
		if (rank != 0) return;
		Check(fs::create_directories(fixture), "fixture must be new");
		for (const auto& name : {"three_d", "upstream", "downstream"}) fs::create_directory(fixture/name);
		WriteThreeDCase(fixture/"three_d");
		WriteOneDCase(fixture/"upstream"); WriteOneDCase(fixture/"downstream");
		WriteUnitDatabase(fixture/"one.ntiga", ranks); WriteGraphCase(fixture, scheme);
	});
	const auto baseline = fixture/"baseline";
	const int base_status = Invoke(comm, Arguments(fixture, baseline, true));
	iga::CollectiveLocalStage(comm, "strong baseline validation", [&] {
		Check(base_status == 0, "healthy baseline failed");
		if (rank == 0) ValidateOutput(baseline, true, scheme);
	});
	std::vector<std::string> stages{"strong ready", "strong begin upstream", "strong begin 3D",
		"strong begin downstream", "strong upstream schedule", "strong upstream input", "strong solve upstream",
		"strong 3D input", "strong solve three_d", "strong measure inlet", "strong measure outlet",
		"strong measure wall", "strong downstream input", "strong solve downstream", "strong trial result",
		"strong convergence input", "strong precommit", "strong prepare upstream", "strong prepare three_d",
		"strong prepare downstream", "strong accepted result", "strong rollback downstream",
		"strong rollback three_d", "strong rollback upstream", "strong relaxation update", "strong step completion",
		"abort-errors"};
	std::size_t observed_objects = 0;
	for (std::size_t i = 0; i < stages.size(); ++i) {
		const bool abort_fault = stages[i] == "abort-errors";
		const auto primary = abort_fault ? "strong precommit" : stages[i];
		const bool postcommit = primary == "strong accepted result" || primary == "strong step completion";
		Snapshot ready, accepted;
		int ready_count = 0, accepted_count = 0, aborted_count = 0, abort_calls = 0, secondary_faults = 0;
		bool injected = false;
		Retained retained;
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		iga::SequentialExecutionHooksForTesting().observe = [&](const char* stage, int step,
			const iga::OneDFlowRuntime& up, const iga::TransientFlowRuntime& flow, const iga::OneDFlowRuntime& down) {
			Check(step == 1, "failure case advanced beyond the first physical step");
			if (std::string_view(stage) == "ready") { ready.Capture(up, flow, down); ++ready_count; }
			if (std::string_view(stage) == "accepted") { accepted.Capture(up, flow, down); ++accepted_count; }
			if (std::string_view(stage) == "aborted") {
				(postcommit ? accepted : ready).Verify(up, flow, down);
				++aborted_count;
			}
		};
		iga::RuntimeConstructionHooksForTesting().after_local_stage = [&](const char* stage) {
			if (std::string_view(stage).find("strong abort ") == 0) {
				++abort_calls;
				if (abort_fault && rank == ranks-1) { ++secondary_faults; throw std::runtime_error("injected abort outcome"); }
			}
			if (rank == ranks-1 && !injected && primary == stage) {
				injected = true;
				if (i == 5) throw 7;
				throw std::bad_alloc();
			}
		};
		const auto output = fixture/("failed-"+std::to_string(i));
		std::ostringstream diagnostic;
		auto* previous = std::cerr.rdbuf(diagnostic.rdbuf());
		const int failed = Invoke(comm, Arguments(fixture, output, true));
		std::cerr.rdbuf(previous);
		iga::RuntimeConstructionHooksForTesting() = {}; iga::SequentialExecutionHooksForTesting() = {};
		observed_objects += retained.count;
		retained.Release(comm, true);
		iga::CollectiveLocalStage(comm, "strong expected failure", [&] {
			Check(failed == 1 && ready_count == 1, "runner accepted fault or missed ready snapshot");
			if (rank == ranks-1) Check(injected, "target stage was not reached");
			Check(abort_calls == (primary == "strong ready" ? 0 : 3), "not all runtime aborts were attempted");
			Check(aborted_count == (primary == "strong ready" ? 0 : 1), "abort did not restore the expected physical snapshot");
			Check(accepted_count == (postcommit ? 1 : 0), "failed precommit partially committed a runtime");
			if (abort_fault && rank == ranks-1) Check(secondary_faults == 3, "not all cleanup faults were injected");
			if (rank == 0) {
				Check(!fs::exists(output), "failed strong step published output");
				if (abort_fault) for (const auto* part : {"primary:", "abort downstream:", "abort three_d:", "abort upstream:"})
					Check(diagnostic.str().find(part) != std::string::npos, "missing primary or cleanup diagnostic");
				std::ofstream log(fixture/("failed-"+std::to_string(i)+".stderr")); log << diagnostic.str(); log.close();
				Check(bool(log), "cannot save failure diagnostic");
			}
		});
		const auto retry = fixture/("retry-"+std::to_string(i));
		iga::RuntimeConstructionHooksForTesting().object_created = [&](PetscObject object) { retained.Add(object); };
		const int retry_status = Invoke(comm, Arguments(fixture, retry, true));
		iga::RuntimeConstructionHooksForTesting() = {};
		observed_objects += retained.count;
		retained.Release(comm, true);
		iga::CollectiveLocalStage(comm, "strong retry comparison", [&] {
			Check(retry_status == 0, "healthy retry failed");
			if (rank == 0) for (const auto& file : {"strong_coupling_history.csv", "strong_coupling_iterations.csv"})
				Check(ReadText(retry/file) == ReadText(baseline/file), "retry differs from healthy baseline");
		});
	}
	int first_converged = 0, current_iteration = 0, accepted_iteration = 0;
	bool vetoed = false;
	iga::SequentialExecutionHooksForTesting().convergence = [&](int step, int iteration, bool& converged) {
		if (step != 1) return;
		current_iteration = iteration;
		if (converged && !first_converged) {
			first_converged = iteration;
			if (rank == ranks-1) { converged = false; vetoed = true; }
		}
	};
	iga::SequentialExecutionHooksForTesting().observe = [&](const char* stage, int step,
		const iga::OneDFlowRuntime&, const iga::TransientFlowRuntime&, const iga::OneDFlowRuntime&) {
		if (step == 1 && std::string_view(stage) == "accepted") accepted_iteration = current_iteration;
	};
	const auto veto_output = fixture/"veto";
	const int veto_status = Invoke(comm, Arguments(fixture, veto_output, true));
	iga::SequentialExecutionHooksForTesting() = {};
	iga::CollectiveLocalStage(comm, "strong convergence veto", [&] {
		Check(veto_status == 0 && first_converged > 0 && accepted_iteration > first_converged,
			"a rank committed before the global convergence decision");
		if (rank == ranks-1) Check(vetoed, "convergence veto did not run");
		if (rank == 0) ValidateOutput(veto_output, true, scheme, first_converged);
	});
	iga::RequireCollectiveSameInt(comm, "strong accepted iteration agreement", accepted_iteration);
	if (rank == 0) std::cout << "sequential strong ranks=" << ranks << " scheme=" << scheme
		<< " faults=" << stages.size() << " retries=" << stages.size() << " veto=1"
		<< " observed_objects=" << observed_objects << " snapshots_restored retained_objects_released\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Check(argc == 2 && ranks == 3, "usage: mpiexec -np 3 sequential_strong_failure_test NEW_OUTPUT_PARENT");
		const fs::path root(argv[1]);
		for (const auto& scheme : {"fixed", "aitken"}) RunCases(PETSC_COMM_WORLD, root/"world", scheme);
		MPI_Comm subgroup = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &subgroup);
		const auto group_root = root/(rank == 0 ? "single" : "pair");
		for (const auto& scheme : {"fixed", "aitken"}) RunCases(subgroup, group_root, scheme);
		MPI_Comm_free(&subgroup);
		PetscFinalize(); return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; MPI_Abort(PETSC_COMM_WORLD, 2); return 2;
	}
}

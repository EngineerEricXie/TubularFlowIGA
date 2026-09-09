#include "CollectivePressureFlowExecution.hpp"
#include "AbortReportAllocation.hpp"

// Reuse the numerical fixture and exact result comparator from the serial
// contract tests, including trial/prepare/finalize state tracking.
#define main PressureFlowSerialContractMain
#include "test_pressure_flow_executor.cpp"
#undef main

namespace {

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

void TestAbortReporting(MPI_Comm comm)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const auto graph = Chain();
	const iga::DomainStepContext step{0, 0.0, 0.01};
	const std::map<std::string, double> initial{{"left", 10.0}, {"right", 20.0}};
	std::size_t total = 0;
	for (std::size_t ordinal = 0; ordinal <= total; ++ordinal) {
		Fixture fixture(graph);
		auto policy = iga::CollectivePressureFlowExecution(comm);
		const auto common = policy.outcome;
		int aborts = 0, reports = 0;
		std::size_t allocations = 0;
		bool report_error = false;
		policy.outcome = [&](const char* stage, std::exception_ptr error) {
			if (std::strcmp(stage, "pressure-flow abort") == 0) {
				++aborts;
				try { common(stage, error); }
				catch (...) {
					if (aborts == 3 && rank == ranks-1) abort_report_allocation::Arm(ordinal);
					throw;
				}
				if (aborts == 3 && rank == ranks-1) abort_report_allocation::Arm(ordinal);
				return;
			}
			if (std::strcmp(stage, "pressure-flow abort reporting") == 0) {
				allocations = abort_report_allocation::allocations;
				report_error = bool(error); ++reports;
				abort_report_allocation::Disarm();
			}
			common(stage, error);
		};
		if (rank == ranks-1) {
			fixture.down->fail_solve = true;
			fixture.up->fail_abort = fixture.mid->fail_abort = true;
		}
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", {}, policy);
		std::string diagnostic;
		try { executor.Advance(step, initial); }
		catch (const std::exception& error) { diagnostic = error.what(); }
		abort_report_allocation::Disarm();
		iga::CollectiveLocalStage(comm, "pressure abort report checks", [&] {
			Require(aborts == 3 && reports == 1, "report ran before all abort outcomes");
			Require(fixture.up->aborts == 1 && fixture.mid->aborts == 1 && fixture.down->aborts == 1,
				"formatting failure skipped a runtime abort");
			Require(fixture.up->commits == 0 && fixture.mid->commits == 0 && fixture.down->commits == 0,
				"failed cleanup committed a domain");
			Require(diagnostic.find("injected fake solve failure") != std::string::npos, "primary failure was lost");
			if (ordinal) {
				Require(abort_report_allocation::injected == (rank == ranks-1), "missed report allocation");
				Require(report_error == (rank == ranks-1), "formatter swallowed allocation failure");
				Require(diagnostic.find("cleanup reporting failed: pressure-flow abort reporting") != std::string::npos,
					"missing coordinated reporting failure");
			} else {
				Require(diagnostic.find("up:") != std::string::npos && diagnostic.find("mid:") != std::string::npos,
					"healthy diagnostic lost an abort failure");
			}
		});
		iga::RequireCollectiveSameText(comm, "pressure abort report agreement", diagnostic);
		if (!ordinal) {
			unsigned long long local = rank == ranks-1 ? allocations : 0, maximum = 0;
			MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm);
			Require(maximum > 0 && maximum < 64, "unexpected report allocation count");
			total = static_cast<std::size_t>(maximum);
		}
		fixture.down->fail_solve = fixture.up->fail_abort = fixture.mid->fail_abort = false;
		const auto retry = executor.Advance(step, initial);
		Fixture reference(graph);
		iga::PressureFlowComponentExecutor serial(*reference.registry, "up", {});
		RequireSameStepResult(retry, serial.Advance(step, initial));
		Require(fixture.up->commits == 1 && fixture.mid->commits == 1 && fixture.down->commits == 1,
			"retry did not commit after reporting failure");
	}
	if (rank == 0) std::cout << "pressure abort reporting ranks=" << ranks
		<< " allocation_faults=" << total << " retries=" << total+1 << " passed\n";
}

void RunCases(MPI_Comm communicator, int seed)
{
	TestAbortReporting(communicator);
	int rank = 0, ranks = 0;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	const int failed_rank = ranks-1;
	const auto graph = Chain();
	const iga::DomainStepContext step{0, 0.0, 0.01};
	const std::map<std::string, double> initial{{"left", 10.0}, {"right", 20.0}};
	const std::vector<std::string> faults{
		"input", "solve", "port", "prepare", "callback", "abort",
		"pressure-flow begin", "pressure-flow pressure input",
		"pressure-flow flow input", "pressure-flow provider validation",
		"pressure-flow port state", "pressure-flow port validation",
		"pressure-flow edge result", "pressure-flow iteration result",
		"pressure-flow relaxation", "pressure-flow rollback",
		"pressure-flow observation reference", "pressure-flow accepted port",
		"pressure-flow before commit", "pressure-flow prepare commit"};
	int retries = 0;
	for (const auto& fault : faults) {
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		controls.method = iga::PressureFlowIterationMethod::Fixed;
		controls.maximum_iterations = 4;
		auto policy = iga::CollectivePressureFlowExecution(communicator);
		const auto common_outcome = policy.outcome;
		bool injected = false;
		bool armed = true;
		policy.outcome = [&](const char* stage, std::exception_ptr error) {
			// The policy sees an outcome only after the local operation returns.
			// Inject a caught allocation exception at that boundary; no production
			// allocator or fault-injection hook is changed.
			if (armed && !injected && rank == failed_rank && fault == stage) {
				try { throw std::bad_alloc(); }
				catch (...) { error = std::current_exception(); }
				injected = true;
			}
			iga::RequireCollectiveSameText(communicator, "executor stage order", stage);
			common_outcome(stage, error);
		};
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls, policy);
		auto pressure = initial;
		if (fault == "pressure-flow relaxation" || fault == "pressure-flow rollback")
			pressure = {{"left", 0.0}, {"right", 0.0}};
		if (rank == failed_rank) {
			if (fault == "input") pressure.erase("left");
			if (fault == "solve" || fault == "abort") fixture.down->fail_solve = true;
			if (fault == "port") fixture.up->nonfinite_state = true;
			if (fault == "prepare") fixture.down->fail_prepare = true;
			if (fault == "abort") {
				fixture.mid->fail_abort = true;
				fixture.up->fail_abort = true;
			}
		}
		std::string diagnostic;
		try {
			executor.Advance(step, pressure, [&](const iga::PressureFlowStepResult&) {
				if (rank == failed_rank && fault == "callback")
					throw std::runtime_error("local precommit observation failed");
			});
		} catch (const std::exception& error) { diagnostic = error.what(); }
		iga::CollectiveLocalStage(communicator, "test expected rejection", [&] {
			Require(!diagnostic.empty(), "fault was accepted");
			Require(fixture.up->commits == 0 && fixture.mid->commits == 0
				&& fixture.down->commits == 0, "failed trial was committed");
			if (fault == "abort") {
				Require(fixture.up->aborts == 1 && fixture.mid->aborts == 1
					&& fixture.down->aborts == 1, "cleanup skipped a runtime");
				Require(diagnostic.find("injected fake solve failure") != std::string::npos
					&& diagnostic.find("up:") != std::string::npos
					&& diagnostic.find("mid:") != std::string::npos,
					"cleanup lost the primary or an abort failure");
			}
		});
		iga::RequireCollectiveSameText(communicator, "test common diagnostic", diagnostic);
		armed = false;
		fixture.up->nonfinite_state = false;
		fixture.down->fail_solve = false;
		fixture.down->fail_prepare = false;
		fixture.mid->fail_abort = false;
		fixture.up->fail_abort = false;
		const auto retry = executor.Advance(step, initial);
		Fixture reference(graph);
		iga::PressureFlowComponentExecutor serial(*reference.registry, "up", controls);
		RequireSameStepResult(retry, serial.Advance(step, initial));
		Require(fixture.up->commits == 1 && fixture.mid->commits == 1
			&& fixture.down->commits == 1, "healthy retry did not commit every domain");
		++retries;
	}
	// One participant rejects the first convergence decision. Everyone must
	// perform the same second sweep even though the local fixed point is exact.
	for (const auto method : {iga::PressureFlowIterationMethod::Fixed,
		iga::PressureFlowIterationMethod::Aitken}) {
		Fixture fixture(graph);
		iga::PressureFlowExecutionControls controls;
		controls.method = method;
		controls.maximum_iterations = 3;
		auto policy = iga::CollectivePressureFlowExecution(communicator);
		const auto all_converged = policy.all_converged;
		int decisions = 0;
		policy.all_converged = [&](bool value) {
			++decisions;
			return all_converged(value && !(rank == failed_rank && decisions == 1));
		};
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", controls, policy);
		const auto result = executor.Advance(step, initial);
		Require(result.iterations.size() == 2 && result.iterations.back().converged,
			"participants did not share the convergence decision");
		Require(fixture.up->rollbacks == 1 && fixture.mid->rollbacks == 1
			&& fixture.down->rollbacks == 1, "global rejection did not roll back every domain");
	}
	// Extra work in the singleton subgroup catches accidental COMM_WORLD use.
	for (int i = 0; i < seed; ++i) {
		Fixture fixture(graph);
		iga::PressureFlowComponentExecutor executor(*fixture.registry, "up", {},
			iga::CollectivePressureFlowExecution(communicator));
		Require(executor.Advance(step, initial).iterations.size() == 1,
			"explicit subgroup execution failed");
	}
	if (rank == 0) std::cout << "pressure executor ranks=" << ranks << " seed=" << seed
		<< " faults=" << faults.size() << " retries=" << retries
		<< " global_convergence=2\n";
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
		Require(ranks == 3, "run with exactly three ranks");
		RunCases(PETSC_COMM_WORLD, 0);
		MPI_Comm subgroup = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &subgroup);
		RunCases(subgroup, rank == 0 ? 2 : 1);
		MPI_Comm_free(&subgroup);
		PetscFinalize();
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
		return 2;
	}
}

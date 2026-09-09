#include "CollectiveSpeciesPressureFlowExecution.hpp"
#include "AbortReportAllocation.hpp"

#define main SpeciesSerialContractMain
#include "test_species_pressure_flow_executor.cpp"
#undef main

namespace {

void Require(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

void SamePorts(const std::map<iga::PortRef, iga::PortState>& first,
	const std::map<iga::PortRef, iga::PortState>& second)
{
	Require(first.size() == second.size(), "port count changed");
	for (const auto& item : first) {
		const auto& a = item.second;
		const auto& b = second.at(item.first);
		Require(a.time_s == b.time_s && a.area_m2 == b.area_m2
			&& a.outward_flow_m3_s == b.outward_flow_m3_s
			&& a.mean_pressure_pa == b.mean_pressure_pa
			&& a.mean_normal_traction_pa == b.mean_normal_traction_pa
			&& a.total_pressure_pa == b.total_pressure_pa
			&& a.concentration == b.concentration
			&& a.outward_species_flux == b.outward_species_flux, "port field changed");
	}
}

void SameResult(const iga::SpeciesPressureFlowStepResult& a,
	const iga::SpeciesPressureFlowStepResult& b)
{
	SamePorts(a.accepted_ports, b.accepted_ports);
	SamePorts(a.transport_ports, b.transport_ports);
	Require(a.donor_ownership == b.donor_ownership
		&& a.transport_domain_order == b.transport_domain_order, "routing changed");
	Require(a.hydraulic_iterations.size() == b.hydraulic_iterations.size(), "hydraulic work changed");
	Require(a.edge_amounts.size() == b.edge_amounts.size()
		&& a.domain_balances.size() == b.domain_balances.size()
		&& a.global_balances.size() == b.global_balances.size(), "accounting coverage changed");
	for (std::size_t i = 0; i < a.edge_amounts.size(); ++i) {
		const auto& x = a.edge_amounts[i]; const auto& y = b.edge_amounts[i];
		Require(x.edge_id == y.edge_id && x.species_id == y.species_id
			&& x.first == y.first && x.second == y.second && x.donor == y.donor
			&& x.first_outward_amount == y.first_outward_amount
			&& x.second_outward_amount == y.second_outward_amount
			&& x.residual == y.residual && x.normalized_residual == y.normalized_residual,
			"edge amount changed");
	}
	for (std::size_t i = 0; i < a.domain_balances.size(); ++i) {
		const auto& x = a.domain_balances[i]; const auto& y = b.domain_balances[i];
		Require(x.domain_id == y.domain_id && x.species_id == y.species_id
			&& x.accounting.initial_mass == y.accounting.initial_mass
			&& x.accounting.final_mass == y.accounting.final_mass
			&& x.accounting.source_amount == y.accounting.source_amount
			&& x.accounting.outward_port_amount == y.accounting.outward_port_amount
			&& x.accounting.residual == y.accounting.residual
			&& x.recomputed_residual == y.recomputed_residual
			&& x.normalized_residual == y.normalized_residual, "domain amount changed");
	}
	for (const auto& item : a.global_balances) {
		const auto& x = item.second; const auto& y = b.global_balances.at(item.first);
		Require(x.initial_mass == y.initial_mass && x.final_mass == y.final_mass
			&& x.source_amount == y.source_amount && x.outward_amount == y.outward_amount
			&& x.residual == y.residual && x.gross_activity == y.gross_activity
			&& x.normalized_residual == y.normalized_residual, "global amount changed");
	}
}

void TestAbortReporting(MPI_Comm comm)
{
	int rank = 0, ranks = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	const iga::DomainStepContext step{0, 0.0, 0.1};
	const std::map<std::string, double> initial{{"edge", 1.0}};
	std::size_t total = 0;
	for (std::size_t ordinal = 0; ordinal <= total; ++ordinal) {
		Fixture fixture;
		auto policy = iga::CollectiveSpeciesPressureFlowExecution(comm);
		const auto common = policy.execution.outcome;
		int aborts = 0, reports = 0;
		std::size_t allocations = 0;
		bool report_error = false;
		policy.execution.outcome = [&](const char* stage, std::exception_ptr error) {
			if (std::strcmp(stage, "species abort") == 0) {
				++aborts;
				try { common(stage, error); }
				catch (...) {
					if (aborts == 2 && rank == ranks-1) abort_report_allocation::Arm(ordinal);
					throw;
				}
				if (aborts == 2 && rank == ranks-1) abort_report_allocation::Arm(ordinal);
				return;
			}
			if (std::strcmp(stage, "species abort reporting") == 0) {
				allocations = abort_report_allocation::allocations;
				report_error = bool(error); ++reports;
				abort_report_allocation::Disarm();
			}
			common(stage, error);
		};
		if (rank == ranks-1) {
			fixture.b->fail_transport = true;
			fixture.a->fail_abort = fixture.b->fail_abort = true;
		}
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls(), policy);
		std::string diagnostic;
		try { executor.Advance(step, initial); }
		catch (const std::exception& error) { diagnostic = error.what(); }
		abort_report_allocation::Disarm();
		iga::CollectiveLocalStage(comm, "species abort report checks", [&] {
			Require(aborts == 2 && reports == 1, "report ran before all abort outcomes");
			Require(fixture.a->aborts == 1 && fixture.b->aborts == 1,
				"formatting failure skipped a runtime abort");
			Require(fixture.a->commits == 0 && fixture.b->commits == 0,
				"failed cleanup committed a domain");
			Require(diagnostic.find("injected transport failure") != std::string::npos, "primary failure was lost");
			if (ordinal) {
				Require(abort_report_allocation::injected == (rank == ranks-1), "missed report allocation");
				Require(report_error == (rank == ranks-1), "formatter swallowed allocation failure");
				Require(diagnostic.find("cleanup reporting failed: species abort reporting") != std::string::npos,
					"missing coordinated reporting failure");
			} else {
				Require(diagnostic.find("a:") != std::string::npos && diagnostic.find("b:") != std::string::npos,
					"healthy diagnostic lost an abort failure");
			}
		});
		iga::RequireCollectiveSameText(comm, "species abort report agreement", diagnostic);
		if (!ordinal) {
			unsigned long long local = rank == ranks-1 ? allocations : 0, maximum = 0;
			MPI_Allreduce(&local, &maximum, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm);
			Require(maximum > 0 && maximum < 64, "unexpected report allocation count");
			total = static_cast<std::size_t>(maximum);
		}
		fixture.b->fail_transport = fixture.a->fail_abort = fixture.b->fail_abort = false;
		const auto retry = executor.Advance(step, initial);
		Fixture reference;
		iga::SpeciesPressureFlowComponentExecutor serial(*reference.registry, "a", Controls());
		SameResult(retry, serial.Advance(step, initial));
		Require(executor.CommittedDonorOwnership() == retry.donor_ownership, "retry did not restore donor ownership");
		Require(fixture.a->commits == 1 && fixture.b->commits == 1,
			"retry did not commit after reporting failure");
	}
	if (rank == 0) std::cout << "species abort reporting ranks=" << ranks
		<< " allocation_faults=" << total << " retries=" << total+1 << " passed\n";
}


void RunCases(MPI_Comm communicator, int seed)
{
	TestAbortReporting(communicator);
	int rank = 0, ranks = 0;
	MPI_Comm_rank(communicator, &rank); MPI_Comm_size(communicator, &ranks);
	const int failed_rank = ranks-1;
	const iga::DomainStepContext step{0, 0.0, 0.1};
	const std::map<std::string, double> initial{{"edge", 1.0}};
	std::vector<std::string> faults{
		"input", "transport", "prepare", "callback", "nonstandard", "accounting", "abort", "port",
		"species step input", "species begin", "species pressure input", "species hydraulic solve",
		"species provider validation", "species flow input", "species hydraulic port",
		"species hydraulic port validation", "species hydraulic edge result",
		"species hydraulic iteration result", "species hydraulic relaxation", "species hydraulic rollback",
		"species hydraulic observation reference", "species accepted hydraulic port",
		"species transport schedule", "species donor port", "species concentration input",
		"species set concentration", "species transport solve", "species native accounting",
		"species accounting catalog", "species amount verification",
		"species transport observation reference", "species final port observation",
		"species accepted transport port", "species before commit", "species prepare commit"};
	if (ranks > 1) faults.push_back("schedule");
	for (const auto& fault : faults) {
		Fixture fixture;
		auto controls = Controls();
		controls.hydraulic.method = iga::PressureFlowIterationMethod::Fixed;
		controls.hydraulic.maximum_iterations = 4;
		auto policy = iga::CollectiveSpeciesPressureFlowExecution(communicator);
		const auto common_outcome = policy.execution.outcome;
		bool armed = true, injected = false;
		policy.execution.outcome = [&](const char* stage, std::exception_ptr error) {
			if (armed && !injected && rank == failed_rank && fault == stage) {
				try { throw std::bad_alloc(); }
				catch (...) { error = std::current_exception(); }
				injected = true;
			}
			iga::RequireCollectiveSameText(communicator, "species stage order", stage);
			common_outcome(stage, error);
		};
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", controls, policy);
		auto pressure = initial;
		if (fault == "species hydraulic relaxation" || fault == "species hydraulic rollback")
			pressure.at("edge") = 0.0;
		if (rank == failed_rank) {
			if (fault == "input") pressure.clear();
			if (fault == "transport" || fault == "abort") fixture.b->fail_transport = true;
			if (fault == "prepare") fixture.b->fail_prepare = true;
			if (fault == "abort") fixture.a->fail_abort = fixture.b->fail_abort = true;
			if (fault == "accounting") fixture.b->corrupt_port_amounts = true;
			if (fault == "port") fixture.a->signed_flow = std::numeric_limits<double>::quiet_NaN();
			if (fault == "schedule") fixture.a->signed_flow = fixture.b->signed_flow = -1.0;
		}
		std::string diagnostic;
		try {
			executor.Advance(step, pressure, [&](const iga::SpeciesPressureFlowStepResult&) {
				if (rank == failed_rank && fault == "callback") throw std::runtime_error("species observation failed");
				if (rank == failed_rank && fault == "nonstandard") throw 42;
			});
		} catch (const std::exception& error) { diagnostic = error.what(); }
		iga::CollectiveLocalStage(communicator, "species expected rejection", [&] {
			Require(!diagnostic.empty(), "fault was accepted");
			Require(fixture.a->commits == 0 && fixture.b->commits == 0, "failed trial was committed");
			Require(executor.CommittedDonorOwnership().empty(), "failed trial committed donor ownership");
			if (fault == "schedule") {
				Require(diagnostic.find("species transport schedule agreement") != std::string::npos,
					"schedule divergence was rejected at the wrong stage");
				Require(fixture.a->transport_solves == 0 && fixture.b->transport_solves == 0,
					"inconsistent schedule reached a transport solve");
			}
			if (fault == "abort") {
				Require(fixture.a->aborts == 1 && fixture.b->aborts == 1, "abort skipped a domain");
				Require(diagnostic.find("injected transport failure") != std::string::npos
					&& diagnostic.find("a:") != std::string::npos
					&& diagnostic.find("b:") != std::string::npos, "cleanup lost failures");
			}
		});
		iga::RequireCollectiveSameText(communicator, "species common diagnostic", diagnostic);
		armed = false;
		fixture.b->fail_transport = fixture.b->fail_prepare = false;
		fixture.a->fail_abort = fixture.b->fail_abort = false;
		fixture.b->corrupt_port_amounts = false;
		fixture.a->signed_flow = fixture.b->signed_flow = 1.0;
		const auto retry = executor.Advance(step, initial);
		Fixture reference;
		iga::SpeciesPressureFlowComponentExecutor serial(*reference.registry, "a", controls);
		SameResult(retry, serial.Advance(step, initial));
		Require(fixture.a->commits == 1 && fixture.b->commits == 1, "retry did not commit all domains");
		Require(executor.CommittedDonorOwnership() == retry.donor_ownership, "retry donor ownership missing");
	}
	for (const auto method : {iga::PressureFlowIterationMethod::Fixed, iga::PressureFlowIterationMethod::Aitken}) {
		Fixture fixture;
		auto controls = Controls(); controls.hydraulic.method = method;
		controls.hydraulic.maximum_iterations = 3;
		auto policy = iga::CollectiveSpeciesPressureFlowExecution(communicator);
		const auto all_converged = policy.execution.all_converged;
		int decisions = 0;
		policy.execution.all_converged = [&](bool value) {
			++decisions;
			return all_converged(value && !(rank == failed_rank && decisions == 1));
		};
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", controls, policy);
		const auto result = executor.Advance(step, initial);
		Require(result.hydraulic_iterations.size() == 2 && result.hydraulic_iterations.back().converged,
			"hydraulic convergence decision differs");
		Require(fixture.a->hydraulic_rollbacks == 1 && fixture.b->hydraulic_rollbacks == 1,
			"global rejection skipped a rollback");
	}
	for (int i = 0; i <= seed; ++i) {
		Fixture fixture;
		fixture.a->signed_flow = fixture.b->signed_flow = -1.0;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls(),
			iga::CollectiveSpeciesPressureFlowExecution(communicator));
		const auto result = executor.Advance(step, initial);
		Require(result.transport_domain_order == std::vector<std::string>{"b", "a"}
			&& fixture.a->last_received == 5.0 && fixture.a->last_received_dye == 15.0,
			"legitimate reverse flow did not use the reverse donor");
	}
	if (rank == 0) std::cout << "species executor ranks=" << ranks << " seed=" << seed
		<< " faults=" << faults.size() << " retries=" << faults.size()
		<< " convergence=2 reverse=" << seed+1 << '\n';
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	try {
		int rank = 0, ranks = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
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

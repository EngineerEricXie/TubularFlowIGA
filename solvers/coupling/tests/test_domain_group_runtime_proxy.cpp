#include "DomainGroupRuntimeProxy.hpp"

#include <petscsys.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

class Runtime final : public iga::CoupledDomainRuntime,
	public iga::StagedFlowTransportDomainRuntime {
public:
	Runtime(std::string id, MPI_Comm communicator, bool fail)
		: id_(std::move(id)), communicator_(communicator), fail_(fail) {}

	const std::string& DomainId() const noexcept override { return id_; }
	iga::DomainKind Kind() const noexcept override { return iga::DomainKind::OneDFlow; }
	const std::vector<iga::CouplingPort>& Ports() const noexcept override { return ports_; }
	void BeginStep(const iga::DomainStepContext& step) override { time_ = step.EndTime(); }
	void SetPortInput(const std::string&, const iga::PortBoundaryData& input) override
	{
		pressure_ = input.mean_pressure_pa.value_or(0.0);
	}
	void SolveTrial() override { Solve(false); }
	iga::PortState GetPortState(const std::string&) const override { return State(); }
	void RollbackTrial() override {}
	void AbortStep() override {}
	void PrepareCommitStep() override {}
	void FinalizeCommitStep() noexcept override {}
	void SolveHydraulicTrial() override { Solve(false); }
	iga::PortState GetHydraulicPortState(const std::string&) const override { return State(); }
	void RollbackHydraulicTrial() override {}
	void SetTransportConcentration(const std::string&, double,
		const std::map<std::string, double>& concentration) override
	{
		concentration_ = concentration;
	}
	void SolveTransportTrial() override { Solve(true); }
	iga::PortState GetTransportPortState(const std::string&) const override { return State(); }
	void RollbackTransportTrial() override {}
	std::map<std::string, iga::SpeciesStepAccounting> GetSpeciesStepAccounting() const override
	{
		iga::SpeciesStepAccounting value;
		value.initial_mass = sum_;
		value.final_mass = sum_+1.0;
		value.outward_port_amount["port"] = -1.0;
		return {{"tracer", value}};
	}

private:
	void Solve(bool transport)
	{
		int rank = 0;
		MPI_Comm_rank(communicator_, &rank);
		double local = rank+1.0;
		MPI_Allreduce(&local, &sum_, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		if (fail_ && !transport) {
			fail_ = false;
			throw std::runtime_error("injected group failure");
		}
	}
	iga::PortState State() const
	{
		iga::PortState state;
		state.time_s = time_;
		state.outward_flow_m3_s = sum_+pressure_;
		state.mean_pressure_pa = pressure_;
		state.concentration = concentration_;
		return state;
	}

	std::string id_;
	MPI_Comm communicator_;
	bool fail_ = false;
	double time_ = 0.0;
	double pressure_ = 0.0;
	double sum_ = 0.0;
	std::map<std::string, double> concentration_;
	std::vector<iga::CouplingPort> ports_;
};

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		if (ranks != 3) throw std::runtime_error("domain group proxy test requires three ranks");
		const int color = rank < 2 ? 0 : 1;
		MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
		auto local_a = rank < 2
			? std::unique_ptr<iga::CoupledDomainRuntime>(new Runtime("a", group, true)) : nullptr;
		auto local_b = rank == 2
			? std::unique_ptr<iga::CoupledDomainRuntime>(new Runtime("b", group, false)) : nullptr;
		iga::StagedDomainGroupRuntimeProxy a(PETSC_COMM_WORLD, 0, "a",
			iga::DomainKind::OneDFlow, {}, std::move(local_a));
		iga::StagedDomainGroupRuntimeProxy b(PETSC_COMM_WORLD, 2, "b",
			iga::DomainKind::OneDFlow, {}, std::move(local_b));
		const iga::DomainStepContext step{0, 0.0, 0.25};
		a.BeginStep(step);
		b.BeginStep(step);
		iga::PortBoundaryData input;
		input.time_s = step.EndTime();
		input.mean_pressure_pa = 7.0;
		a.SetPortInput("port", input);
		b.SetPortInput("port", input);
		bool rejected = false;
		std::vector<iga::CoupledDomainRuntime*> batch{&a, &b};
		try { iga::SolveDomainGroupBatch(PETSC_COMM_WORLD, batch,
			iga::DomainGroupSolveKind::Hydraulic); }
		catch (const std::exception& error) {
			rejected = std::string(error.what()).find("injected group failure") != std::string::npos;
		}
		if (!rejected) throw std::runtime_error("group-local failure was not propagated globally");
		iga::SolveDomainGroupBatch(PETSC_COMM_WORLD, batch,
			iga::DomainGroupSolveKind::Hydraulic);
		const auto state_a = a.GetHydraulicPortState("port");
		const auto state_b = b.GetHydraulicPortState("port");
		if (state_a.outward_flow_m3_s != 10.0 || state_b.outward_flow_m3_s != 8.0)
			throw std::runtime_error("cross-group hydraulic state differs");
		a.SetTransportConcentration("port", step.EndTime(), {{"tracer", 2.5}});
		iga::SolveDomainGroupBatch(PETSC_COMM_WORLD, batch,
			iga::DomainGroupSolveKind::Transport);
		const auto transport = a.GetTransportPortState("port");
		if (transport.concentration.at("tracer") != 2.5)
			throw std::runtime_error("cross-group concentration differs");
		const auto accounting = a.GetSpeciesStepAccounting().at("tracer");
		if (accounting.initial_mass != 3.0 || accounting.final_mass != 4.0
			|| accounting.outward_port_amount.at("port") != -1.0)
			throw std::runtime_error("cross-group accounting differs");
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "domain group runtime proxy tests passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		if (group != MPI_COMM_NULL) MPI_Comm_free(&group);
		MPI_Abort(PETSC_COMM_WORLD, 1);
	}
	PetscFinalize();
	return 0;
}

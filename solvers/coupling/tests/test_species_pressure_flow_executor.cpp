#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "PressureFlowCheckpointControls.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

void RequireRejected(const std::function<void()>& operation)
{
	bool rejected = false;
	try { operation(); }
	catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

iga::CouplingPort Port(const std::string& domain, const std::string& id,
	std::set<iga::PortQuantity> hydraulic_requires)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = domain;
	port.locator_kind = "fake";
	port.locator = id;
	port.provides = {iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure,
		iga::PortQuantity::SpeciesConcentration, iga::PortQuantity::SpeciesFlux};
	port.requires = std::move(hydraulic_requires);
	port.requires.insert(iga::PortQuantity::SpeciesConcentration);
	port.requires.insert(iga::PortQuantity::SpeciesFlux);
	port.species = {"dye", "tracer"};
	return port;
}

iga::SimulationGraph Graph()
{
	auto a = iga::DomainNode{"a", iga::DomainKind::OneDFlow,
		{Port("a", "p", {iga::PortQuantity::MeanPressure})},
		{{"dye", "native_dye_a"}, {"tracer", "native_a"}}};
	auto b = iga::DomainNode{"b", iga::DomainKind::ThreeDBodyFittedFlow,
		{Port("b", "p", {iga::PortQuantity::FlowRate})},
		{{"dye", "native_dye_b"}, {"tracer", "native_b"}}};
	return iga::SimulationGraph({std::move(a), std::move(b)},
		{{"edge", {"a", "p"}, {"b", "p"}, iga::CouplingLaw::PressureFlow, {"dye", "tracer"}}},
		{{"dye", "mol/m3"}, {"tracer", "mol/m3"}});
}

iga::SimulationGraph CycleGraph()
{
	auto a = iga::DomainNode{"a", iga::DomainKind::OneDFlow,
		{Port("a", "ab", {iga::PortQuantity::MeanPressure}),
			Port("a", "ad", {iga::PortQuantity::MeanPressure})},
		{{"dye", "a_dye"}, {"tracer", "a_tracer"}}};
	auto b = iga::DomainNode{"b", iga::DomainKind::ThreeDBodyFittedFlow,
		{Port("b", "ab", {iga::PortQuantity::FlowRate}),
			Port("b", "bc", {iga::PortQuantity::MeanPressure})},
		{{"dye", "b_dye"}, {"tracer", "b_tracer"}}};
	auto c = iga::DomainNode{"c", iga::DomainKind::OneDFlow,
		{Port("c", "bc", {iga::PortQuantity::FlowRate}),
			Port("c", "cd", {iga::PortQuantity::MeanPressure})},
		{{"dye", "c_dye"}, {"tracer", "c_tracer"}}};
	auto d = iga::DomainNode{"d", iga::DomainKind::ThreeDBodyFittedFlow,
		{Port("d", "cd", {iga::PortQuantity::FlowRate}),
			Port("d", "ad", {iga::PortQuantity::FlowRate})},
		{{"dye", "d_dye"}, {"tracer", "d_tracer"}}};
	return iga::SimulationGraph({std::move(a), std::move(b), std::move(c), std::move(d)}, {
		{"ab", {"a", "ab"}, {"b", "ab"}, iga::CouplingLaw::PressureFlow, {"dye", "tracer"}},
		{"bc", {"b", "bc"}, {"c", "bc"}, iga::CouplingLaw::PressureFlow, {"dye", "tracer"}},
		{"cd", {"c", "cd"}, {"d", "cd"}, iga::CouplingLaw::PressureFlow, {"dye", "tracer"}},
		{"ad", {"a", "ad"}, {"d", "ad"}, iga::CouplingLaw::PressureFlow, {"dye", "tracer"}}},
		{{"dye", "mol/m3"}, {"tracer", "mol/m3"}});
}

class FakeStagedRuntime : public iga::CoupledDomainRuntime,
	public iga::StagedFlowTransportDomainRuntime {
public:
	FakeStagedRuntime(std::string id, iga::DomainKind kind, std::vector<iga::CouplingPort> ports,
		double concentration)
		: id_(std::move(id)), kind_(kind), ports_(std::move(ports))
	{
		concentration_.emplace("tracer", concentration);
		concentration_.emplace("dye", concentration+10.0);
	}

	const std::string& DomainId() const noexcept override { return id_; }
	iga::DomainKind Kind() const noexcept override { return kind_; }
	const std::vector<iga::CouplingPort>& Ports() const noexcept override { return ports_; }
	void BeginStep(const iga::DomainStepContext& step) override
	{
		step.Validate();
		if (phase_ != Phase::Ready) throw std::runtime_error("fake begin phase");
		end_time_ = step.EndTime(); phase_ = Phase::Open; ++begins;
	}
	void SetPortInput(const std::string&, const iga::PortBoundaryData& input) override
	{
		if (phase_ != Phase::Open) throw std::runtime_error("fake input phase");
		iga::ValidatePortBoundaryData(input);
	}
	void SolveTrial() override { SolveHydraulicTrial(); SolveTransportTrial(); }
	iga::PortState GetPortState(const std::string& port) const override { return GetTransportPortState(port); }
	void RollbackTrial() override { RollbackTransportTrial(); RollbackHydraulicTrial(); }
	void AbortStep() override
	{
		if (phase_ == Phase::Ready) return;
		phase_ = Phase::Ready; input_.clear(); ++aborts;
		if (fail_abort) throw std::runtime_error("injected abort failure");
	}
	void PrepareCommitStep() override
	{
		if (phase_ != Phase::Transport) throw std::runtime_error("fake prepare phase");
		if (fail_prepare) throw std::runtime_error("injected prepare failure");
		phase_ = Phase::Prepared;
	}
	void FinalizeCommitStep() noexcept override { if (phase_ != Phase::Prepared) std::terminate(); phase_ = Phase::Ready; ++commits; }

	void SolveHydraulicTrial() override
	{
		if (phase_ != Phase::Open) throw std::runtime_error("fake hydraulic phase");
		phase_ = Phase::Hydraulic; ++hydraulic_solves;
	}
	iga::PortState GetHydraulicPortState(const std::string& port) const override
	{
		if (phase_ != Phase::Hydraulic && phase_ != Phase::Transport)
			throw std::runtime_error("fake hydraulic state phase");
		iga::PortState state; state.time_s = end_time_;
		state.outward_flow_m3_s = Flow(port); state.mean_pressure_pa = 1.0;
		return state;
	}
	void RollbackHydraulicTrial() override
	{
		if (phase_ != Phase::Hydraulic) throw std::runtime_error("fake hydraulic rollback phase");
		phase_ = Phase::Open; ++hydraulic_rollbacks;
	}
	void SetTransportConcentration(const std::string&, double time,
		const std::map<std::string, double>& concentration) override
	{
		if (phase_ != Phase::Hydraulic || time != end_time_ || concentration.size() != 2
			|| !concentration.count("tracer") || !concentration.count("dye"))
			throw std::runtime_error("fake concentration input");
		input_ = concentration; last_received = concentration.at("tracer");
		last_received_dye = concentration.at("dye"); ++concentration_inputs;
	}
	void SolveTransportTrial() override
	{
		if (phase_ != Phase::Hydraulic) throw std::runtime_error("fake transport phase");
		if (fail_transport) throw std::runtime_error("injected transport failure");
		if (!input_.empty()) concentration_ = input_;
		phase_ = Phase::Transport; ++transport_solves;
	}
	iga::PortState GetTransportPortState(const std::string& port) const override
	{
		if (phase_ != Phase::Transport) throw std::runtime_error("fake transport state phase");
		iga::PortState state; state.time_s = end_time_;
		state.outward_flow_m3_s = Flow(port); state.mean_pressure_pa = 1.0;
		for (const auto& value : concentration_) {
			state.concentration.emplace(value);
			state.outward_species_flux.emplace(value.first, Flow(port)*value.second);
		}
		return state;
	}
	void RollbackTransportTrial() override { if (phase_ != Phase::Transport) throw std::runtime_error("fake transport rollback"); phase_ = Phase::Hydraulic; input_.clear(); }
	std::map<std::string, iga::SpeciesStepAccounting> GetSpeciesStepAccounting() const override
	{
		if (phase_ != Phase::Transport) throw std::runtime_error("fake accounting phase");
		std::map<std::string, iga::SpeciesStepAccounting> result;
		for (const auto& species : {"dye", "tracer"}) {
			iga::SpeciesStepAccounting value;
			const double amount = Flow()*0.1;
			value.initial_mass = 10.0;
			value.final_mass = 10.0-amount+mass_roundoff;
			value.outward_port_amount.emplace("p", corrupt_port_amounts ? Flow() : amount);
			value.residual = corrupt_port_amounts ? 0.0 : mass_roundoff;
			result.emplace(species, std::move(value));
		}
		return result;
	}

	double signed_flow = 1.0;
	bool fail_transport = false;
	bool fail_prepare = false;
	bool fail_abort = false;
	bool corrupt_port_amounts = false;
	double mass_roundoff = 0.0;
	int begins = 0, aborts = 0, commits = 0, hydraulic_solves = 0;
	int hydraulic_rollbacks = 0, transport_solves = 0, concentration_inputs = 0;
	double last_received = 0.0, last_received_dye = 0.0;
	std::map<std::string, double> port_flows;

private:
	enum class Phase { Ready, Open, Hydraulic, Transport, Prepared };
	double Flow(const std::string& port = "p") const
	{
		const auto found = port_flows.find(port);
		if (found != port_flows.end()) return found->second;
		return id_ == "a" ? signed_flow : -signed_flow;
	}
	std::string id_; iga::DomainKind kind_; std::vector<iga::CouplingPort> ports_;
	std::map<std::string, double> concentration_; double end_time_ = 0.0;
	std::map<std::string, double> input_;
	Phase phase_ = Phase::Ready;
};

struct Fixture {
	Fixture()
	{
		graph = std::make_unique<iga::SimulationGraph>(Graph());
		auto first = std::make_unique<FakeStagedRuntime>("a", iga::DomainKind::OneDFlow,
			graph->Domain("a").ports, 2.0);
		a = first.get();
		auto second = std::make_unique<FakeStagedRuntime>("b", iga::DomainKind::ThreeDBodyFittedFlow,
			graph->Domain("b").ports, 5.0);
		b = second.get();
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		runtimes.push_back(std::move(second)); runtimes.push_back(std::move(first));
		registry = std::make_unique<iga::DomainRuntimeRegistry>(*graph, std::move(runtimes));
	}
	FakeStagedRuntime* a = nullptr; FakeStagedRuntime* b = nullptr;
	std::unique_ptr<iga::SimulationGraph> graph;
	std::unique_ptr<iga::DomainRuntimeRegistry> registry;
};

struct CycleFixture {
	CycleFixture()
	{
		graph = std::make_unique<iga::SimulationGraph>(CycleGraph());
		std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> runtimes;
		for (const auto& id : {"a", "b", "c", "d"}) {
			auto runtime = std::make_unique<FakeStagedRuntime>(id, graph->Domain(id).kind,
				graph->Domain(id).ports, 1.0);
			by_id.emplace(id, runtime.get());
			runtimes.push_back(std::move(runtime));
		}
		by_id.at("a")->port_flows = {{"ab", 1.0}, {"ad", -1.0}};
		by_id.at("b")->port_flows = {{"ab", -1.0}, {"bc", 1.0}};
		by_id.at("c")->port_flows = {{"bc", -1.0}, {"cd", 1.0}};
		by_id.at("d")->port_flows = {{"cd", -1.0}, {"ad", 1.0}};
		registry = std::make_unique<iga::DomainRuntimeRegistry>(*graph, std::move(runtimes));
	}
	std::unique_ptr<iga::SimulationGraph> graph;
	std::map<std::string, FakeStagedRuntime*> by_id;
	std::unique_ptr<iga::DomainRuntimeRegistry> registry;
};

iga::SpeciesPressureFlowExecutionControls Controls()
{
	iga::SpeciesPressureFlowExecutionControls controls;
	controls.amount_tolerances.emplace("tracer", iga::SpeciesAmountTolerance{});
	controls.amount_tolerances.emplace("dye", iga::SpeciesAmountTolerance{});
	return controls;
}

} // namespace

int main()
{
	const iga::DomainStepContext step{0, 0.0, 0.1};
	{
		Fixture original;
		iga::SpeciesPressureFlowComponentExecutor executor(*original.registry, "a", Controls());
		RequireRejected([&] { executor.CaptureCheckpointDonors(); });
		const std::map<std::pair<std::string, std::string>, iga::SpeciesDonor> first{
			{{"edge", "dye"}, iga::SpeciesDonor::First}, {{"edge", "tracer"}, iga::SpeciesDonor::First}};
		const auto accepted = executor.Advance(step, {{"edge", 0.0}}, [&](const auto&) {
			RequireRejected([&] { executor.CaptureCheckpointDonors(); });
			RequireRejected([&] { executor.RestoreCheckpointDonors(first); });
		});
		const auto donors = executor.CaptureCheckpointDonors(); assert(donors == first);
		const auto digest = std::string(64, 'a');
		const iga::CoupledCheckpointEpoch epoch{digest, {}, {digest, digest, digest, 1}, 1, 0.1, 0.1};
		const auto controls = iga::MakePressureFlowCheckpointControls(step, accepted.hydraulic_iterations.back(), donors);
		assert(controls.next_pressure_pa.at("edge") == 1.0);
		const auto encoded = iga::SerializePressureFlowCheckpointControls(controls);
		const auto decoded = iga::ParsePressureFlowCheckpointControls(encoded, *original.graph, epoch);
		assert(iga::SerializePressureFlowCheckpointControls(decoded) == encoded);
		assert(decoded.NextStep().step_index == 1 && decoded.NextStep().start_time_s == 0.1);
		for (std::size_t n = 0; n < encoded.size(); ++n)
			RequireRejected([&] { iga::ParsePressureFlowCheckpointControls(std::string_view(encoded).substr(0, n), *original.graph, epoch); });
		for (int mutation = 0; mutation < 5; ++mutation) {
			auto invalid = controls;
			if (mutation == 0) invalid.next_pressure_pa.clear();
			if (mutation == 1) invalid.next_pressure_pa.emplace("unknown", 0);
			if (mutation == 2) invalid.donors.clear();
			if (mutation == 3) invalid.accepted_step.step_index = 1;
			if (mutation == 4) invalid.accepted_step.dt_s = 0.2;
			const auto bad = iga::SerializePressureFlowCheckpointControls(invalid);
			RequireRejected([&] { iga::ParsePressureFlowCheckpointControls(bad, *original.graph, epoch); });
		}
		Fixture restarted;
		iga::SpeciesPressureFlowComponentExecutor resumed(*restarted.registry, "a", Controls());
		for (int mutation = 0; mutation < 3; ++mutation) {
			auto invalid = donors;
			if (mutation == 0) invalid.erase({"edge", "tracer"});
			if (mutation == 1) invalid[{"unknown", "tracer"}] = iga::SpeciesDonor::First;
			if (mutation == 2) invalid[{"edge", "tracer"}] = static_cast<iga::SpeciesDonor>(99);
			RequireRejected([&] { resumed.RestoreCheckpointDonors(invalid); });
			assert(resumed.CommittedDonorOwnership().empty());
			RequireRejected([&] { resumed.CaptureCheckpointDonors(); });
		}
		resumed.RestoreCheckpointDonors(decoded.donors);
		assert(resumed.CaptureCheckpointDonors() == donors);
		RequireRejected([&] { resumed.RestoreCheckpointDonors(donors); });
		restarted.a->signed_flow = 0.0; restarted.b->signed_flow = 0.0;
		const auto zero = resumed.Advance(decoded.NextStep(), decoded.next_pressure_pa);
		assert(zero.donor_ownership == donors && resumed.CaptureCheckpointDonors() == donors);
		restarted.a->signed_flow = -1.0; restarted.b->signed_flow = -1.0;
		const auto reverse = resumed.Advance({2, 0.2, 0.1}, {{"edge", 1.0}});
		assert(reverse.donor_ownership.at({"edge", "tracer"}) == iga::SpeciesDonor::Second);
		restarted.b->fail_prepare = true;
		RequireRejected([&] { resumed.Advance({3, 0.3, 0.1}, {{"edge", 1.0}}); });
		assert(resumed.CaptureCheckpointDonors() == reverse.donor_ownership);
	}
	{
		Fixture fixture;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		const auto result = executor.Advance(step, {{"edge", 0.0}});
		assert((result.transport_domain_order == std::vector<std::string>{"a", "b"}));
		assert(result.donor_ownership.at({"edge", "tracer"}) == iga::SpeciesDonor::First);
		assert(fixture.b->concentration_inputs == 1 && fixture.b->transport_solves == 1);
		assert(fixture.b->last_received == 2.0);
		assert(fixture.b->last_received_dye == 12.0);
		assert(fixture.a->transport_solves == 1 && fixture.a->commits == 1 && fixture.b->commits == 1);
		assert(result.edge_amounts.front().residual == 0.0);
		assert(result.global_balances.at("tracer").residual == 0.0);
		fixture.a->signed_flow = 0.0;
		fixture.b->signed_flow = 0.0;
		const auto zero = executor.Advance({1, 0.1, 0.1}, {{"edge", 0.0}});
		assert(zero.donor_ownership.at({"edge", "tracer"}) == iga::SpeciesDonor::First);
		const auto committed = executor.CommittedDonorOwnership();
		RequireRejected([&] { executor.Advance({2, 0.2, 0.1}, {{"edge", 0.0}},
			[](const iga::SpeciesPressureFlowStepResult&) { throw std::runtime_error("callback failure"); }); });
		assert(executor.CommittedDonorOwnership() == committed);
		assert(fixture.a->commits == 2 && fixture.b->commits == 2
			&& fixture.a->aborts == 1 && fixture.b->aborts == 1);
	}
	{
		Fixture fixture;
		auto controls = Controls();
		controls.hydraulic.method = iga::PressureFlowIterationMethod::Fixed;
		controls.hydraulic.maximum_iterations = 24;
		controls.hydraulic.pressure_relative_tolerance = 1.0e-6;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", controls);
		const auto result = executor.Advance(step, {{"edge", 0.0}});
		assert(result.hydraulic_iterations.size() > 1
			&& result.hydraulic_iterations.back().converged);
		assert(fixture.a->hydraulic_rollbacks+1 == fixture.a->hydraulic_solves
			&& fixture.a->transport_solves == 1 && fixture.b->transport_solves == 1);
	}
	{
		Fixture fixture;
		fixture.a->signed_flow = -1.0;
		fixture.b->signed_flow = -1.0;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		const auto result = executor.Advance(step, {{"edge", 0.0}});
		assert((result.transport_domain_order == std::vector<std::string>{"b", "a"}));
		assert(result.donor_ownership.at({"edge", "tracer"}) == iga::SpeciesDonor::Second);
		assert(fixture.a->concentration_inputs == 1);
		assert(fixture.a->last_received == 5.0);
	}
	{
		Fixture fixture;
		fixture.a->signed_flow = -1.0;
		fixture.b->signed_flow = -1.0;
		std::vector<std::pair<std::string, std::string>> observed_dependencies;
		iga::SpeciesPressureFlowExecutionSynchronization synchronization;
		synchronization.execution.hydraulic_batches = {{"a"}, {"b"}};
		synchronization.execution.solve_hydraulic_batch = [&](const auto& batch) {
			for (const auto& domain : batch)
				dynamic_cast<iga::StagedFlowTransportDomainRuntime&>(
					fixture.registry->Runtime(domain)).SolveHydraulicTrial();
		};
		synchronization.execution.make_batches = [&](const auto& order, const auto& dependencies) {
			observed_dependencies = dependencies;
			std::vector<std::vector<std::string>> batches;
			for (const auto& domain : order) batches.push_back({domain});
			return batches;
		};
		synchronization.execution.solve_transport_batch = [&](const auto& batch) {
			for (const auto& domain : batch)
				dynamic_cast<iga::StagedFlowTransportDomainRuntime&>(
					fixture.registry->Runtime(domain)).SolveTransportTrial();
		};
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a",
			Controls(), synchronization);
		const auto result = executor.Advance(step, {{"edge", 0.0}});
		assert((result.transport_domain_order == std::vector<std::string>{"b", "a"}));
		assert((observed_dependencies
			== std::vector<std::pair<std::string, std::string>>{{"b", "a"}, {"b", "a"}}));
		assert(fixture.a->transport_solves == 1 && fixture.b->transport_solves == 1);
	}
	{
		Fixture fixture;
		fixture.a->corrupt_port_amounts = true;
		fixture.b->corrupt_port_amounts = true;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		RequireRejected([&] { executor.Advance(step, {{"edge", 0.0}}); });
		assert(fixture.a->commits == 0 && fixture.b->commits == 0
			&& fixture.a->aborts == 1 && fixture.b->aborts == 1);
	}
	{
		Fixture fixture;
		fixture.a->mass_roundoff = 1.0e-12;
		fixture.b->mass_roundoff = 1.0e-12;
		auto controls = Controls();
		controls.amount_tolerances.at("tracer").absolute_tolerance = 1.0e-10;
		controls.amount_tolerances.at("dye").absolute_tolerance = 1.0e-10;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", controls);
		const auto result = executor.Advance(step, {{"edge", 0.0}});
		const auto& global = result.global_balances.at("tracer");
		assert(std::abs(global.residual-2.0e-12) < 1.0e-14 && global.gross_activity > 0.39
			&& global.normalized_residual > 0.0);
	}
	{
		Fixture fixture;
		fixture.b->fail_prepare = true;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		RequireRejected([&] { executor.Advance(step, {{"edge", 0.0}}); });
		assert(fixture.a->commits == 0 && fixture.b->commits == 0
			&& fixture.a->aborts == 1 && fixture.b->aborts == 1);
	}
	{
		Fixture fixture;
		fixture.b->fail_transport = true;
		fixture.a->fail_abort = true;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		bool rejected = false;
		try { executor.Advance(step, {{"edge", 0.0}}); }
		catch (const std::runtime_error& error) {
			assert(std::string(error.what()).find("injected transport failure") != std::string::npos);
			assert(std::string(error.what()).find("injected abort failure") != std::string::npos);
			rejected = true;
		}
		assert(rejected);
		assert(fixture.a->commits == 0 && fixture.b->commits == 0);
	}
	{
		CycleFixture fixture;
		iga::SpeciesPressureFlowComponentExecutor executor(*fixture.registry, "a", Controls());
		RequireRejected([&] { executor.Advance(step,
			{{"ab", 0.0}, {"ad", 0.0}, {"bc", 0.0}, {"cd", 0.0}}); });
		for (const auto& runtime : fixture.by_id) {
			assert(runtime.second->commits == 0 && runtime.second->aborts == 1
				&& runtime.second->transport_solves == 0);
		}
	}
	std::cout << "species pressure-flow executor tests passed\n";
	return 0;
}

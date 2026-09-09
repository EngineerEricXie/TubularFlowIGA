#ifndef IGA_ONE_D_FLOW_DOMAIN_ADAPTER_HPP
#define IGA_ONE_D_FLOW_DOMAIN_ADAPTER_HPP

#include "CoupledDomainRuntime.hpp"
#include "FlowDomainPortMetadata.hpp"
#include "OneDRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

class OneDFlowDomainAdapter : public CoupledDomainRuntime {
public:
	OneDFlowDomainAdapter(std::string domain_id, OneDFlowRuntime& runtime,
		std::vector<CouplingPort> ports, OneDInletPolicy inlet_policy)
		: domain_id_(std::move(domain_id)), runtime_(runtime), ports_(std::move(ports)),
		  inlet_policy_(inlet_policy)
	{
		ValidateOneDFlowDomainMetadata(domain_id_, ports_, inlet_policy_);
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::OneDFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		runtime_.RunLocalAdapterStage("1d flow adapter begin step", [&] {
			step.Validate();
			runtime_.BeginStep(step.start_time_s, step.dt_s);
			if (inlet_policy_ == OneDInletPolicy::ConfiguredOpenLoop)
				runtime_.SetConfiguredOpenLoopInlet();
		});
	}

	void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) override
	{
		runtime_.RunLocalAdapterStage("1d flow adapter input", [&] {
			runtime_.SetPortInput(Port(port_id).locator, input);
		});
	}

	void SolveTrial() override { runtime_.SolveTrial(); }

	PortState GetPortState(const std::string& port_id) const override
	{
		return runtime_.GetPortState(Port(port_id).locator);
	}

	void RollbackTrial() override
	{
		runtime_.RunLocalAdapterStage("1d flow adapter rollback", [&] {
			runtime_.RollbackTrial();
		});
	}
	void AbortStep() override
	{
		runtime_.RunLocalAdapterStage("1d flow adapter abort", [&] {
			runtime_.AbortStep();
		});
	}
	void PrepareCommitStep() override
	{
		runtime_.RunLocalAdapterStage("1d flow adapter prepare commit", [&] {
			runtime_.PrepareCommitStep();
		});
	}
	void FinalizeCommitStep() noexcept override { runtime_.FinalizeCommitStep(); }

private:
	const CouplingPort& Port(const std::string& port_id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(),
			[&](const CouplingPort& port) { return port.id == port_id; });
		if (found == ports_.end())
			throw std::runtime_error("1D domain adapter has no port '"+port_id+"'");
		return *found;
	}

	std::string domain_id_;
	OneDFlowRuntime& runtime_;
	std::vector<CouplingPort> ports_;
	OneDInletPolicy inlet_policy_;
};

// The species-aware adapter owns the staged coupling contract.  The original
// OneDFlowDomainAdapter remains flow-only so established callers continue to
// use OneDFlowRuntime::SolveTrial without an altered numerical route.
struct OneDFlowTransportDomainControls {
	std::optional<double> species_flow_epsilon_m3_s;

	void Validate() const
	{
		if (species_flow_epsilon_m3_s
			&& (*species_flow_epsilon_m3_s < 0.0
				|| !std::isfinite(*species_flow_epsilon_m3_s)))
			throw std::runtime_error(
				"1D staged species flow epsilon must be finite and nonnegative");
	}
};

class OneDFlowTransportDomainAdapter : public CoupledDomainRuntime,
	public StagedFlowTransportDomainRuntime {
public:
	OneDFlowTransportDomainAdapter(std::string domain_id, OneDFlowRuntime& runtime,
		std::vector<CouplingPort> ports, OneDInletPolicy inlet_policy,
		std::map<std::string, std::string> species_bindings,
		OneDFlowTransportDomainControls controls = {})
		: domain_id_(std::move(domain_id)), runtime_(runtime), ports_(std::move(ports)),
		  inlet_policy_(inlet_policy), species_bindings_(std::move(species_bindings)),
		  controls_(std::move(controls))
	{
		ValidateOneDFlowDomainMetadata(domain_id_, ports_, inlet_policy_);
		controls_.Validate();
		std::set<std::string> native_species;
		for (const auto& binding : species_bindings_) {
			if (binding.first.empty() || binding.second.empty()
				|| !FindOneDSpecies(runtime_.Transports(), binding.second))
				throw std::runtime_error("1D staged species binding must name a transported native field");
			if (!native_species.insert(binding.second).second)
				throw std::runtime_error("1D staged species bindings require unique native fields");
		}
		for (const auto& port : ports_)
			for (const auto& species : port.species)
				if (!species_bindings_.count(species))
					throw std::runtime_error("1D staged coupled port species '"+species
						+"' has no native binding");
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::OneDFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter begin step", [&] {
			step.Validate();
			if (runtime_.Configuration().physiology.vasodilation)
				throw std::runtime_error(
					"staged 1D flow/transport rejects enabled concentration-driven vasodilation");
			runtime_.BeginStep(step.start_time_s, step.dt_s);
			if (inlet_policy_ == OneDInletPolicy::ConfiguredOpenLoop)
				runtime_.SetConfiguredOpenLoopInlet();
			step_ = step;
			hydraulic_inputs_.clear();
			concentration_inputs_.clear();
			hydraulic_trial_succeeded_ = false;
			transport_trial_succeeded_ = false;
		});
	}

	void SetPortInput(const std::string& port_id, const PortBoundaryData& input) override
	{
		decltype(hydraulic_inputs_) hydraulics;
		decltype(concentration_inputs_) concentrations;
		runtime_.RunLocalAdapterStage("1d staged adapter input", [&] {
			hydraulics = hydraulic_inputs_;
			concentrations = concentration_inputs_;

			RequireOpen("1D flow/transport input");
			const auto& port = Port(port_id);
			ValidatePortBoundaryData(input);
			ValidateTime(input.time_s);
			if (!input.outward_species_flux.empty())
				throw std::runtime_error("1D peer species flux is executor-owned, not a boundary input");
			const auto hydraulic = HydraulicInput(input);
			if (HydraulicValueCount(hydraulic) > 0) {
				if (hydraulics.count(port_id))
					throw std::runtime_error("1D staged hydraulic input was supplied more than once");
				hydraulics.emplace(port_id, hydraulic);
			}
			if (!input.concentration.empty()) {
				if (!port.requires.count(PortQuantity::SpeciesConcentration))
					throw std::runtime_error("1D staged port does not accept concentration input");
				if (concentrations.count(port_id))
					throw std::runtime_error("1D staged concentration was supplied more than once");
				concentrations.emplace(port_id, input.concentration);
			}
		});
		hydraulic_inputs_.swap(hydraulics);
		concentration_inputs_.swap(concentrations);
	}

	void SolveTrial() override
	{
		SolveHydraulicTrial();
		SolveTransportTrial();
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		return GetTransportPortState(port_id);
	}

	void RollbackTrial() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter rollback preparation", [&] {
			if (!hydraulic_trial_succeeded_)
				throw std::runtime_error("1D staged rollback requires a hydraulic trial");
		});
		if (transport_trial_succeeded_) RollbackTransportTrial();
		RollbackHydraulicTrial();
	}

	void AbortStep() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter abort", [&] {
			runtime_.AbortStep();
			hydraulic_inputs_.clear();
			concentration_inputs_.clear();
			hydraulic_trial_succeeded_ = false;
			transport_trial_succeeded_ = false;
		});
	}

	void PrepareCommitStep() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter prepare commit", [&] {
			if (!hydraulic_trial_succeeded_ || !transport_trial_succeeded_)
				throw std::runtime_error("1D staged prepare requires hydraulic and transport trials");
			runtime_.PrepareCommitStep();
		});
	}

	void FinalizeCommitStep() noexcept override
	{
		runtime_.FinalizeCommitStep();
		hydraulic_inputs_.clear();
		concentration_inputs_.clear();
		hydraulic_trial_succeeded_ = false;
		transport_trial_succeeded_ = false;
	}

	void SolveHydraulicTrial() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter hydraulic preparation", [&] {
			RequireOpen("1D staged hydraulic solve");
			for (const auto& port : ports_) {
				const auto found = hydraulic_inputs_.find(port.id);
				const int supplied = found == hydraulic_inputs_.end() ? 0 : HydraulicValueCount(found->second);
				if (HydraulicRequirementCount(port) != supplied)
					throw std::runtime_error("1D staged hydraulic input is missing or extra for port '"+port.id+"'");
				if (supplied) runtime_.SetPortInput(port.locator, found->second);
			}
		});
		hydraulic_trial_succeeded_ = false;
		runtime_.SolveHydraulicTrial();
		hydraulic_trial_succeeded_ = true;
	}

	PortState GetHydraulicPortState(const std::string& port_id) const override
	{
		if (!hydraulic_trial_succeeded_
			|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::HydraulicSolved)
			throw std::runtime_error("1D staged hydraulic state requires a successful hydraulic trial");
		auto state = runtime_.GetPortState(Port(port_id).locator);
		state.concentration.clear();
		state.outward_species_flux.clear();
		return state;
	}

	void RollbackHydraulicTrial() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter hydraulic rollback preparation", [&] {
			if (!hydraulic_trial_succeeded_
				|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::HydraulicSolved)
				throw std::runtime_error("1D staged hydraulic rollback requires a solved hydraulic trial");
		});
		runtime_.RollbackHydraulicTrial();
		hydraulic_inputs_.clear();
		hydraulic_trial_succeeded_ = false;
	}

	void SetTransportConcentration(const std::string& port_id, double time_s,
		const std::map<std::string, double>& concentration) override
	{
		decltype(concentration_inputs_) candidate;
		runtime_.RunLocalAdapterStage("1d staged adapter concentration", [&] {
			candidate = concentration_inputs_;

			if (!hydraulic_trial_succeeded_
				|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::HydraulicSolved)
				throw std::runtime_error("1D staged concentration requires accepted hydraulics");
			const auto& port = Port(port_id);
			if (!port.requires.count(PortQuantity::SpeciesConcentration)
				|| port.species.empty() || SpeciesKeys(concentration) != port.species)
				throw std::runtime_error("1D staged concentration requires the complete logical species map");
			for (const auto& value : concentration)
				if (!std::isfinite(value.second))
					throw std::runtime_error("1D staged concentration must be finite");
			ValidateTime(time_s);
			if (candidate.count(port_id))
				throw std::runtime_error("1D staged concentration was supplied more than once");
			candidate.emplace(port_id, concentration);
		});
		concentration_inputs_.swap(candidate);
	}

	void SolveTransportTrial() override
	{
		// An illegal phase is not a new scalar attempt: preserve the previous
		// input/acceptance image, including when another group member rejects it.
		runtime_.RunLocalAdapterStage("1d staged adapter transport state", [&] {
			if (!hydraulic_trial_succeeded_
				|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::HydraulicSolved)
				throw std::runtime_error("1D staged transport requires accepted hydraulics");
		});
		try {
			std::map<std::string, double> root;
			std::map<int, std::map<std::string, double>> outlets;
			OneDStagedRootTransportOwnership root_ownership
				= OneDStagedRootTransportOwnership::Legacy;
			runtime_.RunLocalAdapterStage("1d staged adapter transport preparation", [&] {
				for (const auto& port : ports_) {
					if (port.species.empty()
						|| !port.requires.count(PortQuantity::SpeciesConcentration)) continue;
					const auto found = concentration_inputs_.find(port.id);
					const bool supplied = found != concentration_inputs_.end();
					const auto flow = GetHydraulicPortState(port.id).outward_flow_m3_s;
					if (!flow) throw std::runtime_error("1D staged species port omitted outward flow");
					const double epsilon = SpeciesFlowEpsilon();
					const bool configured_root = inlet_policy_ == OneDInletPolicy::ConfiguredOpenLoop
						&& port.locator == "root";
					if (configured_root && supplied)
						throw std::runtime_error(
							"1D configured-open-loop root owns its scheduled concentration");
					const bool receiver = configured_root || *flow < -epsilon
						|| (std::abs(*flow) <= epsilon && supplied);
					if (!configured_root && *flow < -epsilon && !supplied)
						throw std::runtime_error("1D inward species port requires a complete concentration");
					if (*flow > epsilon && supplied)
						throw std::runtime_error("1D outward species port must not impose concentration");
					ValidateFrameOwnership(port, receiver, epsilon);
					if (port.locator == "root") {
						root_ownership = receiver
							? OneDStagedRootTransportOwnership::BoundaryConcentration
							: OneDStagedRootTransportOwnership::InteriorDonor;
					}
					if (!supplied) continue;
					std::map<std::string, double> native;
					for (const auto& logical : port.species)
						native.emplace(species_bindings_.at(logical), found->second.at(logical));
					if (port.locator == "root") root = std::move(native);
					else outlets.emplace(NativeOutletNode(port), std::move(native));
				}
				transport_trial_succeeded_ = false;
			});
			runtime_.SolveStagedTransportTrial(root, outlets, root_ownership);
			transport_trial_succeeded_ = true;
		} catch (...) {
			// Boundary data belongs to one scalar attempt.  Whether rejection
			// happened before native replay or inside it, retries must explicitly
			// provide a fresh complete concentration map against fixed frames.
			concentration_inputs_.clear();
			transport_trial_succeeded_ = false;
			throw;
		}
	}

	PortState GetTransportPortState(const std::string& port_id) const override
	{
		if (!hydraulic_trial_succeeded_ || !transport_trial_succeeded_
			|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::TrialSolved)
			throw std::runtime_error("1D staged transport state requires successful trials");
		auto state = runtime_.GetPortState(Port(port_id).locator);
		LogicalizeState(Port(port_id), state);
		return state;
	}

	void RollbackTransportTrial() override
	{
		runtime_.RunLocalAdapterStage("1d staged adapter transport rollback preparation", [&] {
			if (!hydraulic_trial_succeeded_ || !transport_trial_succeeded_
				|| runtime_.CurrentPhase() != OneDFlowRuntime::Phase::TrialSolved)
				throw std::runtime_error("1D staged transport rollback requires successful trials");
		});
		runtime_.RollbackStagedTransportTrial();
		concentration_inputs_.clear();
		transport_trial_succeeded_ = false;
	}

	std::map<std::string, SpeciesStepAccounting> GetSpeciesStepAccounting() const override
	{
		if (!hydraulic_trial_succeeded_ || !transport_trial_succeeded_)
			throw std::runtime_error("1D staged accounting requires successful trials");
		const auto native = runtime_.GetSpeciesStepAccounting();
		std::map<std::string, SpeciesStepAccounting> result;
		for (const auto& binding : species_bindings_) {
			const auto& value = native.at(binding.second);
			SpeciesStepAccounting accounting;
			accounting.initial_mass = value.initial_mass;
			accounting.final_mass = value.final_mass;
			accounting.source_amount = value.source_amount;
			for (const auto& port : ports_) if (port.species.count(binding.first)) {
				if (port.locator == "root")
					accounting.outward_port_amount.emplace(port.id, value.root_outward_amount);
				else accounting.outward_port_amount.emplace(port.id,
					value.outlet_outward_amount.at(NativeOutletNode(port)));
			}
			accounting.residual = value.balance_residual;
			result.emplace(binding.first, std::move(accounting));
		}
		return result;
	}

private:
	const CouplingPort& Port(const std::string& port_id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(),
			[&](const CouplingPort& port) { return port.id == port_id; });
		if (found == ports_.end()) throw std::runtime_error("1D staged adapter has no port '"+port_id+"'");
		return *found;
	}

	int NativeOutletNode(const CouplingPort& port) const
	{
		const int node_id = ParseOneDOutletNodeLocator(port);
		const auto found = runtime_.Network().node_index.find(node_id);
		if (found == runtime_.Network().node_index.end())
			throw std::runtime_error("1D staged outlet locator is not in the native network");
		return found->second;
	}

	double FrameOutwardFlow(const CouplingPort& port,
		const OneDFlowRuntime::HydraulicFrame& frame) const
	{
		if (port.locator == "root") return -frame.post_flow.inlet_flow;
		const int node = NativeOutletNode(port);
		const auto found = std::find_if(frame.post_flow.outlets.begin(),
			frame.post_flow.outlets.end(),
			[node](const OneDOutletState& outlet) { return outlet.node == node; });
		if (found == frame.post_flow.outlets.end())
			throw std::runtime_error("1D staged frame is missing a configured outlet state");
		return found->flow;
	}

	void ValidateFrameOwnership(const CouplingPort& port, bool receiver,
		double epsilon) const
	{
		const auto& frames = runtime_.HydraulicFrames();
		for (std::size_t index = 0; index < frames.size(); ++index) {
			const auto& frame = frames[index];
			if (port.locator == "root") {
				bool materially_inward = false;
				bool materially_outward = false;
				for (const double native_flow : runtime_.HydraulicFrameRootBranchNativeFlows(index)) {
					materially_inward = materially_inward || native_flow > epsilon;
					materially_outward = materially_outward || native_flow < -epsilon;
				}
				if (materially_inward && materially_outward)
					throw std::runtime_error(
						"1D staged root has mixed material branch directions within a hydraulic frame");
				if (receiver && materially_outward)
					throw std::runtime_error(
						"1D staged root branch changes from receiver to donor within hydraulic frames");
				if (!receiver && materially_inward)
					throw std::runtime_error(
						"1D staged root branch changes from donor to receiver within hydraulic frames");
				continue;
			}
			const double flow = FrameOutwardFlow(port, frame);
			if (receiver && flow > epsilon)
				throw std::runtime_error("1D staged port changes from receiver to donor within hydraulic frames");
			if (!receiver && flow < -epsilon)
				throw std::runtime_error("1D staged port changes from donor to receiver within hydraulic frames");
		}
	}

	double SpeciesFlowEpsilon() const
	{
		return controls_.species_flow_epsilon_m3_s.value_or(
			runtime_.Configuration().coupling.flow_epsilon_m3_s);
	}

	void RequireOpen(const char* operation) const
	{
		if (runtime_.CurrentPhase() != OneDFlowRuntime::Phase::TrialOpen)
			throw std::runtime_error(std::string(operation)+" requires a trial-open step");
	}

	void ValidateTime(double time_s) const
	{
		const double scale = std::max({1.0, std::abs(step_.EndTime()), std::abs(time_s)});
		if (!std::isfinite(time_s) || std::abs(time_s-step_.EndTime()) > 1.0e-12*scale)
			throw std::runtime_error("1D staged input time does not match the active step");
	}

	static int HydraulicRequirementCount(const CouplingPort& port)
	{
		return port.requires.count(PortQuantity::FlowRate)
			+port.requires.count(PortQuantity::MeanPressure);
	}

	static int HydraulicValueCount(const PortBoundaryData& input)
	{
		return static_cast<int>(input.outward_flow_m3_s.has_value())
			+static_cast<int>(input.mean_pressure_pa.has_value());
	}

	static PortBoundaryData HydraulicInput(const PortBoundaryData& input)
	{
		PortBoundaryData result;
		result.time_s = input.time_s;
		result.outward_flow_m3_s = input.outward_flow_m3_s;
		result.mean_pressure_pa = input.mean_pressure_pa;
		return result;
	}

	static std::set<std::string> SpeciesKeys(const std::map<std::string, double>& values)
	{
		std::set<std::string> result;
		for (const auto& value : values) result.insert(value.first);
		return result;
	}

	void LogicalizeState(const CouplingPort& port, PortState& state) const
	{
		std::map<std::string, double> concentration;
		std::map<std::string, double> flux;
		for (const auto& logical : port.species) {
			const auto& native = species_bindings_.at(logical);
			concentration.emplace(logical, state.concentration.at(native));
			flux.emplace(logical, state.outward_species_flux.at(native));
		}
		state.concentration = std::move(concentration);
		state.outward_species_flux = std::move(flux);
		ValidatePortState(state);
	}

	std::string domain_id_;
	OneDFlowRuntime& runtime_;
	std::vector<CouplingPort> ports_;
	OneDInletPolicy inlet_policy_;
	std::map<std::string, std::string> species_bindings_;
	OneDFlowTransportDomainControls controls_;
	DomainStepContext step_;
	std::map<std::string, PortBoundaryData> hydraulic_inputs_;
	std::map<std::string, std::map<std::string, double>> concentration_inputs_;
	bool hydraulic_trial_succeeded_ = false;
	bool transport_trial_succeeded_ = false;
};

} // namespace iga

#endif

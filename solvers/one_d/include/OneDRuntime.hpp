#ifndef IGA_ONE_D_RUNTIME_HPP
#define IGA_ONE_D_RUNTIME_HPP

#include "OneDCoupling.hpp"
#include "CouplingPort.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace iga {

// Owns all mutable native 1d state.  PETSc is deliberately injected for the
// implicit scheme so the common lifecycle and its fast tests remain C++17-only.
class OneDFlowRuntime {
public:
	using ImplicitAdvance = std::function<void(const OneDNetwork&,
		const OneDFlowSystemDefinition&, OneDFlowState&, double, double)>;

	enum class Phase { Uninitialized, Ready, TrialOpen, TrialSolved };

	OneDFlowRuntime(OneDConfiguration configuration, OneDFlowSystemDefinition flow,
		OneDNetwork network, OneDInletState inlet, std::filesystem::path case_directory,
		ImplicitAdvance implicit_advance = {})
		: configuration_(std::move(configuration)), flow_(std::move(flow)),
		  network_(std::move(network)), inlet_(std::move(inlet)),
		  case_directory_(std::move(case_directory)), implicit_advance_(std::move(implicit_advance))
	{
		flow_state_.outlets = ResolveOneDOutlets(configuration_, network_);
		for (const auto& transport : configuration_.transport_systems)
			if (transport.flow_system == flow_.name)
				transports_.push_back(InitializeOneDTransport(configuration_, transport, network_));
	}

	void InitializeOpenLoop(double inlet_flow)
	{
		VascularInletState inlet;
		inlet.has_flow = true;
		inlet.flow_m3_s = inlet_flow;
		InitializeOpenLoop(inlet);
	}

	void InitializeOpenLoop(const VascularInletState& inlet)
	{
		RequirePhase(Phase::Uninitialized, "InitializeOpenLoop");
		ValidateVascularInletState(inlet);
		if (!inlet.has_flow || !Close(inlet.time_s, 0.0))
			throw std::runtime_error("1d initial open-loop inlet requires flow at time zero");
		InitializeFlow(inlet.flow_m3_s);
		last_inlet_ = inlet;
		phase_ = Phase::Ready;
	}

	void InitializeCoupled(const VascularInletState& inlet)
	{
		RequirePhase(Phase::Uninitialized, "InitializeCoupled");
		const double flow = ApplyOneDCoupledInlet(configuration_, transports_, inlet);
		InitializeFlow(flow);
		last_inlet_ = inlet;
		phase_ = Phase::Ready;
	}

	VascularInletState OpenLoopInlet(double time_s, double flow_m3_s) const
	{
		VascularInletState result;
		result.time_s = time_s;
		result.has_flow = true;
		result.flow_m3_s = flow_m3_s;
		for (const auto& transport : transports_)
			for (const auto& species : transport.species)
				result.species.emplace(species.definition.field,
					EvaluateOneDSpeciesInlet(configuration_, species, case_directory_, time_s));
		if (configuration_.physiology.enabled) {
			result.has_hematocrit = true;
			result.hematocrit_percent = configuration_.physiology.hematocrit_percent;
		}
		return result;
	}

	void BeginStep(double time_s, double dt_s)
	{
		RequirePhase(Phase::Ready, "BeginStep");
		if (!(dt_s > 0.0) || !std::isfinite(dt_s))
			throw std::runtime_error("1d BeginStep requires a finite positive dt_s");
		if (!Close(time_s, flow_state_.physical_time))
			throw std::runtime_error("1d BeginStep time_s must equal committed physical time");
		committed_ = Snapshot{configuration_, network_, flow_state_, transports_, last_inlet_};
		trial_time_s_ = time_s;
		trial_dt_s_ = dt_s;
		trial_inlet_.reset();
		trial_outlet_pressure_overrides_.clear();
		trial_coupled_ = false;
		trial_solve_succeeded_ = false;
		phase_ = Phase::TrialOpen;
	}

	// Generic common-port entry point.  Native root flow is root-to-leaf, hence
	// positive native flow is inward at the root and outward at distal leaves.
	void SetPortInput(const std::string& port_id, const PortBoundaryData& input)
	{
		RequirePhase(Phase::TrialOpen, "SetPortInput");
		ValidatePortBoundaryData(input);
		if (!Close(input.time_s, trial_time_s_+trial_dt_s_))
			throw std::runtime_error("1d port input time_s must equal trial end time");
		if (port_id == "root") {
			if (!input.outward_flow_m3_s)
				throw std::runtime_error("native 1d root trial input requires outward_flow_m3_s");
			if (input.mean_pressure_pa || input.mean_normal_traction_pa || input.total_pressure_pa
				|| !input.outward_species_flux.empty())
				throw std::runtime_error(
					"native 1d root supports flow and concentration trial data only");
			VascularInletState inlet;
			inlet.time_s = input.time_s;
			inlet.has_flow = true;
			inlet.flow_m3_s = RootOrientation().ToNative(*input.outward_flow_m3_s);
			inlet.species = input.concentration;
			trial_inlet_ = std::move(inlet);
			trial_coupled_ = true;
			return;
		}
		const int node = OutletNode(port_id);
		if (!input.mean_pressure_pa || input.outward_flow_m3_s
			|| input.mean_normal_traction_pa || input.total_pressure_pa
			|| !input.concentration.empty() || !input.outward_species_flux.empty())
			throw std::runtime_error(
				"native 1d outlet trial input requires only mean_pressure_pa");
		const auto& outlet = flow_state_.outlets.at(OutletIndex(node));
		if (outlet.kind != OneDOutletKind::Pressure)
			throw std::runtime_error(
				"native 1d outlet pressure cannot replace a resistance or RCR closure");
		trial_outlet_pressure_overrides_[node] = *input.mean_pressure_pa;
	}

	// Compatibility path for VCA/replay inputs; coupled perfusate mutations are
	// applied only while solving the trial and are included in the rollback image.
	void SetCoupledInlet(const VascularInletState& inlet)
	{
		RequirePhase(Phase::TrialOpen, "SetCoupledInlet");
		ValidateVascularInletState(inlet);
		if (!Close(inlet.time_s, trial_time_s_+trial_dt_s_))
			throw std::runtime_error("1d coupled inlet time_s must equal trial end time");
		trial_inlet_ = inlet;
		trial_coupled_ = true;
	}

	void SetOpenLoopInlet(const VascularInletState& inlet)
	{
		RequirePhase(Phase::TrialOpen, "SetOpenLoopInlet");
		ValidateVascularInletState(inlet);
		if (!inlet.has_flow || !Close(inlet.time_s, trial_time_s_+trial_dt_s_))
			throw std::runtime_error("1d open-loop inlet must contain end-time flow data");
		trial_inlet_ = inlet;
		trial_coupled_ = false;
	}

	void SolveTrial()
	{
		RequirePhase(Phase::TrialOpen, "SolveTrial");
		if (!trial_inlet_ || !trial_inlet_->has_flow)
			throw std::runtime_error("1d SolveTrial requires a root flow input");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		try {
			for (const auto& override : trial_outlet_pressure_overrides_) {
				auto& outlet = flow_state_.outlets.at(OutletIndex(override.first));
				outlet.pressure = override.second;
				outlet.reference_pressure = override.second;
				outlet.capacitor_pressure = override.second;
			}
			const double inlet_flow = trial_coupled_
				? ApplyOneDCoupledInlet(configuration_, transports_, *trial_inlet_)
				: trial_inlet_->flow_m3_s;
			if (flow_.scheme == OneDFlowScheme::SteadyPoiseuille)
				SolveRigidOneD(network_, flow_, flow_state_, inlet_flow, trial_dt_s_);
			else if (flow_.scheme == OneDFlowScheme::ExplicitRusanov)
				AdvanceExplicitOneD(network_, flow_, flow_state_, inlet_flow, trial_dt_s_);
			else {
				if (!implicit_advance_)
					throw std::runtime_error(
						"1d implicit trial solve requires an injected PETSc advance function");
				implicit_advance_(network_, flow_, flow_state_, inlet_flow, trial_dt_s_);
			}
			for (auto& transport : transports_)
				AdvanceOneDTransport(configuration_, network_, flow_state_, transport,
					case_directory_, trial_time_s_, trial_dt_s_);
			ApplyOneDVasodilation(configuration_, network_, transports_, trial_dt_s_,
				flow_.dynamic_viscosity);
			flow_state_.completed_step = committed_.flow.completed_step+1;
			flow_state_.physical_time = trial_time_s_+trial_dt_s_;
			last_inlet_ = *trial_inlet_;
			trial_solve_succeeded_ = true;
			phase_ = Phase::TrialSolved;
		} catch (...) {
			phase_ = Phase::TrialSolved;
			throw;
		}
	}

	void RollbackTrial()
	{
		RequirePhase(Phase::TrialSolved, "RollbackTrial");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		phase_ = Phase::TrialOpen;
	}

	void CommitStep()
	{
		RequirePhase(Phase::TrialSolved, "CommitStep");
		if (!trial_solve_succeeded_)
			throw std::runtime_error("1d CommitStep requires a successful trial solve");
		committed_ = Snapshot{};
		trial_outlet_pressure_overrides_.clear();
		trial_solve_succeeded_ = false;
		phase_ = Phase::Ready;
	}

	PortState GetPortState(const std::string& port_id) const
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::Ready)
			throw std::runtime_error("1d GetPortState requires a solved or committed state");
		int node = network_.root;
		PortOrientation orientation = RootOrientation();
		bool root = port_id == "root";
		if (!root) {
			node = OutletNode(port_id);
			orientation.native_to_outward_sign = 1;
		}
		PortState state;
		state.time_s = flow_state_.physical_time;
		const int segment_index = root ? OneDSegmentsOutOfNode(network_, node).front()
			: OneDSegmentIntoNode(network_, node);
		const auto& segment = network_.segments.at(static_cast<std::size_t>(segment_index));
		const int cell_index = root ? segment.cell_offset : segment.cell_offset+segment.cells-1;
		const auto cell = static_cast<std::size_t>(cell_index);
		state.area_m2 = flow_state_.area.at(cell);
		const double native_flow = root ? flow_state_.inlet_flow
			: flow_state_.outlets.at(OutletIndex(node)).flow;
		state.outward_flow_m3_s = orientation.ToOutward(native_flow);
		state.mean_pressure_pa = flow_state_.node_pressure.at(static_cast<std::size_t>(node));
		for (const auto& transport : transports_)
			for (const auto& species : transport.species) {
				double concentration = species.concentration.at(cell);
				if (root) {
					const auto boundary = last_inlet_.species.find(species.definition.field);
					if (boundary != last_inlet_.species.end()) concentration = boundary->second;
					else concentration = species.inlet_value;
				}
				state.concentration.emplace(species.definition.field, concentration);
				state.outward_species_flux.emplace(species.definition.field,
					*state.outward_flow_m3_s*concentration);
			}
		ValidatePortState(state);
		return state;
	}

	void RestoreCommittedState(OneDFlowState flow, std::vector<OneDTransportState> transports,
		OneDNetwork network)
	{
		RequirePhase(Phase::Ready, "RestoreCommittedState");
		flow_state_ = std::move(flow);
		transports_ = std::move(transports);
		network_ = std::move(network);
	}

	const OneDConfiguration& Configuration() const { return configuration_; }
	const OneDFlowSystemDefinition& FlowSystem() const { return flow_; }
	const OneDNetwork& Network() const { return network_; }
	const OneDInletState& InletDefinition() const { return inlet_; }
	const OneDFlowState& FlowState() const { return flow_state_; }
	const std::vector<OneDTransportState>& Transports() const { return transports_; }
	const VascularInletState& LastInlet() const { return last_inlet_; }
	Phase CurrentPhase() const { return phase_; }

	static PortOrientation RootOrientation() { return {-1}; }

private:
	struct Snapshot {
		OneDConfiguration configuration;
		OneDNetwork network;
		OneDFlowState flow;
		std::vector<OneDTransportState> transports;
		VascularInletState last_inlet;
	};

	static bool Close(double first, double second)
	{
		return std::abs(first-second) <= 1.0e-12*std::max({1.0, std::abs(first), std::abs(second)});
	}

	void RequirePhase(Phase required, const char* operation) const
	{
		if (phase_ != required) throw std::runtime_error(std::string("illegal 1d runtime transition: ")+operation);
	}

	void InitializeFlow(double inlet_flow)
	{
		if (flow_.model == OneDFlowModel::Rigid)
			SolveRigidOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
		else InitializeCompliantOneDFromRigid(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
	}

	void RestoreCommitted()
	{
		configuration_ = committed_.configuration;
		network_ = committed_.network;
		flow_state_ = committed_.flow;
		transports_ = committed_.transports;
		last_inlet_ = committed_.last_inlet;
	}

	int OutletNode(const std::string& port_id) const
	{
		const std::string prefix = "outlet:";
		if (port_id.compare(0, prefix.size(), prefix) != 0)
			throw std::runtime_error("1d port id must be root or outlet:<node-id>");
		int node_id = 0;
		const auto first = port_id.data()+prefix.size();
		const auto parsed = std::from_chars(first, port_id.data()+port_id.size(), node_id);
		if (first == port_id.data()+port_id.size() || parsed.ec != std::errc{}
			|| parsed.ptr != port_id.data()+port_id.size())
			throw std::runtime_error("1d outlet port node id must be an integer");
		const auto found = network_.node_index.find(node_id);
		if (found == network_.node_index.end())
			throw std::runtime_error("unknown 1d outlet port node");
		if (std::find(network_.outlet_nodes.begin(), network_.outlet_nodes.end(), found->second)
			== network_.outlet_nodes.end())
			throw std::runtime_error("1d port is not an outlet leaf");
		return found->second;
	}

	std::size_t OutletIndex(int node) const
	{
		for (std::size_t index = 0; index < flow_state_.outlets.size(); ++index)
			if (flow_state_.outlets[index].node == node) return index;
		throw std::runtime_error("1d outlet state is missing");
	}

	OneDConfiguration configuration_;
	OneDFlowSystemDefinition flow_;
	OneDNetwork network_;
	OneDInletState inlet_;
	std::filesystem::path case_directory_;
	ImplicitAdvance implicit_advance_;
	OneDFlowState flow_state_;
	std::vector<OneDTransportState> transports_;
	VascularInletState last_inlet_;
	Snapshot committed_;
	std::optional<VascularInletState> trial_inlet_;
	std::map<int, double> trial_outlet_pressure_overrides_;
	bool trial_coupled_ = false;
	bool trial_solve_succeeded_ = false;
	double trial_time_s_ = 0.0;
	double trial_dt_s_ = 0.0;
	Phase phase_ = Phase::Uninitialized;
};

} // namespace iga

#endif

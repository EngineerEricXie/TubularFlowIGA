#ifndef IGA_ONE_D_RUNTIME_HPP
#define IGA_ONE_D_RUNTIME_HPP

#include "OneDCoupling.hpp"
#include "CouplingPort.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

// Owns all mutable native 1d state.  PETSc is deliberately injected for the
// implicit scheme so the common lifecycle and its fast tests remain C++17-only.
class OneDFlowRuntime {
public:
	using ImplicitAdvance = std::function<void(const OneDNetwork&,
		const OneDFlowSystemDefinition&, OneDFlowState&, double, double)>;

	enum class Phase { Uninitialized, Ready, TrialOpen, TrialSolved, CommitPrepared };
	enum class TrialInletMode { None, HeldOpenLoop, Coupled, ConfiguredOpenLoop };

	struct TrialDiagnostics {
		int planned_configured_substeps = 0;
		int attempted_configured_substeps = 0;
		int completed_configured_substeps = 0;
		long long explicit_cfl_substep_delta = 0;
		std::vector<double> configured_open_loop_endpoint_times_s;
		std::vector<double> configured_open_loop_endpoint_flows_m3_s;
		std::vector<double> substep_endpoint_flows_m3_s;
	};

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
		const int configured_substeps = ConfiguredSubsteps(dt_s);
		if (flow_state_.completed_step > configuration_.time.steps-configured_substeps)
			throw std::runtime_error("1d BeginStep would exceed configured step count");
		committed_ = Snapshot{configuration_, network_, flow_state_, transports_, last_inlet_};
		trial_time_s_ = time_s;
		trial_dt_s_ = dt_s;
		trial_inlet_.reset();
		trial_outlet_pressure_overrides_.clear();
		trial_outlet_concentrations_.clear();
		trial_inlet_mode_ = TrialInletMode::None;
		trial_diagnostics_ = {};
		trial_diagnostics_.planned_configured_substeps = configured_substeps;
		trial_configured_open_loop_schedule_.clear();
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
			trial_inlet_mode_ = TrialInletMode::Coupled;
			return;
		}
		const int node = OutletNode(port_id);
		if (!input.mean_pressure_pa || input.outward_flow_m3_s
			|| input.mean_normal_traction_pa || input.total_pressure_pa
			|| !input.outward_species_flux.empty())
			throw std::runtime_error(
				"native 1d outlet trial input requires mean_pressure_pa and optional concentrations");
		for (const auto& concentration : input.concentration)
			if (!FindOneDSpecies(transports_, concentration.first))
				throw std::runtime_error("1d outlet concentration species '"
					+concentration.first+"' is not transported");
		const auto& outlet = flow_state_.outlets.at(OutletIndex(node));
		if (outlet.kind != OneDOutletKind::Pressure)
			throw std::runtime_error(
				"native 1d outlet pressure cannot replace a resistance or RCR closure");
		trial_outlet_pressure_overrides_[node] = *input.mean_pressure_pa;
		trial_outlet_concentrations_[node] = input.concentration;
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
		trial_inlet_mode_ = TrialInletMode::Coupled;
	}

	void SetOpenLoopInlet(const VascularInletState& inlet)
	{
		RequirePhase(Phase::TrialOpen, "SetOpenLoopInlet");
		ValidateVascularInletState(inlet);
		if (!inlet.has_flow || !Close(inlet.time_s, trial_time_s_+trial_dt_s_))
			throw std::runtime_error("1d open-loop inlet must contain end-time flow data");
		trial_inlet_ = inlet;
		trial_inlet_mode_ = TrialInletMode::HeldOpenLoop;
	}

	void SetConfiguredOpenLoopInlet()
	{
		RequirePhase(Phase::TrialOpen, "SetConfiguredOpenLoopInlet");
		trial_inlet_.reset();
		trial_inlet_mode_ = TrialInletMode::ConfiguredOpenLoop;
		trial_configured_open_loop_schedule_.clear();
		for (int substep = 0; substep < trial_diagnostics_.planned_configured_substeps; ++substep) {
			const double endpoint = substep+1 == trial_diagnostics_.planned_configured_substeps
				? trial_time_s_+trial_dt_s_ : trial_time_s_+(substep+1)*configuration_.time.dt;
			trial_configured_open_loop_schedule_.push_back(OpenLoopInlet(endpoint,
				EvaluateOneDInlet(configuration_, inlet_, case_directory_, endpoint, network_.segments.front().area0)));
		}
	}

	void SolveTrial()
	{
		RequirePhase(Phase::TrialOpen, "SolveTrial");
		if (trial_inlet_mode_ != TrialInletMode::ConfiguredOpenLoop
			&& (!trial_inlet_ || !trial_inlet_->has_flow))
			throw std::runtime_error("1d SolveTrial requires a root flow input");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		trial_diagnostics_.attempted_configured_substeps = 0;
		trial_diagnostics_.completed_configured_substeps = 0;
		trial_diagnostics_.explicit_cfl_substep_delta = 0;
		trial_diagnostics_.configured_open_loop_endpoint_times_s.clear();
		trial_diagnostics_.configured_open_loop_endpoint_flows_m3_s.clear();
		trial_diagnostics_.substep_endpoint_flows_m3_s.clear();
		const long long cfl_before = flow_state_.internal_substeps;
		try {
			double coupled_inlet_flow = 0.0;
			if (trial_inlet_mode_ == TrialInletMode::Coupled)
				coupled_inlet_flow = ApplyOneDCoupledInlet(configuration_, transports_, *trial_inlet_);
			for (auto& transport : transports_)
				for (auto& species : transport.species)
					ResetOneDSpeciesStepAccounting(network_, flow_state_, species);
			for (int substep = 0; substep < trial_diagnostics_.planned_configured_substeps; ++substep) {
				++trial_diagnostics_.attempted_configured_substeps;
				const double sub_start = trial_time_s_+substep*configuration_.time.dt;
				const double sub_end = substep+1 == trial_diagnostics_.planned_configured_substeps
					? trial_time_s_+trial_dt_s_ : trial_time_s_+(substep+1)*configuration_.time.dt;
				for (const auto& override : trial_outlet_pressure_overrides_) {
					auto& outlet = flow_state_.outlets.at(OutletIndex(override.first));
					outlet.pressure = override.second;
					outlet.reference_pressure = override.second;
					outlet.capacitor_pressure = override.second;
				}
				VascularInletState inlet = trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop
					? trial_configured_open_loop_schedule_.at(static_cast<std::size_t>(substep)) : *trial_inlet_;
				if (trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop)
				{
					trial_diagnostics_.configured_open_loop_endpoint_times_s.push_back(sub_end);
					trial_diagnostics_.configured_open_loop_endpoint_flows_m3_s.push_back(inlet.flow_m3_s);
				}
				const double inlet_flow = trial_inlet_mode_ == TrialInletMode::Coupled ? coupled_inlet_flow : inlet.flow_m3_s;
				trial_diagnostics_.substep_endpoint_flows_m3_s.push_back(inlet_flow);
				const auto transport_initial_area = flow_state_.area;
				if (flow_.scheme == OneDFlowScheme::SteadyPoiseuille)
					SolveRigidOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				else if (flow_.scheme == OneDFlowScheme::RigidInertance)
					SolveRigidInertanceOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				else if (flow_.scheme == OneDFlowScheme::ExplicitRusanov)
					AdvanceExplicitOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				else {
					if (!implicit_advance_) throw std::runtime_error("1d implicit trial solve requires an injected PETSc advance function");
					implicit_advance_(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				}
				for (auto& transport : transports_)
					AdvanceOneDTransport(configuration_, network_, flow_state_, transport,
						case_directory_, sub_start, configuration_.time.dt,
						&transport_initial_area, &inlet.species,
						&trial_outlet_concentrations_);
				ApplyOneDVasodilation(configuration_, network_, transports_, configuration_.time.dt, flow_.dynamic_viscosity);
				++flow_state_.completed_step;
				flow_state_.physical_time = sub_end;
				last_inlet_ = std::move(inlet);
				++trial_diagnostics_.completed_configured_substeps;
			}
			trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
			for (auto& transport : transports_)
				for (auto& species : transport.species)
					species.step_accounting_valid = true;
			trial_solve_succeeded_ = true;
			phase_ = Phase::TrialSolved;
		} catch (...) {
			trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
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

	void AbortStep()
	{
		if (phase_ == Phase::Ready) return;
		if (phase_ != Phase::TrialOpen && phase_ != Phase::TrialSolved
			&& phase_ != Phase::CommitPrepared)
			throw std::runtime_error("illegal 1d runtime transition: AbortStep");
		RestoreCommitted();
		CloseTrialState();
		trial_inlet_mode_ = TrialInletMode::None;
		trial_diagnostics_ = {};
		trial_time_s_ = 0.0;
		trial_dt_s_ = 0.0;
		phase_ = Phase::Ready;
	}

	void PrepareCommitStep()
	{
		RequirePhase(Phase::TrialSolved, "PrepareCommitStep");
		if (!trial_solve_succeeded_)
			throw std::runtime_error("1d PrepareCommitStep requires a successful trial solve");
		phase_ = Phase::CommitPrepared;
	}

	void FinalizeCommitStep() noexcept
	{
		if (phase_ != Phase::CommitPrepared) std::terminate();
		CloseTrialState();
		phase_ = Phase::Ready;
	}

	void CommitStep()
	{
		PrepareCommitStep();
		FinalizeCommitStep();
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
				bool root_concentration_supplied = false;
				if (root) {
					const auto boundary = last_inlet_.species.find(species.definition.field);
					root_concentration_supplied = boundary != last_inlet_.species.end();
					double native_inward_flow = 0.0;
					double native_outward_flow = 0.0;
					double outward_weighted_concentration = 0.0;
					for (const int root_segment_index : OneDSegmentsOutOfNode(network_, node)) {
						const auto& root_segment = network_.segments.at(
							static_cast<std::size_t>(root_segment_index));
						const auto root_cell = static_cast<std::size_t>(root_segment.cell_offset);
						const double root_flow = flow_state_.flow.at(root_cell);
						if (root_flow > configuration_.coupling.flow_epsilon_m3_s)
							native_inward_flow += root_flow;
						else if (root_flow < -configuration_.coupling.flow_epsilon_m3_s) {
							native_outward_flow -= root_flow;
							outward_weighted_concentration -= root_flow
								*species.concentration.at(root_cell);
						}
					}
					if (native_inward_flow > 0.0 && native_outward_flow > 0.0)
						throw std::runtime_error(
							"1d aggregate root port cannot report mixed-direction branch flow");
					if (!root_concentration_supplied && native_outward_flow > 0.0)
						concentration = outward_weighted_concentration/native_outward_flow;
					else if (root_concentration_supplied) concentration = boundary->second;
					else concentration = species.inlet_value;
				}
				state.concentration.emplace(species.definition.field, concentration);
				double native_flux = 0.0;
				const double dx = segment.length/segment.cells;
				if (species.boundary_flux_valid) {
					if (root) native_flux = species.root_native_flux;
					else native_flux = species.outlet_native_flux.at(node);
				} else if (root) {
					double boundary_concentration = species.inlet_value;
					const auto boundary = last_inlet_.species.find(species.definition.field);
					if (boundary != last_inlet_.species.end())
						boundary_concentration = boundary->second;
					for (const int root_segment_index : OneDSegmentsOutOfNode(network_, node)) {
						const auto& root_segment = network_.segments.at(
							static_cast<std::size_t>(root_segment_index));
						const auto root_cell = static_cast<std::size_t>(root_segment.cell_offset);
						const double root_flow = flow_state_.flow.at(root_cell);
						const double exterior_concentration = !root_concentration_supplied
							&& root_flow < -configuration_.coupling.flow_epsilon_m3_s
							? species.concentration.at(root_cell) : boundary_concentration;
						native_flux += OneDSpeciesFaceFlux(root_flow,
							exterior_concentration, species.concentration.at(root_cell),
							flow_state_.area.at(root_cell), species.definition.diffusivity,
							root_segment.length/root_segment.cells);
					}
				} else native_flux = OneDSpeciesFaceFlux(native_flow,
					species.concentration.at(cell), species.concentration.at(cell),
					flow_state_.area.at(cell), species.definition.diffusivity, dx);
				state.outward_species_flux.emplace(species.definition.field,
					orientation.ToOutward(native_flux));
			}
		ValidatePortState(state);
		return state;
	}

	std::map<std::string, OneDSpeciesStepAccounting> GetSpeciesStepAccounting() const
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::CommitPrepared
			&& phase_ != Phase::Ready)
			throw std::runtime_error("1d species step accounting requires a solved or committed state");
		std::map<std::string, OneDSpeciesStepAccounting> result;
		for (const auto& transport : transports_)
			for (const auto& species : transport.species)
				result.emplace(species.definition.field,
					GetOneDSpeciesStepAccounting(network_, flow_state_, species));
		return result;
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
	const TrialDiagnostics& Diagnostics() const { return trial_diagnostics_; }

	static PortOrientation RootOrientation() { return {-1}; }

private:
	struct Snapshot {
		OneDConfiguration configuration;
		OneDNetwork network;
		OneDFlowState flow;
		std::vector<OneDTransportState> transports;
		VascularInletState last_inlet;
	};
	static_assert(std::is_nothrow_move_assignable<Snapshot>::value,
		"1D trial snapshot cleanup must remain noexcept");

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

	void CloseTrialState() noexcept
	{
		committed_ = Snapshot{};
		trial_inlet_.reset();
		trial_outlet_pressure_overrides_.clear();
		trial_outlet_concentrations_.clear();
		trial_configured_open_loop_schedule_.clear();
		trial_solve_succeeded_ = false;
	}

	int ConfiguredSubsteps(double macro_dt_s) const
	{
		const long double ratio = static_cast<long double>(macro_dt_s)/configuration_.time.dt;
		const long long rounded = std::llround(ratio);
		if (ratio < 1.0L || rounded < 1 || rounded > std::numeric_limits<int>::max()
			|| std::abs(ratio-rounded) > 1.0e-12L*std::max(1.0L, std::abs(ratio)))
			throw std::runtime_error("1d BeginStep macro dt must be an integer multiple of configured dt");
		return static_cast<int>(rounded);
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
	std::vector<VascularInletState> trial_configured_open_loop_schedule_;
	std::map<int, double> trial_outlet_pressure_overrides_;
	std::map<int, std::map<std::string, double>> trial_outlet_concentrations_;
	TrialInletMode trial_inlet_mode_ = TrialInletMode::None;
	bool trial_solve_succeeded_ = false;
	TrialDiagnostics trial_diagnostics_;
	double trial_time_s_ = 0.0;
	double trial_dt_s_ = 0.0;
	Phase phase_ = Phase::Uninitialized;
};

} // namespace iga

#endif

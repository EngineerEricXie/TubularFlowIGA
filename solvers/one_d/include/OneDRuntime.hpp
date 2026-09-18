#ifndef IGA_ONE_D_RUNTIME_HPP
#define IGA_ONE_D_RUNTIME_HPP

#include "OneDCoupling.hpp"
#include "CouplingPort.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <exception>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

struct OneDFlowCheckpointState {
	std::string configuration_identity_sha256;
	OneDFlowState flow;
	std::vector<OneDTransportState> transports;
	std::vector<double> segment_radii_m;
	VascularInletState last_inlet;
	// Physiology Hct/Hb followed by perfusate oxygen Hct/Hb. Coupled inlet
	// application mutates all four; the last inlet alone cannot reconstruct
	// them when a subsequent inlet omits hematocrit.
	std::array<double, 4> blood_state{};
	int accepted_macro_steps = 0;
	double last_macro_start_s = 0.0, last_macro_dt_s = 0.0;
};

// Owns all mutable native 1d state.  PETSc is deliberately injected for the
// implicit scheme so the common lifecycle and its fast tests remain C++17-only.
class OneDFlowRuntime {
public:
	using ImplicitAdvance = std::function<void(const OneDNetwork&,
		const OneDFlowSystemDefinition&, OneDFlowState&, double, double)>;

	// Optional agreement for LOCAL phases of combined and staged trial paths.
	// The embedding layer supplies MPI; default execution stays C++17-only.
	using FailureAgreement = std::function<void(const char*, std::exception_ptr)>;

	enum class Phase { Uninitialized, Ready, TrialOpen, HydraulicSolved, TrialSolved, CommitPrepared };
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

	// A staged hydraulic trial records the exact flow image consumed by each
	// conservative transport substep.  It deliberately owns copies: transport
	// retries must never recompute or mutate the accepted hydraulic trajectory.
	struct HydraulicFrame {
		std::vector<double> pre_flow_area;
		OneDFlowState post_flow;
		VascularInletState inlet;
		double start_time_s = 0.0;
		double dt_s = 0.0;
	};

	OneDFlowRuntime(OneDConfiguration configuration, OneDFlowSystemDefinition flow,
		OneDNetwork network, OneDInletState inlet, std::filesystem::path case_directory,
		ImplicitAdvance implicit_advance = {}, FailureAgreement failure_agreement = {},
		std::string checkpoint_identity_sha256 = {})
		: configuration_(std::move(configuration)), flow_(std::move(flow)),
		  network_(std::move(network)), inlet_(std::move(inlet)),
		  case_directory_(std::move(case_directory)), implicit_advance_(std::move(implicit_advance)),
		  failure_agreement_(std::move(failure_agreement)), checkpoint_identity_sha256_(std::move(checkpoint_identity_sha256))
	{
		// The embedding provider binds the verified configuration, selected
		// system, network and external inputs before constructing this owner.
		// Legacy callers can omit the identity but cannot use checkpoint APIs.
		if (!checkpoint_identity_sha256_.empty() && (checkpoint_identity_sha256_.size() != 64
			|| checkpoint_identity_sha256_.find_first_not_of("0123456789abcdef") != std::string::npos))
			throw std::runtime_error("invalid 1d checkpoint configuration identity");
		flow_state_.outlets = ResolveOneDOutlets(configuration_, network_);
		for (const auto& transport : configuration_.transport_systems)
			if (transport.flow_system == flow_.name)
				transports_.push_back(InitializeOneDTransport(configuration_, transport, network_));
		for (const auto& transport : transports_)
			for (const auto& species : transport.species)
				configured_species_waveforms_.emplace(species.definition.field, species.inlet_waveform);
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
		if (accepted_macro_steps_ == std::numeric_limits<int>::max())
			throw std::runtime_error("1d macro-step counter overflows");
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
		RunLocalTrialStage("1d trial preparation", [&] {
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
		});
		const long long cfl_before = flow_state_.internal_substeps;
		try {
			double coupled_inlet_flow = 0.0;
			RunLocalTrialStage("1d trial accounting", [&] {
				if (trial_inlet_mode_ == TrialInletMode::Coupled)
					coupled_inlet_flow = ApplyOneDCoupledInlet(configuration_, transports_, *trial_inlet_);
				for (auto& transport : transports_)
					for (auto& species : transport.species)
						ResetOneDSpeciesStepAccounting(network_, flow_state_, species);
			});
			for (int substep = 0; substep < trial_diagnostics_.planned_configured_substeps; ++substep) {
				++trial_diagnostics_.attempted_configured_substeps;
				const double sub_start = trial_time_s_+substep*configuration_.time.dt;
				const double sub_end = substep+1 == trial_diagnostics_.planned_configured_substeps
					? trial_time_s_+trial_dt_s_ : trial_time_s_+(substep+1)*configuration_.time.dt;
				VascularInletState inlet;
				double inlet_flow = 0.0;
				std::vector<double> transport_initial_area;
				RunLocalTrialStage("1d substep preparation", [&] {
					for (const auto& override : trial_outlet_pressure_overrides_) {
						auto& outlet = flow_state_.outlets.at(OutletIndex(override.first));
						outlet.pressure = override.second;
						outlet.reference_pressure = override.second;
						outlet.capacitor_pressure = override.second;
					}
					inlet = trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop
						? trial_configured_open_loop_schedule_.at(static_cast<std::size_t>(substep)) : *trial_inlet_;
					if (trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop) {
						trial_diagnostics_.configured_open_loop_endpoint_times_s.push_back(sub_end);
						trial_diagnostics_.configured_open_loop_endpoint_flows_m3_s.push_back(inlet.flow_m3_s);
					}
					inlet_flow = trial_inlet_mode_ == TrialInletMode::Coupled ? coupled_inlet_flow : inlet.flow_m3_s;
					trial_diagnostics_.substep_endpoint_flows_m3_s.push_back(inlet_flow);
					transport_initial_area = flow_state_.area;
					if ((flow_.scheme == OneDFlowScheme::ImplicitPetsc
						|| flow_.scheme == OneDFlowScheme::Petsc) && !implicit_advance_)
						throw std::runtime_error("1d implicit trial solve requires an injected PETSc advance function");
				});
				if (flow_.scheme == OneDFlowScheme::ImplicitPetsc
					|| flow_.scheme == OneDFlowScheme::Petsc) {
					// This callback may contain collectives; never run it inside a local stage.
					implicit_advance_(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				} else RunLocalTrialStage("1d local flow solve", [&] {
					if (flow_.scheme == OneDFlowScheme::SteadyPoiseuille)
						SolveRigidOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
					else if (flow_.scheme == OneDFlowScheme::RigidInertance)
						SolveRigidInertanceOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
					else if (flow_.scheme == OneDFlowScheme::ExplicitRusanov)
						AdvanceExplicitOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				});
				RunLocalTrialStage("1d substep transport", [&] {
					for (auto& transport : transports_)
						AdvanceOneDTransport(configuration_, network_, flow_state_, transport,
							case_directory_, sub_start, configuration_.time.dt,
							&transport_initial_area, &inlet.species, &trial_outlet_concentrations_);
					ApplyOneDVasodilation(configuration_, network_, transports_, configuration_.time.dt, flow_.dynamic_viscosity);
					++flow_state_.completed_step;
					flow_state_.physical_time = sub_end;
					last_inlet_ = std::move(inlet);
					++trial_diagnostics_.completed_configured_substeps;
				});
			}
			RunLocalTrialStage("1d trial completion", [&] {
				trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
				for (auto& transport : transports_)
					for (auto& species : transport.species)
						species.step_accounting_valid = true;
				trial_solve_succeeded_ = true;
				phase_ = Phase::TrialSolved;
			});
		} catch (...) {
			trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
			trial_solve_succeeded_ = false;
			phase_ = Phase::TrialSolved;
			throw;
		}
	}

	// Staged flow/transport support is intentionally separate from SolveTrial:
	// flow-only callers retain the established combined numerical path.
	void SolveHydraulicTrial()
	{
		RunLocalTrialStage("1d hydraulic preparation", [&] {
			RequirePhase(Phase::TrialOpen, "SolveHydraulicTrial");
			if (configuration_.physiology.vasodilation)
				throw std::runtime_error("staged 1d transport does not support concentration-driven vasodilation");
			if (trial_inlet_mode_ != TrialInletMode::ConfiguredOpenLoop
				&& (!trial_inlet_ || !trial_inlet_->has_flow))
				throw std::runtime_error("1d staged hydraulic solve requires a root flow input");
			RestoreCommitted();
			trial_solve_succeeded_ = false;
			hydraulic_frames_.clear();
			trial_diagnostics_.attempted_configured_substeps = 0;
			trial_diagnostics_.completed_configured_substeps = 0;
			trial_diagnostics_.explicit_cfl_substep_delta = 0;
			trial_diagnostics_.configured_open_loop_endpoint_times_s.clear();
			trial_diagnostics_.configured_open_loop_endpoint_flows_m3_s.clear();
			trial_diagnostics_.substep_endpoint_flows_m3_s.clear();
		});
		const long long cfl_before = flow_state_.internal_substeps;
		try {
			const double coupled_inlet_flow = trial_inlet_mode_ == TrialInletMode::Coupled
				? trial_inlet_->flow_m3_s : 0.0;
			for (int substep = 0; substep < trial_diagnostics_.planned_configured_substeps; ++substep) {
				++trial_diagnostics_.attempted_configured_substeps;
				const double sub_start = trial_time_s_+substep*configuration_.time.dt;
				const double sub_end = substep+1 == trial_diagnostics_.planned_configured_substeps
					? trial_time_s_+trial_dt_s_ : trial_time_s_+(substep+1)*configuration_.time.dt;
				HydraulicFrame frame;
				VascularInletState inlet;
				double inlet_flow = 0.0;
				RunLocalTrialStage("1d hydraulic substep preparation", [&] {
					for (const auto& override : trial_outlet_pressure_overrides_) {
						auto& outlet = flow_state_.outlets.at(OutletIndex(override.first));
						outlet.pressure = override.second;
						outlet.reference_pressure = override.second;
						outlet.capacitor_pressure = override.second;
					}
					inlet = trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop
						? trial_configured_open_loop_schedule_.at(static_cast<std::size_t>(substep)) : *trial_inlet_;
					if (trial_inlet_mode_ == TrialInletMode::ConfiguredOpenLoop) {
						trial_diagnostics_.configured_open_loop_endpoint_times_s.push_back(sub_end);
						trial_diagnostics_.configured_open_loop_endpoint_flows_m3_s.push_back(inlet.flow_m3_s);
					}
					inlet_flow = trial_inlet_mode_ == TrialInletMode::Coupled ? coupled_inlet_flow : inlet.flow_m3_s;
					trial_diagnostics_.substep_endpoint_flows_m3_s.push_back(inlet_flow);
					frame.pre_flow_area = flow_state_.area;
					frame.inlet = inlet;
					frame.start_time_s = sub_start;
					frame.dt_s = configuration_.time.dt;
					if ((flow_.scheme == OneDFlowScheme::ImplicitPetsc
						|| flow_.scheme == OneDFlowScheme::Petsc) && !implicit_advance_)
						throw std::runtime_error("1d implicit trial solve requires an injected PETSc advance function");
				});
				if (flow_.scheme == OneDFlowScheme::ImplicitPetsc
					|| flow_.scheme == OneDFlowScheme::Petsc) {
					// Group collectives stay outside the local work callback.
					implicit_advance_(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				} else RunLocalTrialStage("1d hydraulic local solve", [&] {
					if (flow_.scheme == OneDFlowScheme::SteadyPoiseuille)
						SolveRigidOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
					else if (flow_.scheme == OneDFlowScheme::RigidInertance)
						SolveRigidInertanceOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
					else if (flow_.scheme == OneDFlowScheme::ExplicitRusanov)
						AdvanceExplicitOneD(network_, flow_, flow_state_, inlet_flow, configuration_.time.dt);
				});
				RunLocalTrialStage("1d hydraulic frame capture", [&] {
					++flow_state_.completed_step;
					flow_state_.physical_time = sub_end;
					frame.post_flow = flow_state_;
					hydraulic_frames_.push_back(std::move(frame));
					last_inlet_ = std::move(inlet);
					++trial_diagnostics_.completed_configured_substeps;
				});
			}
			RunLocalTrialStage("1d hydraulic completion", [&] {
				trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
				trial_solve_succeeded_ = true;
				phase_ = Phase::HydraulicSolved;
			});
		} catch (...) {
			trial_solve_succeeded_ = false;
			trial_diagnostics_.explicit_cfl_substep_delta = flow_state_.internal_substeps-cfl_before;
			phase_ = Phase::HydraulicSolved;
			throw;
		}
	}

	void SolveStagedTransportTrial(const std::map<std::string, double>& root_concentrations,
		const std::map<int, std::map<std::string, double>>& outlet_concentrations,
		OneDStagedRootTransportOwnership root_ownership
			= OneDStagedRootTransportOwnership::Legacy)
	{
		// Keep accepted hydraulics and the committed scalar image untouched
		// until every group member has finished the entire local replay.
		std::vector<OneDTransportState> trial_transports;
		VascularInletState trial_last_inlet;
		RunLocalTrialStage("1d staged transport", [&] {
			RequirePhase(Phase::HydraulicSolved, "SolveStagedTransportTrial");
			if (!trial_solve_succeeded_ || hydraulic_frames_.empty())
				throw std::runtime_error("1d staged transport requires a successful hydraulic trial");
			trial_transports = committed_.transports;
			OneDFlowState initial_flow = hydraulic_frames_.front().post_flow;
			initial_flow.area = hydraulic_frames_.front().pre_flow_area;
			for (auto& transport : trial_transports)
				for (auto& species : transport.species)
					ResetOneDSpeciesStepAccounting(network_, initial_flow, species);
			for (const auto& frame : hydraulic_frames_) {
				const auto& frame_root_concentrations = trial_inlet_mode_
					== TrialInletMode::ConfiguredOpenLoop ? frame.inlet.species : root_concentrations;
				for (auto& transport : trial_transports)
					AdvanceOneDTransport(configuration_, network_, frame.post_flow, transport,
						case_directory_, frame.start_time_s, frame.dt_s, &frame.pre_flow_area,
						&frame_root_concentrations, &outlet_concentrations, root_ownership);
			}
			for (auto& transport : trial_transports)
				for (auto& species : transport.species)
					species.step_accounting_valid = true;
			trial_last_inlet = hydraulic_frames_.back().inlet;
			if (trial_inlet_mode_ != TrialInletMode::ConfiguredOpenLoop)
				trial_last_inlet.species = root_concentrations;
		});
		static_assert(std::is_nothrow_move_assignable<VascularInletState>::value,
			"staged inlet publication must not throw after group agreement");
		transports_.swap(trial_transports);
		last_inlet_ = std::move(trial_last_inlet);
		staged_root_transport_ownership_ = root_ownership;
		phase_ = Phase::TrialSolved;
	}

	void RollbackHydraulicTrial()
	{
		RunLocalTrialStage("1d hydraulic rollback", [&] {
			RequirePhase(Phase::HydraulicSolved, "RollbackHydraulicTrial");
			RestoreCommitted();
			hydraulic_frames_.clear();
			trial_inlet_.reset();
			trial_outlet_pressure_overrides_.clear();
			trial_outlet_concentrations_.clear();
			trial_inlet_mode_ = TrialInletMode::None;
			trial_solve_succeeded_ = false;
			staged_root_transport_ownership_ = OneDStagedRootTransportOwnership::Legacy;
			phase_ = Phase::TrialOpen;
		});
	}

	void RollbackStagedTransportTrial()
	{
		RunLocalTrialStage("1d staged transport rollback", [&] {
			RequirePhase(Phase::TrialSolved, "RollbackStagedTransportTrial");
			transports_ = committed_.transports;
			last_inlet_ = hydraulic_frames_.back().inlet;
			staged_root_transport_ownership_ = OneDStagedRootTransportOwnership::Legacy;
			phase_ = Phase::HydraulicSolved;
		});
	}

	void RollbackTrial()
	{
		RequirePhase(Phase::TrialSolved, "RollbackTrial");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		staged_root_transport_ownership_ = OneDStagedRootTransportOwnership::Legacy;
		phase_ = Phase::TrialOpen;
	}

	void AbortStep()
	{
		if (phase_ == Phase::Ready) return;
		if (phase_ != Phase::TrialOpen && phase_ != Phase::HydraulicSolved && phase_ != Phase::TrialSolved
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
		++accepted_macro_steps_;
		last_macro_start_s_ = trial_time_s_; last_macro_dt_s_ = trial_dt_s_;
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
		if (phase_ != Phase::HydraulicSolved && phase_ != Phase::TrialSolved && phase_ != Phase::Ready)
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
					const bool explicit_interior_donor = phase_ == Phase::TrialSolved
						&& staged_root_transport_ownership_
							== OneDStagedRootTransportOwnership::InteriorDonor;
					const auto boundary = last_inlet_.species.find(species.definition.field);
					root_concentration_supplied = !explicit_interior_donor
						&& boundary != last_inlet_.species.end();
					double native_inward_flow = 0.0;
					double native_outward_flow = 0.0;
					double outward_weighted_concentration = 0.0;
					for (const int root_segment_index : OneDSegmentsOutOfNode(network_, node)) {
						const auto& root_segment = network_.segments.at(
							static_cast<std::size_t>(root_segment_index));
						const auto root_cell = static_cast<std::size_t>(root_segment.cell_offset);
						const double root_flow = flow_state_.flow.at(root_cell);
						if (explicit_interior_donor) {
							const double weight = std::abs(root_flow);
							native_outward_flow += weight;
							outward_weighted_concentration += weight
								*species.concentration.at(root_cell);
						} else if (root_flow > configuration_.coupling.flow_epsilon_m3_s)
							native_inward_flow += root_flow;
						else if (root_flow < -configuration_.coupling.flow_epsilon_m3_s) {
							native_outward_flow -= root_flow;
							outward_weighted_concentration -= root_flow
								*species.concentration.at(root_cell);
						}
					}
					if (!explicit_interior_donor && native_inward_flow > 0.0 && native_outward_flow > 0.0)
						throw std::runtime_error(
							"1d aggregate root port cannot report mixed-direction branch flow");
					if (explicit_interior_donor) {
						if (native_outward_flow > 0.0)
							concentration = outward_weighted_concentration/native_outward_flow;
					} else if (!root_concentration_supplied && native_outward_flow > 0.0)
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
						const bool interior_donor = (phase_ == Phase::TrialSolved
							&& staged_root_transport_ownership_
								== OneDStagedRootTransportOwnership::InteriorDonor)
							|| (!root_concentration_supplied
								&& root_flow < -configuration_.coupling.flow_epsilon_m3_s);
						const double exterior_concentration = interior_donor
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

	OneDFlowCheckpointState CaptureCheckpointState() const
	{
		RequirePhase(Phase::Ready, "CaptureCheckpointState");
		OneDFlowCheckpointState value;
		value.configuration_identity_sha256 = checkpoint_identity_sha256_;
		value.flow = flow_state_; value.transports = transports_; value.last_inlet = last_inlet_;
		for (const auto& segment : network_.segments) value.segment_radii_m.push_back(segment.radius0);
		value.blood_state = {{configuration_.physiology.hematocrit_percent, configuration_.physiology.hemoglobin_g_dl,
			configuration_.coupling.perfusate.oxygen.hematocrit_percent, configuration_.coupling.perfusate.oxygen.hemoglobin_g_dl}};
		value.accepted_macro_steps = accepted_macro_steps_;
		value.last_macro_start_s = last_macro_start_s_; value.last_macro_dt_s = last_macro_dt_s_;
		ValidateCheckpointState(value); return value;
	}

	// Initialize an unpublished, fresh candidate. No MPI calls are made here;
	// the provider must coordinate all candidate owners before publication.
	void RestoreCheckpointState(OneDFlowCheckpointState value)
	{
		RequirePhase(Phase::Ready, "RestoreCheckpointState");
		if (accepted_macro_steps_ || flow_state_.completed_step)
			throw std::runtime_error("1d checkpoint restore requires a fresh runtime");
		ValidateCheckpointState(value);
		auto configuration = configuration_; auto network = network_;
		for (std::size_t i = 0; i < network.segments.size(); ++i) {
			auto& segment = network.segments[i]; segment.radius0 = value.segment_radii_m[i];
			segment.area0 = OneDPi*segment.radius0*segment.radius0;
			segment.resistance = 8.0*flow_.dynamic_viscosity*segment.length/(OneDPi*std::pow(segment.radius0, 4.0));
			if (!(segment.area0 > 0.0) || !std::isfinite(segment.area0)
				|| !(segment.resistance > 0.0) || !std::isfinite(segment.resistance))
				throw std::runtime_error("1d checkpoint radius produces invalid derived geometry");
		}
		configuration.physiology.hematocrit_percent = value.blood_state[0];
		configuration.physiology.hemoglobin_g_dl = value.blood_state[1];
		configuration.coupling.perfusate.oxygen.hematocrit_percent = value.blood_state[2];
		configuration.coupling.perfusate.oxygen.hemoglobin_g_dl = value.blood_state[3];
		static_assert(std::is_nothrow_move_assignable<OneDConfiguration>::value
			&& std::is_nothrow_move_assignable<OneDNetwork>::value
			&& std::is_nothrow_move_assignable<OneDFlowState>::value
			&& std::is_nothrow_move_assignable<VascularInletState>::value, "1d checkpoint publication must not throw");
		configuration_ = std::move(configuration); network_ = std::move(network);
		flow_state_ = std::move(value.flow); transports_.swap(value.transports); last_inlet_ = std::move(value.last_inlet);
		accepted_macro_steps_ = value.accepted_macro_steps;
		last_macro_start_s_ = value.last_macro_start_s; last_macro_dt_s_ = value.last_macro_dt_s;
	}

	const std::string& CheckpointConfigurationIdentity() const { return checkpoint_identity_sha256_; }
	const OneDConfiguration& Configuration() const { return configuration_; }
	const OneDFlowSystemDefinition& FlowSystem() const { return flow_; }
	const OneDNetwork& Network() const { return network_; }
	const OneDInletState& InletDefinition() const { return inlet_; }
	const OneDFlowState& FlowState() const { return flow_state_; }
	const std::vector<OneDTransportState>& Transports() const { return transports_; }
	const VascularInletState& LastInlet() const { return last_inlet_; }
	Phase CurrentPhase() const { return phase_; }
	const TrialDiagnostics& Diagnostics() const { return trial_diagnostics_; }
	const std::vector<HydraulicFrame>& HydraulicFrames() const { return hydraulic_frames_; }

	std::vector<double> HydraulicFrameRootBranchNativeFlows(std::size_t frame_index) const
	{
		if (frame_index >= hydraulic_frames_.size())
			throw std::runtime_error("1d hydraulic frame index is out of range");
		std::vector<double> result;
		for (const int segment_index : OneDSegmentsOutOfNode(network_, network_.root)) {
			const auto& segment = network_.segments.at(static_cast<std::size_t>(segment_index));
			result.push_back(hydraulic_frames_[frame_index].post_flow.flow.at(
				static_cast<std::size_t>(segment.cell_offset)));
		}
		return result;
	}

	static PortOrientation RootOrientation() { return {-1}; }

	// Embedding adapters use the same agreement for their LOCAL mutations and
	// solve preparation. With agreement enabled, every group member must call
	// stages in the same order. Work must not call a collective runtime method.
	// Read-only port queries remain local and may be used by one rank alone.
	template<class Work>
	void RunLocalAdapterStage(const char* name, Work&& work)
	{
		RunLocalTrialStage(name, std::forward<Work>(work));
	}

private:
	void ValidateCheckpointState(const OneDFlowCheckpointState& value) const
	{
		auto require = [](bool condition, const char* reason) {
			if (!condition) throw std::runtime_error(std::string("1d checkpoint: ")+reason);
		};
		require(!checkpoint_identity_sha256_.empty() && value.configuration_identity_sha256 == checkpoint_identity_sha256_, "configuration identity differs or is unbound");
		require(value.accepted_macro_steps > 0 && std::isfinite(value.last_macro_start_s)
			&& value.last_macro_start_s >= 0.0
			&& std::isfinite(value.last_macro_dt_s) && value.last_macro_dt_s > 0.0, "invalid macro clock");
		const int substeps = ConfiguredSubsteps(value.last_macro_dt_s);
		require(value.flow.completed_step > 0 && value.flow.completed_step <= configuration_.time.steps
			&& static_cast<long long>(value.accepted_macro_steps)*substeps == value.flow.completed_step
			&& std::isfinite(value.flow.physical_time) && value.flow.physical_time == value.last_macro_start_s+value.last_macro_dt_s
			&& Close(value.flow.physical_time, value.accepted_macro_steps*value.last_macro_dt_s)
			&& value.flow.internal_substeps >= 0 && std::isfinite(value.flow.inlet_flow), "inconsistent accepted counters or time");
		auto field = [&](const std::vector<double>& values, std::size_t expected, bool positive) {
			require(values.size() == expected, "field shape differs");
			for (double entry : values) require(std::isfinite(entry) && (!positive || entry > 0.0), "invalid field value");
		};
		field(value.flow.area, network_.cells, true); field(value.flow.flow, network_.cells, false);
		field(value.flow.pressure, network_.cells, false); field(value.flow.node_pressure, network_.nodes.size(), false);
		field(value.flow.segment_flow, network_.segments.size(), false); field(value.segment_radii_m, network_.segments.size(), true);
		for (std::size_t i = 0; i < network_.segments.size(); ++i)
			if (!configuration_.physiology.vasodilation) require(value.segment_radii_m[i] == network_.segments[i].radius0, "fixed radius changed");
		require(value.flow.outlets.size() == flow_state_.outlets.size(), "outlet count differs");
		for (std::size_t i = 0; i < value.flow.outlets.size(); ++i) {
			const auto& actual = value.flow.outlets[i]; const auto& expected = flow_state_.outlets[i];
			require(actual.node == expected.node && actual.kind == expected.kind && actual.resistance == expected.resistance
				&& actual.proximal_resistance == expected.proximal_resistance && actual.distal_resistance == expected.distal_resistance
				&& actual.capacitance == expected.capacitance, "outlet model differs");
			for (double scalar : {actual.pressure, actual.reference_pressure, actual.capacitor_pressure, actual.flow})
				require(std::isfinite(scalar), "nonfinite outlet state");
			if (actual.kind != OneDOutletKind::Pressure) require(actual.reference_pressure == expected.reference_pressure, "outlet reference pressure changed");
		}
		ValidateVascularInletState(value.last_inlet);
		require(value.last_inlet.has_flow && Close(value.last_inlet.time_s, value.flow.physical_time), "last inlet clock differs");
		for (double scalar : value.blood_state) require(std::isfinite(scalar) && scalar >= 0.0, "invalid blood state");
		require(value.blood_state[0] <= 100.0 && value.blood_state[2] <= 100.0, "invalid hematocrit");
		require(value.transports.size() == transports_.size(), "transport system count differs");
		auto outlet_map = [&](const std::map<int, double>& values) {
			require(values.size() == network_.outlet_nodes.size(), "outlet accounting coverage differs");
			for (int node : network_.outlet_nodes) { const auto found = values.find(node);
				require(found != values.end() && std::isfinite(found->second), "invalid outlet accounting entry"); }
		};
		std::set<std::string> species_names;
		for (std::size_t i = 0; i < transports_.size(); ++i) {
			const auto& actual = value.transports[i]; const auto& expected = transports_[i];
			require(actual.name == expected.name && actual.species.size() == expected.species.size(), "transport catalog differs");
			for (std::size_t j = 0; j < actual.species.size(); ++j) {
				const auto& a = actual.species[j]; const auto& b = expected.species[j]; species_names.insert(b.definition.field);
				require(a.definition.field == b.definition.field && a.definition.diffusivity == b.definition.diffusivity
					&& a.definition.reaction_rate == b.definition.reaction_rate && a.definition.volume_source == b.definition.volume_source
					&& a.wall_kind == b.wall_kind && a.wall_value == b.wall_value && a.wall_coefficient == b.wall_coefficient
					&& a.exterior_value == b.exterior_value, "species model differs");
				field(a.concentration, network_.cells, false);
				for (double scalar : {a.inlet_value, a.root_native_flux, a.step_initial_mass, a.step_root_native_amount, a.step_source_amount})
					require(std::isfinite(scalar), "nonfinite species boundary state");
				require(a.inlet_value >= 0.0
					&& (a.inlet_waveform.empty() || a.inlet_waveform == configured_species_waveforms_.at(b.definition.field)), "invalid species inlet state");
				require(a.boundary_flux_valid && a.step_accounting_valid, "species accepted accounting is unavailable");
				outlet_map(a.outlet_native_flux); outlet_map(a.step_outlet_native_amount);
			}
		}
		for (const auto& species : value.last_inlet.species) require(species_names.count(species.first), "unknown last inlet species");
	}
	// Capture errors without converting the work lambda to std::function:
	// such a conversion could allocate before entering the agreement protocol.
	template<class Work>
	void RunLocalTrialStage(const char* name, Work&& work)
	{
		std::exception_ptr error;
		try { std::forward<Work>(work)(); }
		catch (...) { error = std::current_exception(); }
		if (failure_agreement_) failure_agreement_(name, error);
		// Also preserve failure if an embedding callback accidentally returns.
		if (error) std::rethrow_exception(error);
	}

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
		else if (flow_.model == OneDFlowModel::Lumped)
			InitializeLumpedOneD(network_, flow_, flow_state_, inlet_flow);
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
		hydraulic_frames_.clear();
		trial_solve_succeeded_ = false;
		staged_root_transport_ownership_ = OneDStagedRootTransportOwnership::Legacy;
	}

	int ConfiguredSubsteps(double macro_dt_s) const
	{
		const long double ratio = static_cast<long double>(macro_dt_s)/configuration_.time.dt;
		if (!std::isfinite(ratio) || ratio < 1.0L || ratio > std::numeric_limits<int>::max())
			throw std::runtime_error("1d macro dt ratio is outside supported range");
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
	FailureAgreement failure_agreement_;
	std::string checkpoint_identity_sha256_;
	std::map<std::string, std::string> configured_species_waveforms_;
	int accepted_macro_steps_ = 0;
	double last_macro_start_s_ = 0.0, last_macro_dt_s_ = 0.0;
	OneDFlowState flow_state_;
	std::vector<OneDTransportState> transports_;
	VascularInletState last_inlet_;
	Snapshot committed_;
	std::optional<VascularInletState> trial_inlet_;
	std::vector<VascularInletState> trial_configured_open_loop_schedule_;
	std::vector<HydraulicFrame> hydraulic_frames_;
	std::map<int, double> trial_outlet_pressure_overrides_;
	std::map<int, std::map<std::string, double>> trial_outlet_concentrations_;
	TrialInletMode trial_inlet_mode_ = TrialInletMode::None;
	OneDStagedRootTransportOwnership staged_root_transport_ownership_
		= OneDStagedRootTransportOwnership::Legacy;
	bool trial_solve_succeeded_ = false;
	TrialDiagnostics trial_diagnostics_;
	double trial_time_s_ = 0.0;
	double trial_dt_s_ = 0.0;
	Phase phase_ = Phase::Uninitialized;
};

} // namespace iga

#endif

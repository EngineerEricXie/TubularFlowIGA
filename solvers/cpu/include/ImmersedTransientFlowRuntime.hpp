#ifndef IGA_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP

#include "PetscSolverOptions.hpp"

// Fixed-geometry backward-Euler immersed flow transaction.  This is purposefully
// separate from ImmersedStaticFlowRuntime: the fixed entry accepts only an
// identity transition, while the supplied moving entry is target-epoch-local.
#include "MovingCutGeometry.hpp"
#include "ImmersedFlowPort.hpp"
#include "ImmersedNitscheWall.hpp"
#include "ElementAssemblyExecution.hpp"
#include <iostream>

#include <petscksp.h>

#include <algorithm>
#include <array>
#include "PetscPhaseProfile.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ImmersedTransientFlowOptions {
	NavierStokesParameters parameters{1.0, 1.0, 0.0};
	std::vector<int> wall_labels;
	std::vector<ImmersedFlowPortDefinition> ports;
	double wall_gamma0 = 2.0;
	// Optional backward-Euler wall impedance.  Zero is the exact legacy wall.
	double wall_inertial_gamma0 = 0.0;
	NavierStokesBodyForceEvaluator body_force = [](const std::array<double,3>&) { return std::array<double,3>{{0,0,0}}; };
	bool include_pressure_gauge = true;
	PetscInt nonlinear_maximum_iterations = 12, ksp_maximum_iterations = 2000;
	double ksp_relative_tolerance = 1e-10, nonlinear_relative_tolerance = 1e-9, nonlinear_absolute_tolerance = 1e-11;
	double flow_controller_relative_tolerance = 1e-10, flow_controller_absolute_tolerance_m3_s = 1e-15, flow_controller_reference_flow_m3_s = 1e-12;
	double minimum_damping = 1.0/128.0, lu_pivot_shift = 0.0, nonlinear_block_reduction = 0.0;
	std::string solver_options_prefix; // Empty selects immersed_transient_.
	// Work-only control; never changes the frozen physical input identity.
	bool reuse_accepted_line_search_assembly = true;
	bool line_search_residual_only = true;
	bool cache_volume_basis = true;
	// P5: reuse the previous current-epoch LU only as a right preconditioner for
	// the current Jacobian.  A failed quality gate rebuilds from the current
	// Jacobian once; no factor handle crosses a runtime/geometry epoch.
	bool reuse_preconditioner = true;
	PetscInt reused_preconditioner_maximum_iterations = 30;
	double reused_preconditioner_maximum_iteration_factor = 4.0;
	double reused_preconditioner_true_linear_relative_tolerance = 1e-10;
};

struct ImmersedTransientFlowNewtonStep {
	PetscInt iteration = 0, ksp_iterations = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	double residual_norm = 0.0, update_norm = 0.0, linear_relative_residual = 0.0, damping = 0.0;
};

struct ImmersedTransientFlowConservationDiagnostics {
	std::map<int,double> surface_flow_by_boundary_label_m3_s;
	// Material velocities are sampled through immutable canonical provenance on
	// every target surface label.  Q_w is therefore a closed-boundary quantity;
	// the wall-only view is retained solely for u-w leakage reporting.
	std::map<int,double> material_surface_outward_flow_by_boundary_label_m3_s;
	std::map<int,double> material_wall_outward_flow_by_boundary_label_m3_s;
	// Retain pre-cancellation summation data for the algebraic roundoff check.
	double absolute_surface_flux_sum_m3_s = 0.0;
	std::uint64_t surface_flux_term_count = 0;
	double endpoint_volume_divergence_m3_s = 0.0, total_surface_outward_flow_m3_s = 0.0;
	double open_port_outward_flow_m3_s = 0.0, wall_outward_flow_m3_s = 0.0;
	double total_material_surface_outward_flow_m3_s = 0.0;
	double total_material_wall_outward_flow_m3_s = 0.0, wall_relative_leakage_m3_s = 0.0;
	// R_div = int_Omega div(u) - int_boundary u.n: a volume/surface quadrature
	// consistency measurement.  The historical normalized_open_balance keeps
	// this exact meaning and is not a discrete moving-wall continuity residual.
	double divergence_theorem_defect_m3_s = 0.0;
	// R_cont = Q_port + Q_w,wall uses every configured port and only configured
	// material-wall fluxes at the target endpoint.
	double discrete_moving_wall_continuity_defect_m3_s = 0.0;
	double discrete_moving_wall_continuity_normalization_scale_m3_s = 1.0;
	double normalized_open_balance = 0.0, normalized_wall_leakage = 0.0;
	double normalized_discrete_moving_wall_continuity_defect = 0.0;
};

namespace immersed_transient_detail {

// This identity is evaluated from separately reduced aggregates:
// (Q_port + Q_w,material) + (Q_w,fluid - Q_w,material) = Q_total,fluid.
// Its tolerance retains raw quadrature terms that can cancel even within a
// single aggregate, plus the number of accumulated terms. The allowance is roundoff-only;
// it is not a physical continuity tolerance.
inline double MovingWallContinuityIdentityRoundoffTolerance(
	double open_port_outward_flow_m3_s,
	double wall_outward_flow_m3_s,
	double total_material_wall_outward_flow_m3_s,
	double discrete_moving_wall_continuity_defect_m3_s,
	double wall_relative_leakage_m3_s,
	double total_fluid_surface_outward_flow_m3_s,
	double absolute_surface_flux_sum_m3_s = 0.0,
	std::uint64_t surface_flux_term_count = 0)
{
	for(const double value : {open_port_outward_flow_m3_s,wall_outward_flow_m3_s,
		total_material_wall_outward_flow_m3_s,discrete_moving_wall_continuity_defect_m3_s,
		wall_relative_leakage_m3_s,total_fluid_surface_outward_flow_m3_s,absolute_surface_flux_sum_m3_s})
		if(!std::isfinite(value)) return std::numeric_limits<double>::quiet_NaN();
	if(absolute_surface_flux_sum_m3_s<0.0) return std::numeric_limits<double>::quiet_NaN();
	const double accumulated_epsilon=static_cast<double>(surface_flux_term_count)*std::numeric_limits<double>::epsilon();
	if(accumulated_epsilon>=0.5) return std::numeric_limits<double>::quiet_NaN();
	// gamma_n bounds sequential accumulation.  The extra factor accounts for
	// the separate total/subset sums and the rounded absolute sum itself.
	const double gamma=accumulated_epsilon/(1.0-accumulated_epsilon);
	const double raw_term_scale_m3_s=std::max({std::numeric_limits<double>::min(),
		absolute_surface_flux_sum_m3_s,
		std::abs(open_port_outward_flow_m3_s),
		std::abs(wall_outward_flow_m3_s),
		std::abs(total_material_wall_outward_flow_m3_s),
		std::abs(discrete_moving_wall_continuity_defect_m3_s),
		std::abs(wall_relative_leakage_m3_s),
		std::abs(total_fluid_surface_outward_flow_m3_s)});
	return (128.0*std::numeric_limits<double>::epsilon()+4.0*gamma/(1.0-gamma))*raw_term_scale_m3_s;
}

inline bool MovingWallContinuityIdentityReconciles(
	double open_port_outward_flow_m3_s,
	double wall_outward_flow_m3_s,
	double total_material_wall_outward_flow_m3_s,
	double discrete_moving_wall_continuity_defect_m3_s,
	double wall_relative_leakage_m3_s,
	double total_fluid_surface_outward_flow_m3_s,
	double absolute_surface_flux_sum_m3_s = 0.0,
	std::uint64_t surface_flux_term_count = 0)
{
	const double tolerance=MovingWallContinuityIdentityRoundoffTolerance(open_port_outward_flow_m3_s,
		wall_outward_flow_m3_s,total_material_wall_outward_flow_m3_s,
		discrete_moving_wall_continuity_defect_m3_s,wall_relative_leakage_m3_s,total_fluid_surface_outward_flow_m3_s,
		absolute_surface_flux_sum_m3_s,surface_flux_term_count);
	const double reconstructed_total_fluid_surface_outward_flow_m3_s=discrete_moving_wall_continuity_defect_m3_s+wall_relative_leakage_m3_s;
	if(!std::isfinite(tolerance)||!std::isfinite(reconstructed_total_fluid_surface_outward_flow_m3_s)) return false;
	const double residual_m3_s=reconstructed_total_fluid_surface_outward_flow_m3_s-total_fluid_surface_outward_flow_m3_s;
	return std::isfinite(residual_m3_s)&&std::abs(residual_m3_s)<=tolerance;
}

} // namespace immersed_transient_detail

struct ImmersedTransientFlowDiagnostics {
	std::size_t active_nodes = 0, physical_dofs = 0, total_dofs = 0, volume_cells = 0, surface_cells = 0, ghost_faces = 0;
	double target_time_s = 0.0, dt_s = 0.0, pressure_measure = 0.0, pressure_gauge_defect = 0.0;
	double residual_norm = 0.0, true_linear_relative_residual = 0.0;
	PetscInt nonlinear_iterations = 0, ksp_iterations = 0;
	std::size_t identity_history_nodes = 0, committed_history_nodes = 0, extended_history_nodes = 0, missing_history_nodes = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	bool idle = true, trial_active = false, converged = false, prepared = false, committed = true, scalar_diagonal_structure_verified = false;
	std::size_t attempt_count = 0, abort_count = 0, rollback_count = 0, prepare_count = 0, finalize_count = 0, commit_count = 0;
	double last_assembly_seconds = 0.0, last_linear_solve_seconds = 0.0;
	ImmersedNitscheWallDiagnostics wall_penalty;
	std::string geometry_identity_sha256, layout_hash_sha256, committed_state_hash_sha256, trial_state_hash_sha256, history_hash_sha256;
	// Deterministic publication identities deliberately exclude timing and lifetime
	// attempt/abort counters.  They make an exact retry auditable.
	std::string input_hash_sha256, solved_state_hash_sha256, prepared_hash_sha256, attempt_hash_sha256, moving_map_identity_sha256;
	std::size_t attempt_assembly_count = 0;
	std::size_t attempt_full_assembly_count = 0, attempt_residual_only_count = 0;
	bool volume_basis_cache_enabled = false;
	std::size_t volume_basis_cache_build_misses = 0, volume_basis_cache_hits = 0;
	std::size_t volume_basis_cache_misses = 0, volume_basis_cache_bytes = 0;
	std::string volume_basis_cache_key_sha256;
	bool preconditioner_reuse_enabled = false;
	std::size_t preconditioner_builds = 0, preconditioner_reuse_attempts = 0;
	std::size_t preconditioner_reuse_accepts = 0, preconditioner_rebuilds = 0;
	std::size_t preconditioner_rejections_ksp = 0, preconditioner_rejections_iterations = 0;
	std::size_t preconditioner_rejections_true_residual = 0;
	// Measurements describe a particular published state.  Idle mutation clears
	// them instead of allowing a previous committed state to leak into an input.
	struct Port { std::string id; int boundary_label = -1; ImmersedFlowPortControlMode control_mode = ImmersedFlowPortControlMode::Pressure; double target = 0.0, multiplier = 0.0, controller_error = 0.0; PetscInt multiplier_row = -1; ImmersedFlowPortMeasurement measurement{}; bool measurement_valid = false; };
	std::vector<Port> ports;
	std::vector<ImmersedTransientFlowNewtonStep> newton_steps;
};

class ImmersedTransientFlowRuntime {
public:
	ImmersedTransientFlowRuntime(const MovingCutGeometry& geometry, ImmersedTransientFlowOptions options = {},
		std::shared_ptr<PetscSolverOptions> solver_options = {})
		: solver_options_(std::move(solver_options)), geometry_(geometry), domain_(geometry.Domain()), volume_(geometry.Volume()), surface_(geometry.Surface()), ghost_(geometry.Ghost()), options_(std::move(options))
	{
		PhaseScope geometry_phase(ProfilePhase::Geometry);
		// A constructor whose body throws does not run this object's destructor.
		// Keep all work following possible PETSc handle creation inside this guard.
		try {
		ValidateOptions(); diagnostics_.preconditioner_reuse_enabled=options_.reuse_preconditioner;
		ValidateGeometry(); ConfigurePorts(); PreflightCatalogs(); BuildVolumeBasisCache();
		layout_ = ImmersedActiveLayout::Build(domain_, volume_, geometry_.GeometryIdentitySha256(), ControllerIds(), HasGauge());
		for (auto& port : diagnostics_.ports)
			if (port.control_mode == ImmersedFlowPortControlMode::FlowRate)
				port.multiplier_row = static_cast<PetscInt>(layout_.ControllerRow(static_cast<std::uint64_t>(port.boundary_label)));
		if (layout_.NodeIds().empty()) throw std::runtime_error("immersed transient active layout is empty");
		diagnostics_.active_nodes=layout_.NodeIds().size(); diagnostics_.physical_dofs=layout_.NodeFieldRows(); diagnostics_.total_dofs=layout_.Rows();
		diagnostics_.geometry_identity_sha256=geometry_.GeometryIdentitySha256(); diagnostics_.layout_hash_sha256=layout_.HashSha256();
		BuildGaugeWeights(); CreatePetsc();
		if (!solver_options_) solver_options_.reset(new PetscSolverOptions(PETSC_COMM_SELF,
			options_.solver_options_prefix.empty() ? "immersed_transient_" : options_.solver_options_prefix,
			nullptr, {}, {}, "immersed_transient_", false));
		if (solver_options_->Prefix() != (options_.solver_options_prefix.empty() ? "immersed_transient_" : options_.solver_options_prefix))
			throw std::invalid_argument("transient shared solver prefix differs from runtime options");
		solver_options_->Attach(ksp_);
		solver_options_->Attach(jacobian_);
		solver_options_->Call("immersed transient solver options", [&] { return KSPSetFromOptions(ksp_); });
		solver_options_->RecordUsed();
		std::vector<std::array<double,4>> zero(layout_.NodeIds().size());
		committed_global_.reset(new ImmersedGlobalFlowState(geometry_.Evaluation().EvaluatedTimeS(), 0, layout_, std::move(zero), InitialPortMultipliers(), HasGauge(), 0.0));
		SetVectorFromGlobal(*committed_global_, committed_); Check(VecCopy(committed_, state_), "VecCopy initial state");
		SyncPortMultipliers(diagnostics_.ports, committed_); RefreshHashes();
		} catch (...) { Destroy(); throw; }
	}
	PetscKspConfiguration SolverConfiguration() const { return CaptureKspConfiguration(ksp_); }
	~ImmersedTransientFlowRuntime() { Destroy(); }
	ImmersedTransientFlowRuntime(const ImmersedTransientFlowRuntime&) = delete;
	ImmersedTransientFlowRuntime& operator=(const ImmersedTransientFlowRuntime&) = delete;

	const ImmersedActiveLayout& Layout() const noexcept { return layout_; }
	const std::vector<int>& ConfiguredWallLabels() const noexcept { return options_.wall_labels; }
	const std::vector<ImmersedFlowPortDefinition>& ConfiguredPorts() const noexcept { return options_.ports; }
	const ImmersedTransientFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const ImmersedGlobalFlowState& CommittedGlobalState() const { return *committed_global_; }
	std::vector<PetscScalar> CommittedState() const { return Copy(committed_); }
	std::vector<PetscScalar> TrialState() const { return Copy(state_); }
	PetscInt Dof(std::int32_t id, int field) const { if(field<0||field>3) throw std::out_of_range("immersed transient field is invalid"); return static_cast<PetscInt>(4*layout_.LocalNode(id)+field); }
	PetscInt PortMultiplierDof(const std::string& id) const { for(const auto& p:diagnostics_.ports) if(p.id==id) return p.multiplier_row; throw std::out_of_range("immersed transient port id is absent"); }
	PetscInt GaugeDof() const { if(!HasGauge()) throw std::logic_error("immersed transient gauge is absent"); return static_cast<PetscInt>(layout_.GaugeRow()); }
	void FailNextPrepareForTesting() noexcept { fail_next_prepare_ = true; }
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	void SetVolumeProbeForTesting(std::function<void(std::uint64_t)> probe) { volume_probe_for_testing_=std::move(probe); }
	void SetLineSearchNormProbeForTesting(std::function<void(PetscInt,double,double&,double)> probe) { line_search_norm_probe_for_testing_=std::move(probe); }
	void SetNewtonBoundaryProbeForTesting(std::function<void(PetscInt)> probe) { newton_boundary_probe_for_testing_=std::move(probe); }
	class FirstAssemblyFailureScopeForTesting {
	public:
		FirstAssemblyFailureScopeForTesting() noexcept { FirstAssemblyFaultArmedForTesting()=true; }
		~FirstAssemblyFailureScopeForTesting() { FirstAssemblyFaultArmedForTesting()=false; }
		bool Consumed() const noexcept { return !FirstAssemblyFaultArmedForTesting(); }
		FirstAssemblyFailureScopeForTesting(const FirstAssemblyFailureScopeForTesting&)=delete;
	};
	// This fault injector exists only in the moving-runtime focused test build.
	// It changes no production code path and is deliberately not a publication
	// or convergence bypass.
	class NonfiniteMaterialVelocityScopeForTesting {
	public:
		NonfiniteMaterialVelocityScopeForTesting() noexcept { InjectNonfiniteMaterialVelocityForTesting()=true; }
		~NonfiniteMaterialVelocityScopeForTesting() { InjectNonfiniteMaterialVelocityForTesting()=false; }
		NonfiniteMaterialVelocityScopeForTesting(const NonfiniteMaterialVelocityScopeForTesting&) = delete;
		NonfiniteMaterialVelocityScopeForTesting& operator=(const NonfiniteMaterialVelocityScopeForTesting&) = delete;
	};
	class FiniteConservationDefectScopeForTesting {
	public:
		FiniteConservationDefectScopeForTesting() noexcept { InjectFiniteConservationDefectForTesting()=true; }
		~FiniteConservationDefectScopeForTesting() { InjectFiniteConservationDefectForTesting()=false; }
		FiniteConservationDefectScopeForTesting(const FiniteConservationDefectScopeForTesting&) = delete;
		FiniteConservationDefectScopeForTesting& operator=(const FiniteConservationDefectScopeForTesting&) = delete;
	};
	class FiniteDiscreteMovingWallContinuityDefectScopeForTesting {
	public:
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting() noexcept { InjectFiniteDiscreteMovingWallContinuityDefectForTesting()=true; }
		~FiniteDiscreteMovingWallContinuityDefectScopeForTesting() { InjectFiniteDiscreteMovingWallContinuityDefectForTesting()=false; }
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting(const FiniteDiscreteMovingWallContinuityDefectScopeForTesting&) = delete;
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting& operator=(const FiniteDiscreteMovingWallContinuityDefectScopeForTesting&) = delete;
	};
	double MaxSurfaceRelativeVelocityNormForTesting() const { return MaxSurfaceRelativeVelocityNorm(); }
#endif
	// A trial-only warm start for finite-difference and restart workflows.  It
	// deliberately cannot alter the frozen committed state, identity history, or
	// rollback seed.
	void SetTrialState(const std::vector<PetscScalar>& value)
	{
		RequireTrial("set trial state");
		if(value.size()!=layout_.Rows()) throw std::invalid_argument("immersed transient trial state size is invalid");
		for(const auto x:value) if(!std::isfinite(PetscRealPart(x))) throw std::invalid_argument("immersed transient trial state is not finite");
		SetVector(state_,value); ++state_generation_; InvalidateSolved(); RefreshHashes();
	}

	void SetCommittedGlobalState(const ImmersedGlobalFlowState& value)
	{
		RequireIdle("set committed state"); ValidateGlobal(value);
		std::unique_ptr<ImmersedGlobalFlowState> candidate(new ImmersedGlobalFlowState(value));
		SetVectorFromGlobal(*candidate, prepared_); // all throwing PETSc work precedes publication
		Check(VecCopy(prepared_, base_),"VecCopy committed candidate state");
		auto candidate_ports=diagnostics_.ports;
		SyncPortMultipliers(candidate_ports, prepared_);
		InvalidatePortMeasurements(candidate_ports);
		auto candidate_hash=candidate->HashSha256();
		std::swap(committed_,prepared_); std::swap(state_,base_); committed_global_.swap(candidate);
		diagnostics_.ports.swap(candidate_ports);
		diagnostics_.committed_state_hash_sha256.swap(candidate_hash);
		ClearIdlePublicationDiagnostics(); RefreshHashes();
	}
	void InitializeCommittedGlobalState(const ImmersedGlobalFlowState& value) { SetCommittedGlobalState(value); }
	void SetPortControlValue(const std::string& id, double value)
	{
		RequireIdle("set port control"); if(!std::isfinite(value)) throw std::invalid_argument("immersed transient port control is not finite");
		for(std::size_t i=0;i<options_.ports.size();++i) if(options_.ports[i].id==id) { const double old=options_.ports[i].value; options_.ports[i].value=value; try { ValidateAllFlowCompatibility(); } catch(...) { options_.ports[i].value=old; throw; } diagnostics_.ports[i].target=value; InvalidatePortMeasurements(diagnostics_.ports); ClearIdlePublicationDiagnostics(); RefreshHashes(); return; }
		throw std::out_of_range("immersed transient port id is absent");
	}

	// Freeze the exact identity history and warm-start seed.  Geometry movement,
	// extension, an imprecise target time, or a stale committed state is rejected.
	void BeginTrial(double target_time_s, std::uint64_t target_index, double dt_s)
	{
		RequireIdle("begin trial");
		// Build every fallible part of the transaction in local storage.  In
		// particular, a load callback is allowed to throw, and no idle-visible
		// value (including dt) may have changed when that happens.
		auto candidate_parameters=options_.parameters; candidate_parameters.dt=dt_s; ValidateTransientNavierStokesPreflight(candidate_parameters);
		if(committed_global_->Index()==std::numeric_limits<std::uint64_t>::max()) throw std::overflow_error("immersed transient target index overflows uint64");
		ValidateAllFlowCompatibility();
		if(target_index != committed_global_->Index()+1 || target_time_s != CheckedTransientTargetTime(committed_global_->TimeS(),dt_s)
			|| target_time_s != geometry_.Evaluation().EvaluatedTimeS()) throw std::invalid_argument("immersed transient target time/index does not match fixed geometry transition");
		CertifyZeroMaterialWallVelocity();
		std::unique_ptr<ImmersedVelocityHistory> candidate_history(new ImmersedVelocityHistory(BuildIdentityImmersedVelocityHistory(*committed_global_,layout_,target_time_s)));
		auto candidate_trial_ports=diagnostics_.ports, candidate_frozen_ports=diagnostics_.ports;
		auto candidate_forces=FreezeBodyForceCandidate();
		auto candidate_seed=Copy(committed_);
		auto candidate_diagnostics=diagnostics_;
		SetHistoryDiagnostics(candidate_diagnostics,candidate_history.get());
		candidate_diagnostics.target_time_s=target_time_s; candidate_diagnostics.dt_s=dt_s;
		candidate_diagnostics.idle=false; candidate_diagnostics.trial_active=true; candidate_diagnostics.committed=false;
		candidate_diagnostics.converged=false; candidate_diagnostics.prepared=false;
		ResetAttemptWork(candidate_diagnostics); ++candidate_diagnostics.attempt_count;
		candidate_diagnostics.trial_state_hash_sha256=HashVector(committed_);
		candidate_diagnostics.input_hash_sha256=InputHash(candidate_parameters,*candidate_history,target_time_s,target_index,candidate_frozen_ports,candidate_forces,candidate_seed);
		// prepared_ is an idle spare.  Copying the seed there is the final fallible
		// PETSc operation; pointer exchange then publishes it without a VecCopy.
		Check(VecCopy(committed_,prepared_),"VecCopy staged trial seed");
		options_.parameters=candidate_parameters; history_.swap(candidate_history);
		frozen_body_force_.swap(candidate_forces); frozen_seed_.swap(candidate_seed);
		trial_ports_.swap(candidate_trial_ports); frozen_ports_.swap(candidate_frozen_ports);
		std::swap(state_,prepared_); target_time_s_=target_time_s; target_index_=target_index; moving_trial_=false; moving_map_.reset();
		std::swap(diagnostics_,candidate_diagnostics);
		++state_generation_;
	}
	// A target-geometry runtime owns only target-layout PETSc objects.  The
	// caller supplies the immutable old-to-target continuation and pressure seed;
	// no old-layout vector is ever copied into this runtime.
	void BeginMovingTrial(double target_time_s, std::uint64_t target_index, double dt_s,
		const ImmersedVelocityHistory& supplied_history, const ImmersedGlobalFlowState& supplied_seed,
		const ImmersedMovingTrialMapIdentity& map_identity)
	{
		RequireIdle("begin moving trial");
		auto candidate_parameters=options_.parameters; candidate_parameters.dt=dt_s; ValidateTransientNavierStokesPreflight(candidate_parameters);
		ValidateAllFlowCompatibility(); ValidateMovingInputs(target_time_s,target_index,dt_s,supplied_history,supplied_seed,map_identity);
		CertifyFiniteMaterialSurfaceVelocity();
		std::unique_ptr<ImmersedVelocityHistory> candidate_history(new ImmersedVelocityHistory(supplied_history));
		std::unique_ptr<ImmersedMovingTrialMapIdentity> candidate_map(new ImmersedMovingTrialMapIdentity(map_identity));
		auto candidate_trial_ports=diagnostics_.ports, candidate_frozen_ports=diagnostics_.ports;
		auto candidate_forces=FreezeBodyForceCandidate();
		// The exact supplied target-layout seed is copied while publication is
		// still idle.  prepared_ is the target-sized staging vector.
		SetVectorFromGlobal(supplied_seed,prepared_);
		auto candidate_seed=Copy(prepared_);
		auto candidate_diagnostics=diagnostics_;
		SetHistoryDiagnostics(candidate_diagnostics,candidate_history.get());
		candidate_diagnostics.target_time_s=target_time_s; candidate_diagnostics.dt_s=dt_s;
		candidate_diagnostics.idle=false; candidate_diagnostics.trial_active=true; candidate_diagnostics.committed=false;
		candidate_diagnostics.converged=false; candidate_diagnostics.prepared=false;
		ResetAttemptWork(candidate_diagnostics); ++candidate_diagnostics.attempt_count;
		candidate_diagnostics.trial_state_hash_sha256=HashVector(prepared_);
		candidate_diagnostics.moving_map_identity_sha256=map_identity.HashSha256();
		candidate_diagnostics.input_hash_sha256=MovingInputHash(candidate_parameters,*candidate_history,target_time_s,target_index,candidate_frozen_ports,candidate_forces,candidate_seed,map_identity);
		options_.parameters=candidate_parameters; history_.swap(candidate_history); moving_map_.swap(candidate_map);
		frozen_body_force_.swap(candidate_forces); frozen_seed_.swap(candidate_seed);
		trial_ports_.swap(candidate_trial_ports); frozen_ports_.swap(candidate_frozen_ports);
		std::swap(state_,prepared_); target_time_s_=target_time_s; target_index_=target_index; moving_trial_=true;
		std::swap(diagnostics_,candidate_diagnostics);
		++state_generation_;
	}
	void Assemble() { AssembleImpl(NavierStokesAssemblyRequest::ResidualAndJacobian); }
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	void AssembleResidualOnlyForTesting() { AssembleImpl(NavierStokesAssemblyRequest::ResidualOnly); }
#endif
private:
	void AssembleImpl(NavierStokesAssemblyRequest request)
	{
		const bool full=request==NavierStokesAssemblyRequest::ResidualAndJacobian;
		assembly_execution_.RequireCaller();
		PhaseScope assembly_phase(ProfilePhase::Assembly);
		RequireTrial("assemble"); const auto start=std::chrono::steady_clock::now();
		DetailRecord detail("assembly");
		if(detail.Enabled()) {
			detail.String("assembly_request",full?"ResidualAndJacobian":"ResidualOnly");
			detail.Number("step",target_index_); detail.Number("active_dofs",layout_.Rows());
			detail.String("geometry_epoch",geometry_.GeometryIdentitySha256());
			detail.String("input_identity",diagnostics_.input_hash_sha256);
			detail.String("state_identity",HashVector(state_));
			detail.String("history_identity",diagnostics_.history_hash_sha256);
			detail.String("moving_map_identity",diagnostics_.moving_map_identity_sha256);
			detail.Number("state_generation",state_generation_);
			detail.Number("inside_cells",domain_.Diagnostics().inside_count);
			detail.Number("cut_cells",domain_.Diagnostics().cut_count);
			detail.Number("volume_quadrature_points",volume_.Diagnostics().logical_output_points);
			detail.Number("surface_quadrature_points",surface_.Diagnostics().output_points);
			detail.Number("volume_basis_cache_enabled",options_.cache_volume_basis?1:0);
			detail.String("volume_basis_cache_key",diagnostics_.volume_basis_cache_key_sha256);
			detail.Number("volume_basis_cache_build_misses",diagnostics_.volume_basis_cache_build_misses);
			detail.Number("volume_basis_cache_bytes",diagnostics_.volume_basis_cache_bytes);
		}
		if(diagnostics_.converged) InvalidateSolved();
		++diagnostics_.attempt_assembly_count; if(full) ++diagnostics_.attempt_full_assembly_count; else ++diagnostics_.attempt_residual_only_count;
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if(FirstAssemblyFaultArmedForTesting()) {
			FirstAssemblyFaultArmedForTesting()=false;
			throw std::runtime_error("injected first immersed transient assembly failure");
		}
#endif
		// A failed batch can leave earlier MatSetValues contributions pending.
		// PETSc rejects MatZeroEntries until these insertions are assembled.
		// Complete the pending insertion phase, then discard all its values.
		if(assembly_pending_) {
			Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin recovery");
			Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd recovery");
			Check(VecAssemblyBegin(rhs_),"VecAssemblyBegin recovery");
			Check(VecAssemblyEnd(rhs_),"VecAssemblyEnd recovery");
		}
		if(full) Check(MatZeroEntries(jacobian_),"MatZeroEntries");
		Check(VecSet(rhs_,0.0),"VecSet rhs");
		assembly_pending_=true;
		std::size_t candidate_volume_cells=0, candidate_surface_cells=0, candidate_ghost_faces=0;
		ImmersedNitscheWallDiagnostics candidate_wall_penalty;
		std::vector<std::uint64_t> cells;
		for(std::uint64_t cell=0;cell<domain_.Cells().size();++cell) if(Usable(cell)) cells.push_back(cell);
		struct PreparedVolume {
			std::uint64_t cell;
			Element element;
			std::vector<std::array<double,4>> nodal;
		};
        struct ElementSystems {
            NavierStokesSystem volume;
            std::optional<NavierStokesSystem> trace;
            std::optional<ImmersedNitscheWallAssembly> wall;
            double volume_seconds=0.,trace_seconds=0.,wall_seconds=0.;
			std::size_t volume_basis_cache_hits=0,volume_basis_cache_misses=0;
			std::size_t ResultPayloadBytes() const noexcept {
				return volume.ResultPayloadBytes()+(trace?trace->ResultPayloadBytes():0)
					+(wall?wall->system.ResultPayloadBytes():0);
			}
        };
        auto batch_options=assembly_execution_.Options(); batch_options.profile_detail=detail.Enabled();
        double volume_work=0.,trace_work=0.,wall_work=0.;
		std::size_t volume_basis_cache_hits=0,volume_basis_cache_misses=0;
        const auto batch = ForEachElementBatch(cells.size(), batch_options,
            [&](std::size_t index) {
                auto element=domain_.Background().MaterializeElement(cells[index]);
                auto nodal=Gather(element);
                return PreparedVolume{cells[index],std::move(element),std::move(nodal)};
            },
            [this,&detail,request](const PreparedVolume& input, std::size_t) {
                ElementSystems result;
                { DetailTimer timer(result.volume_seconds,detail.Enabled()); result.volume=BuildVolume(input.element,input.nodal,input.cell,request,result.volume_basis_cache_hits,result.volume_basis_cache_misses); }
                if(domain_.Cells()[input.cell].classification==CellClassification::Cut) {
                    const auto& rule=surface_.UsableRule(domain_,input.cell);
                    { DetailTimer timer(result.trace_seconds,detail.Enabled()); result.trace=BuildImmersedConservativeMixedTraceElement(input.element,rule,input.nodal,request); }
                    if(HasWallPoint(rule)) {
                        DetailTimer timer(result.wall_seconds,detail.Enabled());
                        const auto local_history=LocalizeImmersedVelocityHistory(input.element,layout_,*history_,target_time_s_);
                        std::vector<std::array<double,4>> old(local_history.size());
                        for(std::size_t i=0;i<old.size();++i) for(int q=0;q<3;++q) old[i][q]=local_history[i][q];
                        result.wall=BuildImmersedNitscheWallElementFromVolumeSystemMaterialAware(domain_,volume_,surface_,input.cell,input.nodal,old,options_.parameters,options_.wall_labels,result.volume,ghost_,options_.wall_gamma0,options_.wall_inertial_gamma0,MaterialVelocity());
                    }
                }
                return result;
            },
            [&](const PreparedVolume& input, const ElementSystems& systems, std::size_t) {
                volume_work+=systems.volume_seconds; trace_work+=systems.trace_seconds; wall_work+=systems.wall_seconds;
				if(std::numeric_limits<std::size_t>::max()-volume_basis_cache_hits<systems.volume_basis_cache_hits
					||std::numeric_limits<std::size_t>::max()-volume_basis_cache_misses<systems.volume_basis_cache_misses)
					throw std::overflow_error("immersed transient volume basis cache counters overflow");
				volume_basis_cache_hits+=systems.volume_basis_cache_hits;
				volume_basis_cache_misses+=systems.volume_basis_cache_misses;
                // PETSc and diagnostic publication stay on the caller in cell order.
                { auto timer=detail.Time("scatter"); Scatter(input.element.connectivity,systems.volume); } ++candidate_volume_cells;
                if(systems.trace) {
                    { auto timer=detail.Time("scatter"); Scatter(input.element.connectivity,*systems.trace); }
                    if(systems.wall) {
                        AccumulateWallPenaltyDiagnostics(candidate_wall_penalty,systems.wall->diagnostics);
                        { auto timer=detail.Time("scatter"); SubtractAndScatter(input.element.connectivity,systems.wall->system,systems.volume); } ++candidate_surface_cells;
                    }
                    const auto& rule=surface_.UsableRule(domain_,input.cell);
                    auto ports_timer=detail.Time("ports");
                    for(std::size_t p=0;p<options_.ports.size();++p) if(RuleHasLabel(rule,options_.ports[p].boundary_label)) ScatterPort(input.element,input.nodal,rule,p,full);
                }
            });
		detail.Number("volume_worker_work_s",volume_work); detail.Number("trace_worker_work_s",trace_work); detail.Number("wall_worker_work_s",wall_work);
		if(std::numeric_limits<std::size_t>::max()-diagnostics_.volume_basis_cache_hits<volume_basis_cache_hits
			||std::numeric_limits<std::size_t>::max()-diagnostics_.volume_basis_cache_misses<volume_basis_cache_misses)
			throw std::overflow_error("immersed transient lifetime volume basis cache counters overflow");
		diagnostics_.volume_basis_cache_hits+=volume_basis_cache_hits;
		diagnostics_.volume_basis_cache_misses+=volume_basis_cache_misses;
		detail.Number("volume_basis_cache_hits",volume_basis_cache_hits);
		detail.Number("volume_basis_cache_misses",volume_basis_cache_misses);
		detail.Number("prepare_wall_s",batch.prepare_wall_seconds); detail.Number("compute_wall_s",batch.compute_wall_seconds); detail.Number("consume_wall_s",batch.consume_wall_seconds);
		detail.Number("worker_work_s",batch.worker_work_seconds); detail.Number("max_worker_s",batch.maximum_worker_seconds);
		detail.Number("maximum_worker_completion_spread_s",batch.maximum_worker_completion_spread_seconds);
		detail.Number("maximum_resident_result_payload_bytes",batch.maximum_resident_result_payload_bytes);
		detail.Number("team_size",batch.maximum_team_size); detail.Number("batch_capacity",batch_options.capacity);
		if(CurrentPhaseProfile().Enabled())
			std::cout << "element_assembly threads_requested=" << assembly_execution_.Options().threads
				<< " team_size=" << batch.maximum_team_size << " cells=" << batch.items
				<< " batches=" << batch.batches << " maximum_resident_items=" << batch.maximum_resident_items << '\n';
		{ auto timer=detail.Time("ghost"); for(std::size_t f=0;f<ghost_.Faces().size();++f) { const auto block=full?ghost_.AssembleFaceLocal(f,domain_,volume_,[this](std::int32_t node,int q){return Value(state_,Dof(node,q));},options_.parameters.dynamic_viscosity):ghost_.AssembleFaceResidualLocal(f,domain_,volume_,[this](std::int32_t node,int q){return Value(state_,Dof(node,q));},options_.parameters.dynamic_viscosity); ScatterBlock(block.connectivity,block.jacobian,block.negative_residual); ++candidate_ghost_faces; } }
		{ auto timer=detail.Time("mat_vec_assembly");
		if(full) { Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin"); Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd"); } Check(VecAssemblyBegin(rhs_),"VecAssemblyBegin"); Check(VecAssemblyEnd(rhs_),"VecAssemblyEnd");
		}
		{ auto timer=detail.Time("port_measurement_gauge"); MeasurePorts(); if(HasGauge()) InsertGauge(full); }
		{ auto timer=detail.Time("mat_vec_assembly"); if(full) { Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin gauge"); Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd gauge"); } Check(VecAssemblyBegin(rhs_),"VecAssemblyBegin gauge"); Check(VecAssemblyEnd(rhs_),"VecAssemblyEnd gauge"); }
		assembly_pending_=false;
		diagnostics_.volume_cells=candidate_volume_cells; diagnostics_.surface_cells=candidate_surface_cells; diagnostics_.ghost_faces=candidate_ghost_faces;
		diagnostics_.wall_penalty=std::move(candidate_wall_penalty);
		diagnostics_.last_assembly_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
		{ auto timer=detail.Time("hash_diagnostic"); RefreshHashes(); }
		detail.Number("ghost_faces",candidate_ghost_faces);
		residual_generation_=state_generation_; if(full) jacobian_generation_=state_generation_;
		detail.Number("residual_generation",residual_generation_); detail.Number("jacobian_generation",jacobian_generation_);
	}
	struct LinearSolveResult {
		PetscInt iterations=0;
		KSPConvergedReason reason=KSP_CONVERGED_ITERATING;
		double update_norm=0.0;
		double true_relative_residual=std::numeric_limits<double>::infinity();
	};
	LinearSolveResult SolveNewtonLinearSystem(bool reuse_preconditioner)
	{
		DetailRecord detail("linear-solve-attempt");
		detail.Number("reuse_preconditioner",reuse_preconditioner?1:0);
		if(reuse_preconditioner) ++diagnostics_.preconditioner_reuse_attempts;
		else ++diagnostics_.preconditioner_builds;
		Check(KSPSetReusePreconditioner(ksp_,reuse_preconditioner?PETSC_TRUE:PETSC_FALSE),"KSPSetReusePreconditioner");
		Check(KSPSetOperators(ksp_,jacobian_,jacobian_),"KSPSetOperators");
		solver_options_->Call("immersed transient solve", [&] { return SolveProfiledKsp(ksp_,rhs_,update_); });
		solver_options_->RecordUsed();
		LinearSolveResult result;
		Check(KSPGetIterationNumber(ksp_,&result.iterations),"KSPGetIterationNumber");
		Check(KSPGetConvergedReason(ksp_,&result.reason),"KSPGetConvergedReason");
		if(result.reason>0) {
			Check(MatMult(jacobian_,update_,action_output_),"MatMult");
			Check(VecAXPY(action_output_,-1.0,action_input_),"VecAXPY");
			PetscReal linear=0,update=0,rhs_norm=0;
			Check(VecNorm(action_output_,NORM_2,&linear),"VecNorm linear");
			Check(VecNorm(update_,NORM_2,&update),"VecNorm update");
			Check(VecNorm(action_input_,NORM_2,&rhs_norm),"VecNorm linear rhs");
			result.update_norm=update;
			result.true_relative_residual=rhs_norm>0?linear/rhs_norm:linear;
		}
		detail.Number("ksp_iterations",result.iterations);
		detail.Number("ksp_reason",static_cast<int>(result.reason));
		detail.Number("true_linear_relative_residual",result.true_relative_residual);
		return result;
	}
public:
	bool SolveTrial()
	{
		RequireTrial("solve"); ResetAttemptWork(); double initial=-1.0; std::array<double,3> initial_blocks{{0.0,0.0,0.0}};
		// A ready full R/J belongs ONLY to this synchronous SolveTrial invocation.
		// No external setter, new trial, rollback or retry can inherit this flag.
		// Frozen history/map/ports/forces and geometry cannot change in this loop.
		bool accepted_full_ready=false,preconditioner_available=false;
		PetscInt last_fresh_preconditioner_iterations=-1;
		try { for(PetscInt it=0;it<options_.nonlinear_maximum_iterations;++it) {
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
			if(it>0 && newton_boundary_probe_for_testing_) newton_boundary_probe_for_testing_(it);
#endif
			if(!accepted_full_ready) { AssemblyCallScope tag(it==0?"initial":"next-Newton",it); Assemble(); }
			else {
				AssemblyCallScope tag("accepted-full-reuse",it); DetailRecord reuse("assembly-reuse");
				if(reuse.Enabled()) {
					reuse.Number("step",target_index_); reuse.Number("active_dofs",layout_.Rows());
					reuse.Number("state_generation",state_generation_);
					reuse.Number("residual_generation",residual_generation_); reuse.Number("jacobian_generation",jacobian_generation_);
					reuse.String("geometry_epoch",geometry_.GeometryIdentitySha256());
					reuse.String("input_identity",diagnostics_.input_hash_sha256);
					reuse.String("state_identity",HashVector(state_));
					reuse.String("history_identity",diagnostics_.history_hash_sha256);
					reuse.String("moving_map_identity",diagnostics_.moving_map_identity_sha256);
				}
			}
			accepted_full_ready=false; // Consume once; only an accepted full candidate can re-arm.
			if(residual_generation_!=state_generation_ || jacobian_generation_!=state_generation_) throw std::logic_error("immersed transient Newton iteration requires current full R/J");
			PetscReal norm=0; Check(VecNorm(rhs_,NORM_2,&norm),"VecNorm residual"); if(!std::isfinite(norm)) throw std::runtime_error("immersed transient nonlinear residual is not finite"); if(initial<0) { initial=norm; initial_blocks=ResidualBlockNorms(Copy(rhs_)); } diagnostics_.residual_norm=norm; diagnostics_.nonlinear_iterations=it;
			if(norm<=std::max(options_.nonlinear_absolute_tolerance,options_.nonlinear_relative_tolerance*initial) && BlockReductionSatisfied(initial_blocks,Copy(rhs_),initial) && ControllersSatisfied()) { MarkSolved(); return true; }
			Check(VecCopy(rhs_,action_input_),"VecCopy rhs");
			const auto linear_begin=std::chrono::steady_clock::now();
			const bool request_reuse=options_.reuse_preconditioner&&preconditioner_available;
			PetscInt attempted_iterations=0;
			LinearSolveResult linear_result=SolveNewtonLinearSystem(request_reuse);
			attempted_iterations+=linear_result.iterations;
			if(request_reuse) {
				const bool rejected_ksp=linear_result.reason<=0;
				const bool rejected_iterations=linear_result.iterations>options_.reused_preconditioner_maximum_iterations
					||(last_fresh_preconditioner_iterations>0
						&&static_cast<double>(linear_result.iterations)>options_.reused_preconditioner_maximum_iteration_factor*static_cast<double>(last_fresh_preconditioner_iterations));
				const bool rejected_true_residual=!std::isfinite(linear_result.true_relative_residual)
					||linear_result.true_relative_residual>options_.reused_preconditioner_true_linear_relative_tolerance;
				if(rejected_ksp||rejected_iterations||rejected_true_residual) {
					++diagnostics_.preconditioner_rebuilds;
					if(rejected_ksp) ++diagnostics_.preconditioner_rejections_ksp;
					if(rejected_iterations) ++diagnostics_.preconditioner_rejections_iterations;
					if(rejected_true_residual) ++diagnostics_.preconditioner_rejections_true_residual;
					linear_result=SolveNewtonLinearSystem(false);
					attempted_iterations+=linear_result.iterations;
					last_fresh_preconditioner_iterations=linear_result.iterations;
				} else ++diagnostics_.preconditioner_reuse_accepts;
			} else last_fresh_preconditioner_iterations=linear_result.iterations;
			diagnostics_.last_linear_solve_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-linear_begin).count();
			diagnostics_.ksp_iterations+=attempted_iterations;
			diagnostics_.ksp_reason=linear_result.reason;
			if(linear_result.reason<=0) throw std::runtime_error("immersed transient KSP failed");
			preconditioner_available=true;
			ImmersedTransientFlowNewtonStep step; step.iteration=it; step.ksp_iterations=attempted_iterations; step.ksp_reason=diagnostics_.ksp_reason; step.residual_norm=norm; step.update_norm=linear_result.update_norm; step.linear_relative_residual=linear_result.true_relative_residual; diagnostics_.true_linear_relative_residual=step.linear_relative_residual;
			Check(VecCopy(state_,base_),"VecCopy line-search base");
			double d=1; bool accepted=false; while(d>=options_.minimum_damping) {
				accepted_full_ready=false; // Invalidate BEFORE every restore/update, rejected or not.
				Check(VecCopy(base_,state_),"VecCopy line-search restore"); Check(VecAXPY(state_,d,update_),"VecAXPY update"); ++state_generation_; { AssemblyCallScope tag("line-search",it,d); if(options_.line_search_residual_only) AssembleImpl(NavierStokesAssemblyRequest::ResidualOnly); else Assemble(); } PetscReal candidate=0; Check(VecNorm(rhs_,NORM_2,&candidate),"VecNorm candidate");
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
				if(line_search_norm_probe_for_testing_) line_search_norm_probe_for_testing_(it,d,candidate,norm);
#endif
				if(std::isfinite(candidate)&&candidate<norm) { accepted=true; step.damping=d; diagnostics_.residual_norm=candidate; ++diagnostics_.nonlinear_iterations; diagnostics_.newton_steps.push_back(step); if(candidate<=std::max(options_.nonlinear_absolute_tolerance,options_.nonlinear_relative_tolerance*initial) && BlockReductionSatisfied(initial_blocks,Copy(rhs_),initial) && ControllersSatisfied()) { MarkSolved(); return true; } accepted_full_ready=!options_.line_search_residual_only&&options_.reuse_accepted_line_search_assembly; break; }
				d*=.5;
			}
			if(!accepted) { Check(VecCopy(base_,state_),"VecCopy line-search failure restore"); ++state_generation_; }
			if(!accepted) { diagnostics_.newton_steps.push_back(step); throw std::runtime_error("immersed transient backtracking failed"); }
		} throw std::runtime_error("immersed transient nonlinear iteration cap reached"); } catch(...) { accepted_full_ready=false; Rollback(); throw; }
	}
	void PrepareCommit()
	{
		if(!diagnostics_.trial_active||!diagnostics_.converged||diagnostics_.prepared) throw std::logic_error("cannot prepare immersed transient trial");
		if(diagnostics_.solved_state_hash_sha256.empty() || HashVector(state_)!=diagnostics_.solved_state_hash_sha256)
			throw std::logic_error("immersed transient trial state changed after convergence");
		if(fail_next_prepare_) { fail_next_prepare_=false; throw std::runtime_error("injected immersed transient prepare failure"); }
		Check(VecCopy(state_,prepared_),"VecCopy prepare");
		std::unique_ptr<ImmersedGlobalFlowState> global(new ImmersedGlobalFlowState(target_time_s_,target_index_,layout_,ReadCoefficients(prepared_),ReadMultipliers(prepared_),HasGauge(),HasGauge()?Value(prepared_,GaugeDof()):0.0));
		auto ports=trial_ports_; SyncPortMultipliers(ports, prepared_);
		std::string hash=PreparedHash(*global,ports), global_hash=global->HashSha256();
		prepared_global_.swap(global); prepared_ports_.swap(ports); prepared_global_hash_.swap(global_hash);
		diagnostics_.prepared_hash_sha256.swap(hash); diagnostics_.prepared=true; ++diagnostics_.prepare_count;
	}
	void FinalizeCommit() noexcept
	{
		if(!diagnostics_.prepared) return;
		std::swap(committed_,prepared_); committed_global_.swap(prepared_global_); diagnostics_.ports.swap(prepared_ports_);
		diagnostics_.committed_state_hash_sha256.swap(prepared_global_hash_); // already allocated in PrepareCommit
		history_.reset(); moving_map_.reset(); moving_trial_=false; frozen_body_force_.clear(); frozen_seed_.clear(); trial_ports_.clear(); diagnostics_.trial_state_hash_sha256.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.moving_map_identity_sha256.clear(); diagnostics_.prepared=false; diagnostics_.trial_active=false; diagnostics_.idle=true; diagnostics_.committed=true; ++diagnostics_.finalize_count; ++diagnostics_.commit_count;
	}
	void Commit() { PrepareCommit(); FinalizeCommit(); }
	void Rollback()
	{
		if(!diagnostics_.trial_active) return;
		++state_generation_;
		diagnostics_.prepared=false; prepared_global_.reset(); prepared_ports_.clear(); prepared_global_hash_.clear(); diagnostics_.prepared_hash_sha256.clear(); if(!frozen_seed_.empty()) SetVector(state_,frozen_seed_); trial_ports_=frozen_ports_; ResetAttemptWork(); ++diagnostics_.rollback_count; RefreshHashes();
	}
	void AbortTrial()
	{
		if(!diagnostics_.trial_active && !diagnostics_.prepared) return;
		Check(VecCopy(committed_,state_),"VecCopy abort committed state");
		++state_generation_;
		diagnostics_.prepared=false; prepared_global_.reset(); prepared_ports_.clear(); prepared_global_hash_.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.input_hash_sha256.clear(); diagnostics_.moving_map_identity_sha256.clear(); trial_ports_.clear(); history_.reset(); moving_map_.reset(); moving_trial_=false; frozen_body_force_.clear(); frozen_seed_.clear(); diagnostics_.trial_active=false; diagnostics_.idle=true; diagnostics_.committed=true; ResetAttemptWork(); ++diagnostics_.abort_count; RefreshHashes();
	}
	void AbortPrepared() { AbortTrial(); }
	MatInfo JacobianStorageInfo() const { MatInfo info{};Check(MatGetInfo(jacobian_,MAT_LOCAL,&info),"MatGetInfo Jacobian");return info; }
	std::vector<PetscScalar> AssembledNegativeResidual() const { if(residual_generation_!=state_generation_) throw std::logic_error("immersed transient residual does not belong to current state"); return Copy(rhs_); }
	std::vector<PetscScalar> AssembledJacobianAction(const std::vector<PetscScalar>& x) const { if(x.size()!=layout_.Rows()) throw std::invalid_argument("immersed transient Jacobian action size is invalid"); if(jacobian_generation_!=state_generation_) throw std::logic_error("immersed transient Jacobian does not belong to current state"); SetVector(action_input_,x); Check(MatMult(jacobian_,action_input_,action_output_),"MatMult action"); return Copy(action_output_); }
	ImmersedTransientFlowConservationDiagnostics ConservationDiagnostics() const { return MeasureConservation(); }

private:
	static void Check(PetscErrorCode c,const char* op) { if(c) throw std::runtime_error(std::string("PETSc ")+op+" failed: "+std::to_string(static_cast<long long>(c))); }
	void ValidateOptions() const { if(!options_.body_force||!std::isfinite(options_.parameters.density)||!(options_.parameters.density>0)||!std::isfinite(options_.parameters.dynamic_viscosity)||!(options_.parameters.dynamic_viscosity>0)||!std::isfinite(options_.wall_gamma0)||!(options_.wall_gamma0>0)||!std::isfinite(options_.wall_inertial_gamma0)||options_.wall_inertial_gamma0<0||options_.nonlinear_maximum_iterations<=0||options_.ksp_maximum_iterations<=0||!std::isfinite(options_.ksp_relative_tolerance)||!(options_.ksp_relative_tolerance>0)||!std::isfinite(options_.nonlinear_relative_tolerance)||!(options_.nonlinear_relative_tolerance>0)||!std::isfinite(options_.nonlinear_absolute_tolerance)||!(options_.nonlinear_absolute_tolerance>0)||!std::isfinite(options_.flow_controller_relative_tolerance)||!(options_.flow_controller_relative_tolerance>0)||!std::isfinite(options_.flow_controller_absolute_tolerance_m3_s)||!(options_.flow_controller_absolute_tolerance_m3_s>0)||!std::isfinite(options_.flow_controller_reference_flow_m3_s)||!(options_.flow_controller_reference_flow_m3_s>0)||!std::isfinite(options_.minimum_damping)||!(options_.minimum_damping>0)||options_.minimum_damping>1||!std::isfinite(options_.lu_pivot_shift)||options_.lu_pivot_shift<0||!std::isfinite(options_.nonlinear_block_reduction)||options_.nonlinear_block_reduction<0||options_.reused_preconditioner_maximum_iterations<0||!std::isfinite(options_.reused_preconditioner_maximum_iteration_factor)||options_.reused_preconditioner_maximum_iteration_factor<1.0||!std::isfinite(options_.reused_preconditioner_true_linear_relative_tolerance)||!(options_.reused_preconditioner_true_linear_relative_tolerance>0)) throw std::invalid_argument("immersed transient options are invalid"); ValidateImmersedNitscheWallLabels(options_.wall_labels); }
	void ValidateGeometry() const { if(!surface_.Usable()||!ghost_.Usable()) throw std::invalid_argument("immersed transient geometry catalogs are unusable"); ghost_.ValidateBinding(domain_,volume_); }
	void ConfigurePorts() { std::vector<int> labels; for(const auto& p:options_.ports){ValidateImmersedFlowPortDefinition(p); labels.push_back(p.boundary_label);} std::sort(labels.begin(),labels.end()); if(std::adjacent_find(labels.begin(),labels.end())!=labels.end()) throw std::invalid_argument("immersed transient port labels must be unique"); for(std::size_t i=0;i<options_.ports.size();++i){for(std::size_t j=0;j<i;++j)if(options_.ports[i].id==options_.ports[j].id)throw std::invalid_argument("immersed transient port ids must be unique"); if(std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),options_.ports[i].boundary_label)) throw std::invalid_argument("immersed transient wall and port labels overlap"); diagnostics_.ports.push_back({options_.ports[i].id,options_.ports[i].boundary_label,options_.ports[i].control_mode,options_.ports[i].value});}
		for (const auto& entry : surface_.Diagnostics().source_area_by_boundary_id) {
			const int label = static_cast<int>(entry.first);
			if (!std::binary_search(options_.wall_labels.begin(), options_.wall_labels.end(), label)
				&& !std::binary_search(labels.begin(), labels.end(), label))
				throw std::invalid_argument("immersed transient surface label is neither a wall nor a port");
		}
		for (const int label : options_.wall_labels) { const auto i=surface_.Diagnostics().source_area_by_boundary_id.find(static_cast<std::uint32_t>(label)); if(i==surface_.Diagnostics().source_area_by_boundary_id.end()||!(i->second>0)||!std::isfinite(i->second)) throw std::invalid_argument("immersed transient wall label has no positive source area"); }
		for (const auto& port : options_.ports) { const auto i=surface_.Diagnostics().source_area_by_boundary_id.find(static_cast<std::uint32_t>(port.boundary_label)); if(i==surface_.Diagnostics().source_area_by_boundary_id.end()||!(i->second>0)||!std::isfinite(i->second)) throw std::invalid_argument("immersed transient port label has no positive source area"); long double retained=0; for(std::uint64_t c=0;c<surface_.Cells().size();++c) for(const auto& point:surface_.UsableRule(domain_,c).Points()) if(point.boundary_id==port.boundary_label) retained+=point.weight; const long double scale=std::max(std::abs(retained),std::abs(static_cast<long double>(i->second))); if(!(retained>0)||!std::isfinite(static_cast<double>(retained))||std::abs(retained-static_cast<long double>(i->second))>1024.0L*std::numeric_limits<double>::epsilon()*std::max(scale,std::numeric_limits<long double>::denorm_min())) throw std::runtime_error("immersed transient port retained area fails the surface catalog audit tolerance"); }
		if(!HasGauge() && std::none_of(options_.ports.begin(),options_.ports.end(),[](const auto&p){return IsImmersedFlowPressureLike(p.control_mode);} )) throw std::invalid_argument("immersed transient gauge is mandatory without a pressure-like port");
		ValidateAllFlowCompatibility();
	}
	void ValidateAllFlowCompatibility() const { if(options_.ports.empty()||std::any_of(options_.ports.begin(),options_.ports.end(),[](const auto&p){return IsImmersedFlowPressureLike(p.control_mode);})) return; long double sum=0,scale=0; for(const auto&p:options_.ports){sum+=p.value;scale+=std::abs(static_cast<long double>(p.value));} const long double tolerance=options_.flow_controller_absolute_tolerance_m3_s+options_.flow_controller_relative_tolerance*std::max(scale,static_cast<long double>(options_.flow_controller_reference_flow_m3_s)); if(std::abs(sum)>tolerance) throw std::invalid_argument("all flow-controlled immersed ports require compatible net outward flow"); }
	void PreflightCatalogs() const
	{
		for (std::uint64_t cell=0;cell<domain_.Cells().size();++cell) {
			const auto classification=domain_.Cells()[cell].classification;
			if (classification!=CellClassification::Inside && classification!=CellClassification::Cut) continue;
			if (!volume_.Cell(cell).usable) throw std::runtime_error("classified immersed transient cell has an unusable volume rule");
			if (volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) volume_.ValidateUsableRule(domain_,cell);
			else volume_.ValidateUsableCompactRule(domain_,cell);
			const auto expected=VolumePointCount(cell);
			// Catalog validation above certifies an empty rule. The positive-rule
			// iterator deliberately rejects empty rules, so do not enter it here.
			if (!expected) {
				if (classification==CellClassification::Inside) throw std::runtime_error("inside immersed transient cell has an empty volume rule");
				continue;
			}
			std::size_t visited=0;
			ForEachUsableVolumePoint(cell,[&](const VolumeQuadraturePoint&) {
				if (visited==std::numeric_limits<std::size_t>::max()) throw std::overflow_error("immersed transient volume point count overflows");
				++visited;
			});
			if (visited!=expected) throw std::logic_error("immersed transient volume point iterator count is invalid");
			if (classification==CellClassification::Cut && volume_.Cell(cell).diagnostics.estimated_reference_volume>0) {
				if (!ghost_.Covered(cell)) throw std::runtime_error("positive immersed transient cut cell "+std::to_string(cell)+" requires ghost coverage");
				surface_.ValidateUsableRule(domain_,cell);
			}
		}
	}
	void BuildVolumeBasisCache()
	{
		diagnostics_.volume_basis_cache_enabled=options_.cache_volume_basis;
		Sha256 key;
		immersed_transient_detail::AppendString(key,"ImmersedTransientVolumeBasisCache/v1");
		immersed_transient_detail::AppendString(key,geometry_.GeometryIdentitySha256());
		key.AppendLittleEndian32(static_cast<std::uint32_t>(volume_.StorageMode()));
		diagnostics_.volume_basis_cache_key_sha256=key.Hex();
		if(!options_.cache_volume_basis) return;
		volume_basis_cache_.resize(domain_.Cells().size());
		for(std::uint64_t cell=0;cell<domain_.Cells().size();++cell) if(Usable(cell)) {
			const auto element=domain_.Background().MaterializeElement(cell);
			auto& entries=volume_basis_cache_[cell];
			const auto count=VolumePointCount(cell); entries.reserve(count);
			ForEachUsableVolumePoint(cell,[&](const VolumeQuadraturePoint& point) {
				entries.push_back(BuildNavierStokesVolumePointCacheEntry(element,point));
			});
			if(entries.size()!=count) throw std::logic_error("immersed transient volume basis cache point count is invalid");
			if(std::numeric_limits<std::size_t>::max()-diagnostics_.volume_basis_cache_build_misses<count)
				throw std::overflow_error("immersed transient volume basis cache build counter overflows");
			diagnostics_.volume_basis_cache_build_misses+=count;
		}
		auto add_bytes=[this](std::size_t value) {
			if(std::numeric_limits<std::size_t>::max()-diagnostics_.volume_basis_cache_bytes<value)
				throw std::overflow_error("immersed transient volume basis cache byte count overflows");
			diagnostics_.volume_basis_cache_bytes+=value;
		};
		add_bytes(volume_basis_cache_.capacity()*sizeof(std::vector<NavierStokesVolumePointCacheEntry>));
		for(const auto& entries:volume_basis_cache_) {
			add_bytes(entries.capacity()*sizeof(NavierStokesVolumePointCacheEntry));
			for(const auto& entry:entries) add_bytes(entry.OwnedPayloadBytes());
		}
	}
	// Force callbacks are deliberately evaluated only here, at BeginTrial.  The
	// stored cell/point order is the canonical volume-rule order and therefore
	// neither retries nor later callback mutation can change an assembly.
	std::vector<std::vector<std::array<double,3>>> FreezeBodyForceCandidate() const { std::vector<std::vector<std::array<double,3>>> result(domain_.Cells().size()); for(std::uint64_t c=0;c<domain_.Cells().size();++c) if(Usable(c)){const auto element=domain_.Background().MaterializeElement(c); auto& forces=result[c]; const auto count=VolumePointCount(c); forces.reserve(count); ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint& point){const auto force=options_.body_force(EvaluateElementGeometry(element,point.parametric).physical); for(double value:force) if(!std::isfinite(value)) throw std::invalid_argument("immersed transient body force is not finite on the active quadrature"); forces.push_back(force);}); if(forces.size()!=count)throw std::logic_error("immersed transient frozen body-force point count is invalid");} return result; }
	std::vector<std::uint64_t> ControllerIds() const { std::vector<std::uint64_t> ids; for(const auto&p:options_.ports)if(p.control_mode==ImmersedFlowPortControlMode::FlowRate) ids.push_back(static_cast<std::uint64_t>(p.boundary_label)); std::sort(ids.begin(),ids.end()); return ids; }
	bool HasGauge() const { return options_.include_pressure_gauge && std::none_of(options_.ports.begin(),options_.ports.end(),[](const auto&p){return IsImmersedFlowPressureLike(p.control_mode);}); }
	void RequireIdle(const char* what) const { if(!diagnostics_.idle||diagnostics_.trial_active||diagnostics_.prepared) throw std::logic_error(std::string("cannot ")+what+" outside idle state"); }
	void RequireTrial(const char* what) const { if(!diagnostics_.trial_active||diagnostics_.prepared||!history_) throw std::logic_error(std::string("cannot ")+what+" without active transient trial"); }
	void ValidateGlobal(const ImmersedGlobalFlowState& x) const { if(!x.Valid()||x.GeometryIdentity()!=layout_.GeometryIdentity()||x.NodeIds()!=layout_.NodeIds()||x.PortIds()!=layout_.PortIds()||x.HasGaugeMultiplier()!=HasGauge()) throw std::invalid_argument("immersed transient committed global state layout mismatch"); }
	static bool IsCanonicalPositiveZero(double value) noexcept { const double zero=0.0; return std::memcmp(&value,&zero,sizeof(double))==0; }
	void ValidateMovingInputs(double target_time_s, std::uint64_t target_index, double dt_s,
		const ImmersedVelocityHistory& history, const ImmersedGlobalFlowState& seed,
		const ImmersedMovingTrialMapIdentity& map) const
	{
		if (!map.Valid() || !history.Valid() || !seed.Valid()) throw std::invalid_argument("immersed moving trial input identity is invalid");
		if (seed.GeometryIdentity()!=layout_.GeometryIdentity() || seed.NodeIds()!=layout_.NodeIds()
			|| seed.PortIds()!=layout_.PortIds() || seed.HasGaugeMultiplier()!=HasGauge())
			throw std::invalid_argument("immersed moving trial seed does not match target layout");
		if (seed.Index()==std::numeric_limits<std::uint64_t>::max() || target_index!=seed.Index()+1
			|| seed.TimeS()!=history.SourceTimeS() || target_time_s!=CheckedTransientTargetTime(seed.TimeS(),dt_s)
			|| target_time_s!=geometry_.Evaluation().EvaluatedTimeS())
			throw std::invalid_argument("immersed moving trial seed time/index is inconsistent");
		if (!IsCanonicalPositiveZero(seed.GaugeMultiplier()))
			throw std::invalid_argument("immersed moving trial seed gauge must be canonical positive zero");
		if (history.TargetGeometryIdentity()!=layout_.GeometryIdentity() || history.TargetTimeS()!=target_time_s
			|| history.NodeIds()!=layout_.NodeIds()) throw std::invalid_argument("immersed moving trial history does not exactly cover target geometry/layout/time");
		if (map.OldStateHashSha256().empty() || map.OldGeometryIdentity()!=history.SourceGeometryIdentity()
			|| map.NewGeometryIdentity()!=geometry_.GeometryIdentitySha256()
			|| map.TargetPublicationIdentity()!=geometry_.PublicationIdentitySha256()
			|| map.NewLayoutHashSha256()!=layout_.HashSha256() || map.VelocityHistoryHashSha256()!=history.HashSha256()
			|| map.SourceTimeS()!=seed.TimeS() || map.SourceIndex()!=seed.Index()
			|| map.TargetTimeS()!=target_time_s || map.TargetIndex()!=target_index || map.DtS()!=dt_s
			|| !map.ResetsGaugeToCanonicalZero()) throw std::invalid_argument("immersed moving trial map identity does not bind supplied inputs");
		if (map.ControllerIds()!=layout_.PortIds() || map.ControllerValues().size()!=seed.PortMultipliers().size())
			throw std::invalid_argument("immersed moving trial map controller binding is invalid");
		for (std::size_t i=0;i<map.ControllerValues().size();++i)
			if (std::memcmp(&map.ControllerValues()[i],&seed.PortMultipliers()[i],sizeof(double))!=0)
				throw std::invalid_argument("immersed moving trial seed controller values are not map-bound");
	}
	std::size_t VolumePointCount(std::uint64_t cell) const { if(volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) return volume_.Cell(cell).rule.Points().size(); return CompactCutCellVolumeLogicalPointCount(volume_.Cell(cell).compact_rule); }
	template<class Callback> void ForEachUsableVolumePoint(std::uint64_t cell,Callback&& callback) const { if(volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) for(const auto& point:volume_.UsableRule(domain_,cell).Points()) callback(point); else ForEachVolumePoint(volume_.UsableCompactRule(domain_,cell),std::forward<Callback>(callback)); }
	bool Usable(std::uint64_t cell) const { const auto& q=volume_.Cell(cell); const auto c=domain_.Cells()[cell].classification; return q.usable&&(c==CellClassification::Inside||c==CellClassification::Cut)&&VolumePointCount(cell)!=0; }
	static bool RuleHasLabel(const SurfaceQuadratureRule&r,int label){return std::any_of(r.Points().begin(),r.Points().end(),[&](const auto&p){return p.boundary_id==label;});}
	bool HasWallPoint(const SurfaceQuadratureRule&r) const {return std::any_of(r.Points().begin(),r.Points().end(),[this](const auto&p){return std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),p.boundary_id);});}
	std::array<double,3> MaterialVelocityAt(const ImmersedSurfaceQuadraturePointProvenance& provenance) const
	{
		auto value=geometry_.Evaluation().WallVelocity(provenance.canonical_triangle,provenance.canonical_barycentric);
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if (InjectNonfiniteMaterialVelocityForTesting()) value[0]=std::numeric_limits<double>::quiet_NaN();
#endif
		return value;
	}
	ImmersedMaterialWallVelocityEvaluator MaterialVelocity() const { return [this](const SurfaceQuadraturePoint&,const ImmersedSurfaceQuadraturePointProvenance& p){return MaterialVelocityAt(p);}; }
	void CertifyZeroMaterialWallVelocity() const { for(std::uint64_t c=0;c<surface_.Cells().size();++c){const auto&r=surface_.UsableRule(domain_,c);const auto&v=surface_.UsableProvenance(domain_,c);for(std::size_t i=0;i<r.Points().size();++i)if(std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),r.Points()[i].boundary_id)){const auto w=geometry_.Evaluation().WallVelocity(v[i].canonical_triangle,v[i].canonical_barycentric);for(double x:w)if(!std::isfinite(x)||x!=0.0)throw std::invalid_argument("immersed transient fixed geometry requires exactly zero material wall velocity");}} }
	void CertifyFiniteMaterialSurfaceVelocity() const { for(std::uint64_t c=0;c<surface_.Cells().size();++c){const auto&r=surface_.UsableRule(domain_,c);const auto&v=surface_.UsableProvenance(domain_,c);if(r.Points().size()!=v.size())throw std::logic_error("immersed transient material surface provenance size is invalid");for(std::size_t i=0;i<r.Points().size();++i){const auto w=MaterialVelocityAt(v[i]);for(double x:w)if(!std::isfinite(x))throw std::invalid_argument("immersed transient moving material surface velocity is not finite");}} }
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	static bool& FirstAssemblyFaultArmedForTesting() noexcept { static bool value=false; return value; }
	static bool& InjectNonfiniteMaterialVelocityForTesting() noexcept { static bool value=false; return value; }
	static bool& InjectFiniteConservationDefectForTesting() noexcept { static bool value=false; return value; }
	static bool& InjectFiniteDiscreteMovingWallContinuityDefectForTesting() noexcept { static bool value=false; return value; }
#endif
	// This transaction intentionally has no PETSc command-line override path:
	// every effective solver setting below is fixed or an explicit option and is
	// included in InputHash.  In particular, no untracked options prefix exists.
	// Store the right-preconditioned Krylov directions in FGMRES. Applying
	// a shifted, ill-conditioned LU only to the final combined direction
	// can lose accuracy even when GMRES reports a tiny recursive residual.
	void PreallocateJacobian()
	{
		// Final assembly discards unused SeqAIJ slots. Seed the complete stencil
		// graph before the scalar-diagonal audit, rather than compressing an
		// otherwise empty matrix and reallocating it during every first assembly.
		Mat pattern=nullptr;
		try {
			const PetscInt n=static_cast<PetscInt>(layout_.Rows());
			Check(MatCreate(PETSC_COMM_SELF,&pattern),"MatCreate preallocator");
			Check(MatSetType(pattern,MATPREALLOCATOR),"MatSetType preallocator");
			Check(MatSetSizes(pattern,n,n,n,n),"MatSetSizes preallocator");
			Check(MatSetUp(pattern),"MatSetUp preallocator");
			for(std::uint64_t cell=0;cell<domain_.Cells().size();++cell) if(Usable(cell)) {
				const auto element=domain_.Background().MaterializeElement(cell);
				std::vector<PetscInt> rows;
				for(auto node:element.connectivity) for(int field=0;field<4;++field) rows.push_back(Dof(node,field));
				Check(MatSetValues(pattern,rows.size(),rows.data(),rows.size(),rows.data(),nullptr,INSERT_VALUES),"MatSetValues volume pattern");
				if(domain_.Cells()[cell].classification==CellClassification::Cut) {
					const auto& rule=surface_.UsableRule(domain_,cell);
					for(const auto& port:diagnostics_.ports) if(port.multiplier_row>=0 && RuleHasLabel(rule,port.boundary_label)) {
						std::vector<PetscInt> velocity;
						for(auto node:element.connectivity) for(int field=0;field<3;++field) velocity.push_back(Dof(node,field));
						const PetscInt scalar=port.multiplier_row;
						Check(MatSetValues(pattern,1,&scalar,velocity.size(),velocity.data(),nullptr,INSERT_VALUES),"MatSetValues controller row pattern");
						Check(MatSetValues(pattern,velocity.size(),velocity.data(),1,&scalar,nullptr,INSERT_VALUES),"MatSetValues controller column pattern");
					}
				}
			}
			for(const auto& face:ghost_.Faces()) {
				const auto nodes=CubicCartesianSplineFaceConnectivity(domain_,face.minus_cell,face.plus_cell);
				for(int field=0;field<4;++field) {
					std::vector<PetscInt> rows;
					for(auto node:nodes) rows.push_back(Dof(node,field));
					Check(MatSetValues(pattern,rows.size(),rows.data(),rows.size(),rows.data(),nullptr,INSERT_VALUES),"MatSetValues ghost pattern");
				}
			}
			if(HasGauge()) {
				const PetscInt scalar=GaugeDof();std::vector<PetscInt> pressure;
				for(auto node:layout_.NodeIds()) pressure.push_back(Dof(node,3));
				Check(MatSetValues(pattern,1,&scalar,pressure.size(),pressure.data(),nullptr,INSERT_VALUES),"MatSetValues gauge row pattern");
				Check(MatSetValues(pattern,pressure.size(),pressure.data(),1,&scalar,nullptr,INSERT_VALUES),"MatSetValues gauge column pattern");
			}
			for(const auto scalar:AppendedScalarRows()) Check(MatSetValue(pattern,scalar,scalar,0.,INSERT_VALUES),"MatSetValue scalar pattern");
			Check(MatAssemblyBegin(pattern,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin preallocator");
			Check(MatAssemblyEnd(pattern,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd preallocator");
			Check(MatPreallocatorPreallocate(pattern,PETSC_TRUE,jacobian_),"MatPreallocatorPreallocate");
			Check(MatDestroy(&pattern),"MatDestroy preallocator");
		} catch(...) { if(pattern) MatDestroy(&pattern);throw; }
	}
	void CreatePetsc() { const PetscInt n=static_cast<PetscInt>(layout_.Rows()); try { Check(MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,0,nullptr,&jacobian_),"MatCreateSeqAIJ"); PreallocateJacobian(); Check(MatSetOption(jacobian_,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE),"MatSetOption"); Check(MatSetOption(jacobian_,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE),"MatSetOption retain scalar diagonals"); for(const auto row:AppendedScalarRows()) Check(MatSetValue(jacobian_,row,row,0.0,INSERT_VALUES),"MatSetValue scalar structural diagonal"); Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin scalar structural diagonals");Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd scalar structural diagonals"); AuditAppendedScalarDiagonals(); Check(MatSetOption(jacobian_,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE),"MatSetOption ignore zero entries"); Check(VecCreateSeq(PETSC_COMM_SELF,n,&state_),"VecCreateSeq"); Check(VecDuplicate(state_,&committed_),"VecDuplicate");Check(VecDuplicate(state_,&prepared_),"VecDuplicate");Check(VecDuplicate(state_,&base_),"VecDuplicate");Check(VecDuplicate(state_,&rhs_),"VecDuplicate");Check(VecDuplicate(state_,&update_),"VecDuplicate");Check(VecDuplicate(state_,&action_input_),"VecDuplicate");Check(VecDuplicate(state_,&action_output_),"VecDuplicate"); Check(KSPCreate(PETSC_COMM_SELF,&ksp_),"KSPCreate"); Check(KSPSetType(ksp_,KSPFGMRES),"KSPSetType"); Check(KSPGMRESSetRestart(ksp_,30),"KSPGMRESSetRestart"); Check(KSPSetPCSide(ksp_,PC_RIGHT),"KSPSetPCSide"); Check(KSPSetNormType(ksp_,KSP_NORM_UNPRECONDITIONED),"KSPSetNormType"); Check(KSPSetTolerances(ksp_,options_.ksp_relative_tolerance,1e-50,1e5,options_.ksp_maximum_iterations),"KSPSetTolerances"); PC pc=nullptr;Check(KSPGetPC(ksp_,&pc),"KSPGetPC");Check(PCSetType(pc,PCLU),"PCSetType");if(options_.lu_pivot_shift>0){Check(PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO),"PCFactorSetShiftType");Check(PCFactorSetShiftAmount(pc,options_.lu_pivot_shift),"PCFactorSetShiftAmount");} }catch(...){Destroy();throw;} }
	std::vector<PetscInt> AppendedScalarRows() const {std::vector<PetscInt> rows;for(const auto& p:diagnostics_.ports)if(p.multiplier_row>=0)rows.push_back(p.multiplier_row);if(HasGauge())rows.push_back(GaugeDof());return rows;}
	void AuditAppendedScalarDiagonals(){for(const auto row:AppendedScalarRows()){PetscInt count=0;const PetscInt* columns=nullptr;const PetscScalar* values=nullptr;Check(MatGetRow(jacobian_,row,&count,&columns,&values),"MatGetRow scalar structural diagonal");bool found=false,zero=false;for(PetscInt i=0;i<count;++i)if(columns[i]==row){found=true;zero=PetscRealPart(values[i])==0.0;break;}Check(MatRestoreRow(jacobian_,row,&count,&columns,&values),"MatRestoreRow scalar structural diagonal");if(!found||!zero)throw std::runtime_error("immersed transient scalar row lacks an exact-zero structural diagonal");}diagnostics_.scalar_diagonal_structure_verified=true;}
	void Destroy() noexcept {if(ksp_)KSPDestroy(&ksp_);if(action_output_)VecDestroy(&action_output_);if(action_input_)VecDestroy(&action_input_);if(update_)VecDestroy(&update_);if(rhs_)VecDestroy(&rhs_);if(base_)VecDestroy(&base_);if(prepared_)VecDestroy(&prepared_);if(committed_)VecDestroy(&committed_);if(state_)VecDestroy(&state_);if(jacobian_)MatDestroy(&jacobian_);}
	static double Value(Vec v,PetscInt i){const PetscScalar*x=nullptr;Check(VecGetArrayRead(v,&x),"VecGetArrayRead");const double r=PetscRealPart(x[i]);Check(VecRestoreArrayRead(v,&x),"VecRestoreArrayRead");return r;}
	static void SetVector(Vec v,const std::vector<PetscScalar>&x){PetscScalar*y=nullptr;Check(VecGetArray(v,&y),"VecGetArray");for(std::size_t i=0;i<x.size();++i)y[i]=x[i];Check(VecRestoreArray(v,&y),"VecRestoreArray");}
	std::vector<PetscScalar> Copy(Vec v) const {PetscInt n=0;Check(VecGetSize(v,&n),"VecGetSize");std::vector<PetscScalar>r(static_cast<std::size_t>(n));const PetscScalar*x=nullptr;Check(VecGetArrayRead(v,&x),"VecGetArrayRead");std::copy(x,x+n,r.begin());Check(VecRestoreArrayRead(v,&x),"VecRestoreArrayRead");return r;}
	void SetVectorFromGlobal(const ImmersedGlobalFlowState& x,Vec v){std::vector<PetscScalar>a(layout_.Rows(),0);for(std::size_t i=0;i<x.Coefficients().size();++i)for(int f=0;f<4;++f)a[4*i+f]=x.Coefficients()[i][f];for(std::size_t i=0;i<x.PortMultipliers().size();++i)a[layout_.ControllerRow(x.PortIds()[i])]=x.PortMultipliers()[i];if(HasGauge())a[layout_.GaugeRow()]=x.GaugeMultiplier();SetVector(v,a);}
	std::vector<std::array<double,4>> ReadCoefficients(Vec v) const {auto a=Copy(v);std::vector<std::array<double,4>>r(layout_.NodeIds().size());for(std::size_t i=0;i<r.size();++i)for(int q=0;q<4;++q)r[i][q]=PetscRealPart(a[4*i+q]);return r;}
	std::vector<double> ReadMultipliers(Vec v) const {auto a=Copy(v);std::vector<double>r;for(auto id:layout_.PortIds())r.push_back(PetscRealPart(a[layout_.ControllerRow(id)]));return r;}
	std::vector<double> InitialPortMultipliers() const {return std::vector<double>(layout_.PortIds().size(),0.0);}
	std::vector<std::array<double,4>> Gather(const Element&e)const{auto a=Copy(state_);std::vector<std::array<double,4>>r(e.connectivity.size());for(std::size_t i=0;i<r.size();++i)for(int q=0;q<4;++q)r[i][q]=PetscRealPart(a[Dof(e.connectivity[i],q)]);return r;}
	NavierStokesSystem BuildVolume(const Element&e,const std::vector<std::array<double,4>>&n,std::uint64_t cell,NavierStokesAssemblyRequest request,std::size_t& cache_hits,std::size_t& cache_misses)const{
		cache_hits=cache_misses=0;
		const auto& forces=frozen_body_force_.at(cell); std::size_t point=0;
		const NavierStokesBodyForceEvaluator frozen=[&forces,&point](const std::array<double,3>&){if(point>=forces.size())throw std::logic_error("immersed transient frozen body-force point count is invalid");return forces[point++];};
		ValidateTransientNavierStokesPreflight(options_.parameters);
		if (!moving_trial_) {
			if(!history_||!history_->Valid()||!layout_.Valid()||history_->SourceGeometryIdentity()!=layout_.GeometryIdentity()||history_->TargetGeometryIdentity()!=layout_.GeometryIdentity()) throw std::invalid_argument("immersed velocity history does not match fixed geometry layout");
			for(const auto provenance:history_->Provenance()) if(provenance!=ImmersedVelocityHistoryProvenance::Committed) throw std::invalid_argument("fixed-geometry immersed velocity history must be committed");
			const double expected_target=CheckedTransientTargetTime(history_->SourceTimeS(),options_.parameters.dt);
			if(history_->TargetTimeS()!=expected_target||target_time_s_!=expected_target) throw std::invalid_argument("immersed velocity history target time does not match assembly time");
		} else {
			if (!moving_map_ || !moving_map_->Valid() || !history_ || !history_->Valid() || !layout_.Valid()
				|| history_->TargetGeometryIdentity()!=layout_.GeometryIdentity() || history_->TargetTimeS()!=target_time_s_
				|| history_->NodeIds()!=layout_.NodeIds() || moving_map_->VelocityHistoryHashSha256()!=history_->HashSha256()
				|| moving_map_->NewGeometryIdentity()!=geometry_.GeometryIdentitySha256()
				|| moving_map_->TargetPublicationIdentity()!=geometry_.PublicationIdentitySha256()
				|| moving_map_->NewLayoutHashSha256()!=layout_.HashSha256()
				|| target_time_s_!=CheckedTransientTargetTime(history_->SourceTimeS(),options_.parameters.dt))
				throw std::invalid_argument("immersed moving velocity history does not match target geometry/layout/time");
		}
		const auto local=LocalizeImmersedVelocityHistory(e,layout_,*history_,target_time_s_);
		std::vector<std::array<double,4>> old(local.size()); for(std::size_t i=0;i<old.size();++i) for(int q=0;q<3;++q) old[i][q]=local[i][q];
		NavierStokesSystem result;
		if(options_.cache_volume_basis) {
			const auto& entries=volume_basis_cache_.at(cell);
			if(entries.size()!=forces.size()) throw std::logic_error("immersed transient volume basis cache does not match frozen force order");
			cache_hits=entries.size();
			result=BuildNavierStokesElementFromPreparedPoints(e,n,old,options_.parameters,[&entries](const auto& consume){for(const auto& entry:entries)consume(entry.point,entry.basis,entry.physical);},frozen,NavierStokesResolvedMixedForm::Conservative,request);
		} else {
			cache_misses=VolumePointCount(cell);
			result=BuildNavierStokesElementFromPoints(e,n,old,options_.parameters,[this,cell](const auto& consume){ForEachUsableVolumePoint(cell,consume);},frozen,NavierStokesResolvedMixedForm::Conservative,request);
		}
		if(point!=forces.size()) throw std::logic_error("immersed transient frozen body-force point order is invalid");
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if(volume_probe_for_testing_) volume_probe_for_testing_(cell);
#endif
		return result;
	}
	void Scatter(const std::vector<std::int32_t>&nodes,const NavierStokesSystem&s){std::vector<PetscInt>r;for(auto n:nodes)for(int q=0;q<4;++q)r.push_back(Dof(n,q));if(!s.jacobian.empty())Check(MatSetValues(jacobian_,r.size(),r.data(),r.size(),r.data(),s.jacobian.data(),ADD_VALUES),"MatSetValues");Check(VecSetValues(rhs_,r.size(),r.data(),s.negative_residual.data(),ADD_VALUES),"VecSetValues");}
	void ScatterBlock(const std::vector<std::int32_t>&nodes,const std::vector<PetscScalar>&a,const std::vector<PetscScalar>&b){NavierStokesSystem s{a,b};Scatter(nodes,s);}
	void SubtractAndScatter(const std::vector<std::int32_t>&nodes,NavierStokesSystem a,const NavierStokesSystem&b){for(std::size_t i=0;i<a.jacobian.size();++i)a.jacobian[i]-=b.jacobian[i];for(std::size_t i=0;i<a.negative_residual.size();++i)a.negative_residual[i]-=b.negative_residual[i];Scatter(nodes,a);}
	void ScatterPort(const Element&e,const std::vector<std::array<double,4>>&n,const SurfaceQuadratureRule&r,std::size_t i,bool full){const auto& p=options_.ports[i];const auto x=BuildImmersedFlowPortElement(e,r,p.boundary_label,p.control_mode,p.value,n);std::vector<PetscInt> rows;for(auto id:e.connectivity)for(int q=0;q<4;++q)rows.push_back(Dof(id,q));if(IsImmersedFlowPressureLike(p.control_mode)){Check(VecSetValues(rhs_,rows.size(),rows.data(),x.negative_residual.data(),ADD_VALUES),"VecSetValues port");return;}const auto scalar=static_cast<PetscInt>(layout_.ControllerRow(static_cast<std::uint64_t>(p.boundary_label)));const double lambda=Value(state_,scalar);for(std::size_t k=0;k<rows.size();++k){const double c=PetscRealPart(x.flow_coefficient[k]);if(c){if(full){Check(MatSetValue(jacobian_,rows[k],scalar,c,ADD_VALUES),"MatSetValue port col");Check(MatSetValue(jacobian_,scalar,rows[k],c,ADD_VALUES),"MatSetValue port row");}Check(VecSetValue(rhs_,rows[k],-lambda*c,ADD_VALUES),"VecSetValue port");}}}
	void MeasurePorts(){for(std::size_t p=0;p<options_.ports.size();++p){ImmersedFlowPortMeasurement t{};for(std::uint64_t c=0;c<surface_.Cells().size();++c){const auto&r=surface_.UsableRule(domain_,c);if(!RuleHasLabel(r,options_.ports[p].boundary_label))continue;const auto m=MeasureImmersedFlowPortElement(domain_.Background().MaterializeElement(c),r,options_.ports[p].boundary_label,Gather(domain_.Background().MaterializeElement(c)),options_.parameters.dynamic_viscosity);t.area_m2+=m.area_m2;t.outward_flow_m3_s+=m.outward_flow_m3_s;t.mean_pressure_pa+=m.mean_pressure_pa*m.area_m2;t.mean_normal_traction_pa+=m.mean_normal_traction_pa*m.area_m2;t.mean_velocity_squared_m2_s2+=m.mean_velocity_squared_m2_s2*m.area_m2;}if(!(t.area_m2>0))throw std::runtime_error("immersed transient port area is zero");t.mean_pressure_pa/=t.area_m2;t.mean_normal_traction_pa/=t.area_m2;t.mean_velocity_squared_m2_s2/=t.area_m2;auto&d=trial_ports_.at(p);d.measurement=t;d.measurement_valid=true;if(options_.ports[p].control_mode==ImmersedFlowPortControlMode::FlowRate){d.controller_error=options_.ports[p].value-t.outward_flow_m3_s;d.multiplier=Value(state_,d.multiplier_row);Check(VecSetValue(rhs_,d.multiplier_row,d.controller_error,ADD_VALUES),"VecSetValue flow constraint");}}}
	bool ControllersSatisfied()const{for(std::size_t i=0;i<options_.ports.size();++i)if(options_.ports[i].control_mode==ImmersedFlowPortControlMode::FlowRate){const auto&d=trial_ports_.at(i);const double tol=options_.flow_controller_absolute_tolerance_m3_s+options_.flow_controller_relative_tolerance*std::max(std::abs(options_.ports[i].value),options_.flow_controller_reference_flow_m3_s);if(std::abs(d.controller_error)>tol)return false;}return true;}
	void BuildGaugeWeights(){gauge_weights_.assign(layout_.NodeIds().size(),0);for(std::uint64_t c=0;c<domain_.Cells().size();++c)if(Usable(c)){const auto e=domain_.Background().MaterializeElement(c);ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint&p){const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false);const double m=p.weight*b.raw_determinant;diagnostics_.pressure_measure+=m;for(std::size_t a=0;a<e.connectivity.size();++a)gauge_weights_[layout_.LocalNode(e.connectivity[a])]+=b.value[a]*m;});}}
	void SetHistoryDiagnostics(ImmersedTransientFlowDiagnostics& diagnostics,const ImmersedVelocityHistory* history) const {diagnostics.identity_history_nodes=history?history->NodeIds().size():0;diagnostics.committed_history_nodes=0;diagnostics.extended_history_nodes=0;diagnostics.missing_history_nodes=layout_.NodeIds().size();diagnostics.history_hash_sha256.clear();if(!history)return;diagnostics.history_hash_sha256=history->HashSha256();for(const auto p:history->Provenance())if(p==ImmersedVelocityHistoryProvenance::Committed)++diagnostics.committed_history_nodes;else if(p==ImmersedVelocityHistoryProvenance::Extended)++diagnostics.extended_history_nodes;diagnostics.missing_history_nodes=layout_.NodeIds().size()-diagnostics.identity_history_nodes;}
	void RefreshHistoryDiagnostics(){SetHistoryDiagnostics(diagnostics_,history_.get());}
	void InsertGauge(bool full){const PetscInt q=GaugeDof();double g=0;for(std::size_t i=0;i<gauge_weights_.size();++i){const double w=gauge_weights_[i];g+=w*Value(state_,static_cast<PetscInt>(4*i+3));if(full){Check(MatSetValue(jacobian_,static_cast<PetscInt>(4*i+3),q,w,ADD_VALUES),"MatSetValue gauge");Check(MatSetValue(jacobian_,q,static_cast<PetscInt>(4*i+3),w,ADD_VALUES),"MatSetValue gauge");}Check(VecSetValue(rhs_,static_cast<PetscInt>(4*i+3),-w*Value(state_,q),ADD_VALUES),"VecSetValue gauge");}Check(VecSetValue(rhs_,q,-g,ADD_VALUES),"VecSetValue gauge");diagnostics_.pressure_gauge_defect=std::abs(g);}
	double MaxSurfaceRelativeVelocityNorm() const
	{
		double maximum=0.0;
		for(std::uint64_t c=0;c<domain_.Cells().size();++c) {
			if(!Usable(c) || domain_.Cells()[c].classification!=CellClassification::Cut) continue;
			const auto e=domain_.Background().MaterializeElement(c); const auto n=Gather(e);
			const auto& rule=surface_.UsableRule(domain_,c); const auto& provenance=surface_.UsableProvenance(domain_,c);
			if(rule.Points().size()!=provenance.size()) throw std::logic_error("immersed transient relative-velocity surface provenance size is invalid");
			for(std::size_t i=0;i<rule.Points().size();++i) {
				const auto& p=rule.Points()[i]; const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false);
				const auto w=MaterialVelocityAt(provenance[i]); double squared=0.0;
				for(int q=0;q<3;++q) { if(!std::isfinite(w[q])) throw std::runtime_error("immersed transient relative-velocity material surface velocity is not finite"); double u=0.0; for(std::size_t a=0;a<n.size();++a)u+=n[a][q]*b.value[a]; if(!std::isfinite(u))throw std::runtime_error("immersed transient relative-velocity fluid value is not finite"); const double difference=u-w[q]; squared+=difference*difference; }
				if(!std::isfinite(squared)) throw std::runtime_error("immersed transient relative-velocity norm is not finite");
				maximum=std::max(maximum,std::sqrt(squared));
			}
		}
		return maximum;
	}
	ImmersedTransientFlowConservationDiagnostics MeasureConservation()const
	{
		ImmersedTransientFlowConservationDiagnostics d;
		auto add=[](double& total,double value,const char* what) {
			if(!std::isfinite(value) || !std::isfinite(total) || !std::isfinite(total+value))
				throw std::runtime_error(std::string("immersed transient conservation ")+what+" is nonfinite");
			total+=value;
		};
		auto add_raw_flux=[&](double value) {
			if(d.surface_flux_term_count==std::numeric_limits<std::uint64_t>::max())
				throw std::overflow_error("immersed transient conservation flux count overflows");
			add(d.absolute_surface_flux_sum_m3_s,std::abs(value),"absolute surface flux");
			++d.surface_flux_term_count;
		};
		for(std::uint64_t c=0;c<domain_.Cells().size();++c) if(Usable(c)) {
			const auto e=domain_.Background().MaterializeElement(c); const auto n=Gather(e);
			ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint& p) {
				const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false); double div=0;
				for(std::size_t a=0;a<n.size();++a) for(int q=0;q<3;++q) div+=n[a][q]*b.gradient[a][q];
				add(d.endpoint_volume_divergence_m3_s,p.weight*b.raw_determinant*div,"volume divergence");
			});
			if(domain_.Cells()[c].classification!=CellClassification::Cut) continue;
			const auto& rule=surface_.UsableRule(domain_,c); const auto& provenance=surface_.UsableProvenance(domain_,c);
			if(rule.Points().size()!=provenance.size()) throw std::logic_error("immersed transient conservation surface provenance size is invalid");
			for(std::size_t i=0;i<rule.Points().size();++i) {
				const auto& p=rule.Points()[i]; const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false); double fluid=0;
				for(std::size_t a=0;a<n.size();++a) for(int q=0;q<3;++q) fluid+=n[a][q]*b.value[a]*p.normal[q]*p.weight;
				add_raw_flux(fluid);
				add(d.surface_flow_by_boundary_label_m3_s[p.boundary_id],fluid,"fluid boundary flux"); add(d.total_surface_outward_flow_m3_s,fluid,"total fluid boundary flux");
				const auto w=MaterialVelocityAt(provenance[i]); double material=0;
				for(int q=0;q<3;++q) { if(!std::isfinite(w[q])) throw std::runtime_error("immersed transient conservation material surface velocity is not finite"); material+=w[q]*p.normal[q]*p.weight; }
				add(d.material_surface_outward_flow_by_boundary_label_m3_s[p.boundary_id],material,"material surface flux"); add(d.total_material_surface_outward_flow_m3_s,material,"total material surface flux");
				if(std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),p.boundary_id)) {
					add_raw_flux(material);
					add(d.wall_outward_flow_m3_s,fluid,"wall fluid flux");
					add(d.material_wall_outward_flow_by_boundary_label_m3_s[p.boundary_id],material,"material wall flux"); add(d.total_material_wall_outward_flow_m3_s,material,"total material wall flux");
				} else if(std::any_of(options_.ports.begin(),options_.ports.end(),[&](const auto& port){return port.boundary_label==p.boundary_id;})) add(d.open_port_outward_flow_m3_s,fluid,"open port flux");
				else throw std::runtime_error("immersed transient conservation found an unconfigured surface label");
			}
		}
		d.wall_relative_leakage_m3_s=d.wall_outward_flow_m3_s-d.total_material_wall_outward_flow_m3_s;
		d.divergence_theorem_defect_m3_s=d.endpoint_volume_divergence_m3_s-d.total_surface_outward_flow_m3_s;
		d.discrete_moving_wall_continuity_defect_m3_s=d.open_port_outward_flow_m3_s+d.total_material_wall_outward_flow_m3_s;
		d.discrete_moving_wall_continuity_normalization_scale_m3_s=std::max({options_.flow_controller_reference_flow_m3_s,std::abs(d.open_port_outward_flow_m3_s),std::abs(d.total_material_wall_outward_flow_m3_s)});
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if(InjectFiniteConservationDefectForTesting()) d.divergence_theorem_defect_m3_s=1.0;
		if(InjectFiniteDiscreteMovingWallContinuityDefectForTesting()) d.discrete_moving_wall_continuity_defect_m3_s=1.0;
#endif
		if(!std::isfinite(d.total_material_surface_outward_flow_m3_s)||!std::isfinite(d.wall_relative_leakage_m3_s)||!std::isfinite(d.divergence_theorem_defect_m3_s)||!std::isfinite(d.discrete_moving_wall_continuity_defect_m3_s)||!std::isfinite(d.discrete_moving_wall_continuity_normalization_scale_m3_s)||!(d.discrete_moving_wall_continuity_normalization_scale_m3_s>0.0)) throw std::runtime_error("immersed transient conservation defect is nonfinite or has invalid scale");
		auto validate_map=[](const std::map<int,double>& values,double total,const char* what) {
			double sum=0.0; for(const auto& value:values) { if(!std::isfinite(value.second)||!std::isfinite(sum+value.second)) throw std::runtime_error(std::string("immersed transient conservation ")+what+" label flux is nonfinite"); sum+=value.second; }
			if(std::abs(sum-total)>128.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(total))*std::max<std::size_t>(1,values.size())) throw std::logic_error(std::string("immersed transient conservation ")+what+" label fluxes do not reconcile with total");
		};
		validate_map(d.surface_flow_by_boundary_label_m3_s,d.total_surface_outward_flow_m3_s,"fluid surface");
		validate_map(d.material_surface_outward_flow_by_boundary_label_m3_s,d.total_material_surface_outward_flow_m3_s,"material surface");
		validate_map(d.material_wall_outward_flow_by_boundary_label_m3_s,d.total_material_wall_outward_flow_m3_s,"material wall");
		for(const auto& port:options_.ports) if(d.surface_flow_by_boundary_label_m3_s.find(port.boundary_label)==d.surface_flow_by_boundary_label_m3_s.end()) throw std::logic_error("immersed transient conservation configured port is absent from surface flux map");
		for(const auto& item:d.material_wall_outward_flow_by_boundary_label_m3_s) if(!std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),item.first)) throw std::logic_error("immersed transient conservation material-wall map contains an unconfigured label");
		if(!immersed_transient_detail::MovingWallContinuityIdentityReconciles(d.open_port_outward_flow_m3_s,
			d.wall_outward_flow_m3_s,d.total_material_wall_outward_flow_m3_s,
			d.discrete_moving_wall_continuity_defect_m3_s,d.wall_relative_leakage_m3_s,d.total_surface_outward_flow_m3_s,
			d.absolute_surface_flux_sum_m3_s,d.surface_flux_term_count))
			throw std::logic_error("immersed transient moving-wall continuity roundoff identity does not reconcile");
		const double legacy_scale=std::max(options_.flow_controller_reference_flow_m3_s,std::abs(d.open_port_outward_flow_m3_s));
		d.normalized_wall_leakage=std::abs(d.wall_relative_leakage_m3_s)/legacy_scale;
		d.normalized_open_balance=std::abs(d.divergence_theorem_defect_m3_s)/legacy_scale;
		d.normalized_discrete_moving_wall_continuity_defect=std::abs(d.discrete_moving_wall_continuity_defect_m3_s)/d.discrete_moving_wall_continuity_normalization_scale_m3_s;
		if(!std::isfinite(d.normalized_wall_leakage)||!std::isfinite(d.normalized_open_balance)||!std::isfinite(d.normalized_discrete_moving_wall_continuity_defect)) throw std::runtime_error("immersed transient normalized conservation quotient is nonfinite");
		return d;
	}
	void RefreshHashes(){diagnostics_.committed_state_hash_sha256=committed_global_->HashSha256();diagnostics_.trial_state_hash_sha256=diagnostics_.trial_active?HashVector(state_):std::string{};if(!history_){diagnostics_.identity_history_nodes=diagnostics_.committed_history_nodes=diagnostics_.extended_history_nodes=0;diagnostics_.missing_history_nodes=layout_.NodeIds().size();diagnostics_.history_hash_sha256.clear();}}
	void RefreshHashesNoexcept()noexcept{try{RefreshHashes();}catch(...){}}
	std::string HashVector(Vec v)const{Sha256 h;immersed_transient_detail::AppendString(h,"ImmersedTransientTrial/v1");for(auto x:Copy(v))h.AppendNormalizedDouble(PetscRealPart(x));return h.Hex();}
	void SyncPortMultipliers(std::vector<ImmersedTransientFlowDiagnostics::Port>& ports,Vec state) const { for(auto& port:ports) if(port.multiplier_row>=0) port.multiplier=Value(state,port.multiplier_row); }
	void InvalidatePortMeasurements(std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) const { for(auto& port:ports) { port.measurement={}; port.measurement_valid=false; port.controller_error=0.0; } }
	// Explicit Assemble() is permitted as pre-solve diagnostic work.  It is not
	// solve-attempt work; reset at solve/rollback boundaries makes retries exact.
	static void ResetAttemptWork(ImmersedTransientFlowDiagnostics& diagnostics) noexcept { diagnostics.converged=false; diagnostics.solved_state_hash_sha256.clear(); diagnostics.attempt_hash_sha256.clear(); diagnostics.residual_norm=0.0; diagnostics.true_linear_relative_residual=0.0; diagnostics.nonlinear_iterations=0; diagnostics.ksp_iterations=0; diagnostics.ksp_reason=KSP_CONVERGED_ITERATING; diagnostics.last_assembly_seconds=0.0; diagnostics.last_linear_solve_seconds=0.0; diagnostics.attempt_assembly_count=diagnostics.attempt_full_assembly_count=diagnostics.attempt_residual_only_count=0; diagnostics.preconditioner_builds=diagnostics.preconditioner_reuse_attempts=diagnostics.preconditioner_reuse_accepts=diagnostics.preconditioner_rebuilds=0; diagnostics.preconditioner_rejections_ksp=diagnostics.preconditioner_rejections_iterations=diagnostics.preconditioner_rejections_true_residual=0; diagnostics.volume_cells=diagnostics.surface_cells=diagnostics.ghost_faces=0; diagnostics.wall_penalty={}; diagnostics.newton_steps.clear(); }
	static void AccumulateWallPenaltyDiagnostics(ImmersedNitscheWallDiagnostics& aggregate,
		const ImmersedNitscheWallDiagnostics& local)
	{
		// Build a complete candidate first.  Assembly already delays publication
		// of its candidate aggregate, and this keeps the helper equally atomic
		// for any future caller that accumulates directly into a published record.
		auto candidate=aggregate;
		AccumulateWallPenaltyDiagnosticsInto(candidate,local);
		using std::swap;
		swap(aggregate,candidate);
	}
	static void AccumulateWallPenaltyDiagnosticsInto(ImmersedNitscheWallDiagnostics& aggregate,
		const ImmersedNitscheWallDiagnostics& local)
	{
		bool local_selected=false, aggregate_selected=false;
		for(const auto& item:local.by_boundary_id) local_selected=local_selected||item.second.selected_points!=0;
		for(const auto& item:aggregate.by_boundary_id) aggregate_selected=aggregate_selected||item.second.selected_points!=0;
		if(!local_selected) return;
		const std::array<double,17> values{{local.minimum_h_n_m,local.maximum_h_n_m,
			local.minimum_eta,local.maximum_eta,
			local.minimum_eta_mu,local.maximum_eta_mu,local.minimum_eta_t,local.maximum_eta_t,
			local.minimum_eta_mu_h_n_over_mu,local.maximum_eta_mu_h_n_over_mu,
			local.minimum_eta_t_dt_over_rho_h_n,local.maximum_eta_t_dt_over_rho_h_n,
			local.minimum_eta_mu_fraction,local.maximum_eta_mu_fraction,
			local.minimum_eta_t_fraction,local.maximum_eta_t_fraction,local.maximum_gap_norm}};
		for(double value:values) if(!std::isfinite(value)||value<0.0)
			throw std::overflow_error("immersed transient wall penalty diagnostic is invalid");
		if(!(local.minimum_h_n_m>0.0)||!(local.maximum_h_n_m>=local.minimum_h_n_m))
			throw std::overflow_error("immersed transient wall h_n diagnostic is invalid");
		const std::array<double,3> local_fractions{{local.fraction_lower,local.fraction_estimate,local.fraction_upper}};
		for(double value:local_fractions) if(!std::isfinite(value)||value<0.0||value>1.0)
			throw std::overflow_error("immersed transient local wall reference fraction is invalid");
		if(!std::isfinite(local.fraction_ordering_tolerance)||local.fraction_ordering_tolerance<0.0
			||local.fraction_ordering_tolerance>2.0e-12)
			throw std::overflow_error("immersed transient local wall fraction tolerance is invalid");
		if(local.fraction_lower>local.fraction_estimate+local.fraction_ordering_tolerance
			||local.fraction_estimate>local.fraction_upper+local.fraction_ordering_tolerance)
			throw std::overflow_error("immersed transient local wall reference fraction ordering is invalid");
		const auto checked_add=[](double& out,double value,const char* description) {
			if(!std::isfinite(out)||!std::isfinite(value)||value<0.0||!std::isfinite(out+value))
				throw std::overflow_error(description);
			out+=value;
		};
		const auto merge_min=[aggregate_selected](double& out,double value) { out=aggregate_selected?std::min(out,value):value; };
		const auto merge_max=[](double& out,double value) { out=std::max(out,value); };
		merge_min(aggregate.minimum_h_n_m,local.minimum_h_n_m); merge_max(aggregate.maximum_h_n_m,local.maximum_h_n_m);
		merge_min(aggregate.minimum_eta,local.minimum_eta); merge_max(aggregate.maximum_eta,local.maximum_eta);
		merge_min(aggregate.minimum_eta_mu,local.minimum_eta_mu); merge_max(aggregate.maximum_eta_mu,local.maximum_eta_mu);
		merge_min(aggregate.minimum_eta_t,local.minimum_eta_t); merge_max(aggregate.maximum_eta_t,local.maximum_eta_t);
		merge_min(aggregate.minimum_eta_mu_h_n_over_mu,local.minimum_eta_mu_h_n_over_mu); merge_max(aggregate.maximum_eta_mu_h_n_over_mu,local.maximum_eta_mu_h_n_over_mu);
		merge_min(aggregate.minimum_eta_t_dt_over_rho_h_n,local.minimum_eta_t_dt_over_rho_h_n); merge_max(aggregate.maximum_eta_t_dt_over_rho_h_n,local.maximum_eta_t_dt_over_rho_h_n);
		merge_min(aggregate.minimum_eta_mu_fraction,local.minimum_eta_mu_fraction); merge_max(aggregate.maximum_eta_mu_fraction,local.maximum_eta_mu_fraction);
		merge_min(aggregate.minimum_eta_t_fraction,local.minimum_eta_t_fraction); merge_max(aggregate.maximum_eta_t_fraction,local.maximum_eta_t_fraction);
		merge_max(aggregate.maximum_eta_h_n_over_mu,local.maximum_eta_h_n_over_mu);
		merge_max(aggregate.maximum_gap_norm,local.maximum_gap_norm);
		checked_add(aggregate.fraction_lower,local.fraction_lower,"immersed transient wall reference fraction lower total overflows");
		checked_add(aggregate.fraction_estimate,local.fraction_estimate,"immersed transient wall reference fraction estimate total overflows");
		checked_add(aggregate.fraction_upper,local.fraction_upper,"immersed transient wall reference fraction upper total overflows");
		checked_add(aggregate.fraction_ordering_tolerance,local.fraction_ordering_tolerance,
			"immersed transient wall fraction tolerance total overflows");
		const double roundoff=8.0*std::numeric_limits<double>::epsilon()*std::max({1.0,
			aggregate.fraction_lower,aggregate.fraction_estimate,aggregate.fraction_upper});
		checked_add(aggregate.fraction_ordering_tolerance,roundoff,"immersed transient wall fraction roundoff overflows");
		if(aggregate.fraction_lower>aggregate.fraction_estimate+aggregate.fraction_ordering_tolerance
			||aggregate.fraction_estimate>aggregate.fraction_upper+aggregate.fraction_ordering_tolerance)
			throw std::logic_error("immersed transient wall reference fraction totals are unordered");
		aggregate.ghost_covered_policy=local.ghost_covered_policy;
		for(const auto& item:local.by_boundary_id) {
			auto& out=aggregate.by_boundary_id[item.first];
			if(std::numeric_limits<std::size_t>::max()-out.selected_points<item.second.selected_points
				||std::numeric_limits<std::size_t>::max()-out.skipped_points<item.second.skipped_points)
				throw std::overflow_error("immersed transient wall point diagnostic overflows");
			out.selected_points+=item.second.selected_points; out.skipped_points+=item.second.skipped_points;
			for(const double value:{item.second.selected_area_m2,item.second.skipped_area_m2}) if(!std::isfinite(value)||value<0.0)
				throw std::overflow_error("immersed transient wall area diagnostic is invalid");
			if(!std::isfinite(out.selected_area_m2+item.second.selected_area_m2)||!std::isfinite(out.skipped_area_m2+item.second.skipped_area_m2))
				throw std::overflow_error("immersed transient wall area diagnostic overflows");
			out.selected_area_m2+=item.second.selected_area_m2; out.skipped_area_m2+=item.second.skipped_area_m2;
		}
	}
	void ResetAttemptWork() noexcept { ResetAttemptWork(diagnostics_); }
	std::uint64_t state_generation_=0,residual_generation_=std::numeric_limits<std::uint64_t>::max(),jacobian_generation_=std::numeric_limits<std::uint64_t>::max();
	void InvalidateSolved() noexcept { ResetAttemptWork(); }
	void ClearIdlePublicationDiagnostics() noexcept { diagnostics_.input_hash_sha256.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.trial_state_hash_sha256.clear(); diagnostics_.history_hash_sha256.clear(); diagnostics_.moving_map_identity_sha256.clear(); diagnostics_.target_time_s=diagnostics_.dt_s=diagnostics_.pressure_gauge_defect=0.0; diagnostics_.prepared=false; ResetAttemptWork(); }
	void MarkSolved() { RefreshHashes(); diagnostics_.solved_state_hash_sha256=diagnostics_.trial_state_hash_sha256; diagnostics_.converged=true; diagnostics_.attempt_hash_sha256=AttemptHash(); }
	static void AppendPorts(Sha256& h,const std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) { h.AppendLittleEndian64(ports.size()); for(const auto&p:ports){immersed_transient_detail::AppendString(h,p.id);h.AppendLittleEndian32(static_cast<std::uint32_t>(p.boundary_label));h.AppendLittleEndian32(static_cast<std::uint32_t>(p.control_mode));h.AppendNormalizedDouble(p.target);h.AppendNormalizedDouble(p.multiplier);h.AppendNormalizedDouble(p.controller_error);h.AppendLittleEndian64(static_cast<std::uint64_t>(p.multiplier_row));h.AppendNormalizedDouble(p.measurement.area_m2);h.AppendNormalizedDouble(p.measurement.outward_flow_m3_s);h.AppendNormalizedDouble(p.measurement.mean_pressure_pa);h.AppendNormalizedDouble(p.measurement.mean_normal_traction_pa);h.AppendNormalizedDouble(p.measurement.mean_velocity_squared_m2_s2);h.AppendLittleEndian32(p.measurement_valid?1:0);}}
	static void AppendFrozenBodyForce(Sha256& h,const std::vector<std::vector<std::array<double,3>>>& forces) { h.AppendLittleEndian64(forces.size()); for(const auto& cell:forces) { h.AppendLittleEndian64(cell.size()); for(const auto& force:cell) for(double value:force) h.AppendNormalizedDouble(value); } }
	void AppendFrozenBodyForce(Sha256& h) const { AppendFrozenBodyForce(h,frozen_body_force_); }
	void AppendFixedSolverConfiguration(Sha256& h) const { immersed_transient_detail::AppendString(h,"KSPFGMRES"); h.AppendLittleEndian64(30); immersed_transient_detail::AppendString(h,"PCLU"); immersed_transient_detail::AppendString(h,"PC_RIGHT"); immersed_transient_detail::AppendString(h,"KSP_NORM_UNPRECONDITIONED"); h.AppendNormalizedDouble(1e-50); h.AppendNormalizedDouble(1e5);
		PetscOptionEntries configured;
		const auto prefix = "-"+CanonicalPetscOptionName(solver_options_->Prefix());
		for (const auto& entry : CapturePetscOptionEntries(solver_options_->Database()))
			if (CanonicalPetscOptionName(entry.first).compare(0, prefix.size(), prefix)==0) configured["-"+CanonicalPetscOptionName(entry.first).substr(prefix.size())] = entry.second;
		immersed_transient_detail::AppendString(h, SerializePetscOptionEntries(configured));
	}
	std::array<double,3> ResidualBlockNorms(const std::vector<PetscScalar>& residual) const { if(residual.size()!=diagnostics_.total_dofs) throw std::invalid_argument("immersed transient residual block vector size is invalid"); std::array<double,3> norms{{0.0,0.0,0.0}}; for(std::size_t row=0;row<residual.size();++row){const double value=PetscRealPart(residual[row]);const std::size_t block=row>=diagnostics_.physical_dofs?2:((row%4)==3?1:0);norms[block]+=value*value;}for(auto& value:norms)value=std::sqrt(value);return norms; }
	bool BlockReductionSatisfied(const std::array<double,3>& initial,const std::vector<PetscScalar>& residual,double global_initial) const { if(options_.nonlinear_block_reduction==0.0)return true; const auto current=ResidualBlockNorms(residual); const double zero_block_tolerance=options_.nonlinear_absolute_tolerance+options_.nonlinear_relative_tolerance*global_initial; for(std::size_t block=0;block<current.size();++block)if(initial[block]>0.0?current[block]>initial[block]/options_.nonlinear_block_reduction:current[block]>zero_block_tolerance)return false;return true; }
	void AppendPreconditionerReuseConfiguration(Sha256& h) const { h.AppendLittleEndian32(options_.reuse_preconditioner?1:0); h.AppendLittleEndian64(static_cast<std::uint64_t>(options_.reused_preconditioner_maximum_iterations)); h.AppendNormalizedDouble(options_.reused_preconditioner_maximum_iteration_factor); h.AppendNormalizedDouble(options_.reused_preconditioner_true_linear_relative_tolerance); }
	std::string InputHash(const NavierStokesParameters& parameters,const ImmersedVelocityHistory& history,double target_time_s,std::uint64_t target_index,const std::vector<ImmersedTransientFlowDiagnostics::Port>& frozen_ports,const std::vector<std::vector<std::array<double,3>>>& frozen_forces,const std::vector<PetscScalar>& frozen_seed) const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientInput/v7"); immersed_transient_detail::AppendString(h,layout_.HashSha256()); immersed_transient_detail::AppendString(h,geometry_.GeometryIdentitySha256()); h.AppendLittleEndian32(static_cast<std::uint32_t>(volume_.StorageMode())); immersed_transient_detail::AppendString(h,history.HashSha256()); h.AppendNormalizedDouble(committed_global_->TimeS());h.AppendLittleEndian64(committed_global_->Index());h.AppendNormalizedDouble(target_time_s);h.AppendLittleEndian64(target_index);h.AppendNormalizedDouble(parameters.density);h.AppendNormalizedDouble(parameters.dynamic_viscosity);h.AppendNormalizedDouble(parameters.dt);h.AppendNormalizedDouble(options_.wall_gamma0);h.AppendNormalizedDouble(options_.wall_inertial_gamma0);h.AppendLittleEndian32(options_.include_pressure_gauge?1:0);h.AppendLittleEndian64(options_.nonlinear_maximum_iterations);h.AppendLittleEndian64(options_.ksp_maximum_iterations);h.AppendNormalizedDouble(options_.ksp_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_absolute_tolerance);h.AppendNormalizedDouble(options_.flow_controller_relative_tolerance);h.AppendNormalizedDouble(options_.flow_controller_absolute_tolerance_m3_s);h.AppendNormalizedDouble(options_.flow_controller_reference_flow_m3_s);h.AppendNormalizedDouble(options_.minimum_damping);h.AppendNormalizedDouble(options_.lu_pivot_shift);h.AppendNormalizedDouble(options_.nonlinear_block_reduction);AppendPreconditionerReuseConfiguration(h);AppendFixedSolverConfiguration(h);AppendPorts(h,frozen_ports);AppendFrozenBodyForce(h,frozen_forces);for(auto x:frozen_seed)h.AppendNormalizedDouble(PetscRealPart(x));return h.Hex(); }
	std::string MovingInputHash(const NavierStokesParameters& parameters,const ImmersedVelocityHistory& history,double target_time_s,std::uint64_t target_index,const std::vector<ImmersedTransientFlowDiagnostics::Port>& frozen_ports,const std::vector<std::vector<std::array<double,3>>>& frozen_forces,const std::vector<PetscScalar>& frozen_seed,const ImmersedMovingTrialMapIdentity& map) const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientMovingInput/v4"); immersed_transient_detail::AppendString(h,map.HashSha256()); immersed_transient_detail::AppendString(h,layout_.HashSha256()); immersed_transient_detail::AppendString(h,geometry_.GeometryIdentitySha256()); immersed_transient_detail::AppendString(h,geometry_.PublicationIdentitySha256()); h.AppendLittleEndian32(static_cast<std::uint32_t>(volume_.StorageMode())); immersed_transient_detail::AppendString(h,history.HashSha256()); h.AppendNormalizedDouble(target_time_s);h.AppendLittleEndian64(target_index);h.AppendNormalizedDouble(parameters.density);h.AppendNormalizedDouble(parameters.dynamic_viscosity);h.AppendNormalizedDouble(parameters.dt);h.AppendNormalizedDouble(options_.wall_gamma0);h.AppendNormalizedDouble(options_.wall_inertial_gamma0);h.AppendLittleEndian32(options_.include_pressure_gauge?1:0);h.AppendLittleEndian64(options_.nonlinear_maximum_iterations);h.AppendLittleEndian64(options_.ksp_maximum_iterations);h.AppendNormalizedDouble(options_.ksp_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_absolute_tolerance);h.AppendNormalizedDouble(options_.flow_controller_relative_tolerance);h.AppendNormalizedDouble(options_.flow_controller_absolute_tolerance_m3_s);h.AppendNormalizedDouble(options_.flow_controller_reference_flow_m3_s);h.AppendNormalizedDouble(options_.minimum_damping);h.AppendNormalizedDouble(options_.lu_pivot_shift);h.AppendNormalizedDouble(options_.nonlinear_block_reduction);AppendPreconditionerReuseConfiguration(h);AppendFixedSolverConfiguration(h);AppendPorts(h,frozen_ports);AppendFrozenBodyForce(h,frozen_forces);for(auto x:frozen_seed)h.AppendNormalizedDouble(PetscRealPart(x));return h.Hex(); }
	std::string InputHash() const { if(!history_) return {}; return InputHash(options_.parameters,*history_,target_time_s_,target_index_,frozen_ports_,frozen_body_force_,frozen_seed_); }
	std::string AttemptHash() const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientAttempt/v4"); immersed_transient_detail::AppendString(h,diagnostics_.input_hash_sha256); immersed_transient_detail::AppendString(h,HashVector(state_)); AppendPorts(h,trial_ports_); h.AppendLittleEndian64(diagnostics_.attempt_assembly_count);h.AppendLittleEndian64(diagnostics_.attempt_full_assembly_count);h.AppendLittleEndian64(diagnostics_.attempt_residual_only_count);h.AppendLittleEndian64(diagnostics_.preconditioner_builds);h.AppendLittleEndian64(diagnostics_.preconditioner_reuse_attempts);h.AppendLittleEndian64(diagnostics_.preconditioner_reuse_accepts);h.AppendLittleEndian64(diagnostics_.preconditioner_rebuilds);h.AppendLittleEndian64(diagnostics_.preconditioner_rejections_ksp);h.AppendLittleEndian64(diagnostics_.preconditioner_rejections_iterations);h.AppendLittleEndian64(diagnostics_.preconditioner_rejections_true_residual);h.AppendLittleEndian64(diagnostics_.newton_steps.size());for(const auto&s:diagnostics_.newton_steps){h.AppendLittleEndian64(s.iteration);h.AppendLittleEndian64(s.ksp_iterations);h.AppendLittleEndian32(static_cast<std::uint32_t>(s.ksp_reason));h.AppendNormalizedDouble(s.residual_norm);h.AppendNormalizedDouble(s.update_norm);h.AppendNormalizedDouble(s.linear_relative_residual);h.AppendNormalizedDouble(s.damping);} h.AppendLittleEndian64(diagnostics_.ksp_iterations);h.AppendLittleEndian64(diagnostics_.nonlinear_iterations);h.AppendLittleEndian32(static_cast<std::uint32_t>(diagnostics_.ksp_reason));h.AppendNormalizedDouble(diagnostics_.residual_norm);h.AppendNormalizedDouble(diagnostics_.true_linear_relative_residual);return h.Hex(); }
	std::string PreparedHash(const ImmersedGlobalFlowState& global,const std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientPrepared/v2"); immersed_transient_detail::AppendString(h,diagnostics_.attempt_hash_sha256); immersed_transient_detail::AppendString(h,global.HashSha256());AppendPorts(h,ports);return h.Hex(); }
	std::shared_ptr<PetscSolverOptions> solver_options_;
	ElementAssemblyExecution assembly_execution_;
	bool assembly_pending_=false;
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	std::function<void(std::uint64_t)> volume_probe_for_testing_;
	std::function<void(PetscInt,double,double&,double)> line_search_norm_probe_for_testing_;
	std::function<void(PetscInt)> newton_boundary_probe_for_testing_;
#endif
	const MovingCutGeometry& geometry_; const CartesianDomainClassification& domain_; const CutCellVolumeQuadratureCatalog& volume_; const ImmersedSurfaceQuadratureCatalog& surface_; const CutCellGhostPenaltyCatalog& ghost_; ImmersedTransientFlowOptions options_; ImmersedActiveLayout layout_; std::vector<double> gauge_weights_; std::vector<std::vector<NavierStokesVolumePointCacheEntry>> volume_basis_cache_; std::unique_ptr<ImmersedGlobalFlowState> committed_global_,prepared_global_; std::unique_ptr<ImmersedVelocityHistory> history_; std::unique_ptr<ImmersedMovingTrialMapIdentity> moving_map_; std::vector<PetscScalar> frozen_seed_; std::vector<std::vector<std::array<double,3>>> frozen_body_force_; std::vector<ImmersedTransientFlowDiagnostics::Port> frozen_ports_,trial_ports_,prepared_ports_; std::string prepared_global_hash_; double target_time_s_=0; std::uint64_t target_index_=0; Mat jacobian_=nullptr;Vec state_=nullptr,committed_=nullptr,prepared_=nullptr,base_=nullptr,rhs_=nullptr,update_=nullptr,action_input_=nullptr,action_output_=nullptr;KSP ksp_=nullptr;ImmersedTransientFlowDiagnostics diagnostics_{};bool moving_trial_=false,fail_next_prepare_=false;
};

} // namespace iga

#endif

#ifndef IGA_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP

// Fixed-geometry backward-Euler immersed flow transaction.  This is purposefully
// separate from ImmersedStaticFlowRuntime: a moving evaluation is an immutable
// input and this PR accepts only an identity transition.
#include "MovingCutGeometry.hpp"
#include "ImmersedFlowPort.hpp"
#include "ImmersedNitscheWall.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
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
	NavierStokesBodyForceEvaluator body_force = [](const std::array<double,3>&) { return std::array<double,3>{{0,0,0}}; };
	bool include_pressure_gauge = true;
	PetscInt nonlinear_maximum_iterations = 12, ksp_maximum_iterations = 2000;
	double ksp_relative_tolerance = 1e-10, nonlinear_relative_tolerance = 1e-9, nonlinear_absolute_tolerance = 1e-11;
	double flow_controller_relative_tolerance = 1e-10, flow_controller_absolute_tolerance_m3_s = 1e-15, flow_controller_reference_flow_m3_s = 1e-12;
	double minimum_damping = 1.0/128.0, lu_pivot_shift = 0.0, nonlinear_block_reduction = 0.0;
};

struct ImmersedTransientFlowNewtonStep {
	PetscInt iteration = 0, ksp_iterations = 0;
	KSPConvergedReason ksp_reason = KSP_CONVERGED_ITERATING;
	double residual_norm = 0.0, update_norm = 0.0, linear_relative_residual = 0.0, damping = 0.0;
};

struct ImmersedTransientFlowConservationDiagnostics {
	std::map<int,double> surface_flow_by_boundary_label_m3_s;
	double endpoint_volume_divergence_m3_s = 0.0, total_surface_outward_flow_m3_s = 0.0;
	double open_port_outward_flow_m3_s = 0.0, wall_outward_flow_m3_s = 0.0, wall_relative_leakage_m3_s = 0.0;
	double normalized_open_balance = 0.0, normalized_wall_leakage = 0.0;
};

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
	std::string geometry_identity_sha256, layout_hash_sha256, committed_state_hash_sha256, trial_state_hash_sha256, history_hash_sha256;
	// Deterministic publication identities deliberately exclude timing and lifetime
	// attempt/abort counters.  They make an exact retry auditable.
	std::string input_hash_sha256, solved_state_hash_sha256, prepared_hash_sha256, attempt_hash_sha256;
	std::size_t attempt_assembly_count = 0;
	// Measurements describe a particular published state.  Idle mutation clears
	// them instead of allowing a previous committed state to leak into an input.
	struct Port { std::string id; int boundary_label = -1; ImmersedFlowPortControlMode control_mode = ImmersedFlowPortControlMode::Pressure; double target = 0.0, multiplier = 0.0, controller_error = 0.0; PetscInt multiplier_row = -1; ImmersedFlowPortMeasurement measurement{}; bool measurement_valid = false; };
	std::vector<Port> ports;
	std::vector<ImmersedTransientFlowNewtonStep> newton_steps;
};

class ImmersedTransientFlowRuntime {
public:
	ImmersedTransientFlowRuntime(const MovingCutGeometry& geometry, ImmersedTransientFlowOptions options = {})
		: geometry_(geometry), domain_(geometry.Domain()), volume_(geometry.Volume()), surface_(geometry.Surface()), ghost_(geometry.Ghost()), options_(std::move(options))
	{
		// A constructor whose body throws does not run this object's destructor.
		// Keep all work following possible PETSc handle creation inside this guard.
		try {
		ValidateOptions(); ValidateGeometry(); ConfigurePorts(); PreflightCatalogs();
		layout_ = ImmersedActiveLayout::Build(domain_, volume_, geometry_.GeometryIdentitySha256(), ControllerIds(), HasGauge());
		for (auto& port : diagnostics_.ports)
			if (port.control_mode == ImmersedFlowPortControlMode::FlowRate)
				port.multiplier_row = static_cast<PetscInt>(layout_.ControllerRow(static_cast<std::uint64_t>(port.boundary_label)));
		if (layout_.NodeIds().empty()) throw std::runtime_error("immersed transient active layout is empty");
		diagnostics_.active_nodes=layout_.NodeIds().size(); diagnostics_.physical_dofs=layout_.NodeFieldRows(); diagnostics_.total_dofs=layout_.Rows();
		diagnostics_.geometry_identity_sha256=geometry_.GeometryIdentitySha256(); diagnostics_.layout_hash_sha256=layout_.HashSha256();
		BuildGaugeWeights(); CreatePetsc();
		std::vector<std::array<double,4>> zero(layout_.NodeIds().size());
		committed_global_.reset(new ImmersedGlobalFlowState(geometry_.Evaluation().EvaluatedTimeS(), 0, layout_, std::move(zero), InitialPortMultipliers(), HasGauge(), 0.0));
		SetVectorFromGlobal(*committed_global_, committed_); Check(VecCopy(committed_, state_), "VecCopy initial state");
		SyncPortMultipliers(diagnostics_.ports, committed_); RefreshHashes();
		} catch (...) { Destroy(); throw; }
	}
	~ImmersedTransientFlowRuntime() { Destroy(); }
	ImmersedTransientFlowRuntime(const ImmersedTransientFlowRuntime&) = delete;
	ImmersedTransientFlowRuntime& operator=(const ImmersedTransientFlowRuntime&) = delete;

	const ImmersedActiveLayout& Layout() const noexcept { return layout_; }
	const ImmersedTransientFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const ImmersedGlobalFlowState& CommittedGlobalState() const { return *committed_global_; }
	std::vector<PetscScalar> CommittedState() const { return Copy(committed_); }
	std::vector<PetscScalar> TrialState() const { return Copy(state_); }
	PetscInt Dof(std::int32_t id, int field) const { if(field<0||field>3) throw std::out_of_range("immersed transient field is invalid"); return static_cast<PetscInt>(4*layout_.LocalNode(id)+field); }
	PetscInt PortMultiplierDof(const std::string& id) const { for(const auto& p:diagnostics_.ports) if(p.id==id) return p.multiplier_row; throw std::out_of_range("immersed transient port id is absent"); }
	PetscInt GaugeDof() const { if(!HasGauge()) throw std::logic_error("immersed transient gauge is absent"); return static_cast<PetscInt>(layout_.GaugeRow()); }
	void FailNextPrepareForTesting() noexcept { fail_next_prepare_ = true; }
	// A trial-only warm start for finite-difference and restart workflows.  It
	// deliberately cannot alter the frozen committed state, identity history, or
	// rollback seed.
	void SetTrialState(const std::vector<PetscScalar>& value)
	{
		RequireTrial("set trial state");
		if(value.size()!=layout_.Rows()) throw std::invalid_argument("immersed transient trial state size is invalid");
		for(const auto x:value) if(!std::isfinite(PetscRealPart(x))) throw std::invalid_argument("immersed transient trial state is not finite");
		SetVector(state_,value); InvalidateSolved(); RefreshHashes();
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
		std::swap(state_,prepared_); target_time_s_=target_time_s; target_index_=target_index;
		std::swap(diagnostics_,candidate_diagnostics);
	}
	void Assemble()
	{
		RequireTrial("assemble"); const auto start=std::chrono::steady_clock::now();
		if(diagnostics_.converged) InvalidateSolved();
		++diagnostics_.attempt_assembly_count;
		Check(MatZeroEntries(jacobian_),"MatZeroEntries"); Check(VecSet(rhs_,0.0),"VecSet rhs"); diagnostics_.volume_cells=diagnostics_.surface_cells=diagnostics_.ghost_faces=0;
		for(std::uint64_t cell=0;cell<domain_.Cells().size();++cell) if(Usable(cell)) {
			const auto element=domain_.Background().MaterializeElement(cell); const auto nodal=Gather(element); const auto volume_system=BuildVolume(element,nodal,cell);
			Scatter(element.connectivity,volume_system); ++diagnostics_.volume_cells;
			if(domain_.Cells()[cell].classification==CellClassification::Cut) {
				const auto& rule=surface_.UsableRule(domain_,cell); Scatter(element.connectivity,BuildImmersedConservativeMixedTraceElement(element,rule,nodal));
				if(HasWallPoint(rule)) { const auto local_history=LocalizeImmersedVelocityHistory(element,layout_,*history_,target_time_s_); std::vector<std::array<double,4>> old(local_history.size()); for(std::size_t i=0;i<old.size();++i) for(int q=0;q<3;++q) old[i][q]=local_history[i][q];
					auto wall=BuildImmersedNitscheWallElementFromVolumeSystemMaterialAware(domain_,volume_,surface_,cell,nodal,old,options_.parameters,options_.wall_labels,volume_system,ghost_,options_.wall_gamma0,MaterialVelocity());
					SubtractAndScatter(element.connectivity,wall.system,volume_system); ++diagnostics_.surface_cells; }
				for(std::size_t p=0;p<options_.ports.size();++p) if(RuleHasLabel(rule,options_.ports[p].boundary_label)) ScatterPort(element,nodal,rule,p);
			}
		}
		for(std::size_t f=0;f<ghost_.Faces().size();++f) { const auto block=ghost_.AssembleFaceLocal(f,domain_,volume_,[this](std::int32_t node,int q){return Value(state_,Dof(node,q));},options_.parameters.dynamic_viscosity); ScatterBlock(block.connectivity,block.jacobian,block.negative_residual); ++diagnostics_.ghost_faces; }
		Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin"); Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd"); Check(VecAssemblyBegin(rhs_),"VecAssemblyBegin"); Check(VecAssemblyEnd(rhs_),"VecAssemblyEnd");
		MeasurePorts(); if(HasGauge()) InsertGauge(); Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin gauge"); Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd gauge"); Check(VecAssemblyBegin(rhs_),"VecAssemblyBegin gauge"); Check(VecAssemblyEnd(rhs_),"VecAssemblyEnd gauge");
		diagnostics_.last_assembly_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); RefreshHashes();
	}
	bool SolveTrial()
	{
		RequireTrial("solve"); ResetAttemptWork(); double initial=-1.0; std::array<double,3> initial_blocks{{0.0,0.0,0.0}};
		try { for(PetscInt it=0;it<options_.nonlinear_maximum_iterations;++it) { Assemble(); PetscReal norm=0; Check(VecNorm(rhs_,NORM_2,&norm),"VecNorm residual"); if(!std::isfinite(norm)) throw std::runtime_error("immersed transient nonlinear residual is not finite"); if(initial<0) { initial=norm; initial_blocks=ResidualBlockNorms(Copy(rhs_)); } diagnostics_.residual_norm=norm; diagnostics_.nonlinear_iterations=it;
			if(norm<=std::max(options_.nonlinear_absolute_tolerance,options_.nonlinear_relative_tolerance*initial) && BlockReductionSatisfied(initial_blocks,Copy(rhs_),initial) && ControllersSatisfied()) { MarkSolved(); return true; }
			Check(VecCopy(rhs_,action_input_),"VecCopy rhs"); Check(KSPSetOperators(ksp_,jacobian_,jacobian_),"KSPSetOperators"); const auto begin=std::chrono::steady_clock::now(); Check(KSPSolve(ksp_,rhs_,update_),"KSPSolve"); diagnostics_.last_linear_solve_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
			PetscInt ki=0; Check(KSPGetIterationNumber(ksp_,&ki),"KSPGetIterationNumber"); diagnostics_.ksp_iterations+=ki; Check(KSPGetConvergedReason(ksp_,&diagnostics_.ksp_reason),"KSPGetConvergedReason"); if(diagnostics_.ksp_reason<=0) throw std::runtime_error("immersed transient KSP failed");
			Check(MatMult(jacobian_,update_,action_output_),"MatMult"); Check(VecAXPY(action_output_,-1.0,action_input_),"VecAXPY"); PetscReal linear=0,up=0; Check(VecNorm(action_output_,NORM_2,&linear),"VecNorm linear"); Check(VecNorm(update_,NORM_2,&up),"VecNorm update");
			ImmersedTransientFlowNewtonStep step; step.iteration=it; step.ksp_iterations=ki; step.ksp_reason=diagnostics_.ksp_reason; step.residual_norm=norm; step.update_norm=up; step.linear_relative_residual=norm>0?linear/norm:linear; diagnostics_.true_linear_relative_residual=step.linear_relative_residual;
			Check(VecCopy(state_,base_),"VecCopy line-search base");
			double d=1; bool accepted=false; while(d>=options_.minimum_damping) {
				Check(VecCopy(base_,state_),"VecCopy line-search restore"); Check(VecAXPY(state_,d,update_),"VecAXPY update"); Assemble(); PetscReal candidate=0; Check(VecNorm(rhs_,NORM_2,&candidate),"VecNorm candidate");
				if(std::isfinite(candidate)&&candidate<norm) { accepted=true; step.damping=d; diagnostics_.residual_norm=candidate; ++diagnostics_.nonlinear_iterations; diagnostics_.newton_steps.push_back(step); if(candidate<=std::max(options_.nonlinear_absolute_tolerance,options_.nonlinear_relative_tolerance*initial) && BlockReductionSatisfied(initial_blocks,Copy(rhs_),initial) && ControllersSatisfied()) { MarkSolved(); return true; } break; }
				d*=.5;
			}
			if(!accepted) Check(VecCopy(base_,state_),"VecCopy line-search failure restore");
			if(!accepted) { diagnostics_.newton_steps.push_back(step); throw std::runtime_error("immersed transient backtracking failed"); }
		} throw std::runtime_error("immersed transient nonlinear iteration cap reached"); } catch(...) { Rollback(); throw; }
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
		history_.reset(); frozen_body_force_.clear(); frozen_seed_.clear(); trial_ports_.clear(); diagnostics_.trial_state_hash_sha256.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.prepared=false; diagnostics_.trial_active=false; diagnostics_.idle=true; diagnostics_.committed=true; diagnostics_.converged=false; ++diagnostics_.finalize_count; ++diagnostics_.commit_count;
	}
	void Commit() { PrepareCommit(); FinalizeCommit(); }
	void Rollback()
	{
		if(!diagnostics_.trial_active) return;
		diagnostics_.prepared=false; prepared_global_.reset(); prepared_ports_.clear(); prepared_global_hash_.clear(); diagnostics_.prepared_hash_sha256.clear(); if(!frozen_seed_.empty()) SetVector(state_,frozen_seed_); trial_ports_=frozen_ports_; ResetAttemptWork(); ++diagnostics_.rollback_count; RefreshHashes();
	}
	void AbortTrial()
	{
		if(!diagnostics_.trial_active && !diagnostics_.prepared) return;
		Check(VecCopy(committed_,state_),"VecCopy abort committed state");
		diagnostics_.prepared=false; prepared_global_.reset(); prepared_ports_.clear(); prepared_global_hash_.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.input_hash_sha256.clear(); trial_ports_.clear(); history_.reset(); frozen_body_force_.clear(); frozen_seed_.clear(); diagnostics_.trial_active=false; diagnostics_.idle=true; diagnostics_.committed=true; ResetAttemptWork(); ++diagnostics_.abort_count; RefreshHashes();
	}
	void AbortPrepared() { AbortTrial(); }
	std::vector<PetscScalar> AssembledNegativeResidual() const { return Copy(rhs_); }
	std::vector<PetscScalar> AssembledJacobianAction(const std::vector<PetscScalar>& x) const { if(x.size()!=layout_.Rows()) throw std::invalid_argument("immersed transient Jacobian action size is invalid"); SetVector(action_input_,x); Check(MatMult(jacobian_,action_input_,action_output_),"MatMult action"); return Copy(action_output_); }
	ImmersedTransientFlowConservationDiagnostics ConservationDiagnostics() const { return MeasureConservation(); }

private:
	static void Check(PetscErrorCode c,const char* op) { if(c) throw std::runtime_error(std::string("PETSc ")+op+" failed: "+std::to_string(static_cast<long long>(c))); }
	void ValidateOptions() const { if(!options_.body_force||!std::isfinite(options_.parameters.density)||!(options_.parameters.density>0)||!std::isfinite(options_.parameters.dynamic_viscosity)||!(options_.parameters.dynamic_viscosity>0)||!std::isfinite(options_.wall_gamma0)||!(options_.wall_gamma0>0)||options_.nonlinear_maximum_iterations<=0||options_.ksp_maximum_iterations<=0||!std::isfinite(options_.ksp_relative_tolerance)||!(options_.ksp_relative_tolerance>0)||!std::isfinite(options_.nonlinear_relative_tolerance)||!(options_.nonlinear_relative_tolerance>0)||!std::isfinite(options_.nonlinear_absolute_tolerance)||!(options_.nonlinear_absolute_tolerance>0)||!std::isfinite(options_.flow_controller_relative_tolerance)||!(options_.flow_controller_relative_tolerance>0)||!std::isfinite(options_.flow_controller_absolute_tolerance_m3_s)||!(options_.flow_controller_absolute_tolerance_m3_s>0)||!std::isfinite(options_.flow_controller_reference_flow_m3_s)||!(options_.flow_controller_reference_flow_m3_s>0)||!std::isfinite(options_.minimum_damping)||!(options_.minimum_damping>0)||options_.minimum_damping>1||!std::isfinite(options_.lu_pivot_shift)||options_.lu_pivot_shift<0||!std::isfinite(options_.nonlinear_block_reduction)||options_.nonlinear_block_reduction<0) throw std::invalid_argument("immersed transient options are invalid"); ValidateImmersedNitscheWallLabels(options_.wall_labels); }
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
	void PreflightCatalogs() const { for(std::uint64_t c=0;c<domain_.Cells().size();++c){const auto classification=domain_.Cells()[c].classification; if(classification!=CellClassification::Inside&&classification!=CellClassification::Cut)continue; if(!volume_.Cell(c).usable)throw std::runtime_error("classified immersed transient cell has an unusable volume rule"); if(volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) volume_.ValidateUsableRule(domain_,c); else volume_.ValidateUsableCompactRule(domain_,c); const auto expected=VolumePointCount(c); std::size_t visited=0; ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint&){if(visited==std::numeric_limits<std::size_t>::max())throw std::overflow_error("immersed transient volume point count overflows");++visited;}); if(visited!=expected)throw std::logic_error("immersed transient volume point iterator count is invalid"); if(visited==0){if(classification==CellClassification::Inside)throw std::runtime_error("inside immersed transient cell has an empty volume rule");continue;} if(classification==CellClassification::Cut && volume_.Cell(c).diagnostics.estimated_reference_volume>0){if(!ghost_.Covered(c))throw std::runtime_error("positive immersed transient cut cell "+std::to_string(c)+" requires ghost coverage");surface_.ValidateUsableRule(domain_,c);}} }
	// Force callbacks are deliberately evaluated only here, at BeginTrial.  The
	// stored cell/point order is the canonical volume-rule order and therefore
	// neither retries nor later callback mutation can change an assembly.
	std::vector<std::vector<std::array<double,3>>> FreezeBodyForceCandidate() const { std::vector<std::vector<std::array<double,3>>> result(domain_.Cells().size()); for(std::uint64_t c=0;c<domain_.Cells().size();++c) if(Usable(c)){const auto element=domain_.Background().MaterializeElement(c); auto& forces=result[c]; const auto count=VolumePointCount(c); forces.reserve(count); ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint& point){const auto force=options_.body_force(EvaluateElementGeometry(element,point.parametric).physical); for(double value:force) if(!std::isfinite(value)) throw std::invalid_argument("immersed transient body force is not finite on the active quadrature"); forces.push_back(force);}); if(forces.size()!=count)throw std::logic_error("immersed transient frozen body-force point count is invalid");} return result; }
	std::vector<std::uint64_t> ControllerIds() const { std::vector<std::uint64_t> ids; for(const auto&p:options_.ports)if(p.control_mode==ImmersedFlowPortControlMode::FlowRate) ids.push_back(static_cast<std::uint64_t>(p.boundary_label)); std::sort(ids.begin(),ids.end()); return ids; }
	bool HasGauge() const { return options_.include_pressure_gauge && std::none_of(options_.ports.begin(),options_.ports.end(),[](const auto&p){return IsImmersedFlowPressureLike(p.control_mode);}); }
	void RequireIdle(const char* what) const { if(!diagnostics_.idle||diagnostics_.trial_active||diagnostics_.prepared) throw std::logic_error(std::string("cannot ")+what+" outside idle state"); }
	void RequireTrial(const char* what) const { if(!diagnostics_.trial_active||diagnostics_.prepared||!history_) throw std::logic_error(std::string("cannot ")+what+" without active transient trial"); }
	void ValidateGlobal(const ImmersedGlobalFlowState& x) const { if(!x.Valid()||x.GeometryIdentity()!=layout_.GeometryIdentity()||x.NodeIds()!=layout_.NodeIds()||x.PortIds()!=layout_.PortIds()||x.HasGaugeMultiplier()!=HasGauge()) throw std::invalid_argument("immersed transient committed global state layout mismatch"); }
	std::size_t VolumePointCount(std::uint64_t cell) const { if(volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) return volume_.Cell(cell).rule.Points().size(); return CompactCutCellVolumeLogicalPointCount(volume_.Cell(cell).compact_rule); }
	template<class Callback> void ForEachUsableVolumePoint(std::uint64_t cell,Callback&& callback) const { if(volume_.StorageMode()==CutCellVolumeQuadratureStorageMode::Expanded) for(const auto& point:volume_.UsableRule(domain_,cell).Points()) callback(point); else ForEachVolumePoint(volume_.UsableCompactRule(domain_,cell),std::forward<Callback>(callback)); }
	bool Usable(std::uint64_t cell) const { const auto& q=volume_.Cell(cell); const auto c=domain_.Cells()[cell].classification; return q.usable&&(c==CellClassification::Inside||c==CellClassification::Cut)&&VolumePointCount(cell)!=0; }
	static bool RuleHasLabel(const SurfaceQuadratureRule&r,int label){return std::any_of(r.Points().begin(),r.Points().end(),[&](const auto&p){return p.boundary_id==label;});}
	bool HasWallPoint(const SurfaceQuadratureRule&r) const {return std::any_of(r.Points().begin(),r.Points().end(),[this](const auto&p){return std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),p.boundary_id);});}
	ImmersedMaterialWallVelocityEvaluator MaterialVelocity() const { return [this](const SurfaceQuadraturePoint&,const ImmersedSurfaceQuadraturePointProvenance& p){return geometry_.Evaluation().WallVelocity(p.canonical_triangle,p.canonical_barycentric);}; }
	void CertifyZeroMaterialWallVelocity() const { for(std::uint64_t c=0;c<surface_.Cells().size();++c){const auto&r=surface_.UsableRule(domain_,c);const auto&v=surface_.UsableProvenance(domain_,c);for(std::size_t i=0;i<r.Points().size();++i)if(std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),r.Points()[i].boundary_id)){const auto w=geometry_.Evaluation().WallVelocity(v[i].canonical_triangle,v[i].canonical_barycentric);for(double x:w)if(!std::isfinite(x)||x!=0.0)throw std::invalid_argument("immersed transient fixed geometry requires exactly zero material wall velocity");}} }
	// This transaction intentionally has no PETSc command-line override path:
	// every effective solver setting below is fixed or an explicit option and is
	// included in InputHash.  In particular, no untracked options prefix exists.
	void CreatePetsc() { const PetscInt n=static_cast<PetscInt>(layout_.Rows()); try { Check(MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,300,nullptr,&jacobian_),"MatCreateSeqAIJ"); Check(MatSetOption(jacobian_,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),"MatSetOption"); Check(MatSetOption(jacobian_,MAT_IGNORE_ZERO_ENTRIES,PETSC_FALSE),"MatSetOption retain scalar diagonals"); for(const auto row:AppendedScalarRows()) Check(MatSetValue(jacobian_,row,row,0.0,INSERT_VALUES),"MatSetValue scalar structural diagonal"); Check(MatAssemblyBegin(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin scalar structural diagonals");Check(MatAssemblyEnd(jacobian_,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd scalar structural diagonals"); AuditAppendedScalarDiagonals(); Check(MatSetOption(jacobian_,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE),"MatSetOption ignore zero entries"); Check(VecCreateSeq(PETSC_COMM_SELF,n,&state_),"VecCreateSeq"); Check(VecDuplicate(state_,&committed_),"VecDuplicate");Check(VecDuplicate(state_,&prepared_),"VecDuplicate");Check(VecDuplicate(state_,&base_),"VecDuplicate");Check(VecDuplicate(state_,&rhs_),"VecDuplicate");Check(VecDuplicate(state_,&update_),"VecDuplicate");Check(VecDuplicate(state_,&action_input_),"VecDuplicate");Check(VecDuplicate(state_,&action_output_),"VecDuplicate"); Check(KSPCreate(PETSC_COMM_SELF,&ksp_),"KSPCreate"); Check(KSPSetType(ksp_,KSPGMRES),"KSPSetType"); Check(KSPGMRESSetRestart(ksp_,30),"KSPGMRESSetRestart"); Check(KSPSetPCSide(ksp_,PC_RIGHT),"KSPSetPCSide"); Check(KSPSetNormType(ksp_,KSP_NORM_UNPRECONDITIONED),"KSPSetNormType"); Check(KSPSetTolerances(ksp_,options_.ksp_relative_tolerance,1e-50,1e5,options_.ksp_maximum_iterations),"KSPSetTolerances"); PC pc=nullptr;Check(KSPGetPC(ksp_,&pc),"KSPGetPC");Check(PCSetType(pc,PCLU),"PCSetType");if(options_.lu_pivot_shift>0){Check(PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO),"PCFactorSetShiftType");Check(PCFactorSetShiftAmount(pc,options_.lu_pivot_shift),"PCFactorSetShiftAmount");} }catch(...){Destroy();throw;} }
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
	NavierStokesSystem BuildVolume(const Element&e,const std::vector<std::array<double,4>>&n,std::uint64_t cell)const{
		const auto& forces=frozen_body_force_.at(cell); std::size_t point=0;
		const NavierStokesBodyForceEvaluator frozen=[&forces,&point](const std::array<double,3>&){if(point>=forces.size())throw std::logic_error("immersed transient frozen body-force point count is invalid");return forces[point++];};
		ValidateTransientNavierStokesPreflight(options_.parameters);
		if(!history_||!history_->Valid()||!layout_.Valid()||history_->SourceGeometryIdentity()!=layout_.GeometryIdentity()||history_->TargetGeometryIdentity()!=layout_.GeometryIdentity()) throw std::invalid_argument("immersed velocity history does not match fixed geometry layout");
		for(const auto provenance:history_->Provenance()) if(provenance!=ImmersedVelocityHistoryProvenance::Committed) throw std::invalid_argument("fixed-geometry immersed velocity history must be committed");
		const double expected_target=CheckedTransientTargetTime(history_->SourceTimeS(),options_.parameters.dt);
		if(history_->TargetTimeS()!=expected_target||target_time_s_!=expected_target) throw std::invalid_argument("immersed velocity history target time does not match assembly time");
		const auto local=LocalizeImmersedVelocityHistory(e,layout_,*history_,target_time_s_);
		std::vector<std::array<double,4>> old(local.size()); for(std::size_t i=0;i<old.size();++i) for(int q=0;q<3;++q) old[i][q]=local[i][q];
		const auto result=BuildNavierStokesElementFromPoints(e,n,old,options_.parameters,[this,cell](const auto& consume){ForEachUsableVolumePoint(cell,consume);},frozen,NavierStokesResolvedMixedForm::Conservative);
		if(point!=forces.size()) throw std::logic_error("immersed transient frozen body-force point order is invalid");
		return result;
	}
	void Scatter(const std::vector<std::int32_t>&nodes,const NavierStokesSystem&s){std::vector<PetscInt>r;for(auto n:nodes)for(int q=0;q<4;++q)r.push_back(Dof(n,q));Check(MatSetValues(jacobian_,r.size(),r.data(),r.size(),r.data(),s.jacobian.data(),ADD_VALUES),"MatSetValues");Check(VecSetValues(rhs_,r.size(),r.data(),s.negative_residual.data(),ADD_VALUES),"VecSetValues");}
	void ScatterBlock(const std::vector<std::int32_t>&nodes,const std::vector<PetscScalar>&a,const std::vector<PetscScalar>&b){NavierStokesSystem s{a,b};Scatter(nodes,s);}
	void SubtractAndScatter(const std::vector<std::int32_t>&nodes,NavierStokesSystem a,const NavierStokesSystem&b){for(std::size_t i=0;i<a.jacobian.size();++i)a.jacobian[i]-=b.jacobian[i];for(std::size_t i=0;i<a.negative_residual.size();++i)a.negative_residual[i]-=b.negative_residual[i];Scatter(nodes,a);}
	void ScatterPort(const Element&e,const std::vector<std::array<double,4>>&n,const SurfaceQuadratureRule&r,std::size_t i){const auto& p=options_.ports[i];const auto x=BuildImmersedFlowPortElement(e,r,p.boundary_label,p.control_mode,p.value,n);std::vector<PetscInt> rows;for(auto id:e.connectivity)for(int q=0;q<4;++q)rows.push_back(Dof(id,q));if(IsImmersedFlowPressureLike(p.control_mode)){Check(VecSetValues(rhs_,rows.size(),rows.data(),x.negative_residual.data(),ADD_VALUES),"VecSetValues port");return;}const auto scalar=static_cast<PetscInt>(layout_.ControllerRow(static_cast<std::uint64_t>(p.boundary_label)));const double lambda=Value(state_,scalar);for(std::size_t k=0;k<rows.size();++k){const double c=PetscRealPart(x.flow_coefficient[k]);if(c){Check(MatSetValue(jacobian_,rows[k],scalar,c,ADD_VALUES),"MatSetValue port col");Check(MatSetValue(jacobian_,scalar,rows[k],c,ADD_VALUES),"MatSetValue port row");Check(VecSetValue(rhs_,rows[k],-lambda*c,ADD_VALUES),"VecSetValue port");}}}
	void MeasurePorts(){for(std::size_t p=0;p<options_.ports.size();++p){ImmersedFlowPortMeasurement t{};for(std::uint64_t c=0;c<surface_.Cells().size();++c){const auto&r=surface_.UsableRule(domain_,c);if(!RuleHasLabel(r,options_.ports[p].boundary_label))continue;const auto m=MeasureImmersedFlowPortElement(domain_.Background().MaterializeElement(c),r,options_.ports[p].boundary_label,Gather(domain_.Background().MaterializeElement(c)),options_.parameters.dynamic_viscosity);t.area_m2+=m.area_m2;t.outward_flow_m3_s+=m.outward_flow_m3_s;t.mean_pressure_pa+=m.mean_pressure_pa*m.area_m2;t.mean_normal_traction_pa+=m.mean_normal_traction_pa*m.area_m2;t.mean_velocity_squared_m2_s2+=m.mean_velocity_squared_m2_s2*m.area_m2;}if(!(t.area_m2>0))throw std::runtime_error("immersed transient port area is zero");t.mean_pressure_pa/=t.area_m2;t.mean_normal_traction_pa/=t.area_m2;t.mean_velocity_squared_m2_s2/=t.area_m2;auto&d=trial_ports_.at(p);d.measurement=t;d.measurement_valid=true;if(options_.ports[p].control_mode==ImmersedFlowPortControlMode::FlowRate){d.controller_error=options_.ports[p].value-t.outward_flow_m3_s;d.multiplier=Value(state_,d.multiplier_row);Check(VecSetValue(rhs_,d.multiplier_row,d.controller_error,ADD_VALUES),"VecSetValue flow constraint");}}}
	bool ControllersSatisfied()const{for(std::size_t i=0;i<options_.ports.size();++i)if(options_.ports[i].control_mode==ImmersedFlowPortControlMode::FlowRate){const auto&d=trial_ports_.at(i);const double tol=options_.flow_controller_absolute_tolerance_m3_s+options_.flow_controller_relative_tolerance*std::max(std::abs(options_.ports[i].value),options_.flow_controller_reference_flow_m3_s);if(std::abs(d.controller_error)>tol)return false;}return true;}
	void BuildGaugeWeights(){gauge_weights_.assign(layout_.NodeIds().size(),0);for(std::uint64_t c=0;c<domain_.Cells().size();++c)if(Usable(c)){const auto e=domain_.Background().MaterializeElement(c);ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint&p){const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false);const double m=p.weight*b.raw_determinant;diagnostics_.pressure_measure+=m;for(std::size_t a=0;a<e.connectivity.size();++a)gauge_weights_[layout_.LocalNode(e.connectivity[a])]+=b.value[a]*m;});}}
	void SetHistoryDiagnostics(ImmersedTransientFlowDiagnostics& diagnostics,const ImmersedVelocityHistory* history) const {diagnostics.identity_history_nodes=history?history->NodeIds().size():0;diagnostics.committed_history_nodes=0;diagnostics.extended_history_nodes=0;diagnostics.missing_history_nodes=layout_.NodeIds().size();diagnostics.history_hash_sha256.clear();if(!history)return;diagnostics.history_hash_sha256=history->HashSha256();for(const auto p:history->Provenance())if(p==ImmersedVelocityHistoryProvenance::Committed)++diagnostics.committed_history_nodes;else if(p==ImmersedVelocityHistoryProvenance::Extended)++diagnostics.extended_history_nodes;diagnostics.missing_history_nodes=layout_.NodeIds().size()-diagnostics.identity_history_nodes;}
	void RefreshHistoryDiagnostics(){SetHistoryDiagnostics(diagnostics_,history_.get());}
	void InsertGauge(){const PetscInt q=GaugeDof();double g=0;for(std::size_t i=0;i<gauge_weights_.size();++i){const double w=gauge_weights_[i];g+=w*Value(state_,static_cast<PetscInt>(4*i+3));Check(MatSetValue(jacobian_,static_cast<PetscInt>(4*i+3),q,w,ADD_VALUES),"MatSetValue gauge");Check(MatSetValue(jacobian_,q,static_cast<PetscInt>(4*i+3),w,ADD_VALUES),"MatSetValue gauge");Check(VecSetValue(rhs_,static_cast<PetscInt>(4*i+3),-w*Value(state_,q),ADD_VALUES),"VecSetValue gauge");}Check(VecSetValue(rhs_,q,-g,ADD_VALUES),"VecSetValue gauge");diagnostics_.pressure_gauge_defect=std::abs(g);}
	ImmersedTransientFlowConservationDiagnostics MeasureConservation()const{ImmersedTransientFlowConservationDiagnostics d;for(std::uint64_t c=0;c<domain_.Cells().size();++c)if(Usable(c)){const auto e=domain_.Background().MaterializeElement(c);const auto n=Gather(e);ForEachUsableVolumePoint(c,[&](const VolumeQuadraturePoint&p){const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false);double div=0;for(std::size_t a=0;a<n.size();++a)for(int q=0;q<3;++q)div+=n[a][q]*b.gradient[a][q];d.endpoint_volume_divergence_m3_s+=p.weight*b.raw_determinant*div;});if(domain_.Cells()[c].classification==CellClassification::Cut)for(const auto&p:surface_.UsableRule(domain_,c).Points()){const auto b=EvaluateBasis(e,p.parametric[0],p.parametric[1],p.parametric[2],false);double f=0;for(std::size_t a=0;a<n.size();++a)for(int q=0;q<3;++q)f+=n[a][q]*b.value[a]*p.normal[q]*p.weight;d.surface_flow_by_boundary_label_m3_s[p.boundary_id]+=f;d.total_surface_outward_flow_m3_s+=f;if(std::binary_search(options_.wall_labels.begin(),options_.wall_labels.end(),p.boundary_id))d.wall_outward_flow_m3_s+=f;else if(std::any_of(options_.ports.begin(),options_.ports.end(),[&](const auto&port){return port.boundary_label==p.boundary_id;}))d.open_port_outward_flow_m3_s+=f;else throw std::runtime_error("immersed transient conservation found an unconfigured surface label");}}const double scale=std::max(options_.flow_controller_reference_flow_m3_s,std::abs(d.open_port_outward_flow_m3_s));d.wall_relative_leakage_m3_s=d.wall_outward_flow_m3_s;d.normalized_wall_leakage=std::abs(d.wall_outward_flow_m3_s)/scale;d.normalized_open_balance=std::abs(d.endpoint_volume_divergence_m3_s-d.total_surface_outward_flow_m3_s)/scale;return d;}
	void RefreshHashes(){diagnostics_.committed_state_hash_sha256=committed_global_->HashSha256();diagnostics_.trial_state_hash_sha256=diagnostics_.trial_active?HashVector(state_):std::string{};if(!history_){diagnostics_.identity_history_nodes=diagnostics_.committed_history_nodes=diagnostics_.extended_history_nodes=0;diagnostics_.missing_history_nodes=layout_.NodeIds().size();diagnostics_.history_hash_sha256.clear();}}
	void RefreshHashesNoexcept()noexcept{try{RefreshHashes();}catch(...){}}
	std::string HashVector(Vec v)const{Sha256 h;immersed_transient_detail::AppendString(h,"ImmersedTransientTrial/v1");for(auto x:Copy(v))h.AppendNormalizedDouble(PetscRealPart(x));return h.Hex();}
	void SyncPortMultipliers(std::vector<ImmersedTransientFlowDiagnostics::Port>& ports,Vec state) const { for(auto& port:ports) if(port.multiplier_row>=0) port.multiplier=Value(state,port.multiplier_row); }
	void InvalidatePortMeasurements(std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) const { for(auto& port:ports) { port.measurement={}; port.measurement_valid=false; port.controller_error=0.0; } }
	// Explicit Assemble() is permitted as pre-solve diagnostic work.  It is not
	// solve-attempt work; reset at solve/rollback boundaries makes retries exact.
	static void ResetAttemptWork(ImmersedTransientFlowDiagnostics& diagnostics) noexcept { diagnostics.converged=false; diagnostics.solved_state_hash_sha256.clear(); diagnostics.attempt_hash_sha256.clear(); diagnostics.residual_norm=0.0; diagnostics.true_linear_relative_residual=0.0; diagnostics.nonlinear_iterations=0; diagnostics.ksp_iterations=0; diagnostics.ksp_reason=KSP_CONVERGED_ITERATING; diagnostics.last_assembly_seconds=0.0; diagnostics.last_linear_solve_seconds=0.0; diagnostics.attempt_assembly_count=0; diagnostics.volume_cells=diagnostics.surface_cells=diagnostics.ghost_faces=0; diagnostics.newton_steps.clear(); }
	void ResetAttemptWork() noexcept { ResetAttemptWork(diagnostics_); }
	void InvalidateSolved() noexcept { ResetAttemptWork(); }
	void ClearIdlePublicationDiagnostics() noexcept { diagnostics_.input_hash_sha256.clear(); diagnostics_.prepared_hash_sha256.clear(); diagnostics_.trial_state_hash_sha256.clear(); diagnostics_.history_hash_sha256.clear(); diagnostics_.target_time_s=diagnostics_.dt_s=diagnostics_.pressure_gauge_defect=0.0; diagnostics_.prepared=false; ResetAttemptWork(); }
	void MarkSolved() { RefreshHashes(); diagnostics_.solved_state_hash_sha256=diagnostics_.trial_state_hash_sha256; diagnostics_.converged=true; diagnostics_.attempt_hash_sha256=AttemptHash(); }
	static void AppendPorts(Sha256& h,const std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) { h.AppendLittleEndian64(ports.size()); for(const auto&p:ports){immersed_transient_detail::AppendString(h,p.id);h.AppendLittleEndian32(static_cast<std::uint32_t>(p.boundary_label));h.AppendLittleEndian32(static_cast<std::uint32_t>(p.control_mode));h.AppendNormalizedDouble(p.target);h.AppendNormalizedDouble(p.multiplier);h.AppendNormalizedDouble(p.controller_error);h.AppendLittleEndian64(static_cast<std::uint64_t>(p.multiplier_row));h.AppendNormalizedDouble(p.measurement.area_m2);h.AppendNormalizedDouble(p.measurement.outward_flow_m3_s);h.AppendNormalizedDouble(p.measurement.mean_pressure_pa);h.AppendNormalizedDouble(p.measurement.mean_normal_traction_pa);h.AppendNormalizedDouble(p.measurement.mean_velocity_squared_m2_s2);h.AppendLittleEndian32(p.measurement_valid?1:0);}}
	static void AppendFrozenBodyForce(Sha256& h,const std::vector<std::vector<std::array<double,3>>>& forces) { h.AppendLittleEndian64(forces.size()); for(const auto& cell:forces) { h.AppendLittleEndian64(cell.size()); for(const auto& force:cell) for(double value:force) h.AppendNormalizedDouble(value); } }
	void AppendFrozenBodyForce(Sha256& h) const { AppendFrozenBodyForce(h,frozen_body_force_); }
	static void AppendFixedSolverConfiguration(Sha256& h) { immersed_transient_detail::AppendString(h,"KSPGMRES"); h.AppendLittleEndian64(30); immersed_transient_detail::AppendString(h,"PCLU"); immersed_transient_detail::AppendString(h,"PC_RIGHT"); immersed_transient_detail::AppendString(h,"KSP_NORM_UNPRECONDITIONED"); h.AppendNormalizedDouble(1e-50); h.AppendNormalizedDouble(1e5); }
	std::array<double,3> ResidualBlockNorms(const std::vector<PetscScalar>& residual) const { if(residual.size()!=diagnostics_.total_dofs) throw std::invalid_argument("immersed transient residual block vector size is invalid"); std::array<double,3> norms{{0.0,0.0,0.0}}; for(std::size_t row=0;row<residual.size();++row){const double value=PetscRealPart(residual[row]);const std::size_t block=row>=diagnostics_.physical_dofs?2:((row%4)==3?1:0);norms[block]+=value*value;}for(auto& value:norms)value=std::sqrt(value);return norms; }
	bool BlockReductionSatisfied(const std::array<double,3>& initial,const std::vector<PetscScalar>& residual,double global_initial) const { if(options_.nonlinear_block_reduction==0.0)return true; const auto current=ResidualBlockNorms(residual); const double zero_block_tolerance=options_.nonlinear_absolute_tolerance+options_.nonlinear_relative_tolerance*global_initial; for(std::size_t block=0;block<current.size();++block)if(initial[block]>0.0?current[block]>initial[block]/options_.nonlinear_block_reduction:current[block]>zero_block_tolerance)return false;return true; }
	std::string InputHash(const NavierStokesParameters& parameters,const ImmersedVelocityHistory& history,double target_time_s,std::uint64_t target_index,const std::vector<ImmersedTransientFlowDiagnostics::Port>& frozen_ports,const std::vector<std::vector<std::array<double,3>>>& frozen_forces,const std::vector<PetscScalar>& frozen_seed) const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientInput/v4"); immersed_transient_detail::AppendString(h,layout_.HashSha256()); immersed_transient_detail::AppendString(h,geometry_.GeometryIdentitySha256()); h.AppendLittleEndian32(static_cast<std::uint32_t>(volume_.StorageMode())); immersed_transient_detail::AppendString(h,history.HashSha256()); h.AppendNormalizedDouble(committed_global_->TimeS());h.AppendLittleEndian64(committed_global_->Index());h.AppendNormalizedDouble(target_time_s);h.AppendLittleEndian64(target_index);h.AppendNormalizedDouble(parameters.density);h.AppendNormalizedDouble(parameters.dynamic_viscosity);h.AppendNormalizedDouble(parameters.dt);h.AppendNormalizedDouble(options_.wall_gamma0);h.AppendLittleEndian32(options_.include_pressure_gauge?1:0);h.AppendLittleEndian64(options_.nonlinear_maximum_iterations);h.AppendLittleEndian64(options_.ksp_maximum_iterations);h.AppendNormalizedDouble(options_.ksp_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_relative_tolerance);h.AppendNormalizedDouble(options_.nonlinear_absolute_tolerance);h.AppendNormalizedDouble(options_.flow_controller_relative_tolerance);h.AppendNormalizedDouble(options_.flow_controller_absolute_tolerance_m3_s);h.AppendNormalizedDouble(options_.flow_controller_reference_flow_m3_s);h.AppendNormalizedDouble(options_.minimum_damping);h.AppendNormalizedDouble(options_.lu_pivot_shift);h.AppendNormalizedDouble(options_.nonlinear_block_reduction);AppendFixedSolverConfiguration(h);AppendPorts(h,frozen_ports);AppendFrozenBodyForce(h,frozen_forces);for(auto x:frozen_seed)h.AppendNormalizedDouble(PetscRealPart(x));return h.Hex(); }
	std::string InputHash() const { if(!history_) return {}; return InputHash(options_.parameters,*history_,target_time_s_,target_index_,frozen_ports_,frozen_body_force_,frozen_seed_); }
	std::string AttemptHash() const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientAttempt/v2"); immersed_transient_detail::AppendString(h,diagnostics_.input_hash_sha256); immersed_transient_detail::AppendString(h,HashVector(state_)); AppendPorts(h,trial_ports_); h.AppendLittleEndian64(diagnostics_.attempt_assembly_count);h.AppendLittleEndian64(diagnostics_.newton_steps.size());for(const auto&s:diagnostics_.newton_steps){h.AppendLittleEndian64(s.iteration);h.AppendLittleEndian64(s.ksp_iterations);h.AppendLittleEndian32(static_cast<std::uint32_t>(s.ksp_reason));h.AppendNormalizedDouble(s.residual_norm);h.AppendNormalizedDouble(s.update_norm);h.AppendNormalizedDouble(s.linear_relative_residual);h.AppendNormalizedDouble(s.damping);} h.AppendLittleEndian64(diagnostics_.ksp_iterations);h.AppendLittleEndian64(diagnostics_.nonlinear_iterations);h.AppendLittleEndian32(static_cast<std::uint32_t>(diagnostics_.ksp_reason));h.AppendNormalizedDouble(diagnostics_.residual_norm);h.AppendNormalizedDouble(diagnostics_.true_linear_relative_residual);return h.Hex(); }
	std::string PreparedHash(const ImmersedGlobalFlowState& global,const std::vector<ImmersedTransientFlowDiagnostics::Port>& ports) const { Sha256 h; immersed_transient_detail::AppendString(h,"ImmersedTransientPrepared/v2"); immersed_transient_detail::AppendString(h,diagnostics_.attempt_hash_sha256); immersed_transient_detail::AppendString(h,global.HashSha256());AppendPorts(h,ports);return h.Hex(); }
	const MovingCutGeometry& geometry_; const CartesianDomainClassification& domain_; const CutCellVolumeQuadratureCatalog& volume_; const ImmersedSurfaceQuadratureCatalog& surface_; const CutCellGhostPenaltyCatalog& ghost_; ImmersedTransientFlowOptions options_; ImmersedActiveLayout layout_; std::vector<double> gauge_weights_; std::unique_ptr<ImmersedGlobalFlowState> committed_global_,prepared_global_; std::unique_ptr<ImmersedVelocityHistory> history_; std::vector<PetscScalar> frozen_seed_; std::vector<std::vector<std::array<double,3>>> frozen_body_force_; std::vector<ImmersedTransientFlowDiagnostics::Port> frozen_ports_,trial_ports_,prepared_ports_; std::string prepared_global_hash_; double target_time_s_=0; std::uint64_t target_index_=0; Mat jacobian_=nullptr;Vec state_=nullptr,committed_=nullptr,prepared_=nullptr,base_=nullptr,rhs_=nullptr,update_=nullptr,action_input_=nullptr,action_output_=nullptr;KSP ksp_=nullptr;ImmersedTransientFlowDiagnostics diagnostics_{};bool fail_next_prepare_=false;
};

} // namespace iga

#endif

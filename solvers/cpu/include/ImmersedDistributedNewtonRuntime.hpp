#ifndef IGA_IMMERSED_DISTRIBUTED_NEWTON_RUNTIME_HPP
#define IGA_IMMERSED_DISTRIBUTED_NEWTON_RUNTIME_HPP

#include "ImmersedStaticDistributedOperator.hpp"
#include "PetscSolverOptions.hpp"
#include "PetscPhaseProfile.hpp"

namespace iga {

// Borrowed source database. A null database preserves the global-source API.
// A supplied database must outlive this numerical runtime.
struct ImmersedNewtonOptionsSource { PetscOptions value=nullptr; };

// Shared MPI Newton solve with owned committed/prepared states. Public MPI
// operations are collective on the borrowed communicator. FinalizeCommit is
// the nonthrowing publication half of a collectively successful preparation.
template<class Operator> class ImmersedDistributedNewtonRuntime {
public:
	template<class... Arguments>
	ImmersedDistributedNewtonRuntime(MPI_Comm communicator,const char* prefix,Arguments&&... arguments)
		: ImmersedDistributedNewtonRuntime(communicator,ImmersedNewtonOptionsSource{},prefix,std::forward<Arguments>(arguments)...)
	{}
	template<class... Arguments>
	ImmersedDistributedNewtonRuntime(MPI_Comm communicator,ImmersedNewtonOptionsSource source,const char* prefix,Arguments&&... arguments)
		: communicator_(communicator)
	{
		op_ = AllocateCollectiveRuntime<Operator>(communicator_,communicator_,std::forward<Arguments>(arguments)...);
		const auto& options = op_->Options();
		try {
			std::string checked_prefix, inherited_prefix;
			CollectiveLocalStage(communicator_,"immersed Newton prefix preflight",[&] {
				if (!prefix || !*prefix) throw std::invalid_argument("immersed Newton options prefix is empty");
				inherited_prefix = prefix;
				checked_prefix = options.solver_options_prefix.empty() ? inherited_prefix : options.solver_options_prefix;
			});
			RequireCollectiveSameText(communicator_,"immersed Newton prefix agreement",checked_prefix);
			CollectiveLocalStage(communicator_,"distributed static diagnostics storage",[&] { diagnostics_ = op_->Diagnostics(); });
			for (auto entry : {&committed_,&prepared_,&update_,&linear_rhs_,&linear_action_})
				Check("distributed static vector create",VecDuplicate(op_->Assembly().State(),entry));
			Check("distributed static committed clear",VecSet(committed_,0.0));
			RequireCollectivePetscOptions(communicator_,source.value);
			Check("distributed static KSP create",KSPCreate(communicator_,&solver_));
			Check("distributed static KSP type",KSPSetType(solver_,KSPGMRES));
			// A shifted LU is only a preconditioner. Stop on the physical residual,
			// rather than allowing a small preconditioned norm to hide its error.
			Check("distributed static PC side",KSPSetPCSide(solver_,PC_RIGHT));
			Check("distributed static KSP norm",KSPSetNormType(solver_,KSP_NORM_UNPRECONDITIONED));
			Check("distributed static KSP tolerances",KSPSetTolerances(solver_,options.ksp_relative_tolerance,PETSC_DEFAULT,PETSC_DEFAULT,options.ksp_maximum_iterations));
			PC pc = nullptr;
			Check("distributed static PC lookup",KSPGetPC(solver_,&pc));
			Check("distributed static PC type",PCSetType(pc,PCLU));
			// Keep the same pivoting backend for single- and multi-rank reference
			// solves. A caller may select another PC/backend through the prefix.
			Check("distributed static LU backend",PCFactorSetMatSolverType(pc,MATSOLVERMUMPS));
			if (options.lu_pivot_shift > 0.0) {
				Check("distributed static LU shift type",PCFactorSetShiftType(pc,MAT_SHIFT_NONZERO));
				Check("distributed static LU shift amount",PCFactorSetShiftAmount(pc,options.lu_pivot_shift));
			}
			solver_options_ = AllocateCollectiveRuntime<PetscSolverOptions>(communicator_, communicator_, checked_prefix,
				source.value, PetscOptionEntries{}, std::set<std::string>{}, inherited_prefix, false);
			solver_options_->Attach(solver_);
			solver_options_->Call("distributed static KSP options", [&] { return KSPSetFromOptions(solver_); });
			solver_options_->RecordUsed();
		} catch (...) { Release(); throw; }
	}
	~ImmersedDistributedNewtonRuntime() { Release(); }
	ImmersedDistributedNewtonRuntime(const ImmersedDistributedNewtonRuntime&) = delete;
	ImmersedDistributedNewtonRuntime& operator=(const ImmersedDistributedNewtonRuntime&) = delete;
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	PetscKspConfiguration SolverConfiguration() const { RequireOpen(); return CaptureKspConfiguration(solver_); }
	MPI_Comm Communicator() const noexcept { return communicator_; }
	const std::vector<ImmersedFlowPortDefinition>& PortDefinitions() const noexcept { return op_->Options().ports; }
	void SetPortControlValue(const std::string& id,double value)
	{
		RequireStateAgreement();
		CollectiveLocalStage(communicator_,"distributed static port trial guard",[&] {
			RequireOpen();
			if (diagnostics_.trial_active) throw std::logic_error("cannot change immersed port control during an active trial");
		});
		op_->SetPortControlValue(id,value);
		for (auto& port : diagnostics_.ports) if (port.id == id) port.target = value;
		diagnostics_.converged = false;
	}
	const ImmersedStaticFlowSetup& Topology() const noexcept { return op_->Topology(); }
	std::size_t OwnedStencilCount() const noexcept { return op_->Assembly().OwnedStencils().size(); }
	std::size_t RequiredStateRows() const noexcept { return op_->Assembly().RequiredRows().size(); }
	PetscInt RowBegin() const noexcept { return op_->Assembly().RowBegin(); }
	PetscInt RowEnd() const noexcept { return op_->Assembly().RowEnd(); }
	// Borrowed vectors are for read-only inspection; input uses the owned setter.
	Vec State() const { return op_->Assembly().State(); }
	Vec CommittedState() const { RequireOpen(); return committed_; }
	ImmersedStaticFlowConservationDiagnostics ConservationDiagnostics() const { return op_->ConservationDiagnostics(); }

	void SetCommittedOwnedState(const std::vector<PetscScalar>& values)
	{
		RequireStateAgreement();
		CollectiveLocalStage(communicator_,"distributed static owned input",[&] {
			RequireOpen();
			if (diagnostics_.trial_active || values.size() != static_cast<std::size_t>(RowEnd()-RowBegin()))
				throw std::invalid_argument("invalid static owned state size or active trial");
			for (auto value : values) if (!std::isfinite(PetscRealPart(value))) throw std::invalid_argument("nonfinite static owned state");
		});
		CollectiveLocalStage(communicator_,"distributed static input copy",[&] {
			PetscScalar* target = nullptr;
			if (VecGetArray(prepared_,&target)) throw std::runtime_error("cannot acquire static input vector");
			std::copy(values.begin(),values.end(),target);
			if (VecRestoreArray(prepared_,&target)) throw std::runtime_error("cannot restore static input vector");
		});
		try { Check("distributed static input state",VecCopy(prepared_,State())); }
		catch (...) { Check("distributed static input rollback",VecCopy(committed_,State())); throw; }
		std::swap(committed_,prepared_);
		diagnostics_.converged = false;
	}
	void Assemble()
	{
		op_->Assemble();
		CollectiveLocalStage(communicator_,"distributed static assembly diagnostics",[&] {
			if (candidate_assembly_ && fail_next_candidate_) {
				fail_next_candidate_ = false;
				throw std::runtime_error("injected distributed Newton candidate failure");
			}
			const auto& source = op_->Diagnostics();
			auto candidate = diagnostics_;
			candidate.volume_cells = source.volume_cells; candidate.surface_cells = source.surface_cells; candidate.ghost_faces = source.ghost_faces;
			candidate.pressure_measure = source.pressure_measure; candidate.constant_pressure_defect = source.constant_pressure_defect;
			candidate.ports = source.ports; candidate.wall_selected_points = source.wall_selected_points;
			candidate.wall_selected_points_by_label = source.wall_selected_points_by_label;
			candidate.scalar_diagonal_structure_verified = source.scalar_diagonal_structure_verified;
			candidate.last_assembly_seconds = source.last_assembly_seconds; candidate.aggregate_assembly_seconds = source.aggregate_assembly_seconds;
			diagnostics_ = std::move(candidate);
		});
	}
	bool SolveTrial()
	{
		RequireStateAgreement();
		CollectiveLocalStage(communicator_,"distributed static trial preflight",[&] {
			RequireOpen(); if (diagnostics_.trial_active) throw std::logic_error("static trial is already active");
		});
		Check("distributed static begin trial",VecCopy(committed_,State()));
		diagnostics_.trial_active = true; diagnostics_.committed = false; diagnostics_.prepared = false; diagnostics_.converged = false;
		diagnostics_.nonlinear_iterations = diagnostics_.ksp_iterations = 0; diagnostics_.ksp_reason = KSP_CONVERGED_ITERATING;
		diagnostics_.damping = diagnostics_.residual_norm = 0.0; diagnostics_.newton_steps.clear();
		try {
			double initial = -1.0; std::array<double,3> initial_blocks{};
			for (PetscInt iteration = 0; iteration < Options().nonlinear_maximum_iterations; ++iteration) {
				Assemble(); const double residual = Norm(op_->Assembly().Residual());
				if (initial < 0.0) { initial = residual; initial_blocks = BlockNorms(); }
				diagnostics_.residual_norm = residual; diagnostics_.nonlinear_iterations = iteration;
				if (Converged(residual,initial,initial_blocks)) { diagnostics_.converged = true; return true; }
				ImmersedStaticFlowNewtonStep step; step.iteration = iteration; step.residual_norm = residual;
				Check("distributed static original linear RHS",VecCopy(op_->Assembly().Residual(),linear_rhs_));
				if (!matrix_options_bound_) {
					// Matrix() becomes available after the first assembly. Bind
					// before setup creates a backend factor, not in construction.
					solver_options_->Attach(op_->Assembly().Matrix());
					matrix_options_bound_=true;
				}
				Check("distributed static operators",KSPSetOperators(solver_,op_->Assembly().Matrix(),op_->Assembly().Matrix()));
				RequireKspFactorBackend(solver_,op_->Assembly().Matrix(),communicator_);
				const auto start = std::chrono::steady_clock::now();
				solver_options_->Call("distributed static solve", [&] { return SolveProfiledKsp(solver_,op_->Assembly().Residual(),update_); });
				solver_options_->RecordUsed();
				diagnostics_.last_linear_solve_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
				diagnostics_.aggregate_linear_solve_seconds += diagnostics_.last_linear_solve_seconds;
				Check("distributed static KSP reason",KSPGetConvergedReason(solver_,&step.ksp_reason));
				Check("distributed static KSP iterations",KSPGetIterationNumber(solver_,&step.ksp_iterations));
				PetscReal reported_residual = 0.0;
				Check("distributed static KSP residual",KSPGetResidualNorm(solver_,&reported_residual));
				step.ksp_residual_norm = reported_residual;
				diagnostics_.ksp_reason = step.ksp_reason; diagnostics_.ksp_iterations += step.ksp_iterations;
				CollectiveLocalStage(communicator_,"distributed static linear convergence",[&] {
					if (!std::isfinite(step.ksp_residual_norm)) throw std::runtime_error("nonfinite static KSP residual");
					if (step.ksp_reason <= 0) {
						std::string message="static immersed-flow KSP failed with reason "+std::to_string(static_cast<int>(step.ksp_reason));
						// These local getters do not enter MPI. Preserve the original
						// KSP reason even if the additional PC lookup is unavailable.
						PC pc=nullptr;PCFailedReason failure=PC_NOERROR;
						if (!KSPGetPC(solver_,&pc) && !PCGetFailedReason(pc,&failure)) {
							message+="; PC failure reason "+std::to_string(static_cast<int>(failure));
							switch (failure) {
							case PC_FACTOR_STRUCT_ZEROPIVOT: message+=" (structural zero pivot)";break;
							case PC_FACTOR_NUMERIC_ZEROPIVOT: message+=" (numerical zero pivot)";break;
							case PC_FACTOR_OUTMEMORY: message+=" (factor memory exhausted)";break;
							default: break;
							}
						}
						throw std::runtime_error(message);
					}
				});
				step.update_norm = Norm(update_);
				Check("distributed static true linear action",MatMult(op_->Assembly().Matrix(),update_,linear_action_));
				Check("distributed static true linear residual",VecAXPY(linear_action_,-1.0,linear_rhs_));
				step.linear_residual_norm = Norm(linear_action_); step.linear_relative_residual = step.linear_residual_norm/residual;
				bool accepted = false;
				for (double damping = 1.0; damping >= Options().minimum_damping; damping *= 0.5) {
					Check("distributed static candidate",VecAXPY(State(),damping,update_));
					const auto candidate_start = std::chrono::steady_clock::now();
					AssembleCandidate();
					diagnostics_.last_line_search_candidate_assembly_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-candidate_start).count();
					diagnostics_.aggregate_line_search_candidate_assembly_seconds += diagnostics_.last_line_search_candidate_assembly_seconds;
					const double candidate = Norm(op_->Assembly().Residual());
					if (candidate < residual) {
						accepted = true; step.damping = damping; step.candidate_residual_norm = candidate;
						diagnostics_.damping = damping; diagnostics_.residual_norm = candidate; ++diagnostics_.nonlinear_iterations;
						RecordStep(step);
						if (Converged(candidate,initial,initial_blocks)) { diagnostics_.converged = true; return true; }
						break;
					}
					Check("distributed static candidate undo",VecAXPY(State(),-damping,update_));
				}
				if (!accepted) { RecordStep(step); throw std::runtime_error("static immersed-flow backtracking failed"); }
			}
			throw std::runtime_error("static immersed-flow nonlinear iteration cap reached");
		} catch (...) { Rollback(); throw; }
	}
	void PrepareCommit()
	{
		RequireStateAgreement();
		CollectiveLocalStage(communicator_,"distributed static prepare preflight",[&] {
			RequireOpen();
			if (!diagnostics_.trial_active || !diagnostics_.converged || diagnostics_.prepared)
				throw std::logic_error("cannot prepare an unconverged or already prepared static-flow trial");
			if (fail_next_prepare_) { fail_next_prepare_ = false; throw std::runtime_error("injected distributed static prepare failure"); }
		});
		Check("distributed static prepare copy",VecCopy(State(),prepared_));
		diagnostics_.prepared = true; ++diagnostics_.prepare_count;
	}
	void FinalizeCommit() noexcept
	{
		if (!diagnostics_.prepared) return;
		std::swap(committed_,prepared_); diagnostics_.prepared = false; diagnostics_.trial_active = false;
		diagnostics_.committed = true; ++diagnostics_.commit_count; ++diagnostics_.finalize_count;
	}
	void Commit() { PrepareCommit(); FinalizeCommit(); }
	void AbortPrepared() noexcept { diagnostics_.prepared = false; }
	void Rollback()
	{
		RequireStateAgreement();
		CollectiveLocalStage(communicator_,"distributed static rollback preflight",[&] { RequireOpen(); });
		if (!diagnostics_.trial_active) return;
		Check("distributed static rollback",VecCopy(committed_,State()));
		diagnostics_.prepared = false; diagnostics_.trial_active = false; diagnostics_.committed = true; diagnostics_.converged = false;
		++diagnostics_.rollback_count;
	}
	void FailNextCandidateForTesting() noexcept { fail_next_candidate_ = true; }
	void FailNextPrepareForTesting() noexcept { fail_next_prepare_ = true; }
	void Close()
	{
		Release(); std::exception_ptr error;
		try { op_->Close(); } catch (...) { error = std::current_exception(); }
		try { cleanup_.Check(communicator_,"distributed static close"); } catch (...) { if (!error) error = std::current_exception(); }
		if (error) std::rethrow_exception(error);
	}
protected:
	Operator& FlowOperator() noexcept { return *op_; }
	const Operator& FlowOperator() const noexcept { return *op_; }
private:
	void AssembleCandidate()
	{
		candidate_assembly_ = true;
		try { Assemble(); } catch (...) { candidate_assembly_ = false; throw; }
		candidate_assembly_ = false;
	}
	const auto& Options() const noexcept { return op_->Options(); }
	void RequireOpen() const { if (!solver_ || !committed_) throw std::logic_error("distributed static runtime is closed"); }
	void RequireStateAgreement() const
	{
		const int local = (diagnostics_.trial_active ? 1 : 0) | (diagnostics_.prepared ? 2 : 0) | (diagnostics_.converged ? 4 : 0) | (solver_ ? 8 : 0);
		int minimum = 0,maximum = 0;
		MPI_Allreduce(&local,&minimum,1,MPI_INT,MPI_MIN,communicator_); MPI_Allreduce(&local,&maximum,1,MPI_INT,MPI_MAX,communicator_);
		CollectiveLocalStage(communicator_,"distributed static state agreement",[&] { if (minimum != maximum) throw std::logic_error("inconsistent static trial state"); });
	}
	double Norm(Vec vector) const
	{
		PetscReal value = 0.0; Check("distributed static vector norm",VecNorm(vector,NORM_2,&value));
		CollectiveLocalStage(communicator_,"distributed static finite norm",[&] { if (!std::isfinite(value)) throw std::runtime_error("nonfinite static vector norm"); });
		return value;
	}
	std::array<double,3> BlockNorms() const
	{
		std::array<double,3> local{},global{};
		CollectiveLocalStage(communicator_,"distributed static residual blocks",[&] {
			PetscReadArray values; values.Acquire(op_->Assembly().Residual());
			for (PetscInt row = RowBegin(); row < RowEnd(); ++row) {
				const auto block = static_cast<std::size_t>(row) >= diagnostics_.physical_dofs ? 2 : row%4 == 3 ? 1 : 0;
				const double value = PetscRealPart(values.Data()[row-RowBegin()]); local[block] += value*value;
			}
			values.Restore();
		});
		MPI_Allreduce(local.data(),global.data(),3,MPI_DOUBLE,MPI_SUM,communicator_);
		CollectiveLocalStage(communicator_,"distributed static finite blocks",[&] {
			for (auto& value : global) {
				value = std::sqrt(value);
				if (!std::isfinite(value)) throw std::runtime_error("nonfinite static residual block norm");
			}
		});
		return global;
	}
	bool Converged(double residual,double initial,const std::array<double,3>& initial_blocks) const
	{
		if (residual > std::max(Options().nonlinear_absolute_tolerance,Options().nonlinear_relative_tolerance*initial)) return false;
		for (const auto& port : diagnostics_.ports) if (port.control_mode == ImmersedFlowPortControlMode::FlowRate)
			if (port.absolute_flow_residual_m3_s > port.flow_tolerance_m3_s) return false;
		if (Options().nonlinear_block_reduction > 0.0) {
			const auto current = BlockNorms();
			const double zero = Options().nonlinear_absolute_tolerance+Options().nonlinear_relative_tolerance*initial;
			for (std::size_t block = 0; block < current.size(); ++block)
				if (current[block] > (initial_blocks[block] > 0.0 ? initial_blocks[block]/Options().nonlinear_block_reduction : zero)) return false;
		}
		return true;
	}
	void RecordStep(const ImmersedStaticFlowNewtonStep& step)
	{
		CollectiveLocalStage(communicator_,"distributed static Newton history",[&] { diagnostics_.newton_steps.push_back(step); });
	}
	void Check(const char* stage,PetscErrorCode code) const { RequireCollectivePetscSuccess(communicator_,stage,code); }
	void Release() noexcept
	{
		candidate_assembly_ = false;
		diagnostics_.prepared = false; diagnostics_.trial_active = false; diagnostics_.converged = false;
		if (solver_) { cleanup_.Observe("static solver destroy",KSPDestroy(&solver_)); solver_ = nullptr; }
		for (auto entry : {&linear_action_,&linear_rhs_,&update_,&prepared_,&committed_})
			if (*entry) { cleanup_.Observe("static vector destroy",VecDestroy(entry)); *entry = nullptr; }
	}
	MPI_Comm communicator_;
	std::unique_ptr<Operator> op_;
	Vec committed_ = nullptr,prepared_ = nullptr,update_ = nullptr,linear_rhs_ = nullptr,linear_action_ = nullptr;
	std::unique_ptr<PetscSolverOptions> solver_options_;
	KSP solver_ = nullptr;
	ImmersedStaticFlowDiagnostics diagnostics_{};
	bool fail_next_prepare_ = false,fail_next_candidate_ = false,candidate_assembly_ = false;
	bool matrix_options_bound_ = false;
	RuntimeCleanupResult cleanup_;
};

} // namespace iga
#endif

#ifndef IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP
#define IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP

#include "ImmersedDistributedNewtonRuntime.hpp"
#include "ImmersedMovingTransientDistributedOperator.hpp"

namespace iga {

struct ImmersedMovingDistributedClock {
	double time_s=0,target_time_s=0,dt_s=0;
	std::uint64_t index=0,target_index=0;
	bool trial_active=false;
};

// Owns accepted and candidate geometry epochs. Numerical fields remain in
// distributed PETSc vectors; replicated metadata contains no global field.
// Construction, mutation and Close are collective on the borrowed communicator.
// FinalizeCommit is a local, nonthrowing publication on every member, after a
// successful collective PrepareCommit. Destruction must precede MPI_Finalize.
class ImmersedMovingTransientDistributedRuntime
{
	using Base=ImmersedDistributedNewtonRuntime<ImmersedMovingTransientDistributedOperator>;
	class Numerical : public Base {
	public:
		using Base::Base;
		using Base::FlowOperator;
	};
	struct Epoch {
		std::unique_ptr<MovingCutGeometry> geometry;
		std::unique_ptr<Numerical> flow;
		Epoch(MPI_Comm comm,std::unique_ptr<MovingCutGeometry> input,const ImmersedTransientFlowOptions& options,PetscOptions source)
			: geometry(std::move(input))
		{
			flow=AllocateCollectiveRuntime<Numerical>(comm,comm,ImmersedNewtonOptionsSource{source},"immersed_moving_",*geometry,options,ImmersedWorkPartition::FixedBackground);
		}
		const ImmersedActiveLayout& Layout() const noexcept
		{
			return flow->FlowOperator().Layout();
		}
		void Close()
		{
			flow->Close();
		}
	};
public:
	ImmersedMovingTransientDistributedRuntime(MPI_Comm comm,std::unique_ptr<MovingCutGeometry> geometry,
		const ImmersedTransientFlowOptions& options,std::uint64_t initial_index=0)
		: comm_(comm),options_(PrepareRuntimeConstructionInput<ImmersedTransientFlowOptions>(comm,"moving runtime options",options))
	{
		std::string index_identity,prefix;
		CollectiveLocalStage(comm_,"moving initial epoch preflight",[&] {
			if (!geometry) throw std::invalid_argument("moving initial geometry is absent");
			const double time=geometry->Evaluation().EvaluatedTimeS();
			if (!std::isfinite(time) || time<0) throw std::invalid_argument("invalid moving initial time");
			clock_.time_s=clock_.target_time_s=time;clock_.index=clock_.target_index=initial_index;
			index_identity=std::to_string(initial_index);
			prefix=options_.solver_options_prefix.empty()?"immersed_moving_":options_.solver_options_prefix;
		});
		RequireCollectiveSameText(comm_,"moving initial index agreement",index_identity);
		options_snapshot_=AllocateCollectiveRuntime<PetscSolverOptions>(comm_,comm_,prefix,nullptr,
			PetscOptionEntries{},std::set<std::string>{},"immersed_moving_",false);
		committed_=AllocateCollectiveRuntime<Epoch>(comm_,comm_,std::move(geometry),options_,options_snapshot_->Database());
		options_snapshot_->RecordUsed();
	}
	ImmersedMovingTransientDistributedRuntime(const ImmersedMovingTransientDistributedRuntime&)=delete;
	ImmersedMovingTransientDistributedRuntime& operator=(const ImmersedMovingTransientDistributedRuntime&)=delete;
	MPI_Comm Communicator() const noexcept
	{
		return comm_;
	}
	const ImmersedMovingDistributedClock& Clock() const noexcept
	{
		return clock_;
	}
	const MovingCutGeometry& CommittedGeometry() const
	{
		RequireLocalOpen();return *committed_->geometry;
	}
	const ImmersedActiveLayout& CommittedLayout() const
	{
		RequireLocalOpen();return committed_->Layout();
	}
	Vec CommittedState() const
	{
		RequireLocalOpen();return committed_->flow->CommittedState();
	}
	const MovingCutGeometry& TrialGeometry() const
	{
		RequireLocalTrial();return *trial_->geometry;
	}
	const ImmersedActiveLayout& TrialLayout() const
	{
		RequireLocalTrial();return trial_->Layout();
	}
	Vec TrialState() const
	{
		RequireLocalTrial();return trial_->flow->State();
	}
	const ImmersedStaticFlowDiagnostics& TrialDiagnostics() const
	{
		RequireLocalTrial();return trial_->flow->Diagnostics();
	}
	const ImmersedDistributedVelocityHistory& TrialHistory() const
	{
		RequireLocalTrial();return trial_->flow->FlowOperator().Inputs().History();
	}
	void SetCommittedOwnedState(const std::vector<PetscScalar>& values)
	{
		RequireIdle();committed_->flow->SetCommittedOwnedState(values);
	}
	void BeginTrial(std::unique_ptr<MovingCutGeometry> geometry,std::uint64_t target_index,
		std::uint32_t extension_layers=1,ImmersedVelocityExtensionOptions extension_options={})
	{
		RequireIdle();
		double target=0,dt=0;
		CollectiveLocalStage(comm_,"moving target epoch preflight",[&] {
			if (!geometry) throw std::invalid_argument("moving target geometry is absent");
			const auto& material=geometry->Evaluation();target=material.EvaluatedTimeS();dt=material.DtS();
			if (material.StepStartS()!=clock_.time_s || target!=material.StepEndS()
				|| target!=CheckedTransientTargetTime(clock_.time_s,dt)
				|| clock_.index==std::numeric_limits<std::uint64_t>::max() || target_index!=clock_.index+1)
				throw std::invalid_argument("moving target does not follow the accepted clock");
		});
		// Retired resources survive logical publication and are destroyed only
		// here (or in Close), while all members are in a collective operation.
		if (retired_)
		{
			retired_->Close();retired_.reset();
		}
		auto candidate=AllocateCollectiveRuntime<Epoch>(comm_,comm_,std::move(geometry),options_,options_snapshot_->Database());
		auto extension=AllocateCollectiveRuntime<DistributedImmersedVelocityExtension>(comm_,comm_,
			*committed_->geometry,committed_->Layout(),*candidate->geometry,candidate->Layout(),extension_layers,extension_options,options_snapshot_->Database());
		const auto seed=extension->BuildOwnedTargetSeed(committed_->flow->CommittedState(),committed_->Layout(),
			candidate->flow->State(),candidate->Layout(),clock_.time_s,clock_.index,target,target_index,dt);
		candidate->flow->SetCommittedOwnedState(seed);
		candidate->flow->FlowOperator().Freeze(*extension,committed_->flow->CommittedState(),committed_->Layout(),
			clock_.time_s,clock_.index,target,target_index,dt);
		// Frozen inputs own their values. No extension or old-geometry reference
		// is needed by Newton, preparation or local publication.
		extension->Close();extension.reset();
		options_snapshot_->RecordUsed();
		trial_=std::move(candidate);clock_.target_time_s=target;clock_.target_index=target_index;
		clock_.dt_s=dt;clock_.trial_active=true;
	}
	void Assemble()
	{
		RequireTrial();trial_->flow->Assemble();
	}
	bool SolveTrial()
	{
		RequireTrial();const bool converged=trial_->flow->SolveTrial();
		options_snapshot_->RecordUsed();return converged;
	}
	void PrepareCommit()
	{
		RequireTrial();trial_->flow->PrepareCommit();
	}
	void FinalizeCommit() noexcept
	{
		if (closed_ || !trial_ || retired_ || !trial_->flow->Diagnostics().prepared) return;
		trial_->flow->FinalizeCommit();trial_->flow->FlowOperator().ReleaseTrial();
		retired_.swap(committed_);committed_.swap(trial_);
		clock_.time_s=clock_.target_time_s;clock_.index=clock_.target_index;
		clock_.trial_active=false;clock_.dt_s=0;
	}
	void Commit()
	{
		PrepareCommit();FinalizeCommit();
	}
	void AbortPrepared() noexcept
	{
		if (trial_) trial_->flow->AbortPrepared();
	}
	// A failed Newton attempt rolls back to the mapped seed and retains frozen
	// history for retry. AbortTrial instead discards the entire candidate epoch.
	void Rollback()
	{
		RequireTrial();trial_->flow->Rollback();
	}
	void AbortTrial()
	{
		CollectiveLocalStage(comm_,"moving abort guard",[&] { RequireLocalOpen(); });
		if (!trial_) return;
		auto discarded=std::move(trial_);
		clock_.trial_active=false;clock_.dt_s=0;clock_.target_time_s=clock_.time_s;clock_.target_index=clock_.index;
		discarded->Close();
	}
	void FailNextPrepareForTesting()
	{
		RequireLocalTrial();trial_->flow->FailNextPrepareForTesting();
	}
	void FailNextCandidateForTesting()
	{
		RequireLocalTrial();trial_->flow->FailNextCandidateForTesting();
	}
	void Close()
	{
		closed_=true;clock_.trial_active=false;clock_.dt_s=0;
		std::exception_ptr error;
		for (auto* epoch:{&trial_,&retired_,&committed_}) if (*epoch) {
			try { (*epoch)->Close(); } catch (...) { if (!error) error=std::current_exception(); }
			epoch->reset();
		}
		if (options_snapshot_) {
			try { options_snapshot_->RecordUsed(); } catch (...) { if (!error) error=std::current_exception(); }
			options_snapshot_.reset();
		}
		if (error) std::rethrow_exception(error);
	}
private:
	void RequireLocalOpen() const
	{
		if (closed_ || !committed_) throw std::logic_error("moving distributed runtime is closed");
	}
	void RequireLocalTrial() const
	{
		RequireLocalOpen();if (!trial_) throw std::logic_error("moving distributed runtime has no trial");
	}
	void RequireIdle() const
	{
		CollectiveLocalStage(comm_,"moving idle guard",[&] {
			RequireLocalOpen();if (trial_) throw std::logic_error("moving distributed runtime already has a trial");
		});
	}
	void RequireTrial() const
	{
		CollectiveLocalStage(comm_,"moving trial guard",[&] {
			RequireLocalTrial();if (trial_->flow->Diagnostics().prepared) throw std::logic_error("moving distributed trial is already prepared");
		});
	}
	MPI_Comm comm_;
	ImmersedTransientFlowOptions options_;
	ImmersedMovingDistributedClock clock_;
	// Declared before epochs so their borrowed options source outlives them.
	std::unique_ptr<PetscSolverOptions> options_snapshot_;
	std::unique_ptr<Epoch> committed_,trial_,retired_;
	bool closed_=false;
};

} // namespace iga
#endif

#ifndef IGA_IMMERSED_TRANSIENT_DISTRIBUTED_RUNTIME_HPP
#define IGA_IMMERSED_TRANSIENT_DISTRIBUTED_RUNTIME_HPP

#include "ImmersedDistributedNewtonRuntime.hpp"
#include "ImmersedTransientDistributedOperator.hpp"

namespace iga {

struct ImmersedDistributedTransientClock {
	double time_s = 0,target_time_s = 0,dt_s = 0;
	std::uint64_t index = 0,target_index = 0;
	bool trial_active = false;
};

// Fixed, stationary geometry for every accepted time. Numerical fields and
// frozen history remain distributed; only scalar clock metadata is replicated.
// MPI methods are collective. FinalizeCommit and AbortPrepared are noexcept
// publication operations called on all members after coordinated preparation.
class ImmersedTransientDistributedRuntime : private ImmersedDistributedNewtonRuntime<ImmersedTransientDistributedOperator> {
	using Base = ImmersedDistributedNewtonRuntime<ImmersedTransientDistributedOperator>;
public:
	ImmersedTransientDistributedRuntime(MPI_Comm communicator,
		const CartesianDomainClassification& domain,const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,const CutCellGhostPenaltyCatalog& ghost,
		const std::string& geometry_identity,const ImmersedTransientFlowOptions& options,
		ImmersedWorkPartition partition = ImmersedWorkPartition::CellCount,
		double initial_time = 0,std::uint64_t initial_index = 0)
		: Base(communicator,"immersed_transient_",domain,volume,surface,ghost,geometry_identity,options,partition)
	{
		AgreeClock(initial_time,initial_index);
		clock_.time_s = clock_.target_time_s = initial_time;
		clock_.index = clock_.target_index = initial_index;
	}
	using Base::Communicator;
	// Diagnostics describes the numerical Newton phase; Clock owns the broader
	// frozen-trial lifetime, including a retry after a failed Newton attempt.
	using Base::Diagnostics;
	using Base::SolverConfiguration;
	using Base::RowBegin;
	using Base::RowEnd;
	using Base::State;
	using Base::CommittedState;
	using Base::OwnedStencilCount;
	using Base::RequiredStateRows;
	using Base::PortDefinitions;
	using Base::FailNextPrepareForTesting;
	using Base::FailNextCandidateForTesting;
	const ImmersedDistributedTransientClock& Clock() const noexcept { return clock_; }
	const ImmersedActiveLayout& Layout() const noexcept { return FlowOperator().Layout(); }
	const ImmersedDistributedVelocityHistory& History() const noexcept { return FlowOperator().Inputs().History(); }

	void SetCommittedOwnedState(const std::vector<PetscScalar>& values,double time,std::uint64_t index)
	{
		RequireIdle("transient committed input guard"); AgreeClock(time,index);
		Base::SetCommittedOwnedState(values);
		clock_.time_s = clock_.target_time_s = time; clock_.index = clock_.target_index = index; clock_.dt_s = 0;
	}
	void SetPortControlValue(const std::string& id,double value)
	{
		RequireIdle("transient port input guard"); Base::SetPortControlValue(id,value);
	}
	void BeginTrial(double target_time,std::uint64_t target_index,double dt)
	{
		RequireIdle("transient begin trial guard");
		std::string signature;
		CollectiveLocalStage(Communicator(),"transient target clock",[&] {
			if (!std::isfinite(dt) || !(dt > 0) || target_time != CheckedTransientTargetTime(clock_.time_s,dt)
				|| clock_.index == std::numeric_limits<std::uint64_t>::max() || target_index != clock_.index+1)
				throw std::invalid_argument("transient target time or index does not follow committed clock");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10)
				<< clock_.time_s << ':' << clock_.index << ':' << target_time << ':' << target_index << ':' << dt;
			signature = text.str();
		});
		RequireCollectiveSameText(Communicator(),"transient target clock agreement",signature);
		FlowOperator().Freeze(CommittedState(),Layout(),clock_.time_s,clock_.index,target_time,target_index,dt);
		clock_.target_time_s = target_time; clock_.target_index = target_index; clock_.dt_s = dt; clock_.trial_active = true;
	}
	void Assemble()
	{
		RequireTrial("transient assembly guard"); Base::Assemble();
	}
	bool SolveTrial()
	{
		RequireTrial("transient solve guard");
		// Base restores the committed field on failure. Frozen history/force and
		// this target clock survive, so the same trial can be retried exactly.
		return Base::SolveTrial();
	}
	void PrepareCommit()
	{
		RequireTrial("transient prepare guard"); Base::PrepareCommit();
	}
	void FinalizeCommit() noexcept
	{
		if (closed_ || !clock_.trial_active || !Diagnostics().prepared || !History().Active()) return;
		Base::FinalizeCommit();
		clock_.time_s = clock_.target_time_s; clock_.index = clock_.target_index;
		clock_.trial_active = false; clock_.dt_s = 0;
		FlowOperator().ReleaseTrial();
	}
	void Commit()
	{
		PrepareCommit(); FinalizeCommit();
	}
	void AbortPrepared() noexcept
	{
		Base::AbortPrepared();
	}
	void Rollback()
	{
		RequireOpen("transient rollback guard"); Base::Rollback();
	}
	void AbortTrial()
	{
		RequireOpen("transient abort guard"); Base::Rollback();
		FlowOperator().ReleaseTrial(); clock_.trial_active = false; clock_.dt_s = 0;
		clock_.target_time_s = clock_.time_s; clock_.target_index = clock_.index;
	}
	ImmersedStaticFlowConservationDiagnostics ConservationDiagnostics() const
	{
		RequireOpen("transient conservation guard"); return Base::ConservationDiagnostics();
	}
	void Close()
	{
		closed_ = true; clock_.trial_active = false; clock_.dt_s = 0;
		Base::Close();
	}
private:
	void AgreeClock(double time,std::uint64_t index) const
	{
		std::string signature;
		CollectiveLocalStage(Communicator(),"transient source clock",[&] {
			if (!std::isfinite(time) || time < 0) throw std::invalid_argument("invalid transient committed time");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10) << time << ':' << index; signature = text.str();
		});
		RequireCollectiveSameText(Communicator(),"transient source clock agreement",signature);
	}
	void RequireOpen(const char* stage) const
	{
		CollectiveLocalStage(Communicator(),stage,[&] { if (closed_) throw std::logic_error("transient runtime is closed"); });
	}
	void RequireIdle(const char* stage) const
	{
		CollectiveLocalStage(Communicator(),stage,[&] {
			if (closed_ || clock_.trial_active) throw std::logic_error("transient runtime is not idle and open");
		});
	}
	void RequireTrial(const char* stage) const
	{
		CollectiveLocalStage(Communicator(),stage,[&] {
			if (closed_ || !clock_.trial_active || !History().Active() || Diagnostics().prepared)
				throw std::logic_error("transient runtime requires an unprepared frozen trial");
		});
	}
	ImmersedDistributedTransientClock clock_;
	bool closed_ = false;
};

} // namespace iga
#endif

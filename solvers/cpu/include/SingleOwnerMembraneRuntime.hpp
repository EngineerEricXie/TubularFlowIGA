#ifndef IGA_SINGLE_OWNER_MEMBRANE_RUNTIME_HPP
#define IGA_SINGLE_OWNER_MEMBRANE_RUNTIME_HPP

#include "SingleOwnerSurfaceTraction.hpp"
#include "SingleOwnerSurfaceKinematics.hpp"
#include "PretensionedMembrane.hpp"
#include <memory>
#include <type_traits>

namespace iga {
// Borrowed communicator; all mutating methods are collective. Matching fluid
// and structure publication ownership is required. The numerical membrane
// exists only on owner, while every rank retains its local publication state.
class SingleOwnerMembraneRuntime {
public:
	static_assert(std::is_nothrow_swappable<SurfaceKinematics>::value && std::is_nothrow_move_assignable<SurfaceKinematics>::value,
		"membrane publication commit and cleanup must not throw");
	SingleOwnerMembraneRuntime(const SingleOwnerMembraneRuntime&)=delete;
	SingleOwnerMembraneRuntime& operator=(const SingleOwnerMembraneRuntime&)=delete;
	SingleOwnerMembraneRuntime(MPI_Comm comm,DistributedSurfaceLayout partition,int owner,
		DistributedSurfaceInterface structure,DistributedSurfaceInterface fluid,
		PretensionedMembraneMaterial material,std::vector<std::uint64_t> clamps,
		PretensionedMembraneOptions options={},PointIdentityLimits limits={})
		: comm_(comm),partition_(std::move(partition)),owner_(owner),structure_(std::move(structure)),
		fluid_(std::move(fluid)),options_(options),limits_(limits)
	{
		MPI_Comm_rank(comm_,&rank_);
		std::string config;
		CollectiveLocalStage(comm_,"single owner membrane configuration",[&] {
			ValidateDistributedSurfaceInterface(structure_);ValidateDistributedSurfaceInterface(fluid_);
			if(structure_.reference_mesh_identity_sha256!=partition_.reference_mesh_identity_sha256
				||fluid_.reference_mesh_identity_sha256!=partition_.reference_mesh_identity_sha256
				||structure_.boundary_labels!=fluid_.boundary_labels)
				throw std::invalid_argument("single owner membrane interface geometry mismatch");
			if(fluid_.provides!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::TractionOnStructure}
				||fluid_.requires!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity})
				throw std::invalid_argument("single owner membrane requires fluid traction role");
			Sha256 hash;
			for(const auto& interface:{structure_,fluid_})distributed_surface_detail::AppendString(hash,BuildDistributedSurfaceInterfaceIdentitySha256(interface));
			for(double value:{material.areal_mass_kg_per_m2,material.damping_kg_per_m2_s,material.pretension_n_per_m,material.foundation_n_per_m3})hash.AppendNormalizedDouble(value);
			hash.AppendLittleEndian64(options.maximum_nodes);hash.AppendLittleEndian64(clamps.size());
			for(auto id:clamps)hash.AppendLittleEndian64(id);
			config=hash.Hex();
		});
		RequireCollectiveSameText(comm_,"single owner membrane model agreement",config);
		layout_=GatherSurfaceLayoutAtOwner(comm_,structure_.id,partition_,owner_,options_.maximum_nodes,limits_);
		CollectiveLocalStage(comm_,"single owner membrane model construction",[&] {
			if(rank_==owner_)membrane_=std::make_unique<PretensionedMembrane>(*layout_,structure_,fluid_.id,material,std::move(clamps),options_);
		});
	}

#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
	enum class FailurePoint { None,AfterSolve,AfterPrepare };
	void SetFailureForTesting(FailurePoint point,int rank) noexcept
	{ failure_point_=point;failure_rank_=rank; }
#endif

	void SolveTrial(const FsiTrialContext& context,const SurfaceTraction& traction,
		const SurfaceFieldStamp& expected_stamp,const std::string& expected_projection)
	{
		CollectiveLocalStage(comm_,"single owner membrane trial preflight",[&] {
			ValidateFsiTrialContext(context);
			if(phase_!=Phase::Idle||context.start_time_s!=committed_time_
				||committed_step_==std::numeric_limits<std::uint64_t>::max()||context.step!=committed_step_+1
				||expected_stamp.time_s!=context.EndTime()||expected_stamp.step!=context.step
				||expected_stamp.coupling_iteration!=context.coupling_iteration)
				throw std::runtime_error("single owner membrane trial context or phase mismatch");
		});
		try {
			const auto input=GatherSurfaceTractionAtOwner(comm_,partition_,traction,fluid_.id,expected_stamp,expected_projection,owner_,options_.maximum_nodes,limits_);
			std::string identity;
			CollectiveLocalStage(comm_,"single owner membrane numerical trial",[&] {
				if(rank_!=owner_)return;
				if(BuildDistributedSurfacePartitionIdentitySha256(input->layout)!=BuildDistributedSurfacePartitionIdentitySha256(*layout_))
					throw std::runtime_error("single owner membrane numerical layout changed");
				PretensionedMembraneTrialContext numerical;
				numerical.step=context.step;numerical.start_time_s=context.start_time_s;numerical.dt_s=context.dt_s;numerical.coupling_iteration=context.coupling_iteration;
				numerical.expected_traction_stamp=input->traction.stamp;
				trial_=membrane_->SolveTrial(numerical,input->traction);
				identity=BuildSurfaceKinematicsIdentitySha256(trial_->kinematics,*layout_);
			});
#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
			InjectFailureForTesting(FailurePoint::AfterSolve);
#endif
			auto output=DistributeSingleOwnerSurfaceKinematics(comm_,partition_,owner_,layout_?&*layout_:nullptr,
				trial_?&trial_->kinematics:nullptr,identity,structure_.id,context,limits_);
			publication_=std::move(output);context_=context;phase_=Phase::Solved;
		} catch(...) { ClearTrial();throw; }
	}

	const SurfaceKinematics& TrialKinematics() const
	{
		if(phase_==Phase::Idle)throw std::runtime_error("single owner membrane has no trial publication");
		return publication_;
	}
	const SurfaceKinematics& CommittedKinematics() const
	{
		if(!has_committed_publication_)throw std::runtime_error("single owner membrane has no committed publication");
		return committed_publication_;
	}
	void AbortTrial()
	{
		CollectiveLocalStage(comm_,"single owner membrane abort",[&] { ClearTrial(); });
	}
	void PrepareCommit()
	{
		CollectiveLocalStage(comm_,"single owner membrane prepare preflight",[&] {
			if(phase_!=Phase::Solved)throw std::runtime_error("single owner membrane prepare requires solved trial");
		});
		try {
			CollectiveLocalStage(comm_,"single owner membrane prepare",[&] {
				prepared_publication_=publication_;
				if(rank_==owner_)trial_=membrane_->PrepareTrial(*trial_);
			});
#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
			InjectFailureForTesting(FailurePoint::AfterPrepare);
#endif
			phase_=Phase::Prepared;
		} catch(...) { ClearTrial();throw; }
	}
	void Commit()
	{
		CollectiveLocalStage(comm_,"single owner membrane finalize preflight",[&] {
			CoordinatorRequireFinalizeAllowed();
		});
		CoordinatorFinalizeCommitNoexcept();
	}
private:
	friend class DistributedFsiCommitCoordinator;
	void CoordinatorRequireCommitContext(MPI_Comm comm,const FsiTrialContext& context) const
	{
		int comparison=MPI_UNEQUAL;MPI_Comm_compare(comm_,comm,&comparison);
		if((comparison!=MPI_IDENT&&comparison!=MPI_CONGRUENT)||phase_!=Phase::Solved
			||context.step!=context_.step||context.start_time_s!=context_.start_time_s
			||context.dt_s!=context_.dt_s||context.coupling_iteration!=context_.coupling_iteration)
			throw std::runtime_error("single owner membrane paired commit context mismatch");
	}
	void CoordinatorRequireFinalizeAllowed() const
	{
		if(phase_!=Phase::Prepared)throw std::runtime_error("single owner membrane finalize requires preparation");
		if(rank_==owner_)membrane_->RequireFinalizeAllowed(*trial_);
	}
	void CoordinatorAbortNoexcept() noexcept { ClearTrial(); }
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		// Every fallible check completed collectively before the first mutation.
		if(rank_==owner_)membrane_->FinalizePreparedTrialNoexcept(std::move(*trial_));
		using std::swap;swap(committed_publication_,prepared_publication_);
		committed_time_=context_.EndTime();committed_step_=context_.step;has_committed_publication_=true;
		ClearTrial();
	}
#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
	void InjectFailureForTesting(FailurePoint point)
	{
		CollectiveLocalStage(comm_,"single owner membrane injected failure",[&] {
			if(failure_point_==point&&rank_==failure_rank_)
				throw std::runtime_error("injected single owner membrane failure");
		});
	}
	FailurePoint failure_point_=FailurePoint::None;
	int failure_rank_=-1;
#endif
	enum class Phase { Idle,Solved,Prepared };
	void ClearTrial() noexcept
	{
		if(membrane_)membrane_->AbortTrial();
		trial_.reset();
		publication_=SurfaceKinematics{};prepared_publication_=SurfaceKinematics{};phase_=Phase::Idle;
	}
	MPI_Comm comm_;DistributedSurfaceLayout partition_;int owner_=0,rank_=0;
	DistributedSurfaceInterface structure_,fluid_;PretensionedMembraneOptions options_;PointIdentityLimits limits_;
	std::optional<DistributedSurfaceLayout> layout_;std::unique_ptr<PretensionedMembrane> membrane_;
	std::optional<PretensionedMembraneTrial> trial_;SurfaceKinematics publication_,prepared_publication_,committed_publication_;
	FsiTrialContext context_;Phase phase_=Phase::Idle;double committed_time_=0.;std::uint64_t committed_step_=0;
	bool has_committed_publication_=false;
};
} // namespace iga
#endif

#ifndef IGA_SINGLE_OWNER_MEMBRANE_RUNTIME_HPP
#define IGA_SINGLE_OWNER_MEMBRANE_RUNTIME_HPP

#include "SingleOwnerSurfaceTraction.hpp"
#include "SingleOwnerSurfaceKinematics.hpp"
#include "MembraneCheckpoint.hpp"
#include <memory>
#include <type_traits>

namespace iga {
// Borrowed communicator; all mutating methods are collective. Matching fluid
// and structure publication ownership is required. The numerical membrane
// exists only on owner, while every rank retains its local publication state.
struct SingleOwnerMembraneCheckpoint {
	std::string metadata,numerical;
};

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
		fluid_(std::move(fluid)),material_(material),clamps_(std::move(clamps)),options_(options),limits_(limits)
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
			hash.AppendLittleEndian64(options.maximum_nodes);hash.AppendLittleEndian64(clamps_.size());
			for(auto id:clamps_)hash.AppendLittleEndian64(id);
			config=hash.Hex();
		});
		RequireCollectiveSameText(comm_,"single owner membrane model agreement",config);
		configuration_identity_=std::move(config);
		layout_=GatherSurfaceLayoutAtOwner(comm_,structure_.id,partition_,owner_,options_.maximum_nodes,limits_);
		CollectiveLocalStage(comm_,"single owner membrane model construction",[&] {
			if(rank_==owner_)membrane_=std::make_unique<PretensionedMembrane>(*layout_,structure_,fluid_.id,material_,clamps_,options_);
		});
	}

#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
	enum class FailurePoint { None,AfterSolve,AfterPrepare,BeforeRestoreFinalize };
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

	// Collective, owner-only bounded payload. The bundle must authenticate the
	// metadata hash; that metadata in turn authenticates the numerical bytes.
	std::optional<SingleOwnerMembraneCheckpoint> CaptureCheckpoint() const
	{
		std::optional<SingleOwnerMembraneCheckpoint> result;
		CollectiveLocalStage(comm_,"single owner membrane checkpoint capture",[&] {
			if(phase_!=Phase::Idle||!has_committed_publication_)
				throw std::runtime_error("membrane checkpoint requires idle accepted publication");
			if(rank_!=owner_)return;
			SingleOwnerMembraneCheckpoint payload;
			payload.numerical=SerializeMembraneCheckpoint(*membrane_,committed_step_,committed_time_,options_.maximum_nodes);
			checkpoint_metadata::Writer output;
			output.Text("IGA_SINGLE_OWNER_MEMBRANE/1");output.Text(configuration_identity_);
			output.Text(CheckpointBytesIdentity(payload.numerical));output.Text(membrane_->CommittedStateIdentitySha256());
			output.Unsigned(committed_context_.step);output.Real(committed_context_.start_time_s);
			output.Real(committed_context_.dt_s);output.Unsigned(committed_context_.coupling_iteration);
			output.Text(committed_publication_.stamp.producer_state_identity_sha256);
			payload.metadata=output.Bytes();result=std::move(payload);
		});
		return result;
	}

	// Restore only into a newly constructed runtime. All candidate allocation,
	// validation and publication redistribution precede the no-throw local swap.
	// A paired bundle should construct both fresh runtimes before publishing them.
	void RestoreCheckpoint(const SingleOwnerMembraneCheckpoint* owner_payload,
		const std::string& expected_metadata_identity)
	{
		CollectiveLocalStage(comm_,"single owner membrane restore preflight",[&] {
			if(phase_!=Phase::Idle||has_committed_publication_||committed_step_!=0
				||!IsLowercaseSha256(expected_metadata_identity)||(rank_==owner_)!=(owner_payload!=nullptr))
				throw std::runtime_error("invalid fresh membrane restore authority or phase");
		});
		RequireCollectiveSameText(comm_,"membrane restore metadata authority",expected_metadata_identity);
		std::unique_ptr<PretensionedMembrane> candidate;
		SurfaceKinematics source;std::string source_identity;FsiTrialContext accepted;
		CollectiveLocalStage(comm_,"single owner membrane restore candidate",[&] {
			if(rank_!=owner_)return;
			using checkpoint_metadata::Require;
			Require(owner_payload->metadata.size()<=checkpoint_metadata::maximum_bytes
				&&owner_payload->numerical.size()<=checkpoint_metadata::maximum_bytes,"membrane restore payload exceeds limit");
			Require(CheckpointBytesIdentity(owner_payload->metadata)==expected_metadata_identity,"membrane metadata authority differs");
			checkpoint_metadata::Reader input(owner_payload->metadata);
			Require(input.Text()=="IGA_SINGLE_OWNER_MEMBRANE/1","unsupported membrane runtime checkpoint");
			Require(input.Text()==configuration_identity_,"membrane runtime configuration differs");
			Require(input.Text()==CheckpointBytesIdentity(owner_payload->numerical),"membrane numerical payload hash differs");
			const auto state_identity=input.Text();accepted.step=input.Unsigned();accepted.start_time_s=input.Real();
			accepted.dt_s=input.Real();accepted.coupling_iteration=input.Unsigned();
			const auto producer=input.Text();input.Finish();ValidateFsiTrialContext(accepted);
			Require(accepted.step>0&&accepted.start_time_s>=0&&IsLowercaseSha256(producer),"invalid membrane accepted publication");
			auto state=ParseMembraneCheckpoint(owner_payload->numerical,*membrane_,state_identity,
				accepted.step,accepted.EndTime(),options_.maximum_nodes);
			candidate=std::make_unique<PretensionedMembrane>(*layout_,structure_,fluid_.id,material_,clamps_,options_,std::move(state));
			source.interface=structure_.id;source.stamp.time_s=accepted.EndTime();source.stamp.step=accepted.step;
			source.stamp.coupling_iteration=accepted.coupling_iteration;source.stamp.producer_state_identity_sha256=producer;
			source.stamp.reference_mesh_identity_sha256=layout_->reference_mesh_identity_sha256;
			source.stamp.layout_identity_sha256=layout_->layout_identity_sha256;
			source.stamp.partition_identity_sha256=BuildDistributedSurfacePartitionIdentitySha256(*layout_);
			const auto& restored=candidate->CommittedState();const auto& normals=candidate->ReferenceVertexNormals();
			source.displacement_m.resize(normals.size());source.velocity_m_per_s.resize(normals.size());
			for(std::size_t node=0;node<normals.size();++node)for(int axis=0;axis<3;++axis) {
				source.displacement_m[node][axis]=restored.displacement_m[node]*normals[node][axis];
				source.velocity_m_per_s[node][axis]=restored.velocity_m_per_s[node]*normals[node][axis];
			}
			source_identity=BuildSurfaceKinematicsIdentitySha256(source,*layout_);
		});
		MPI_Bcast(&accepted.step,1,MPI_UINT64_T,owner_,comm_);
		MPI_Bcast(&accepted.start_time_s,1,MPI_DOUBLE,owner_,comm_);
		MPI_Bcast(&accepted.dt_s,1,MPI_DOUBLE,owner_,comm_);
		MPI_Bcast(&accepted.coupling_iteration,1,MPI_UINT64_T,owner_,comm_);
		auto publication=DistributeSingleOwnerSurfaceKinematics(comm_,partition_,owner_,layout_?&*layout_:nullptr,
			rank_==owner_?&source:nullptr,source_identity,structure_.id,accepted,limits_);
#ifdef IGA_SINGLE_OWNER_MEMBRANE_TESTING
		InjectFailureForTesting(FailurePoint::BeforeRestoreFinalize);
#endif
		if(rank_==owner_)membrane_.swap(candidate);
		using std::swap;swap(committed_publication_,publication);
		committed_context_=accepted;committed_step_=accepted.step;committed_time_=accepted.EndTime();has_committed_publication_=true;
	}

	MPI_Comm Communicator() const noexcept { return comm_; }
	const FsiTrialContext& CommittedContext() const
	{
		if(!has_committed_publication_)throw std::logic_error("membrane committed publication is absent");
		return committed_context_;
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
	static std::string CheckpointBytesIdentity(std::string_view bytes)
	{ Sha256 hash;hash.Append(bytes.data(),bytes.size());return hash.Hex(); }
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
		committed_context_=context_;committed_time_=context_.EndTime();committed_step_=context_.step;has_committed_publication_=true;
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
	DistributedSurfaceInterface structure_,fluid_;PretensionedMembraneMaterial material_;std::vector<std::uint64_t> clamps_;
	std::string configuration_identity_;PretensionedMembraneOptions options_;PointIdentityLimits limits_;
	std::optional<DistributedSurfaceLayout> layout_;std::unique_ptr<PretensionedMembrane> membrane_;
	std::optional<PretensionedMembraneTrial> trial_;SurfaceKinematics publication_,prepared_publication_,committed_publication_;
	FsiTrialContext context_,committed_context_;Phase phase_=Phase::Idle;double committed_time_=0.;std::uint64_t committed_step_=0;
	bool has_committed_publication_=false;
};
} // namespace iga
#endif

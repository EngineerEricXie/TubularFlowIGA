#ifndef IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP
#define IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP

#include "ImmersedDistributedNewtonRuntime.hpp"
#include "ImmersedMovingTransientDistributedOperator.hpp"
#include "DistributedFluidSurfaceTraction.hpp"
#include "OwnedCheckpointVector.hpp"
#include "MovingConservationCheckpoint.hpp"

namespace iga {

struct ImmersedMovingDistributedClock {
	double time_s=0,target_time_s=0,dt_s=0;
	std::uint64_t index=0,target_index=0;
	bool trial_active=false;
};

// In-memory restore candidate. The outer bundle authenticates configuration,
// geometry and file payloads before supplying this state. Fields remain owned.
struct ImmersedMovingAcceptedCheckpoint {
	OwnedCheckpointVector field;
	std::map<std::string,double> port_control_values;
	std::string geometry_identity_sha256,publication_identity_sha256,layout_identity_sha256;
	std::string material_identity_sha256;
	std::uint64_t step=0;
	double time_s=0;
	ImmersedMovingDistributedConservation conservation;
	std::optional<MaterialSurfaceKinematics> previous_material,current_material;
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
	const std::vector<ImmersedFlowPortDefinition>& PortDefinitions() const noexcept
	{
		return options_.ports;
	}
	PetscInt CommittedRowBegin() const noexcept
	{
		return committed_?committed_->flow->RowBegin():0;
	}
	PetscInt CommittedRowEnd() const noexcept
	{
		return committed_?committed_->flow->RowEnd():0;
	}
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept
	{
		if(trial_)return trial_->flow->Diagnostics();
		if(committed_)return committed_->flow->Diagnostics();
		return closed_diagnostics_;
	}
	PetscKspConfiguration CommittedSolverConfiguration() const
	{
		RequireLocalOpen();return committed_->flow->SolverConfiguration();
	}
	std::size_t CommittedOwnedStencilCount() const
	{
		RequireLocalOpen();return committed_->flow->OwnedStencilCount();
	}
	std::size_t CommittedRequiredStateRows() const
	{
		RequireLocalOpen();return committed_->flow->RequiredStateRows();
	}
	void SetPortControlValue(const std::string& id,double value)
	{
		RequireIdle();committed_->flow->SetPortControlValue(id,value);
		for(auto& port:options_.ports)if(port.id==id)port.value=value;
		// Controls configure the next epoch; the retained transition record
		// still describes the unchanged accepted fields and geometry.
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
	ImmersedMovingAcceptedCheckpoint CaptureAcceptedCheckpoint() const
	{
		RequireIdle();ImmersedMovingAcceptedCheckpoint result;
		CollectiveLocalStage(comm_,"moving accepted checkpoint capture",[&] {
			if(!committed_conservation_||!committed_previous_material_||clock_.index==0)
				throw std::logic_error("moving checkpoint requires an accepted transition");
			result.field=CaptureOwnedCheckpointVector(committed_->flow->CommittedState());
			for(const auto& port:options_.ports)result.port_control_values.emplace(port.id,port.value);
			result.geometry_identity_sha256=committed_->geometry->GeometryIdentitySha256();
			result.publication_identity_sha256=committed_->geometry->PublicationIdentitySha256();
			result.layout_identity_sha256=committed_->Layout().HashSha256();
			result.material_identity_sha256=committed_->geometry->Evaluation().ContentIdentitySha256();
			result.step=clock_.index;result.time_s=clock_.time_s;result.conservation=*committed_conservation_;
			result.previous_material.emplace(*committed_previous_material_);result.current_material.emplace(committed_->geometry->Evaluation());
		});
		return result;
	}
	// Fresh target geometry and initial_index must already represent the saved
	// accepted epoch. Shard redistribution belongs to the outer checkpoint reader.
	void RestoreAcceptedCheckpoint(const ImmersedMovingAcceptedCheckpoint& state)
	{
		RequireIdle();std::unique_ptr<ImmersedMovingDistributedConservation> conservation;
		std::unique_ptr<MaterialSurfaceKinematics> previous_material;
		std::vector<PetscScalar> values;std::string conservation_identity;
		CollectiveLocalStage(comm_,"moving accepted checkpoint restore preflight",[&] {
			if(has_accepted_transition_||committed_conservation_||clock_.index==0||!state.previous_material||!state.current_material
				||state.step!=clock_.index||state.time_s!=clock_.time_s
				||state.geometry_identity_sha256!=committed_->geometry->GeometryIdentitySha256()
				||state.publication_identity_sha256!=committed_->geometry->PublicationIdentitySha256()
				||state.layout_identity_sha256!=committed_->Layout().HashSha256()
				||state.material_identity_sha256!=committed_->geometry->Evaluation().ContentIdentitySha256()
				||state.conservation.target_index!=state.step||state.conservation.target_time_s!=state.time_s
				||state.conservation.target_geometry_identity_sha256!=state.geometry_identity_sha256
				||state.conservation.target_publication_identity_sha256!=state.publication_identity_sha256)
				throw std::invalid_argument("moving checkpoint target epoch or identity differs");
			if(state.port_control_values.size()!=options_.ports.size())throw std::invalid_argument("moving checkpoint port catalog differs");
			for(const auto& port:options_.ports) {
				const auto found=state.port_control_values.find(port.id);
				if(found==state.port_control_values.end()||!std::isfinite(found->second)||found->second!=port.value)
					throw std::invalid_argument("moving checkpoint port control differs from fresh target");
			}
			state.previous_material->Validate();state.current_material->Validate();
			if(state.current_material->ContentIdentitySha256()!=state.material_identity_sha256
				||(!committed_->geometry->PreviousMaterialIdentitySha256().empty()
					&&state.previous_material->ContentIdentitySha256()!=committed_->geometry->PreviousMaterialIdentitySha256())
				||state.previous_material->EvaluatedTimeS()!=state.conservation.source_time_s)
				throw std::invalid_argument("moving checkpoint material history differs");
			previous_material=std::make_unique<MaterialSurfaceKinematics>(*state.previous_material);
			const auto metadata=SerializeMovingConservationCheckpoint(state.conservation);
			checkpoint_metadata::Writer controls;controls.Reals(state.port_control_values);
			Sha256 hash;distributed_surface_detail::AppendString(hash,metadata);
			distributed_surface_detail::AppendString(hash,controls.Bytes());conservation_identity=hash.Hex();
			ValidateOwnedCheckpointVector(committed_->flow->CommittedState(),state.field);
			conservation=std::make_unique<ImmersedMovingDistributedConservation>(state.conservation);
			values.assign(state.field.values.begin(),state.field.values.end());
		});
		RequireCollectiveSameText(comm_,"moving restore metadata agreement",conservation_identity);
		// The numerical setter stages owned values and completes fallible PETSc
		// copies before swapping its accepted vector. Metadata publication follows.
		committed_->flow->SetCommittedOwnedState(values);
		committed_conservation_.swap(conservation);committed_previous_material_.swap(previous_material);has_accepted_transition_=true;
	}
	void SetCommittedOwnedState(const std::vector<PetscScalar>& values)
	{
		RequireIdle();committed_->flow->SetCommittedOwnedState(values);committed_conservation_.reset();
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
	std::vector<FluidSurfaceElementState> CaptureTrialPatchState(const MaterialSurfacePatchMap& map)
	{
		CollectiveLocalStage(comm_, "moving trial patch capture guard", [&] {
			RequireLocalTrial();
			if (!trial_->flow->Diagnostics().converged || trial_->flow->Diagnostics().prepared)
				throw std::logic_error("moving patch capture requires a converged unprepared trial");
		});
		return trial_->flow->FlowOperator().CaptureOwnedPatchState(map);
	}

	// Publish only coefficients captured from this runtime's converged trial.
	// The FSI transaction supplies its context and independently expected material.
	SurfaceTraction BuildTrialSurfaceTraction(const MaterialSurfacePatchMap& map,
		const DistributedSurfaceInterface& surface, const FsiTrialContext& context,
		const std::string& expected_material, int projection_owner,
		std::size_t maximum_nodes=4096, PointIdentityLimits limits={})
	{
		CollectiveLocalStage(comm_, "moving traction transaction guard", [&] {
			RequireLocalTrial();ValidateFsiTrialContext(context);
			if (context.step!=clock_.target_index || context.start_time_s!=clock_.time_s
				|| context.EndTime()!=clock_.target_time_s
				|| expected_material!=trial_->geometry->Evaluation().ContentIdentitySha256())
				throw std::invalid_argument("moving traction transaction differs from active trial");
			// The validated context binds both endpoints. Material DtS is their
			// floating-point difference, which need not equal the nominal step.
		});
		const auto state=CaptureTrialPatchState(map);
		const auto& geometry=*trial_->geometry;
		std::vector<std::uint64_t> owned_cells;
		std::string producer;
		CollectiveLocalStage(comm_, "moving traction producer capture", [&] {
			int rank=0,size=1;MPI_Comm_rank(comm_,&rank);MPI_Comm_size(comm_,&size);
			const auto count=geometry.Domain().Background().ElementCount();
			const auto quotient=count/static_cast<std::uint64_t>(size), remainder=count%static_cast<std::uint64_t>(size);
			const auto first=quotient*rank+std::min(remainder,static_cast<std::uint64_t>(rank));
			const auto last=quotient*(rank+1)+std::min(remainder,static_cast<std::uint64_t>(rank+1));
			if(last-first>owned_cells.max_size())throw std::overflow_error("moving traction ownership exceeds storage");
			owned_cells.reserve(static_cast<std::size_t>(last-first));
			for(auto cell=first;cell<last;++cell)owned_cells.push_back(cell);
			producer=BuildFluidSurfaceTractionStateIdentitySha256(geometry.Domain(),geometry.Surface(),
				geometry.Evaluation(),map,options_.parameters.dynamic_viscosity,state);
		});
		return BuildDistributedFluidSurfaceTraction(comm_,geometry.Domain(),geometry.Surface(),geometry.Evaluation(),
			map,surface,options_.parameters.dynamic_viscosity,owned_cells,state,context,expected_material,
			producer,projection_owner,maximum_nodes,limits);
	}

	ImmersedMovingDistributedConservation ConservationDiagnostics() const
	{
		CollectiveLocalStage(comm_,"moving conservation guard",[&] {
			RequireLocalOpen();
			if(!trial_&&!committed_conservation_) throw std::logic_error("moving conservation needs a transition");
		});
		if(trial_&&!trial_->flow->Diagnostics().prepared) return BuildConservation();
		ImmersedMovingDistributedConservation result;
		CollectiveLocalStage(comm_,"moving conservation record copy",[&] {
			result=*(trial_?trial_conservation_:committed_conservation_);
		});
		return result;
	}
	void PrepareCommit()
	{
		RequireTrial();
		auto value=BuildConservation();
		std::unique_ptr<ImmersedMovingDistributedConservation> candidate;
		std::unique_ptr<MaterialSurfaceKinematics> previous_material;
		CollectiveLocalStage(comm_,"moving conservation preparation",[&] {
			candidate=std::make_unique<ImmersedMovingDistributedConservation>(std::move(value));
			previous_material=std::make_unique<MaterialSurfaceKinematics>(committed_->geometry->Evaluation());
		});
		trial_->flow->PrepareCommit();trial_conservation_.swap(candidate);trial_previous_material_.swap(previous_material);
	}
	void FinalizeCommit() noexcept
	{
		if (closed_ || !trial_ || retired_ || !trial_->flow->Diagnostics().prepared) return;
		trial_->flow->FinalizeCommit();trial_->flow->FlowOperator().ReleaseTrial();
		retired_.swap(committed_);committed_.swap(trial_);
		committed_conservation_.swap(trial_conservation_);
		committed_previous_material_.swap(trial_previous_material_);trial_previous_material_.reset();
		clock_.time_s=clock_.target_time_s;clock_.index=clock_.target_index;
		clock_.trial_active=false;clock_.dt_s=0;has_accepted_transition_=true;
	}
	void Commit()
	{
		PrepareCommit();FinalizeCommit();
	}
	void AbortPrepared() noexcept
	{
		if (trial_) { trial_->flow->AbortPrepared();trial_conservation_.reset();trial_previous_material_.reset(); }
	}
	// A failed Newton attempt rolls back to the mapped seed and retains frozen
	// history for retry. AbortTrial instead discards the entire candidate epoch.
	void Rollback()
	{
		RequireTrial();trial_->flow->Rollback();
	}
	// Local logical rollback for a collectively agreed paired-FSI failure.
	// Call on every member. PETSc resources survive until BeginTrial or Close,
	// so no allocation, numerical operation or MPI call occurs in this gate.
	void AbortTrialDeferredNoexcept() noexcept
	{
		if (closed_ || !trial_) return;
		// BeginTrial clears retired_ before installing a trial, and successful
		// FinalizeCommit removes trial_; both pointers cannot be populated here.
		retired_.swap(trial_);
		trial_conservation_.reset();trial_previous_material_.reset();
		clock_.trial_active=false;clock_.dt_s=0;
		clock_.target_time_s=clock_.time_s;clock_.target_index=clock_.index;
	}

	void AbortTrial()
	{
		CollectiveLocalStage(comm_,"moving abort guard",[&] { RequireLocalOpen(); });
		if (!trial_) return;
		auto discarded=std::move(trial_);trial_conservation_.reset();trial_previous_material_.reset();
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
		trial_conservation_.reset();trial_previous_material_.reset();committed_conservation_.reset();committed_previous_material_.reset();
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
	ImmersedMovingDistributedConservation BuildConservation() const
	{
		ImmersedMovingDistributedConservation result;
		result.endpoint=trial_->flow->FlowOperator().MaterialConservationDiagnostics();
		CollectiveLocalStage(comm_,"moving transition conservation",[&] {
			result.source_geometry_identity_sha256=committed_->geometry->GeometryIdentitySha256();
			result.target_geometry_identity_sha256=trial_->geometry->GeometryIdentitySha256();
			result.source_publication_identity_sha256=committed_->geometry->PublicationIdentitySha256();
			result.target_publication_identity_sha256=trial_->geometry->PublicationIdentitySha256();
			result.source_time_s=clock_.time_s;result.target_time_s=clock_.target_time_s;result.dt_s=clock_.dt_s;
			result.source_index=clock_.index;result.target_index=clock_.target_index;
			result.source_audited_volume_m3=committed_->geometry->Diagnostics().closed_surface_physical_volume_m3;
			result.target_audited_volume_m3=trial_->geometry->Diagnostics().closed_surface_physical_volume_m3;
			for(double volume:{result.source_audited_volume_m3,result.target_audited_volume_m3})
				if(!std::isfinite(volume)||volume<0) throw std::runtime_error("moving audited volume is invalid");
			if(!(result.dt_s>0)||!std::isfinite(result.dt_s)) throw std::logic_error("moving conservation interval is invalid");
			const double delta=result.target_audited_volume_m3-result.source_audited_volume_m3;
			result.backward_euler_volume_rate_m3_s=delta/result.dt_s;
			const auto& endpoint=result.endpoint;
			result.reynolds_defect_m3_s=result.backward_euler_volume_rate_m3_s-endpoint.total_material_surface_outward_flow_m3_s;
			result.moving_mass_defect_m3_s=result.backward_euler_volume_rate_m3_s+endpoint.total_surface_outward_flow_m3_s-endpoint.total_material_surface_outward_flow_m3_s;
			result.normalization_scale_m3_s=std::max({options_.flow_controller_reference_flow_m3_s,
				std::abs(result.backward_euler_volume_rate_m3_s),std::abs(endpoint.total_surface_outward_flow_m3_s),std::abs(endpoint.total_material_surface_outward_flow_m3_s)});
			if(!(result.normalization_scale_m3_s>0)||!std::isfinite(result.normalization_scale_m3_s))
				throw std::runtime_error("moving conservation scale is invalid");
			result.normalized_divergence_theorem_defect=std::abs(endpoint.divergence_theorem_defect_m3_s)/result.normalization_scale_m3_s;
			result.normalized_reynolds_defect=std::abs(result.reynolds_defect_m3_s)/result.normalization_scale_m3_s;
			result.normalized_moving_mass_defect=std::abs(result.moving_mass_defect_m3_s)/result.normalization_scale_m3_s;
			result.normalized_wall_relative_leakage=std::abs(endpoint.wall_relative_leakage_m3_s)/result.normalization_scale_m3_s;
			for(double value:{delta,result.backward_euler_volume_rate_m3_s,result.reynolds_defect_m3_s,result.moving_mass_defect_m3_s,
				result.normalized_divergence_theorem_defect,result.normalized_reynolds_defect,result.normalized_moving_mass_defect,result.normalized_wall_relative_leakage})
				if(!std::isfinite(value)) throw std::runtime_error("moving transition conservation is nonfinite");
		});
		return result;
	}
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
	std::unique_ptr<ImmersedMovingDistributedConservation> committed_conservation_,trial_conservation_;
	std::unique_ptr<MaterialSurfaceKinematics> committed_previous_material_,trial_previous_material_;
	ImmersedStaticFlowDiagnostics closed_diagnostics_;
	bool closed_=false,has_accepted_transition_=false;
};

} // namespace iga
#endif

#ifndef IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP
#define IGA_IMMERSED_MOVING_TRANSIENT_DISTRIBUTED_RUNTIME_HPP

#include "ImmersedDistributedNewtonRuntime.hpp"
#include "ImmersedMovingTransientDistributedOperator.hpp"
#include "DistributedFluidSurfaceTraction.hpp"

namespace iga {

struct ImmersedMovingDistributedClock {
	double time_s=0,target_time_s=0,dt_s=0;
	std::uint64_t index=0,target_index=0;
	bool trial_active=false;
};

struct ImmersedMovingDistributedConservation {
	ImmersedTransientFlowConservationDiagnostics endpoint;
	std::string source_geometry_identity_sha256,target_geometry_identity_sha256;
	std::string source_publication_identity_sha256,target_publication_identity_sha256;
	double source_time_s=0,target_time_s=0,dt_s=0;
	std::uint64_t source_index=0,target_index=0;
	double source_audited_volume_m3=0,target_audited_volume_m3=0;
	double backward_euler_volume_rate_m3_s=0,reynolds_defect_m3_s=0,moving_mass_defect_m3_s=0;
	double normalization_scale_m3_s=1;
	double normalized_divergence_theorem_defect=0,normalized_reynolds_defect=0;
	double normalized_moving_mass_defect=0,normalized_wall_relative_leakage=0;
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
		CollectiveLocalStage(comm_,"moving conservation preparation",[&] {
			candidate=std::make_unique<ImmersedMovingDistributedConservation>(std::move(value));
		});
		trial_->flow->PrepareCommit();trial_conservation_.swap(candidate);
	}
	void FinalizeCommit() noexcept
	{
		if (closed_ || !trial_ || retired_ || !trial_->flow->Diagnostics().prepared) return;
		trial_->flow->FinalizeCommit();trial_->flow->FlowOperator().ReleaseTrial();
		retired_.swap(committed_);committed_.swap(trial_);
		committed_conservation_.swap(trial_conservation_);
		clock_.time_s=clock_.target_time_s;clock_.index=clock_.target_index;
		clock_.trial_active=false;clock_.dt_s=0;
	}
	void Commit()
	{
		PrepareCommit();FinalizeCommit();
	}
	void AbortPrepared() noexcept
	{
		if (trial_) { trial_->flow->AbortPrepared();trial_conservation_.reset(); }
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
		trial_conservation_.reset();
		clock_.trial_active=false;clock_.dt_s=0;
		clock_.target_time_s=clock_.time_s;clock_.target_index=clock_.index;
	}

	void AbortTrial()
	{
		CollectiveLocalStage(comm_,"moving abort guard",[&] { RequireLocalOpen(); });
		if (!trial_) return;
		auto discarded=std::move(trial_);trial_conservation_.reset();
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
		trial_conservation_.reset();committed_conservation_.reset();
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
	ImmersedStaticFlowDiagnostics closed_diagnostics_;
	bool closed_=false;
};

} // namespace iga
#endif

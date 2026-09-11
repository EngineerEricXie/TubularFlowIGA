#ifndef IGA_IMMERSED_MOVING_DISTRIBUTED_FSI_RUNTIME_HPP
#define IGA_IMMERSED_MOVING_DISTRIBUTED_FSI_RUNTIME_HPP

#include "DistributedMaterialSurfacePatchKinematics.hpp"
#include "ImmersedMovingTransientDistributedRuntime.hpp"
#include "MovingFsiPublicationCheckpoint.hpp"
#include <functional>

namespace iga {

// The runtime and patch map are borrowed and must be mutated only through this
// adapter during a transaction. Geometry construction is local; the adapter
// coordinates failures before entering any distributed numerical operation.
class ImmersedMovingDistributedFsiRuntime
{
public:
	using GeometryProvider=std::function<std::unique_ptr<MovingCutGeometry>(
		const MovingCutGeometry&,const MaterialSurfaceKinematics&)>;
	struct Policy {
		double divergence_theorem=std::numeric_limits<double>::quiet_NaN();
		double reynolds=std::numeric_limits<double>::quiet_NaN();
		double moving_mass=std::numeric_limits<double>::quiet_NaN();
		double wall_relative_leakage=std::numeric_limits<double>::quiet_NaN();
		double discrete_continuity=std::numeric_limits<double>::quiet_NaN();
		std::uint32_t extension_layers=2;
		std::size_t maximum_patch_nodes=4096,maximum_material_vertices=1000000;
	};
	ImmersedMovingDistributedFsiRuntime(const ImmersedMovingDistributedFsiRuntime&)=delete;
	ImmersedMovingDistributedFsiRuntime& operator=(const ImmersedMovingDistributedFsiRuntime&)=delete;
	ImmersedMovingDistributedFsiRuntime(ImmersedMovingTransientDistributedRuntime& runtime,
		const MaterialSurfacePatchMap& map,const DistributedSurfaceInterface& fluid,
		const GeometryProvider& provider,int projection_owner,Policy policy,PointIdentityLimits limits={})
		: runtime_(runtime),map_(map),owner_(projection_owner),policy_(policy),limits_(limits)
	{
		std::string identity;
		CollectiveLocalStage(Communicator(),"moving FSI adapter configuration",[&] {
			if(runtime_.Clock().trial_active || !provider || !policy_.extension_layers
				|| !policy_.maximum_patch_nodes || !policy_.maximum_material_vertices)
				throw std::invalid_argument("moving FSI adapter requires idle runtime and bounded geometry policy");
			ValidateDistributedSurfaceInterface(fluid);
			if(fluid.reference_mesh_identity_sha256!=map_.Layout().reference_mesh_identity_sha256
				|| fluid.boundary_labels!=map_.Interface().boundary_labels
				|| fluid.provides!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::TractionOnStructure}
				|| fluid.requires!=std::vector<SurfaceFieldQuantity>{SurfaceFieldQuantity::Displacement,SurfaceFieldQuantity::Velocity})
				throw std::invalid_argument("moving FSI adapter fluid interface mismatch");
			int ranks=0;MPI_Comm_size(Communicator(),&ranks);
			if(owner_<0 || owner_>=ranks)throw std::invalid_argument("moving FSI projection owner is outside communicator");
			Sha256 hash;
			distributed_surface_detail::AppendString(hash,BuildDistributedSurfaceInterfaceIdentitySha256(fluid));
			for(double value:Limits()) {
				if(!std::isfinite(value) || value<0)throw std::invalid_argument("moving FSI conservation limits must be finite and nonnegative");
				hash.AppendNormalizedDouble(value);
			}
			hash.AppendLittleEndian64(policy_.extension_layers);
			hash.AppendLittleEndian64(policy_.maximum_patch_nodes);hash.AppendLittleEndian64(policy_.maximum_material_vertices);
			distributed_surface_detail::AppendString(hash,map_.ReferenceIdentitySha256());
			configuration_identity_=hash.Hex();hash.AppendLittleEndian64(owner_);
			identity=hash.Hex();fluid_=fluid;provider_=provider;
		});
		RequireCollectiveSameText(Communicator(),"moving FSI adapter policy agreement",identity);
	}
	MPI_Comm Communicator() const noexcept { return runtime_.Communicator(); }
	bool UsesFlowRuntime(const ImmersedMovingTransientDistributedRuntime& flow) const noexcept { return &runtime_==&flow; }
	void SolveTrial(const FsiTrialContext& context,const SurfaceKinematics& kinematics,
		const SurfaceFieldStamp& expected_stamp,const std::string& expected_committed_material)
	{
		CollectiveLocalStage(Communicator(),"moving FSI trial preflight",[&] {
			ValidateFsiTrialContext(context);
			(void)runtime_.CommittedGeometry();
			const auto& clock=runtime_.Clock();
			if(trial_ || clock.trial_active || context.start_time_s!=clock.time_s
				|| clock.index==std::numeric_limits<std::uint64_t>::max() || context.step!=clock.index+1)
				throw std::runtime_error("moving FSI trial context or phase mismatch");
		});
		try {
			const auto composition=ComposeDistributedMaterialSurfacePatch(Communicator(),map_,
				runtime_.CommittedGeometry().Evaluation(),kinematics,expected_stamp,expected_committed_material,
				context,policy_.maximum_patch_nodes,policy_.maximum_material_vertices,limits_);
			std::unique_ptr<MovingCutGeometry> geometry;
			CollectiveLocalStage(Communicator(),"moving FSI target construction",[&] {
				geometry=provider_(runtime_.CommittedGeometry(),composition.target);
				if(!geometry || geometry->Evaluation().ContentIdentitySha256()!=composition.target.ContentIdentitySha256())
					throw std::runtime_error("moving FSI provider changed composed material epoch");
			});
			runtime_.BeginTrial(std::move(geometry),context.step,policy_.extension_layers);
			const bool solved=runtime_.SolveTrial();
			CollectiveLocalStage(Communicator(),"moving FSI convergence",[&] {
				if(!solved)throw std::runtime_error("moving FSI trial did not converge");
			});
			CheckConservation();
			auto traction=runtime_.BuildTrialSurfaceTraction(map_,fluid_,context,
				composition.target.ContentIdentitySha256(),owner_,policy_.maximum_patch_nodes,limits_);
			std::unique_ptr<Publication> candidate;
			CollectiveLocalStage(Communicator(),"moving FSI publication preparation",[&] {
				candidate=std::make_unique<Publication>(Publication{context,std::move(traction),composition.composition_identity_sha256});
			});
			trial_.swap(candidate);
		} catch(...) { CoordinatorAbortNoexcept();throw; }
	}
	const SurfaceTraction& TrialTraction() const
	{
		if(!trial_)throw std::logic_error("moving FSI trial publication is absent");
		return trial_->traction;
	}
	const SurfaceTraction& CommittedTraction() const
	{
		if(!committed_)throw std::logic_error("moving FSI committed publication is absent");
		return committed_->traction;
	}
	const FsiTrialContext& CommittedContext() const
	{
		if(!committed_)throw std::logic_error("moving FSI committed publication is absent");
		return committed_->context;
	}
	std::string CaptureCheckpoint() const
	{
		std::string bytes;
		CollectiveLocalStage(Communicator(),"moving FSI checkpoint capture",[&] {
			if(trial_||prepared_||runtime_.Clock().trial_active||!committed_)
				throw std::runtime_error("FSI checkpoint requires idle accepted publication");
			const auto& geometry=runtime_.CommittedGeometry();
			if(committed_->context.step!=runtime_.Clock().index||committed_->context.EndTime()!=runtime_.Clock().time_s)
				throw std::runtime_error("FSI publication differs from accepted fluid clock");
			bytes=SerializeMovingFsiPublicationCheckpoint({committed_->context,committed_->traction,committed_->composition_identity,
				geometry.Evaluation().ContentIdentitySha256(),geometry.GeometryIdentitySha256()},map_.Layout(),configuration_identity_,policy_.maximum_patch_nodes);
		});
		return bytes;
	}
	// The enclosing bundle supplies the authenticated local payload identity.
	// Restore into a fresh adapter over an already restored, unpublished fluid
	// candidate. Caller publishes the fluid/structure owners only as a pair.
	void RestoreCheckpoint(std::string_view bytes,const std::string& expected_payload_identity)
	{
		std::unique_ptr<Publication> candidate;std::string common;
		CollectiveLocalStage(Communicator(),"moving FSI restore phase",[&] {
			if(trial_||prepared_||committed_||runtime_.Clock().trial_active)
				throw std::runtime_error("FSI restore requires fresh idle adapter");
		});
		const auto conservation=runtime_.ConservationDiagnostics();
		CollectiveLocalStage(Communicator(),"moving FSI publication restore",[&] {
			const auto& geometry=runtime_.CommittedGeometry();
			auto state=ParseMovingFsiPublicationCheckpoint(bytes,expected_payload_identity,map_.Layout(),fluid_.id,configuration_identity_,
				geometry.Evaluation().ContentIdentitySha256(),geometry.GeometryIdentitySha256(),runtime_.Clock().index,runtime_.Clock().time_s,policy_.maximum_patch_nodes);
			if(state.context.start_time_s!=conservation.source_time_s||state.context.step!=conservation.target_index)
				throw std::runtime_error("FSI restored context differs from fluid history");
			Sha256 hash;
			for(const auto* id:{&state.composition_identity,&state.material_identity,&state.geometry_identity})distributed_surface_detail::AppendString(hash,*id);
			hash.AppendLittleEndian64(state.context.step);hash.AppendNormalizedDouble(state.context.start_time_s);
			hash.AppendNormalizedDouble(state.context.dt_s);hash.AppendLittleEndian64(state.context.coupling_iteration);common=hash.Hex();
			candidate=std::make_unique<Publication>(Publication{state.context,std::move(state.traction),std::move(state.composition_identity)});
		});
		RequireCollectiveSameText(Communicator(),"moving FSI restored publication agreement",common);
		committed_.swap(candidate);
	}
	void AbortTrial()
	{
		CollectiveLocalStage(Communicator(),"moving FSI abort",[&] { CoordinatorAbortNoexcept(); });
	}
	void PrepareCommit()
	{
		CollectiveLocalStage(Communicator(),"moving FSI prepare preflight",[&] {
			if(!trial_ || prepared_)throw std::logic_error("moving FSI prepare requires solved trial");
		});
		CheckConservation();runtime_.PrepareCommit();prepared_=true;
	}
private:
	friend class DistributedFsiCommitCoordinator;
	struct Publication { FsiTrialContext context;SurfaceTraction traction;std::string composition_identity; };
	std::array<double,5> Limits() const noexcept
	{
		return {{policy_.divergence_theorem,policy_.reynolds,policy_.moving_mass,
			policy_.wall_relative_leakage,policy_.discrete_continuity}};
	}
	void CheckConservation() const
	{
		const auto value=runtime_.ConservationDiagnostics();
		CollectiveLocalStage(Communicator(),"moving FSI physical conservation",[&] {
			const std::array<double,5> actual{{value.normalized_divergence_theorem_defect,value.normalized_reynolds_defect,
				value.normalized_moving_mass_defect,value.normalized_wall_relative_leakage,
				value.endpoint.normalized_discrete_moving_wall_continuity_defect}};
			const auto limits=Limits();
			for(std::size_t i=0;i<actual.size();++i)
				if(!std::isfinite(actual[i]) || actual[i]>limits[i]) {
					std::ostringstream message;message<<std::setprecision(17)
						<<"moving FSI conservation gate "<<i<<" failed: "<<actual[i]<<" > "<<limits[i]
						<<" normalization_m3_s="<<value.normalization_scale_m3_s
						<<" divergence_defect_m3_s="<<value.endpoint.divergence_theorem_defect_m3_s
						<<" moving_mass_defect_m3_s="<<value.moving_mass_defect_m3_s
						<<" wall_relative_leakage_m3_s="<<value.endpoint.wall_relative_leakage_m3_s;
					throw std::runtime_error(message.str());
				}
		});
	}
	void CoordinatorRequireCommitContext(MPI_Comm comm,const FsiTrialContext& context) const
	{
		int comparison=MPI_UNEQUAL;MPI_Comm_compare(Communicator(),comm,&comparison);
		if((comparison!=MPI_IDENT && comparison!=MPI_CONGRUENT) || !trial_ || prepared_
			|| context.step!=trial_->context.step || context.start_time_s!=trial_->context.start_time_s
			|| context.dt_s!=trial_->context.dt_s || context.coupling_iteration!=trial_->context.coupling_iteration)
			throw std::runtime_error("moving FSI paired commit context mismatch");
	}
	void CoordinatorRequireFinalizeAllowed() const
	{
		if(!trial_ || !prepared_ || !runtime_.Clock().trial_active || !runtime_.TrialDiagnostics().prepared)
			throw std::runtime_error("moving FSI finalize requires prepared trial");
	}
	void CoordinatorAbortNoexcept() noexcept
	{
		runtime_.AbortTrialDeferredNoexcept();trial_.reset();prepared_=false;
	}
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		runtime_.FinalizeCommit();committed_.swap(trial_);trial_.reset();prepared_=false;
	}
	ImmersedMovingTransientDistributedRuntime& runtime_;
	const MaterialSurfacePatchMap& map_;
	DistributedSurfaceInterface fluid_;GeometryProvider provider_;
	std::string configuration_identity_;
	int owner_;Policy policy_;PointIdentityLimits limits_;
	std::unique_ptr<Publication> trial_,committed_;
	bool prepared_=false;
};

} // namespace iga
#endif

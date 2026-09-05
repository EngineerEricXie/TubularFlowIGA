#ifndef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_FSI_RUNTIME_HPP
#define IGA_MOVING_IMMERSED_TRANSIENT_FLOW_FSI_RUNTIME_HPP

// Rank-local fluid-side FSI adapter.  This owns the moving epoch transaction
// and the surface-publication transaction together, but deliberately owns no
// transport/coordinator and exposes no mutable numerical owner alias.
#include "FluidSurfaceTraction.hpp"
#include "FsiDomainRuntime.hpp"
#include "MaterialSurfacePatchKinematics.hpp"
#include "MovingImmersedTransientFlowRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

class StrongFluidStructureCouplingAccess;

class MovingImmersedTransientFlowFsiRuntime final : public FsiFluidDomainRuntime {
public:
	static_assert(std::is_nothrow_swappable<std::optional<SurfaceTraction>>::value,
		"solved traction publication must swap without throwing");
	static_assert(std::is_nothrow_swappable<std::unique_ptr<MaterialSurfaceKinematics>>::value,
		"solved material publication must swap without throwing");
	static_assert(std::is_nothrow_swappable<std::optional<ImmersedTransientFlowDiagnostics>>::value,
		"solved diagnostics publication must swap without throwing");
	MovingImmersedTransientFlowFsiRuntime(const MovingImmersedTransientFlowFsiRuntime&) = delete;
	MovingImmersedTransientFlowFsiRuntime& operator=(const MovingImmersedTransientFlowFsiRuntime&) = delete;
	MovingImmersedTransientFlowFsiRuntime(MovingImmersedTransientFlowFsiRuntime&&) = delete;
	MovingImmersedTransientFlowFsiRuntime& operator=(MovingImmersedTransientFlowFsiRuntime&&) = delete;

	MovingImmersedTransientFlowFsiRuntime(std::string domain_id, std::string subsystem_id,
		FsiCouplingEdge edge, DistributedSurfaceInterface fluid_surface,
		DistributedSurfaceInterface structure_surface, DistributedSurfaceLayout fluid_layout,
		DistributedSurfaceLayout structure_layout, MaterialSurfaceKinematics initial_full,
		MaterialSurfacePatchMap patch_map, MovingImmersedTransientFlowOptions moving_options,
		FluidSurfaceTractionProjectionOptions traction_options = {})
		: domain_id_(std::move(domain_id)), subsystem_id_(std::move(subsystem_id)), edge_(std::move(edge)),
		  catalog_{std::move(fluid_surface)}, structure_surface_(std::move(structure_surface)),
		  fluid_layout_(std::move(fluid_layout)), structure_layout_(std::move(structure_layout)),
		  patch_map_(std::move(patch_map)), committed_full_(new MaterialSurfaceKinematics(initial_full)),
		  moving_options_viscosity_(moving_options.flow.parameters.dynamic_viscosity),
		  moving_(*committed_full_, std::move(moving_options)), traction_options_(traction_options),
		  lifecycle_(domain_id_, subsystem_id_, edge_, catalog_.front(), structure_surface_,
			fluid_layout_, structure_layout_)
	{
		ValidateConstruction();
	}

	const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept override
	{ return catalog_; }
	const FsiTrialLifecycle& Lifecycle() const noexcept { return lifecycle_; }
	const FsiCouplingEdge& StrongCouplingEdge() const noexcept { return edge_; }
	const DistributedSurfaceLayout& StrongCouplingLayout() const noexcept { return fluid_layout_; }
	const MovingImmersedTransientFlowDiagnostics& MovingDiagnostics() const noexcept { return moving_.Diagnostics(); }
	const MovingCutGeometry& CommittedGeometry() const noexcept { return moving_.CommittedGeometry(); }
	const ImmersedGlobalFlowState& CommittedGlobalState() const { return moving_.CommittedGlobalState(); }
	ImmersedTransientFlowDiagnostics TrialFlowDiagnostics() const
	{
		if (!trial_flow_diagnostics_.has_value() || !lifecycle_.HasTrialOutput())
			throw std::logic_error("fluid FSI trial flow diagnostics are unavailable before a solved trial");
		return *trial_flow_diagnostics_;
	}
	std::string TrialCompositionIdentitySha256() const
	{
		if (trial_composition_identity_sha256_.empty() || !lifecycle_.HasTrialOutput())
			throw std::logic_error("fluid FSI composition identity is unavailable before a solved trial");
		return trial_composition_identity_sha256_;
	}
	std::string CommittedCompositionIdentitySha256() const
	{
		if (committed_composition_identity_sha256_.empty() || !lifecycle_.HasCommittedOutput())
			throw std::logic_error("fluid FSI composition identity is unavailable before commit");
		return committed_composition_identity_sha256_;
	}
	MaterialSurfaceKinematics TrialMaterialKinematicsSnapshot() const
	{
		if (!trial_target_ || !lifecycle_.HasTrialOutput())
			throw std::logic_error("fluid FSI target material state is unavailable before a solved trial");
		return *trial_target_;
	}
	MaterialSurfaceKinematics CommittedMaterialKinematicsSnapshot() const { return *committed_full_; }
	std::optional<SurfaceTraction> CommittedSurfaceTractionSnapshot() const
	{ return committed_traction_; }
	FluidSurfaceTractionDiagnostics TrialSurfaceTractionDiagnostics() const
	{
		if (!trial_traction_diagnostics_.has_value() || !lifecycle_.HasTrialOutput())
			throw std::logic_error("fluid FSI traction diagnostics are unavailable before a solved trial");
		return *trial_traction_diagnostics_;
	}
	FluidSurfaceTractionDiagnostics CommittedSurfaceTractionDiagnostics() const
	{
		if (!committed_traction_diagnostics_.has_value() || !lifecycle_.HasCommittedOutput())
			throw std::logic_error("fluid FSI traction diagnostics are unavailable before commit");
		return *committed_traction_diagnostics_;
	}
	MovingImmersedTransientFlowConservationDiagnostics ConservationDiagnostics() const
	{
		// A target-endpoint measurement becomes externally meaningful only with
		// the exact solved traction that describes that same trial.  Once the
		// trial is gone, retain only the committed transition measurement.
		if (lifecycle_.HasTrialOutput()) return moving_.ConservationDiagnostics();
		if (lifecycle_.HasCommittedOutput() && moving_.Idle()) return moving_.ConservationDiagnostics();
		throw std::logic_error("fluid FSI conservation is unavailable before a solved trial or committed transition");
	}

	void BeginMacroStep(const DomainStepContext& step)
	{
		step.Validate();
		if (!moving_.Idle()) throw std::runtime_error("fluid FSI runtime retains an unexpected moving trial");
		if (lifecycle_.HasActiveStep()) throw std::runtime_error("fluid FSI runtime already has an active macro step");
		const auto& global = moving_.CommittedGlobalState();
		if (global.Index() == std::numeric_limits<std::uint64_t>::max()
			|| static_cast<std::uint64_t>(step.step_index) != global.Index()+1
			|| step.start_time_s != global.TimeS()
			|| committed_full_->EvaluatedTimeS() != step.start_time_s
			|| committed_full_->StepEndS() != step.start_time_s
			|| moving_.CommittedGeometry().Kinematics().ContentIdentitySha256() != committed_full_->ContentIdentitySha256())
			throw std::runtime_error("fluid FSI macro step is not the exact next committed moving epoch");
		lifecycle_.BeginStep(step);
	}
	void BeginStep(const DomainStepContext& step) { BeginMacroStep(step); }

	void BeginCouplingIteration(std::uint64_t iteration,
		const SurfaceFieldStamp& expected_kinematics, const SurfaceFieldStampEnvelope& traction_envelope)
	{
		if (!moving_.Idle() || input_.has_value() || trial_traction_.has_value())
			throw std::runtime_error("fluid FSI runtime retains an unexpected trial publication");
		SurfaceFieldStamp input_copy = expected_kinematics;
		SurfaceFieldStampEnvelope output_copy = traction_envelope;
		lifecycle_.BeginIteration(iteration, input_copy, output_copy);
		using std::swap;
		swap(expected_kinematics_, input_copy);
	}
	void BeginIteration(std::uint64_t iteration, const SurfaceFieldStamp& expected_kinematics,
		const SurfaceFieldStampEnvelope& traction_envelope)
	{ BeginCouplingIteration(iteration, expected_kinematics, traction_envelope); }

	void SetSurfaceKinematics(const std::string& interface_id,
		const SurfaceKinematics& kinematics) override
	{
		if (interface_id != edge_.fluid.interface_id)
			throw std::runtime_error("unknown fluid FSI surface interface");
		// Do this before constructing or replacing any snapshot.  In particular,
		// a duplicate submission while InputReady must preserve the accepted input
		// and its ability to drive the pending solve.
		if (lifecycle_.Phase() != FsiTrialPhase::IterationAwaitingInput || input_.has_value())
			throw std::runtime_error("fluid FSI kinematics are unavailable outside the exact awaiting iteration");
		SurfaceKinematics candidate = kinematics;
		ValidateFsiFluidKinematicsInput(edge_, candidate, structure_layout_, expected_kinematics_);
		// Complete all fallible ownership work locally.  MarkInput validates before
		// it changes its phase; after it succeeds the optional exchange is noexcept
		// and makes the input snapshot/lifecycle gate visible together.
		std::optional<SurfaceKinematics> accepted;
		accepted.emplace(std::move(candidate));
		lifecycle_.MarkInput(accepted->interface, accepted->stamp);
		input_.swap(accepted);
	}

	void SolveFluidTrial()
	{
		lifecycle_.RequireSolveAllowed();
		if (!input_.has_value()) throw std::runtime_error("fluid FSI runtime has no accepted kinematics snapshot");
		if (!moving_.Idle()) throw std::runtime_error("fluid FSI runtime already owns a moving target epoch");
		const auto composed = MaterialSurfacePatchKinematics::ComposeTarget(patch_map_, *committed_full_,
			*input_, lifecycle_.Context());
		try {
			moving_.BeginTrial(composed.target, lifecycle_.Context().step, lifecycle_.Context().dt_s);
			moving_.Assemble();
			if (!moving_.SolveTrial()) return; // preserve diagnostics; caller must reject or abort.
			const auto state = GatherRequiredTrialSurfaceState();
			SurfaceFieldStamp actual = MakeActualOutputStamp(state);
			auto result = BuildFluidSurfaceTraction(moving_.TrialGeometry().Domain(), moving_.TrialGeometry().Surface(),
				composed.target, catalog_.front(), fluid_layout_, patch_map_, actual,
				DynamicViscosity(), state, traction_options_);
			// Complete every allocating copy before the lifecycle becomes observable.
			// Once MarkSolved succeeds, publication is only noexcept exchanges and
			// scalar state, so no half-published output can escape this adapter.
			std::optional<SurfaceTraction> staged_traction;
			std::optional<FluidSurfaceTractionDiagnostics> staged_traction_diagnostics;
			std::unique_ptr<MaterialSurfaceKinematics> staged_target(new MaterialSurfaceKinematics(composed.target));
			std::optional<ImmersedTransientFlowDiagnostics> staged_flow_diagnostics;
			std::string staged_composition_identity_sha256;
			staged_traction.emplace(std::move(result.traction));
			staged_traction_diagnostics.emplace(result.diagnostics);
			staged_flow_diagnostics.emplace(moving_.TrialDiagnostics());
			staged_composition_identity_sha256 = composed.composition_identity_sha256;
			if (staged_composition_identity_sha256.empty())
				throw std::logic_error("fluid FSI composition identity is absent");
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
			if (fail_next_late_publication_) {
				fail_next_late_publication_ = false;
				throw std::runtime_error("injected fluid FSI late-publication failure");
			}
#endif
			// This is the last fallible publication. The staged payload is complete
			// before the lifecycle can advertise any solved output.
			lifecycle_.MarkSolved(staged_traction->interface, staged_traction->stamp);
			trial_traction_.swap(staged_traction);
			trial_traction_diagnostics_.swap(staged_traction_diagnostics);
			trial_target_.swap(staged_target);
			trial_flow_diagnostics_.swap(staged_flow_diagnostics);
			trial_composition_identity_sha256_.swap(staged_composition_identity_sha256);
		}
		catch (...) {
			moving_.AbortTrial();
			trial_traction_.reset(); trial_traction_diagnostics_.reset(); trial_target_.reset(); trial_flow_diagnostics_.reset();
			trial_composition_identity_sha256_.clear();
			throw;
		}
	}
	void SolveTrial() { SolveFluidTrial(); }

	SurfaceTraction GetSurfaceTraction(const std::string& interface_id) const override
	{
		if (interface_id != edge_.fluid.interface_id)
			throw std::runtime_error("unknown fluid FSI surface interface");
		if (!trial_traction_.has_value()) throw std::runtime_error("fluid FSI runtime has no solved traction");
		lifecycle_.RequireTrialOutput(trial_traction_->interface, trial_traction_->stamp);
		return *trial_traction_;
	}

	SurfaceTraction GetCommittedSurfaceTraction(const std::string& interface_id,
		const SurfaceFieldStamp& stamp) const
	{
		if (interface_id != edge_.fluid.interface_id)
			throw std::runtime_error("unknown fluid FSI surface interface");
		if (!committed_traction_.has_value()) throw std::runtime_error("fluid FSI traction has not been committed");
		lifecycle_.RequireCommittedOutput(catalog_.front().id, stamp);
		return *committed_traction_;
	}

	void RejectCouplingIteration()
	{
		if (!lifecycle_.HasActiveStep()) throw std::runtime_error("fluid FSI reject requires an active step");
		DiscardTrialNoexcept();
		lifecycle_.RejectIteration();
	}
	void RollbackTrial() { RejectCouplingIteration(); }

	void PrepareCommitStep()
	{
		if (!trial_traction_.has_value() || !trial_traction_diagnostics_.has_value() || !trial_target_ || moving_.Idle() || moving_.Prepared())
			throw std::runtime_error("fluid FSI prepare requires a solved moving trial and traction");
		lifecycle_.RequireTrialOutput(trial_traction_->interface, trial_traction_->stamp);
		// All allocations/copies and both potentially throwing preparations occur
		// before either prepared state is made visible to a caller.
		SurfaceTraction traction_copy = *trial_traction_;
		FluidSurfaceTractionDiagnostics diagnostics_copy = *trial_traction_diagnostics_;
		std::unique_ptr<MaterialSurfaceKinematics> target_copy(new MaterialSurfaceKinematics(*trial_target_));
		try {
			moving_.PrepareCommit();
			lifecycle_.PrepareCommit();
			prepared_traction_.emplace(std::move(traction_copy));
			prepared_traction_diagnostics_.emplace(diagnostics_copy);
			prepared_target_.swap(target_copy);
		}
		catch (...) {
			moving_.AbortPrepared();
			prepared_traction_.reset(); prepared_traction_diagnostics_.reset(); prepared_target_.reset();
			trial_traction_.reset(); trial_traction_diagnostics_.reset(); trial_target_.reset(); trial_flow_diagnostics_.reset(); input_.reset();
			trial_composition_identity_sha256_.clear();
			if (lifecycle_.HasActiveStep()) lifecycle_.RejectIteration();
			throw;
		}
	}
	void PrepareCommit() { PrepareCommitStep(); }

	void FinalizeCommitStep()
	{
		if (!prepared_traction_.has_value() || !prepared_target_)
			throw std::runtime_error("fluid FSI finalize requires prepared publication snapshots");
		if (!moving_.Prepared()) throw std::runtime_error("fluid FSI moving trial is not prepared");
		lifecycle_.RequireFinalizeAllowed();
		// Preconditions above complete every throwing check.  The coordinated
		// ownership handoff below consists only of noexcept finalizes/swaps.
		moving_.FinalizeCommit();
		using std::swap;
		swap(committed_traction_, prepared_traction_);
		swap(committed_traction_diagnostics_, prepared_traction_diagnostics_);
		committed_full_.swap(prepared_target_);
		committed_composition_identity_sha256_.swap(trial_composition_identity_sha256_);
		lifecycle_.FinalizePreparedCommitNoexcept();
		trial_traction_.reset(); trial_traction_diagnostics_.reset(); trial_target_.reset(); trial_flow_diagnostics_.reset(); input_.reset();
		prepared_traction_.reset(); prepared_traction_diagnostics_.reset(); prepared_target_.reset();
	}
	void FinalizeCommit() { FinalizeCommitStep(); }

	void AbortStep() noexcept
	{
		DiscardTrialNoexcept();
		lifecycle_.AbortStep();
	}

	// These are read only while the prepared snapshots remain abortable.  They
	// let a strong coordinator stage its complete result before finalization.
	std::string CoordinatorPreparedCommittedTractionIdentitySha256() const
	{
		if (!prepared_traction_.has_value()) throw std::runtime_error("fluid FSI prepared traction is unavailable");
		return BuildSurfaceTractionIdentitySha256(*prepared_traction_, fluid_layout_);
	}
	std::string CoordinatorPreparedCommittedCompositionIdentitySha256() const
	{
		if (!prepared_target_ || trial_composition_identity_sha256_.empty())
			throw std::runtime_error("fluid FSI prepared composition is unavailable");
		return trial_composition_identity_sha256_;
	}

#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	void FailNextPrepareForTesting() { moving_.FailNextPrepareForTesting(); }
	void FailNextLatePublicationForTesting() noexcept { fail_next_late_publication_ = true; }
	double MaxTrialSurfaceRelativeVelocityNormForTesting() const
	{ return moving_.MaxTrialSurfaceRelativeVelocityNormForTesting(); }
#endif

private:
	friend class StrongFluidStructureCouplingAccess;
	void CoordinatorRequireFinalizeAllowed() const
	{
		if (!prepared_traction_.has_value() || !prepared_target_ || !moving_.Prepared())
			throw std::runtime_error("fluid FSI coordinator finalize requires prepared state");
		lifecycle_.RequireFinalizeAllowed();
	}
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		moving_.FinalizeCommit(); using std::swap;
		swap(committed_traction_, prepared_traction_); swap(committed_traction_diagnostics_, prepared_traction_diagnostics_);
		committed_full_.swap(prepared_target_); committed_composition_identity_sha256_.swap(trial_composition_identity_sha256_);
		lifecycle_.FinalizePreparedCommitNoexcept(); trial_traction_.reset(); trial_traction_diagnostics_.reset(); trial_target_.reset(); trial_flow_diagnostics_.reset(); input_.reset(); prepared_traction_.reset(); prepared_traction_diagnostics_.reset(); prepared_target_.reset();
	}
	void ValidateConstruction() const
	{
		ValidateFsiCouplingEdge(edge_);
		ValidateFsiFluidSurfaceInterface(catalog_.front());
		ValidateFsiStructureSurfaceInterface(structure_surface_);
		ValidateDistributedSurfaceLayout(fluid_layout_);
		ValidateDistributedSurfaceLayout(structure_layout_);
		committed_full_->Validate(); patch_map_.FullReference().Validate();
		if (!(catalog_.front().id == edge_.fluid) || !(structure_surface_.id == edge_.structure)
			|| catalog_.front().id.domain_id != domain_id_ || catalog_.front().id.subsystem_id != subsystem_id_)
			throw std::runtime_error("fluid FSI runtime surfaces do not bind the exact local edge endpoint");
		if (!(patch_map_.Interface().id == edge_.structure)
			|| patch_map_.Layout().layout_identity_sha256 != structure_layout_.layout_identity_sha256
			|| BuildDistributedSurfacePartitionIdentitySha256(patch_map_.Layout())
				!= BuildDistributedSurfacePartitionIdentitySha256(structure_layout_)
			|| fluid_layout_.layout_identity_sha256 != patch_map_.Layout().layout_identity_sha256
			|| BuildDistributedSurfacePartitionIdentitySha256(fluid_layout_)
				!= BuildDistributedSurfacePartitionIdentitySha256(patch_map_.Layout()))
			throw std::runtime_error("fluid FSI runtime layouts do not bind the immutable material patch map");
		if (committed_full_->MaterialIdentitySha256() != patch_map_.FullReference().MaterialIdentitySha256()
			|| committed_full_->TopologyIdentitySha256() != patch_map_.FullReference().TopologyIdentitySha256()
			|| moving_.CommittedGeometry().Kinematics().ContentIdentitySha256() != committed_full_->ContentIdentitySha256())
			throw std::runtime_error("fluid FSI runtime initial full material state does not bind the patch map/moving geometry");
		if (fluid_layout_.partition_count != 1 || fluid_layout_.partition_rank != 0
			|| structure_layout_.partition_count != 1 || structure_layout_.partition_rank != 0)
			throw std::runtime_error("fluid FSI runtime is bounded to one complete rank-local partition");
		if (!std::isfinite(traction_options_.conservation_absolute_force_tolerance_n)
			|| !std::isfinite(traction_options_.conservation_absolute_moment_tolerance_n_m)
			|| !std::isfinite(traction_options_.conservation_relative_tolerance))
			throw std::runtime_error("fluid FSI traction projection options are invalid");
	}

	double DynamicViscosity() const
	{
		const double result = moving_options_viscosity_;
		if (!(result > 0.0) || !std::isfinite(result))
			throw std::runtime_error("fluid FSI moving-flow dynamic viscosity is invalid");
		return result;
	}

	std::vector<FluidSurfaceElementState> GatherRequiredTrialSurfaceState() const
	{
		const auto& geometry = moving_.TrialGeometry(); const auto& domain = geometry.Domain();
		const auto& catalog = geometry.Surface(); const auto& layout = moving_.TrialLayout();
		const auto values = moving_.TrialState();
		if (values.size() < 4*layout.NodeIds().size())
			throw std::runtime_error("fluid FSI trial vector does not cover active IGA nodes");
		std::vector<bool> selected(geometry.Kinematics().CanonicalTriangleProvenance().size(), false);
		for (std::size_t i = 0; i < patch_map_.LayoutTriangleToSourceTriangles().size(); ++i)
			selected[patch_map_.CanonicalTriangleForLayoutTriangle(i)] = true;
		std::vector<FluidSurfaceElementState> result;
		for (std::uint64_t cell = 0; cell < domain.Cells().size(); ++cell) {
			if (domain.Cells()[static_cast<std::size_t>(cell)].classification != CellClassification::Cut) continue;
			const auto& provenance = catalog.UsableProvenance(domain, cell);
			if (!std::any_of(provenance.begin(), provenance.end(), [&selected](const auto& point) {
				return point.canonical_triangle < selected.size() && selected[point.canonical_triangle];
			})) continue;
			const auto element = domain.Background().MaterializeElement(cell);
			FluidSurfaceElementState state; state.cell_id = cell; state.nodal_state.resize(element.connectivity.size());
			for (std::size_t node = 0; node < element.connectivity.size(); ++node) {
				const auto local = layout.LocalNode(element.connectivity[node]);
				for (int component = 0; component < 4; ++component)
					state.nodal_state[node][component] = PetscRealPart(values[4*local+component]);
			}
			result.push_back(std::move(state));
		}
		if (result.empty()) throw std::runtime_error("fluid FSI target has no mapped patch cut cells");
		return result;
	}

	SurfaceFieldStamp MakeActualOutputStamp(const std::vector<FluidSurfaceElementState>& state) const
	{
		const auto& context = lifecycle_.Context();
		SurfaceFieldStamp result;
		result.time_s = context.EndTime(); result.step = context.step; result.coupling_iteration = context.coupling_iteration;
		result.reference_mesh_identity_sha256 = fluid_layout_.reference_mesh_identity_sha256;
		result.layout_identity_sha256 = fluid_layout_.layout_identity_sha256;
		result.partition_identity_sha256 = BuildDistributedSurfacePartitionIdentitySha256(fluid_layout_);
		result.producer_state_identity_sha256 = BuildFluidSurfaceTractionStateIdentitySha256(
			moving_.TrialGeometry().Domain(), moving_.TrialGeometry().Surface(), moving_.TrialGeometry().Kinematics(),
			patch_map_, DynamicViscosity(), state);
		return result;
	}

	void DiscardTrialNoexcept() noexcept
	{
		moving_.AbortTrial(); input_.reset(); trial_traction_.reset(); trial_traction_diagnostics_.reset(); trial_target_.reset(); trial_flow_diagnostics_.reset();
		trial_composition_identity_sha256_.clear();
		prepared_traction_.reset(); prepared_traction_diagnostics_.reset(); prepared_target_.reset();
	}

	std::string domain_id_, subsystem_id_;
	FsiCouplingEdge edge_;
	std::vector<DistributedSurfaceInterface> catalog_;
	DistributedSurfaceInterface structure_surface_;
	DistributedSurfaceLayout fluid_layout_, structure_layout_;
	MaterialSurfacePatchMap patch_map_;
	std::unique_ptr<MaterialSurfaceKinematics> committed_full_;
	double moving_options_viscosity_ = 0.0;
	MovingImmersedTransientFlowRuntime moving_;
	FluidSurfaceTractionProjectionOptions traction_options_;
	FsiTrialLifecycle lifecycle_;
	SurfaceFieldStamp expected_kinematics_;
	std::optional<SurfaceKinematics> input_;
	std::optional<SurfaceTraction> trial_traction_, prepared_traction_, committed_traction_;
	std::optional<FluidSurfaceTractionDiagnostics> trial_traction_diagnostics_, prepared_traction_diagnostics_, committed_traction_diagnostics_;
	std::unique_ptr<MaterialSurfaceKinematics> trial_target_;
	std::optional<ImmersedTransientFlowDiagnostics> trial_flow_diagnostics_;
	std::unique_ptr<MaterialSurfaceKinematics> prepared_target_;
	std::string trial_composition_identity_sha256_, committed_composition_identity_sha256_;
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	bool fail_next_late_publication_ = false;
#endif
};

} // namespace iga

#endif

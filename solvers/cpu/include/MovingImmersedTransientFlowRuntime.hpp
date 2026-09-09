#ifndef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_HPP

// Transactional owner for a sequence of moving immersed-flow epochs.  An epoch
// intentionally owns a complete geometry and a complete fixed-geometry flow
// runtime: PETSc objects are never rebound from one cut layout to another.
#include "ImmersedTransientFlowRuntime.hpp"
#include "ImmersedVelocityExtension.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct MovingImmersedTransientFlowOptions {
	CubicCartesianGridSpec grid;
	MovingCutGeometryOptions geometry;
	std::uint32_t extension_layers = 3;
	ImmersedVelocityExtensionOptions extension;
	ImmersedTransientFlowOptions flow;
};

// This is a measurement of a single old -> target endpoint transition.  All
// signs are outward: positive volume rate means expansion and positive flux
// leaves the target fluid domain.
struct MovingImmersedTransientFlowConservationDiagnostics {
	std::string source_geometry_identity_sha256, source_publication_identity_sha256;
	std::string target_geometry_identity_sha256, target_publication_identity_sha256, transition_identity_sha256;
	double source_time_s = 0.0, target_time_s = 0.0, dt_s = 0.0;
	std::uint64_t source_index = 0, target_index = 0;
	double source_audited_volume_m3 = 0.0, target_audited_volume_m3 = 0.0;
	double backward_euler_volume_rate_m3_s = 0.0;
	std::map<int,double> fluid_surface_outward_flow_by_boundary_label_m3_s;
	std::map<int,double> material_surface_outward_flow_by_boundary_label_m3_s;
	std::map<int,double> material_wall_outward_flow_by_boundary_label_m3_s;
	// Exact target-endpoint inner primitives retained for transition consumers.
	double endpoint_volume_divergence_m3_s = 0.0;
	double total_fluid_surface_outward_flow_m3_s = 0.0, total_material_surface_outward_flow_m3_s = 0.0;
	double open_port_outward_flow_m3_s = 0.0, wall_outward_flow_m3_s = 0.0, total_material_wall_outward_flow_m3_s = 0.0;
	double divergence_theorem_defect_m3_s = 0.0, reynolds_defect_m3_s = 0.0;
	double moving_mass_defect_m3_s = 0.0, wall_relative_leakage_m3_s = 0.0;
	double discrete_moving_wall_continuity_defect_m3_s = 0.0;
	double discrete_moving_wall_continuity_normalization_scale_m3_s = 1.0;
	// This is the historical denominator for normalized_open_balance and
	// normalized_wall_leakage.  Retain it explicitly so those values remain
	// reproducible after the richer moving-transition normalization was added.
	double legacy_normalization_scale_m3_s = 1.0;
	// One explicit scale keeps all reported normalized rates comparable.  It is
	// max(reference flow, |G_BE|, |Q_u|, |Q_w|), hence never zero.
	double normalization_scale_m3_s = 1.0;
	double normalized_divergence_theorem_defect = 0.0, normalized_reynolds_defect = 0.0;
	double normalized_moving_mass_defect = 0.0, normalized_wall_relative_leakage = 0.0;
	double normalized_open_balance = 0.0, normalized_wall_leakage = 0.0;
	double normalized_discrete_moving_wall_continuity_defect = 0.0;
};

struct MovingImmersedTransientFlowDiagnostics {
	bool idle = true, trial_active = false, prepared = false;
	std::size_t begin_count = 0, rollback_count = 0, abort_count = 0, prepare_count = 0, finalize_count = 0;
	std::string committed_geometry_identity_sha256, committed_publication_identity_sha256;
	std::string trial_geometry_identity_sha256, trial_publication_identity_sha256;
	std::string source_geometry_identity_sha256, source_publication_identity_sha256;
	std::string target_geometry_identity_sha256, target_publication_identity_sha256;
	std::string extension_hash_sha256, scalar_extension_hash_sha256, map_identity_sha256;
	// This versioned owner identity is complete before the inner trial is
	// published, so PrepareCommit has no outer allocation or hashing work.
	std::string transition_identity_sha256;
	// Immutable geometry/input audit fields are bound before trial publication.
	double source_time_s = 0.0, target_time_s = 0.0, dt_s = 0.0;
	std::uint64_t source_index = 0, target_index = 0;
	double source_audited_volume_m3 = 0.0, target_audited_volume_m3 = 0.0;
};

class MovingImmersedTransientFlowRuntime {
public:
	explicit MovingImmersedTransientFlowRuntime(MaterialSurfaceKinematics initial,
		MovingImmersedTransientFlowOptions options)
		: options_(std::move(options))
	{
		PhaseScope geometry_phase(ProfilePhase::Geometry);
		auto geometry = MovingCutGeometry::Build(options_.grid, std::move(initial), options_.geometry);
		std::unique_ptr<Epoch> epoch(new Epoch(std::move(geometry), options_.flow));
		committed_.swap(epoch); RefreshCommittedDiagnostics();
	}
	~MovingImmersedTransientFlowRuntime() = default;
	MovingImmersedTransientFlowRuntime(const MovingImmersedTransientFlowRuntime&) = delete;
	MovingImmersedTransientFlowRuntime& operator=(const MovingImmersedTransientFlowRuntime&) = delete;

	const MovingCutGeometry& CommittedGeometry() const noexcept { return *committed_->geometry; }
	const ImmersedActiveLayout& CommittedLayout() const noexcept { return committed_->runtime->Layout(); }
	const ImmersedGlobalFlowState& CommittedGlobalState() const { return committed_->runtime->CommittedGlobalState(); }
	std::vector<PetscScalar> CommittedState() const { return committed_->runtime->CommittedState(); }
	const std::vector<int>& ConfiguredWallLabels() const noexcept { return committed_->runtime->ConfiguredWallLabels(); }
	const std::vector<ImmersedFlowPortDefinition>& ConfiguredPorts() const noexcept { return committed_->runtime->ConfiguredPorts(); }
	const ImmersedTransientFlowDiagnostics& CommittedDiagnostics() const noexcept { return committed_->runtime->Diagnostics(); }
	const MovingImmersedTransientFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	bool Idle() const noexcept { return !trial_; }
	bool Prepared() const noexcept { return trial_ && diagnostics_.prepared; }
	const MovingCutGeometry& TrialGeometry() const { RequireTrial("access trial geometry"); return *trial_->geometry; }
	const ImmersedActiveLayout& TrialLayout() const { RequireTrial("access trial layout"); return trial_->runtime->Layout(); }
	std::vector<PetscScalar> TrialState() const { RequireTrial("access trial state"); return trial_->runtime->TrialState(); }
	const ImmersedTransientFlowDiagnostics& TrialDiagnostics() const { RequireTrial("access trial diagnostics"); return trial_->runtime->Diagnostics(); }
	// For an active (including prepared) trial this is a fresh read-only
	// measurement of its target state.  After finalization the exact record
	// prepared for that committed transition is retained.
	MovingImmersedTransientFlowConservationDiagnostics ConservationDiagnostics() const
	{
		if (trial_) return BuildConservation(*trial_, diagnostics_);
		if (!committed_conservation_) throw std::logic_error("moving immersed transient conservation is unavailable before a committed transition");
		return *committed_conservation_;
	}

	void SetCommittedGlobalState(const ImmersedGlobalFlowState& value)
	{
		RequireIdle("set committed state"); committed_->runtime->SetCommittedGlobalState(value); committed_conservation_.reset(); RefreshCommittedDiagnostics();
	}
	void InitializeCommittedGlobalState(const ImmersedGlobalFlowState& value) { SetCommittedGlobalState(value); }
	void SetPortControlValue(const std::string& id, double value)
	{
		RequireIdle("set port control");
		// Validate/mutate the live runtime first.  Its mutation is transactional;
		// only then persist the same option for every future epoch.
		committed_->runtime->SetPortControlValue(id, value);
		for (auto& port : options_.flow.ports) if (port.id == id) { port.value = value; committed_conservation_.reset(); return; }
		throw std::logic_error("moving immersed transient port option is absent");
	}

	void BeginTrial(MaterialSurfaceKinematics target, std::uint64_t target_index, double dt_s)
	{
		PhaseScope geometry_phase(ProfilePhase::Geometry);
		RequireIdle("begin trial");
		const auto& old_state = committed_->runtime->CommittedGlobalState();
		if (old_state.Index() == std::numeric_limits<std::uint64_t>::max())
			throw std::overflow_error("moving immersed transient target index overflows uint64");
		const double target_time_s = target.EvaluatedTimeS();
		if (target_index != old_state.Index()+1 || target_time_s != CheckedTransientTargetTime(old_state.TimeS(), dt_s))
			throw std::invalid_argument("moving immersed transient target time/index is not the exact next epoch");
		target.ValidateNextEpoch(committed_->geometry->Kinematics(), target_time_s, dt_s);

		// Everything through BeginMovingTrial is local.  A failure destroys the
		// target PETSc runtime before its referenced target geometry.
		TestingThrow(FaultStage::Geometry);
		auto target_geometry = MovingCutGeometry::Build(options_.grid, std::move(target), options_.geometry, committed_->geometry.get());
		std::unique_ptr<Epoch> candidate(new Epoch(std::move(target_geometry), options_.flow));
		// This seam is deliberately after target PETSc creation.  It proves an
		// unfinished target epoch cannot retain old-layout handles on failure.
		TestingThrow(FaultStage::InnerRuntime);
		TestingThrow(FaultStage::Extension);
		ImmersedVelocityExtension extension = ImmersedVelocityExtension::Build(*committed_->geometry,
			committed_->runtime->Layout(), old_state, *candidate->geometry, candidate->runtime->Layout(),
			options_.extension_layers, options_.extension);
		TestingThrow(FaultStage::Scalar);
		std::vector<double> old_pressure(old_state.Coefficients().size());
		for (std::size_t i=0; i<old_pressure.size(); ++i) old_pressure[i] = old_state.Coefficients()[i][3];
		const ImmersedScalarExtension scalar = extension.ExtendScalar(old_pressure);
		TestingThrow(FaultStage::Seed);
		const auto seed = BuildSeed(*candidate->runtime, extension, scalar, old_state);
		TestingThrow(FaultStage::Map);
		const auto map = BuildMap(*candidate, committed_->runtime->Layout(), extension, scalar, old_state, seed, dt_s, target_index);
		TestingThrow(FaultStage::BodyForce);
		TestingThrow(FaultStage::MovingBegin);
		candidate->runtime->BeginMovingTrial(target_time_s, target_index, dt_s, extension.TargetHistory(), seed, map);
		// The inner runtime freezes ports, body force, solver configuration, seed,
		// time, and map as part of BeginMovingTrial.  Bind its effective input
		// identity only after that succeeds; all outer hashing remains local until
		// the noexcept pointer/string publication below.
		auto outer = BuildTrialDiagnostics(*candidate, extension, scalar, map);
		// Pointer and string swaps are noexcept; publication is the first visible
		// transaction mutation.
		trial_.swap(candidate); std::swap(diagnostics_, outer); diagnostics_.idle=false; diagnostics_.trial_active=true; ++diagnostics_.begin_count;
	}

	void Assemble() { RequireTrial("assemble"); trial_->runtime->Assemble(); }
	bool SolveTrial() { RequireTrial("solve"); return trial_->runtime->SolveTrial(); }
	void PrepareCommit()
	{
		RequireTrial("prepare commit"); if (diagnostics_.prepared) throw std::logic_error("moving immersed transient trial is already prepared");
		// This potentially allocating/throwing audit is deliberately complete
		// before the inner PrepareCommit transaction is published.
		std::unique_ptr<MovingImmersedTransientFlowConservationDiagnostics> conservation(
			new MovingImmersedTransientFlowConservationDiagnostics(BuildConservation(*trial_, diagnostics_)));
		trial_->runtime->PrepareCommit();
		// No allocation, hashing, or copying after the inner prepared publication.
		trial_conservation_.swap(conservation);
		diagnostics_.prepared=true; ++diagnostics_.prepare_count;
	}
	void FinalizeCommit() noexcept
	{
		if (!trial_ || !diagnostics_.prepared) std::terminate();
		trial_->runtime->FinalizeCommit();
		committed_.swap(trial_); trial_.reset();
		committed_conservation_.swap(trial_conservation_);
		diagnostics_.idle=true; diagnostics_.trial_active=false; diagnostics_.prepared=false; ++diagnostics_.finalize_count;
		RefreshCommittedDiagnosticsNoexcept();
	}
	void Commit() { PrepareCommit(); FinalizeCommit(); }
	void Rollback()
	{
		RequireTrial("rollback"); trial_->runtime->Rollback(); diagnostics_.prepared=false; ++diagnostics_.rollback_count;
	}
	// Unlike the inner fixed-geometry abort this is an owner discard: no target
	// PETSc vector is copied and no old-layout state is touched.
	void AbortTrial() noexcept
	{
		if (!trial_) return;
		trial_.reset(); diagnostics_.idle=true; diagnostics_.trial_active=false; diagnostics_.prepared=false; ++diagnostics_.abort_count;
		trial_conservation_.reset();
		diagnostics_.trial_geometry_identity_sha256.clear(); diagnostics_.trial_publication_identity_sha256.clear();
		diagnostics_.source_geometry_identity_sha256.clear(); diagnostics_.source_publication_identity_sha256.clear(); diagnostics_.target_geometry_identity_sha256.clear(); diagnostics_.target_publication_identity_sha256.clear();
		diagnostics_.extension_hash_sha256.clear(); diagnostics_.scalar_extension_hash_sha256.clear(); diagnostics_.map_identity_sha256.clear(); diagnostics_.transition_identity_sha256.clear();
		diagnostics_.source_time_s=diagnostics_.target_time_s=diagnostics_.dt_s=0.0; diagnostics_.source_index=diagnostics_.target_index=0;
		diagnostics_.source_audited_volume_m3=diagnostics_.target_audited_volume_m3=0.0;
	}
	void AbortPrepared() noexcept { AbortTrial(); }

#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	// Test-only forwarding seam for the inner transactional PrepareCommit
	// failure path; this is not a production convergence or publication bypass.
	void FailNextPrepareForTesting()
	{
		RequireTrial("inject prepare failure"); trial_->runtime->FailNextPrepareForTesting();
	}
	double MaxTrialSurfaceRelativeVelocityNormForTesting() const
	{
		RequireTrial("measure target surface relative velocity"); return trial_->runtime->MaxSurfaceRelativeVelocityNormForTesting();
	}
	class FiniteConservationDefectScopeForTesting {
	public:
		FiniteConservationDefectScopeForTesting() noexcept { InjectFiniteConservationDefectForTesting()=true; }
		~FiniteConservationDefectScopeForTesting() { InjectFiniteConservationDefectForTesting()=false; }
		FiniteConservationDefectScopeForTesting(const FiniteConservationDefectScopeForTesting&) = delete;
		FiniteConservationDefectScopeForTesting& operator=(const FiniteConservationDefectScopeForTesting&) = delete;
	};
	class FiniteDiscreteMovingWallContinuityDefectScopeForTesting {
	public:
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting() noexcept { InjectFiniteDiscreteMovingWallContinuityDefectForTesting()=true; }
		~FiniteDiscreteMovingWallContinuityDefectScopeForTesting() { InjectFiniteDiscreteMovingWallContinuityDefectForTesting()=false; }
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting(const FiniteDiscreteMovingWallContinuityDefectScopeForTesting&) = delete;
		FiniteDiscreteMovingWallContinuityDefectScopeForTesting& operator=(const FiniteDiscreteMovingWallContinuityDefectScopeForTesting&) = delete;
	};
	enum class TestingFault : std::uint8_t { Geometry, InnerRuntime, Extension, Scalar, Seed, Map, BodyForce, MovingBegin, Count };
	class TestingFaultScope {
	public:
		explicit TestingFaultScope(TestingFault fault) noexcept : fault_(fault) { TestingFaults()[static_cast<std::size_t>(fault_)]=true; }
		~TestingFaultScope() { TestingFaults().fill(false); }
		TestingFaultScope(const TestingFaultScope&) = delete;
		TestingFaultScope& operator=(const TestingFaultScope&) = delete;
	private: TestingFault fault_;
	};
#endif

private:
	enum class FaultStage : std::uint8_t { Geometry, InnerRuntime, Extension, Scalar, Seed, Map, BodyForce, MovingBegin, Count };
	// Declaration order is intentional: runtime is destroyed before geometry.
	struct Epoch {
		std::unique_ptr<MovingCutGeometry> geometry;
		std::unique_ptr<ImmersedTransientFlowRuntime> runtime;
		Epoch(std::unique_ptr<MovingCutGeometry> value, const ImmersedTransientFlowOptions& options)
			: geometry(std::move(value)), runtime(new ImmersedTransientFlowRuntime(*geometry, options)) {}
	};
	static std::vector<double> CarryControllers(const ImmersedGlobalFlowState& old_state, const ImmersedActiveLayout& target)
	{
		std::map<std::uint64_t,double> values;
		for (std::size_t i=0; i<old_state.PortIds().size(); ++i) values.emplace(old_state.PortIds()[i], old_state.PortMultipliers()[i]);
		std::vector<double> result; result.reserve(target.PortIds().size());
		for (const auto id : target.PortIds()) { const auto found=values.find(id); if(found==values.end()) throw std::invalid_argument("moving immersed transient controller ID changed across epochs"); result.push_back(found->second); }
		return result;
	}
	static ImmersedGlobalFlowState BuildSeed(const ImmersedTransientFlowRuntime& target, const ImmersedVelocityExtension& extension,
		const ImmersedScalarExtension& scalar, const ImmersedGlobalFlowState& old_state)
	{
		const auto& layout=target.Layout();
		if (extension.TargetHistory().NodeIds()!=layout.NodeIds() || scalar.target_node_ids!=layout.NodeIds() || scalar.target_values.size()!=layout.NodeIds().size())
			throw std::invalid_argument("moving immersed transient continuation does not exactly cover target layout");
		std::vector<std::array<double,4>> fields(layout.NodeIds().size());
		for (std::size_t i=0; i<fields.size(); ++i) { for (int q=0; q<3; ++q) fields[i][q]=extension.TargetHistory().Velocities()[i][q]; fields[i][3]=scalar.target_values[i]; }
		return ImmersedGlobalFlowState(old_state.TimeS(), old_state.Index(), layout, std::move(fields), CarryControllers(old_state,layout), layout.HasGaugeRow(), 0.0);
	}
	static ImmersedMovingTrialMapIdentity BuildMap(const Epoch& target, const ImmersedActiveLayout& old_layout, const ImmersedVelocityExtension& extension,
		const ImmersedScalarExtension& scalar, const ImmersedGlobalFlowState& old_state, const ImmersedGlobalFlowState& seed,
		double dt_s, std::uint64_t target_index)
	{
		return ImmersedMovingTrialMapIdentity::Create(old_state.HashSha256(), old_state.GeometryIdentity(), target.geometry->GeometryIdentitySha256(),
			target.geometry->PublicationIdentitySha256(), old_layout.HashSha256(),
			target.runtime->Layout().HashSha256(), extension.HashSha256(), extension.OperatorHashSha256(), extension.TargetHistory().HashSha256(), scalar.hash_sha256,
			target.runtime->Layout().PortIds(), seed.PortMultipliers(), true, old_state.TimeS(), old_state.Index(), target.geometry->Evaluation().EvaluatedTimeS(), target_index, dt_s);
	}
	MovingImmersedTransientFlowDiagnostics BuildTrialDiagnostics(const Epoch& candidate, const ImmersedVelocityExtension& extension,
		const ImmersedScalarExtension& scalar, const ImmersedMovingTrialMapIdentity& map) const
	{
		MovingImmersedTransientFlowDiagnostics result=diagnostics_;
		result.source_geometry_identity_sha256=committed_->geometry->GeometryIdentitySha256(); result.source_publication_identity_sha256=committed_->geometry->PublicationIdentitySha256();
		result.trial_geometry_identity_sha256=candidate.geometry->GeometryIdentitySha256(); result.trial_publication_identity_sha256=candidate.geometry->PublicationIdentitySha256();
		result.target_geometry_identity_sha256=result.trial_geometry_identity_sha256; result.target_publication_identity_sha256=result.trial_publication_identity_sha256;
		result.extension_hash_sha256=extension.HashSha256(); result.scalar_extension_hash_sha256=scalar.hash_sha256; result.map_identity_sha256=map.HashSha256();
		result.source_time_s=map.SourceTimeS(); result.source_index=map.SourceIndex(); result.target_time_s=map.TargetTimeS(); result.target_index=map.TargetIndex(); result.dt_s=map.DtS();
		result.source_audited_volume_m3=AuditedVolume(*committed_->geometry); result.target_audited_volume_m3=AuditedVolume(*candidate.geometry);
		Sha256 h; immersed_transient_detail::AppendString(h,"MovingImmersedTransientTransition/v3");
		const auto& inner_input=candidate.runtime->Diagnostics().input_hash_sha256;
		if(inner_input.empty()) throw std::logic_error("moving immersed transient inner input identity is absent after begin");
		const std::array<const std::string*,13> inputs{{&committed_->geometry->GeometryIdentitySha256(), &committed_->geometry->PublicationIdentitySha256(), &committed_->runtime->Layout().HashSha256(),
			&committed_->runtime->CommittedGlobalState().HashSha256(), &result.trial_geometry_identity_sha256, &result.trial_publication_identity_sha256,
			&candidate.runtime->Layout().HashSha256(), &result.extension_hash_sha256, &extension.OperatorHashSha256(), &extension.TargetHistory().HashSha256(), &result.scalar_extension_hash_sha256, &result.map_identity_sha256, &inner_input}};
		for (const auto* v : inputs) immersed_transient_detail::AppendString(h,*v);
		h.AppendNormalizedDouble(result.source_time_s); h.AppendLittleEndian64(result.source_index); h.AppendNormalizedDouble(result.target_time_s); h.AppendLittleEndian64(result.target_index); h.AppendNormalizedDouble(result.dt_s);
		h.AppendNormalizedDouble(result.source_audited_volume_m3); h.AppendNormalizedDouble(result.target_audited_volume_m3);
		result.transition_identity_sha256=h.Hex(); return result;
	}
	static double AuditedVolume(const MovingCutGeometry& geometry)
	{
		const double value=geometry.Diagnostics().closed_surface_physical_volume_m3;
		if(!std::isfinite(value) || value<0.0) throw std::runtime_error("moving immersed transient audited geometry volume is invalid");
		return value;
	}
	MovingImmersedTransientFlowConservationDiagnostics BuildConservation(const Epoch& target,
		const MovingImmersedTransientFlowDiagnostics& identity) const
	{
		if(!std::isfinite(identity.dt_s) || !(identity.dt_s>0.0) || !std::isfinite(identity.source_time_s)
			|| !std::isfinite(identity.target_time_s) || identity.source_time_s!=committed_->runtime->CommittedGlobalState().TimeS()
			|| identity.source_index!=committed_->runtime->CommittedGlobalState().Index()
			|| identity.source_index==std::numeric_limits<std::uint64_t>::max() || identity.target_index!=identity.source_index+1
			|| identity.target_time_s!=target.geometry->Evaluation().EvaluatedTimeS()
			|| identity.source_geometry_identity_sha256!=committed_->geometry->GeometryIdentitySha256()
			|| identity.target_geometry_identity_sha256!=target.geometry->GeometryIdentitySha256()
			|| identity.source_publication_identity_sha256!=committed_->geometry->PublicationIdentitySha256()
			|| identity.target_publication_identity_sha256!=target.geometry->PublicationIdentitySha256()
			|| identity.source_audited_volume_m3!=AuditedVolume(*committed_->geometry)
			|| identity.target_audited_volume_m3!=AuditedVolume(*target.geometry))
			throw std::logic_error("moving immersed transient conservation transition identity is inconsistent");
		const auto& inner_diagnostics=target.runtime->Diagnostics();
		if(inner_diagnostics.target_time_s!=identity.target_time_s || inner_diagnostics.dt_s!=identity.dt_s)
			throw std::logic_error("moving immersed transient conservation inner target identity is inconsistent");
		const auto inner=target.runtime->ConservationDiagnostics();
		MovingImmersedTransientFlowConservationDiagnostics result;
		result.source_geometry_identity_sha256=identity.source_geometry_identity_sha256; result.source_publication_identity_sha256=identity.source_publication_identity_sha256;
		result.target_geometry_identity_sha256=identity.target_geometry_identity_sha256; result.target_publication_identity_sha256=identity.target_publication_identity_sha256; result.transition_identity_sha256=identity.transition_identity_sha256;
		result.source_time_s=identity.source_time_s; result.source_index=identity.source_index; result.target_time_s=identity.target_time_s; result.target_index=identity.target_index; result.dt_s=identity.dt_s;
		result.source_audited_volume_m3=identity.source_audited_volume_m3; result.target_audited_volume_m3=identity.target_audited_volume_m3;
		const double volume_delta=result.target_audited_volume_m3-result.source_audited_volume_m3;
		if(!std::isfinite(volume_delta)) throw std::overflow_error("moving immersed transient audited volume delta overflows");
		result.backward_euler_volume_rate_m3_s=volume_delta/result.dt_s;
		result.fluid_surface_outward_flow_by_boundary_label_m3_s=inner.surface_flow_by_boundary_label_m3_s;
		result.material_surface_outward_flow_by_boundary_label_m3_s=inner.material_surface_outward_flow_by_boundary_label_m3_s;
		result.material_wall_outward_flow_by_boundary_label_m3_s=inner.material_wall_outward_flow_by_boundary_label_m3_s;
		result.endpoint_volume_divergence_m3_s=inner.endpoint_volume_divergence_m3_s;
		result.total_fluid_surface_outward_flow_m3_s=inner.total_surface_outward_flow_m3_s;
		result.total_material_surface_outward_flow_m3_s=inner.total_material_surface_outward_flow_m3_s;
		result.open_port_outward_flow_m3_s=inner.open_port_outward_flow_m3_s;
		result.wall_outward_flow_m3_s=inner.wall_outward_flow_m3_s;
		result.total_material_wall_outward_flow_m3_s=inner.total_material_wall_outward_flow_m3_s;
		result.divergence_theorem_defect_m3_s=inner.divergence_theorem_defect_m3_s;
		result.wall_relative_leakage_m3_s=inner.wall_relative_leakage_m3_s;
		result.discrete_moving_wall_continuity_defect_m3_s=inner.discrete_moving_wall_continuity_defect_m3_s;
		result.discrete_moving_wall_continuity_normalization_scale_m3_s=inner.discrete_moving_wall_continuity_normalization_scale_m3_s;
		result.legacy_normalization_scale_m3_s=std::max(options_.flow.flow_controller_reference_flow_m3_s,
			std::abs(result.open_port_outward_flow_m3_s));
		result.normalized_open_balance=inner.normalized_open_balance;
		result.normalized_wall_leakage=inner.normalized_wall_leakage;
		result.normalized_discrete_moving_wall_continuity_defect=inner.normalized_discrete_moving_wall_continuity_defect;
		result.reynolds_defect_m3_s=result.backward_euler_volume_rate_m3_s-result.total_material_surface_outward_flow_m3_s;
		result.moving_mass_defect_m3_s=result.backward_euler_volume_rate_m3_s+result.total_fluid_surface_outward_flow_m3_s-result.total_material_surface_outward_flow_m3_s;
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if(InjectFiniteConservationDefectForTesting()) result.reynolds_defect_m3_s=1.0;
		if(InjectFiniteDiscreteMovingWallContinuityDefectForTesting()) result.discrete_moving_wall_continuity_defect_m3_s=1.0;
#endif
		result.normalization_scale_m3_s=std::max({options_.flow.flow_controller_reference_flow_m3_s,std::abs(result.backward_euler_volume_rate_m3_s),std::abs(result.total_fluid_surface_outward_flow_m3_s),std::abs(result.total_material_surface_outward_flow_m3_s)});
		for(const double value : {result.endpoint_volume_divergence_m3_s,result.total_fluid_surface_outward_flow_m3_s,result.open_port_outward_flow_m3_s,result.wall_outward_flow_m3_s,result.total_material_surface_outward_flow_m3_s,result.total_material_wall_outward_flow_m3_s,result.divergence_theorem_defect_m3_s,result.wall_relative_leakage_m3_s,result.discrete_moving_wall_continuity_defect_m3_s,result.discrete_moving_wall_continuity_normalization_scale_m3_s,result.legacy_normalization_scale_m3_s,result.reynolds_defect_m3_s,result.moving_mass_defect_m3_s,result.normalization_scale_m3_s}) if(!std::isfinite(value)) throw std::runtime_error("moving immersed transient conservation value is nonfinite");
		auto validate_map=[](const std::map<int,double>& values,double total,const char* what) {
			double sum=0.0; for(const auto& value:values) { if(!std::isfinite(value.second)||!std::isfinite(sum+value.second)) throw std::runtime_error(std::string("moving immersed transient ")+what+" label flux is nonfinite"); sum+=value.second; }
			if(std::abs(sum-total)>128.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(total))*std::max<std::size_t>(1,values.size())) throw std::logic_error(std::string("moving immersed transient ")+what+" label fluxes do not reconcile with total");
		};
		validate_map(result.fluid_surface_outward_flow_by_boundary_label_m3_s,result.total_fluid_surface_outward_flow_m3_s,"fluid surface");
		validate_map(result.material_surface_outward_flow_by_boundary_label_m3_s,result.total_material_surface_outward_flow_m3_s,"material surface");
		validate_map(result.material_wall_outward_flow_by_boundary_label_m3_s,result.total_material_wall_outward_flow_m3_s,"material wall");
		if(!immersed_transient_detail::MovingWallContinuityIdentityReconciles(result.open_port_outward_flow_m3_s,
			result.wall_outward_flow_m3_s,result.total_material_wall_outward_flow_m3_s,
			result.discrete_moving_wall_continuity_defect_m3_s,result.wall_relative_leakage_m3_s,result.total_fluid_surface_outward_flow_m3_s))
			throw std::logic_error("moving immersed transient moving-wall continuity roundoff identity does not reconcile");
		if(!(result.normalization_scale_m3_s>0.0)) throw std::runtime_error("moving immersed transient conservation normalization scale is invalid");
		if(!(result.legacy_normalization_scale_m3_s>0.0)) throw std::runtime_error("moving immersed transient legacy normalization scale is invalid");
		if(!(result.discrete_moving_wall_continuity_normalization_scale_m3_s>0.0)) throw std::runtime_error("moving immersed transient discrete continuity normalization scale is invalid");
		result.normalized_divergence_theorem_defect=std::abs(result.divergence_theorem_defect_m3_s)/result.normalization_scale_m3_s;
		result.normalized_reynolds_defect=std::abs(result.reynolds_defect_m3_s)/result.normalization_scale_m3_s;
		result.normalized_moving_mass_defect=std::abs(result.moving_mass_defect_m3_s)/result.normalization_scale_m3_s;
		result.normalized_wall_relative_leakage=std::abs(result.wall_relative_leakage_m3_s)/result.normalization_scale_m3_s;
		if(!std::isfinite(result.normalized_divergence_theorem_defect)||!std::isfinite(result.normalized_reynolds_defect)||!std::isfinite(result.normalized_moving_mass_defect)||!std::isfinite(result.normalized_wall_relative_leakage)||!std::isfinite(result.normalized_open_balance)||!std::isfinite(result.normalized_wall_leakage)||!std::isfinite(result.normalized_discrete_moving_wall_continuity_defect)) throw std::runtime_error("moving immersed transient normalized conservation quotient is nonfinite");
		return result;
	}
	void RequireIdle(const char* action) const { if (trial_) throw std::logic_error(std::string("cannot ")+action+" while a moving immersed transient trial exists"); }
	void RequireTrial(const char* action) const { if (!trial_) throw std::logic_error(std::string("cannot ")+action+" without a moving immersed transient trial"); }
	void RefreshCommittedDiagnostics()
	{
		diagnostics_.committed_geometry_identity_sha256=committed_->geometry->GeometryIdentitySha256();
		diagnostics_.committed_publication_identity_sha256=committed_->geometry->PublicationIdentitySha256();
	}
	void RefreshCommittedDiagnosticsNoexcept() noexcept
	{
		// The target identities are already held by the previous trial diagnostics,
		// so move them into committed fields without allocation.
		diagnostics_.committed_geometry_identity_sha256.swap(diagnostics_.trial_geometry_identity_sha256);
		diagnostics_.committed_publication_identity_sha256.swap(diagnostics_.trial_publication_identity_sha256);
		diagnostics_.trial_geometry_identity_sha256.clear(); diagnostics_.trial_publication_identity_sha256.clear();
	}
	static void TestingThrow(FaultStage fault)
	{
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
		if (TestingFaults()[static_cast<std::size_t>(fault)]) throw std::runtime_error("injected moving immersed transient begin failure");
#else
		(void)fault;
#endif
	}
#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	static std::array<bool,static_cast<std::size_t>(FaultStage::Count)>& TestingFaults()
	{ static std::array<bool,static_cast<std::size_t>(FaultStage::Count)> faults{}; return faults; }
	static bool& InjectFiniteConservationDefectForTesting() noexcept { static bool value=false; return value; }
	static bool& InjectFiniteDiscreteMovingWallContinuityDefectForTesting() noexcept { static bool value=false; return value; }
#endif
	MovingImmersedTransientFlowOptions options_;
	std::unique_ptr<Epoch> committed_, trial_;
	std::unique_ptr<MovingImmersedTransientFlowConservationDiagnostics> committed_conservation_, trial_conservation_;
	MovingImmersedTransientFlowDiagnostics diagnostics_;
};

} // namespace iga

#endif

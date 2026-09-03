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

struct MovingImmersedTransientFlowDiagnostics {
	bool idle = true, trial_active = false, prepared = false;
	std::size_t begin_count = 0, rollback_count = 0, abort_count = 0, prepare_count = 0, finalize_count = 0;
	std::string committed_geometry_identity_sha256, committed_publication_identity_sha256;
	std::string trial_geometry_identity_sha256, trial_publication_identity_sha256;
	std::string extension_hash_sha256, scalar_extension_hash_sha256, map_identity_sha256;
	// This versioned owner identity is complete before the inner trial is
	// published, so PrepareCommit has no outer allocation or hashing work.
	std::string transition_identity_sha256;
};

class MovingImmersedTransientFlowRuntime {
public:
	explicit MovingImmersedTransientFlowRuntime(PrescribedSurfaceMotion::Evaluation initial,
		MovingImmersedTransientFlowOptions options)
		: options_(std::move(options))
	{
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
	const ImmersedTransientFlowDiagnostics& CommittedDiagnostics() const noexcept { return committed_->runtime->Diagnostics(); }
	const MovingImmersedTransientFlowDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	bool Idle() const noexcept { return !trial_; }
	bool Prepared() const noexcept { return trial_ && diagnostics_.prepared; }
	const MovingCutGeometry& TrialGeometry() const { RequireTrial("access trial geometry"); return *trial_->geometry; }
	const ImmersedActiveLayout& TrialLayout() const { RequireTrial("access trial layout"); return trial_->runtime->Layout(); }
	std::vector<PetscScalar> TrialState() const { RequireTrial("access trial state"); return trial_->runtime->TrialState(); }
	const ImmersedTransientFlowDiagnostics& TrialDiagnostics() const { RequireTrial("access trial diagnostics"); return trial_->runtime->Diagnostics(); }

	void SetCommittedGlobalState(const ImmersedGlobalFlowState& value)
	{
		RequireIdle("set committed state"); committed_->runtime->SetCommittedGlobalState(value); RefreshCommittedDiagnostics();
	}
	void InitializeCommittedGlobalState(const ImmersedGlobalFlowState& value) { SetCommittedGlobalState(value); }
	void SetPortControlValue(const std::string& id, double value)
	{
		RequireIdle("set port control");
		// Validate/mutate the live runtime first.  Its mutation is transactional;
		// only then persist the same option for every future epoch.
		committed_->runtime->SetPortControlValue(id, value);
		for (auto& port : options_.flow.ports) if (port.id == id) { port.value = value; return; }
		throw std::logic_error("moving immersed transient port option is absent");
	}

	void BeginTrial(PrescribedSurfaceMotion::Evaluation target, std::uint64_t target_index, double dt_s)
	{
		RequireIdle("begin trial");
		const auto& old_state = committed_->runtime->CommittedGlobalState();
		if (old_state.Index() == std::numeric_limits<std::uint64_t>::max())
			throw std::overflow_error("moving immersed transient target index overflows uint64");
		const double target_time_s = target.EvaluatedTimeS();
		if (target_index != old_state.Index()+1 || target_time_s != CheckedTransientTargetTime(old_state.TimeS(), dt_s))
			throw std::invalid_argument("moving immersed transient target time/index is not the exact next epoch");

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
		trial_->runtime->PrepareCommit();
		// No allocation, hashing, or copying after the inner prepared publication.
		diagnostics_.prepared=true; ++diagnostics_.prepare_count;
	}
	void FinalizeCommit() noexcept
	{
		if (!trial_ || !diagnostics_.prepared) std::terminate();
		trial_->runtime->FinalizeCommit();
		committed_.swap(trial_); trial_.reset();
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
		diagnostics_.trial_geometry_identity_sha256.clear(); diagnostics_.trial_publication_identity_sha256.clear();
		diagnostics_.extension_hash_sha256.clear(); diagnostics_.scalar_extension_hash_sha256.clear(); diagnostics_.map_identity_sha256.clear(); diagnostics_.transition_identity_sha256.clear();
	}
	void AbortPrepared() noexcept { AbortTrial(); }

#ifdef IGA_MOVING_IMMERSED_TRANSIENT_FLOW_RUNTIME_TESTING
	// Test-only forwarding seam for the inner transactional PrepareCommit
	// failure path; this is not a production convergence or publication bypass.
	void FailNextPrepareForTesting()
	{
		RequireTrial("inject prepare failure"); trial_->runtime->FailNextPrepareForTesting();
	}
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
		result.trial_geometry_identity_sha256=candidate.geometry->GeometryIdentitySha256(); result.trial_publication_identity_sha256=candidate.geometry->PublicationIdentitySha256();
		result.extension_hash_sha256=extension.HashSha256(); result.scalar_extension_hash_sha256=scalar.hash_sha256; result.map_identity_sha256=map.HashSha256();
		Sha256 h; immersed_transient_detail::AppendString(h,"MovingImmersedTransientTransition/v2");
		const auto& inner_input=candidate.runtime->Diagnostics().input_hash_sha256;
		if(inner_input.empty()) throw std::logic_error("moving immersed transient inner input identity is absent after begin");
		const std::array<const std::string*,13> inputs{{&committed_->geometry->GeometryIdentitySha256(), &committed_->geometry->PublicationIdentitySha256(), &committed_->runtime->Layout().HashSha256(),
			&committed_->runtime->CommittedGlobalState().HashSha256(), &result.trial_geometry_identity_sha256, &result.trial_publication_identity_sha256,
			&candidate.runtime->Layout().HashSha256(), &result.extension_hash_sha256, &extension.OperatorHashSha256(), &extension.TargetHistory().HashSha256(), &result.scalar_extension_hash_sha256, &result.map_identity_sha256, &inner_input}};
		for (const auto* v : inputs) immersed_transient_detail::AppendString(h,*v);
		result.transition_identity_sha256=h.Hex(); return result;
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
#endif
	MovingImmersedTransientFlowOptions options_;
	std::unique_ptr<Epoch> committed_, trial_;
	MovingImmersedTransientFlowDiagnostics diagnostics_;
};

} // namespace iga

#endif

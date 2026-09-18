#include "StrongFluidStructureCoupling.hpp"

#include <cassert>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string H(const std::string& value)
{
	iga::Sha256 hash;
	hash.Append(value.data(), value.size());
	return hash.Hex();
}

iga::DistributedSurfaceLayout L()
{
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = H("mesh");
	layout.global_node_count = 3;
	layout.partition_count = 1;
	layout.owned_global_node_ids = {0, 1, 2};
	layout.reference_positions = {{0, {{0, 0, 0}}}, {1, {{1, 0, 0}}}, {2, {{0, 1, 0}}}};
	layout.reference_triangles = {{{0, 1, 2}}};
	// Uneven weights: sum .36; nodes zero and two deliberately dominate the RMS.
	layout.owned_reference_lumped_areas_m2 = {.03, .09, .24};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return layout;
}

iga::FsiCouplingEdge E()
{
	return {"edge", {"fluid", "f", "wall"}, {"structure", "s", "wall"},
		iga::FsiCouplingLaw::FluidStructureTractionKinematics};
}

iga::SurfaceFieldStamp S(const iga::DistributedSurfaceLayout& layout, double time, int step,
	int iteration, const std::string& producer)
{
	return {time, static_cast<std::uint64_t>(step), static_cast<std::uint64_t>(iteration),
		layout.reference_mesh_identity_sha256, layout.layout_identity_sha256,
		iga::BuildDistributedSurfacePartitionIdentitySha256(layout), H(producer)};
}

bool InEnvelope(const iga::SurfaceFieldStamp& stamp, const iga::SurfaceFieldStampEnvelope& envelope)
{
	return stamp.time_s == envelope.time_s && stamp.step == envelope.step
		&& stamp.coupling_iteration == envelope.coupling_iteration
		&& stamp.reference_mesh_identity_sha256 == envelope.reference_mesh_identity_sha256
		&& stamp.layout_identity_sha256 == envelope.layout_identity_sha256
		&& stamp.partition_identity_sha256 == envelope.partition_identity_sha256
		&& iga::IsLowercaseSha256(stamp.producer_state_identity_sha256);
}

struct F {
	iga::DistributedSurfaceLayout l;
	iga::FsiCouplingEdge e;
	iga::DomainStepContext step{};
	iga::SurfaceKinematics in{};
	iga::SurfaceTraction tr{};
	iga::SurfaceFieldStamp expected_input{};
	iga::SurfaceFieldStampEnvelope expected_output{};
	std::vector<iga::SurfaceKinematics> inputs;
	std::optional<iga::SurfaceTraction> prepared_traction, committed_traction;
	std::string prepared_composition_identity_sha256, committed_composition_identity_sha256;
	int it = 0, commits = 0, prepares = 0, prevalidations = 0;
	bool active = false, awaiting = false, solved = false, prepared = false;
	bool fail_solve = false, fail_prepare = false, fail_pre = false;

	F(iga::DistributedSurfaceLayout layout, iga::FsiCouplingEdge edge) : l(std::move(layout)), e(std::move(edge)) {}
	const iga::FsiCouplingEdge& StrongCouplingEdge() const noexcept { return e; }
	const iga::DistributedSurfaceLayout& StrongCouplingLayout() const noexcept { return l; }
	void BeginMacroStep(const iga::DomainStepContext& value) { if (active) throw std::runtime_error("active"); step = value; active = true; }
	void BeginCouplingIteration(std::uint64_t iteration, const iga::SurfaceFieldStamp& expected,
		const iga::SurfaceFieldStampEnvelope& output)
	{ it = static_cast<int>(iteration); expected_input = expected; expected_output = output; awaiting = true; }
	void SetSurfaceKinematics(const std::string& id, const iga::SurfaceKinematics& value)
	{
		if (id != e.fluid.interface_id || !awaiting
			|| iga::BuildSurfaceFieldStampIdentitySha256(value.stamp, l) != iga::BuildSurfaceFieldStampIdentitySha256(expected_input, l))
			throw std::runtime_error("fluid input stamp");
		in = value; inputs.push_back(value); awaiting = false;
	}
	void SolveFluidTrial()
	{
		if (fail_solve || awaiting) throw std::runtime_error("fluid solve");
		tr.interface = e.fluid; tr.stamp = S(l, step.start_time_s+step.dt_s, step.step_index, it, "fluid"+std::to_string(it));
		if (!InEnvelope(tr.stamp, expected_output)) throw std::runtime_error("fluid output envelope");
		tr.traction_on_structure_pa.assign(3, {{0, 0, 0}}); tr.consistent_nodal_force_n.assign(3, {{0, 0, 0}});
		for (int i = 0; i < 3; ++i) tr.consistent_nodal_force_n[i][2] = in.displacement_m[i][2];
		tr.projection_identity_sha256 = H("projection"); solved = true;
	}
	iga::SurfaceTraction GetSurfaceTraction(const std::string& id) const { if (id != e.fluid.interface_id || !solved) throw std::runtime_error("traction"); return tr; }
	void RejectCouplingIteration() { solved = false; awaiting = false; }
	void PrepareCommitStep()
	{
		++prepares;
		if (!solved || fail_prepare) throw std::runtime_error("fluid prepare");
		prepared_traction = tr;
		prepared_composition_identity_sha256 = H("fluid composition"+std::to_string(it));
		prepared = true;
	}
	std::string CoordinatorPreparedCommittedTractionIdentitySha256() const
	{
		if (!prepared || !prepared_traction) throw std::runtime_error("fluid identity");
		return iga::BuildSurfaceTractionIdentitySha256(*prepared_traction, l);
	}
	std::string CoordinatorPreparedCommittedCompositionIdentitySha256() const
	{
		if (!prepared || prepared_composition_identity_sha256.empty()) throw std::runtime_error("fluid composition");
		return prepared_composition_identity_sha256;
	}
	const iga::SurfaceTraction& CommittedTraction() const
	{ if (!committed_traction) throw std::runtime_error("fluid committed traction"); return *committed_traction; }
	const std::string& CommittedCompositionIdentitySha256() const
	{ if (committed_composition_identity_sha256.empty()) throw std::runtime_error("fluid committed composition"); return committed_composition_identity_sha256; }
	bool HasCommittedPublication() const noexcept
	{ return committed_traction.has_value() || !committed_composition_identity_sha256.empty(); }
	void AbortStep() noexcept
	{
		active = awaiting = solved = prepared = false;
		prepared_traction.reset(); prepared_composition_identity_sha256.clear();
	}
	void CoordinatorRequireFinalizeAllowed()
	{
		++prevalidations;
		if (!prepared || !prepared_traction || prepared_composition_identity_sha256.empty() || fail_pre)
			throw std::runtime_error("fluid prevalidate");
	}
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		++commits;
		using std::swap;
		swap(committed_traction, prepared_traction);
		committed_composition_identity_sha256.swap(prepared_composition_identity_sha256);
		AbortStep();
	}
};

// Keep the accepted fluid geometry distinct from the raw membrane state.
struct MovingHistoryFluid : F {
    using F::F;
    std::array<std::array<double,3>,3> history{};
    std::vector<std::array<double,3>> StrongCouplingCommittedDisplacementM() const
    { return {history.begin(),history.end()}; }
    void SetSurfaceKinematics(const std::string& id,const iga::SurfaceKinematics& value)
    {
        for(std::size_t i=0;i<history.size();++i) for(int c=0;c<3;++c)
            if(std::abs(value.displacement_m[i][c]-history[i][c]-step.dt_s*value.velocity_m_per_s[i][c])>1e-12)
                throw std::runtime_error("moving interface violated its own backward-Euler history");
        F::SetSurfaceKinematics(id,value);
    }
    void CoordinatorFinalizeCommitNoexcept() noexcept
    {
        for(std::size_t i=0;i<history.size();++i) history[i]=in.displacement_m[i];
        F::CoordinatorFinalizeCommitNoexcept();
    }
};

struct M {
	iga::DistributedSurfaceLayout l;
	iga::FsiCouplingEdge e;
	iga::StrongCouplingStructureSnapshot snap;
	iga::DomainStepContext step{};
	iga::SurfaceTraction tr{};
	iga::SurfaceKinematics raw{};
	iga::SurfaceFieldStamp expected_input{};
	iga::SurfaceFieldStampEnvelope expected_output{};
	std::optional<iga::SurfaceKinematics> prepared_kinematics, committed_kinematics;
	std::optional<iga::StrongCouplingStructureSnapshot> prepared_snapshot;
	double b, c;
	int it = 0, commits = 0, prepares = 0, prevalidations = 0;
	bool active = false, awaiting = false, solved = false, prepared = false;
	bool fail_solve = false, fail_prepare = false, fail_pre = false;

	M(iga::DistributedSurfaceLayout layout, iga::FsiCouplingEdge edge, double offset, double slope)
		: l(std::move(layout)), e(std::move(edge)), b(offset), c(slope)
	{
		snap.displacement_m = {0, 0, 0}; snap.velocity_m_per_s = {0, 0, 0};
		snap.immutable_reference_normals.assign(3, {{0, 0, 1}}); snap.clamped = {false, true, false};
		snap.committed_state_identity_sha256 = H("state"); snap.model_identity_sha256 = H("model");
	}
	const iga::FsiCouplingEdge& StrongCouplingEdge() const noexcept { return e; }
	const iga::DistributedSurfaceLayout& StrongCouplingLayout() const noexcept { return l; }
	iga::StrongCouplingStructureSnapshot StrongCouplingCommittedSnapshot() const { return snap; }
	void BeginMacroStep(const iga::DomainStepContext& value) { if (active) throw std::runtime_error("active"); step = value; active = true; }
	void BeginCouplingIteration(std::uint64_t iteration, const iga::SurfaceFieldStamp& expected,
		const iga::SurfaceFieldStampEnvelope& output)
	{ it = static_cast<int>(iteration); expected_input = expected; expected_output = output; awaiting = true; }
	void SetSurfaceTraction(const std::string& id, const iga::SurfaceTraction& value)
	{
		if (id != e.structure.interface_id || !awaiting
			|| iga::BuildSurfaceFieldStampIdentitySha256(value.stamp, l) != iga::BuildSurfaceFieldStampIdentitySha256(expected_input, l))
			throw std::runtime_error("structure input stamp");
		tr = value; awaiting = false;
	}
	void SolveMembraneTrial()
	{
		if (fail_solve || awaiting) throw std::runtime_error("structure solve");
		raw.interface = e.structure; raw.stamp = S(l, step.start_time_s+step.dt_s, step.step_index, it, "structure"+std::to_string(it));
		if (!InEnvelope(raw.stamp, expected_output)) throw std::runtime_error("structure output envelope");
		raw.displacement_m.assign(3, {{0, 0, 0}}); raw.velocity_m_per_s.assign(3, {{0, 0, 0}});
		for (int i = 0; i < 3; ++i) {
			const double displacement = snap.clamped[i] ? 0.0 : b+c*tr.consistent_nodal_force_n[i][2];
			raw.displacement_m[i][2] = displacement;
			raw.velocity_m_per_s[i][2] = (displacement-snap.displacement_m[i])/step.dt_s;
		}
		solved = true;
	}
	iga::SurfaceKinematics GetSurfaceKinematics(const std::string& id) const { if (id != e.structure.interface_id || !solved) throw std::runtime_error("kinematics"); return raw; }
	void RejectCouplingIteration() { solved = false; awaiting = false; }
	std::string StateIdentity(const iga::StrongCouplingStructureSnapshot& value) const
	{
		iga::Sha256 hash;
		const std::string kind = "strong coupling mock structure state/v1";
		hash.Append(kind.data(), kind.size());
		for (std::size_t i = 0; i < value.displacement_m.size(); ++i) {
			hash.AppendNormalizedDouble(value.displacement_m[i]);
			hash.AppendNormalizedDouble(value.velocity_m_per_s[i]);
		}
		return hash.Hex();
	}
	void PrepareCommitStep()
	{
		++prepares;
		if (!solved || fail_prepare) throw std::runtime_error("structure prepare");
		prepared_kinematics = raw;
		prepared_snapshot = snap;
		for (int i = 0; i < 3; ++i) {
			prepared_snapshot->displacement_m[i] = raw.displacement_m[i][2];
			prepared_snapshot->velocity_m_per_s[i] = raw.velocity_m_per_s[i][2];
		}
		prepared_snapshot->committed_state_identity_sha256 = StateIdentity(*prepared_snapshot);
		prepared = true;
	}
	std::string CoordinatorPreparedCommittedKinematicsIdentitySha256() const
	{
		if (!prepared || !prepared_kinematics) throw std::runtime_error("structure identity");
		return iga::BuildSurfaceKinematicsIdentitySha256(*prepared_kinematics, l);
	}
	std::string CoordinatorPreparedCommittedStateIdentitySha256() const
	{
		if (!prepared || !prepared_snapshot) throw std::runtime_error("structure state");
		return prepared_snapshot->committed_state_identity_sha256;
	}
	const iga::SurfaceKinematics& CommittedKinematics() const
	{ if (!committed_kinematics) throw std::runtime_error("structure committed kinematics"); return *committed_kinematics; }
	const std::string& CommittedStateIdentitySha256() const
	{ return snap.committed_state_identity_sha256; }
	bool HasCommittedPublication() const noexcept { return committed_kinematics.has_value(); }
	void AbortStep() noexcept
	{
		active = awaiting = solved = prepared = false;
		prepared_kinematics.reset(); prepared_snapshot.reset();
	}
	void CoordinatorRequireFinalizeAllowed()
	{
		++prevalidations;
		if (!prepared || !prepared_kinematics || !prepared_snapshot || fail_pre)
			throw std::runtime_error("structure prevalidate");
	}
	void CoordinatorFinalizeCommitNoexcept() noexcept
	{
		++commits;
		using std::swap;
		swap(committed_kinematics, prepared_kinematics);
		swap(snap, *prepared_snapshot);
		AbortStep();
	}
};

iga::StrongFluidStructureCouplingOptions O(int iterations = 32)
{
	iga::StrongFluidStructureCouplingOptions options;
	options.maximum_iterations = iterations; options.absolute_displacement_tolerance_m = 1.e-9;
	options.relative_displacement_tolerance = 1.e-9; options.reference_displacement_scale_m = 1.0;
	options.aitken_controls.initial_relaxation = .5; options.aitken_controls.minimum_relaxation = .05;
	options.aitken_controls.maximum_relaxation = .75;
	return options;
}

template <class X> void R(X&& action)
{
	bool failed = false; try { action(); } catch (const std::runtime_error&) { failed = true; } assert(failed);
}

} // namespace

int main()
{
	const auto layout = L(); const auto edge = E(); const iga::DomainStepContext step{1, 0, .1};
    {
        MovingHistoryFluid fluid(layout,edge); M membrane(layout,edge,1,-3);
        auto options=O(); options.absolute_displacement_tolerance_m=.01;
        options.aitken_controls.initial_relaxation=.15; options.aitken_controls.maximum_relaxation=.2;
        iga::StrongFluidStructureCoupling<MovingHistoryFluid,M> coordinator(fluid,membrane,edge,layout,options);
        const auto first=coordinator.Execute(step);
        assert(first.converged && first.history.back().area_weighted_rms_residual_m>1e-8);
        assert(std::abs(membrane.snap.displacement_m[0]-fluid.history[0][2])>1e-8);
        const auto input_count=fluid.inputs.size();
        const auto second=coordinator.Execute({2,.1,.1});
        assert(second.converged && second.iterations>1 && fluid.inputs.size()>input_count+1);
        assert(fluid.commits==2 && membrane.commits==2);
    }
	// r=(1,0,1), A=(.03,.09,.24): sqrt((.03+.24)/.36)=sqrt(3/4).
	{
		F fluid(layout, edge); M membrane(layout, edge, 1, 0); auto options = O();
		options.absolute_displacement_tolerance_m = 2.5;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, options);
		auto result = coordinator.Execute(step); const auto& d = result.history.front();
		assert(result.iterations == 1 && d.converged && !d.aitken_proposal_applied && !d.relaxation && !d.aitken_status);
		assert(std::abs(d.area_weighted_rms_residual_m-std::sqrt(.75)) < 1.e-14);
		assert(d.max_residual_m == 1.0 && d.displacement_scale_m == 1.0
			&& std::abs(d.convergence_threshold_m-2.500000001) < 1.e-15);
		assert(result.accepted_relaxed_kinematics_identity_sha256 == d.current_kinematics_identity_sha256);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 1, -3); iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		auto result = coordinator.Execute(step); assert(result.converged && result.iterations == 3 && fluid.commits == 1 && membrane.commits == 1);
		const auto& first = result.history[0]; const auto& second = result.history[1];
		assert(first.aitken_proposal_applied && *first.relaxation == .5 && *first.aitken_status == iga::AitkenRelaxationStatus::Initial);
		assert(second.aitken_proposal_applied && *second.relaxation == .25 && *second.unclamped_relaxation == .25
			&& *second.aitken_status == iga::AitkenRelaxationStatus::Dynamic && second.aitken_had_previous_residual);
		assert(!first.aitken_proposal_identity_sha256.empty() && !first.aitken_control_state_identity_sha256.empty());
		assert(fluid.inputs[1].displacement_m[0][2] == .5 && fluid.inputs[1].velocity_m_per_s[0][2] == 5.0);
		assert(fluid.inputs[1].displacement_m[1][2] == 0.0 && fluid.inputs[1].velocity_m_per_s[1][2] == 0.0);
		assert(first.current_kinematics_identity_sha256 != first.relaxed_kinematics_identity_sha256);
		assert(fluid.inputs[0].stamp.producer_state_identity_sha256 != fluid.inputs[1].stamp.producer_state_identity_sha256);
		assert(result.final_committed_fluid_traction_identity_sha256 == result.accepted_fluid_traction_identity_sha256);
		assert(result.final_committed_structure_kinematics_identity_sha256 == result.accepted_raw_kinematics_identity_sha256);
		assert(!result.final_committed_structure_state_identity_sha256.empty() && !result.final_committed_fluid_composition_identity_sha256.empty());
		assert(result.final_committed_fluid_traction_identity_sha256
			== iga::BuildSurfaceTractionIdentitySha256(fluid.CommittedTraction(), layout));
		assert(result.final_committed_fluid_composition_identity_sha256 == fluid.CommittedCompositionIdentitySha256());
		assert(result.final_committed_structure_kinematics_identity_sha256
			== iga::BuildSurfaceKinematicsIdentitySha256(membrane.CommittedKinematics(), layout));
		assert(result.final_committed_structure_state_identity_sha256 == membrane.CommittedStateIdentitySha256());
		assert(result.final_committed_structure_state_identity_sha256
			== membrane.StateIdentity(membrane.StrongCouplingCommittedSnapshot()));
		const auto replay = coordinator.Execute({2, .1, .1}); assert(replay.identity_sha256 != result.identity_sha256);
		F replay_fluid(layout, edge); M replay_membrane(layout, edge, 1, -3);
		iga::StrongFluidStructureCoupling<F, M> replay_coordinator(replay_fluid, replay_membrane, edge, layout, O());
		const auto exact_replay = replay_coordinator.Execute(step);
		assert(exact_replay.identity_sha256 == result.identity_sha256
			&& exact_replay.history[1].aitken_proposal_identity_sha256 == result.history[1].aitken_proposal_identity_sha256);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); auto options = O(); options.fail_before_finalize_for_testing = true;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, options);
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active && !fluid.HasCommittedPublication() && !membrane.HasCommittedPublication());
		assert(coordinator.Execute(step).converged && fluid.commits == 1 && membrane.commits == 1);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); fluid.fail_pre = true;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active && !fluid.HasCommittedPublication() && !membrane.HasCommittedPublication());
		assert(fluid.prevalidations == 1 && membrane.prevalidations == 0);
		fluid.fail_pre = false; assert(coordinator.Execute(step).converged);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); membrane.fail_pre = true;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active && !fluid.HasCommittedPublication() && !membrane.HasCommittedPublication());
		assert(fluid.prevalidations == 1 && membrane.prevalidations == 1);
		membrane.fail_pre = false; assert(coordinator.Execute(step).converged && fluid.commits == 1 && membrane.commits == 1);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); fluid.fail_prepare = true;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active && !fluid.HasCommittedPublication() && !membrane.HasCommittedPublication());
		assert(fluid.prepares == 1 && membrane.prepares == 0 && fluid.prevalidations == 0 && membrane.prevalidations == 0);
		fluid.fail_prepare = false; assert(coordinator.Execute(step).converged && fluid.commits == 1 && membrane.commits == 1);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); membrane.fail_prepare = true;
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active && !fluid.HasCommittedPublication() && !membrane.HasCommittedPublication());
		assert(fluid.prepares == 1 && membrane.prepares == 1 && fluid.prevalidations == 0 && membrane.prevalidations == 0);
		membrane.fail_prepare = false; assert(coordinator.Execute(step).converged && fluid.commits == 1 && membrane.commits == 1);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 1, 1); iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O(2));
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0 && !fluid.active && !membrane.active);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); membrane.snap.immutable_reference_normals[0] = {{0, 0, 2}};
		iga::StrongFluidStructureCoupling<F, M> coordinator(fluid, membrane, edge, layout, O());
		R([&] { coordinator.Execute(step); }); assert(fluid.commits == 0 && membrane.commits == 0);
	}
	{
		F fluid(layout, edge); M membrane(layout, edge, 0, 0); auto loose = O(); auto tight = O(); tight.absolute_displacement_tolerance_m = 2.e-9;
		iga::StrongFluidStructureCoupling<F, M> first(fluid, membrane, edge, layout, loose); const auto a = first.Execute(step);
		F fluid_two(layout, edge); M membrane_two(layout, edge, 0, 0);
		iga::StrongFluidStructureCoupling<F, M> second(fluid_two, membrane_two, edge, layout, tight); const auto b = second.Execute(step);
		assert(a.options_identity_sha256 != b.options_identity_sha256 && a.identity_sha256 != b.identity_sha256);
	}
	return 0;
}

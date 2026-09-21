#ifndef IGA_STRONG_FLUID_STRUCTURE_COUPLING_HPP
#define IGA_STRONG_FLUID_STRUCTURE_COUPLING_HPP

// Bounded single-partition Dirichlet--Neumann fixed-point coordinator.  It
// owns coupling control, never a fluid or structure numerical state.  Runtime
// adapters supply snapshots and the exact lifecycle calls named below.
#include "DynamicWeightedAitkenRelaxation.hpp"
#include "FsiCouplingEdge.hpp"
#include "FsiDomainRuntime.hpp"
#include "PhaseProfile.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace iga {

namespace strong_fsi_detail {
template <class Runtime, class = void> struct HasCommittedFluidDisplacement : std::false_type {};
template <class Runtime> struct HasCommittedFluidDisplacement<Runtime,
    std::void_t<decltype(std::declval<const Runtime&>().StrongCouplingCommittedDisplacementM())>> : std::true_type {};
template <class Runtime, class = void> struct HasTrialFlowDiagnostics : std::false_type {};
template <class Runtime> struct HasTrialFlowDiagnostics<Runtime,
	std::void_t<decltype(std::declval<const Runtime&>().TrialFlowDiagnostics())>> : std::true_type {};
template <class Runtime, class = void> struct HasTrialCouplingAcceptanceDiagnostics : std::false_type {};
template <class Runtime> struct HasTrialCouplingAcceptanceDiagnostics<Runtime,
	std::void_t<decltype(std::declval<const Runtime&>().TrialCouplingAcceptanceDiagnostics())>> : std::true_type {};
}

class StrongFluidStructureCouplingAccess {
private:
	template <class FluidRuntime, class StructureRuntime>
	friend class StrongFluidStructureCoupling;
	template <class Runtime>
	static void PrevalidateFinalize(Runtime& runtime)
	{ runtime.CoordinatorRequireFinalizeAllowed(); }
	template <class Runtime>
	static void FinalizeNoexcept(Runtime& runtime) noexcept
	{
		static_assert(noexcept(runtime.CoordinatorFinalizeCommitNoexcept()),
			"strong FSI runtime coordinator finalization must be noexcept");
		runtime.CoordinatorFinalizeCommitNoexcept();
	}
};

struct StrongCouplingStructureSnapshot {
	std::vector<double> displacement_m;
	std::vector<double> velocity_m_per_s;
	std::vector<std::array<double, 3>> immutable_reference_normals;
	std::vector<bool> clamped;
	std::string committed_state_identity_sha256;
	std::string model_identity_sha256;
};

struct StrongCouplingAcceptanceDiagnostics {
	double residual = 0.0;
	double interface_power_defect_w = 0.0;
};

struct StrongFluidStructureCouplingOptions {
	std::uint64_t maximum_iterations = 32;
	double absolute_displacement_tolerance_m = 1.e-10;
	double relative_displacement_tolerance = 1.e-6;
	double reference_displacement_scale_m = 1.e-6;
	// Infinity preserves legacy displacement-only behavior. Finite values make
	// each quantity an independent acceptance gate.
	double absolute_nodal_force_residual_tolerance_n = std::numeric_limits<double>::infinity();
	double absolute_traction_residual_tolerance_pa = std::numeric_limits<double>::infinity();
	double fluid_residual_tolerance = std::numeric_limits<double>::infinity();
	double structure_residual_tolerance_n = std::numeric_limits<double>::infinity();
	double interface_power_defect_tolerance_w = std::numeric_limits<double>::infinity();
	AitkenRelaxationControls aitken_controls{};
	// This is intentionally a narrow test seam: it runs after every allocating
	// result/preparation operation and before either runtime can commit.
	bool fail_before_finalize_for_testing = false;
};

struct StrongFluidStructureIterationDiagnostics {
	std::uint64_t iteration = 0;
	double area_weighted_rms_residual_m = 0.0;
	double max_residual_m = 0.0;
	double displacement_scale_m = 0.0;
	double convergence_threshold_m = 0.0;
	bool converged = false;
	std::optional<double> nodal_force_residual_n;
	std::optional<double> traction_residual_pa;
	std::optional<double> fluid_acceptance_residual;
	std::optional<double> structure_acceptance_residual_n;
	std::optional<double> interface_power_defect_w;
	bool displacement_gate_passed = false;
	bool traction_gate_passed = false;
	bool nodal_force_gate_passed = false;
	bool fluid_gate_passed = false;
	bool structure_gate_passed = false;
	bool interface_power_gate_passed = false;
	bool aitken_proposal_applied = false;
	std::optional<double> relaxation;
	std::optional<double> unclamped_relaxation;
	std::optional<AitkenRelaxationStatus> aitken_status;
	std::optional<double> aitken_global_numerator;
	std::optional<double> aitken_global_denominator;
	std::optional<double> aitken_residual_scale;
	bool aitken_had_previous_residual = false;
	std::string fluid_traction_identity_sha256;
	std::string raw_kinematics_identity_sha256;
	std::string current_kinematics_identity_sha256;
	std::string relaxed_kinematics_identity_sha256;
	std::string aitken_proposal_identity_sha256;
	std::string aitken_control_state_identity_sha256;
	// Present only when the concrete fluid adapter exposes real trial flow
	// diagnostics.  They make a production benchmark's aggregate solve work
	// auditable without changing the mock adapter contract.
	std::optional<std::uint64_t> fluid_nonlinear_iterations;
	std::optional<std::uint64_t> fluid_ksp_iterations;
	std::optional<double> fluid_residual_norm;
	std::optional<double> fluid_linear_relative_residual;
};

struct StrongFluidStructureCouplingResult {
	bool converged = false;
	std::uint64_t iterations = 0;
	std::string status;
	std::string identity_sha256;
	std::string accepted_fluid_traction_identity_sha256;
	std::string accepted_raw_kinematics_identity_sha256;
	std::string accepted_relaxed_kinematics_identity_sha256;
	std::string final_committed_fluid_traction_identity_sha256;
	std::string final_committed_structure_kinematics_identity_sha256;
	std::string final_committed_structure_state_identity_sha256;
	std::string final_committed_fluid_composition_identity_sha256;
	std::string options_identity_sha256;
	std::vector<StrongFluidStructureIterationDiagnostics> history;
};

// FluidRuntime and StructureRuntime are intentionally structural contracts.
// This keeps deterministic tests free of PETSc/the real immersed solve while
// the production adapters expose the same narrow operations.
template <class FluidRuntime, class StructureRuntime>
class StrongFluidStructureCoupling final {
public:
	StrongFluidStructureCoupling(const StrongFluidStructureCoupling&) = delete;
	StrongFluidStructureCoupling& operator=(const StrongFluidStructureCoupling&) = delete;
	StrongFluidStructureCoupling(StrongFluidStructureCoupling&&) = delete;
	StrongFluidStructureCoupling& operator=(StrongFluidStructureCoupling&&) = delete;

	StrongFluidStructureCoupling(FluidRuntime& fluid, StructureRuntime& structure,
		FsiCouplingEdge edge, DistributedSurfaceLayout layout,
		StrongFluidStructureCouplingOptions options = {})
		: fluid_(fluid), structure_(structure), edge_(std::move(edge)), layout_(std::move(layout)),
		  partition_identity_sha256_(BuildDistributedSurfacePartitionIdentitySha256(layout_)), options_(options),
		  aitken_(partition_identity_sha256_, layout_.owned_reference_lumped_areas_m2,
			layout_.owned_reference_lumped_areas_m2.empty() ? 0.0 : AreaTotal(), options_.aitken_controls)
	{
		ValidateConstruction();
	}

	const FsiCouplingEdge& Edge() const noexcept { return edge_; }
	const DistributedSurfaceLayout& Layout() const noexcept { return layout_; }
	const DynamicWeightedAitkenRelaxation& Aitken() const noexcept { return aitken_; }

	StrongFluidStructureCouplingResult Execute(const DomainStepContext& step)
	{
		PhaseScope coupling_phase(ProfilePhase::Coupling);
		static_assert(std::is_nothrow_move_constructible<StrongFluidStructureCouplingResult>::value,
			"strong FSI result return must not throw after paired finalization");
		if (active_) throw std::runtime_error("strong FSI coordinator already has an active macro step");
		step.Validate();
		active_ = true;
		aitken_.Reset();
		try {
			StrongFluidStructureCouplingResult result;
			std::optional<SurfaceTraction> previous_traction;
			const StrongCouplingStructureSnapshot committed = structure_.StrongCouplingCommittedSnapshot();
			ValidateSnapshot(committed);
			SurfaceKinematics current = BuildPredictor(step, committed);
			fluid_.BeginMacroStep(step); structure_.BeginMacroStep(step);
			for (std::uint64_t iteration = 0; iteration < options_.maximum_iterations; ++iteration) {
				fluid_.BeginCouplingIteration(iteration, current.stamp, Envelope(step, iteration));
				fluid_.SetSurfaceKinematics(edge_.fluid.interface_id, current);
				fluid_.SolveFluidTrial();
				const SurfaceTraction traction = fluid_.GetSurfaceTraction(edge_.fluid.interface_id);
				structure_.BeginCouplingIteration(iteration, traction.stamp, Envelope(step, iteration));
				structure_.SetSurfaceTraction(edge_.structure.interface_id, traction);
				structure_.SolveMembraneTrial();
				const SurfaceKinematics raw = structure_.GetSurfaceKinematics(edge_.structure.interface_id);
				const std::vector<double> residual = ScalarResidual(raw, current, committed.immutable_reference_normals);
				const double raw_scale = ScalarMagnitude(raw, committed.immutable_reference_normals);
				const double current_scale = ScalarMagnitude(current, committed.immutable_reference_normals);
				const double scale = std::max(options_.reference_displacement_scale_m, std::max(raw_scale, current_scale));
				const double rms = WeightedRms(residual), maximum = MaxMagnitude(residual);
				StrongFluidStructureIterationDiagnostics diagnostic;
				diagnostic.iteration = iteration; diagnostic.area_weighted_rms_residual_m = rms;
				diagnostic.max_residual_m = maximum; diagnostic.displacement_scale_m = scale;
				diagnostic.convergence_threshold_m = options_.absolute_displacement_tolerance_m
					+options_.relative_displacement_tolerance*scale;
				diagnostic.displacement_gate_passed = rms <= diagnostic.convergence_threshold_m;
				if (std::isinf(options_.absolute_nodal_force_residual_tolerance_n)) {
					diagnostic.nodal_force_gate_passed = true;
				} else if (previous_traction) {
					diagnostic.nodal_force_residual_n = NodalForceResidual(traction,*previous_traction);
					diagnostic.nodal_force_gate_passed = *diagnostic.nodal_force_residual_n
						<= options_.absolute_nodal_force_residual_tolerance_n;
				}
				if(std::isinf(options_.absolute_traction_residual_tolerance_pa)){
					diagnostic.traction_gate_passed=true;
				}else if(previous_traction){
					diagnostic.traction_residual_pa=TractionResidual(traction,*previous_traction);
					diagnostic.traction_gate_passed=*diagnostic.traction_residual_pa
						<=options_.absolute_traction_residual_tolerance_pa;
				}
				diagnostic.fluid_gate_passed = std::isinf(options_.fluid_residual_tolerance);
				diagnostic.structure_gate_passed = std::isinf(options_.structure_residual_tolerance_n);
				diagnostic.interface_power_gate_passed = std::isinf(options_.interface_power_defect_tolerance_w);
				if constexpr (strong_fsi_detail::HasTrialCouplingAcceptanceDiagnostics<FluidRuntime>::value) {
					const auto gates=fluid_.TrialCouplingAcceptanceDiagnostics();
					if(!std::isfinite(gates.residual)||gates.residual<0.0
						||!std::isfinite(gates.interface_power_defect_w)||gates.interface_power_defect_w<0.0)
						throw std::runtime_error("strong FSI fluid acceptance diagnostics are invalid");
					diagnostic.fluid_acceptance_residual=gates.residual;
					diagnostic.interface_power_defect_w=gates.interface_power_defect_w;
					diagnostic.fluid_gate_passed=gates.residual<=options_.fluid_residual_tolerance;
					diagnostic.interface_power_gate_passed=gates.interface_power_defect_w
						<=options_.interface_power_defect_tolerance_w;
				}
				if constexpr (strong_fsi_detail::HasTrialCouplingAcceptanceDiagnostics<StructureRuntime>::value) {
					const auto gates=structure_.TrialCouplingAcceptanceDiagnostics();
					if(!std::isfinite(gates.residual)||gates.residual<0.0)
						throw std::runtime_error("strong FSI structure acceptance diagnostics are invalid");
					diagnostic.structure_acceptance_residual_n=gates.residual;
					diagnostic.structure_gate_passed=gates.residual<=options_.structure_residual_tolerance_n;
				}
				diagnostic.fluid_traction_identity_sha256 = BuildSurfaceTractionIdentitySha256(traction, layout_);
				diagnostic.raw_kinematics_identity_sha256 = BuildSurfaceKinematicsIdentitySha256(raw, layout_);
				diagnostic.current_kinematics_identity_sha256 = BuildSurfaceKinematicsIdentitySha256(current, layout_);
				if constexpr (strong_fsi_detail::HasTrialFlowDiagnostics<FluidRuntime>::value) {
					const auto flow_diagnostics = fluid_.TrialFlowDiagnostics();
					diagnostic.fluid_nonlinear_iterations = static_cast<std::uint64_t>(flow_diagnostics.nonlinear_iterations);
					diagnostic.fluid_ksp_iterations = static_cast<std::uint64_t>(flow_diagnostics.ksp_iterations);
					diagnostic.fluid_residual_norm = flow_diagnostics.residual_norm;
					diagnostic.fluid_linear_relative_residual = flow_diagnostics.true_linear_relative_residual;
				}
                if(CurrentPhaseProfile().Enabled()) {
                    std::printf("strong_fsi_iteration step=%llu iteration=%llu residual_rms_m=%.17g threshold_m=%.17g converged=%d\n",
                        static_cast<unsigned long long>(step.step_index),static_cast<unsigned long long>(iteration),
                        rms,diagnostic.convergence_threshold_m,
						diagnostic.displacement_gate_passed&&diagnostic.traction_gate_passed
						&&diagnostic.nodal_force_gate_passed
						&&diagnostic.fluid_gate_passed&&diagnostic.structure_gate_passed
						&&diagnostic.interface_power_gate_passed ? 1 : 0);
                    std::fflush(stdout);
                }
				if (diagnostic.displacement_gate_passed&&diagnostic.traction_gate_passed
					&&diagnostic.nodal_force_gate_passed
					&&diagnostic.fluid_gate_passed&&diagnostic.structure_gate_passed
					&&diagnostic.interface_power_gate_passed) {
					diagnostic.relaxed_kinematics_identity_sha256 = BuildSurfaceKinematicsIdentitySha256(current, layout_);
					diagnostic.converged = true; result.history.push_back(std::move(diagnostic));
					result.converged = true; result.iterations = iteration+1; result.status = "converged";
					const auto& accepted = result.history.back();
					result.accepted_fluid_traction_identity_sha256 = accepted.fluid_traction_identity_sha256;
					result.accepted_raw_kinematics_identity_sha256 = accepted.raw_kinematics_identity_sha256;
					result.accepted_relaxed_kinematics_identity_sha256 = accepted.relaxed_kinematics_identity_sha256;
					PrepareAndFinalizeAtomically(step, result);
					active_ = false;
					return result;
				}
				previous_traction=traction;
				fluid_.RejectCouplingIteration(); structure_.RejectCouplingIteration();
				const auto proposal = aitken_.Propose(ScalarDisplacement(current, committed.immutable_reference_normals),
					residual, options_.reference_displacement_scale_m, partition_identity_sha256_);
				aitken_.AcceptApplied(proposal, residual, proposal.relaxation_factor, partition_identity_sha256_);
				diagnostic.aitken_proposal_applied = true; diagnostic.relaxation = proposal.relaxation_factor;
				diagnostic.unclamped_relaxation = proposal.unclamped_relaxation_factor;
				diagnostic.aitken_status = proposal.status; diagnostic.aitken_global_numerator = proposal.global_numerator;
				diagnostic.aitken_global_denominator = proposal.global_denominator; diagnostic.aitken_residual_scale = proposal.residual_scale;
				diagnostic.aitken_had_previous_residual = proposal.has_previous_residual;
				diagnostic.aitken_proposal_identity_sha256 = proposal.proposal_identity_sha256;
				diagnostic.aitken_control_state_identity_sha256 = proposal.control_state_identity_sha256;
				current = BuildRelaxed(step, iteration+1, committed, raw, proposal);
				diagnostic.relaxed_kinematics_identity_sha256 = BuildSurfaceKinematicsIdentitySha256(current, layout_);
				result.history.push_back(std::move(diagnostic));
			}
			throw std::runtime_error("strong FSI fixed-point iteration did not converge");
		}
		catch (...) { fluid_.AbortStep(); structure_.AbortStep(); aitken_.Reset(); active_ = false; throw; }
	}

private:
	static double Dot(const std::array<double, 3>& a, const std::array<double, 3>& b) noexcept
	{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
	double AreaTotal() const
	{ double total = 0.0; for (const double area : layout_.owned_reference_lumped_areas_m2) total += area; return total; }
	void ValidateConstruction() const
	{
		ValidateFsiCouplingEdge(edge_); ValidateDistributedSurfaceLayout(layout_);
		if (layout_.partition_count != 1 || layout_.partition_rank != 0)
			throw std::runtime_error("strong FSI coordinator is bounded to one partition");
		if (!(fluid_.StrongCouplingEdge() == edge_) || !(structure_.StrongCouplingEdge() == edge_)
			|| fluid_.StrongCouplingLayout().layout_identity_sha256 != layout_.layout_identity_sha256
			|| structure_.StrongCouplingLayout().layout_identity_sha256 != layout_.layout_identity_sha256
			|| BuildDistributedSurfacePartitionIdentitySha256(fluid_.StrongCouplingLayout()) != partition_identity_sha256_
			|| BuildDistributedSurfacePartitionIdentitySha256(structure_.StrongCouplingLayout()) != partition_identity_sha256_)
			throw std::runtime_error("strong FSI coordinator does not bind exact edge/layout/partition runtimes");
		if (options_.maximum_iterations == 0 || !(options_.absolute_displacement_tolerance_m >= 0.0)
			|| !(options_.relative_displacement_tolerance >= 0.0)
			|| !(options_.reference_displacement_scale_m > 0.0)
			|| !std::isfinite(options_.absolute_displacement_tolerance_m)
			|| !std::isfinite(options_.relative_displacement_tolerance)
			|| !std::isfinite(options_.reference_displacement_scale_m)
			|| !(options_.absolute_nodal_force_residual_tolerance_n>=0.0)
			|| !(options_.absolute_traction_residual_tolerance_pa>=0.0)
			|| !(options_.fluid_residual_tolerance>=0.0)
			|| !(options_.structure_residual_tolerance_n>=0.0)
			|| !(options_.interface_power_defect_tolerance_w>=0.0))
			throw std::runtime_error("strong FSI convergence options are invalid");
		if(!std::isinf(options_.fluid_residual_tolerance)
			&&!strong_fsi_detail::HasTrialCouplingAcceptanceDiagnostics<FluidRuntime>::value)
			throw std::runtime_error("strong FSI finite fluid gate requires runtime diagnostics");
		if(!std::isinf(options_.interface_power_defect_tolerance_w)
			&&!strong_fsi_detail::HasTrialCouplingAcceptanceDiagnostics<FluidRuntime>::value)
			throw std::runtime_error("strong FSI finite power gate requires fluid diagnostics");
		if(!std::isinf(options_.structure_residual_tolerance_n)
			&&!strong_fsi_detail::HasTrialCouplingAcceptanceDiagnostics<StructureRuntime>::value)
			throw std::runtime_error("strong FSI finite structure gate requires runtime diagnostics");
	}
	void ValidateSnapshot(const StrongCouplingStructureSnapshot& snapshot) const
	{
		const std::size_t count = layout_.owned_global_node_ids.size();
		if (snapshot.displacement_m.size() != count || snapshot.velocity_m_per_s.size() != count
			|| snapshot.immutable_reference_normals.size() != count || snapshot.clamped.size() != count
			|| !IsLowercaseSha256(snapshot.committed_state_identity_sha256)
			|| !IsLowercaseSha256(snapshot.model_identity_sha256)) throw std::runtime_error("strong FSI structure snapshot is invalid");
		for (std::size_t i = 0; i < count; ++i) {
			const double length = std::sqrt(Dot(snapshot.immutable_reference_normals[i], snapshot.immutable_reference_normals[i]));
			if (!std::isfinite(snapshot.displacement_m[i]) || !std::isfinite(snapshot.velocity_m_per_s[i])
				|| !std::isfinite(length) || std::abs(length-1.0) > 128.0*std::numeric_limits<double>::epsilon()
				|| (snapshot.clamped[i] && (snapshot.displacement_m[i] != 0.0 || snapshot.velocity_m_per_s[i] != 0.0)))
				throw std::runtime_error("strong FSI structure snapshot violates scalar/clamp contract");
		}
	}
	SurfaceFieldStampEnvelope Envelope(const DomainStepContext& step, std::uint64_t iteration) const
	{ return {step.start_time_s+step.dt_s, static_cast<std::uint64_t>(step.step_index), iteration,
		layout_.reference_mesh_identity_sha256, layout_.layout_identity_sha256, partition_identity_sha256_}; }
	std::string KinematicsProducerIdentity(const char* kind, const DomainStepContext& step, std::uint64_t iteration,
		const StrongCouplingStructureSnapshot& state, const std::string& raw_identity,
		const DynamicWeightedAitkenProposal* proposal, const std::vector<double>& scalar) const
	{
		Sha256 hash; distributed_surface_detail::AppendString(hash, kind);
		distributed_surface_detail::AppendString(hash, BuildFsiCouplingEdgeIdentitySha256(edge_));
		distributed_surface_detail::AppendString(hash, partition_identity_sha256_);
		distributed_surface_detail::AppendString(hash, state.committed_state_identity_sha256);
		distributed_surface_detail::AppendString(hash, state.model_identity_sha256);
		distributed_surface_detail::AppendString(hash, raw_identity); hash.AppendLittleEndian64(static_cast<std::uint64_t>(step.step_index));
		hash.AppendNormalizedDouble(step.start_time_s); hash.AppendNormalizedDouble(step.dt_s); hash.AppendLittleEndian64(iteration);
		if (proposal) { distributed_surface_detail::AppendString(hash, proposal->proposal_identity_sha256); distributed_surface_detail::AppendString(hash, proposal->control_state_identity_sha256); hash.AppendNormalizedDouble(proposal->relaxation_factor); }
		for (const double value : scalar) hash.AppendNormalizedDouble(value);
		return hash.Hex();
	}
    void ApplyFluidKinematicHistory(SurfaceKinematics& value, const DomainStepContext& step) const
    {
        if constexpr (strong_fsi_detail::HasCommittedFluidDisplacement<FluidRuntime>::value) {
            // The fluid commits the relaxed iterate; the membrane commits its
            // raw solution. A finite accepted residual must not change the
            // fluid's backward-Euler displacement/velocity history next step.
            const auto history=fluid_.StrongCouplingCommittedDisplacementM();
            if(history.size()!=value.displacement_m.size())
                throw std::runtime_error("committed fluid interface history has incompatible size");
            Sha256 hash;
            distributed_surface_detail::AppendString(hash,"StrongFsiFluidKinematicHistory/v1");
            distributed_surface_detail::AppendString(hash,value.stamp.producer_state_identity_sha256);
            for(std::size_t i=0;i<history.size();++i) for(int c=0;c<3;++c) {
                if(!std::isfinite(history[i][c]))
                    throw std::runtime_error("committed fluid interface history is nonfinite");
                value.velocity_m_per_s[i][c]=(value.displacement_m[i][c]-history[i][c])/step.dt_s;
                hash.AppendNormalizedDouble(history[i][c]);
                hash.AppendNormalizedDouble(value.velocity_m_per_s[i][c]);
            }
            value.stamp.producer_state_identity_sha256=hash.Hex();
        }
    }
	SurfaceKinematics BuildPredictor(const DomainStepContext& step, const StrongCouplingStructureSnapshot& state) const
	{
		SurfaceKinematics result; result.interface = edge_.structure; const auto envelope = Envelope(step, 0);
		result.stamp = {envelope.time_s,envelope.step,envelope.coupling_iteration,envelope.reference_mesh_identity_sha256,envelope.layout_identity_sha256,envelope.partition_identity_sha256,{}};
		std::vector<double> scalar(state.displacement_m.size()); result.displacement_m.resize(scalar.size()); result.velocity_m_per_s.resize(scalar.size());
		for (std::size_t i=0;i<scalar.size();++i) { scalar[i]=state.clamped[i] ? 0.0 : state.displacement_m[i]+step.dt_s*state.velocity_m_per_s[i]; const double v=state.clamped[i]?0.0:state.velocity_m_per_s[i]; for(int c=0;c<3;++c){result.displacement_m[i][c]=scalar[i]*state.immutable_reference_normals[i][c]; result.velocity_m_per_s[i][c]=v*state.immutable_reference_normals[i][c];} }
		result.stamp.producer_state_identity_sha256=KinematicsProducerIdentity("StrongFsiPredictor/v1",step,0,state,"",nullptr,scalar); ApplyFluidKinematicHistory(result,step); ValidateSurfaceKinematics(result,layout_); return result;
	}
	std::vector<double> ScalarDisplacement(const SurfaceKinematics& value, const std::vector<std::array<double,3>>& normals) const
	{ ValidateSurfaceKinematics(value,layout_); std::vector<double> r(value.displacement_m.size()); for(std::size_t i=0;i<r.size();++i) r[i]=Dot(value.displacement_m[i],normals[i]); return r; }
	std::vector<double> ScalarResidual(const SurfaceKinematics& raw, const SurfaceKinematics& current, const std::vector<std::array<double,3>>& normals) const
	{ const auto a=ScalarDisplacement(raw,normals), b=ScalarDisplacement(current,normals); std::vector<double> r(a.size()); for(std::size_t i=0;i<r.size();++i) r[i]=a[i]-b[i]; return r; }
	double ScalarMagnitude(const SurfaceKinematics& value, const std::vector<std::array<double,3>>& normals) const
	{ double r=0.; for(const double x:ScalarDisplacement(value,normals)) r=std::max(r,std::abs(x)); return r; }
	double WeightedRms(const std::vector<double>& values) const
	{ double sum=0.; for(std::size_t i=0;i<values.size();++i) sum+=layout_.owned_reference_lumped_areas_m2[i]*values[i]*values[i]; return std::sqrt(sum/AreaTotal()); }
	double NodalForceResidual(const SurfaceTraction& first,const SurfaceTraction& second)const
	{
		ValidateSurfaceTraction(first,layout_);ValidateSurfaceTraction(second,layout_);
		double sum=0.0;const std::size_t count=first.consistent_nodal_force_n.size();
		for(std::size_t node=0;node<count;++node)for(int component=0;component<3;++component){
			const double difference=first.consistent_nodal_force_n[node][component]
				-second.consistent_nodal_force_n[node][component];sum+=difference*difference;}
		return std::sqrt(sum/static_cast<double>(count));
	}
	double TractionResidual(const SurfaceTraction& first,const SurfaceTraction& second)const
	{
		ValidateSurfaceTraction(first,layout_);ValidateSurfaceTraction(second,layout_);
		double sum=0.0;
		for(std::size_t node=0;node<first.traction_on_structure_pa.size();++node){
			double squared=0.0;for(int component=0;component<3;++component){
				const double difference=first.traction_on_structure_pa[node][component]
					-second.traction_on_structure_pa[node][component];squared+=difference*difference;}
			sum+=layout_.owned_reference_lumped_areas_m2[node]*squared;}
		return std::sqrt(sum/AreaTotal());
	}
	static double MaxMagnitude(const std::vector<double>& values) { double r=0.; for(double x:values) r=std::max(r,std::abs(x)); return r; }
	SurfaceKinematics BuildRelaxed(const DomainStepContext& step, std::uint64_t iteration, const StrongCouplingStructureSnapshot& state,
		const SurfaceKinematics& raw, const DynamicWeightedAitkenProposal& proposal) const
	{
		SurfaceKinematics result; result.interface=edge_.structure; const auto envelope=Envelope(step,iteration); result.stamp={envelope.time_s,envelope.step,envelope.coupling_iteration,envelope.reference_mesh_identity_sha256,envelope.layout_identity_sha256,envelope.partition_identity_sha256,{}};
		result.displacement_m.resize(proposal.next.size()); result.velocity_m_per_s.resize(proposal.next.size()); const std::string raw_identity=BuildSurfaceKinematicsIdentitySha256(raw,layout_);
		for(std::size_t i=0;i<proposal.next.size();++i){ const double d=state.clamped[i]?0.:proposal.next[i]; const double v=state.clamped[i]?0.:(d-state.displacement_m[i])/step.dt_s; for(int c=0;c<3;++c){result.displacement_m[i][c]=d*state.immutable_reference_normals[i][c]; result.velocity_m_per_s[i][c]=v*state.immutable_reference_normals[i][c];} }
		result.stamp.producer_state_identity_sha256=KinematicsProducerIdentity("StrongFsiRelaxedIterate/v1",step,iteration,state,raw_identity,&proposal,proposal.next); ApplyFluidKinematicHistory(result,step); ValidateSurfaceKinematics(result,layout_); return result;
	}
	void PrepareAndFinalizeAtomically(const DomainStepContext& step, StrongFluidStructureCouplingResult& result)
	{
		try {
			fluid_.PrepareCommitStep(); structure_.PrepareCommitStep();
			// Finish all allocation/string hashing and obtain the exact prepared
			// identities while both runtimes can still be aborted.
			result.final_committed_fluid_traction_identity_sha256 = fluid_.CoordinatorPreparedCommittedTractionIdentitySha256();
			result.final_committed_structure_kinematics_identity_sha256 = structure_.CoordinatorPreparedCommittedKinematicsIdentitySha256();
			result.final_committed_structure_state_identity_sha256 = structure_.CoordinatorPreparedCommittedStateIdentitySha256();
			result.final_committed_fluid_composition_identity_sha256 = fluid_.CoordinatorPreparedCommittedCompositionIdentitySha256();
			result.options_identity_sha256 = OptionsIdentity();
			result.identity_sha256 = ResultIdentity(step, result);
			if (options_.fail_before_finalize_for_testing) {
				options_.fail_before_finalize_for_testing = false;
				throw std::runtime_error("injected strong FSI pre-finalize failure");
			}
			StrongFluidStructureCouplingAccess::PrevalidateFinalize(fluid_);
			StrongFluidStructureCouplingAccess::PrevalidateFinalize(structure_);
		}
		catch (...) { fluid_.AbortStep(); structure_.AbortStep(); throw; }
		// No allocation, validation, catchable call, or potentially throwing move
		// is permitted below this line: the first handoff is irreversible.
		StrongFluidStructureCouplingAccess::FinalizeNoexcept(fluid_); StrongFluidStructureCouplingAccess::FinalizeNoexcept(structure_);
	}
	std::string OptionsIdentity() const
	{
		Sha256 hash; distributed_surface_detail::AppendString(hash,"StrongFluidStructureCouplingOptions/v3");
		hash.AppendLittleEndian64(options_.maximum_iterations); hash.AppendNormalizedDouble(options_.absolute_displacement_tolerance_m);
		hash.AppendNormalizedDouble(options_.relative_displacement_tolerance); hash.AppendNormalizedDouble(options_.reference_displacement_scale_m);
		AppendGate(hash,options_.absolute_nodal_force_residual_tolerance_n);
		AppendGate(hash,options_.absolute_traction_residual_tolerance_pa);
		AppendGate(hash,options_.fluid_residual_tolerance);
		AppendGate(hash,options_.structure_residual_tolerance_n);
		AppendGate(hash,options_.interface_power_defect_tolerance_w);
		hash.AppendNormalizedDouble(options_.aitken_controls.initial_relaxation); hash.AppendNormalizedDouble(options_.aitken_controls.minimum_relaxation);
		hash.AppendNormalizedDouble(options_.aitken_controls.maximum_relaxation); hash.AppendNormalizedDouble(options_.aitken_controls.scaled_difference_threshold);
		return hash.Hex();
	}
	std::string ResultIdentity(const DomainStepContext& step, const StrongFluidStructureCouplingResult& result) const
	{ Sha256 hash; distributed_surface_detail::AppendString(hash,"StrongFluidStructureCouplingResult/v5"); distributed_surface_detail::AppendString(hash,BuildFsiCouplingEdgeIdentitySha256(edge_)); distributed_surface_detail::AppendString(hash,partition_identity_sha256_); distributed_surface_detail::AppendString(hash,result.options_identity_sha256); hash.AppendLittleEndian64(static_cast<std::uint64_t>(step.step_index)); hash.AppendNormalizedDouble(step.start_time_s); hash.AppendNormalizedDouble(step.dt_s); hash.AppendLittleEndian64(result.iterations); hash.AppendLittleEndian32(result.converged ? 1U : 0U); distributed_surface_detail::AppendString(hash,result.status); for(const auto& d:result.history){hash.AppendLittleEndian64(d.iteration);hash.AppendNormalizedDouble(d.area_weighted_rms_residual_m);hash.AppendNormalizedDouble(d.max_residual_m);hash.AppendNormalizedDouble(d.displacement_scale_m);hash.AppendNormalizedDouble(d.convergence_threshold_m);AppendOptional(hash,d.nodal_force_residual_n);AppendOptional(hash,d.traction_residual_pa);AppendOptional(hash,d.fluid_acceptance_residual);AppendOptional(hash,d.structure_acceptance_residual_n);AppendOptional(hash,d.interface_power_defect_w);hash.AppendLittleEndian32(d.displacement_gate_passed?1U:0U);hash.AppendLittleEndian32(d.traction_gate_passed?1U:0U);hash.AppendLittleEndian32(d.nodal_force_gate_passed?1U:0U);hash.AppendLittleEndian32(d.fluid_gate_passed?1U:0U);hash.AppendLittleEndian32(d.structure_gate_passed?1U:0U);hash.AppendLittleEndian32(d.interface_power_gate_passed?1U:0U);hash.AppendLittleEndian32(d.converged ? 1U : 0U);hash.AppendLittleEndian32(d.aitken_proposal_applied ? 1U : 0U); AppendOptional(hash,d.relaxation);AppendOptional(hash,d.unclamped_relaxation);AppendOptionalStatus(hash,d.aitken_status);AppendOptional(hash,d.aitken_global_numerator);AppendOptional(hash,d.aitken_global_denominator);AppendOptional(hash,d.aitken_residual_scale);hash.AppendLittleEndian32(d.aitken_had_previous_residual ? 1U : 0U); AppendOptionalUint64(hash,d.fluid_nonlinear_iterations);AppendOptionalUint64(hash,d.fluid_ksp_iterations);AppendOptional(hash,d.fluid_residual_norm);AppendOptional(hash,d.fluid_linear_relative_residual);distributed_surface_detail::AppendString(hash,d.fluid_traction_identity_sha256);distributed_surface_detail::AppendString(hash,d.raw_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,d.current_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,d.relaxed_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,d.aitken_proposal_identity_sha256);distributed_surface_detail::AppendString(hash,d.aitken_control_state_identity_sha256);} distributed_surface_detail::AppendString(hash,result.accepted_fluid_traction_identity_sha256);distributed_surface_detail::AppendString(hash,result.accepted_raw_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,result.accepted_relaxed_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,result.final_committed_fluid_traction_identity_sha256);distributed_surface_detail::AppendString(hash,result.final_committed_structure_kinematics_identity_sha256);distributed_surface_detail::AppendString(hash,result.final_committed_structure_state_identity_sha256);distributed_surface_detail::AppendString(hash,result.final_committed_fluid_composition_identity_sha256); return hash.Hex(); }
	static void AppendOptional(Sha256& hash, const std::optional<double>& value)
	{ hash.AppendLittleEndian32(value ? 1U : 0U); if (value) hash.AppendNormalizedDouble(*value); }
	static void AppendGate(Sha256& hash,double value)
	{hash.AppendLittleEndian32(std::isinf(value)?0U:1U);if(!std::isinf(value))hash.AppendNormalizedDouble(value);}
	static void AppendOptionalUint64(Sha256& hash, const std::optional<std::uint64_t>& value)
	{ hash.AppendLittleEndian32(value ? 1U : 0U); if (value) hash.AppendLittleEndian64(*value); }
	static void AppendOptionalStatus(Sha256& hash, const std::optional<AitkenRelaxationStatus>& value)
	{ hash.AppendLittleEndian32(value ? 1U : 0U); if (value) hash.AppendLittleEndian32(static_cast<std::uint32_t>(*value)); }

	FluidRuntime& fluid_; StructureRuntime& structure_; FsiCouplingEdge edge_; DistributedSurfaceLayout layout_; std::string partition_identity_sha256_; StrongFluidStructureCouplingOptions options_; DynamicWeightedAitkenRelaxation aitken_; bool active_=false;
};

} // namespace iga

#endif

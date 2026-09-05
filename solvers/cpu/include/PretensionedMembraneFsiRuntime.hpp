#ifndef IGA_PRETENSIONED_MEMBRANE_FSI_RUNTIME_HPP
#define IGA_PRETENSIONED_MEMBRANE_FSI_RUNTIME_HPP

// Rank-local structure-side FSI adapter.  It deliberately contains neither
// MPI transport nor a coupling coordinator: a caller supplies the exact field
// stamps for each already-scheduled exchange.
#include "FsiDomainRuntime.hpp"
#include "PretensionedMembrane.hpp"

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

class PretensionedMembraneFsiRuntime final : public FsiStructureDomainRuntime {
public:
	PretensionedMembraneFsiRuntime(const PretensionedMembraneFsiRuntime&) = delete;
	PretensionedMembraneFsiRuntime& operator=(const PretensionedMembraneFsiRuntime&) = delete;
	PretensionedMembraneFsiRuntime(PretensionedMembraneFsiRuntime&&) = delete;
	PretensionedMembraneFsiRuntime& operator=(PretensionedMembraneFsiRuntime&&) = delete;

	PretensionedMembraneFsiRuntime(std::string domain_id, std::string subsystem_id,
		FsiCouplingEdge edge, DistributedSurfaceInterface structure_surface,
		DistributedSurfaceInterface fluid_surface, DistributedSurfaceLayout structure_layout,
		DistributedSurfaceLayout fluid_layout, PretensionedMembraneMaterial material,
		std::vector<std::uint64_t> clamped_global_node_ids,
		PretensionedMembraneOptions options = PretensionedMembraneOptions{},
		PretensionedMembraneState initial_state = PretensionedMembraneState{})
		: domain_id_(std::move(domain_id)), subsystem_id_(std::move(subsystem_id)),
		  edge_(std::move(edge)), catalog_{std::move(structure_surface)},
		  fluid_surface_(std::move(fluid_surface)), structure_layout_(std::move(structure_layout)),
		  fluid_layout_(std::move(fluid_layout)), membrane_(structure_layout_, catalog_.front(),
			fluid_surface_.id, material, std::move(clamped_global_node_ids), options,
			std::move(initial_state)), lifecycle_(domain_id_, subsystem_id_, edge_, catalog_.front(),
			fluid_surface_, structure_layout_, fluid_layout_)
	{
		// The membrane and lifecycle repeat the important checks independently:
		// neither numerical ownership nor publication availability can be bound
		// to a different endpoint/layout than the other.
		if (!(catalog_.front().id == edge_.structure) || !(fluid_surface_.id == edge_.fluid))
			throw std::runtime_error("membrane FSI runtime surfaces do not match the exact edge endpoints");
		if (catalog_.front().id.domain_id != domain_id_
			|| catalog_.front().id.subsystem_id != subsystem_id_)
			throw std::runtime_error("membrane FSI runtime local surface does not match its domain/subsystem");
		if (structure_layout_.partition_count != 1 || structure_layout_.partition_rank != 0)
			throw std::runtime_error("membrane FSI runtime inherits the single-partition restriction");
	}

	const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept override
	{ return catalog_; }
	const FsiTrialLifecycle& Lifecycle() const noexcept { return lifecycle_; }
	PretensionedMembraneState CommittedStateSnapshot() const { return membrane_.CommittedState(); }
	std::string ModelIdentitySha256() const { return membrane_.ModelIdentitySha256(); }
	std::string CommittedStateIdentitySha256() const { return membrane_.CommittedStateIdentitySha256(); }

	void BeginMacroStep(const DomainStepContext& step)
	{
		if (trial_.has_value()) throw std::runtime_error("membrane FSI runtime retains an unexpected trial");
		lifecycle_.BeginStep(step);
	}
	void BeginStep(const DomainStepContext& step) { BeginMacroStep(step); }

	void BeginCouplingIteration(std::uint64_t iteration,
		const SurfaceFieldStamp& expected_traction, const SurfaceFieldStampEnvelope& output_envelope)
	{
		if (trial_.has_value()) throw std::runtime_error("membrane FSI runtime retains an unexpected trial");
		// Snapshot all caller-owned fields before changing lifecycle availability.
		// FsiTrialLifecycle repeats the authoritative context validation and its
		// own transition is transactional.
		SurfaceFieldStamp traction_copy = expected_traction;
		SurfaceFieldStampEnvelope envelope_copy = output_envelope;
		lifecycle_.BeginIteration(iteration, traction_copy, envelope_copy);
		using std::swap;
		swap(expected_traction_, traction_copy);
	}
	void BeginIteration(std::uint64_t iteration, const SurfaceFieldStamp& expected_traction,
		const SurfaceFieldStampEnvelope& output_envelope)
	{ BeginCouplingIteration(iteration, expected_traction, output_envelope); }

	void SetSurfaceTraction(const std::string& interface_id,
		const SurfaceTraction& traction) override
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown membrane structure surface interface");
		// Snapshot before changing availability.  No caller-owned mutable field is
		// retained, even if a later validation error is caught by the coordinator.
		SurfaceTraction candidate = traction;
		ValidateFsiStructureTractionInput(edge_, candidate, structure_layout_, expected_traction_);
		lifecycle_.MarkInput(candidate.interface, candidate.stamp);
		traction_ = std::move(candidate);
	}

	void SolveMembraneTrial()
	{
		lifecycle_.RequireSolveAllowed();
		if (!traction_.has_value()) throw std::runtime_error("membrane FSI runtime has no accepted traction snapshot");
		PretensionedMembraneTrial candidate = membrane_.SolveTrial(MembraneContext(), *traction_);
		try {
			lifecycle_.MarkSolved(candidate.kinematics.interface, candidate.kinematics.stamp);
		}
		catch (...) {
			membrane_.AbortTrial();
			throw std::runtime_error("membrane FSI kinematics stamp does not match the declared exact output");
		}
		trial_.emplace(std::move(candidate));
	}
	void SolveTrial() { SolveMembraneTrial(); }

	SurfaceKinematics GetSurfaceKinematics(const std::string& interface_id) const override
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown membrane structure surface interface");
		if (!trial_.has_value()) throw std::runtime_error("membrane FSI runtime has no solved trial");
		lifecycle_.RequireTrialOutput(trial_->kinematics.interface, trial_->kinematics.stamp);
		return trial_->kinematics;
	}

	SurfaceKinematics GetCommittedSurfaceKinematics(const std::string& interface_id,
		const SurfaceFieldStamp& stamp) const
	{
		if (interface_id != edge_.structure.interface_id)
			throw std::runtime_error("unknown membrane structure surface interface");
		lifecycle_.RequireCommittedOutput(catalog_.front().id, stamp);
		return committed_kinematics_;
	}

	void RejectCouplingIteration()
	{
		if (!lifecycle_.HasActiveStep()) throw std::runtime_error("membrane FSI runtime reject requires an active step");
		DiscardTrialNoexcept();
		lifecycle_.RejectIteration();
	}
	void RollbackTrial() { RejectCouplingIteration(); }

	void PrepareCommitStep()
	{
		if (!trial_.has_value()) throw std::runtime_error("membrane FSI runtime prepare requires a solved trial");
		lifecycle_.RequireTrialOutput(trial_->kinematics.interface, trial_->kinematics.stamp);
		// Copy the prospective public snapshot before either owner is prepared.
		SurfaceKinematics candidate_kinematics = trial_->kinematics;
		try {
			lifecycle_.PrepareCommit();
			PretensionedMembraneTrial prepared = membrane_.PrepareTrial(*trial_);
			trial_.emplace(std::move(prepared));
			prepared_kinematics_ = std::move(candidate_kinematics);
		}
		catch (...) {
			membrane_.AbortTrial();
			trial_.reset();
			traction_.reset();
			if (lifecycle_.HasActiveStep()) lifecycle_.RejectIteration();
			throw;
		}
	}
	void PrepareCommit() { PrepareCommitStep(); }

	void FinalizeCommitStep()
	{
		if (!trial_.has_value()) throw std::runtime_error("membrane FSI runtime finalize requires a prepared trial");
		// Both checks perform every remaining fallible integrity operation before
		// either owner changes.  The paired exchanges below are noexcept.
		membrane_.RequireFinalizeAllowed(*trial_);
		lifecycle_.RequireFinalizeAllowed();
		membrane_.FinalizePreparedTrialNoexcept(std::move(*trial_));
		using std::swap;
		swap(committed_kinematics_, prepared_kinematics_);
		lifecycle_.FinalizePreparedCommitNoexcept();
		trial_.reset();
		traction_.reset();
	}
	void FinalizeCommit() { FinalizeCommitStep(); }

	void AbortStep() noexcept
	{
		membrane_.AbortTrial();
		trial_.reset();
		traction_.reset();
		lifecycle_.AbortStep();
	}

private:
	PretensionedMembraneTrialContext MembraneContext() const
	{
		const auto& context = lifecycle_.Context();
		PretensionedMembraneTrialContext result;
		result.step = context.step;
		result.start_time_s = context.start_time_s;
		result.dt_s = context.dt_s;
		result.coupling_iteration = context.coupling_iteration;
		result.expected_traction_stamp = expected_traction_;
		return result;
	}
	void DiscardTrialNoexcept() noexcept
	{
		membrane_.AbortTrial();
		trial_.reset();
		traction_.reset();
	}

	std::string domain_id_, subsystem_id_;
	FsiCouplingEdge edge_;
	std::vector<DistributedSurfaceInterface> catalog_;
	DistributedSurfaceInterface fluid_surface_;
	DistributedSurfaceLayout structure_layout_, fluid_layout_;
	PretensionedMembrane membrane_;
	FsiTrialLifecycle lifecycle_;
	SurfaceFieldStamp expected_traction_;
	std::optional<SurfaceTraction> traction_;
	std::optional<PretensionedMembraneTrial> trial_;
	SurfaceKinematics prepared_kinematics_, committed_kinematics_;
};

} // namespace iga

#endif

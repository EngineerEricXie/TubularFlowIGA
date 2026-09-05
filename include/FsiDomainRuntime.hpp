#ifndef IGA_FSI_DOMAIN_RUNTIME_HPP
#define IGA_FSI_DOMAIN_RUNTIME_HPP

// Field-valued FSI capabilities deliberately sit beside CoupledDomainRuntime.
// FsiTrialLifecycle is a small owner that an FSI runtime embeds; it is not a
// second commit owner and does not alter scalar port/executor semantics.
#include "CoupledDomainRuntime.hpp"
#include "FsiCouplingEdge.hpp"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class FsiTrialPhase : std::uint8_t {
	Idle,
	StepActive,
	IterationAwaitingInput,
	InputReady,
	Solved,
	Prepared
};

struct FsiTrialContext {
	std::uint64_t step = 0;
	double start_time_s = 0.0;
	double dt_s = 0.0;
	std::uint64_t coupling_iteration = 0;

	double EndTime() const { return start_time_s+dt_s; }
};

inline void ValidateFsiTrialContext(const FsiTrialContext& context)
{
	if (!std::isfinite(context.start_time_s) || !std::isfinite(context.dt_s)
		|| !(context.dt_s > 0.0) || !std::isfinite(context.EndTime()))
		throw std::runtime_error("FSI trial context requires finite time and positive finite dt");
}

// Owns only field-trial availability. The host runtime remains responsible
// for numerical state and its existing BeginStep/SolveTrial/commit methods.
class FsiTrialLifecycle {
public:
	FsiTrialLifecycle(std::string domain_id, std::string subsystem_id,
		FsiCouplingEdge edge, DistributedSurfaceInterface local_surface,
		DistributedSurfaceInterface peer_surface, DistributedSurfaceLayout local_layout,
		DistributedSurfaceLayout peer_layout)
		: domain_id_(std::move(domain_id)), subsystem_id_(std::move(subsystem_id)),
		  edge_(std::move(edge)), local_surface_(std::move(local_surface)),
		  peer_surface_(std::move(peer_surface)), local_layout_(std::move(local_layout)),
		  peer_layout_(std::move(peer_layout))
	{
		if (domain_id_.empty() || subsystem_id_.empty())
			throw std::runtime_error("FSI lifecycle domain and subsystem IDs must be nonempty");
		const auto& fluid_surface = edge_.fluid == local_surface_.id ? local_surface_ : peer_surface_;
		const auto& structure_surface = edge_.structure == local_surface_.id ? local_surface_ : peer_surface_;
		const auto& fluid_layout = edge_.fluid == local_surface_.id ? local_layout_ : peer_layout_;
		const auto& structure_layout = edge_.structure == local_surface_.id ? local_layout_ : peer_layout_;
		ValidateFsiCouplingEdgeLayouts(edge_, fluid_surface, structure_surface,
			fluid_layout, structure_layout);
		if (local_surface_.id.domain_id != domain_id_
			|| local_surface_.id.subsystem_id != subsystem_id_)
			throw std::runtime_error("FSI lifecycle local surface does not bind its domain/subsystem");
		if (!(local_surface_.id == edge_.fluid) && !(local_surface_.id == edge_.structure))
			throw std::runtime_error("FSI lifecycle local surface is not an edge endpoint");
		if (peer_surface_.id == local_surface_.id)
			throw std::runtime_error("FSI lifecycle peer surface must be the other edge endpoint");
	}

	FsiTrialPhase Phase() const noexcept { return phase_; }
	bool HasActiveStep() const noexcept { return phase_ != FsiTrialPhase::Idle; }
	bool HasTrialOutput() const noexcept
	{
		return phase_ == FsiTrialPhase::Solved || phase_ == FsiTrialPhase::Prepared;
	}
	bool HasCommittedOutput() const noexcept { return has_committed_output_; }
	const FsiTrialContext& Context() const
	{
		if (!HasActiveStep()) throw std::runtime_error("FSI lifecycle has no active step context");
		return context_;
	}

	void BeginStep(const DomainStepContext& step)
	{
		step.Validate();
		if (phase_ != FsiTrialPhase::Idle)
			throw std::runtime_error("FSI lifecycle step is already active");
		context_.step = static_cast<std::uint64_t>(step.step_index);
		context_.start_time_s = step.start_time_s;
		context_.dt_s = step.dt_s;
		context_.coupling_iteration = 0;
		phase_ = FsiTrialPhase::StepActive;
	}

	void BeginIteration(std::uint64_t coupling_iteration,
		const SurfaceFieldStamp& expected_input, const SurfaceFieldStamp& expected_output)
	{
		if (phase_ != FsiTrialPhase::StepActive)
			throw std::runtime_error("FSI lifecycle iteration requires an active step with no prior trial");
		context_.coupling_iteration = coupling_iteration;
		ValidateFsiTrialContext(context_);
		ValidateExpectedStamp(expected_input, peer_layout_);
		ValidateExpectedStamp(expected_output, local_layout_);
		expected_input_ = expected_input;
		expected_output_ = expected_output;
		has_expected_stamps_ = true;
		has_input_ = false;
		has_output_ = false;
		phase_ = FsiTrialPhase::IterationAwaitingInput;
	}

	void MarkInput(const SurfaceInterfaceRef& producer, const SurfaceFieldStamp& stamp)
	{
		if (phase_ != FsiTrialPhase::IterationAwaitingInput || !has_expected_stamps_)
			throw std::runtime_error("FSI input is unavailable outside its exact awaiting iteration");
		if (!(producer == peer_surface_.id))
			throw std::runtime_error("FSI input belongs to the wrong producer endpoint");
		ValidateExactStamp(stamp, expected_input_, peer_layout_);
		has_input_ = true;
		phase_ = FsiTrialPhase::InputReady;
	}

	void RequireSolveAllowed() const
	{
		if (phase_ != FsiTrialPhase::InputReady || !has_input_)
			throw std::runtime_error("FSI solve requires its exact iteration input");
	}

	void MarkSolved(const SurfaceInterfaceRef& producer, const SurfaceFieldStamp& stamp)
	{
		RequireSolveAllowed();
		if (!(producer == local_surface_.id))
			throw std::runtime_error("FSI output belongs to the wrong local endpoint");
		ValidateExactStamp(stamp, expected_output_, local_layout_);
		has_output_ = true;
		phase_ = FsiTrialPhase::Solved;
	}

	void RequireTrialOutput(const SurfaceInterfaceRef& producer,
		const SurfaceFieldStamp& stamp) const
	{
		if (!HasTrialOutput() || !has_output_)
			throw std::runtime_error("FSI output is unavailable before a successful exact trial solve");
		if (!(producer == local_surface_.id))
			throw std::runtime_error("FSI output belongs to the wrong local endpoint");
		ValidateExactStamp(stamp, expected_output_, local_layout_);
	}

	void PrepareCommit()
	{
		if (phase_ != FsiTrialPhase::Solved || !has_output_)
			throw std::runtime_error("FSI prepare requires a solved exact trial");
		phase_ = FsiTrialPhase::Prepared;
	}

	void FinalizeCommit()
	{
		if (phase_ != FsiTrialPhase::Prepared || !has_output_)
			throw std::runtime_error("FSI finalize requires a prepared exact trial");
		committed_output_ = expected_output_;
		has_committed_output_ = true;
		ClearTrial();
		context_ = FsiTrialContext{};
		phase_ = FsiTrialPhase::Idle;
	}

	void RequireCommittedOutput(const SurfaceInterfaceRef& producer,
		const SurfaceFieldStamp& stamp) const
	{
		if (!has_committed_output_) throw std::runtime_error("FSI output has not been committed");
		if (!(producer == local_surface_.id))
			throw std::runtime_error("FSI committed output belongs to the wrong local endpoint");
		ValidateStoredExactStamp(stamp, committed_output_, local_layout_);
	}

	void RejectIteration()
	{
		if (!HasActiveStep()) throw std::runtime_error("FSI reject requires an active step");
		ClearTrial();
		phase_ = FsiTrialPhase::StepActive;
	}
	void RollbackTrial() { RejectIteration(); }
	void AbortStep() noexcept
	{
		ClearTrial();
		context_ = FsiTrialContext{};
		phase_ = FsiTrialPhase::Idle;
	}

	std::string ContextIdentitySha256() const
	{
		if (!HasActiveStep()) throw std::runtime_error("FSI lifecycle has no active context identity");
		ValidateFsiTrialContext(context_);
		Sha256 hash;
		distributed_surface_detail::AppendString(hash, "FsiTrialContext/v1");
		distributed_surface_detail::AppendString(hash, domain_id_);
		distributed_surface_detail::AppendString(hash, subsystem_id_);
		distributed_surface_detail::AppendString(hash, BuildFsiCouplingEdgeIdentitySha256(edge_));
		distributed_surface_detail::AppendString(hash, local_surface_.id.domain_id);
		distributed_surface_detail::AppendString(hash, local_surface_.id.subsystem_id);
		distributed_surface_detail::AppendString(hash, local_surface_.id.interface_id);
		distributed_surface_detail::AppendString(hash, peer_surface_.id.domain_id);
		distributed_surface_detail::AppendString(hash, peer_surface_.id.subsystem_id);
		distributed_surface_detail::AppendString(hash, peer_surface_.id.interface_id);
		distributed_surface_detail::AppendString(hash, local_layout_.layout_identity_sha256);
		distributed_surface_detail::AppendString(hash, BuildDistributedSurfacePartitionIdentitySha256(local_layout_));
		distributed_surface_detail::AppendString(hash, peer_layout_.layout_identity_sha256);
		distributed_surface_detail::AppendString(hash, BuildDistributedSurfacePartitionIdentitySha256(peer_layout_));
		hash.AppendLittleEndian64(context_.step);
		hash.AppendNormalizedDouble(context_.start_time_s);
		hash.AppendNormalizedDouble(context_.dt_s);
		hash.AppendLittleEndian64(context_.coupling_iteration);
		return hash.Hex();
	}

private:
	void ValidateExpectedStamp(const SurfaceFieldStamp& stamp,
		const DistributedSurfaceLayout& layout) const
	{
		ValidateSurfaceFieldStamp(stamp, layout);
		if (stamp.step != context_.step || stamp.coupling_iteration != context_.coupling_iteration
			|| stamp.time_s != context_.EndTime())
			throw std::runtime_error("FSI field stamp does not match the active macro step, time, or iteration");
	}
	void ValidateExactStamp(const SurfaceFieldStamp& value, const SurfaceFieldStamp& expected,
		const DistributedSurfaceLayout& layout) const
	{
		ValidateExpectedStamp(value, layout);
		ValidateExpectedStamp(expected, layout);
		if (BuildSurfaceFieldStampIdentitySha256(value, layout)
			!= BuildSurfaceFieldStampIdentitySha256(expected, layout))
			throw std::runtime_error("FSI field stamp is stale or belongs to a different trial");
	}
	void ValidateStoredExactStamp(const SurfaceFieldStamp& value, const SurfaceFieldStamp& expected,
		const DistributedSurfaceLayout& layout) const
	{
		ValidateSurfaceFieldStamp(value, layout);
		ValidateSurfaceFieldStamp(expected, layout);
		if (BuildSurfaceFieldStampIdentitySha256(value, layout)
			!= BuildSurfaceFieldStampIdentitySha256(expected, layout))
			throw std::runtime_error("FSI committed field stamp is stale or belongs to a different output");
	}
	void ClearTrial() noexcept
	{
		has_expected_stamps_ = false;
		has_input_ = false;
		has_output_ = false;
		expected_input_ = SurfaceFieldStamp{};
		expected_output_ = SurfaceFieldStamp{};
	}

	std::string domain_id_;
	std::string subsystem_id_;
	FsiCouplingEdge edge_;
	DistributedSurfaceInterface local_surface_;
	DistributedSurfaceInterface peer_surface_;
	DistributedSurfaceLayout local_layout_;
	DistributedSurfaceLayout peer_layout_;
	FsiTrialContext context_;
	SurfaceFieldStamp expected_input_;
	SurfaceFieldStamp expected_output_;
	SurfaceFieldStamp committed_output_;
	FsiTrialPhase phase_ = FsiTrialPhase::Idle;
	bool has_expected_stamps_ = false;
	bool has_input_ = false;
	bool has_output_ = false;
	bool has_committed_output_ = false;
};

class FsiFluidDomainRuntime {
public:
	virtual ~FsiFluidDomainRuntime() = default;
	virtual const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept = 0;
	virtual void SetSurfaceKinematics(const std::string& interface_id,
		const SurfaceKinematics& kinematics) = 0;
	virtual SurfaceTraction GetSurfaceTraction(const std::string& interface_id) const = 0;
};

class FsiStructureDomainRuntime {
public:
	virtual ~FsiStructureDomainRuntime() = default;
	virtual const std::vector<DistributedSurfaceInterface>& SurfaceInterfaces() const noexcept = 0;
	virtual void SetSurfaceTraction(const std::string& interface_id,
		const SurfaceTraction& traction) = 0;
	virtual SurfaceKinematics GetSurfaceKinematics(const std::string& interface_id) const = 0;
};

} // namespace iga

#endif

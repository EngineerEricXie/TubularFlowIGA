#ifndef IGA_THREE_D_IMMERSED_MOVING_DISTRIBUTED_FLOW_DOMAIN_HPP
#define IGA_THREE_D_IMMERSED_MOVING_DISTRIBUTED_FLOW_DOMAIN_HPP

#include "ThreeDImmersedDistributedDomain.hpp"
#include "../solvers/cpu/include/ImmersedMovingTransientDistributedRuntime.hpp"
#include <functional>

namespace iga {

// Explicit case policy: the reference flow in the numerical options must have
// the physical scale of the case. No acceptance tolerance is inferred here.
struct ImmersedMovingConservationLimits {
	double divergence_theorem,reynolds,moving_mass,wall_relative_leakage,discrete_continuity;
};

inline void ValidateImmersedMovingGraphInterval(const MaterialSurfaceKinematics& geometry,const DomainStepContext& step)
{
	if(geometry.StepStartS()!=step.start_time_s || geometry.StepEndS()!=step.EndTime()
		|| geometry.EvaluatedTimeS()!=step.EndTime())
		throw std::invalid_argument("moving graph geometry does not match its step");
	// DtS is the difference of these validated endpoints. Requiring that
	// subtraction to equal the nominal dt rejects ordinary decimal steps.
}

// Borrowed runtime and immutable geometry provider outlive the graph domain.
// The provider performs LOCAL geometry construction only; the backend wraps it
// in collective failure coordination. Numerical work is never run by it.
class ImmersedMovingDistributedGraphBackend
{
public:
	using GeometryProvider=std::function<std::unique_ptr<MovingCutGeometry>(
		const MovingCutGeometry&,const DomainStepContext&)>;
	ImmersedMovingDistributedGraphBackend(ImmersedMovingTransientDistributedRuntime& runtime,
		const GeometryProvider& provider,std::uint32_t layers,ImmersedMovingConservationLimits limits)
		: runtime_(runtime),layers_(layers),limits_(limits)
	{
		std::string policy;
		CollectiveLocalStage(Communicator(),"moving graph geometry and policy",[&] {
			if(!provider)throw std::invalid_argument("moving graph geometry provider is absent");
			provider_=provider;
			std::ostringstream text;text.exceptions(std::ios::badbit|std::ios::failbit);
			text<<std::setprecision(std::numeric_limits<double>::max_digits10)<<layers_;
			for(double value:{limits_.divergence_theorem,limits_.reynolds,limits_.moving_mass,
				limits_.wall_relative_leakage,limits_.discrete_continuity}) {
				if(!std::isfinite(value)||value<0)throw std::invalid_argument("moving conservation limits must be finite and nonnegative");
				text<<':'<<value;
			}
			policy=text.str();
		});
		RequireCollectiveSameText(Communicator(),"moving graph policy agreement",policy);
	}
	MPI_Comm Communicator() const noexcept
	{
		return runtime_.Communicator();
	}
	const ImmersedMovingDistributedClock& Clock() const noexcept
	{
		return runtime_.Clock();
	}
	const std::vector<ImmersedFlowPortDefinition>& PortDefinitions() const noexcept
	{
		return runtime_.PortDefinitions();
	}
	const ImmersedStaticFlowDiagnostics& Diagnostics() const noexcept
	{
		return runtime_.Diagnostics();
	}
	// Graph snapshots always refer to accepted rows, even while the target
	// active layout owns a different number of rows on this rank.
	PetscInt RowBegin() const noexcept
	{
		return runtime_.CommittedRowBegin();
	}
	PetscInt RowEnd() const noexcept
	{
		return runtime_.CommittedRowEnd();
	}
	Vec CommittedState() const
	{
		return runtime_.CommittedState();
	}
	void SetPortControlValue(const std::string& id,double value)
	{
		runtime_.SetPortControlValue(id,value);
	}
	void BeginTrial(const DomainStepContext& step)
	{
		std::unique_ptr<MovingCutGeometry> geometry;
		CollectiveLocalStage(Communicator(),"moving graph target geometry",[&] {
			geometry=provider_(runtime_.CommittedGeometry(),step);
			if(!geometry)throw std::invalid_argument("moving graph geometry is absent");
			ValidateImmersedMovingGraphInterval(geometry->Evaluation(),step);
		});
		runtime_.BeginTrial(std::move(geometry),static_cast<std::uint64_t>(step.step_index)+1,layers_);
	}
	bool SolveTrial()
	{
		if(!runtime_.SolveTrial())return false;
		CheckConservation();return true;
	}
	void PrepareCommit()
	{
		CheckConservation();runtime_.PrepareCommit();
	}
	void FinalizeCommit() noexcept
	{
		runtime_.FinalizeCommit();
	}
	void AbortPrepared() noexcept
	{
		runtime_.AbortPrepared();
	}
	void AbortTrial()
	{
		runtime_.AbortTrial();
	}
private:
	void CheckConservation() const
	{
		const auto value=runtime_.ConservationDiagnostics();
		CollectiveLocalStage(Communicator(),"moving graph physical conservation",[&] {
			const std::array<double,5> actual{{value.normalized_divergence_theorem_defect,value.normalized_reynolds_defect,
				value.normalized_moving_mass_defect,value.normalized_wall_relative_leakage,value.endpoint.normalized_discrete_moving_wall_continuity_defect}};
			const std::array<double,5> limit{{limits_.divergence_theorem,limits_.reynolds,limits_.moving_mass,limits_.wall_relative_leakage,limits_.discrete_continuity}};
			for(std::size_t i=0;i<actual.size();++i)if(!std::isfinite(actual[i])||actual[i]>limit[i]) {
				std::ostringstream message;message<<std::setprecision(17)<<"moving graph conservation gate "<<i<<" failed: "<<actual[i]<<" > "<<limit[i];
				throw std::runtime_error(message.str());
			}
		});
	}
	ImmersedMovingTransientDistributedRuntime& runtime_;
	GeometryProvider provider_;
	std::uint32_t layers_;
	ImmersedMovingConservationLimits limits_;
};

struct ImmersedMovingDomainLifecycle {
	static const char* Name() noexcept
	{
		return "moving_backward_euler";
	}
	static void ValidateInitial(const ImmersedMovingDistributedGraphBackend& runtime,double time)
	{
		if(runtime.Clock().time_s!=time||runtime.Clock().index!=0||runtime.Clock().trial_active)
			throw std::invalid_argument("moving graph requires a fresh matching backend clock");
	}
	static void ValidateStep(const ImmersedMovingDistributedGraphBackend& runtime,const DomainStepContext& step)
	{
		if(runtime.Clock().trial_active||runtime.Clock().time_s!=step.start_time_s
			||runtime.Clock().index!=static_cast<std::uint64_t>(step.step_index))
			throw std::invalid_argument("moving graph and backend clocks differ");
	}
	static void Begin(ImmersedMovingDistributedGraphBackend& runtime,const DomainStepContext& step)
	{
		runtime.BeginTrial(step);
	}
	static void Abort(ImmersedMovingDistributedGraphBackend& runtime)
	{
		if(runtime.Clock().trial_active)runtime.AbortTrial();
	}
};

class ThreeDImmersedMovingDistributedFlowDomain final
	: public ThreeDImmersedDistributedDomain<ImmersedMovingDistributedGraphBackend,ImmersedMovingDomainLifecycle> {
public:
	using ThreeDImmersedDistributedDomain::ThreeDImmersedDistributedDomain;
};

} // namespace iga
#endif

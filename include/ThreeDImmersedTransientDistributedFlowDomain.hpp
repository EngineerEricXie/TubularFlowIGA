#ifndef IGA_THREE_D_IMMERSED_TRANSIENT_DISTRIBUTED_FLOW_DOMAIN_HPP
#define IGA_THREE_D_IMMERSED_TRANSIENT_DISTRIBUTED_FLOW_DOMAIN_HPP

#include "ThreeDImmersedDistributedDomain.hpp"
#include "../solvers/cpu/include/ImmersedTransientDistributedRuntime.hpp"

namespace iga {

struct ImmersedTransientDomainLifecycle {
	static const char* Name() noexcept { return "backward_euler"; }
	static void ValidateInitial(const ImmersedTransientDistributedRuntime& runtime,double time)
	{
		if (runtime.Clock().time_s != time || runtime.Clock().index != 0 || runtime.Clock().trial_active)
			throw std::invalid_argument("transient graph adapter requires a fresh matching backend clock");
	}
	static void ValidateStep(const ImmersedTransientDistributedRuntime& runtime,const DomainStepContext& step)
	{
		if (runtime.Clock().trial_active || runtime.Clock().time_s != step.start_time_s
			|| runtime.Clock().index != static_cast<std::uint64_t>(step.step_index))
			throw std::invalid_argument("transient graph and backend committed clocks differ");
	}
	static void Begin(ImmersedTransientDistributedRuntime& runtime,const DomainStepContext& step)
	{
		// All controlled boundaries have been applied. Each coupling iteration
		// freezes the same accepted velocity, never the preceding trial field.
		runtime.BeginTrial(step.EndTime(),static_cast<std::uint64_t>(step.step_index)+1,step.dt_s);
	}
	static void Abort(ImmersedTransientDistributedRuntime& runtime)
	{
		// Newton's Rollback deliberately retains frozen inputs. Coupling retry
		// must discard them so new controls and a new force snapshot can be set.
		if (runtime.Clock().trial_active) runtime.AbortTrial();
	}
};

class ThreeDImmersedTransientDistributedFlowDomain final
	: public ThreeDImmersedDistributedDomain<ImmersedTransientDistributedRuntime,ImmersedTransientDomainLifecycle> {
public:
	using ThreeDImmersedDistributedDomain::ThreeDImmersedDistributedDomain;
};

} // namespace iga
#endif

#ifndef IGA_THREE_D_IMMERSED_DISTRIBUTED_FLOW_DOMAIN_HPP
#define IGA_THREE_D_IMMERSED_DISTRIBUTED_FLOW_DOMAIN_HPP

#include "ThreeDImmersedDistributedDomain.hpp"
#include "../solvers/cpu/include/ImmersedStaticDistributedRuntime.hpp"

namespace iga {

struct ImmersedSteadyDomainLifecycle {
	static const char* Name() noexcept { return "steady"; }
	static void ValidateInitial(const ImmersedStaticDistributedRuntime&,double) noexcept {}
	static void ValidateStep(const ImmersedStaticDistributedRuntime&,const DomainStepContext&) noexcept {}
	static void Begin(ImmersedStaticDistributedRuntime&,const DomainStepContext&) noexcept {}
	static void Abort(ImmersedStaticDistributedRuntime& runtime)
	{
		if (runtime.Diagnostics().trial_active) runtime.Rollback();
	}
};

// Coupling time leaves the steady operator's dt=0 contract unchanged.
class ThreeDImmersedDistributedFlowDomain final
	: public ThreeDImmersedDistributedDomain<ImmersedStaticDistributedRuntime,ImmersedSteadyDomainLifecycle> {
public:
	using ThreeDImmersedDistributedDomain::ThreeDImmersedDistributedDomain;
};

} // namespace iga
#endif

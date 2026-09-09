#ifndef IGA_SEQUENTIAL_EXECUTION_HPP
#define IGA_SEQUENTIAL_EXECUTION_HPP

#include "RuntimeConstruction.hpp"
#ifdef IGA_SEQUENTIAL_EXECUTION_TESTING
#include <functional>
#endif

namespace iga {

// Exactly one runtime operation, whose internal collectives must coordinate
// their own errors. Only its returned outcome is agreed here. Do not put
// allocating preparation or multiple collective operations in the callback.
template<class Work>
void SequentialRuntimeStage(MPI_Comm communicator, const char* stage, Work&& work)
{
	std::exception_ptr error;
	try { std::forward<Work>(work)(); }
	catch (...) { error = std::current_exception(); }
	RuntimeConstructionStage(communicator, stage, [&] {
		if (error) std::rethrow_exception(error);
	});
}

class OneDFlowRuntime;
class TransientFlowRuntime;

#ifdef IGA_SEQUENTIAL_EXECUTION_TESTING
struct SequentialExecutionHooks {
	// Both hooks execute inside a coordinated LOCAL stage. Observers must not
	// enter MPI; local vector reads may be used to check committed snapshots.
	std::function<void(int, int, bool&)> convergence;
	std::function<void(const char*, int, const OneDFlowRuntime&,
		const TransientFlowRuntime&, const OneDFlowRuntime&)> observe;
};

inline SequentialExecutionHooks& SequentialExecutionHooksForTesting()
{
	static thread_local SequentialExecutionHooks hooks;
	return hooks;
}
#endif

inline bool SequentialStrongConverged(MPI_Comm communicator, int step, int iteration, bool converged)
{
	RuntimeConstructionStage(communicator, "strong convergence input", [&] {
#ifdef IGA_SEQUENTIAL_EXECUTION_TESTING
		const auto& hook = SequentialExecutionHooksForTesting().convergence;
		if (hook) hook(step, iteration, converged);
#else
		(void)step; (void)iteration;
#endif
	});
	const int local = converged ? 1 : 0;
	int global = 0;
	{
		PhaseScope communication_phase(ProfilePhase::Communication);
		MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
	}
	return global != 0;
}

inline void ObserveSequentialStrong(const char* stage, int step,
	const OneDFlowRuntime& upstream, const TransientFlowRuntime& flow, const OneDFlowRuntime& downstream)
{
#ifdef IGA_SEQUENTIAL_EXECUTION_TESTING
	const auto& hook = SequentialExecutionHooksForTesting().observe;
	if (hook) hook(stage, step, upstream, flow, downstream);
#else
	(void)stage; (void)step; (void)upstream; (void)flow; (void)downstream;
#endif
}

} // namespace iga

#endif

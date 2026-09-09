#ifndef IGA_COLLECTIVE_PRESSURE_FLOW_EXECUTION_HPP
#define IGA_COLLECTIVE_PRESSURE_FLOW_EXECUTION_HPP

#include "CollectiveFailure.hpp"
#include "PressureFlowComponentExecutor.hpp"

namespace iga {

// Construct this policy inside a coordinated local configuration stage since
// std::function storage may allocate. The communicator is borrowed.
inline PressureFlowExecutionSynchronization CollectivePressureFlowExecution(MPI_Comm communicator)
{
	return {
		[communicator](const char* stage, std::exception_ptr error) {
			CollectiveLocalStage(communicator, stage, [&] {
				if (error) std::rethrow_exception(error);
			});
		},
		[communicator](bool converged) {
			const int local = converged ? 1 : 0;
			int global = 0;
			PhaseScope communication_phase(ProfilePhase::Communication);
			MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
			return global != 0;
		}
	};
}

} // namespace iga

#endif

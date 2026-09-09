#ifndef IGA_COLLECTIVE_SPECIES_PRESSURE_FLOW_EXECUTION_HPP
#define IGA_COLLECTIVE_SPECIES_PRESSURE_FLOW_EXECUTION_HPP

#include "CollectivePressureFlowExecution.hpp"
#include "SpeciesPressureFlowComponentExecutor.hpp"

namespace iga {

// Build within a coordinated local construction stage. The communicator is
// borrowed; graph/control agreement remains the caller's precondition.
inline SpeciesPressureFlowExecutionSynchronization CollectiveSpeciesPressureFlowExecution(MPI_Comm communicator)
{
	return {CollectivePressureFlowExecution(communicator),
		[communicator](const char* stage, std::string_view schedule) {
			RequireCollectiveSameText(communicator, stage, schedule);
		}};
}

} // namespace iga

#endif

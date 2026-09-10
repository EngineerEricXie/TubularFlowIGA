#ifndef IGA_DISTRIBUTED_FSI_COMMIT_COORDINATOR_HPP
#define IGA_DISTRIBUTED_FSI_COMMIT_COORDINATOR_HPP

#include "CollectiveFailure.hpp"
#include "FsiDomainRuntime.hpp"
#include "Sha256.hpp"

namespace iga {

// Runtime contract: PrepareCommit is collective on this communicator and
// either returns everywhere or throws everywhere. All Coordinator* methods
// are LOCAL, including abort and finalize. Every fallible result/publication
// operation must finish before this call; nothing may throw after the gate.
// This coordinates live-rank exceptions, not MPI process loss.
class DistributedFsiCommitCoordinator
{
public:
	template <class FluidRuntime, class StructureRuntime>
	static void Commit(MPI_Comm comm, const FsiTrialContext& context,
		FluidRuntime& fluid, StructureRuntime& structure)
	{
		static_assert(noexcept(fluid.CoordinatorAbortNoexcept()) && noexcept(structure.CoordinatorAbortNoexcept()),
			"paired FSI abort must be local and noexcept");
		static_assert(noexcept(fluid.CoordinatorFinalizeCommitNoexcept()) && noexcept(structure.CoordinatorFinalizeCommitNoexcept()),
			"paired FSI finalize must be local and noexcept");
		try {
			std::string identity;
			CollectiveLocalStage(comm,"paired FSI commit context",[&] {
				ValidateFsiTrialContext(context);
				fluid.CoordinatorRequireCommitContext(comm,context);
				structure.CoordinatorRequireCommitContext(comm,context);
				Sha256 hash;hash.AppendLittleEndian64(context.step);
				hash.AppendNormalizedDouble(context.start_time_s);hash.AppendNormalizedDouble(context.dt_s);
				hash.AppendLittleEndian64(context.coupling_iteration);identity=hash.Hex();
			});
			RequireCollectiveSameText(comm,"paired FSI commit epoch agreement",identity);
			fluid.PrepareCommit();
			structure.PrepareCommit();
			CollectiveLocalStage(comm,"paired FSI finalize gate",[&] {
				fluid.CoordinatorRequireFinalizeAllowed();
				structure.CoordinatorRequireFinalizeAllowed();
			});
		} catch (...) {
			fluid.CoordinatorAbortNoexcept();
			structure.CoordinatorAbortNoexcept();
			throw;
		}
		// No collective or fallible operation may separate these state changes.
		fluid.CoordinatorFinalizeCommitNoexcept();
		structure.CoordinatorFinalizeCommitNoexcept();
	}
};

} // namespace iga
#endif

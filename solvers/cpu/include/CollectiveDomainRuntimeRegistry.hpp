#ifndef IGA_COLLECTIVE_DOMAIN_RUNTIME_REGISTRY_HPP
#define IGA_COLLECTIVE_DOMAIN_RUNTIME_REGISTRY_HPP

#include "DomainRuntimeRegistry.hpp"
#include "RuntimeConstruction.hpp"

namespace iga {

// Both registry storage and its local index are prepared collectively before
// runtime ownership moves. If construction fails, the input owners remain with
// the caller on every rank and can be released in the caller's common order.
inline std::unique_ptr<DomainRuntimeRegistry> CreateCollectiveDomainRuntimeRegistry(
	MPI_Comm communicator, const SimulationGraph& graph,
	std::vector<std::unique_ptr<CoupledDomainRuntime>>& runtimes)
{
	DomainRuntimeRegistry::ConstructionOutcome synchronize;
	RuntimeConstructionStage(communicator, "registry synchronization ready", [&] {
		synchronize = [communicator](std::exception_ptr error) {
			RuntimeConstructionStage(communicator, "registry index ready", [&] {
				if (error) std::rethrow_exception(error);
			});
		};
	});
	return AllocateCollectiveRuntime<DomainRuntimeRegistry>(communicator,
		graph, std::move(runtimes), synchronize);
}

} // namespace iga

#endif

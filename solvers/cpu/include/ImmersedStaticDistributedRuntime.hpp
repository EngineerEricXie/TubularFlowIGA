#ifndef IGA_IMMERSED_STATIC_DISTRIBUTED_RUNTIME_HPP
#define IGA_IMMERSED_STATIC_DISTRIBUTED_RUNTIME_HPP

#include "ImmersedDistributedNewtonRuntime.hpp"

namespace iga {

// Preserve the steady constructor/API and its verified PETSc options prefix.
class ImmersedStaticDistributedRuntime : public ImmersedDistributedNewtonRuntime<ImmersedStaticDistributedOperator> {
public:
	ImmersedStaticDistributedRuntime(MPI_Comm communicator,
		const CartesianDomainClassification& domain,const CutCellVolumeQuadratureCatalog& volume,
		const ImmersedSurfaceQuadratureCatalog& surface,const CutCellGhostPenaltyCatalog& ghost,
		const ImmersedStaticFlowOptions& options,ImmersedWorkPartition partition = ImmersedWorkPartition::CellCount)
		: ImmersedDistributedNewtonRuntime(communicator,"immersed_static_",domain,volume,surface,ghost,options,partition)
	{}
};

} // namespace iga
#endif

#ifndef IGA_MULTIDOMAIN_RUNNER_HPP
#define IGA_MULTIDOMAIN_RUNNER_HPP

#include <petscsys.h>

namespace iga {

// Run one complete graph collectively on a borrowed communicator. PETSc must
// already be initialized. The caller keeps communicator valid throughout this
// synchronous call; all graph runtimes are destroyed before it returns.
// This does not initialize/finalize MPI and does not schedule separate domain
// groups inside the graph. Output paths must be distinct for independent runs.
int RunMultidomainFlow(int argc, char** argv, MPI_Comm communicator);

// Sequential positional/schema-v5 runner, with the same lifetime contract.
int RunSequentialFlow(int argc, char** argv, MPI_Comm communicator);

} // namespace iga

#endif

#ifndef IGA_RUNTIME_CLEANUP_HPP
#define IGA_RUNTIME_CLEANUP_HPP

#include "CollectiveFailure.hpp"

namespace iga {

#ifdef IGA_RUNTIME_CLEANUP_TESTING
using RuntimeCleanupProbe = PetscErrorCode (*)(const char*, PetscErrorCode) noexcept;
inline RuntimeCleanupProbe& RuntimeCleanupProbeForTesting()
{
	static thread_local RuntimeCleanupProbe probe = nullptr;
	return probe;
}
#endif

// Retain the first returned error without allocating or stopping the remaining
// releases. Only Check enters the agreement protocol, after all destroys return.
struct RuntimeCleanupResult
{
	PetscErrorCode error = 0;
	const char* operation = nullptr;

	void Observe(const char* name, PetscErrorCode code) noexcept
	{
#ifdef IGA_RUNTIME_CLEANUP_TESTING
		if (const auto probe = RuntimeCleanupProbeForTesting()) code = probe(name, code);
#endif
		if (!error && code) { error = code; operation = name; }
	}

	void Check(MPI_Comm communicator, const char* stage) const
	{
		CollectiveLocalStage(communicator, stage, [&] {
			if (error) throw std::runtime_error(std::string(operation)
				+" returned PETSc error "+std::to_string(error));
		});
	}
};

} // namespace iga
#endif

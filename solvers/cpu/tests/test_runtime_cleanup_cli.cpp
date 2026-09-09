#define IGA_RUNTIME_CLEANUP_TESTING
#define main NativeCleanupMain
#ifdef IGA_CLEANUP_GRAPH
#include "../../coupling/src/iga_1d_3d_bifurcation.cpp"
#elif defined(IGA_CLEANUP_SEQUENTIAL)
#include "../../coupling/src/iga_1d_3d_explicit.cpp"
#else
#include "../src/iga_navier_stokes.cpp"
#endif
#undef main

namespace {
const char* cleanup_target = nullptr;
bool selected_rank = false, injected = false;

PetscErrorCode CleanupFailure(const char* operation, PetscErrorCode code) noexcept
{
	if (selected_rank && !injected && std::strcmp(operation, cleanup_target) == 0) {
		injected = true;
		return PETSC_ERR_USER;
	}
	return code;
}
}

int main(int argc, char** argv)
{
	const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
	const char* target_rank = std::getenv("IGA_TEST_CLEANUP_RANK");
	cleanup_target = std::getenv("IGA_TEST_CLEANUP_TARGET");
	if (!rank || !target_rank || !cleanup_target) return 3;
	selected_rank = std::strcmp(rank, target_rank) == 0;
	iga::RuntimeCleanupProbeForTesting() = CleanupFailure;
	const int status = NativeCleanupMain(argc, argv);
	iga::RuntimeCleanupProbeForTesting() = nullptr;
	std::cerr << "cleanup_test rank=" << rank << " injected=" << injected
		<< " native_status=" << status << '\n';
	return status;
}

#include <hdf5.h>
#include <dlfcn.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

// Test-only Linux interposition at the C++ caller's H5Fcreate boundary.
// Throw before entering HDF5: no exception crosses a library C stack frame.
// Only the exact selected output path is affected, once per process.
extern "C" hid_t H5Fcreate(const char* name, unsigned flags, hid_t creation, hid_t access)
{
	using Create = hid_t (*)(const char*, unsigned, hid_t, hid_t);
	static const auto create = reinterpret_cast<Create>(dlsym(RTLD_NEXT, "H5Fcreate"));
	static bool injected = false;
	if (!create) return -1;
	const auto target = std::getenv("TUBULARFLOWIGA_TEST_HDF_CREATE_PATH");
	const auto mode = std::getenv("TUBULARFLOWIGA_TEST_HDF_CREATE_MODE");
	if (!injected && target && mode && std::strcmp(name, target) == 0) {
		injected = true;
		std::fprintf(stderr, "[hdf-create-injection] %s\n", mode);
		if (std::strcmp(mode, "bad-alloc") == 0) throw std::bad_alloc();
		if (std::strcmp(mode, "nonstandard") == 0) throw 7;
		if (std::strcmp(mode, "error") == 0) return -1;
	}
	return create(name, flags, creation, access);
}

#include <hdf5.h>
#include <dlfcn.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

// Linux test-only preload library. Fail once, only for the exact output file
// selected by the controller; retain the ID so destructor cleanup can retry.
extern "C" herr_t H5Fclose(hid_t id)
{
	using Close = herr_t (*)(hid_t);
	static const auto close = reinterpret_cast<Close>(dlsym(RTLD_NEXT, "H5Fclose"));
	static bool injected = false;
	if (!close) return -1;
	const auto target = std::getenv("TUBULARFLOWIGA_TEST_HDF_CLOSE");
	if (target && !injected) {
		const auto size = H5Fget_name(id, nullptr, 0);
		if (size >= 0) {
			std::vector<char> name(static_cast<std::size_t>(size)+1);
			if (H5Fget_name(id, name.data(), name.size()) == size
				&& std::strcmp(name.data(), target) == 0) {
				injected = true;
				std::fputs("[vtkhdf-close-injection] H5Fclose returned failure\n", stderr);
				return -1;
			}
		}
	}
	return close(id);
}

#include <dlfcn.h>
#include <unistd.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Test only: after returning the real selected file contents, turn its EOF
// into EIO. A streambuf copy must not silently accept the complete prefix.
extern "C" ssize_t read(int descriptor, void* buffer, size_t size)
{
	using Read = ssize_t (*)(int, void*, size_t);
	static const auto real_read = reinterpret_cast<Read>(dlsym(RTLD_NEXT, "read"));
	static bool injected = false;
	const char* target = std::getenv("IGA_TEST_READ_PATH");
	bool selected = false;
	if (target && !injected) {
		char link[64]{}, path[4096]{};
		std::snprintf(link, sizeof(link), "/proc/self/fd/%d", descriptor);
		const auto count = readlink(link, path, sizeof(path)-1);
		selected = count >= 0 && std::strcmp(path, target) == 0;
	}
	const auto count = real_read(descriptor, buffer, size);
	if (selected && count == 0) {
		injected = true;
		std::fputs("[text-read-injection] EOF replaced with EIO\n", stderr);
		errno = EIO;
		return -1;
	}
	return count;
}

#include <dlfcn.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Linux test-only library. Close the selected real file, then return an I/O
// failure once. Exact path matching leaves MPI and diagnostic files unaffected.
extern "C" int fclose(FILE* stream)
{
	using Close = int (*)(FILE*);
	static const auto close = reinterpret_cast<Close>(dlsym(RTLD_NEXT, "fclose"));
	static bool injected = false;
	if (!close) { errno = EIO; return EOF; }
	const char* target = std::getenv("TUBULARFLOWIGA_TEST_TEXT_CLOSE");
	bool selected = false;
	if (target && !injected) {
		char descriptor[64]{};
		char path[4096]{};
		std::snprintf(descriptor, sizeof(descriptor), "/proc/self/fd/%d", fileno(stream));
		const auto size = readlink(descriptor, path, sizeof(path)-1);
		selected = size >= 0 && std::strcmp(path, target) == 0;
	}
	const int status = close(stream);
	if (selected) {
		injected = true;
		std::fputs("[text-close-injection] fclose returned failure\n", stderr);
		errno = EIO;
		return EOF;
	}
	return status;
}

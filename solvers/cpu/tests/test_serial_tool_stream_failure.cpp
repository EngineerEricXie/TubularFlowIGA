// Invoke the native CLI with a fault in its actual numeric serializer.
#include "StringStreamFailure.hpp"
#include <cstdlib>
#define main SerialToolNativeMain
#if defined(IGA_TEST_WOMERSLEY_STREAM)
#include "../src/iga_womersley_reference.cpp"
#else
#include "../src/iga_transport_validate.cpp"
#endif
#undef main

int main(int argc, char** argv)
{
	const char* mode = std::getenv("IGA_TEST_FORMAT_MODE");
	const char* prefix = std::getenv("IGA_TEST_FORMAT_PREFIX");
	if (!mode || !prefix) return SerialToolNativeMain(argc, argv);
	bool injected = false;
	int status = 0;
	{
		iga_test::ScopedStringStreamFailure fault(true, std::atoi(mode), prefix);
		try { status = SerialToolNativeMain(argc, argv); }
		catch (...) { status = 1; }
		injected = fault.Injected();
	}
	std::cerr << "serial_formatter injected=" << injected << " native_status=" << status << '\n';
	return injected ? status : 3;
}

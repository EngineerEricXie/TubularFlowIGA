#define main FsiExporterNativeMain
#include "../src/phase8_compliant_channel_fsi_paraview.cpp"
#undef main
#include "StringStreamFailure.hpp"

int main(int argc, char** argv)
{
	try {
		Require(argc == 2, "provide exported text path");
		const auto text = ReadAll(argv[1]);
		Require(!text.empty(), "empty exported text");
		for (int mode = 0; mode < 3; ++mode) {
			bool rejected = false, injected = false;
			{
				iga_test::ScopedStringStreamFailure fault(true, mode);
				try { Values([](auto& stream) { stream << 1.25; }); }
				catch (...) { rejected = true; }
				injected = fault.Injected();
			}
			Require(rejected && injected, "export Values accepted formatter failure");
			Require(Values([](auto& stream) { stream << 1.25; }) == "1.25", "retry changed values");
		}
		std::cout << "fsi export helper faults=3 retries=3 passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

#include "../../cpu/tests/StringStreamFailure.hpp"

#define main OneDNativeMain
#include "../src/iga_1d.cpp"
#undef main

// The real CLI owns MPI initialization/finalization and its output stages.
// Only the profile filename's numeric formatting is instrumented by a facet.
int main(int argc, char** argv)
{
	if (argc < 4) return 2;
	const int mode = std::stoi(argv[1]), step = std::stoi(argv[2]);
	std::vector<char*> arguments{argv[0]};
	arguments.insert(arguments.end(), argv+3, argv+argc);
	arguments.push_back(nullptr);
	iga_test::ScopedStringStreamFailure fault(mode >= 0, mode, "profile_1d_", step);
	const int result = OneDNativeMain(static_cast<int>(arguments.size())-1, arguments.data());
	if (fault.Injected()) std::cerr << "filename formatting injected mode=" << mode << " step=" << step << '\n';
	return result;
}

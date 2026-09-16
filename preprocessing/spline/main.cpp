#include <iostream>
#include "kernel.h"
#include <exception>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <string>

#include <omp.h>

using namespace std;

namespace {

void PrintUsage()
{
	cerr << "usage: spline CASE_DIR/ [--no-legacy-text] [--legacy-vtk] [--threads N]\n";
}

int ParseThreadCount(const char* text)
{
	errno = 0;
	char* end = nullptr;
	const long value = strtol(text, &end, 10);
	if (errno != 0 || end == text || *end != '\0' || value <= 0 || value > INT_MAX)
	{
		throw invalid_argument("--threads requires a positive integer");
	}
	return static_cast<int>(value);
}

}

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		PrintUsage();
		return 2;
	}
	try
	{
		bool legacy_text = true;
		bool legacy_vtk = false;
		for (int argument = 2; argument < argc; ++argument)
		{
			const string option(argv[argument]);
			if (option == "--no-legacy-text") legacy_text = false;
			else if (option == "--legacy-vtk") legacy_vtk = true;
			else if (option == "--threads")
			{
				if (++argument == argc)
				{
					PrintUsage();
					return 2;
				}
				try
				{
					omp_set_num_threads(ParseThreadCount(argv[argument]));
				}
				catch (const invalid_argument& error)
				{
					cerr << "spline: " << error.what() << "\n";
					PrintUsage();
					return 2;
				}
			}
			else
			{
				PrintUsage();
				return 2;
			}
		}
		kernel app;
		app.run(argv[1], legacy_text, legacy_vtk);
		cout << "DONE!\n";
		return 0;
	}
	catch (const exception& error)
	{
		cerr << "spline: " << error.what() << "\n";
		return 1;
	}
}

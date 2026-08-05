#include <iostream>
#include "kernel.h"
#include <exception>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2 || argc > 3)
	{
		cerr << "usage: spline CASE_DIR/ [--no-legacy-text]\n";
		return 2;
	}
	try
	{
		const bool legacy_text = argc == 2;
		if (!legacy_text && string(argv[2]) != "--no-legacy-text")
		{
			cerr << "usage: spline CASE_DIR/ [--no-legacy-text]\n";
			return 2;
		}
		kernel app;
		app.run(argv[1], legacy_text);
		cout << "DONE!\n";
		return 0;
	}
	catch (const exception& error)
	{
		cerr << "spline: " << error.what() << "\n";
		return 1;
	}
}

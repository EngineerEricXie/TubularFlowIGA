#include <iostream>
#include "kernel.h"
#include <exception>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "usage: spline CASE_DIR/ [--no-legacy-text] [--legacy-vtk]\n";
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
			else
			{
				cerr << "usage: spline CASE_DIR/ [--no-legacy-text] [--legacy-vtk]\n";
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

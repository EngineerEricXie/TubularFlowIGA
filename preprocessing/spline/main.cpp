#include <iostream>
#include "kernel.h"
#include <exception>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		cerr << "usage: spline CASE_DIR/\n";
		return 2;
	}
	try
	{
		kernel app;
		app.run(argv[1]);
		cout << "DONE!\n";
		return 0;
	}
	catch (const exception& error)
	{
		cerr << "spline: " << error.what() << "\n";
		return 1;
	}
}

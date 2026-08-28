#include "VtkOutput.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

int main()
{
	assert(iga::ResolveVisualizationFormat(
		iga::VisualizationFormat::Automatic, true)
		== iga::VisualizationFormat::BezierVtkHdf);
	assert(iga::ResolveVisualizationFormat(
		iga::VisualizationFormat::Automatic, false)
		== iga::VisualizationFormat::Vtu);
	assert(iga::ResolveVisualizationFormat(
		iga::VisualizationFormat::Vtu, true)
		== iga::VisualizationFormat::Vtu);
	const auto directory = fs::temp_directory_path()/"tubularflowiga-vtk-output-test";
	fs::create_directories(directory);
	const auto mesh = directory/"mesh.vtk";
	{
		std::ofstream output(mesh);
		output << "# vtk DataFile Version 3.0\ntest\nASCII\n"
			<< "DATASET UNSTRUCTURED_GRID\n"
			<< "POINTS 4 double\n0 0 0 1 0 0 0 1 0 0 0 1\n"
			<< "CELLS 1 5\n4 0 1 2 3\n"
			<< "CELL_TYPES 1\n10\n";
	}
	const auto grid = iga::ReadLegacyUnstructuredGrid(mesh);
	assert(grid.points.size() == 12);
	assert(grid.offsets.size() == 1 && grid.offsets[0] == 4);
	const auto vtu = directory/"field.step000001.vtu";
	iga::WriteVtu(mesh, vtu,
		{{"scalar", 1, {1.0, 2.0, 3.0, 4.0}},
		 {"vector", 3, {1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
			0.0, 0.0, 1.0, 1.0, 1.0, 1.0}}}, 0.25);
	std::ifstream input(vtu);
	std::ostringstream contents;
	contents << input.rdbuf();
	assert(contents.str().find("Name=\"scalar\"") != std::string::npos);
	assert(contents.str().find("Name=\"TimeValue\"") != std::string::npos);
	const auto pvd = directory/"field.pvd";
	iga::WritePvd(pvd, {{0.25, vtu}});
	std::ifstream pvd_input(pvd);
	std::ostringstream pvd_contents;
	pvd_contents << pvd_input.rdbuf();
	assert(pvd_contents.str().find("field.step000001.vtu") != std::string::npos);
	fs::remove_all(directory);
}

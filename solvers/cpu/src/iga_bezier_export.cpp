#include "BezierVisualization.hpp"
#include "IgaDatabase.hpp"
#include "TemporalVtkHdf.hpp"

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	if (argc < 2 || argc > 3) {
		std::cerr << "usage: iga_bezier_export DATABASE.ntiga [OUTPUT.vtkhdf]\n";
		return 2;
	}
	try {
		const fs::path database_path(argv[1]);
		const fs::path output = argc == 3
			? fs::path(argv[2])
			: database_path.parent_path()/(database_path.stem().string()+".vtkhdf");
		iga::Database database(database_path.string());
		const auto mesh = iga::BuildSourceCoordinateBezierVisualizationMesh(
			database, false);
		const auto report = iga::BezierGeometryReportPath(output);
		iga::WriteBezierGeometryReport(report, mesh.validation);
		iga::RequireValidBezierGeometry(mesh.validation);
		iga::TemporalVtkHdfWriter writer(output, mesh);
		writer.Append(0.0, {});
		std::cout << "wrote cubic Bezier visualization " << output
			<< " with " << mesh.types.size() << " elements and "
			<< mesh.points.size() << " unique points\n"
			<< "geometry report " << report << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "iga_bezier_export: " << error.what() << '\n';
		return 1;
	}
}

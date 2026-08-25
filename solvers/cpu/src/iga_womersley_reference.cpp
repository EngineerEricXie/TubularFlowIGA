#include "IgaDatabase.hpp"
#include "WomersleyReference.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	try {
		if (argc != 5 && argc != 6)
			throw std::runtime_error(
				"usage: iga_womersley_reference DATABASE.ntiga CONTROL_MESH.vtk "
				"WOMERSLEY.json OUTPUT_DIR [MANIFEST.csv]");
		iga::Database database(argv[1]);
		auto points = iga::ReadVtkPointCoordinates(argv[2], database.header().nodes);
		iga::RescalePointCoordinatesLikeSpline(points);
		const auto configuration = iga::ReadWomersleyReferenceConfiguration(argv[3]);
		const fs::path output_directory(argv[4]);
		const fs::path manifest_name(argc == 6 ? argv[5] : "velocity_series.csv");
		if (manifest_name.is_absolute() || manifest_name.has_parent_path())
			throw std::runtime_error("MANIFEST.csv must be a simple relative file name");
		fs::create_directories(output_directory);
		std::ofstream manifest(output_directory/manifest_name);
		if (!manifest) throw std::runtime_error("cannot create Womersley manifest");
		manifest << std::setprecision(17) << "time,file\n";
		for (std::size_t index = 0; index < configuration.sample_times.size(); ++index) {
			std::ostringstream name;
			name << configuration.file_prefix << ".step"
				<< std::setw(6) << std::setfill('0') << index+1 << ".txt";
			std::ofstream output(output_directory/name.str());
			if (!output) throw std::runtime_error(
				"cannot create Womersley velocity field: "+name.str());
			output << std::setprecision(17);
			for (const auto& point : points) {
				const auto velocity = iga::WomersleyVelocity(
					configuration, point, configuration.sample_times[index]);
				output << velocity[0] << ' ' << velocity[1] << ' ' << velocity[2] << '\n';
			}
			if (!output) throw std::runtime_error(
				"cannot write Womersley velocity field: "+name.str());
			manifest << configuration.sample_times[index] << ',' << name.str() << '\n';
		}
		if (!manifest) throw std::runtime_error("cannot write Womersley manifest");
		std::cout << "womersley_nodes=" << points.size()
			<< " snapshots=" << configuration.sample_times.size()
			<< " manifest=" << (output_directory/manifest_name).string() << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "iga_womersley_reference: " << error.what() << '\n';
		return 1;
	}
}

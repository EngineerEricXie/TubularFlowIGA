#ifndef IGA_VTK_OUTPUT_HPP
#define IGA_VTK_OUTPUT_HPP

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class VisualizationFormat { Automatic, Vtu, BezierVtkHdf };

inline VisualizationFormat ParseVisualizationFormat(const std::string& value)
{
	if (value == "auto") return VisualizationFormat::Automatic;
	if (value == "vtu") return VisualizationFormat::Vtu;
	if (value == "vtkhdf" || value == "bezier-vtkhdf")
		return VisualizationFormat::BezierVtkHdf;
	throw std::runtime_error(
		"--visualization-format must be auto, vtu, or vtkhdf");
}

inline VisualizationFormat ResolveVisualizationFormat(
	VisualizationFormat requested, bool transient)
{
	if (requested != VisualizationFormat::Automatic) return requested;
	return transient ? VisualizationFormat::BezierVtkHdf
		: VisualizationFormat::Vtu;
}

struct VtkPointArray {
	std::string name;
	int components = 1;
	std::vector<double> values;
};

struct LegacyUnstructuredGrid {
	std::vector<double> points;
	std::vector<std::int64_t> connectivity;
	std::vector<std::int64_t> offsets;
	std::vector<unsigned int> types;
};

inline std::string EscapeVtkXml(const std::string& value)
{
	std::string escaped;
	for (const char character : value) {
		if (character == '&') escaped += "&amp;";
		else if (character == '<') escaped += "&lt;";
		else if (character == '>') escaped += "&gt;";
		else if (character == '"') escaped += "&quot;";
		else escaped.push_back(character);
	}
	return escaped;
}

inline LegacyUnstructuredGrid ReadLegacyUnstructuredGrid(
	const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open VTK mesh: "+path.string());
	LegacyUnstructuredGrid grid;
	std::string token;
	std::size_t point_count = 0, cell_count = 0;
	while (input >> token) {
		if (token == "POINTS") {
			std::string type;
			input >> point_count >> type;
			grid.points.resize(3*point_count);
			for (auto& value : grid.points)
				if (!(input >> value)) throw std::runtime_error("truncated VTK points: "+path.string());
		} else if (token == "CELLS") {
			std::size_t entries = 0;
			input >> cell_count >> entries;
			grid.offsets.reserve(cell_count);
			for (std::size_t cell = 0; cell < cell_count; ++cell) {
				std::size_t nodes = 0;
				if (!(input >> nodes)) throw std::runtime_error("truncated VTK cells: "+path.string());
				for (std::size_t node = 0; node < nodes; ++node) {
					std::int64_t index = -1;
					if (!(input >> index) || index < 0
						|| static_cast<std::size_t>(index) >= point_count)
						throw std::runtime_error("invalid VTK cell connectivity: "+path.string());
					grid.connectivity.push_back(index);
				}
				grid.offsets.push_back(static_cast<std::int64_t>(grid.connectivity.size()));
			}
		} else if (token == "CELL_TYPES") {
			std::size_t count = 0;
			input >> count;
			grid.types.resize(count);
			for (auto& value : grid.types)
				if (!(input >> value)) throw std::runtime_error("truncated VTK cell types: "+path.string());
			break;
		}
	}
	if (point_count == 0 || cell_count == 0 || grid.offsets.size() != cell_count
		|| grid.types.size() != cell_count)
		throw std::runtime_error("VTK mesh is not a complete ASCII unstructured grid: "+path.string());
	return grid;
}

inline void WriteVtu(const std::filesystem::path& mesh_path,
	const std::filesystem::path& output_path,
	const std::vector<VtkPointArray>& arrays, double physical_time)
{
	const auto grid = ReadLegacyUnstructuredGrid(mesh_path);
	const auto points = grid.points.size()/3;
	for (const auto& array : arrays)
		if (array.components < 1
			|| array.values.size() != points*static_cast<std::size_t>(array.components))
			throw std::runtime_error("VTU point-array size mismatch for '"+array.name+"'");
	if (!output_path.parent_path().empty())
		std::filesystem::create_directories(output_path.parent_path());
	std::ofstream output(output_path);
	if (!output) throw std::runtime_error("cannot create VTU output: "+output_path.string());
	output << std::setprecision(17)
		<< "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <UnstructuredGrid>\n"
		<< "    <FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"
		<< physical_time << "</DataArray></FieldData>\n"
		<< "    <Piece NumberOfPoints=\"" << points << "\" NumberOfCells=\""
		<< grid.offsets.size() << "\">\n"
		<< "      <PointData>\n";
	for (const auto& array : arrays) {
		output << "        <DataArray type=\"Float64\" Name=\""
			<< EscapeVtkXml(array.name) << "\" NumberOfComponents=\""
			<< array.components << "\" format=\"ascii\">\n          ";
		for (const double value : array.values) output << value << ' ';
		output << "\n        </DataArray>\n";
	}
	output << "      </PointData>\n"
		<< "      <Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n        ";
	for (const double value : grid.points) output << value << ' ';
	output << "\n      </DataArray></Points>\n"
		<< "      <Cells>\n"
		<< "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
	for (const auto value : grid.connectivity) output << value << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
	for (const auto value : grid.offsets) output << value << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
	for (const auto value : grid.types) output << value << ' ';
	output << "\n        </DataArray>\n"
		<< "      </Cells>\n"
		<< "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
	if (!output) throw std::runtime_error("cannot write VTU output: "+output_path.string());
}

inline std::filesystem::path VtuStepPath(const std::filesystem::path& base, int step)
{
	std::ostringstream name;
	name << base.stem().string() << ".step" << std::setw(6) << std::setfill('0')
		<< step << ".vtu";
	return base.parent_path()/name.str();
}

inline std::filesystem::path VtuFinalPath(const std::filesystem::path& base)
{
	return base.parent_path()/(base.stem().string()+".vtu");
}

inline std::filesystem::path PvdPath(const std::filesystem::path& base)
{
	return base.parent_path()/(base.stem().string()+".pvd");
}

inline void WritePvd(const std::filesystem::path& path,
	const std::vector<std::pair<double, std::filesystem::path>>& snapshots)
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create PVD output: "+path.string());
	output << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <Collection>\n";
	for (const auto& snapshot : snapshots)
		output << "    <DataSet timestep=\"" << std::setprecision(17) << snapshot.first
			<< "\" group=\"\" part=\"0\" file=\""
			<< EscapeVtkXml(snapshot.second.filename().string()) << "\"/>\n";
	output << "  </Collection>\n</VTKFile>\n";
	if (!output) throw std::runtime_error("cannot write PVD output: "+path.string());
}

} // namespace iga

#endif

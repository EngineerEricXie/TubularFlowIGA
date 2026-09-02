#include "BezierVisualization.hpp"

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

iga::Element MakeElement(std::uint64_t id, double x_offset,
	std::int32_t node_offset, bool share_structured_nodes)
{
	iga::Element element;
	element.id = id;
	element.owner = 0;
	element.connectivity.resize(64);
	element.extraction.resize(64);
	std::size_t point = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++point) {
				element.bezier_points[point] = {{
					x_offset+i/3.0, j/3.0, k/3.0}};
				const auto node = share_structured_nodes
					? (static_cast<std::int32_t>(id)*3+i)+7*j+28*k
					: node_offset+static_cast<std::int32_t>(point);
				element.connectivity[point] = node;
				element.extraction[point].fill(0.0);
				element.extraction[point][point] = 1.0;
			}
	return element;
}

void WriteDatabase(const fs::path& path, const std::vector<iga::Element>& elements,
	std::uint64_t nodes, const iga::GeometryTransform& transform = {})
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	assert(output);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, std::uint32_t{1});
	iga::Write(output, static_cast<std::uint64_t>(elements.size()));
	iga::Write(output, nodes);
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto rank_index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	for (const double value : transform.source_origin) iga::Write(output, value);
	iga::Write(output, transform.source_units_per_normalized_unit);
	iga::Write(output, transform.source_length_scale_to_m);
	const auto offsets_position = output.tellp();
	for (std::size_t index = 0; index <= elements.size(); ++index)
		iga::Write(output, std::uint64_t{0});
	for (const auto& element : elements) iga::Write(output, element.owner);
	std::vector<std::uint64_t> offsets;
	for (const auto& element : elements) {
		offsets.push_back(static_cast<std::uint64_t>(output.tellp()));
		iga::Write(output, element.id);
		iga::Write(output, element.type);
		iga::Write(output, element.owner);
		iga::Write(output, static_cast<std::uint32_t>(element.connectivity.size()));
		output.write(reinterpret_cast<const char*>(element.boundary_labels.data()),
			static_cast<std::streamsize>(sizeof(element.boundary_labels)));
		output.write(reinterpret_cast<const char*>(element.connectivity.data()),
			static_cast<std::streamsize>(element.connectivity.size()*sizeof(std::int32_t)));
		for (const auto& row : element.extraction) {
			std::uint8_t nonzeros = 0;
			for (const auto coefficient : row) if (coefficient != 0.0) ++nonzeros;
			iga::Write(output, nonzeros);
			for (std::uint8_t column = 0; column < row.size(); ++column)
				if (row[column] != 0.0) {
					iga::Write(output, column);
					iga::Write(output, row[column]);
				}
		}
		output.write(reinterpret_cast<const char*>(element.bezier_points.data()),
			static_cast<std::streamsize>(sizeof(element.bezier_points)));
	}
	offsets.push_back(static_cast<std::uint64_t>(output.tellp()));
	const auto rank_index_offset = static_cast<std::uint64_t>(output.tellp());
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, static_cast<std::uint64_t>(elements.size()));
	for (std::uint64_t index = 0; index < elements.size(); ++index)
		iga::Write(output, index);
	output.seekp(offsets_position);
	output.write(reinterpret_cast<const char*>(offsets.data()),
		static_cast<std::streamsize>(offsets.size()*sizeof(std::uint64_t)));
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	assert(output);
}

} // namespace

int main()
{
	const auto directory = fs::temp_directory_path()/"tubularflowiga-bezier-visualization-test";
	fs::create_directories(directory);
	const auto conforming_path = directory/"conforming.ntiga";
	WriteDatabase(conforming_path,
		{MakeElement(0, 0.0, 0, true), MakeElement(1, 1.0, 0, true)}, 112,
		{{{10.0, 20.0, 30.0}}, 2.5, 0.001});
	iga::Database conforming_database(conforming_path.string());
	const auto mesh = iga::BuildBezierVisualizationMesh(conforming_database);
	assert(mesh.validation.local_points == 128);
	assert(mesh.validation.unique_points == 112);
	assert(mesh.validation.shared_point_references == 16);
	assert(mesh.validation.coincident_unmerged_points == 0);
	assert(mesh.validation.overlapping_element_pairs == 0);
	assert(mesh.points.size() == 112);
	assert(mesh.connectivity.size() == 128);
	assert(mesh.offsets == std::vector<std::int64_t>({0, 64, 128}));
	assert(mesh.connectivity[0] == 0);
	assert(mesh.connectivity[1] == 3);
	assert(mesh.connectivity[2] == 15);
	const auto source_mesh = iga::BuildSourceCoordinateBezierVisualizationMesh(
		conforming_database);
	assert(source_mesh.points.size() == mesh.points.size());
	const std::array<double, 3> source_origin{{10.0, 20.0, 30.0}};
	for (std::size_t point = 0; point < mesh.points.size(); ++point)
		for (int direction = 0; direction < 3; ++direction)
			assert(source_mesh.points[point][direction]
				== source_origin[direction]
					+2.5*mesh.points[point][direction]);
	assert(source_mesh.connectivity == mesh.connectivity);
	const auto expected_source_jacobian
		= mesh.validation.minimum_jacobian*2.5*2.5*2.5;
	assert(std::abs(source_mesh.validation.minimum_jacobian-expected_source_jacobian)
		<= 1.0e-14*std::abs(expected_source_jacobian));
	const auto report_path = directory/"geometry.json";
	iga::WriteBezierGeometryReport(report_path, mesh.validation);
	std::ifstream report_input(report_path);
	const std::string report((std::istreambuf_iterator<char>(report_input)),
		std::istreambuf_iterator<char>());
	assert(report.find("\"valid\": true") != std::string::npos);
	assert(report.find("\"unique_points\": 112") != std::string::npos);
	std::vector<double> control_values(112);
	for (std::size_t node = 0; node < control_values.size(); ++node)
		control_values[node] = static_cast<double>(node)+0.25;
	const auto arrays = iga::ExtractBezierPointArrays(
		mesh, {{"scalar", 1, control_values}});
	assert(arrays.size() == 1 && arrays[0].values.size() == 112);
	for (std::size_t point = 0; point < mesh.points.size(); ++point) {
		const auto entry = mesh.signature_offsets[point];
		assert(mesh.signature_offsets[point+1] == entry+1);
		assert(arrays[0].values[point]
			== control_values[static_cast<std::size_t>(mesh.signature_nodes[entry])]);
	}

	const auto coincident_path = directory/"coincident.ntiga";
	WriteDatabase(coincident_path,
		{MakeElement(0, 0.0, 0, false), MakeElement(1, 1.0, 64, false)}, 128);
	iga::Database coincident_database(coincident_path.string());
	const auto coincident = iga::BuildBezierVisualizationMesh(
		coincident_database, false);
	assert(coincident.validation.unique_points == 128);
	assert(coincident.validation.coincident_unmerged_points == 16);
	assert(coincident.validation.overlapping_element_pairs == 0);
	assert((coincident.validation.first_coincident_element_pair
		== std::array<std::uint64_t, 2>{{0, 1}}));

	const auto overlap_path = directory/"overlap.ntiga";
	WriteDatabase(overlap_path,
		{MakeElement(0, 0.0, 0, false), MakeElement(1, 0.5, 64, false)}, 128);
	iga::Database overlap_database(overlap_path.string());
	const auto overlap = iga::BuildBezierVisualizationMesh(overlap_database, false);
	assert(overlap.validation.coincident_unmerged_points == 0);
	assert(overlap.validation.overlap_candidates >= 1);
	assert(overlap.validation.overlapping_element_pairs >= 1);
	bool rejected = false;
	try {
		iga::RequireValidBezierGeometry(overlap.validation);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);

	fs::remove_all(directory);
}

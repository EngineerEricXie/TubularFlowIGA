#include "BezierVisualization.hpp"
#include "TemporalVtkHdf.hpp"

#include <hdf5.h>

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

namespace {

iga::BezierVisualizationMesh MakeMesh()
{
	iga::BezierVisualizationMesh mesh;
	std::size_t point = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++point) {
				mesh.points.push_back({{i/3.0, j/3.0, k/3.0}});
				mesh.signature_nodes.push_back(static_cast<std::int32_t>(point));
				mesh.signature_coefficients.push_back(1.0);
				mesh.signature_offsets.push_back(point+1);
			}
	for (const auto tensor : iga::detail::VtkCubicHexTensorIndices())
		mesh.connectivity.push_back(tensor);
	mesh.offsets = {0, 64};
	mesh.types = {iga::kVtkBezierHexahedron};
	mesh.higher_order_degrees = {{{3, 3, 3}}};
	mesh.element_ids = {17};
	mesh.element_owners = {0};
	mesh.validation.elements = 1;
	mesh.validation.local_points = 64;
	mesh.validation.unique_points = 64;
	return mesh;
}

std::vector<hsize_t> Dimensions(hid_t file, const char* path)
{
	const auto dataset = H5Dopen2(file, path, H5P_DEFAULT);
	assert(dataset >= 0);
	const auto space = H5Dget_space(dataset);
	assert(space >= 0);
	const auto rank = H5Sget_simple_extent_ndims(space);
	assert(rank > 0);
	std::vector<hsize_t> dimensions(static_cast<std::size_t>(rank));
	assert(H5Sget_simple_extent_dims(space, dimensions.data(), nullptr) >= 0);
	H5Sclose(space);
	H5Dclose(dataset);
	return dimensions;
}

template <class T>
std::vector<T> Read(hid_t file, const char* path, hid_t type)
{
	const auto dataset = H5Dopen2(file, path, H5P_DEFAULT);
	assert(dataset >= 0);
	const auto dimensions = Dimensions(file, path);
	std::size_t values = 1;
	for (const auto dimension : dimensions) values *= static_cast<std::size_t>(dimension);
	std::vector<T> result(values);
	assert(H5Dread(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, result.data()) >= 0);
	H5Dclose(dataset);
	return result;
}

} // namespace

int main()
{
	const auto directory = fs::temp_directory_path()/"tubularflowiga-temporal-vtkhdf-test";
	fs::create_directories(directory);
	const auto path = directory/"bezier.vtkhdf";
	const auto mesh = MakeMesh();
	std::vector<double> scalar(64);
	std::vector<double> velocity(64*3);
	for (std::size_t point = 0; point < 64; ++point) {
		scalar[point] = static_cast<double>(point);
		velocity[3*point] = static_cast<double>(point);
		velocity[3*point+1] = 2.0*point;
		velocity[3*point+2] = -static_cast<double>(point);
	}
	{
		iga::TemporalVtkHdfWriter writer(path, mesh);
		writer.Append(0.0, {{"scalar", 1, scalar}, {"velocity", 3, velocity}});
		for (auto& value : scalar) value += 10.0;
		writer.Append(0.5, {{"scalar", 1, scalar}, {"velocity", 3, velocity}});
	}
	{
		iga::TemporalVtkHdfWriter writer(path, mesh, true);
		for (auto& value : scalar) value += 10.0;
		writer.Append(1.0, {{"scalar", 1, scalar}, {"velocity", 3, velocity}});
	}
	const auto file = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	assert(file >= 0);
	const auto group = H5Gopen2(file, "/VTKHDF", H5P_DEFAULT);
	assert(group >= 0);
	const auto type_attribute = H5Aopen(group, "Type", H5P_DEFAULT);
	assert(type_attribute >= 0);
	const auto type = H5Aget_type(type_attribute);
	std::vector<char> type_value(H5Tget_size(type)+1, '\0');
	assert(H5Aread(type_attribute, type, type_value.data()) >= 0);
	assert(std::string(type_value.data()) == "UnstructuredGrid");
	H5Tclose(type);
	H5Aclose(type_attribute);
	const auto version_attribute = H5Aopen(group, "Version", H5P_DEFAULT);
	assert(version_attribute >= 0);
	std::array<std::int64_t, 2> version{};
	assert(H5Aread(version_attribute, H5T_NATIVE_INT64, version.data()) >= 0);
	assert((version == std::array<std::int64_t, 2>{{2, 1}}));
	H5Aclose(version_attribute);
	H5Gclose(group);
	assert((Dimensions(file, "/VTKHDF/Points") == std::vector<hsize_t>{64, 3}));
	assert((Dimensions(file, "/VTKHDF/Connectivity") == std::vector<hsize_t>{64}));
	assert((Dimensions(file, "/VTKHDF/PointData/scalar") == std::vector<hsize_t>{192}));
	assert((Dimensions(file, "/VTKHDF/PointData/velocity") == std::vector<hsize_t>{192, 3}));
	assert((Dimensions(file, "/VTKHDF/Steps/CellOffsets") == std::vector<hsize_t>{3, 1}));
	const auto types = Read<std::uint8_t>(file, "/VTKHDF/Types", H5T_NATIVE_UINT8);
	assert(types == std::vector<std::uint8_t>{iga::kVtkBezierHexahedron});
	const auto times = Read<double>(file, "/VTKHDF/Steps/Values", H5T_NATIVE_DOUBLE);
	assert((times == std::vector<double>{0.0, 0.5, 1.0}));
	const auto offsets = Read<std::int64_t>(file,
		"/VTKHDF/Steps/PointDataOffsets/scalar", H5T_NATIVE_INT64);
	assert((offsets == std::vector<std::int64_t>{0, 64, 128}));
	const auto values = Read<double>(file, "/VTKHDF/PointData/scalar", H5T_NATIVE_DOUBLE);
	assert(values[0] == 0.0 && values[64] == 10.0 && values[128] == 20.0);
	H5Fclose(file);
	if (std::getenv("TUBULARFLOWIGA_KEEP_TEST_OUTPUT"))
		std::cout << path << '\n';
	else
		fs::remove_all(directory);
}

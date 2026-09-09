#include "TemporalVtkHdf.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class Fault { None, Flush, Group, File, Dataset, Count };
Fault fault = Fault::None;
hid_t leaked_dataset = -1;
int injected = 0;

bool Fail(Fault candidate)
{
	if (fault != candidate) return false;
	fault = Fault::None;
	++injected;
	return true;
}

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

template <class Function>
void Reject(Function&& function, const std::string& diagnostic)
{
	try { function(); }
	catch (const std::exception& error) {
		Require(std::string(error.what()).find(diagnostic) != std::string::npos,
			"unexpected diagnostic: "+std::string(error.what()));
		return;
	}
	throw std::runtime_error("missing rejection: "+diagnostic);
}

iga::BezierVisualizationMesh MakeMesh()
{
	iga::BezierVisualizationMesh mesh;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				mesh.signature_nodes.push_back(static_cast<std::int32_t>(mesh.points.size()));
				mesh.points.push_back({{i/3.0, j/3.0, k/3.0}});
				mesh.signature_coefficients.push_back(1.0);
				mesh.signature_offsets.push_back(mesh.points.size());
			}
	for (const auto index : iga::detail::VtkCubicHexTensorIndices())
		mesh.connectivity.push_back(index);
	mesh.offsets = {0, 64};
	mesh.types = {iga::kVtkBezierHexahedron};
	mesh.higher_order_degrees = {{{3, 3, 3}}};
	mesh.element_ids = {17};
	mesh.element_owners = {0};
	return mesh;
}

void Verify(const std::filesystem::path& path)
{
	auto file = iga::hdf_detail::RequireHandle(
		H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "reopen failed");
	auto values = iga::hdf_detail::RequireHandle(
		H5Dopen2(file.get(), "/VTKHDF/PointData/scalar", H5P_DEFAULT), H5Dclose, "dataset missing");
	auto space = iga::hdf_detail::RequireHandle(H5Dget_space(values.get()), H5Sclose, "space missing");
	Require(H5Sget_simple_extent_npoints(space.get()) == 64, "wrong field size");
	std::vector<double> data(64);
	Require(H5Dread(values.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, data.data()) >= 0, "field read failed");
	for (std::size_t i = 0; i < data.size(); ++i)
		Require(data[i] == static_cast<double>(i), "field changed after finalization retry");
	auto times = iga::hdf_detail::RequireHandle(
		H5Dopen2(file.get(), "/VTKHDF/Steps/Values", H5P_DEFAULT), H5Dclose, "time missing");
	auto time_space = iga::hdf_detail::RequireHandle(H5Dget_space(times.get()), H5Sclose, "time space missing");
	Require(H5Sget_simple_extent_npoints(time_space.get()) == 1, "wrong time count");
	double time = -1.0;
	Require(H5Dread(times.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, &time) >= 0 && time == 0.25, "time changed");
}

} // namespace

// GNU/compatible linker wrappers inject returned errors while retaining real
// HDF5 handles. No filesystem or production-code fault hooks are required.
extern "C" {
herr_t __real_H5Fflush(hid_t, H5F_scope_t);
herr_t __real_H5Gclose(hid_t);
herr_t __real_H5Fclose(hid_t);
herr_t __real_H5Dclose(hid_t);
ssize_t __real_H5Fget_obj_count(hid_t, unsigned);

herr_t __wrap_H5Fflush(hid_t id, H5F_scope_t scope)
{
	return Fail(Fault::Flush) ? -1 : __real_H5Fflush(id, scope);
}
herr_t __wrap_H5Gclose(hid_t id)
{
	return Fail(Fault::Group) ? -1 : __real_H5Gclose(id);
}
herr_t __wrap_H5Fclose(hid_t id)
{
	return Fail(Fault::File) ? -1 : __real_H5Fclose(id);
}
herr_t __wrap_H5Dclose(hid_t id)
{
	if (Fail(Fault::Dataset)) {
		leaked_dataset = id;
		return -1;
	}
	return __real_H5Dclose(id);
}
ssize_t __wrap_H5Fget_obj_count(hid_t id, unsigned types)
{
	return Fail(Fault::Count) ? -1 : __real_H5Fget_obj_count(id, types);
}
}

int main(int argc, char** argv)
{
	try {
		Require(argc == 2, "usage: vtkhdf_close_failure_test NEW_OUTPUT_DIRECTORY");
		const std::filesystem::path directory(argv[1]);
		Require(std::filesystem::create_directory(directory), "test directory must be new");
		const auto mesh = MakeMesh();
		std::vector<double> scalar(64);
		for (std::size_t i = 0; i < scalar.size(); ++i) scalar[i] = static_cast<double>(i);
		const std::vector<iga::VtkPointArray> arrays{{"scalar", 1, scalar}};
		const auto initial_objects = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
		const std::vector<std::pair<Fault, std::string>> cases{
			{Fault::Flush, "cannot finalize VTKHDF output"},
			{Fault::Group, "cannot close VTKHDF root group"},
			{Fault::File, "cannot close VTKHDF output"},
			{Fault::Dataset, "only the file handle"},
			{Fault::Count, "only the file handle"}
		};
		for (std::size_t i = 0; i < cases.size(); ++i) {
			const auto path = directory/(std::to_string(i)+".vtkhdf");
			iga::TemporalVtkHdfWriter writer(path, mesh);
			if (cases[i].first == Fault::Dataset) fault = Fault::Dataset;
			writer.Append(0.25, arrays);
			if (cases[i].first != Fault::Dataset) fault = cases[i].first;
			Reject([&] { writer.Close(); }, cases[i].second);
			Require(fault == Fault::None, "fault was not reached");
			Reject([&] { writer.Append(0.5, arrays); }, "closing or closed");
			if (leaked_dataset >= 0) {
				// The test owns recovery of the intentionally orphaned handle;
				// Close must not claim success until it has been released.
				Require(__real_H5Dclose(leaked_dataset) >= 0, "orphan recovery failed");
				leaked_dataset = -1;
			}
			writer.Close();
			writer.Close();
			Reject([&] { writer.Append(0.5, arrays); }, "closing or closed");
			Require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == initial_objects,
				"successful Close retained file objects");
			Verify(path);
			std::cout << "fault=" << i << " rejection, retry, roundtrip: PASS\n";
		}
		// Another reader of the same physical file must remain usable and must
		// not prevent this writer from closing its own file identifier.
		{
			const auto path = directory/"reader.vtkhdf";
			iga::TemporalVtkHdfWriter writer(path, mesh);
			writer.Append(0.25, arrays);
			auto reader = iga::hdf_detail::RequireHandle(
				H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "reader open failed");
			auto group = iga::hdf_detail::OpenGroup(reader.get(), "VTKHDF");
			writer.Close();
			Require(H5Iis_valid(group.get()) > 0, "writer invalidated independent reader");
			Verify(path);
		}
		Require(injected == 5, "wrong injection count");
		Require(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == initial_objects, "file objects leaked");
		std::cout << "faults=5 retries=5 roundtrips=6 independent_reader=PASS\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

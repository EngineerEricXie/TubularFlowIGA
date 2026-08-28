#ifndef IGA_TEMPORAL_VTK_HDF_HPP
#define IGA_TEMPORAL_VTK_HDF_HPP

#include "BezierVisualization.hpp"

#include <hdf5.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

namespace hdf_detail {

class Handle {
public:
	Handle() = default;
	Handle(hid_t value, herr_t (*closer)(hid_t)) : value_(value), closer_(closer) {}
	Handle(const Handle&) = delete;
	Handle& operator=(const Handle&) = delete;
	Handle(Handle&& other) noexcept : value_(other.value_), closer_(other.closer_)
	{
		other.value_ = -1;
	}
	Handle& operator=(Handle&& other) noexcept
	{
		if (this != &other) {
			Close();
			value_ = other.value_;
			closer_ = other.closer_;
			other.value_ = -1;
		}
		return *this;
	}
	~Handle() { Close(); }
	hid_t get() const { return value_; }
	explicit operator bool() const { return value_ >= 0; }
private:
	void Close()
	{
		if (value_ >= 0 && closer_) closer_(value_);
		value_ = -1;
	}
	hid_t value_ = -1;
	herr_t (*closer_)(hid_t) = nullptr;
};

inline void Require(herr_t status, const std::string& message)
{
	if (status < 0) throw std::runtime_error(message);
}

inline Handle RequireHandle(hid_t value, herr_t (*closer)(hid_t),
	const std::string& message)
{
	if (value < 0) throw std::runtime_error(message);
	return Handle(value, closer);
}

inline Handle CreateGroup(hid_t parent, const std::string& name)
{
	return RequireHandle(H5Gcreate2(parent, name.c_str(), H5P_DEFAULT,
		H5P_DEFAULT, H5P_DEFAULT), H5Gclose, "cannot create HDF5 group "+name);
}

inline Handle OpenGroup(hid_t parent, const std::string& name)
{
	return RequireHandle(H5Gopen2(parent, name.c_str(), H5P_DEFAULT), H5Gclose,
		"cannot open HDF5 group "+name);
}

inline void WriteStringAttribute(hid_t object, const std::string& name,
	const std::string& value, std::size_t minimum_size = 0)
{
	const auto size = std::max<std::size_t>({1, minimum_size, value.size()});
	auto type = RequireHandle(H5Tcopy(H5T_C_S1), H5Tclose,
		"cannot copy HDF5 string type");
	Require(H5Tset_size(type.get(), size), "cannot size HDF5 string attribute");
	Require(H5Tset_strpad(type.get(), H5T_STR_NULLPAD),
		"cannot configure HDF5 string padding");
	auto space = RequireHandle(H5Screate(H5S_SCALAR), H5Sclose,
		"cannot create HDF5 scalar space");
	auto attribute = RequireHandle(H5Acreate2(object, name.c_str(), type.get(),
		space.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
		"cannot create HDF5 attribute "+name);
	std::vector<char> buffer(size, '\0');
	std::copy(value.begin(), value.end(), buffer.begin());
	Require(H5Awrite(attribute.get(), type.get(), buffer.data()),
		"cannot write HDF5 attribute "+name);
}

inline std::string ReadStringAttribute(hid_t object, const std::string& name)
{
	auto attribute = RequireHandle(H5Aopen(object, name.c_str(), H5P_DEFAULT),
		H5Aclose, "cannot open HDF5 attribute "+name);
	auto type = RequireHandle(H5Aget_type(attribute.get()), H5Tclose,
		"cannot inspect HDF5 attribute "+name);
	const auto size = H5Tget_size(type.get());
	std::vector<char> buffer(size+1, '\0');
	Require(H5Aread(attribute.get(), type.get(), buffer.data()),
		"cannot read HDF5 attribute "+name);
	return std::string(buffer.data());
}

template <class T>
inline hid_t NativeType();
template <> inline hid_t NativeType<double>() { return H5T_NATIVE_DOUBLE; }
template <> inline hid_t NativeType<std::int32_t>() { return H5T_NATIVE_INT32; }
template <> inline hid_t NativeType<std::int64_t>() { return H5T_NATIVE_INT64; }
template <> inline hid_t NativeType<std::uint8_t>() { return H5T_NATIVE_UINT8; }
template <> inline hid_t NativeType<std::uint64_t>() { return H5T_NATIVE_UINT64; }

template <class T>
inline void WriteScalarAttribute(hid_t object, const std::string& name, T value)
{
	auto space = RequireHandle(H5Screate(H5S_SCALAR), H5Sclose,
		"cannot create scalar attribute space");
	auto attribute = RequireHandle(H5Acreate2(object, name.c_str(), NativeType<T>(),
		space.get(), H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
		"cannot create HDF5 attribute "+name);
	Require(H5Awrite(attribute.get(), NativeType<T>(), &value),
		"cannot write HDF5 attribute "+name);
}

template <class T>
inline T ReadScalarAttribute(hid_t object, const std::string& name)
{
	auto attribute = RequireHandle(H5Aopen(object, name.c_str(), H5P_DEFAULT),
		H5Aclose, "cannot open HDF5 attribute "+name);
	T value{};
	Require(H5Aread(attribute.get(), NativeType<T>(), &value),
		"cannot read HDF5 attribute "+name);
	return value;
}

template <class T>
inline void ReplaceScalarAttribute(hid_t object, const std::string& name, T value)
{
	auto attribute = RequireHandle(H5Aopen(object, name.c_str(), H5P_DEFAULT),
		H5Aclose, "cannot open HDF5 attribute "+name);
	Require(H5Awrite(attribute.get(), NativeType<T>(), &value),
		"cannot update HDF5 attribute "+name);
}

inline Handle DatasetProperties(const std::vector<hsize_t>& dimensions,
	int compression)
{
	auto properties = RequireHandle(H5Pcreate(H5P_DATASET_CREATE), H5Pclose,
		"cannot create HDF5 dataset properties");
	std::vector<hsize_t> chunk = dimensions;
	if (chunk.empty()) chunk.push_back(1);
	chunk[0] = std::max<hsize_t>(1, std::min<hsize_t>(chunk[0], 65536));
	for (auto& value : chunk) value = std::max<hsize_t>(1, value);
	Require(H5Pset_chunk(properties.get(), static_cast<int>(chunk.size()), chunk.data()),
		"cannot configure HDF5 chunks");
	if (compression > 0)
		Require(H5Pset_deflate(properties.get(), static_cast<unsigned int>(compression)),
			"cannot configure HDF5 compression");
	return properties;
}

template <class T>
inline Handle WriteFixedDataset(hid_t parent, const std::string& name,
	const T* values, const std::vector<hsize_t>& dimensions, int compression)
{
	auto space = RequireHandle(H5Screate_simple(static_cast<int>(dimensions.size()),
		dimensions.data(), nullptr), H5Sclose, "cannot create HDF5 dataset space");
	auto properties = DatasetProperties(dimensions, compression);
	auto dataset = RequireHandle(H5Dcreate2(parent, name.c_str(), NativeType<T>(),
		space.get(), H5P_DEFAULT, properties.get(), H5P_DEFAULT), H5Dclose,
		"cannot create HDF5 dataset "+name);
	Require(H5Dwrite(dataset.get(), NativeType<T>(), H5S_ALL, H5S_ALL,
		H5P_DEFAULT, values), "cannot write HDF5 dataset "+name);
	return dataset;
}

template <class T>
inline Handle CreateExpandableDataset(hid_t parent, const std::string& name,
	int components, int compression, bool force_component_dimension = false)
{
	if (components < 1) throw std::runtime_error("invalid HDF5 component count");
	std::vector<hsize_t> dimensions(
		components == 1 && !force_component_dimension ? 1 : 2, 0);
	std::vector<hsize_t> maximum(dimensions.size(), H5S_UNLIMITED);
	std::vector<hsize_t> chunk(dimensions.size(), 1);
	chunk[0] = 4096;
	if (dimensions.size() == 2) {
		dimensions[1] = static_cast<hsize_t>(components);
		maximum[1] = dimensions[1];
		chunk[1] = dimensions[1];
	}
	auto space = RequireHandle(H5Screate_simple(static_cast<int>(dimensions.size()),
		dimensions.data(), maximum.data()), H5Sclose,
		"cannot create expandable HDF5 space");
	auto properties = RequireHandle(H5Pcreate(H5P_DATASET_CREATE), H5Pclose,
		"cannot create HDF5 dataset properties");
	Require(H5Pset_chunk(properties.get(), static_cast<int>(chunk.size()), chunk.data()),
		"cannot configure expandable HDF5 chunks");
	if (compression > 0)
		Require(H5Pset_deflate(properties.get(), static_cast<unsigned int>(compression)),
			"cannot configure expandable HDF5 compression");
	return RequireHandle(H5Dcreate2(parent, name.c_str(), NativeType<T>(), space.get(),
		H5P_DEFAULT, properties.get(), H5P_DEFAULT), H5Dclose,
		"cannot create expandable HDF5 dataset "+name);
}

inline std::vector<hsize_t> DatasetDimensions(hid_t dataset)
{
	auto space = RequireHandle(H5Dget_space(dataset), H5Sclose,
		"cannot inspect HDF5 dataset space");
	const auto rank = H5Sget_simple_extent_ndims(space.get());
	if (rank < 1) throw std::runtime_error("invalid HDF5 dataset rank");
	std::vector<hsize_t> result(static_cast<std::size_t>(rank));
	Require(H5Sget_simple_extent_dims(space.get(), result.data(), nullptr),
		"cannot inspect HDF5 dataset dimensions");
	return result;
}

inline void ResizeRows(hid_t dataset, hsize_t rows)
{
	auto dimensions = DatasetDimensions(dataset);
	dimensions[0] = rows;
	Require(H5Dset_extent(dataset, dimensions.data()),
		"cannot resize HDF5 dataset");
}

template <class T>
inline void WriteRows(hid_t dataset, hsize_t row_offset, hsize_t rows,
	int components, const T* values)
{
	auto dimensions = DatasetDimensions(dataset);
	const auto required = row_offset+rows;
	if (dimensions[0] < required) ResizeRows(dataset, required);
	auto file_space = RequireHandle(H5Dget_space(dataset), H5Sclose,
		"cannot open HDF5 file space");
	std::vector<hsize_t> start(dimensions.size(), 0);
	std::vector<hsize_t> count(dimensions.size(), 1);
	start[0] = row_offset;
	count[0] = rows;
	if (dimensions.size() == 2) count[1] = static_cast<hsize_t>(components);
	Require(H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET,
		start.data(), nullptr, count.data(), nullptr),
		"cannot select HDF5 output rows");
	auto memory_space = RequireHandle(H5Screate_simple(static_cast<int>(count.size()),
		count.data(), nullptr), H5Sclose, "cannot create HDF5 memory space");
	if (H5Dwrite(dataset, NativeType<T>(), memory_space.get(), file_space.get(),
		H5P_DEFAULT, values) < 0) {
		std::array<char, 512> name{};
		H5Iget_name(dataset, name.data(), name.size());
		throw std::runtime_error("cannot write HDF5 output rows: "+std::string(name.data()));
	}
}

template <class T>
inline T ReadRowValue(hid_t dataset, hsize_t row)
{
	auto file_space = RequireHandle(H5Dget_space(dataset), H5Sclose,
		"cannot open HDF5 file space");
	const hsize_t start[1] = {row};
	const hsize_t count[1] = {1};
	Require(H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET,
		start, nullptr, count, nullptr), "cannot select HDF5 row");
	auto memory_space = RequireHandle(H5Screate_simple(1, count, nullptr), H5Sclose,
		"cannot create HDF5 row space");
	T value{};
	Require(H5Dread(dataset, NativeType<T>(), memory_space.get(), file_space.get(),
		H5P_DEFAULT, &value), "cannot read HDF5 row");
	return value;
}

inline std::uint64_t HashBytes(std::uint64_t hash, const void* data, std::size_t bytes)
{
	const auto* characters = static_cast<const unsigned char*>(data);
	for (std::size_t index = 0; index < bytes; ++index) {
		hash ^= characters[index];
		hash *= 1099511628211ULL;
	}
	return hash;
}

inline std::uint64_t GeometryHash(const BezierVisualizationMesh& mesh)
{
	std::uint64_t hash = 1469598103934665603ULL;
	for (const auto& point : mesh.points)
		hash = HashBytes(hash, point.data(), sizeof(point));
	hash = HashBytes(hash, mesh.connectivity.data(),
		mesh.connectivity.size()*sizeof(mesh.connectivity.front()));
	hash = HashBytes(hash, mesh.offsets.data(),
		mesh.offsets.size()*sizeof(mesh.offsets.front()));
	return hash;
}

inline std::string ArraySchema(const std::vector<VtkPointArray>& arrays)
{
	std::set<std::string> names;
	std::ostringstream result;
	for (const auto& array : arrays) {
		if (array.name.empty() || array.name.find('/') != std::string::npos
			|| array.name.find('\n') != std::string::npos
			|| array.name.find('\t') != std::string::npos)
			throw std::runtime_error("invalid VTKHDF point-array name: "+array.name);
		if (!names.insert(array.name).second)
			throw std::runtime_error("duplicate VTKHDF point-array name: "+array.name);
		result << array.name << '\t' << array.components << '\n';
	}
	return result.str();
}

} // namespace hdf_detail

class TemporalVtkHdfWriter {
public:
	TemporalVtkHdfWriter(const std::filesystem::path& path,
		const BezierVisualizationMesh& mesh, bool resume = false, int compression = 4)
		: path_(path), mesh_(mesh), compression_(compression)
	{
		if (compression_ < 0 || compression_ > 9)
			throw std::runtime_error("VTKHDF compression must be between 0 and 9");
		H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
		if (!path_.parent_path().empty())
			std::filesystem::create_directories(path_.parent_path());
		if (resume && std::filesystem::exists(path_)) {
			file_ = hdf_detail::RequireHandle(H5Fopen(path_.string().c_str(),
				H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose,
				"cannot open VTKHDF output "+path_.string());
			OpenExisting();
		} else {
			file_ = hdf_detail::RequireHandle(H5Fcreate(path_.string().c_str(),
				H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose,
				"cannot create VTKHDF output "+path_.string());
			CreateFile();
		}
	}

	void Append(double physical_time, const std::vector<VtkPointArray>& control_arrays)
	{
		if (!std::isfinite(physical_time))
			throw std::runtime_error("VTKHDF physical time is not finite");
		const auto arrays = ExtractBezierPointArrays(mesh_, control_arrays);
		EnsureArraySchema(arrays);
		auto steps = hdf_detail::OpenGroup(root_.get(), "Steps");
		auto values = hdf_detail::RequireHandle(H5Dopen2(steps.get(), "Values", H5P_DEFAULT),
			H5Dclose, "cannot open VTKHDF time values");
		auto step_count = static_cast<std::uint64_t>(
			hdf_detail::ReadScalarAttribute<std::int64_t>(steps.get(), "NSteps"));
		bool replace = false;
		if (step_count > 0) {
			const auto last = hdf_detail::ReadRowValue<double>(values.get(), step_count-1);
			const auto tolerance = 1e-12*std::max({1.0, std::abs(last), std::abs(physical_time)});
			if (physical_time < last-tolerance)
				throw std::runtime_error("VTKHDF time steps must be monotone");
			replace = std::abs(physical_time-last) <= tolerance;
		}
		const auto step = replace ? step_count-1 : step_count;
		auto point_data = hdf_detail::OpenGroup(root_.get(), "PointData");
		for (const auto& array : arrays) {
			auto dataset = hdf_detail::RequireHandle(H5Dopen2(point_data.get(),
				array.name.c_str(), H5P_DEFAULT), H5Dclose,
				"cannot open VTKHDF point array "+array.name);
			hdf_detail::WriteRows(dataset.get(),
				static_cast<hsize_t>(step*mesh_.points.size()), mesh_.points.size(),
				array.components, array.values.data());
		}
		if (!replace) {
			hdf_detail::WriteRows(values.get(), step, 1, 1, &physical_time);
			AppendStepMetadata(steps.get(), step, arrays);
			hdf_detail::ReplaceScalarAttribute<std::int64_t>(steps.get(), "NSteps",
				static_cast<std::int64_t>(step_count+1));
		}
		hdf_detail::Require(H5Fflush(file_.get(), H5F_SCOPE_GLOBAL),
			"cannot flush VTKHDF output "+path_.string());
	}

	const BezierGeometryValidation& validation() const { return mesh_.validation; }
	const std::filesystem::path& path() const { return path_; }

private:
	void CreateFile()
	{
		root_ = hdf_detail::CreateGroup(file_.get(), "VTKHDF");
		const std::array<std::int64_t, 2> version{{2, 1}};
		auto version_space = hdf_detail::RequireHandle(H5Screate_simple(1,
			std::array<hsize_t, 1>{{2}}.data(), nullptr), H5Sclose,
			"cannot create VTKHDF version space");
		auto version_attribute = hdf_detail::RequireHandle(H5Acreate2(root_.get(),
			"Version", H5T_NATIVE_INT64, version_space.get(), H5P_DEFAULT, H5P_DEFAULT),
			H5Aclose, "cannot create VTKHDF version");
		hdf_detail::Require(H5Awrite(version_attribute.get(), H5T_NATIVE_INT64,
			version.data()), "cannot write VTKHDF version");
		hdf_detail::WriteStringAttribute(root_.get(), "Type", "UnstructuredGrid");
		const std::int64_t points = static_cast<std::int64_t>(mesh_.points.size());
		const std::int64_t cells = static_cast<std::int64_t>(mesh_.types.size());
		const std::int64_t connectivity = static_cast<std::int64_t>(mesh_.connectivity.size());
		hdf_detail::WriteFixedDataset(root_.get(), "NumberOfPoints", &points, {1}, 0);
		hdf_detail::WriteFixedDataset(root_.get(), "NumberOfCells", &cells, {1}, 0);
		hdf_detail::WriteFixedDataset(root_.get(), "NumberOfConnectivityIds",
			&connectivity, {1}, 0);
		hdf_detail::WriteFixedDataset(root_.get(), "Points", mesh_.points.front().data(),
			{mesh_.points.size(), 3}, compression_);
		hdf_detail::WriteFixedDataset(root_.get(), "Connectivity", mesh_.connectivity.data(),
			{mesh_.connectivity.size()}, compression_);
		hdf_detail::WriteFixedDataset(root_.get(), "Offsets", mesh_.offsets.data(),
			{mesh_.offsets.size()}, compression_);
		hdf_detail::WriteFixedDataset(root_.get(), "Types", mesh_.types.data(),
			{mesh_.types.size()}, compression_);
		auto point_data = hdf_detail::CreateGroup(root_.get(), "PointData");
		(void)point_data;
		auto cell_data = hdf_detail::CreateGroup(root_.get(), "CellData");
		auto degrees = hdf_detail::WriteFixedDataset(cell_data.get(), "HigherOrderDegrees",
			mesh_.higher_order_degrees.front().data(), {mesh_.higher_order_degrees.size(), 3},
			compression_);
		hdf_detail::WriteStringAttribute(degrees.get(), "Attribute", "HigherOrderDegrees");
		auto element_ids = hdf_detail::WriteFixedDataset(cell_data.get(), "element_id",
			mesh_.element_ids.data(), {mesh_.element_ids.size()}, compression_);
		hdf_detail::WriteStringAttribute(element_ids.get(), "Attribute", "GlobalIds");
		hdf_detail::WriteFixedDataset(cell_data.get(), "metis_owner",
			mesh_.element_owners.data(), {mesh_.element_owners.size()}, compression_);
		auto steps = hdf_detail::CreateGroup(root_.get(), "Steps");
		hdf_detail::WriteScalarAttribute<std::int64_t>(steps.get(), "NSteps", 0);
		hdf_detail::CreateExpandableDataset<double>(steps.get(), "Values", 1, 0);
		hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(), "PartOffsets", 1, 0);
		hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(), "NumberOfParts", 1, 0);
		hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(), "PointOffsets", 1, 0);
		hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(),
			"CellOffsets", 1, 0, true);
		hdf_detail::CreateExpandableDataset<std::int64_t>(steps.get(),
			"ConnectivityIdOffsets", 1, 0, true);
		hdf_detail::CreateGroup(steps.get(), "PointDataOffsets");
		auto metadata = hdf_detail::CreateGroup(file_.get(), "TubularFlowIGA");
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(), "GeometryHash",
			hdf_detail::GeometryHash(mesh_));
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(), "UniquePoints",
			mesh_.points.size());
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(), "Elements",
			mesh_.types.size());
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(),
			"SharedPointReferences", mesh_.validation.shared_point_references);
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(),
			"SignatureCoordinateRepairs", mesh_.validation.signature_coordinate_repairs);
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(),
			"OverlapCandidates", mesh_.validation.overlap_candidates);
		hdf_detail::WriteScalarAttribute<std::uint64_t>(metadata.get(),
			"OverlappingElementPairs", mesh_.validation.overlapping_element_pairs);
	}

	void OpenExisting()
	{
		root_ = hdf_detail::OpenGroup(file_.get(), "VTKHDF");
		if (hdf_detail::ReadStringAttribute(root_.get(), "Type") != "UnstructuredGrid")
			throw std::runtime_error("existing VTKHDF output is not an UnstructuredGrid");
		auto metadata = hdf_detail::OpenGroup(file_.get(), "TubularFlowIGA");
		if (hdf_detail::ReadScalarAttribute<std::uint64_t>(metadata.get(), "GeometryHash")
			!= hdf_detail::GeometryHash(mesh_))
			throw std::runtime_error("existing VTKHDF geometry does not match the current .ntiga database");
		RecoverInterruptedAppend();
	}

	void RecoverInterruptedAppend()
	{
		auto steps = hdf_detail::OpenGroup(root_.get(), "Steps");
		const auto count = static_cast<hsize_t>(
			hdf_detail::ReadScalarAttribute<std::int64_t>(steps.get(), "NSteps"));
		for (const auto* name : {"Values", "PartOffsets", "NumberOfParts", "PointOffsets",
			"CellOffsets", "ConnectivityIdOffsets"}) {
			auto dataset = hdf_detail::RequireHandle(H5Dopen2(steps.get(), name, H5P_DEFAULT),
				H5Dclose, std::string("cannot open VTKHDF step array ")+name);
			hdf_detail::ResizeRows(dataset.get(), count);
		}
		auto metadata = hdf_detail::OpenGroup(file_.get(), "TubularFlowIGA");
		if (H5Aexists(metadata.get(), "PointArraySchema") <= 0) return;
		const auto schema = hdf_detail::ReadStringAttribute(metadata.get(), "PointArraySchema");
		std::istringstream input(schema);
		std::string line;
		auto point_data = hdf_detail::OpenGroup(root_.get(), "PointData");
		auto offsets = hdf_detail::OpenGroup(steps.get(), "PointDataOffsets");
		while (std::getline(input, line)) {
			const auto separator = line.find('\t');
			if (separator == std::string::npos) continue;
			const auto name = line.substr(0, separator);
			auto values = hdf_detail::RequireHandle(H5Dopen2(point_data.get(), name.c_str(),
				H5P_DEFAULT), H5Dclose, "cannot open VTKHDF point array "+name);
			hdf_detail::ResizeRows(values.get(), count*mesh_.points.size());
			auto positions = hdf_detail::RequireHandle(H5Dopen2(offsets.get(), name.c_str(),
				H5P_DEFAULT), H5Dclose, "cannot open VTKHDF point offset "+name);
			hdf_detail::ResizeRows(positions.get(), count);
		}
	}

	void EnsureArraySchema(const std::vector<VtkPointArray>& arrays)
	{
		const auto schema = hdf_detail::ArraySchema(arrays);
		auto metadata = hdf_detail::OpenGroup(file_.get(), "TubularFlowIGA");
		if (H5Aexists(metadata.get(), "PointArraySchema") > 0) {
			if (hdf_detail::ReadStringAttribute(metadata.get(), "PointArraySchema") != schema)
				throw std::runtime_error("VTKHDF point-array schema changed while appending");
			return;
		}
		hdf_detail::WriteStringAttribute(metadata.get(), "PointArraySchema", schema);
		auto point_data = hdf_detail::OpenGroup(root_.get(), "PointData");
		auto steps = hdf_detail::OpenGroup(root_.get(), "Steps");
		auto offsets = hdf_detail::OpenGroup(steps.get(), "PointDataOffsets");
		for (const auto& array : arrays) {
			auto dataset = hdf_detail::CreateExpandableDataset<double>(point_data.get(),
				array.name, array.components, compression_);
			if (array.components == 3 && array.name == "velocity")
				hdf_detail::WriteStringAttribute(dataset.get(), "Attribute", "Vectors");
			hdf_detail::CreateExpandableDataset<std::int64_t>(offsets.get(),
				array.name, 1, 0);
		}
	}

	void AppendStepMetadata(hid_t steps, std::uint64_t step,
		const std::vector<VtkPointArray>& arrays)
	{
		const std::int64_t zero = 0;
		const std::int64_t one = 1;
		for (const auto* name : {"PartOffsets", "PointOffsets", "CellOffsets",
			"ConnectivityIdOffsets"}) {
			auto dataset = hdf_detail::RequireHandle(H5Dopen2(steps, name, H5P_DEFAULT),
				H5Dclose, std::string("cannot open VTKHDF step array ")+name);
			hdf_detail::WriteRows(dataset.get(), step, 1, 1, &zero);
		}
		auto number_of_parts = hdf_detail::RequireHandle(H5Dopen2(steps,
			"NumberOfParts", H5P_DEFAULT), H5Dclose,
			"cannot open VTKHDF NumberOfParts");
		hdf_detail::WriteRows(number_of_parts.get(), step, 1, 1, &one);
		auto offsets = hdf_detail::OpenGroup(steps, "PointDataOffsets");
		const auto point_offset = static_cast<std::int64_t>(step*mesh_.points.size());
		for (const auto& array : arrays) {
			auto dataset = hdf_detail::RequireHandle(H5Dopen2(offsets.get(),
				array.name.c_str(), H5P_DEFAULT), H5Dclose,
				"cannot open VTKHDF point-data offset "+array.name);
			hdf_detail::WriteRows(dataset.get(), step, 1, 1, &point_offset);
		}
	}

	std::filesystem::path path_;
	const BezierVisualizationMesh& mesh_;
	int compression_ = 4;
	hdf_detail::Handle file_;
	hdf_detail::Handle root_;
};

inline std::filesystem::path VtkHdfPath(const std::filesystem::path& base)
{
	return base.parent_path()/(base.stem().string()+".vtkhdf");
}

} // namespace iga

#endif

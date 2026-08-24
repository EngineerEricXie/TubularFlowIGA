#ifndef IGA_CUDA_DEVICE_MESH_HPP
#define IGA_CUDA_DEVICE_MESH_HPP

#include "CudaRuntime.hpp"
#include "IgaDatabase.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga::cuda {

struct FlatMesh {
	std::uint64_t nodes = 0;
	std::vector<int> element_offsets;
	std::vector<int> connectivity;
	std::vector<int> extraction_offsets;
	std::vector<std::uint8_t> extraction_columns;
	std::vector<double> extraction_values;
	std::vector<double> bezier_points;
	std::vector<int> boundary_labels;
	int maximum_basis = 0;

	explicit FlatMesh(iga::Database& database)
	{
		nodes = database.header().nodes;
		if (nodes > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
			throw std::runtime_error("CUDA solver currently requires fewer than INT_MAX nodes");
		const auto count = database.header().elements;
		if (count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
			throw std::runtime_error("CUDA solver currently requires fewer than INT_MAX elements");
		element_offsets.reserve(static_cast<std::size_t>(count) + 1);
		element_offsets.push_back(0);
		extraction_offsets.push_back(0);
		bezier_points.reserve(static_cast<std::size_t>(count) * 64 * 3);
		for (std::uint64_t index = 0; index < count; ++index) {
			auto element = database.Load(index);
			maximum_basis = std::max(maximum_basis, static_cast<int>(element.connectivity.size()));
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				connectivity.push_back(element.connectivity[a]);
				for (std::uint8_t column = 0; column < 64; ++column) {
					const auto coefficient = element.extraction[a][column];
					if (coefficient == 0.0) continue;
					extraction_columns.push_back(column);
					extraction_values.push_back(coefficient);
				}
				extraction_offsets.push_back(CheckedInt(extraction_values.size(), "extraction entries"));
			}
			element_offsets.push_back(CheckedInt(connectivity.size(), "element basis entries"));
			for (const auto& point : element.bezier_points)
				for (double coordinate : point) bezier_points.push_back(coordinate);
			boundary_labels.insert(boundary_labels.end(),
				element.boundary_labels.begin(), element.boundary_labels.end());
		}
		if (maximum_basis > 80)
			throw std::runtime_error("CUDA kernels support at most 80 basis functions per element");
	}

	std::size_t elements() const { return element_offsets.size() - 1; }

private:
	static int CheckedInt(std::size_t value, const char* name)
	{
		if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
			throw std::runtime_error(std::string(name) + " exceed INT_MAX");
		return static_cast<int>(value);
	}
};

struct DeviceMeshView {
	int nodes = 0;
	int elements = 0;
	const int* element_offsets = nullptr;
	const int* connectivity = nullptr;
	const int* extraction_offsets = nullptr;
	const std::uint8_t* extraction_columns = nullptr;
	const double* extraction_values = nullptr;
	const double* bezier_points = nullptr;
};

class DeviceMesh {
public:
	explicit DeviceMesh(const FlatMesh& host)
		: element_offsets_(host.element_offsets.size()),
		  connectivity_(host.connectivity.size()),
		  extraction_offsets_(host.extraction_offsets.size()),
		  extraction_columns_(host.extraction_columns.size()),
		  extraction_values_(host.extraction_values.size()),
		  bezier_points_(host.bezier_points.size())
	{
		element_offsets_.CopyFromHost(host.element_offsets.data(), host.element_offsets.size());
		connectivity_.CopyFromHost(host.connectivity.data(), host.connectivity.size());
		extraction_offsets_.CopyFromHost(host.extraction_offsets.data(), host.extraction_offsets.size());
		extraction_columns_.CopyFromHost(host.extraction_columns.data(), host.extraction_columns.size());
		extraction_values_.CopyFromHost(host.extraction_values.data(), host.extraction_values.size());
		bezier_points_.CopyFromHost(host.bezier_points.data(), host.bezier_points.size());
		view_ = {static_cast<int>(host.nodes), static_cast<int>(host.elements()),
			element_offsets_.data(), connectivity_.data(), extraction_offsets_.data(),
			extraction_columns_.data(), extraction_values_.data(), bezier_points_.data()};
	}

	const DeviceMeshView& view() const { return view_; }
	std::size_t bytes() const
	{
		return element_offsets_.bytes() + connectivity_.bytes() + extraction_offsets_.bytes()
			+ extraction_columns_.bytes() + extraction_values_.bytes() + bezier_points_.bytes();
	}

private:
	DeviceBuffer<int> element_offsets_, connectivity_, extraction_offsets_;
	DeviceBuffer<std::uint8_t> extraction_columns_;
	DeviceBuffer<double> extraction_values_, bezier_points_;
	DeviceMeshView view_;
};

} // namespace iga::cuda

#endif

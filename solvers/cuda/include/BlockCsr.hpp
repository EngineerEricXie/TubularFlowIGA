#ifndef IGA_CUDA_BLOCK_CSR_HPP
#define IGA_CUDA_BLOCK_CSR_HPP

#include "CudaRuntime.hpp"
#include "DeviceMesh.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga::cuda {

struct BlockPattern {
	std::vector<int> row_offsets;
	std::vector<int> columns;
	std::vector<int> diagonal;

	explicit BlockPattern(const FlatMesh& mesh)
	{
		std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(mesh.nodes));
		for (std::size_t element = 0; element < mesh.elements(); ++element) {
			const int begin = mesh.element_offsets[element];
			const int end = mesh.element_offsets[element + 1];
			for (int a = begin; a < end; ++a) {
				auto& row = adjacency[static_cast<std::size_t>(mesh.connectivity[a])];
				row.insert(row.end(), mesh.connectivity.begin() + begin, mesh.connectivity.begin() + end);
			}
		}
		row_offsets.reserve(static_cast<std::size_t>(mesh.nodes) + 1);
		diagonal.resize(static_cast<std::size_t>(mesh.nodes), -1);
		row_offsets.push_back(0);
		for (std::size_t node = 0; node < adjacency.size(); ++node) {
			auto& row = adjacency[node];
			std::sort(row.begin(), row.end());
			row.erase(std::unique(row.begin(), row.end()), row.end());
			if (!std::binary_search(row.begin(), row.end(), static_cast<int>(node)))
				throw std::runtime_error("block pattern is missing a diagonal entry");
			if (columns.size() + row.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::runtime_error("block CSR has more than INT_MAX blocks");
			const auto position = std::lower_bound(row.begin(), row.end(), static_cast<int>(node)) - row.begin();
			diagonal[node] = static_cast<int>(columns.size() + static_cast<std::size_t>(position));
			columns.insert(columns.end(), row.begin(), row.end());
			row_offsets.push_back(static_cast<int>(columns.size()));
		}
	}
};

struct DevicePatternView {
	int nodes = 0;
	int blocks = 0;
	const int* row_offsets = nullptr;
	const int* columns = nullptr;
	const int* diagonal = nullptr;
};

class DevicePattern {
public:
	explicit DevicePattern(const BlockPattern& host)
		: row_offsets_(host.row_offsets.size()), columns_(host.columns.size()), diagonal_(host.diagonal.size())
	{
		row_offsets_.CopyFromHost(host.row_offsets.data(), host.row_offsets.size());
		columns_.CopyFromHost(host.columns.data(), host.columns.size());
		diagonal_.CopyFromHost(host.diagonal.data(), host.diagonal.size());
		view_ = {static_cast<int>(host.diagonal.size()), static_cast<int>(host.columns.size()),
			row_offsets_.data(), columns_.data(), diagonal_.data()};
	}

	const DevicePatternView& view() const { return view_; }
	std::size_t bytes() const { return row_offsets_.bytes() + columns_.bytes() + diagonal_.bytes(); }

private:
	DeviceBuffer<int> row_offsets_, columns_, diagonal_;
	DevicePatternView view_;
};

template <int Fields>
class BlockMatrix {
public:
	explicit BlockMatrix(const DevicePattern& pattern)
		: pattern_(pattern.view()), values_(static_cast<std::size_t>(pattern_.blocks) * Fields * Fields) {}

	void Clear() { values_.Clear(); }
	const DevicePatternView& pattern() const { return pattern_; }
	double* values() { return values_.data(); }
	const double* values() const { return values_.data(); }
	std::size_t bytes() const { return values_.bytes(); }

private:
	DevicePatternView pattern_;
	DeviceBuffer<double> values_;
};

} // namespace iga::cuda

#endif

#ifndef IGA_CUDA_AGGREGATION_AMG_HPP
#define IGA_CUDA_AGGREGATION_AMG_HPP

#include "ScalarKrylov.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

// Plain (unsmoothed) aggregation algebraic multigrid for a scalar CSR matrix
// whose unknowns come in blocks of `block` components per node (3 for the
// velocity block, 1 for a pressure Schur complement). The hierarchy is built
// on the host from the first matrix values: nodes are aggregated along strong
// connections, prolongation is piecewise constant, and every coarse operator
// A_c = P' A P is a sum of fine entries, so Update() refreshes the hierarchy
// from new values on the device without touching the host again. Smoothing
// is damped Jacobi; the coarsest operator (a few hundred unknowns) is
// inverted densely on the host at every Update() and applied with one GEMV.
// Cycle() applies V-cycles as a preconditioner and falls back to Jacobi
// scaling if a cycle produces a non-finite vector.

namespace iga::cuda {

struct AggregationAmgOptions {
	int pre_sweeps = 1;
	int post_sweeps = 1;
	double omega = 0.6;
	// Node j is a strong neighbour of i when |a_ij| >= theta * max_k |a_ik|.
	double strength_threshold = 0.25;
	// Plain aggregation under-corrects; scaling the coarse correction helps.
	double over_correction = 1.3;
	int coarsest_size = 300;
	int maximum_levels = 25;
};

__global__ void SumFineEntriesKernel(int fine_entries, const int* coarse_entry, const double* fine, double* coarse)
{
	const int entry = blockIdx.x*blockDim.x+threadIdx.x;
	if (entry < fine_entries) atomicAdd(coarse+coarse_entry[entry], fine[entry]);
}

__global__ void ResidualKernel(DevicePatternView pattern, const double* values, const double* b, const double* x,
	double* r)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= pattern.nodes) return;
	double sum = 0.0;
	for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry)
		sum += values[entry]*x[pattern.columns[entry]];
	r[row] = b[row]-sum;
}

__global__ void DampedJacobiKernel(int size, double omega, const double* inverse, const double* r, double* x)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size) x[row] += omega*inverse[row]*r[row];
}

__global__ void RestrictKernel(int fine_size, const int* aggregate, const double* r, double* coarse)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < fine_size) atomicAdd(coarse+aggregate[row], r[row]);
}

__global__ void ProlongKernel(int fine_size, const int* aggregate, double scale, const double* coarse, double* x)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < fine_size) x[row] += scale*coarse[aggregate[row]];
}

__global__ void CountNonFiniteKernel(int size, const double* x, unsigned int* count)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size && !isfinite(x[row])) atomicAdd(count, 1u);
}

__global__ void ScaleIntoKernel(int size, const double* inverse, const double* b, double* x)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row < size) x[row] = inverse[row]*b[row];
}

class AggregationAmg {
public:
	// `fine` and `host_pattern` describe the same matrix; `fine` stays owned
	// by the caller and must outlive the hierarchy.
	AggregationAmg(const BlockPattern& host_pattern, const BlockMatrix<1>& fine, int block,
		AggregationAmgOptions options = {})
		: fine_(fine), host_(host_pattern), block_(block), options_(options),
		  size_(static_cast<int>(host_pattern.diagonal.size())),
		  fine_inverse_(host_pattern.diagonal.size()), fine_r_(host_pattern.diagonal.size()), flags_(1)
	{
		Check(cublasSetPointerMode(blas_, CUBLAS_POINTER_MODE_HOST), "cublasSetPointerMode AMG");
		if (block <= 0 || size_ == 0 || size_%block != 0)
			throw std::invalid_argument("aggregation AMG needs a positive block size dividing the system");
		if (options.pre_sweeps < 0 || options.post_sweeps < 0 || !(options.omega > 0.0))
			throw std::invalid_argument("aggregation AMG smoother options are invalid");
	}

	int Levels() const { return static_cast<int>(levels_.size())+1; }
	int CoarsestSize() const { return levels_.empty() ? size_ : levels_.back()->coarse_size; }
	int FallbackCount() const { return fallbacks_; }

	// Build the hierarchy on the first call, then recompute every coarse
	// operator and Jacobi inverse from the current fine values.
	void Update()
	{
		if (!built_) Build();
		const DevicePatternView* pattern = &fine_.pattern();
		const double* values = fine_.values();
		BuildInverse(*pattern, values, fine_inverse_.data());
		for (auto& level : levels_) {
			level->matrix->Clear();
			const int entries = pattern->blocks;
			SumFineEntriesKernel<<<(entries+255)/256,256>>>(entries, level->coarse_entry.data(), values,
				level->matrix->values());
			pattern = &level->matrix->pattern();
			values = level->matrix->values();
			BuildInverse(*pattern, values, level->inverse.data());
		}
		CheckKernel("aggregation AMG update");
		InvertCoarsest(*pattern, values);
	}

	// x = M^-1 b with `cycles` V-cycles from a zero start.
	void Cycle(const double* b, double* x, int cycles)
	{
		if (!built_) throw std::logic_error("aggregation AMG cycles need Update() first");
		Check(cudaMemset(x, 0, static_cast<std::size_t>(size_)*sizeof(double)), "clear AMG solution");
		for (int cycle = 0; cycle < cycles; ++cycle)
			Visit(0, fine_.pattern(), fine_.values(), fine_inverse_.data(), b, x);
		flags_.Clear();
		CountNonFiniteKernel<<<(size_+255)/256,256>>>(size_, x, flags_.data());
		unsigned int bad = 0;
		flags_.CopyToHost(&bad, 1);
		if (bad) {
			++fallbacks_;
			ScaleIntoKernel<<<(size_+255)/256,256>>>(size_, fine_inverse_.data(), b, x);
		}
		CheckKernel("aggregation AMG cycle");
	}

private:
	struct Level {
		Level(int fine_size, int coarse_size, const std::vector<int>& host_aggregate, const BlockPattern& coarse,
			const std::vector<int>& host_coarse_entry)
			: fine_size(fine_size), coarse_size(coarse_size), aggregate(host_aggregate.size()),
			  coarse_entry(host_coarse_entry.size()), pattern(std::make_unique<DevicePattern>(coarse)),
			  matrix(std::make_unique<BlockMatrix<1>>(*pattern)), inverse(coarse_size), x(coarse_size), b(coarse_size),
			  r(fine_size)
		{
			aggregate.CopyFromHost(host_aggregate.data(), host_aggregate.size());
			coarse_entry.CopyFromHost(host_coarse_entry.data(), host_coarse_entry.size());
		}

		int fine_size, coarse_size;
		DeviceBuffer<int> aggregate, coarse_entry;
		std::unique_ptr<DevicePattern> pattern;
		std::unique_ptr<BlockMatrix<1>> matrix;
		DeviceBuffer<double> inverse, x, b, r;
	};

	void Build()
	{
		std::vector<double> values(static_cast<std::size_t>(host_.columns.size()));
		Check(cudaMemcpy(values.data(), fine_.values(), values.size()*sizeof(double), cudaMemcpyDeviceToHost),
			"copy AMG fine values");
		BlockPattern pattern = host_;
		int size = size_;
		for (int level = 0; level < options_.maximum_levels && size > options_.coarsest_size; ++level) {
			std::vector<int> aggregate;
			const int coarse_size = Aggregate(pattern, values, aggregate);
			if (coarse_size*4 > size*3) break;  // too little coarsening to pay off
			BlockPattern coarse;
			std::vector<int> coarse_entry;
			CoarsePattern(pattern, aggregate, coarse_size, coarse, coarse_entry);
			std::vector<double> coarse_values(coarse.columns.size(), 0.0);
			for (std::size_t entry = 0; entry < values.size(); ++entry)
				coarse_values[static_cast<std::size_t>(coarse_entry[entry])] += values[entry];
			levels_.push_back(std::make_unique<Level>(size, coarse_size, aggregate, coarse, coarse_entry));
			pattern = std::move(coarse);
			values = std::move(coarse_values);
			size = coarse_size;
		}
		coarsest_size_ = size;
		coarsest_host_pattern_ = pattern;
		coarsest_host_values_.resize(pattern.columns.size());
		coarsest_dense_.resize(static_cast<std::size_t>(size)*size);
		coarsest_inverse_ = std::make_unique<DeviceBuffer<double>>(coarsest_dense_.size());
		built_ = true;
	}

	// Greedy aggregation along strong connections of the node graph; returns
	// the coarse unknown count and the coarse unknown of every fine unknown.
	int Aggregate(const BlockPattern& pattern, const std::vector<double>& values, std::vector<int>& aggregate) const
	{
		const int size = static_cast<int>(pattern.diagonal.size()), nodes = size/block_;
		// Strong neighbours per node from the first component's row.
		std::vector<std::vector<int>> strong(static_cast<std::size_t>(nodes));
		for (int node = 0; node < nodes; ++node) {
			const int row = node*block_;
			double maximum = 0.0;
			for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry)
				if (pattern.columns[entry]/block_ != node)
					maximum = std::max(maximum, std::abs(values[static_cast<std::size_t>(entry)]));
			for (int entry = pattern.row_offsets[row]; entry < pattern.row_offsets[row+1]; ++entry) {
				const int neighbour = pattern.columns[entry]/block_;
				if (neighbour != node && pattern.columns[entry]%block_ == 0
					&& std::abs(values[static_cast<std::size_t>(entry)]) >= options_.strength_threshold*maximum)
					strong[node].push_back(neighbour);
			}
		}
		std::vector<int> node_aggregate(nodes, -1);
		int count = 0;
		for (int node = 0; node < nodes; ++node) {
			if (node_aggregate[node] >= 0) continue;
			bool free_neighbour = false;
			for (const int neighbour : strong[node]) free_neighbour |= node_aggregate[neighbour] < 0;
			if (!free_neighbour) {
				// Join a strongly connected aggregate rather than stay a singleton.
				for (const int neighbour : strong[node])
					if (node_aggregate[neighbour] >= 0) { node_aggregate[node] = node_aggregate[neighbour]; break; }
				if (node_aggregate[node] < 0) node_aggregate[node] = count++;
				continue;
			}
			node_aggregate[node] = count;
			for (const int neighbour : strong[node])
				if (node_aggregate[neighbour] < 0) node_aggregate[neighbour] = count;
			++count;
		}
		aggregate.resize(size);
		for (int row = 0; row < size; ++row) aggregate[row] = node_aggregate[row/block_]*block_+row%block_;
		return count*block_;
	}

	static void CoarsePattern(const BlockPattern& fine, const std::vector<int>& aggregate, int coarse_size,
		BlockPattern& coarse, std::vector<int>& coarse_entry)
	{
		const int size = static_cast<int>(fine.diagonal.size());
		std::vector<std::vector<int>> rows(static_cast<std::size_t>(coarse_size));
		for (int row = 0; row < size; ++row)
			for (int entry = fine.row_offsets[row]; entry < fine.row_offsets[row+1]; ++entry)
				rows[aggregate[row]].push_back(aggregate[fine.columns[entry]]);
		coarse.row_offsets.assign(1, 0);
		coarse.diagonal.assign(static_cast<std::size_t>(coarse_size), -1);
		for (int row = 0; row < coarse_size; ++row) {
			auto& columns = rows[row];
			std::sort(columns.begin(), columns.end());
			columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
			coarse.diagonal[row] = static_cast<int>(coarse.columns.size()
				+(std::lower_bound(columns.begin(), columns.end(), row)-columns.begin()));
			coarse.columns.insert(coarse.columns.end(), columns.begin(), columns.end());
			coarse.row_offsets.push_back(static_cast<int>(coarse.columns.size()));
		}
		coarse_entry.resize(fine.columns.size());
		for (int row = 0; row < size; ++row) {
			const int coarse_row = aggregate[row];
			const auto begin = coarse.columns.begin()+coarse.row_offsets[coarse_row];
			const auto end = coarse.columns.begin()+coarse.row_offsets[coarse_row+1];
			for (int entry = fine.row_offsets[row]; entry < fine.row_offsets[row+1]; ++entry)
				coarse_entry[entry] = static_cast<int>(std::lower_bound(begin, end, aggregate[fine.columns[entry]])
					-coarse.columns.begin());
		}
	}

	void BuildInverse(const DevicePatternView& pattern, const double* values, double* inverse)
	{
		flags_.Clear();
		BuildScalarJacobiKernel<<<(pattern.nodes+255)/256,256>>>(pattern, values, inverse, flags_.data());
		unsigned int count = 0;
		flags_.CopyToHost(&count, 1);
		if (count) throw std::runtime_error("aggregation AMG met a zero diagonal");
	}

	void Smooth(const DevicePatternView& pattern, const double* values, const double* inverse, const double* b,
		double* x, double* r, int sweeps)
	{
		const int blocks = (pattern.nodes+255)/256;
		for (int sweep = 0; sweep < sweeps; ++sweep) {
			ResidualKernel<<<blocks,256>>>(pattern, values, b, x, r);
			DampedJacobiKernel<<<blocks,256>>>(pattern.nodes, options_.omega, inverse, r, x);
		}
	}

	void Visit(std::size_t depth, const DevicePatternView& pattern, const double* values, const double* inverse,
		const double* b, double* x)
	{
		const int size = pattern.nodes, blocks = (size+255)/256;
		if (depth == levels_.size()) {
			// Coarsest level: x += A_c^-1 b through the dense inverse (row-major on
			// the device, so the transposed column-major GEMV applies it).
			const double one = 1.0;
			Check(cublasDgemv(blas_, CUBLAS_OP_T, size, size, &one, coarsest_inverse_->data(), size, b, 1, &one, x, 1),
				"cublasDgemv AMG coarsest");
			return;
		}
		Level& level = *levels_[depth];
		double* r = depth == 0 ? fine_r_.data() : level.r.data();
		Smooth(pattern, values, inverse, b, x, r, options_.pre_sweeps);
		ResidualKernel<<<blocks,256>>>(pattern, values, b, x, r);
		level.b.Clear();
		level.x.Clear();
		RestrictKernel<<<blocks,256>>>(size, level.aggregate.data(), r, level.b.data());
		Visit(depth+1, level.matrix->pattern(), level.matrix->values(), level.inverse.data(), level.b.data(),
			level.x.data());
		ProlongKernel<<<blocks,256>>>(size, level.aggregate.data(), options_.over_correction, level.x.data(), x);
		Smooth(pattern, values, inverse, b, x, r, options_.post_sweeps);
	}

	// Dense inverse of the coarsest operator by Gauss-Jordan elimination with
	// partial pivoting; a singular coarsest operator falls back to its Jacobi
	// inverse so the cycle stays defined.
	void InvertCoarsest(const DevicePatternView& pattern, const double* values)
	{
		const int n = coarsest_size_;
		Check(cudaMemcpy(coarsest_host_values_.data(), values, coarsest_host_values_.size()*sizeof(double),
			cudaMemcpyDeviceToHost), "copy AMG coarsest values");
		std::vector<double>& a = coarsest_dense_;
		std::fill(a.begin(), a.end(), 0.0);
		for (int row = 0; row < n; ++row)
			for (int entry = coarsest_host_pattern_.row_offsets[row]; entry < coarsest_host_pattern_.row_offsets[row+1]; ++entry)
				a[static_cast<std::size_t>(row)*n+coarsest_host_pattern_.columns[entry]] = coarsest_host_values_[static_cast<std::size_t>(entry)];
		std::vector<double> inverse(a.size(), 0.0);
		for (int i = 0; i < n; ++i) inverse[static_cast<std::size_t>(i)*n+i] = 1.0;
		bool singular = false;
		for (int column = 0; column < n && !singular; ++column) {
			int pivot = column;
			for (int row = column+1; row < n; ++row)
				if (std::abs(a[static_cast<std::size_t>(row)*n+column]) > std::abs(a[static_cast<std::size_t>(pivot)*n+column])) pivot = row;
			const double magnitude = std::abs(a[static_cast<std::size_t>(pivot)*n+column]);
			if (!(magnitude > 0.0) || !std::isfinite(magnitude)) { singular = true; break; }
			if (pivot != column)
				for (int k = 0; k < n; ++k) {
					std::swap(a[static_cast<std::size_t>(column)*n+k], a[static_cast<std::size_t>(pivot)*n+k]);
					std::swap(inverse[static_cast<std::size_t>(column)*n+k], inverse[static_cast<std::size_t>(pivot)*n+k]);
				}
			const double scale = 1.0/a[static_cast<std::size_t>(column)*n+column];
			for (int k = 0; k < n; ++k) {
				a[static_cast<std::size_t>(column)*n+k] *= scale;
				inverse[static_cast<std::size_t>(column)*n+k] *= scale;
			}
			for (int row = 0; row < n; ++row) {
				if (row == column) continue;
				const double factor = a[static_cast<std::size_t>(row)*n+column];
				if (factor == 0.0) continue;
				for (int k = 0; k < n; ++k) {
					a[static_cast<std::size_t>(row)*n+k] -= factor*a[static_cast<std::size_t>(column)*n+k];
					inverse[static_cast<std::size_t>(row)*n+k] -= factor*inverse[static_cast<std::size_t>(column)*n+k];
				}
			}
		}
		if (singular) {
			std::fill(inverse.begin(), inverse.end(), 0.0);
			for (int row = 0; row < n; ++row) {
				const double diagonal = coarsest_host_values_[static_cast<std::size_t>(coarsest_host_pattern_.diagonal[row])];
				inverse[static_cast<std::size_t>(row)*n+row] = diagonal != 0.0 ? 1.0/diagonal : 0.0;
			}
			++singular_coarsest_;
		}
		(void)pattern;
		coarsest_inverse_->CopyFromHost(inverse.data(), inverse.size());
	}

	const BlockMatrix<1>& fine_;
	BlockPattern host_;
	int block_;
	AggregationAmgOptions options_;
	int size_;
	bool built_ = false;
	int fallbacks_ = 0;
	std::vector<std::unique_ptr<Level>> levels_;
	int coarsest_size_ = 0, singular_coarsest_ = 0;
	BlockPattern coarsest_host_pattern_;
	std::vector<double> coarsest_host_values_, coarsest_dense_;
	std::unique_ptr<DeviceBuffer<double>> coarsest_inverse_;
	BlasHandle blas_;
	DeviceBuffer<double> fine_inverse_, fine_r_;
	DeviceBuffer<unsigned int> flags_;
};

} // namespace iga::cuda

#endif

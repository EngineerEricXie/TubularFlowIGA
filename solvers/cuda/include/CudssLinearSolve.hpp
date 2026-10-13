#ifndef IGA_CUDA_CUDSS_LINEAR_SOLVE_HPP
#define IGA_CUDA_CUDSS_LINEAR_SOLVE_HPP

#include "CudssSolver.hpp"
#include "NativeTetCsrMatrix.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace iga::cuda {

// A NativeTetLinearSolve that factorizes host CSR systems with cuDSS. The
// solid tangent and the harmonic mesh-motion stiffness keep their pattern
// across calls, so the device pattern, matrix, and cuDSS analysis are reused
// whenever the pattern matches the previous call; only values are uploaded.
class CudssLinearSolve {
public:
	std::vector<double> operator()(const NativeTetCsrMatrix& matrix, const std::vector<double>& right_hand_side)
	{
		if (matrix.row_offsets.size() != matrix.size+1 || right_hand_side.size() != matrix.size)
			throw std::invalid_argument("cuDSS linear solve received an inconsistent system");
		if (!cache_ || !cache_->SamePattern(matrix)) cache_.emplace(matrix);
		cache_->matrix.CopyValues(matrix.values);
		DeviceBuffer<double>& rhs = cache_->right_hand_side;
		rhs.CopyFromHost(right_hand_side.data(), right_hand_side.size());
		cache_->solver.Solve(rhs.data());
		std::vector<double> solution(matrix.size);
		rhs.CopyToHost(solution.data(), solution.size());
		return solution;
	}

private:
	// Device-side copy of one host CSR pattern.
	struct DeviceCsr {
		explicit DeviceCsr(const NativeTetCsrMatrix& host)
			: columns(host.columns.begin(), host.columns.end()), pattern(HostPattern(host)), values(pattern)
		{}

		static BlockPattern HostPattern(const NativeTetCsrMatrix& host)
		{
			BlockPattern pattern;
			pattern.row_offsets.reserve(host.row_offsets.size());
			for (const auto offset : host.row_offsets) pattern.row_offsets.push_back(static_cast<int>(offset));
			pattern.columns.reserve(host.columns.size());
			for (const auto column : host.columns) pattern.columns.push_back(static_cast<int>(column));
			pattern.diagonal.assign(host.size, -1);
			for (std::size_t row = 0; row < host.size; ++row)
				for (auto entry = host.row_offsets[row]; entry < host.row_offsets[row+1]; ++entry)
					if (host.columns[entry] == row) pattern.diagonal[row] = static_cast<int>(entry);
			for (const int diagonal : pattern.diagonal)
				if (diagonal < 0) throw std::runtime_error("cuDSS linear solve needs a diagonal entry in every row");
			return pattern;
		}

		void CopyValues(const std::vector<double>& host_values)
		{
			Check(cudaMemcpy(values.values(), host_values.data(), host_values.size()*sizeof(double),
				cudaMemcpyHostToDevice), "copy CSR values to cuDSS matrix");
		}

		std::vector<std::uint32_t> columns;
		DevicePattern pattern;
		BlockMatrix<1> values;
	};

	struct Cache {
		explicit Cache(const NativeTetCsrMatrix& host)
			: matrix(host), solver(matrix.values), right_hand_side(host.size), row_offsets(host.row_offsets) {}

		bool SamePattern(const NativeTetCsrMatrix& host) const
		{
			return host.row_offsets == row_offsets && host.columns == matrix.columns;
		}

		DeviceCsr matrix;
		CudssSolver solver;
		DeviceBuffer<double> right_hand_side;
		std::vector<std::size_t> row_offsets;
	};

	std::optional<Cache> cache_;
};

} // namespace iga::cuda

#endif

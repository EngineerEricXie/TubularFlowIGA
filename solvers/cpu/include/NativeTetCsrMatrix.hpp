#ifndef IGA_NATIVE_TET_CSR_MATRIX_HPP
#define IGA_NATIVE_TET_CSR_MATRIX_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

// Marks a constrained unknown in a reduced-index map.
constexpr std::size_t kNativeTetConstrainedDof = std::numeric_limits<std::size_t>::max();

// Square sparse matrix in CSR form (columns sorted within each row) whose
// pattern is the vertex coupling of a tetrahedral mesh. It holds only the
// reduced system: unknown `components*node+component` is row
// `reduced_index[...]`, or absent when that entry is kNativeTetConstrainedDof.
struct NativeTetCsrMatrix
{
	std::size_t size = 0;
	std::vector<std::size_t> row_offsets;
	std::vector<std::uint32_t> columns;
	std::vector<double> values;

	double& At(std::size_t row, std::size_t column)
	{
		const auto begin = columns.begin()+static_cast<std::ptrdiff_t>(row_offsets[row]);
		const auto end = columns.begin()+static_cast<std::ptrdiff_t>(row_offsets[row+1]);
		const auto found = std::lower_bound(begin, end, static_cast<std::uint32_t>(column));
		if (found == end || *found != column)
			throw std::logic_error("native tetrahedral CSR entry is outside the pattern");
		return values[static_cast<std::size_t>(found-columns.begin())];
	}

	void Clear() { std::fill(values.begin(), values.end(), 0.0); }
};

inline NativeTetCsrMatrix BuildNativeTetCsrPattern(const NativeTetMesh& mesh, std::size_t components,
	const std::vector<std::size_t>& reduced_index, std::size_t size)
{
	if (components == 0 || reduced_index.size() != components*mesh.points.size()
		|| size > std::numeric_limits<std::uint32_t>::max())
		throw std::invalid_argument("native tetrahedral CSR pattern input is invalid");
	std::vector<std::vector<std::uint32_t>> rows(size);
	for (const auto& cell : mesh.cells)
		for (const auto row_node : cell.nodes)
			for (std::size_t i = 0; i < components; ++i) {
				const auto row = reduced_index[components*row_node+i];
				if (row == kNativeTetConstrainedDof) continue;
				for (const auto column_node : cell.nodes)
					for (std::size_t k = 0; k < components; ++k) {
						const auto column = reduced_index[components*column_node+k];
						if (column != kNativeTetConstrainedDof)
							rows[row].push_back(static_cast<std::uint32_t>(column));
					}
			}
	NativeTetCsrMatrix matrix;
	matrix.size = size;
	matrix.row_offsets.assign(size+1, 0);
	for (std::size_t row = 0; row < size; ++row) {
		auto& entries = rows[row];
		std::sort(entries.begin(), entries.end());
		entries.erase(std::unique(entries.begin(), entries.end()), entries.end());
		matrix.row_offsets[row+1] = matrix.row_offsets[row]+entries.size();
	}
	matrix.columns.reserve(matrix.row_offsets.back());
	for (const auto& entries : rows) matrix.columns.insert(matrix.columns.end(), entries.begin(), entries.end());
	matrix.values.assign(matrix.columns.size(), 0.0);
	return matrix;
}

// Row-major dense copy for the small-problem fallback solvers.
inline std::vector<double> DenseNativeTetCsrMatrix(const NativeTetCsrMatrix& matrix)
{
	std::vector<double> dense(matrix.size*matrix.size, 0.0);
	for (std::size_t row = 0; row < matrix.size; ++row)
		for (std::size_t entry = matrix.row_offsets[row]; entry < matrix.row_offsets[row+1]; ++entry)
			dense[row*matrix.size+matrix.columns[entry]] = matrix.values[entry];
	return dense;
}

// Solves A x = b for a reduced CSR system. Solvers that are given one use it
// (for example PETSc or cuDSS); an empty function keeps their dense fallback.
using NativeTetLinearSolve = std::function<std::vector<double>(const NativeTetCsrMatrix&, const std::vector<double>&)>;

} // namespace iga

#endif

#ifndef IGA_CARTESIAN_BACKGROUND_HPP
#define IGA_CARTESIAN_BACKGROUND_HPP

#include "IgaDatabase.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace iga {

struct CubicCartesianGridSpec {
	std::array<double, 3> lower_m{};
	std::array<double, 3> upper_m{};
	std::array<std::uint32_t, 3> cells{};
};

class CubicCartesianBackground {
public:
	enum NeighborSlot : std::size_t { ZMinus = 0, YMinus = 1, XPlus = 2,
		YPlus = 3, XMinus = 4, ZPlus = 5 };
	struct CellDescriptor {
		std::array<std::uint32_t, 3> index{};
		std::array<double, 3> lower_m{};
		std::array<double, 3> upper_m{};
		std::array<std::uint64_t, 6> neighbor{{kNoNeighbor, kNoNeighbor, kNoNeighbor,
			kNoNeighbor, kNoNeighbor, kNoNeighbor}};
	};
	static constexpr std::uint64_t kNoNeighbor = std::numeric_limits<std::uint64_t>::max();

	explicit CubicCartesianBackground(CubicCartesianGridSpec spec) : spec_(spec)
	{
		for (int axis = 0; axis < 3; ++axis) {
			if (!std::isfinite(spec_.lower_m[axis]) || !std::isfinite(spec_.upper_m[axis])
				|| !(spec_.upper_m[axis] > spec_.lower_m[axis]) || !spec_.cells[axis])
				throw std::invalid_argument("Cartesian background bounds and cell counts must be finite and positive");
			h_[axis] = (spec_.upper_m[axis]-spec_.lower_m[axis])/spec_.cells[axis];
			if (!std::isfinite(h_[axis]) || !(h_[axis] > 0.0))
				throw std::invalid_argument("Cartesian background spacing must be finite and positive");
		}
		std::array<std::uint64_t, 3> basis{};
		for (int axis = 0; axis < 3; ++axis) basis[axis] = static_cast<std::uint64_t>(spec_.cells[axis])+3;
		const auto count = CheckedProduct(spec_.cells[0], spec_.cells[1], spec_.cells[2]);
		const auto nodes = CheckedProduct(basis[0], basis[1], basis[2]);
		if (!count || !nodes || count > std::numeric_limits<std::uint64_t>::max()/64
			|| nodes > static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max()))
			throw std::overflow_error("Cartesian background count exceeds supported id or int32 connectivity range");
		element_count_ = count;
		node_count_ = nodes;
		for (int axis = 0; axis < 3; ++axis) BuildAxis(axis);
	}

	std::uint64_t ElementCount() const noexcept { return element_count_; }
	std::uint64_t NodeCount() const noexcept { return node_count_; }
	const CubicCartesianGridSpec& Spec() const noexcept { return spec_; }
	std::array<double, 3> Greville(std::uint32_t i, std::uint32_t j, std::uint32_t k) const
	{
		const std::array<std::uint32_t, 3> index{{i, j, k}};
		std::array<double, 3> result{};
		for (int axis = 0; axis < 3; ++axis) {
			const auto basis_count = static_cast<std::uint64_t>(spec_.cells[axis])+3;
			if (static_cast<std::uint64_t>(index[axis]) >= basis_count)
				throw std::out_of_range("Greville index is out of range");
			result[axis] = spec_.lower_m[axis] + h_[axis]
				*(knots_[axis][index[axis]+1]+knots_[axis][index[axis]+2]+knots_[axis][index[axis]+3])/3.0;
		}
		return result;
	}

	CellDescriptor Cell(std::uint64_t id) const
	{
		if (id >= element_count_) throw std::out_of_range("Cartesian background element id is out of range");
		CellDescriptor result;
		result.index = {{static_cast<std::uint32_t>(id%spec_.cells[0]),
			static_cast<std::uint32_t>((id/spec_.cells[0])%spec_.cells[1]),
			static_cast<std::uint32_t>(id/(static_cast<std::uint64_t>(spec_.cells[0])*spec_.cells[1]))}};
		for (int axis = 0; axis < 3; ++axis) {
			result.lower_m[axis] = spec_.lower_m[axis]+h_[axis]*result.index[axis];
			result.upper_m[axis] = result.lower_m[axis]+h_[axis];
		}
		if (result.index[2]) result.neighbor[ZMinus] = id-static_cast<std::uint64_t>(spec_.cells[0])*spec_.cells[1];
		if (result.index[1]) result.neighbor[YMinus] = id-spec_.cells[0];
		if (result.index[0]+1 < spec_.cells[0]) result.neighbor[XPlus] = id+1;
		if (result.index[1]+1 < spec_.cells[1]) result.neighbor[YPlus] = id+spec_.cells[0];
		if (result.index[0]) result.neighbor[XMinus] = id-1;
		if (result.index[2]+1 < spec_.cells[2]) result.neighbor[ZPlus] = id+static_cast<std::uint64_t>(spec_.cells[0])*spec_.cells[1];
		return result;
	}

	Element MaterializeElement(std::uint64_t id) const
	{
		if (id >= element_count_) throw std::out_of_range("Cartesian background element id is out of range");
		const auto ex = static_cast<std::uint32_t>(id%spec_.cells[0]);
		const auto ey = static_cast<std::uint32_t>((id/spec_.cells[0])%spec_.cells[1]);
		const auto ez = static_cast<std::uint32_t>(id/(static_cast<std::uint64_t>(spec_.cells[0])*spec_.cells[1]));
		Element element;
		element.id = id;
		element.connectivity.resize(64);
		element.extraction.resize(64);
		std::size_t row = 0;
		for (std::uint32_t k = 0; k < 4; ++k)
			for (std::uint32_t j = 0; j < 4; ++j)
				for (std::uint32_t i = 0; i < 4; ++i, ++row) {
					element.connectivity[row] = static_cast<std::int32_t>((ex+i)
						+ (spec_.cells[0]+3)*((ey+j) + (spec_.cells[1]+3)*(ez+k)));
					for (std::uint32_t c = 0; c < 4; ++c)
						for (std::uint32_t b = 0; b < 4; ++b)
							for (std::uint32_t a = 0; a < 4; ++a)
								element.extraction[row][a+4*(b+4*c)] = extraction_[0][ex][i][a]
									*extraction_[1][ey][j][b]*extraction_[2][ez][k][c];
				}
		std::size_t point = 0;
		for (int c = 0; c < 4; ++c)
			for (int b = 0; b < 4; ++b)
				for (int a = 0; a < 4; ++a, ++point)
					element.bezier_points[point] = {{spec_.lower_m[0]+h_[0]*(ex+a/3.0),
						spec_.lower_m[1]+h_[1]*(ey+b/3.0), spec_.lower_m[2]+h_[2]*(ez+c/3.0)}};
		return element;
	}

private:
	static std::uint64_t CheckedProduct(std::uint64_t first, std::uint64_t second, std::uint64_t third)
	{
		if (first > std::numeric_limits<std::uint64_t>::max()/second) throw std::overflow_error("Cartesian background product overflows");
		const auto product = first*second;
		if (product > std::numeric_limits<std::uint64_t>::max()/third) throw std::overflow_error("Cartesian background product overflows");
		return product*third;
	}
	void BuildAxis(int axis)
	{
		const auto cells = spec_.cells[axis];
		if (cells <= 5) {
			extraction_[axis] = ExactDenseAxisExtraction(cells);
			BuildKnots(axis);
			return;
		}
		const auto template_extraction = ExactDenseAxisExtraction(6);
		extraction_[axis].resize(cells);
		for (std::uint32_t cell = 0; cell < cells; ++cell) {
			const auto template_cell = cell < 3 ? cell
				: (cell+3 >= cells ? 6-(cells-cell) : 3);
			extraction_[axis][cell] = template_extraction[template_cell];
		}
		BuildKnots(axis);
	}

	void BuildKnots(int axis)
	{
		const auto cells = spec_.cells[axis];
		std::vector<double> knots(4, 0.0);
		for (std::uint32_t i = 1; i < cells; ++i) knots.push_back(static_cast<double>(i));
		knots.insert(knots.end(), 4, static_cast<double>(cells));
		knots_[axis] = knots;
	}

	static std::vector<std::array<std::array<double, 4>, 4>> ExactDenseAxisExtraction(
		std::uint32_t cells)
	{
		std::vector<double> knots(4, 0.0);
		for (std::uint32_t i = 1; i < cells; ++i) knots.push_back(static_cast<double>(i));
		knots.insert(knots.end(), 4, static_cast<double>(cells));
		std::vector<std::vector<double>> controls(cells+3, std::vector<double>(cells+3, 0.0));
		for (std::uint32_t i = 0; i < cells+3; ++i) controls[i][i] = 1.0;
		for (std::uint32_t value = 1; value < cells; ++value)
			for (int repeat = 0; repeat < 2; ++repeat) InsertKnot(knots, controls, value);
		std::vector<std::array<std::array<double, 4>, 4>> result(cells);
		for (std::uint32_t cell = 0; cell < cells; ++cell)
			for (std::uint32_t i = 0; i < 4; ++i)
				for (std::uint32_t a = 0; a < 4; ++a)
					result[cell][i][a] = controls[3*cell+a][cell+i];
		return result;
	}

	static void InsertKnot(std::vector<double>& knots, std::vector<std::vector<double>>& control,
		double value)
	{
		const int degree = 3;
		const int n = static_cast<int>(control.size())-1;
		int span = degree;
		while (span+1 < static_cast<int>(knots.size()) && knots[span+1] <= value) ++span;
		int multiplicity = 0;
		for (double knot : knots) if (knot == value) ++multiplicity;
		std::vector<std::vector<double>> next(control.size()+1, std::vector<double>(control.front().size()));
		for (int i = 0; i <= span-degree; ++i) next[i] = control[i];
		for (int i = span-multiplicity; i <= n; ++i) next[i+1] = control[i];
		for (int i = span-degree+1; i <= span-multiplicity; ++i) {
			const auto alpha = (value-knots[i])/(knots[i+degree]-knots[i]);
			for (std::size_t j = 0; j < next[i].size(); ++j)
				next[i][j] = alpha*control[i][j] + (1.0-alpha)*control[i-1][j];
		}
		knots.insert(knots.begin()+span+1, value);
		control = std::move(next);
	}

	CubicCartesianGridSpec spec_;
	std::array<double, 3> h_{};
	std::array<std::vector<double>, 3> knots_;
	std::array<std::vector<std::array<std::array<double, 4>, 4>>, 3> extraction_;
	std::uint64_t element_count_ = 0;
	std::uint64_t node_count_ = 0;
};

} // namespace iga

#endif

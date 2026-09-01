#include "CartesianBackground.hpp"
#include "ElementGeometry.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

double Cox(const std::vector<double>& knots, int i, int degree, double value)
{
	if (!degree) return (knots[i] <= value && value < knots[i+1]) ? 1.0 : 0.0;
	double result = 0.0;
	if (knots[i+degree] > knots[i]) result += (value-knots[i])/(knots[i+degree]-knots[i])
		*Cox(knots, i, degree-1, value);
	if (knots[i+degree+1] > knots[i+1]) result += (knots[i+degree+1]-value)
		/(knots[i+degree+1]-knots[i+1])*Cox(knots, i+1, degree-1, value);
	return result;
}

double CoxFirst(const std::vector<double>& knots, int i, int degree, double value)
{
	return degree*((knots[i+degree] > knots[i] ? Cox(knots, i, degree-1, value)/(knots[i+degree]-knots[i]) : 0.0)
		-(knots[i+degree+1] > knots[i+1] ? Cox(knots, i+1, degree-1, value)/(knots[i+degree+1]-knots[i+1]) : 0.0));
}

double CoxSecond(const std::vector<double>& knots, int i, int degree, double value)
{
	return degree*((knots[i+degree] > knots[i] ? CoxFirst(knots, i, degree-1, value)/(knots[i+degree]-knots[i]) : 0.0)
		-(knots[i+degree+1] > knots[i+1] ? CoxFirst(knots, i+1, degree-1, value)/(knots[i+degree+1]-knots[i+1]) : 0.0));
}

void Bernstein(double x, std::array<double, 4>& value, std::array<double, 4>& first,
	std::array<double, 4>& second)
{
	value = {{std::pow(1-x, 3), 3*std::pow(1-x, 2)*x, 3*(1-x)*x*x, x*x*x}};
	first = {{-3*std::pow(1-x, 2), 3-12*x+9*x*x, 6*x-9*x*x, 3*x*x}};
	second = {{6*(1-x), -12+18*x, 6-18*x, 6*x}};
}

void ExtractedOneD(const iga::Element& element, int row, double x,
	double& value, double& first, double& second)
{
	std::array<double, 4> b{}, d{}, dd{};
	Bernstein(x, b, d, dd);
	value = first = second = 0.0;
	for (int a = 0; a < 4; ++a) {
		const auto coefficient = element.extraction[row][a];
		value += coefficient*b[a]; first += coefficient*d[a]; second += coefficient*dd[a];
	}
}

void RequireAxisExtraction(int axis)
{
	std::array<std::uint32_t, 3> cells{{1, 1, 1}};
	cells[axis] = 7;
	iga::CubicCartesianBackground grid({{{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, cells});
	std::vector<double> knots(4, 0.0);
	for (int i = 1; i < 7; ++i) knots.push_back(i);
	knots.insert(knots.end(), 4, 7.0);
	for (std::uint32_t cell = 0; cell < 7; ++cell) {
		const auto element = grid.MaterializeElement(cell);
		for (const auto xi : std::array<double, 3>{{0.19, 0.47, 0.81}})
			for (int local = 0; local < 4; ++local) {
				const int row = axis == 0 ? local : axis == 1 ? 4*local : 16*local;
				std::array<double, 4> b{}, d{}, dd{};
				Bernstein(xi, b, d, dd);
				double value = 0.0, first = 0.0, second = 0.0;
				for (int a = 0; a < 4; ++a) {
					const int column = axis == 0 ? a : axis == 1 ? 4*a : 16*a;
					const auto q = element.extraction[row][column];
					value += q*b[a]; first += q*d[a]; second += q*dd[a];
				}
				assert(std::abs(value-Cox(knots, cell+local, 3, cell+xi)) < 2e-13);
				assert(std::abs(first-CoxFirst(knots, cell+local, 3, cell+xi)) < 2e-12);
				assert(std::abs(second-CoxSecond(knots, cell+local, 3, cell+xi)) < 2e-11);
			}
	}
}

int main()
{
	bool oversized_rejected = false;
	try { iga::CubicCartesianBackground oversized({{{0.0,0.0,0.0}},{{1.0,1.0,1.0}},{{std::numeric_limits<std::uint32_t>::max(),1,1}}}); }
	catch (const std::overflow_error&) { oversized_rejected = true; }
	assert(oversized_rejected);
	iga::CubicCartesianGridSpec one{{{0.0, 0.0, 0.0}}, {{1.0, 2.0, 3.0}}, {{1, 1, 1}}};
	iga::CubicCartesianBackground grid(one);
	assert(grid.ElementCount() == 1 && grid.NodeCount() == 64);
	const auto element = grid.MaterializeElement(0);
	for (std::size_t i = 0; i < 64; ++i) {
		assert(element.connectivity[i] == static_cast<std::int32_t>(i));
		assert(element.boundary_labels[i/16] == -1);
		for (std::size_t j = 0; j < 64; ++j)
			assert(std::abs(element.extraction[i][j]-(i == j ? 1.0 : 0.0)) < 1e-14);
	}
	for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
		const auto point = element.bezier_points[i+4*(j+4*k)];
		assert(std::abs(point[0]-i/3.0) < 1e-14 && std::abs(point[1]-2.0*j/3.0) < 1e-14
			&& std::abs(point[2]-k) < 1e-14);
	}
	const auto geometry = iga::EvaluateElementGeometry(element, {{0.31, 0.52, 0.73}});
	assert(std::abs(geometry.raw_determinant-6.0) < 2e-13);
	RequireAxisExtraction(0); RequireAxisExtraction(1); RequireAxisExtraction(2);
	iga::CubicCartesianBackground topology({{{0.0, 0.0, 0.0}}, {{3.0, 3.0, 3.0}}, {{3, 3, 3}}});
	const auto center = topology.Cell(13);
	assert((center.lower_m == std::array<double, 3>{{1.0, 1.0, 1.0}}
		&& center.upper_m == std::array<double, 3>{{2.0, 2.0, 2.0}}));
	assert((center.neighbor == std::array<std::uint64_t, 6>{{4, 10, 14, 16, 12, 22}}));
	const auto corner = topology.Cell(0);
	assert(corner.neighbor[iga::CubicCartesianBackground::ZMinus] == iga::CubicCartesianBackground::kNoNeighbor
		&& corner.neighbor[iga::CubicCartesianBackground::YMinus] == iga::CubicCartesianBackground::kNoNeighbor
		&& corner.neighbor[iga::CubicCartesianBackground::XMinus] == iga::CubicCartesianBackground::kNoNeighbor);
	iga::CubicCartesianGridSpec multi{{{-1.0, 2.0, 4.0}}, {{3.0, 5.0, 10.0}}, {{2, 3, 1}}};
	iga::CubicCartesianBackground background(multi);
	assert(background.ElementCount() == 6 && background.NodeCount() == 120);
	const auto first = background.MaterializeElement(0);
	const auto adjacent = background.MaterializeElement(1);
	int shared = 0;
	for (auto left : first.connectivity)
		for (auto right : adjacent.connectivity)
			if (left == right) ++shared;
	assert(shared == 48);
	const auto cell = background.Cell(1);
	assert((cell.index == std::array<std::uint32_t, 3>{{1, 0, 0}}));
	assert(cell.neighbor[4] == 0 && cell.neighbor[2] == iga::CubicCartesianBackground::kNoNeighbor);
	const std::vector<double> knots_x{{0, 0, 0, 0, 1, 2, 2, 2, 2}};
	const std::vector<double> knots_y{{0, 0, 0, 0, 1, 2, 3, 3, 3, 3}};
	const std::vector<double> knots_z{{0, 0, 0, 0, 1, 1, 1, 1}};
	for (const auto& point : std::array<std::array<double, 3>, 2>{{{{0.23, 0.44, 0.61}}, {{0.71, 0.12, 0.83}}}}) {
		std::array<double, 4> bx{}, by{}, bz{}, dx{}, dy{}, dz{}, ddx{}, ddy{}, ddz{};
		Bernstein(point[0], bx, dx, ddx); Bernstein(point[1], by, dy, ddy); Bernstein(point[2], bz, dz, ddz);
		double sum = 0.0, dsum = 0.0, ddsum = 0.0;
		for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
			const auto row = i+4*(j+4*k);
			double value = 0.0, first = 0.0, second = 0.0;
			for (int c = 0; c < 4; ++c) for (int b = 0; b < 4; ++b) for (int a = 0; a < 4; ++a) {
				const auto coefficient = adjacent.extraction[row][a+4*(b+4*c)];
				value += coefficient*bx[a]*by[b]*bz[c];
				first += coefficient*dx[a]*by[b]*bz[c];
				second += coefficient*ddx[a]*by[b]*bz[c];
			}
			const auto oracle = Cox(knots_x, i+1, 3, 1+point[0])*Cox(knots_y, j, 3, point[1])*Cox(knots_z, k, 3, point[2]);
			const auto oracle_first = CoxFirst(knots_x, i+1, 3, 1+point[0])
				*Cox(knots_y, j, 3, point[1])*Cox(knots_z, k, 3, point[2]);
			const auto oracle_second = CoxSecond(knots_x, i+1, 3, 1+point[0])
				*Cox(knots_y, j, 3, point[1])*Cox(knots_z, k, 3, point[2]);
			assert(value >= -2e-13 && std::abs(value-oracle) < 2e-13);
			assert(std::abs(first-oracle_first) < 2e-12);
			assert(std::abs(second-oracle_second) < 2e-11);
			sum += value; dsum += first; ddsum += second;
		}
		assert(std::abs(sum-1.0) < 2e-13 && std::abs(dsum) < 2e-13 && std::abs(ddsum) < 2e-12);
	}
	assert(std::abs(background.Greville(0, 0, 0)[0]+1.0) < 2e-13);
	assert(std::abs(background.Greville(4, 5, 3)[1]-5.0) < 2e-13);
	for (int axis = 0; axis < 3; ++axis) {
		bool out_of_range = false;
		try { background.Greville(axis == 0 ? 5U : 0U, axis == 1 ? 6U : 0U,
			axis == 2 ? 4U : 0U); } catch (const std::out_of_range&) { out_of_range = true; }
		assert(out_of_range);
		out_of_range = false;
		try { background.Greville(axis == 0 ? std::numeric_limits<std::uint32_t>::max() : 0U,
			axis == 1 ? std::numeric_limits<std::uint32_t>::max() : 0U,
			axis == 2 ? std::numeric_limits<std::uint32_t>::max() : 0U); }
		catch (const std::out_of_range&) { out_of_range = true; }
		assert(out_of_range);
	}
	constexpr std::array<double, 4> weights{{0.17392742260172695, 0.32607257707327305,
		0.32607257707327305, 0.17392742260172695}};
	double volume = 0.0;
	for (std::uint64_t id = 0; id < background.ElementCount(); ++id) {
		const auto item = background.MaterializeElement(id);
		for (double wz : weights) for (double wy : weights) for (double wx : weights)
			volume += wx*wy*wz*iga::EvaluateElementGeometry(item, {{0.5, 0.5, 0.5}}).raw_determinant;
	}
	assert(std::abs(volume-72.0) < 1e-6);
	for (const auto cells : std::array<std::uint32_t, 4>{{2, 3, 5, 7}}) {
		iga::CubicCartesianBackground axis_grid({{{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}},
			{{cells, 1, 1}}});
		for (std::uint32_t interface = 0; interface+1 < cells; ++interface) {
			const auto left = axis_grid.MaterializeElement(interface);
			const auto right = axis_grid.MaterializeElement(interface+1);
			for (std::uint32_t global = interface > 2 ? interface-2 : 0;
				global <= std::min(cells+2, interface+4); ++global) {
				double left_value = 0.0, left_first = 0.0, left_second = 0.0;
				double right_value = 0.0, right_first = 0.0, right_second = 0.0;
				for (int row = 0; row < 4; ++row) {
					if (left.connectivity[row] == static_cast<std::int32_t>(global))
						ExtractedOneD(left, row, 1.0, left_value, left_first, left_second);
					if (right.connectivity[row] == static_cast<std::int32_t>(global))
						ExtractedOneD(right, row, 0.0, right_value, right_first, right_second);
				}
				assert(std::abs(left_value-right_value) < 2e-13);
				assert(std::abs(left_first-right_first) < 2e-12);
				assert(std::abs(left_second-right_second) < 2e-11);
			}
		}
	}
	for (const auto id : std::array<std::uint64_t, 2>{{0, 5}}) {
		const auto item = background.MaterializeElement(id);
		const auto descriptor = background.Cell(id);
		for (const auto& point : std::array<std::array<double, 3>, 2>{{{{0.21, 0.43, 0.67}}, {{0.74, 0.28, 0.39}}}}) {
			std::array<double, 4> bx{}, by{}, bz{}, dx{}, dy{}, dz{}, ddx{}, ddy{}, ddz{};
			Bernstein(point[0], bx, dx, ddx); Bernstein(point[1], by, dy, ddy); Bernstein(point[2], bz, dz, ddz);
			std::array<double, 3> mapped{};
			std::array<std::array<double, 3>, 3> derivative{};
			std::array<std::array<std::array<double, 3>, 3>, 3> second{};
			for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
				const auto row = i+4*(j+4*k);
				const auto node = static_cast<std::uint32_t>(item.connectivity[row]);
				const auto gi = node%5, gj = (node/5)%6, gk = node/30;
				const auto greville = background.Greville(gi, gj, gk);
				double n = 0.0;
				std::array<double, 3> d{};
				std::array<std::array<double, 3>, 3> dd{};
				for (int c = 0; c < 4; ++c) for (int b = 0; b < 4; ++b) for (int a = 0; a < 4; ++a) {
					const auto q = item.extraction[row][a+4*(b+4*c)];
					n += q*bx[a]*by[b]*bz[c];
					d[0] += q*dx[a]*by[b]*bz[c]; d[1] += q*bx[a]*dy[b]*bz[c]; d[2] += q*bx[a]*by[b]*dz[c];
					dd[0][0] += q*ddx[a]*by[b]*bz[c]; dd[1][1] += q*bx[a]*ddy[b]*bz[c]; dd[2][2] += q*bx[a]*by[b]*ddz[c];
					dd[0][1] += q*dx[a]*dy[b]*bz[c]; dd[0][2] += q*dx[a]*by[b]*dz[c]; dd[1][2] += q*bx[a]*dy[b]*dz[c];
				}
				for (int physical = 0; physical < 3; ++physical) {
					mapped[physical] += n*greville[physical];
					for (int direction = 0; direction < 3; ++direction) derivative[direction][physical] += d[direction]*greville[physical];
					for (int first_direction = 0; first_direction < 3; ++first_direction)
						for (int second_direction = first_direction; second_direction < 3; ++second_direction)
							second[first_direction][second_direction][physical] += dd[first_direction][second_direction]*greville[physical];
				}
			}
			const std::array<double, 3> spacing{{2.0, 1.0, 6.0}};
			for (int physical = 0; physical < 3; ++physical) {
				assert(std::abs(mapped[physical]-(descriptor.lower_m[physical]+spacing[physical]*point[physical])) < 2e-12);
				for (int direction = 0; direction < 3; ++direction)
					assert(std::abs(derivative[direction][physical]-(direction == physical ? spacing[physical] : 0.0)) < 2e-12);
				for (int first_direction = 0; first_direction < 3; ++first_direction)
					for (int second_direction = first_direction; second_direction < 3; ++second_direction)
						assert(std::abs(second[first_direction][second_direction][physical]) < 2e-11);
			}
		}
	}
	assert(background.MaterializeElement(5).id == 5);
	const auto first_permutation = background.MaterializeElement(4);
	const auto ignored = background.MaterializeElement(2);
	(void)ignored;
	const auto second_permutation = background.MaterializeElement(4);
	assert(first_permutation.connectivity == second_permutation.connectivity
		&& first_permutation.extraction == second_permutation.extraction
		&& first_permutation.bezier_points == second_permutation.bezier_points);
	for (const auto& item : std::array<iga::Element, 2>{{background.MaterializeElement(0), background.MaterializeElement(5)}})
		for (const auto label : item.boundary_labels) assert(label == -1);
	bool rejected = false;
	try { iga::CubicCartesianBackground invalid({{{1.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}}, {{1, 1, 1}}}); }
	catch (const std::invalid_argument&) { rejected = true; }
	assert(rejected);
	for (int axis = 0; axis < 3; ++axis) {
		auto invalid = multi;
		invalid.cells[axis] = 0;
		rejected = false;
		try { iga::CubicCartesianBackground value(invalid); } catch (const std::invalid_argument&) { rejected = true; }
		assert(rejected);
		invalid = multi; invalid.lower_m[axis] = std::numeric_limits<double>::quiet_NaN();
		rejected = false;
		try { iga::CubicCartesianBackground value(invalid); } catch (const std::invalid_argument&) { rejected = true; }
		assert(rejected);
		invalid = multi; invalid.upper_m[axis] = std::numeric_limits<double>::infinity();
		rejected = false;
		try { iga::CubicCartesianBackground value(invalid); } catch (const std::invalid_argument&) { rejected = true; }
		assert(rejected);
		invalid = multi; invalid.upper_m[axis] = invalid.lower_m[axis];
		rejected = false;
		try { iga::CubicCartesianBackground value(invalid); } catch (const std::invalid_argument&) { rejected = true; }
		assert(rejected);
		invalid = multi; invalid.upper_m[axis] = invalid.lower_m[axis]-1.0;
		rejected = false;
		try { iga::CubicCartesianBackground value(invalid); } catch (const std::invalid_argument&) { rejected = true; }
		assert(rejected);
	}
	for (const auto id : std::array<std::uint64_t, 2>{{background.ElementCount(), std::numeric_limits<std::uint64_t>::max()}}) {
		rejected = false;
		try { background.MaterializeElement(id); } catch (const std::out_of_range&) { rejected = true; }
		assert(rejected);
	}
	for (const auto id : std::array<std::uint64_t, 2>{{background.ElementCount(), std::numeric_limits<std::uint64_t>::max()}}) {
		rejected = false;
		try { background.Cell(id); } catch (const std::out_of_range&) { rejected = true; }
		assert(rejected);
	}
	rejected = false;
	try { iga::CubicCartesianBackground overflow({{{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}},
		{{std::numeric_limits<std::uint32_t>::max(), 1, 1}}}); }
	catch (const std::overflow_error&) { rejected = true; }
	assert(rejected);
	iga::CubicCartesianBackground large({{{0.0, 0.0, 0.0}}, {{1.0, 1.0, 1.0}},
		{{4096, 1, 1}}});
	assert(large.ElementCount() == 4096 && large.NodeCount() == 4099*4*4);
	assert(large.MaterializeElement(2048).connectivity.size() == 64);
	std::cout << "Cartesian cubic background tests passed\n";
}

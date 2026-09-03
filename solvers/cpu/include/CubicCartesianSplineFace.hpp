#ifndef IGA_CUBIC_CARTESIAN_SPLINE_FACE_HPP
#define IGA_CUBIC_CARTESIAN_SPLINE_FACE_HPP

// Shared, deterministic trace data for two adjacent cubic Cartesian elements.
#include "CartesianDomainClassification.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace iga {

struct CubicCartesianSplineFaceQuadraturePoint {
	double tangential_0 = 0.0, tangential_1 = 0.0, weight = 0.0;
};

inline const std::array<CubicCartesianSplineFaceQuadraturePoint, 16>& CubicCartesianSplineFaceTangentialQuadrature()
{
	static const std::array<CubicCartesianSplineFaceQuadraturePoint, 16> points = [] {
		constexpr std::array<double, 4> x{{-.8611363115940526, -.3399810435848563,
			.3399810435848563, .8611363115940526}};
		constexpr std::array<double, 4> w{{.3478548451374538, .6521451548625461,
			.6521451548625461, .3478548451374538}};
		std::array<CubicCartesianSplineFaceQuadraturePoint, 16> result{};
		std::size_t k = 0;
		for (std::size_t a = 0; a < 4; ++a)
			for (std::size_t b = 0; b < 4; ++b)
				result[k++] = {(x[a]+1.0)/2.0, (x[b]+1.0)/2.0, w[a]*w[b]};
		return result;
	}();
	return points;
}

inline std::vector<std::int32_t> CubicCartesianSplineFaceConnectivity(
	const CartesianDomainClassification& domain, std::uint64_t minus_cell, std::uint64_t plus_cell)
{
	auto result = domain.Background().MaterializeElement(minus_cell).connectivity;
	const auto plus = domain.Background().MaterializeElement(plus_cell);
	result.insert(result.end(), plus.connectivity.begin(), plus.connectivity.end());
	std::sort(result.begin(), result.end());
	result.erase(std::unique(result.begin(), result.end()), result.end());
	if (result.size() > 80) throw std::runtime_error("ghost penalty adjacent cubic face union exceeds 80 nodes");
	return result;
}

inline std::array<double, 4> CubicCartesianSplineBezierDerivative(std::uint32_t order, double t)
{
	if (order > 3 || !std::isfinite(t)) throw std::invalid_argument("cubic spline derivative arguments are invalid");
	const double u = 1.0-t;
	if (!order) return {{u*u*u, 3.0*u*u*t, 3.0*u*t*t, t*t*t}};
	if (order == 1) return {{-3.0*u*u, 3.0*u*u-6.0*u*t, 6.0*u*t-3.0*t*t, 3.0*t*t}};
	if (order == 2) return {{6.0*u, -12.0+18.0*t, 6.0-18.0*t, 6.0*t}};
	return {{-6.0, 18.0, -18.0, 6.0}};
}

inline std::vector<double> CubicCartesianSplineFaceNormalDerivativeJump(const Element& minus,
	const Element& plus, std::uint8_t axis, std::uint32_t order, double tangential_0,
	double tangential_1, double h_normal_m, const std::vector<std::int32_t>& nodes)
{
	if (axis > 2 || order > 3 || !std::isfinite(tangential_0) || !std::isfinite(tangential_1)
		|| !std::isfinite(h_normal_m) || !(h_normal_m > 0.0))
		throw std::invalid_argument("cubic Cartesian spline face jump arguments are invalid");
	const auto trace = [&](const Element& element, bool lower) {
		std::vector<double> result(nodes.size());
		const auto normal = CubicCartesianSplineBezierDerivative(order, lower ? 1.0 : 0.0);
		const auto b0 = CubicCartesianSplineBezierDerivative(0, tangential_0);
		const auto b1 = CubicCartesianSplineBezierDerivative(0, tangential_1);
		const double scale = std::pow(h_normal_m, -static_cast<int>(order));
		for (std::size_t row = 0; row < element.connectivity.size(); ++row) {
			double value = 0.0;
			for (int c = 0; c < 4; ++c) for (int b = 0; b < 4; ++b) for (int a = 0; a < 4; ++a) {
				const int index = a+4*(b+4*c), normal_index = axis == 0 ? a : axis == 1 ? b : c;
				const int tangent0_index = axis == 0 ? b : a, tangent1_index = axis == 2 ? b : c;
				value += element.extraction[row][index]*normal[normal_index]*b0[tangent0_index]*b1[tangent1_index];
			}
			const auto position = static_cast<std::size_t>(std::lower_bound(nodes.begin(), nodes.end(), element.connectivity[row])-nodes.begin());
			if (position == nodes.size() || nodes[position] != element.connectivity[row])
				throw std::invalid_argument("cubic Cartesian spline face connectivity does not cover element node");
			result[position] += (lower ? 1.0 : -1.0)*scale*value;
		}
		return result;
	};
	auto result = trace(minus, true);
	const auto other = trace(plus, false);
	for (std::size_t i = 0; i < result.size(); ++i) result[i] += other[i];
	return result;
}

} // namespace iga

#endif

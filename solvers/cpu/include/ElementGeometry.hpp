#ifndef IGA_ELEMENT_GEOMETRY_HPP
#define IGA_ELEMENT_GEOMETRY_HPP

#include "IgaDatabase.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace iga {

struct ElementGeometryPoint {
	std::array<double, 3> physical{};
	std::array<std::array<double, 3>, 3> jacobian{};
	double raw_determinant = 0.0;
};

inline double ElementGeometryDeterminant(const std::array<std::array<double, 3>, 3>& matrix)
{
	return matrix[0][0]*matrix[1][1]*matrix[2][2]
		+ matrix[0][1]*matrix[1][2]*matrix[2][0]
		+ matrix[0][2]*matrix[1][0]*matrix[2][1]
		- matrix[0][2]*matrix[1][1]*matrix[2][0]
		- matrix[0][0]*matrix[1][2]*matrix[2][1]
		- matrix[0][1]*matrix[1][0]*matrix[2][2];
}

inline ElementGeometryPoint EvaluateElementGeometry(const Element& element,
	const std::array<double, 3>& parametric)
{
	const double u = parametric[0];
	const double v = parametric[1];
	const double w = parametric[2];
	const double b[3][4] = {
		{std::pow(1.0-u, 3), 3.0*std::pow(1.0-u, 2)*u, 3.0*(1.0-u)*u*u, u*u*u},
		{std::pow(1.0-v, 3), 3.0*std::pow(1.0-v, 2)*v, 3.0*(1.0-v)*v*v, v*v*v},
		{std::pow(1.0-w, 3), 3.0*std::pow(1.0-w, 2)*w, 3.0*(1.0-w)*w*w, w*w*w}
	};
	const double db[3][4] = {
		{-3.0*std::pow(1.0-u, 2), 3.0-12.0*u+9.0*u*u, 3.0*(2.0-3.0*u)*u, 3.0*u*u},
		{-3.0*std::pow(1.0-v, 2), 3.0-12.0*v+9.0*v*v, 3.0*(2.0-3.0*v)*v, 3.0*v*v},
		{-3.0*std::pow(1.0-w, 2), 3.0-12.0*w+9.0*w*w, 3.0*(2.0-3.0*w)*w, 3.0*w*w}
	};
	ElementGeometryPoint result;
	std::size_t point = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++point) {
				const double value = b[0][i]*b[1][j]*b[2][k];
				const double derivative[3] = {db[0][i]*b[1][j]*b[2][k],
					b[0][i]*db[1][j]*b[2][k], b[0][i]*b[1][j]*db[2][k]};
				for (int physical = 0; physical < 3; ++physical) {
					result.physical[physical] += element.bezier_points[point][physical]*value;
					for (int parameter = 0; parameter < 3; ++parameter)
						result.jacobian[physical][parameter] +=
							element.bezier_points[point][physical]*derivative[parameter];
				}
			}
	result.raw_determinant = ElementGeometryDeterminant(result.jacobian);
	return result;
}

inline std::array<std::array<double, 3>, 3> InvertElementJacobian(
	const std::array<std::array<double, 3>, 3>& matrix, double determinant)
{
	if (!std::isfinite(determinant) || std::abs(determinant) < 1.0e-14)
		throw std::runtime_error("singular element Jacobian");
	std::array<std::array<double, 3>, 3> inverse{};
	inverse[0][0] = (matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])/determinant;
	inverse[0][1] = (matrix[2][1]*matrix[0][2]-matrix[0][1]*matrix[2][2])/determinant;
	inverse[0][2] = (matrix[0][1]*matrix[1][2]-matrix[0][2]*matrix[1][1])/determinant;
	inverse[1][0] = (matrix[1][2]*matrix[2][0]-matrix[1][0]*matrix[2][2])/determinant;
	inverse[1][1] = (matrix[0][0]*matrix[2][2]-matrix[0][2]*matrix[2][0])/determinant;
	inverse[1][2] = (matrix[0][2]*matrix[1][0]-matrix[0][0]*matrix[1][2])/determinant;
	inverse[2][0] = (matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0])/determinant;
	inverse[2][1] = (matrix[0][1]*matrix[2][0]-matrix[0][0]*matrix[2][1])/determinant;
	inverse[2][2] = (matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0])/determinant;
	return inverse;
}

} // namespace iga

#endif

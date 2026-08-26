#ifndef IGA_BOUNDARY_FLOW_HPP
#define IGA_BOUNDARY_FLOW_HPP

#include "IgaDatabase.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace iga {

struct BoundaryBasisValues {
	std::vector<double> value;
	std::vector<std::array<double, 3>> gradient;
	std::array<double, 3> physical_coordinate{};
	std::array<std::array<double, 3>, 3> inverse_jacobian{};
	double determinant = 0.0;
};

inline BoundaryBasisValues EvaluateBoundaryBasis(
	const Element& element, double u, double v, double w)
{
	const double b[3][4] = {
		{std::pow(1.0-u,3), 3.0*std::pow(1.0-u,2)*u, 3.0*(1.0-u)*u*u, u*u*u},
		{std::pow(1.0-v,3), 3.0*std::pow(1.0-v,2)*v, 3.0*(1.0-v)*v*v, v*v*v},
		{std::pow(1.0-w,3), 3.0*std::pow(1.0-w,2)*w, 3.0*(1.0-w)*w*w, w*w*w}
	};
	const double db[3][4] = {
		{-3.0*std::pow(1.0-u,2), 3.0-12.0*u+9.0*u*u, 3.0*(2.0-3.0*u)*u, 3.0*u*u},
		{-3.0*std::pow(1.0-v,2), 3.0-12.0*v+9.0*v*v, 3.0*(2.0-3.0*v)*v, 3.0*v*v},
		{-3.0*std::pow(1.0-w,2), 3.0-12.0*w+9.0*w*w, 3.0*(2.0-3.0*w)*w, 3.0*w*w}
	};
	std::array<double, 64> bernstein{};
	std::array<std::array<double, 3>, 64> parametric{};
	std::size_t location = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++location) {
				bernstein[location] = b[0][i]*b[1][j]*b[2][k];
				parametric[location] = {db[0][i]*b[1][j]*b[2][k],
					b[0][i]*db[1][j]*b[2][k], b[0][i]*b[1][j]*db[2][k]};
			}
	double jacobian[3][3]{};
	for (std::size_t p = 0; p < 64; ++p)
		for (int physical = 0; physical < 3; ++physical)
			for (int parameter = 0; parameter < 3; ++parameter)
				jacobian[physical][parameter] +=
					element.bezier_points[p][physical]*parametric[p][parameter];
	const double determinant =
		jacobian[0][0]*jacobian[1][1]*jacobian[2][2]
		+jacobian[0][1]*jacobian[1][2]*jacobian[2][0]
		+jacobian[0][2]*jacobian[1][0]*jacobian[2][1]
		-jacobian[0][2]*jacobian[1][1]*jacobian[2][0]
		-jacobian[0][0]*jacobian[1][2]*jacobian[2][1]
		-jacobian[0][1]*jacobian[1][0]*jacobian[2][2];
	if (!std::isfinite(determinant) || !(determinant > 0.0))
		throw std::runtime_error("non-positive element Jacobian on boundary face");
	double inverse[3][3]{};
	inverse[0][0] = (jacobian[1][1]*jacobian[2][2]-jacobian[1][2]*jacobian[2][1])/determinant;
	inverse[0][1] = (jacobian[2][1]*jacobian[0][2]-jacobian[0][1]*jacobian[2][2])/determinant;
	inverse[0][2] = (jacobian[0][1]*jacobian[1][2]-jacobian[1][1]*jacobian[0][2])/determinant;
	inverse[1][0] = (jacobian[2][0]*jacobian[1][2]-jacobian[1][0]*jacobian[2][2])/determinant;
	inverse[1][1] = (jacobian[0][0]*jacobian[2][2]-jacobian[2][0]*jacobian[0][2])/determinant;
	inverse[1][2] = (jacobian[1][0]*jacobian[0][2]-jacobian[0][0]*jacobian[1][2])/determinant;
	inverse[2][0] = (jacobian[1][0]*jacobian[2][1]-jacobian[2][0]*jacobian[1][1])/determinant;
	inverse[2][1] = (jacobian[2][0]*jacobian[0][1]-jacobian[0][0]*jacobian[2][1])/determinant;
	inverse[2][2] = (jacobian[0][0]*jacobian[1][1]-jacobian[1][0]*jacobian[0][1])/determinant;
	BoundaryBasisValues result;
	result.value.assign(element.connectivity.size(), 0.0);
	result.gradient.assign(element.connectivity.size(), {0.0, 0.0, 0.0});
	result.inverse_jacobian = {{{inverse[0][0], inverse[0][1], inverse[0][2]},
		{inverse[1][0], inverse[1][1], inverse[1][2]},
		{inverse[2][0], inverse[2][1], inverse[2][2]}}};
	result.determinant = 0.125*determinant;
	for (std::size_t p = 0; p < 64; ++p)
		for (int component = 0; component < 3; ++component)
			result.physical_coordinate[component] +=
				element.bezier_points[p][component]*bernstein[p];
	for (std::size_t a = 0; a < element.extraction.size(); ++a)
		for (std::size_t p = 0; p < 64; ++p) {
			result.value[a] += element.extraction[a][p]*bernstein[p];
			for (int physical = 0; physical < 3; ++physical)
				for (int parameter = 0; parameter < 3; ++parameter)
					result.gradient[a][physical] += element.extraction[a][p]
						*parametric[p][parameter]*inverse[parameter][physical];
		}
	return result;
}

inline double IntegrateVolumeDivergence(const Element& element,
	const std::vector<std::array<double, 4>>& nodal_state)
{
	if (nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("volume-divergence state does not match element connectivity");
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	double integral = 0.0;
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				const auto basis = EvaluateBoundaryBasis(
					element, points[qx], points[qy], points[qz]);
				double divergence = 0.0;
				for (std::size_t a = 0; a < nodal_state.size(); ++a)
					for (int component = 0; component < 3; ++component)
						divergence += nodal_state[a][component]
							*basis.gradient[a][component];
				integral += weights[qx]*weights[qy]*weights[qz]
					*basis.determinant*divergence;
			}
	return integral;
}

inline double IntegrateBoundaryFlow(const Element& element, std::size_t face,
	const std::vector<std::array<double, 4>>& nodal_state)
{
	if (face >= 6) throw std::runtime_error("boundary face index must be below six");
	if (nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("boundary-flow state does not match element connectivity");
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
	constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
	constexpr double outward_sign[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
	double flow = 0.0;
	for (std::size_t qi = 0; qi < 4; ++qi)
		for (std::size_t qj = 0; qj < 4; ++qj) {
			std::array<double, 3> coordinate{};
			coordinate[fixed_axis[face]] = fixed_value[face];
			coordinate[varying_axes[face][0]] = points[qi];
			coordinate[varying_axes[face][1]] = points[qj];
			const auto basis = EvaluateBoundaryBasis(
				element, coordinate[0], coordinate[1], coordinate[2]);
			std::array<double, 3> velocity{};
			for (std::size_t a = 0; a < nodal_state.size(); ++a)
				for (int component = 0; component < 3; ++component)
					velocity[component] += basis.value[a]*nodal_state[a][component];
			double normal_velocity = 0.0;
			for (int component = 0; component < 3; ++component)
				normal_velocity += velocity[component]
					*basis.inverse_jacobian[fixed_axis[face]][component];
			flow += outward_sign[face]*weights[qi]*weights[qj]*2.0
				*basis.determinant*normal_velocity;
		}
	return flow;
}

inline std::array<double, 2> IntegrateBoundaryScalarAndArea(
	const Element& element, std::size_t face,
	const std::vector<std::array<double, 4>>& nodal_state)
{
	if (face >= 6) throw std::runtime_error("boundary face index must be below six");
	if (nodal_state.size() != element.connectivity.size())
		throw std::runtime_error("boundary scalar state does not match element connectivity");
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
	constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
	std::array<double, 2> result{};
	for (std::size_t qi = 0; qi < 4; ++qi)
		for (std::size_t qj = 0; qj < 4; ++qj) {
			std::array<double, 3> coordinate{};
			coordinate[fixed_axis[face]] = fixed_value[face];
			coordinate[varying_axes[face][0]] = points[qi];
			coordinate[varying_axes[face][1]] = points[qj];
			const auto basis = EvaluateBoundaryBasis(
				element, coordinate[0], coordinate[1], coordinate[2]);
			double scalar = 0.0;
			for (std::size_t a = 0; a < nodal_state.size(); ++a)
				scalar += basis.value[a]*nodal_state[a][3];
			const double weight = weights[qi]*weights[qj]*2.0*basis.determinant;
			result[0] += weight*scalar;
			result[1] += weight;
		}
	return result;
}

inline double IntegrateBoundarySpeciesFlux(const Element& element, std::size_t face,
	const std::vector<std::array<double, 4>>& nodal_flow,
	const std::vector<double>& nodal_species)
{
	if (face >= 6) throw std::runtime_error("boundary face index must be below six");
	if (nodal_flow.size() != element.connectivity.size()
		|| nodal_species.size() != element.connectivity.size())
		throw std::runtime_error("boundary species state does not match element connectivity");
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
	constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
	constexpr double outward_sign[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
	double flux = 0.0;
	for (std::size_t qi = 0; qi < 4; ++qi)
		for (std::size_t qj = 0; qj < 4; ++qj) {
			std::array<double, 3> coordinate{};
			coordinate[fixed_axis[face]] = fixed_value[face];
			coordinate[varying_axes[face][0]] = points[qi];
			coordinate[varying_axes[face][1]] = points[qj];
			const auto basis = EvaluateBoundaryBasis(
				element, coordinate[0], coordinate[1], coordinate[2]);
			std::array<double, 3> velocity{};
			double concentration = 0.0;
			for (std::size_t a = 0; a < nodal_flow.size(); ++a) {
				concentration += basis.value[a]*nodal_species[a];
				for (int component = 0; component < 3; ++component)
					velocity[component] += basis.value[a]*nodal_flow[a][component];
			}
			double normal_velocity = 0.0;
			for (int component = 0; component < 3; ++component)
				normal_velocity += velocity[component]
					*basis.inverse_jacobian[fixed_axis[face]][component];
			flux += outward_sign[face]*weights[qi]*weights[qj]*2.0
				*basis.determinant*normal_velocity*concentration;
		}
	return flux;
}

inline std::vector<double> IntegrateBoundaryPressureTraction(
	const Element& element, std::size_t face, double pressure)
{
	if (face >= 6) throw std::runtime_error("pressure-traction face index must be below six");
	if (!std::isfinite(pressure))
		throw std::runtime_error("pressure traction must be finite");
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
		0.6521451548625461, 0.3478548451374539}};
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
	constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
	constexpr double outward_sign[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
	std::vector<double> result(4*element.connectivity.size(), 0.0);
	for (std::size_t qi = 0; qi < 4; ++qi)
		for (std::size_t qj = 0; qj < 4; ++qj) {
			std::array<double, 3> coordinate{};
			coordinate[fixed_axis[face]] = fixed_value[face];
			coordinate[varying_axes[face][0]] = points[qi];
			coordinate[varying_axes[face][1]] = points[qj];
			const auto basis = EvaluateBoundaryBasis(
				element, coordinate[0], coordinate[1], coordinate[2]);
			for (std::size_t a = 0; a < element.connectivity.size(); ++a)
				for (int component = 0; component < 3; ++component)
					result[4*a+component] -= pressure*outward_sign[face]
						*weights[qi]*weights[qj]*2.0*basis.determinant
						*basis.inverse_jacobian[fixed_axis[face]][component]
						*basis.value[a];
		}
	return result;
}

} // namespace iga

#endif

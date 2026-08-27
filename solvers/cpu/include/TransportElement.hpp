#ifndef TRANSPORT_ELEMENT_HPP
#define TRANSPORT_ELEMENT_HPP

#include "CaseInput.hpp"
#include "IgaDatabase.hpp"

#include <petscsys.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct BasisValues {
	std::vector<double> value;
	std::vector<std::array<double, 3>> gradient;
	std::vector<std::array<std::array<double, 3>, 3>> hessian;
	std::array<std::array<double, 3>, 3> inverse_jacobian{};
	double determinant = 0.0;
};

inline double Determinant(const double a[3][3])
{
	return a[0][0] * a[1][1] * a[2][2] + a[0][1] * a[1][2] * a[2][0] + a[0][2] * a[1][0] * a[2][1]
		- a[0][2] * a[1][1] * a[2][0] - a[0][0] * a[1][2] * a[2][1] - a[0][1] * a[1][0] * a[2][2];
}

inline void Inverse(const double a[3][3], double inverse[3][3])
{
	const auto determinant = Determinant(a);
	if (!std::isfinite(determinant) || std::abs(determinant) < 1e-14)
		throw std::runtime_error("singular element Jacobian");
	inverse[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / determinant;
	inverse[0][1] = (a[2][1] * a[0][2] - a[0][1] * a[2][2]) / determinant;
	inverse[0][2] = (a[0][1] * a[1][2] - a[1][1] * a[0][2]) / determinant;
	inverse[1][0] = (a[2][0] * a[1][2] - a[1][0] * a[2][2]) / determinant;
	inverse[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / determinant;
	inverse[1][2] = (a[1][0] * a[0][2] - a[0][0] * a[1][2]) / determinant;
	inverse[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / determinant;
	inverse[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / determinant;
	inverse[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / determinant;
}

inline double ElementJacobianDeterminant(const Element& element, double u, double v, double w)
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
	double jacobian[3][3]{};
	std::size_t p = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++p) {
				const double derivative[3] = {db[0][i]*b[1][j]*b[2][k], b[0][i]*db[1][j]*b[2][k], b[0][i]*b[1][j]*db[2][k]};
				for (int physical = 0; physical < 3; ++physical)
					for (int parameter = 0; parameter < 3; ++parameter)
						jacobian[physical][parameter] += element.bezier_points[p][physical] * derivative[parameter];
			}
	return 0.125 * Determinant(jacobian);
}

struct GeometryQuality {
	std::uint64_t bad_elements = 0;
	std::uint64_t bad_samples = 0;
	std::uint64_t first_bad_element = std::numeric_limits<std::uint64_t>::max();
	double minimum_determinant = std::numeric_limits<double>::infinity();
};

template <class OwnsElement>
inline GeometryQuality InspectGeometry(const std::vector<Element>& elements,
	OwnsElement&& owns_element, MPI_Comm communicator)
{
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
	GeometryQuality local, global;
	for (const auto& element : elements) {
		if (!owns_element(element)) continue;
		bool bad = false;
		for (double w : points) for (double v : points) for (double u : points) {
			const auto determinant = ElementJacobianDeterminant(element, u, v, w);
			local.minimum_determinant = std::min(local.minimum_determinant, determinant);
			if (!std::isfinite(determinant) || determinant <= 0.0) { ++local.bad_samples; bad = true; }
		}
		if (bad) { ++local.bad_elements; local.first_bad_element = std::min(local.first_bad_element, element.id); }
	}
	MPI_Allreduce(&local.bad_elements, &global.bad_elements, 1, MPI_UINT64_T, MPI_SUM, communicator);
	MPI_Allreduce(&local.bad_samples, &global.bad_samples, 1, MPI_UINT64_T, MPI_SUM, communicator);
	MPI_Allreduce(&local.first_bad_element, &global.first_bad_element, 1, MPI_UINT64_T, MPI_MIN, communicator);
	MPI_Allreduce(&local.minimum_determinant, &global.minimum_determinant, 1, MPI_DOUBLE, MPI_MIN, communicator);
	return global;
}

inline GeometryQuality InspectGeometry(const std::vector<Element>& elements,
	int rank, MPI_Comm communicator)
{
	return InspectGeometry(elements,
		[rank](const Element& element) { return element.owner == rank; }, communicator);
}

template <class OwnsElement>
inline void RequireValidGeometry(const std::vector<Element>& elements,
	OwnsElement&& owns_element, MPI_Comm communicator)
{
	const auto quality = InspectGeometry(
		elements, std::forward<OwnsElement>(owns_element), communicator);
	if (quality.bad_elements)
		throw std::runtime_error("invalid mesh: " + std::to_string(quality.bad_elements)
			+ " elements have non-positive Jacobians; first element "
			+ std::to_string(quality.first_bad_element) + ", minimum detJ "
			+ std::to_string(quality.minimum_determinant));
}

inline void RequireValidGeometry(const std::vector<Element>& elements, int rank, MPI_Comm communicator)
{
	RequireValidGeometry(elements,
		[rank](const Element& element) { return element.owner == rank; }, communicator);
}

inline BasisValues EvaluateBasis(const Element& element, double u, double v, double w, bool with_hessian = false)
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
	const double d2b[3][4] = {
		{6.0*(1.0-u), -12.0+18.0*u, 6.0-18.0*u, 6.0*u},
		{6.0*(1.0-v), -12.0+18.0*v, 6.0-18.0*v, 6.0*v},
		{6.0*(1.0-w), -12.0+18.0*w, 6.0-18.0*w, 6.0*w}
	};
	std::array<double, 64> bernstein{};
	std::array<std::array<double, 3>, 64> parametric{};
	std::array<std::array<std::array<double, 3>, 3>, 64> parametric_second{};
	std::size_t location = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++location) {
				bernstein[location] = b[0][i] * b[1][j] * b[2][k];
				parametric[location] = {db[0][i]*b[1][j]*b[2][k], b[0][i]*db[1][j]*b[2][k], b[0][i]*b[1][j]*db[2][k]};
				parametric_second[location][0] = {d2b[0][i]*b[1][j]*b[2][k], db[0][i]*db[1][j]*b[2][k], db[0][i]*b[1][j]*db[2][k]};
				parametric_second[location][1] = {db[0][i]*db[1][j]*b[2][k], b[0][i]*d2b[1][j]*b[2][k], b[0][i]*db[1][j]*db[2][k]};
				parametric_second[location][2] = {db[0][i]*b[1][j]*db[2][k], b[0][i]*db[1][j]*db[2][k], b[0][i]*b[1][j]*d2b[2][k]};
			}
	double jacobian[3][3]{};
	for (std::size_t p = 0; p < 64; ++p)
		for (int physical = 0; physical < 3; ++physical)
			for (int parameter = 0; parameter < 3; ++parameter)
				jacobian[physical][parameter] += element.bezier_points[p][physical] * parametric[p][parameter];
	double inverse[3][3]{};
	Inverse(jacobian, inverse);
	std::array<std::array<double, 3>, 64> physical_gradient{};
	for (std::size_t p = 0; p < 64; ++p)
		for (int physical = 0; physical < 3; ++physical)
			for (int parameter = 0; parameter < 3; ++parameter)
				physical_gradient[p][physical] += parametric[p][parameter] * inverse[parameter][physical];
	BasisValues result;
	result.value.assign(element.connectivity.size(), 0.0);
	result.gradient.assign(element.connectivity.size(), {0.0, 0.0, 0.0});
	result.inverse_jacobian = {{{inverse[0][0], inverse[0][1], inverse[0][2]},
		{inverse[1][0], inverse[1][1], inverse[1][2]}, {inverse[2][0], inverse[2][1], inverse[2][2]}}};
	result.determinant = 0.125 * Determinant(jacobian);
	if (!(result.determinant > 0.0)) throw std::runtime_error("non-positive element Jacobian");
	for (std::size_t a = 0; a < element.extraction.size(); ++a)
		for (std::size_t p = 0; p < 64; ++p) {
			result.value[a] += element.extraction[a][p] * bernstein[p];
			for (int d = 0; d < 3; ++d) result.gradient[a][d] += element.extraction[a][p] * physical_gradient[p][d];
		}
	if (with_hessian) {
		result.hessian.assign(element.connectivity.size(), {});
		double geometry_second[3][3][3]{};
		for (int physical = 0; physical < 3; ++physical)
			for (std::size_t p = 0; p < 64; ++p)
				for (int a = 0; a < 3; ++a)
					for (int b_index = 0; b_index < 3; ++b_index)
						geometry_second[physical][a][b_index] += element.bezier_points[p][physical] * parametric_second[p][a][b_index];
		double inverse_second[3][3][3]{};
		for (int c = 0; c < 3; ++c)
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					for (int physical = 0; physical < 3; ++physical)
						for (int a = 0; a < 3; ++a)
							for (int b_index = 0; b_index < 3; ++b_index)
								inverse_second[c][i][j] -= geometry_second[physical][a][b_index] * inverse[a][i] * inverse[b_index][j] * inverse[c][physical];
		for (std::size_t basis_index = 0; basis_index < element.connectivity.size(); ++basis_index)
			for (std::size_t p = 0; p < 64; ++p) {
				const auto coefficient = element.extraction[basis_index][p];
				if (coefficient == 0.0) continue;
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j) {
						double value = 0.0;
						for (int a = 0; a < 3; ++a) {
							for (int b_index = 0; b_index < 3; ++b_index)
								value += parametric_second[p][a][b_index] * inverse[a][i] * inverse[b_index][j];
							value += parametric[p][a] * inverse_second[a][i][j];
						}
						result.hessian[basis_index][i][j] += coefficient * value;
					}
			}
	}
	return result;
}

struct TransportMatrices {
	std::vector<PetscScalar> left;
	std::vector<PetscScalar> previous;
};

inline TransportMatrices BuildTransportElement(const Element& element,
	const std::vector<std::array<double, 3>>& nodal_velocity, const TransportParameters& p)
{
	constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
	constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461, 0.6521451548625461, 0.3478548451374539}};
	const auto nen = element.connectivity.size();
	const auto ndof = 2 * nen;
	TransportMatrices matrices{std::vector<PetscScalar>(ndof * ndof, 0.0), std::vector<PetscScalar>(ndof * ndof, 0.0)};
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz]);
				const auto measure = weights[qx] * weights[qy] * weights[qz] * basis.determinant;
				std::array<double, 3> velocity{};
				for (std::size_t a = 0; a < nen; ++a)
					for (int d = 0; d < 3; ++d)
						velocity[d] += basis.value[a] * nodal_velocity.at(static_cast<std::size_t>(element.connectivity[a]))[d];
				double inverse_length = 0.0;
				for (std::size_t a = 0; a < nen; ++a)
					inverse_length += std::abs(velocity[0]*basis.gradient[a][0] + velocity[1]*basis.gradient[a][1] + velocity[2]*basis.gradient[a][2]);
				const auto tau_space = inverse_length > 0.0 ? 1.0 / inverse_length : 0.0;
				const auto tau_time = p.dt / 2.0;
				const auto tau = tau_space > 0.0 ? 1.0 / std::sqrt(1.0/(tau_space*tau_space) + 1.0/(tau_time*tau_time)) : 0.0;
				std::vector<double> test_plus(nen);
				for (std::size_t a = 0; a < nen; ++a)
					test_plus[a] = basis.value[a] + tau * (velocity[0]*basis.gradient[a][0] + velocity[1]*basis.gradient[a][1] + velocity[2]*basis.gradient[a][2]);
				for (std::size_t a = 0; a < nen; ++a)
					for (std::size_t b = 0; b < nen; ++b) {
						const auto r0 = 2*a, rp = 2*a+1, c0 = 2*b, cp = 2*b+1;
						const auto gradient_dot = basis.gradient[a][0]*basis.gradient[b][0] + basis.gradient[a][1]*basis.gradient[b][1] + basis.gradient[a][2]*basis.gradient[b][2];
						const auto advection = velocity[0]*basis.gradient[b][0] + velocity[1]*basis.gradient[b][1] + velocity[2]*basis.gradient[b][2];
						matrices.left[r0*ndof+c0] += ((1.0+p.dt*(p.kplus+p.kminus))*basis.value[a]*basis.value[b] + p.dt*p.diffusion*gradient_dot)*measure;
						matrices.left[r0*ndof+cp] += -p.detach_plus*p.dt*basis.value[a]*basis.value[b]*measure;
						matrices.left[rp*ndof+c0] += -p.kplus*p.dt*test_plus[a]*basis.value[b]*measure;
						matrices.left[rp*ndof+cp] += ((1.0+p.dt*p.detach_plus)*test_plus[a]*basis.value[b] + p.dt*test_plus[a]*advection
							+ p.dt*p.vplus*p.artificial_diffusion*test_plus[a]*gradient_dot)*measure;
						matrices.previous[r0*ndof+c0] += basis.value[a]*basis.value[b]*measure;
						matrices.previous[rp*ndof+cp] += test_plus[a]*basis.value[b]*measure;
					}
			}
	return matrices;
}

} // namespace iga

#endif

#ifndef IGA_QUADRATURE_HPP
#define IGA_QUADRATURE_HPP

#include "ElementGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

struct VolumeQuadraturePoint {
	std::array<double, 3> parametric{};
	double weight = 0.0; // Dimensionless dxi deta dzeta on [0,1]^3.
};

struct SurfaceQuadraturePoint {
	std::array<double, 3> parametric{};
	std::array<double, 3> physical{};
	std::array<double, 3> normal{}; // Unit geometric outward normal when known.
	double weight = 0.0; // Already-physical dA in m^2.
	int boundary_id = -1;
};

class VolumeQuadratureRule {
public:
	explicit VolumeQuadratureRule(std::vector<VolumeQuadraturePoint> points = {})
		: points_(std::move(points)) {}

	const std::vector<VolumeQuadraturePoint>& Points() const noexcept { return points_; }

private:
	std::vector<VolumeQuadraturePoint> points_;
};

class SurfaceQuadratureRule {
public:
	explicit SurfaceQuadratureRule(std::vector<SurfaceQuadraturePoint> points = {})
		: points_(std::move(points)) {}

	const std::vector<SurfaceQuadraturePoint>& Points() const noexcept { return points_; }

private:
	std::vector<SurfaceQuadraturePoint> points_;
};

class VolumeQuadratureProvider {
public:
	virtual ~VolumeQuadratureProvider() = default;
	virtual const VolumeQuadratureRule& Rule() const noexcept = 0;
};

class SurfaceQuadratureProvider {
public:
	virtual ~SurfaceQuadratureProvider() = default;
	virtual const SurfaceQuadratureRule& Rule() const noexcept = 0;
};

inline constexpr std::array<double, 4> kGaussFourPoints{{0.06943184420297371,
	0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
inline constexpr std::array<double, 4> kGaussFourWeights{{0.3478548451374539,
	0.6521451548625461, 0.6521451548625461, 0.3478548451374539}};

inline bool QuadratureFinite(const std::array<double, 3>& value)
{
	return std::isfinite(value[0]) && std::isfinite(value[1]) && std::isfinite(value[2]);
}

inline bool QuadratureClose(double first, double second)
{
	return std::abs(first-second) <= 1.0e-11*std::max({1.0, std::abs(first), std::abs(second)});
}

inline void ValidateVolumeQuadratureRule(const Element& element,
	const VolumeQuadratureRule& rule)
{
	if (rule.Points().empty())
		throw std::runtime_error("active volume quadrature rule must not be empty");
	for (const auto& point : rule.Points()) {
		if (!QuadratureFinite(point.parametric)
			|| point.parametric[0] < 0.0 || point.parametric[0] > 1.0
			|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0
			|| point.parametric[2] < 0.0 || point.parametric[2] > 1.0)
			throw std::runtime_error("volume quadrature parametric coordinate is not finite in [0,1]^3");
		if (!std::isfinite(point.weight) || !(point.weight > 0.0))
			throw std::runtime_error("volume quadrature reference weight must be finite and positive");
		const auto geometry = EvaluateElementGeometry(element, point.parametric);
		if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0))
			throw std::runtime_error("volume quadrature point has a non-positive raw element Jacobian");
	}
}

inline int SurfaceFaceForParametric(const std::array<double, 3>& coordinate)
{
	if (QuadratureClose(coordinate[2], 0.0)) return 0;
	if (QuadratureClose(coordinate[1], 0.0)) return 1;
	if (QuadratureClose(coordinate[0], 1.0)) return 2;
	if (QuadratureClose(coordinate[1], 1.0)) return 3;
	if (QuadratureClose(coordinate[0], 0.0)) return 4;
	if (QuadratureClose(coordinate[2], 1.0)) return 5;
	throw std::runtime_error("surface quadrature point is not on a reference-cube face");
}

inline std::array<double, 3> SurfaceOutwardNormal(const ElementGeometryPoint& geometry,
	int face)
{
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	constexpr double outward_sign[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
	const auto inverse = InvertElementJacobian(geometry.jacobian, geometry.raw_determinant);
	std::array<double, 3> normal{};
	for (int physical = 0; physical < 3; ++physical)
		normal[physical] = outward_sign[face]*inverse[fixed_axis[face]][physical];
	const double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]
		+ normal[2]*normal[2]);
	if (!std::isfinite(norm) || !(norm > 0.0))
		throw std::runtime_error("surface quadrature normal is singular");
	for (double& value : normal) value /= norm;
	return normal;
}

inline double SurfaceJacobianMeasure(const ElementGeometryPoint& geometry, int face)
{
	constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
	const auto inverse = InvertElementJacobian(geometry.jacobian, geometry.raw_determinant);
	double inverse_normal_squared = 0.0;
	for (int physical = 0; physical < 3; ++physical)
		inverse_normal_squared += inverse[fixed_axis[face]][physical]
			*inverse[fixed_axis[face]][physical];
	return geometry.raw_determinant*std::sqrt(inverse_normal_squared);
}

inline void ValidateSurfaceQuadratureRule(const Element& element,
	const SurfaceQuadratureRule& rule)
{
	for (const auto& point : rule.Points()) {
		if (!QuadratureFinite(point.parametric)
			|| point.parametric[0] < 0.0 || point.parametric[0] > 1.0
			|| point.parametric[1] < 0.0 || point.parametric[1] > 1.0
			|| point.parametric[2] < 0.0 || point.parametric[2] > 1.0)
			throw std::runtime_error("surface quadrature parametric coordinate is not finite in [0,1]^3");
		if (!QuadratureFinite(point.physical))
			throw std::runtime_error("surface quadrature physical coordinate must be finite");
		if (!QuadratureFinite(point.normal))
			throw std::runtime_error("surface quadrature normal must be finite");
		if (!std::isfinite(point.weight) || !(point.weight > 0.0))
			throw std::runtime_error("surface quadrature physical weight must be finite and positive");
		if (point.boundary_id < 0)
			throw std::runtime_error("surface quadrature boundary id must be nonnegative");
		const auto geometry = EvaluateElementGeometry(element, point.parametric);
		if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0))
			throw std::runtime_error("surface quadrature point has a non-positive raw element Jacobian");
		for (int physical = 0; physical < 3; ++physical)
			if (!QuadratureClose(point.physical[physical], geometry.physical[physical]))
				throw std::runtime_error("surface quadrature physical point is inconsistent with element mapping");
		const double normal_norm = std::sqrt(point.normal[0]*point.normal[0]
			+ point.normal[1]*point.normal[1]
			+ point.normal[2]*point.normal[2]);
		if (!QuadratureClose(normal_norm, 1.0))
			throw std::runtime_error("surface quadrature normal is not approximately unit length");
	}
}

inline void ValidateBodyFittedSurfaceQuadratureRule(const Element& element,
	const SurfaceQuadratureRule& rule)
{
	ValidateSurfaceQuadratureRule(element, rule);
	for (const auto& point : rule.Points()) {
		const auto geometry = EvaluateElementGeometry(element, point.parametric);
		const auto expected = SurfaceOutwardNormal(geometry, SurfaceFaceForParametric(point.parametric));
		double alignment = 0.0;
		for (int physical = 0; physical < 3; ++physical)
			alignment += point.normal[physical]*expected[physical];
		if (alignment < 1.0-1.0e-10)
			throw std::runtime_error("surface quadrature normal is not geometrically outward");
	}
}

class FullCell4x4x4VolumeQuadratureProvider final : public VolumeQuadratureProvider {
public:
	explicit FullCell4x4x4VolumeQuadratureProvider(const Element& element)
		: rule_(BuildRule())
	{
		ValidateVolumeQuadratureRule(element, rule_);
	}

	const VolumeQuadratureRule& Rule() const noexcept override { return rule_; }

private:
	static VolumeQuadratureRule BuildRule()
	{
		std::vector<VolumeQuadraturePoint> points;
		points.reserve(64);
		for (std::size_t qz = 0; qz < 4; ++qz)
			for (std::size_t qy = 0; qy < 4; ++qy)
				for (std::size_t qx = 0; qx < 4; ++qx)
					points.push_back({{kGaussFourPoints[qx], kGaussFourPoints[qy],
						kGaussFourPoints[qz]}, kGaussFourWeights[qx]*kGaussFourWeights[qy]
							*kGaussFourWeights[qz]/8.0});
		return VolumeQuadratureRule(std::move(points));
	}

	VolumeQuadratureRule rule_;
};

class BodyFittedSurface4x4QuadratureProvider final : public SurfaceQuadratureProvider {
public:
	explicit BodyFittedSurface4x4QuadratureProvider(const Element& element)
		: rule_(BuildRule(element))
	{
		ValidateBodyFittedSurfaceQuadratureRule(element, rule_);
	}

	const SurfaceQuadratureRule& Rule() const noexcept override { return rule_; }

private:
	static SurfaceQuadratureRule BuildRule(const Element& element)
	{
		constexpr int fixed_axis[6] = {2, 1, 0, 1, 0, 2};
		constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2},
			{0, 2}, {1, 2}, {0, 1}};
		constexpr double fixed_value[6] = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
		std::vector<SurfaceQuadraturePoint> points;
		for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
			const int boundary_id = element.boundary_labels[face];
			if (boundary_id < 0) continue;
			for (std::size_t qi = 0; qi < 4; ++qi)
				for (std::size_t qj = 0; qj < 4; ++qj) {
					std::array<double, 3> coordinate{};
					coordinate[fixed_axis[face]] = fixed_value[face];
					coordinate[varying_axes[face][0]] = kGaussFourPoints[qi];
					coordinate[varying_axes[face][1]] = kGaussFourPoints[qj];
					const auto geometry = EvaluateElementGeometry(element, coordinate);
					const double physical_weight = kGaussFourWeights[qi]*kGaussFourWeights[qj]
						*SurfaceJacobianMeasure(geometry, static_cast<int>(face))/4.0;
					points.push_back({coordinate, geometry.physical,
						SurfaceOutwardNormal(geometry, static_cast<int>(face)),
						physical_weight, boundary_id});
				}
		}
		return SurfaceQuadratureRule(std::move(points));
	}

	SurfaceQuadratureRule rule_;
};

} // namespace iga

#endif

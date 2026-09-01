#include "Quadrature.hpp"
#include "TransportElement.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace {

bool Close(double first, double second)
{
	return std::abs(first-second) <= 2.0e-12*std::max({1.0, std::abs(first), std::abs(second)});
}

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
}

template <class Map>
iga::Element MakeElement(Map&& map)
{
	iga::Element element;
	element.connectivity.resize(64);
	element.extraction.resize(64);
	std::size_t point = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++point) {
				element.connectivity[point] = static_cast<std::int32_t>(point);
				element.extraction[point].fill(0.0);
				element.extraction[point][point] = 1.0;
				element.bezier_points[point] = map(i/3.0, j/3.0, k/3.0);
			}
	return element;
}

double SurfaceArea(const iga::SurfaceQuadratureRule& rule, int boundary_id)
{
	double area = 0.0;
	for (const auto& point : rule.Points())
		if (point.boundary_id == boundary_id) area += point.weight;
	return area;
}

} // namespace

int main()
{
	const auto unit = MakeElement([](double u, double v, double w) {
		return std::array<double, 3>{{u, v, w}};
	});
	iga::FullCell4x4x4VolumeQuadratureProvider volume_provider(unit);
	const iga::VolumeQuadratureProvider& volume_interface = volume_provider;
	const auto& volume = volume_interface.Rule();
	assert(&volume == &volume_provider.Rule());
	assert(volume.Points().size() == 64);
	double reference_volume = 0.0;
	for (std::size_t qz = 0; qz < 4; ++qz)
		for (std::size_t qy = 0; qy < 4; ++qy)
			for (std::size_t qx = 0; qx < 4; ++qx) {
				const auto index = qz*16+qy*4+qx;
				const auto& point = volume.Points()[index];
				assert(Close(point.parametric[0], iga::kGaussFourPoints[qx])
					&& Close(point.parametric[1], iga::kGaussFourPoints[qy])
					&& Close(point.parametric[2], iga::kGaussFourPoints[qz]));
				reference_volume += point.weight;
			}
	assert(Close(reference_volume, 1.0));
	const auto coordinate = std::array<double, 3>{{0.2, 0.3, 0.4}};
	const auto geometry = iga::EvaluateElementGeometry(unit, coordinate);
	const auto basis = iga::EvaluateBasis(unit, coordinate[0], coordinate[1], coordinate[2]);
	assert(Close(geometry.raw_determinant, 1.0));
	assert(Close(basis.raw_determinant, geometry.raw_determinant));
	assert(Close(basis.determinant, geometry.raw_determinant/8.0));
	assert(Close(iga::ElementJacobianDeterminant(unit, coordinate[0], coordinate[1], coordinate[2]),
		geometry.raw_determinant/8.0));

	auto labelled_unit = unit;
	for (int face = 0; face < 6; ++face) labelled_unit.boundary_labels[face] = 20+face;
	iga::BodyFittedSurface4x4QuadratureProvider unit_surface_provider(labelled_unit);
	const iga::SurfaceQuadratureProvider& surface_interface = unit_surface_provider;
	const auto& unit_surface = surface_interface.Rule();
	assert(&unit_surface == &unit_surface_provider.Rule());
	assert(unit_surface.Points().size() == 96);
	const std::array<std::array<double, 3>, 6> unit_normals{{
		{{0.0, 0.0, -1.0}}, {{0.0, -1.0, 0.0}}, {{1.0, 0.0, 0.0}},
		{{0.0, 1.0, 0.0}}, {{-1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}};
	for (int face = 0; face < 6; ++face) {
		assert(Close(SurfaceArea(unit_surface, 20+face), 1.0));
		const auto& first = unit_surface.Points()[static_cast<std::size_t>(face)*16];
		assert(first.boundary_id == 20+face);
		for (int component = 0; component < 3; ++component)
			assert(Close(first.normal[component], unit_normals[face][component]));
	}
	const auto& face_zero_first = unit_surface.Points().front();
	assert(Close(face_zero_first.parametric[0], iga::kGaussFourPoints[0])
		&& Close(face_zero_first.parametric[1], iga::kGaussFourPoints[0])
		&& Close(face_zero_first.parametric[2], 0.0));
	const auto& face_one_first = unit_surface.Points()[16];
	assert(Close(face_one_first.parametric[0], iga::kGaussFourPoints[0])
		&& Close(face_one_first.parametric[1], 0.0)
		&& Close(face_one_first.parametric[2], iga::kGaussFourPoints[0]));

	const std::array<std::array<double, 3>, 3> affine{{
		{{2.0, 0.3, 0.2}}, {{0.1, 3.0, 0.4}}, {{0.2, 0.3, 4.0}}}};
	auto skewed = MakeElement([&](double u, double v, double w) {
		return std::array<double, 3>{{affine[0][0]*u+affine[0][1]*v+affine[0][2]*w,
			affine[1][0]*u+affine[1][1]*v+affine[1][2]*w,
			affine[2][0]*u+affine[2][1]*v+affine[2][2]*w}};
	});
	for (int face = 0; face < 6; ++face) skewed.boundary_labels[face] = 40+face;
	const double known_affine_determinant = 23.55;
	const auto skew_coordinate = std::array<double, 3>{{0.21, 0.37, 0.64}};
	const auto skew_geometry = iga::EvaluateElementGeometry(skewed, skew_coordinate);
	const auto skew_basis = iga::EvaluateBasis(skewed, skew_coordinate[0],
		skew_coordinate[1], skew_coordinate[2]);
	assert(Close(skew_geometry.raw_determinant, known_affine_determinant));
	assert(Close(skew_basis.raw_determinant, known_affine_determinant));
	assert(Close(skew_basis.raw_determinant, skew_geometry.raw_determinant));
	iga::FullCell4x4x4VolumeQuadratureProvider skewed_volume_provider(skewed);
	double integrated_skewed_volume = 0.0;
	for (const auto& point : skewed_volume_provider.Rule().Points())
		integrated_skewed_volume += point.weight
			*iga::EvaluateElementGeometry(skewed, point.parametric).raw_determinant;
	assert(Close(integrated_skewed_volume, known_affine_determinant));
	iga::BodyFittedSurface4x4QuadratureProvider skewed_provider(skewed);
	const auto& skewed_rule = skewed_provider.Rule();
	for (int face = 0; face < 6; ++face) {
		const auto geometry_point = iga::EvaluateElementGeometry(skewed,
			skewed_rule.Points()[static_cast<std::size_t>(face)*16].parametric);
		assert(Close(SurfaceArea(skewed_rule, 40+face),
			iga::SurfaceJacobianMeasure(geometry_point, face)));
	}

	auto curved = MakeElement([](double u, double v, double w) {
		return std::array<double, 3>{{u+0.08*u*v, v+0.05*v*w, w+0.04*u*w}};
	});
	for (int face = 0; face < 6; ++face) curved.boundary_labels[face] = 60+face;
	iga::BodyFittedSurface4x4QuadratureProvider curved_provider(curved);
	for (const auto& point : curved_provider.Rule().Points()) {
		const double norm = std::sqrt(point.normal[0]*point.normal[0]
			+ point.normal[1]*point.normal[1]
			+ point.normal[2]*point.normal[2]);
		assert(Close(norm, 1.0) && point.weight > 0.0);
	}

	RequireRejected([&] { iga::ValidateVolumeQuadratureRule(unit, iga::VolumeQuadratureRule{}); });
	RequireRejected([&] {
		iga::ValidateVolumeQuadratureRule(unit,
			iga::VolumeQuadratureRule({{{{1.1, 0.5, 0.5}}, 1.0}}));
	});
	RequireRejected([&] {
		iga::ValidateVolumeQuadratureRule(unit,
			iga::VolumeQuadratureRule({{{{0.5, 0.5, 0.5}}, -1.0}}));
	});
	auto inverted = unit;
	for (auto& point : inverted.bezier_points) point[0] = 1.0-point[0];
	RequireRejected([&] { iga::FullCell4x4x4VolumeQuadratureProvider invalid(inverted); });
	auto invalid_surface_points = unit_surface.Points();
	invalid_surface_points.front().physical[0] += 0.1;
	RequireRejected([&] {
		iga::ValidateSurfaceQuadratureRule(labelled_unit,
			iga::SurfaceQuadratureRule(std::move(invalid_surface_points)));
	});
	invalid_surface_points = unit_surface.Points();
	invalid_surface_points.front().normal = {{0.0, 0.0, 0.0}};
	RequireRejected([&] {
		iga::ValidateSurfaceQuadratureRule(labelled_unit,
			iga::SurfaceQuadratureRule(std::move(invalid_surface_points)));
	});
	invalid_surface_points = unit_surface.Points();
	invalid_surface_points.front().normal[2] *= -1.0;
	iga::ValidateSurfaceQuadratureRule(labelled_unit,
		iga::SurfaceQuadratureRule(invalid_surface_points));
	RequireRejected([&] {
		iga::ValidateBodyFittedSurfaceQuadratureRule(labelled_unit,
			iga::SurfaceQuadratureRule(std::move(invalid_surface_points)));
	});
	invalid_surface_points = unit_surface.Points();
	invalid_surface_points.front().weight = 0.0;
	RequireRejected([&] {
		iga::ValidateSurfaceQuadratureRule(labelled_unit,
			iga::SurfaceQuadratureRule(std::move(invalid_surface_points)));
	});
	const auto interior_coordinate = std::array<double, 3>{{0.5, 0.5, 0.5}};
	const auto interior_geometry = iga::EvaluateElementGeometry(unit, interior_coordinate);
	iga::SurfaceQuadratureRule immersed({{interior_coordinate, interior_geometry.physical,
		{{1.0, 0.0, 0.0}}, 0.25, 99}});
	iga::ValidateSurfaceQuadratureRule(unit, immersed);
	RequireRejected([&] { iga::ValidateBodyFittedSurfaceQuadratureRule(unit, immersed); });
	invalid_surface_points = unit_surface.Points();
	invalid_surface_points.front().boundary_id = -1;
	RequireRejected([&] {
		iga::ValidateSurfaceQuadratureRule(labelled_unit,
			iga::SurfaceQuadratureRule(std::move(invalid_surface_points)));
	});

	std::cout << "CPU quadrature provider tests passed\n";
}

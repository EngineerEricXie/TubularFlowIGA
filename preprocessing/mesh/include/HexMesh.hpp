#pragma once

#include "Geometry.hpp"

#include <array>
#include <filesystem>
#include <limits>
#include <vector>

namespace tubular {

using Hex = std::array<int,8>;

struct QualityResult
{
	double minimum_determinant = std::numeric_limits<double>::infinity();
	double minimum_scaled_jacobian = std::numeric_limits<double>::infinity();
	std::size_t bad_elements = 0;
	int first_bad_element = -1;
};

QualityResult EvaluateHexQuality(
	const std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<int>* subset = nullptr);

double MovePointSafely(
	std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<std::vector<int>>& incident_elements,
	int point_index,
	const Vec3& candidate,
	double required_scaled_jacobian);

std::vector<std::vector<int>> BuildIncidentElements(
	std::size_t point_count,
	const std::vector<Hex>& elements);

} // namespace tubular

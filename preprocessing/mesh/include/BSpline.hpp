#pragma once

#include "Geometry.hpp"

#include <vector>

namespace tubular {

struct CurveSample
{
	Vec3 point;
	double diameter = 0.0;
	Vec3 tangent;
};

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	double segment_length,
	int mode);

} // namespace tubular

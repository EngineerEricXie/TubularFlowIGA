#pragma once

#include "Geometry.hpp"

#include <string>
#include <vector>

namespace tubular {

struct CurveSample
{
	Vec3 point;
	double diameter = 0.0;
	Vec3 tangent;
};

struct BranchSamplingOptions
{
	double target_spacing = 0.0;
	double max_spacing_over_diameter = 1.0;
	double max_turn_degrees = 12.0;
	double max_diameter_change_fraction = 0.15;
	double upstream_clearance_over_diameter = 1.0;
	double downstream_clearance_over_diameter = 1.5;
};

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	const BranchSamplingOptions& options,
	int mode,
	const std::string& context);

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	double segment_length,
	int mode);

} // namespace tubular

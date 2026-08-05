#include "BSpline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace tubular {
namespace {

using Value = std::array<double,4>;

Value Add(Value a, const Value& b)
{
	for (int i=0; i<4; ++i) a[i] += b[i];
	return a;
}

Value Scale(Value a, double s)
{
	for (double& value : a) value *= s;
	return a;
}

Value Evaluate(const std::vector<double>& knots, const std::vector<Value>& control, int degree, double u)
{
	const int n = static_cast<int>(control.size())-1;
	if (n < degree || static_cast<int>(knots.size()) != n+degree+2)
		throw std::runtime_error("invalid B-spline representation");
	int span = degree;
	if (u >= knots[n+1]) span = n;
	else {
		for (int i=degree; i<=n; ++i)
			if (u >= knots[i] && u < knots[i+1]) { span = i; break; }
	}
	std::vector<Value> d(degree+1);
	for (int j=0; j<=degree; ++j) d[j] = control[span-degree+j];
	for (int r=1; r<=degree; ++r) {
		for (int j=degree; j>=r; --j) {
			const int i = span-degree+j;
			const double denominator = knots[i+degree-r+1]-knots[i];
			const double alpha = denominator > 0.0 ? (u-knots[i])/denominator : 0.0;
			d[j] = Add(Scale(d[j-1], 1.0-alpha), Scale(d[j], alpha));
		}
	}
	return d[degree];
}

std::vector<double> Colon(double start, double step, double end)
{
	std::vector<double> values;
	if (!std::isfinite(start) || !std::isfinite(step) || !std::isfinite(end) || step <= 0.0)
		return values;
	for (std::size_t i=0;; ++i) {
		const double value = start+static_cast<double>(i)*step;
		if (value > end+1.0e-12*std::max(1.0, std::abs(end))) break;
		values.push_back(std::min(value, end));
		if (i > 100000000) throw std::runtime_error("B-spline sampling generated too many points");
	}
	return values;
}

} // namespace

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	double segment_length,
	int mode)
{
	if (points.size() != diameters.size() || points.size() < 2)
		throw std::runtime_error("SampleBranch requires matching point and diameter arrays");
	double length = 0.0;
	for (std::size_t i=1; i<points.size(); ++i) length += Norm(points[i]-points[i-1]);
	if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(segment_length) || segment_length <= 0.0)
		throw std::runtime_error("SampleBranch requires positive curve and segment lengths");
	for (double value : diameters)
		if (!std::isfinite(value) || value <= 0.0) throw std::runtime_error("SampleBranch received an invalid diameter");

	const double step = segment_length/length;
	std::vector<double> samples;
	if (mode == 1) {
		samples.push_back(0.0);
		auto middle = Colon(1.5*diameters.front()/length, step, 1.0-diameters.back()/length);
		samples.insert(samples.end(), middle.begin(), middle.end());
		samples.push_back(1.0);
	} else if (mode == 2) {
		samples.push_back(0.0);
		auto rest = Colon(1.5*diameters.front()/length, step, 1.0);
		samples.insert(samples.end(), rest.begin(), rest.end());
		if (1.0-samples.back() > 0.05) samples.push_back(1.0);
	} else if (mode == 3) {
		samples = Colon(0.0, step, 1.0-diameters.back()/length);
		samples.push_back(1.0);
	} else if (mode == 4) {
		samples = Colon(0.0, step, 1.0);
		if (samples.empty() || 1.0-samples.back() > 0.05) samples.push_back(1.0);
	} else {
		throw std::runtime_error("unknown branch sampling mode");
	}
	if (samples.size() < 4) samples = {0.0, 0.45, 0.67, 1.0};

	const int count = static_cast<int>(points.size());
	const int degree = count > 3 ? 3 : count-1;
	std::vector<double> chord(count, 0.0);
	for (int i=1; i<count; ++i) chord[i] = chord[i-1]+Norm(points[i]-points[i-1])/length;
	chord.back() = 1.0;
	std::vector<double> knots(degree+1, 0.0);
	if (count > 3)
		for (int i=2; i<=count-3; ++i) knots.push_back(chord[i]);
	knots.insert(knots.end(), degree+1, 1.0);

	std::vector<Value> control(count);
	for (int i=0; i<count; ++i)
		control[i] = {points[i].x, points[i].y, points[i].z, diameters[i]};

	std::vector<Value> derivative_control;
	std::vector<double> derivative_knots;
	if (degree > 0) {
		derivative_control.resize(count-1);
		for (int i=0; i<count-1; ++i) {
			const double denominator = knots[i+degree+1]-knots[i+1];
			if (denominator <= 0.0) throw std::runtime_error("degenerate B-spline derivative knot span");
			for (int component=0; component<4; ++component)
				derivative_control[i][component] = degree*(control[i+1][component]-control[i][component])/denominator;
		}
		derivative_knots.assign(knots.begin()+1, knots.end()-1);
	}

	std::vector<CurveSample> result;
	result.reserve(samples.size());
	for (double u : samples) {
		u = std::max(0.0, std::min(1.0, u));
		const Value value = Evaluate(knots, control, degree, u);
		const Value tangent = Evaluate(derivative_knots, derivative_control, degree-1, u);
		result.push_back({{value[0], value[1], value[2]}, value[3], {tangent[0], tangent[1], tangent[2]}});
	}
	return result;
}

} // namespace tubular

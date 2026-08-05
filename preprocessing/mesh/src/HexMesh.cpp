#include "HexMesh.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tubular {
namespace {

QualityResult EvaluateOne(const std::vector<Vec3>& points, const std::vector<Hex>& elements, const std::vector<int>& indices)
{
	static const std::array<double,6> samples{
		-1.0, -0.861136311594053, -0.339981043584856,
		0.339981043584856, 0.861136311594053, 1.0};
	static const std::array<std::array<int,3>,8> signs{{
		{{-1,-1,-1}}, {{1,-1,-1}}, {{1,1,-1}}, {{-1,1,-1}},
		{{-1,-1,1}}, {{1,-1,1}}, {{1,1,1}}, {{-1,1,1}}
	}};
	QualityResult result;
	for (int element_index : indices) {
		if (element_index < 0 || static_cast<std::size_t>(element_index) >= elements.size())
			throw std::runtime_error("hex quality subset contains an invalid element index");
		const auto& element = elements[element_index];
		for (int node : element)
			if (node < 0 || static_cast<std::size_t>(node) >= points.size())
				throw std::runtime_error("hex element contains an invalid point index");
		bool bad = false;
		for (double r : samples) for (double s : samples) for (double t : samples) {
			std::array<Vec3,3> jacobian{};
			for (int a=0; a<8; ++a) {
				const double dr = 0.125*signs[a][0]*(1.0+signs[a][1]*s)*(1.0+signs[a][2]*t);
				const double ds = 0.125*signs[a][1]*(1.0+signs[a][0]*r)*(1.0+signs[a][2]*t);
				const double dt = 0.125*signs[a][2]*(1.0+signs[a][0]*r)*(1.0+signs[a][1]*s);
				const Vec3& p = points[element[a]];
				jacobian[0] += p*dr;
				jacobian[1] += p*ds;
				jacobian[2] += p*dt;
			}
			const double determinant = Determinant(jacobian);
			const double scale = Norm(jacobian[0])*Norm(jacobian[1])*Norm(jacobian[2]);
			const double scaled = std::isfinite(scale) && scale > 1.0e-15
				? determinant/scale : -std::numeric_limits<double>::infinity();
			result.minimum_determinant = std::min(result.minimum_determinant, determinant);
			result.minimum_scaled_jacobian = std::min(result.minimum_scaled_jacobian, scaled);
			if (!std::isfinite(determinant) || determinant <= 0.0) bad = true;
		}
		if (bad) {
			++result.bad_elements;
			if (result.first_bad_element < 0) result.first_bad_element = element_index;
		}
	}
	return result;
}

} // namespace

QualityResult EvaluateHexQuality(
	const std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<int>* subset)
{
	std::vector<int> all;
	if (!subset) {
		all.resize(elements.size());
		for (std::size_t i=0; i<elements.size(); ++i) all[i] = static_cast<int>(i);
		subset = &all;
	}
	return EvaluateOne(points, elements, *subset);
}

std::vector<std::vector<int>> BuildIncidentElements(
	std::size_t point_count,
	const std::vector<Hex>& elements)
{
	std::vector<std::vector<int>> incident(point_count);
	for (std::size_t e=0; e<elements.size(); ++e)
		for (int point : elements[e]) {
			if (point < 0 || static_cast<std::size_t>(point) >= point_count)
				throw std::runtime_error("cannot build incidence for invalid connectivity");
			incident[point].push_back(static_cast<int>(e));
		}
	for (auto& list : incident) {
		std::sort(list.begin(), list.end());
		list.erase(std::unique(list.begin(), list.end()), list.end());
	}
	return incident;
}

double MovePointSafely(
	std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<std::vector<int>>& incident_elements,
	int point_index,
	const Vec3& candidate,
	double required_scaled_jacobian)
{
	if (point_index < 0 || static_cast<std::size_t>(point_index) >= points.size())
		throw std::runtime_error("point move index is out of range");
	if (!IsFinite(candidate)) throw std::runtime_error("point move candidate is not finite");
	const auto& incident = incident_elements.at(point_index);
	if (incident.empty()) throw std::runtime_error("point move has no incident volume element");
	const Vec3 original = points[point_index];
	static const std::array<double,8> alphas{{1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0}};
	for (double alpha : alphas) {
		points[point_index] = original+(candidate-original)*alpha;
		const auto quality = EvaluateHexQuality(points, elements, &incident);
		if (quality.bad_elements == 0 && quality.minimum_determinant > 0.0
			&& quality.minimum_scaled_jacobian >= required_scaled_jacobian)
			return alpha;
	}
	points[point_index] = original;
	throw std::runtime_error("no valid position found during point backtracking");
}

} // namespace tubular

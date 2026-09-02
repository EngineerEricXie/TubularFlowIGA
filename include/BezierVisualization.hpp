#ifndef IGA_BEZIER_VISUALIZATION_HPP
#define IGA_BEZIER_VISUALIZATION_HPP

#include "IgaDatabase.hpp"
#include "VtkOutput.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

constexpr std::uint8_t kVtkBezierHexahedron = 79;

struct BezierGeometryValidation {
	std::uint64_t elements = 0;
	std::uint64_t local_points = 0;
	std::uint64_t unique_points = 0;
	std::uint64_t shared_point_references = 0;
	std::uint64_t signature_coordinate_repairs = 0;
	std::uint64_t coincident_unmerged_points = 0;
	std::uint64_t collapsed_element_points = 0;
	std::uint64_t nonpositive_jacobian_elements = 0;
	std::uint64_t overlap_candidates = 0;
	std::uint64_t overlapping_element_pairs = 0;
	std::uint64_t first_bad_element = std::numeric_limits<std::uint64_t>::max();
	std::array<std::uint64_t, 2> first_overlap_pair{{
		std::numeric_limits<std::uint64_t>::max(),
		std::numeric_limits<std::uint64_t>::max()}};
	std::array<std::uint64_t, 2> first_coincident_element_pair{{
		std::numeric_limits<std::uint64_t>::max(),
		std::numeric_limits<std::uint64_t>::max()}};
	double minimum_jacobian = std::numeric_limits<double>::infinity();
	double maximum_signature_coordinate_difference = 0.0;
};

struct BezierVisualizationMesh {
	std::vector<std::array<double, 3>> points;
	std::vector<std::int64_t> connectivity;
	std::vector<std::int64_t> offsets;
	std::vector<std::uint8_t> types;
	std::vector<std::array<std::int32_t, 3>> higher_order_degrees;
	std::vector<std::uint64_t> element_ids;
	std::vector<std::int32_t> element_owners;
	std::vector<std::uint64_t> signature_offsets{0};
	std::vector<std::int32_t> signature_nodes;
	std::vector<double> signature_coefficients;
	BezierGeometryValidation validation;
};

inline std::filesystem::path BezierGeometryReportPath(
	const std::filesystem::path& output)
{
	auto result = output;
	result.replace_extension(".bezier_geometry.json");
	return result;
}

inline bool IsValidBezierGeometry(const BezierGeometryValidation& validation)
{
	return validation.collapsed_element_points == 0
		&& validation.nonpositive_jacobian_elements == 0
		&& validation.coincident_unmerged_points == 0
		&& validation.overlapping_element_pairs == 0;
}

inline void RequireValidBezierGeometry(const BezierGeometryValidation& validation)
{
	if (IsValidBezierGeometry(validation)) return;
	throw std::runtime_error("invalid Bezier visualization geometry: collapsed_points="
		+std::to_string(validation.collapsed_element_points)
		+" nonpositive_jacobian_elements="
		+std::to_string(validation.nonpositive_jacobian_elements)
		+" coincident_unmerged_points="
		+std::to_string(validation.coincident_unmerged_points)
		+" overlapping_element_pairs="
		+std::to_string(validation.overlapping_element_pairs));
}

inline void WriteBezierGeometryReport(const std::filesystem::path& path,
	const BezierGeometryValidation& validation)
{
	if (!path.parent_path().empty())
		std::filesystem::create_directories(path.parent_path());
	std::ofstream output(path);
	if (!output) throw std::runtime_error(
		"cannot create Bezier geometry report: "+path.string());
	const auto missing = std::numeric_limits<std::uint64_t>::max();
	output << std::setprecision(17)
		<< "{\n"
		<< "  \"valid\": " << (IsValidBezierGeometry(validation) ? "true" : "false") << ",\n"
		<< "  \"elements\": " << validation.elements << ",\n"
		<< "  \"local_points\": " << validation.local_points << ",\n"
		<< "  \"unique_points\": " << validation.unique_points << ",\n"
		<< "  \"shared_point_references\": " << validation.shared_point_references << ",\n"
		<< "  \"signature_coordinate_repairs\": "
		<< validation.signature_coordinate_repairs << ",\n"
		<< "  \"maximum_signature_coordinate_difference\": "
		<< validation.maximum_signature_coordinate_difference << ",\n"
		<< "  \"coincident_unmerged_points\": " << validation.coincident_unmerged_points << ",\n"
		<< "  \"collapsed_element_points\": " << validation.collapsed_element_points << ",\n"
		<< "  \"nonpositive_jacobian_elements\": "
		<< validation.nonpositive_jacobian_elements << ",\n"
		<< "  \"minimum_sampled_jacobian\": " << validation.minimum_jacobian << ",\n"
		<< "  \"overlap_candidates\": " << validation.overlap_candidates << ",\n"
		<< "  \"overlapping_element_pairs\": "
		<< validation.overlapping_element_pairs << ",\n"
		<< "  \"first_bad_element\": ";
	if (validation.first_bad_element == missing) output << "null";
	else output << validation.first_bad_element;
	output << ",\n  \"first_overlap_pair\": ";
	if (validation.first_overlap_pair[0] == missing) output << "null";
	else output << '[' << validation.first_overlap_pair[0] << ", "
		<< validation.first_overlap_pair[1] << ']';
	output << ",\n  \"first_coincident_element_pair\": ";
	if (validation.first_coincident_element_pair[0] == missing) output << "null";
	else output << '[' << validation.first_coincident_element_pair[0] << ", "
		<< validation.first_coincident_element_pair[1] << ']';
	output << "\n}\n";
	if (!output) throw std::runtime_error(
		"cannot write Bezier geometry report: "+path.string());
}

namespace detail {

using QuantizedSignature = std::vector<std::pair<std::int32_t, std::int64_t>>;

inline std::int64_t QuantizeCoefficient(double value, double scale)
{
	if (!std::isfinite(value)
		|| std::abs(value) > static_cast<double>(
			std::numeric_limits<std::int64_t>::max())/scale)
		throw std::runtime_error("Bezier extraction coefficient cannot be quantized");
	return static_cast<std::int64_t>(std::llround(value*scale));
}

inline std::uint64_t SignatureHash(const QuantizedSignature& signature)
{
	std::uint64_t result = 1469598103934665603ULL;
	for (const auto& entry : signature) {
		for (const auto value : {static_cast<std::uint64_t>(entry.first),
			static_cast<std::uint64_t>(entry.second)}) {
			result ^= value;
			result *= 1099511628211ULL;
		}
	}
	return result;
}

inline bool SignatureMatches(const BezierVisualizationMesh& mesh,
	std::size_t point, const QuantizedSignature& signature,
	double coefficient_quantization)
{
	const auto begin = mesh.signature_offsets[point];
	const auto end = mesh.signature_offsets[point+1];
	if (end-begin != signature.size()) return false;
	for (std::size_t index = 0; index < signature.size(); ++index) {
		const auto entry = begin+index;
		if (mesh.signature_nodes[entry] != signature[index].first
			|| QuantizeCoefficient(
				mesh.signature_coefficients[entry], coefficient_quantization)
				!= signature[index].second)
			return false;
	}
	return true;
}

struct Bounds {
	std::array<double, 3> minimum{{
		std::numeric_limits<double>::infinity(),
		std::numeric_limits<double>::infinity(),
		std::numeric_limits<double>::infinity()}};
	std::array<double, 3> maximum{{
		-std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity()}};
};

struct CoordinateKey {
	std::array<std::int64_t, 3> value{};
	bool operator==(const CoordinateKey& other) const { return value == other.value; }
};

struct CoordinateKeyHash {
	std::size_t operator()(const CoordinateKey& key) const
	{
		std::size_t result = 1469598103934665603ULL;
		for (const auto value : key.value) {
			result ^= static_cast<std::uint64_t>(value);
			result *= 1099511628211ULL;
		}
		return result;
	}
};

inline CoordinateKey CoordinateBucket(
	const std::array<double, 3>& point, double tolerance)
{
	CoordinateKey key;
	for (int direction = 0; direction < 3; ++direction)
		key.value[direction] = static_cast<std::int64_t>(
			std::llround(point[direction]/tolerance));
	return key;
}

inline double DistanceSquared(const std::array<double, 3>& left,
	const std::array<double, 3>& right)
{
	double result = 0.0;
	for (int direction = 0; direction < 3; ++direction) {
		const auto difference = left[direction]-right[direction];
		result += difference*difference;
	}
	return result;
}

inline Bounds ElementBounds(const Element& element)
{
	Bounds bounds;
	for (const auto& point : element.bezier_points)
		for (int direction = 0; direction < 3; ++direction) {
			bounds.minimum[direction] = std::min(
				bounds.minimum[direction], point[direction]);
			bounds.maximum[direction] = std::max(
				bounds.maximum[direction], point[direction]);
		}
	return bounds;
}

inline bool BoundsOverlap(const Bounds& left, const Bounds& right, double tolerance)
{
	for (int direction = 0; direction < 3; ++direction)
		if (std::min(left.maximum[direction], right.maximum[direction])
			-std::max(left.minimum[direction], right.minimum[direction]) <= tolerance)
			return false;
	return true;
}

inline void Bernstein(double coordinate, double values[4], double derivatives[4])
{
	const auto one_minus = 1.0-coordinate;
	values[0] = one_minus*one_minus*one_minus;
	values[1] = 3.0*one_minus*one_minus*coordinate;
	values[2] = 3.0*one_minus*coordinate*coordinate;
	values[3] = coordinate*coordinate*coordinate;
	derivatives[0] = -3.0*one_minus*one_minus;
	derivatives[1] = 3.0-12.0*coordinate+9.0*coordinate*coordinate;
	derivatives[2] = 3.0*(2.0-3.0*coordinate)*coordinate;
	derivatives[3] = 3.0*coordinate*coordinate;
}

inline std::array<double, 3> MapBezier(const Element& element,
	const std::array<double, 3>& parameter,
	std::array<std::array<double, 3>, 3>* jacobian = nullptr)
{
	double basis[3][4]{};
	double derivative[3][4]{};
	for (int direction = 0; direction < 3; ++direction)
		Bernstein(parameter[direction], basis[direction], derivative[direction]);
	std::array<double, 3> point{};
	if (jacobian) *jacobian = {};
	std::size_t location = 0;
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i, ++location) {
				const auto value = basis[0][i]*basis[1][j]*basis[2][k];
				const std::array<double, 3> gradient{{
					derivative[0][i]*basis[1][j]*basis[2][k],
					basis[0][i]*derivative[1][j]*basis[2][k],
					basis[0][i]*basis[1][j]*derivative[2][k]}};
				for (int physical = 0; physical < 3; ++physical) {
					point[physical] += element.bezier_points[location][physical]*value;
					if (jacobian)
						for (int direction = 0; direction < 3; ++direction)
							(*jacobian)[physical][direction]
								+= element.bezier_points[location][physical]
									*gradient[direction];
				}
			}
	return point;
}

inline double Determinant(const std::array<std::array<double, 3>, 3>& matrix)
{
	return matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])
		-matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])
		+matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0]);
}

inline bool Solve3x3(const std::array<std::array<double, 3>, 3>& matrix,
	const std::array<double, 3>& right, std::array<double, 3>& solution)
{
	const auto determinant = Determinant(matrix);
	if (std::abs(determinant) <= 1e-18) return false;
	for (int column = 0; column < 3; ++column) {
		auto replaced = matrix;
		for (int row = 0; row < 3; ++row) replaced[row][column] = right[row];
		solution[column] = Determinant(replaced)/determinant;
	}
	return true;
}

inline bool StrictlyInside(const Element& element, const Bounds& bounds,
	const std::array<double, 3>& target, double tolerance)
{
	std::array<double, 3> parameter{};
	for (int direction = 0; direction < 3; ++direction) {
		const auto width = bounds.maximum[direction]-bounds.minimum[direction];
		parameter[direction] = width > tolerance
			? (target[direction]-bounds.minimum[direction])/width : 0.5;
		parameter[direction] = std::max(0.0, std::min(1.0, parameter[direction]));
	}
	for (int iteration = 0; iteration < 24; ++iteration) {
		std::array<std::array<double, 3>, 3> jacobian{};
		const auto point = MapBezier(element, parameter, &jacobian);
		std::array<double, 3> residual{};
		double residual_squared = 0.0;
		for (int direction = 0; direction < 3; ++direction) {
			residual[direction] = point[direction]-target[direction];
			residual_squared += residual[direction]*residual[direction];
		}
		if (residual_squared <= tolerance*tolerance) {
			constexpr double interior_tolerance = 1e-7;
			return parameter[0] > interior_tolerance && parameter[0] < 1.0-interior_tolerance
				&& parameter[1] > interior_tolerance && parameter[1] < 1.0-interior_tolerance
				&& parameter[2] > interior_tolerance && parameter[2] < 1.0-interior_tolerance;
		}
		std::array<double, 3> update{};
		if (!Solve3x3(jacobian, residual, update)) return false;
		for (int direction = 0; direction < 3; ++direction)
			parameter[direction] -= update[direction];
		if (parameter[0] < -0.25 || parameter[0] > 1.25
			|| parameter[1] < -0.25 || parameter[1] > 1.25
			|| parameter[2] < -0.25 || parameter[2] > 1.25)
			return false;
	}
	return false;
}

inline bool SampledVolumesOverlap(const Element& left, const Bounds& left_bounds,
	const Element& right, const Bounds& right_bounds, double tolerance)
{
	constexpr std::array<double, 3> samples{{0.2, 0.5, 0.8}};
	for (const auto u : samples)
		for (const auto v : samples)
			for (const auto w : samples) {
				const std::array<double, 3> parameter{{u, v, w}};
				if (StrictlyInside(right, right_bounds,
					MapBezier(left, parameter), tolerance)) return true;
				if (StrictlyInside(left, left_bounds,
					MapBezier(right, parameter), tolerance)) return true;
			}
	return false;
}

inline std::size_t SharedPointCount(const std::array<std::int64_t, 64>& left,
	const std::array<std::int64_t, 64>& right)
{
	std::size_t result = 0;
	for (const auto left_point : left)
		if (std::find(right.begin(), right.end(), left_point) != right.end()) ++result;
	return result;
}

inline constexpr std::array<int, 64> VtkCubicHexTensorIndices()
{
	return {{
		0, 3, 15, 12, 48, 51, 63, 60,
		1, 2, 7, 11, 13, 14, 4, 8,
		49, 50, 55, 59, 61, 62, 52, 56,
		16, 32, 19, 35, 31, 47, 28, 44,
		20, 36, 24, 40, 23, 39, 27, 43,
		17, 18, 33, 34, 29, 30, 45, 46,
		5, 6, 9, 10, 53, 54, 57, 58,
		21, 22, 25, 26, 37, 38, 41, 42}};
}

} // namespace detail

inline BezierVisualizationMesh BuildBezierVisualizationMesh(
	Database& database, bool require_valid_geometry = true)
{
	BezierVisualizationMesh mesh;
	const auto elements = database.header().elements;
	mesh.validation.elements = elements;
	mesh.validation.local_points = elements*kBezierPointCount;
	mesh.connectivity.reserve(static_cast<std::size_t>(elements)*kBezierPointCount);
	mesh.offsets.reserve(static_cast<std::size_t>(elements)+1);
	mesh.offsets.push_back(0);
	mesh.types.reserve(static_cast<std::size_t>(elements));
	mesh.higher_order_degrees.reserve(static_cast<std::size_t>(elements));
	mesh.element_ids.reserve(static_cast<std::size_t>(elements));
	mesh.element_owners.reserve(static_cast<std::size_t>(elements));
	std::vector<detail::Bounds> element_bounds(static_cast<std::size_t>(elements));
	std::vector<std::array<std::int64_t, 64>> tensor_connectivity(
		static_cast<std::size_t>(elements));
	detail::Bounds global_bounds;
	for (std::uint64_t index = 0; index < elements; ++index) {
		const auto element = database.Load(index);
		for (const auto& point : element.bezier_points)
			for (const auto coordinate : point)
				if (!std::isfinite(coordinate))
					throw std::runtime_error("Bezier geometry contains a non-finite coordinate");
		element_bounds[static_cast<std::size_t>(index)] = detail::ElementBounds(element);
		for (int direction = 0; direction < 3; ++direction) {
			global_bounds.minimum[direction] = std::min(global_bounds.minimum[direction],
				element_bounds[static_cast<std::size_t>(index)].minimum[direction]);
			global_bounds.maximum[direction] = std::max(global_bounds.maximum[direction],
				element_bounds[static_cast<std::size_t>(index)].maximum[direction]);
		}
	}
	double diagonal_squared = 0.0;
	double maximum_coordinate = 0.0;
	for (int direction = 0; direction < 3; ++direction) {
		const auto width = global_bounds.maximum[direction]-global_bounds.minimum[direction];
		diagonal_squared += width*width;
		maximum_coordinate = std::max({maximum_coordinate,
			std::abs(global_bounds.minimum[direction]),
			std::abs(global_bounds.maximum[direction])});
	}
	const auto geometry_tolerance = 1e-11
		*std::max(std::sqrt(diagonal_squared), 1e-12);
	// The preprocessing cache intentionally preserves the legacy %.6g values.
	// Permit the final displayed digit to differ across two element-local sums.
	const auto point_match_tolerance = std::max(geometry_tolerance,
		2e-6*std::max(1.0, maximum_coordinate));
	const auto point_match_tolerance_squared
		= point_match_tolerance*point_match_tolerance;
	std::unordered_map<std::uint64_t, std::vector<std::int64_t>> registry;
	std::unordered_map<detail::CoordinateKey, std::vector<std::int64_t>,
		detail::CoordinateKeyHash> coordinate_registry;
	std::vector<std::uint64_t> point_elements;
	constexpr double coefficient_quantization = 1e12;
	constexpr std::array<double, 4> quadrature{{
		0.06943184420297371, 0.33000947820757187,
		0.6699905217924281, 0.9305681557970262}};
	const auto vtk_order = detail::VtkCubicHexTensorIndices();
	for (std::uint64_t index = 0; index < elements; ++index) {
		const auto element = database.Load(index);
		auto& local_ids = tensor_connectivity[static_cast<std::size_t>(index)];
		for (std::size_t point = 0; point < kBezierPointCount; ++point) {
			detail::QuantizedSignature key;
			std::vector<std::pair<std::int32_t, double>> signature;
			for (std::size_t row = 0; row < element.extraction.size(); ++row) {
				const auto coefficient = element.extraction[row][point];
				if (coefficient == 0.0) continue;
				const auto node = element.connectivity[row];
				key.push_back({node, detail::QuantizeCoefficient(
					coefficient, coefficient_quantization)});
				signature.push_back({node, coefficient});
			}
			if (key.empty())
				throw std::runtime_error("Bezier extraction contains an empty column");
			std::sort(key.begin(), key.end());
			std::sort(signature.begin(), signature.end(),
				[](const auto& left, const auto& right) { return left.first < right.first; });
			const auto signature_hash = detail::SignatureHash(key);
			const auto found = registry.find(signature_hash);
			std::int64_t matching_point = -1;
			if (found != registry.end())
				for (const auto candidate : found->second)
					if (detail::SignatureMatches(mesh,
						static_cast<std::size_t>(candidate), key,
						coefficient_quantization)) {
						matching_point = candidate;
						break;
					}
			if (matching_point >= 0) {
				const auto global = matching_point;
				const auto distance_squared = detail::DistanceSquared(
					mesh.points[static_cast<std::size_t>(global)], element.bezier_points[point]);
				if (distance_squared > point_match_tolerance_squared) {
					std::ostringstream message;
					message << std::setprecision(17)
						<< "matching Bezier extraction signatures have inconsistent coordinates: first_element="
						<< point_elements[static_cast<std::size_t>(global)]
						<< " current_element=" << element.id
						<< " current_local_point=" << point
						<< " first_coordinate=["
						<< mesh.points[static_cast<std::size_t>(global)][0] << ','
						<< mesh.points[static_cast<std::size_t>(global)][1] << ','
						<< mesh.points[static_cast<std::size_t>(global)][2] << ']'
						<< " current_coordinate=[" << element.bezier_points[point][0] << ','
						<< element.bezier_points[point][1] << ','
						<< element.bezier_points[point][2] << ']'
						<< " distance=" << std::sqrt(distance_squared);
					throw std::runtime_error(message.str());
				}
				if (distance_squared > geometry_tolerance*geometry_tolerance) {
					++mesh.validation.signature_coordinate_repairs;
					mesh.validation.maximum_signature_coordinate_difference = std::max(
						mesh.validation.maximum_signature_coordinate_difference,
						std::sqrt(distance_squared));
				}
				local_ids[point] = global;
				++mesh.validation.shared_point_references;
				continue;
			}
			const auto bucket = detail::CoordinateBucket(
				element.bezier_points[point], point_match_tolerance);
			std::int64_t coincident = -1;
			for (int dz = -1; dz <= 1 && coincident < 0; ++dz)
				for (int dy = -1; dy <= 1 && coincident < 0; ++dy)
					for (int dx = -1; dx <= 1 && coincident < 0; ++dx) {
						auto neighbor = bucket;
						neighbor.value[0] += dx;
						neighbor.value[1] += dy;
						neighbor.value[2] += dz;
						const auto same_bucket = coordinate_registry.find(neighbor);
						if (same_bucket == coordinate_registry.end()) continue;
						for (const auto candidate : same_bucket->second)
							if (detail::DistanceSquared(
								mesh.points[static_cast<std::size_t>(candidate)],
								element.bezier_points[point]) <= point_match_tolerance_squared) {
								coincident = candidate;
								break;
							}
					}
			if (coincident >= 0) {
				++mesh.validation.coincident_unmerged_points;
				mesh.validation.first_bad_element = std::min(
					mesh.validation.first_bad_element, element.id);
				if (mesh.validation.first_coincident_element_pair[0]
					== std::numeric_limits<std::uint64_t>::max())
					mesh.validation.first_coincident_element_pair = {{
						point_elements[static_cast<std::size_t>(coincident)], element.id}};
			}
			const auto global = static_cast<std::int64_t>(mesh.points.size());
			registry[signature_hash].push_back(global);
			coordinate_registry[bucket].push_back(global);
			mesh.points.push_back(element.bezier_points[point]);
			point_elements.push_back(element.id);
			for (const auto& entry : signature) {
				mesh.signature_nodes.push_back(entry.first);
				mesh.signature_coefficients.push_back(entry.second);
			}
			mesh.signature_offsets.push_back(mesh.signature_nodes.size());
			local_ids[point] = global;
		}
		auto sorted = local_ids;
		std::sort(sorted.begin(), sorted.end());
		const auto unique_end = std::unique(sorted.begin(), sorted.end());
		if (unique_end != sorted.end()) {
			++mesh.validation.collapsed_element_points;
			mesh.validation.first_bad_element = std::min(
				mesh.validation.first_bad_element, element.id);
		}
		auto visualization_element = element;
		for (std::size_t point = 0; point < kBezierPointCount; ++point)
			visualization_element.bezier_points[point]
				= mesh.points[static_cast<std::size_t>(local_ids[point])];
		element_bounds[static_cast<std::size_t>(index)]
			= detail::ElementBounds(visualization_element);
		bool bad_jacobian = false;
		for (const auto u : quadrature)
			for (const auto v : quadrature)
				for (const auto w : quadrature) {
					std::array<std::array<double, 3>, 3> jacobian{};
					detail::MapBezier(visualization_element, {{u, v, w}}, &jacobian);
					const auto determinant = detail::Determinant(jacobian);
					mesh.validation.minimum_jacobian = std::min(
						mesh.validation.minimum_jacobian, determinant);
					if (!(determinant > 0.0)) bad_jacobian = true;
				}
		if (bad_jacobian) {
			++mesh.validation.nonpositive_jacobian_elements;
			mesh.validation.first_bad_element = std::min(
				mesh.validation.first_bad_element, element.id);
		}
		for (const auto tensor : vtk_order)
			mesh.connectivity.push_back(local_ids[static_cast<std::size_t>(tensor)]);
		mesh.offsets.push_back(static_cast<std::int64_t>(mesh.connectivity.size()));
		mesh.types.push_back(kVtkBezierHexahedron);
		mesh.higher_order_degrees.push_back({{3, 3, 3}});
		mesh.element_ids.push_back(element.id);
		mesh.element_owners.push_back(element.owner);
	}
	mesh.validation.unique_points = mesh.points.size();
	std::vector<std::uint64_t> order(static_cast<std::size_t>(elements));
	std::iota(order.begin(), order.end(), std::uint64_t{0});
	int sweep_axis = 0;
	for (int direction = 1; direction < 3; ++direction)
		if (global_bounds.maximum[direction]-global_bounds.minimum[direction]
			> global_bounds.maximum[sweep_axis]-global_bounds.minimum[sweep_axis])
			sweep_axis = direction;
	std::sort(order.begin(), order.end(), [&](const auto left, const auto right) {
		return element_bounds[static_cast<std::size_t>(left)].minimum[sweep_axis]
			< element_bounds[static_cast<std::size_t>(right)].minimum[sweep_axis];
	});
	for (std::size_t position = 0; position < order.size(); ++position) {
		const auto left_index = order[position];
		const auto& left_bounds = element_bounds[static_cast<std::size_t>(left_index)];
		for (std::size_t next = position+1; next < order.size(); ++next) {
			const auto right_index = order[next];
			const auto& right_bounds = element_bounds[static_cast<std::size_t>(right_index)];
			if (right_bounds.minimum[sweep_axis]
				>= left_bounds.maximum[sweep_axis]-geometry_tolerance) break;
			if (!detail::BoundsOverlap(left_bounds, right_bounds, geometry_tolerance)) continue;
			if (detail::SharedPointCount(
				tensor_connectivity[static_cast<std::size_t>(left_index)],
				tensor_connectivity[static_cast<std::size_t>(right_index)]) >= 16) continue;
			++mesh.validation.overlap_candidates;
			auto left = database.Load(left_index);
			auto right = database.Load(right_index);
			for (std::size_t point = 0; point < kBezierPointCount; ++point) {
				left.bezier_points[point] = mesh.points[static_cast<std::size_t>(
					tensor_connectivity[static_cast<std::size_t>(left_index)][point])];
				right.bezier_points[point] = mesh.points[static_cast<std::size_t>(
					tensor_connectivity[static_cast<std::size_t>(right_index)][point])];
			}
			if (detail::SampledVolumesOverlap(left, left_bounds, right, right_bounds,
				geometry_tolerance)) {
				++mesh.validation.overlapping_element_pairs;
				if (mesh.validation.first_overlap_pair[0]
					== std::numeric_limits<std::uint64_t>::max())
					mesh.validation.first_overlap_pair = {{left.id, right.id}};
			}
		}
	}
	if (require_valid_geometry) RequireValidBezierGeometry(mesh.validation);
	return mesh;
}

inline void TransformBezierVisualizationToSourceCoordinates(
	BezierVisualizationMesh& mesh, const GeometryTransform& transform)
{
	const auto scale = transform.source_units_per_normalized_unit;
	if (!(scale > 0.0) || !std::isfinite(scale)
		|| !std::all_of(transform.source_origin.begin(), transform.source_origin.end(),
			[](double value) { return std::isfinite(value); }))
		throw std::runtime_error("invalid Bezier visualization geometry transform");
	for (auto& point : mesh.points)
		for (int direction = 0; direction < 3; ++direction)
			point[direction] = transform.source_origin[direction]+scale*point[direction];
	mesh.validation.maximum_signature_coordinate_difference *= scale;
	if (std::isfinite(mesh.validation.minimum_jacobian))
		mesh.validation.minimum_jacobian *= scale*scale*scale;
}

inline BezierVisualizationMesh BuildSourceCoordinateBezierVisualizationMesh(
	Database& database, bool require_valid_geometry = true)
{
	auto mesh = BuildBezierVisualizationMesh(database, false);
	TransformBezierVisualizationToSourceCoordinates(
		mesh, database.header().geometry_transform);
	if (require_valid_geometry) RequireValidBezierGeometry(mesh.validation);
	return mesh;
}

inline std::vector<VtkPointArray> ExtractBezierPointArrays(
	const BezierVisualizationMesh& mesh,
	const std::vector<VtkPointArray>& control_arrays)
{
	std::vector<VtkPointArray> result;
	result.reserve(control_arrays.size());
	for (const auto& control : control_arrays) {
		if (control.components < 1)
			throw std::runtime_error("Bezier point array has no components: "+control.name);
		VtkPointArray output{control.name, control.components,
			std::vector<double>(mesh.points.size()*static_cast<std::size_t>(control.components), 0.0)};
		for (std::size_t point = 0; point < mesh.points.size(); ++point)
			for (auto entry = mesh.signature_offsets[point];
				entry < mesh.signature_offsets[point+1]; ++entry) {
				const auto node = static_cast<std::size_t>(mesh.signature_nodes[entry]);
				const auto coefficient = mesh.signature_coefficients[entry];
				const auto required = (node+1)*static_cast<std::size_t>(control.components);
				if (control.values.size() < required)
					throw std::runtime_error("control-point array is shorter than Bezier extraction connectivity: "+control.name);
				for (int component = 0; component < control.components; ++component)
					output.values[point*static_cast<std::size_t>(control.components)
						+static_cast<std::size_t>(component)] += coefficient
						*control.values[node*static_cast<std::size_t>(control.components)
							+static_cast<std::size_t>(component)];
			}
		result.push_back(std::move(output));
	}
	return result;
}

} // namespace iga

#endif

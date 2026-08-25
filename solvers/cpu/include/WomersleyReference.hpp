#ifndef IGA_WOMERSLEY_REFERENCE_HPP
#define IGA_WOMERSLEY_REFERENCE_HPP

#include "CaseConfig.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct WomersleyReferenceConfiguration {
	int schema_version = 1;
	std::array<double, 3> axis_origin{{0.0, 0.0, 0.0}};
	std::array<double, 3> axis_direction{{0.0, 0.0, 1.0}};
	double radius = 0.0;
	double density = 0.0;
	double dynamic_viscosity = 0.0;
	double period = 0.0;
	double mean_pressure_gradient = 0.0;
	std::vector<double> cosine_pressure_gradient;
	std::vector<double> sine_pressure_gradient;
	std::vector<double> sample_times;
	std::string file_prefix = "womersley";
};

inline std::array<double, 3> ParseWomersleyVector(
	const config_detail::JsonValue& value, const std::string& context)
{
	const auto& array = config_detail::RequireArray(value, context);
	if (array.size() != 3)
		throw std::runtime_error("womersley reference: "+context+" requires three values");
	return {{config_detail::RequireNumber(array[0], context+"[0]"),
		config_detail::RequireNumber(array[1], context+"[1]"),
		config_detail::RequireNumber(array[2], context+"[2]")}};
}

inline std::vector<double> ParseWomersleyArray(
	const config_detail::JsonValue& value, const std::string& context)
{
	const auto& array = config_detail::RequireArray(value, context);
	std::vector<double> result;
	result.reserve(array.size());
	for (std::size_t i = 0; i < array.size(); ++i)
		result.push_back(config_detail::RequireNumber(
			array[i], context+"["+std::to_string(i)+"]"));
	return result;
}

inline WomersleyReferenceConfiguration ParseWomersleyReferenceConfiguration(
	const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "womersley reference");
	RequireKnownKeys(root, {"schema_version", "axis_origin", "axis_direction",
		"radius", "density", "dynamic_viscosity", "period",
		"mean_pressure_gradient", "cosine_pressure_gradient",
		"sine_pressure_gradient", "sample_times", "file_prefix"},
		"womersley reference");
	auto required = [&](const std::string& key) -> const JsonValue& {
		const auto* value = Find(root, key);
		if (!value) throw std::runtime_error("womersley reference requires '"+key+"'");
		return *value;
	};
	WomersleyReferenceConfiguration result;
	result.schema_version = RequireInteger(required("schema_version"), "schema_version");
	result.axis_origin = ParseWomersleyVector(required("axis_origin"), "axis_origin");
	result.axis_direction = ParseWomersleyVector(
		required("axis_direction"), "axis_direction");
	result.radius = RequireNumber(required("radius"), "radius");
	result.density = RequireNumber(required("density"), "density");
	result.dynamic_viscosity = RequireNumber(
		required("dynamic_viscosity"), "dynamic_viscosity");
	result.period = RequireNumber(required("period"), "period");
	result.mean_pressure_gradient = RequireNumber(
		required("mean_pressure_gradient"), "mean_pressure_gradient");
	result.cosine_pressure_gradient = ParseWomersleyArray(
		required("cosine_pressure_gradient"), "cosine_pressure_gradient");
	result.sine_pressure_gradient = ParseWomersleyArray(
		required("sine_pressure_gradient"), "sine_pressure_gradient");
	result.sample_times = ParseWomersleyArray(required("sample_times"), "sample_times");
	if (const auto* prefix = Find(root, "file_prefix"))
		result.file_prefix = RequireString(*prefix, "file_prefix");
	const double axis_norm = std::sqrt(
		result.axis_direction[0]*result.axis_direction[0]
		+result.axis_direction[1]*result.axis_direction[1]
		+result.axis_direction[2]*result.axis_direction[2]);
	if (result.schema_version != 1 || !(axis_norm > 0.0)
		|| !(result.radius > 0.0) || !(result.density > 0.0)
		|| !(result.dynamic_viscosity > 0.0) || !(result.period > 0.0)
		|| result.sample_times.empty() || result.file_prefix.empty())
		throw std::runtime_error("womersley reference configuration contains invalid values");
	for (auto& component : result.axis_direction) component /= axis_norm;
	const auto harmonics = std::max(result.cosine_pressure_gradient.size(),
		result.sine_pressure_gradient.size());
	result.cosine_pressure_gradient.resize(harmonics, 0.0);
	result.sine_pressure_gradient.resize(harmonics, 0.0);
	for (std::size_t i = 0; i < result.sample_times.size(); ++i) {
		if (result.sample_times[i] < 0.0
			|| (i > 0 && !(result.sample_times[i] > result.sample_times[i-1])))
			throw std::runtime_error(
				"womersley sample_times must be nonnegative and strictly increasing");
	}
	const std::filesystem::path prefix(result.file_prefix);
	if (prefix.is_absolute() || prefix.has_parent_path())
		throw std::runtime_error("womersley file_prefix must be a simple relative name");
	return result;
}

inline WomersleyReferenceConfiguration ReadWomersleyReferenceConfiguration(
	const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error(
		"cannot open Womersley reference configuration: "+path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error(
		"cannot read Womersley reference configuration: "+path.string());
	return ParseWomersleyReferenceConfiguration(contents.str());
}

inline std::complex<long double> ComplexBesselJ0(
	const std::complex<long double>& argument)
{
	const auto factor = -argument*argument/4.0L;
	std::complex<long double> term(1.0L, 0.0L);
	std::complex<long double> sum = term;
	for (int order = 1; order <= 10000; ++order) {
		term *= factor/static_cast<long double>(order*order);
		sum += term;
		if (order > std::abs(argument)
			&& std::abs(term) <= 8.0L*std::numeric_limits<long double>::epsilon()
				*std::max(1.0L, std::abs(sum)))
			return sum;
	}
	throw std::runtime_error("complex Bessel J0 series did not converge");
}

inline double WomersleyAxialVelocity(
	const WomersleyReferenceConfiguration& configuration,
	double radial_distance, double physical_time)
{
	if (!std::isfinite(radial_distance) || radial_distance < 0.0)
		throw std::runtime_error("Womersley radial distance must be finite and nonnegative");
	// Legacy VTK control meshes store coordinates to six decimal places.
	const double radial_tolerance = 2e-6*std::max(1.0, configuration.radius);
	if (radial_distance > configuration.radius+radial_tolerance)
		throw std::runtime_error("mesh point radial distance "
			+std::to_string(radial_distance)
			+" lies outside configured Womersley radius "
			+std::to_string(configuration.radius));
	const auto radius_fraction = std::min(radial_distance/configuration.radius, 1.0);
	double velocity = configuration.mean_pressure_gradient
		*(configuration.radius*configuration.radius-radial_distance*radial_distance)
		/(4.0*configuration.dynamic_viscosity);
	constexpr long double pi = 3.141592653589793238462643383279502884L;
	const long double fundamental = 2.0L*pi/configuration.period;
	const std::complex<long double> imaginary(0.0L, 1.0L);
	for (std::size_t index = 0;
		index < configuration.cosine_pressure_gradient.size(); ++index) {
		const long double harmonic = static_cast<long double>(index+1);
		const long double omega = harmonic*fundamental;
		const long double alpha = configuration.radius*std::sqrt(
			omega*configuration.density/configuration.dynamic_viscosity);
		const std::complex<long double> wall_argument(
			alpha/std::sqrt(2.0L), -alpha/std::sqrt(2.0L));
		const auto denominator = ComplexBesselJ0(wall_argument);
		if (std::abs(denominator) <= std::numeric_limits<long double>::min())
			throw std::runtime_error("Womersley Bessel denominator is zero");
		const auto profile = 1.0L-ComplexBesselJ0(
			wall_argument*static_cast<long double>(radius_fraction))/denominator;
		const std::complex<long double> gradient(
			configuration.cosine_pressure_gradient[index],
			-configuration.sine_pressure_gradient[index]);
		const auto amplitude = gradient/(imaginary*omega
			*static_cast<long double>(configuration.density))
			*profile;
		velocity += static_cast<double>(std::real(amplitude*std::exp(
			imaginary*omega*static_cast<long double>(physical_time))));
	}
	return velocity;
}

inline std::array<double, 3> WomersleyVelocity(
	const WomersleyReferenceConfiguration& configuration,
	const std::array<double, 3>& point, double physical_time)
{
	std::array<double, 3> displacement{};
	double axial_coordinate = 0.0;
	for (int component = 0; component < 3; ++component) {
		displacement[component] = point[component]-configuration.axis_origin[component];
		axial_coordinate += displacement[component]*configuration.axis_direction[component];
	}
	double radial_squared = 0.0;
	for (int component = 0; component < 3; ++component) {
		const auto radial = displacement[component]
			-axial_coordinate*configuration.axis_direction[component];
		radial_squared += radial*radial;
	}
	const auto axial_velocity = WomersleyAxialVelocity(
		configuration, std::sqrt(std::max(0.0, radial_squared)), physical_time);
	return {{axial_velocity*configuration.axis_direction[0],
		axial_velocity*configuration.axis_direction[1],
		axial_velocity*configuration.axis_direction[2]}};
}

inline std::vector<std::array<double, 3>> ReadVtkPointCoordinates(
	const std::filesystem::path& path, std::uint64_t expected_nodes)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open VTK mesh: "+path.string());
	std::string token;
	std::uint64_t nodes = 0;
	while (input >> token)
		if (token == "POINTS") {
			std::string type;
			if (!(input >> nodes >> type))
				throw std::runtime_error("invalid VTK POINTS record");
			break;
		}
	if (nodes != expected_nodes)
		throw std::runtime_error("VTK POINTS count does not match database nodes");
	std::vector<std::array<double, 3>> points(static_cast<std::size_t>(nodes));
	for (auto& point : points)
		if (!(input >> point[0] >> point[1] >> point[2]))
			throw std::runtime_error("VTK POINTS array is truncated");
	return points;
}

inline void RescalePointCoordinatesLikeSpline(
	std::vector<std::array<double, 3>>& points)
{
	if (points.empty()) throw std::runtime_error("cannot rescale an empty point set");
	std::array<double, 3> minimum = points.front();
	std::array<double, 3> maximum = points.front();
	for (const auto& point : points)
		for (int component = 0; component < 3; ++component) {
			minimum[component] = std::min(minimum[component], point[component]);
			maximum[component] = std::max(maximum[component], point[component]);
		}
	double minimum_extent = std::numeric_limits<double>::infinity();
	for (int component = 0; component < 3; ++component)
		minimum_extent = std::min(minimum_extent,
			maximum[component]-minimum[component]);
	if (!(minimum_extent > 0.0) || !std::isfinite(minimum_extent))
		throw std::runtime_error("control mesh has a non-positive coordinate extent");
	for (auto& point : points)
		for (int component = 0; component < 3; ++component)
			point[component] = (point[component]-minimum[component])/minimum_extent;
}

} // namespace iga

#endif

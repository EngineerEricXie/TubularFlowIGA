#ifndef IGA_MESH_CONFIG_HPP
#define IGA_MESH_CONFIG_HPP

#include "CaseConfig.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace iga {

struct ThreeDGeometryDefinition
{
	std::string kind;
	std::string file;
	double length_scale_to_m = 1.0;
};

struct MeshSmoothingDefinition
{
	int iterations = 0;
	double bifurcation_ratio = 0.0;
	double noise_ratio = 0.0;
};

struct MeshCenterlineDefinition
{
	double target_spacing = 0.0;
	double max_spacing_over_diameter = 1.0;
	double max_turn_degrees = 12.0;
	double max_diameter_change_fraction = 0.15;
	double maximum_curvature_radius_product = 0.8;
};

struct MeshJunctionDefinition
{
	double max_spacing_over_diameter = 0.25;
	double upstream_clearance_over_diameter = 1.0;
	double downstream_clearance_over_diameter = 1.5;
	double minimum_angle_degrees = 10.0;
	double maximum_radius_ratio = 8.0;
	int optimization_iterations = 4;
};

struct MeshQualityDefinition
{
	double minimum_scaled_jacobian = 0.1;
	bool check_self_intersection = true;
	double collision_safety_factor = 1.0;
};

struct MeshDefinition
{
	MeshSmoothingDefinition smoothing;
	MeshCenterlineDefinition centerline;
	MeshJunctionDefinition junction;
	MeshQualityDefinition quality;
};

struct ThreeDMeshCaseConfiguration
{
	int schema_version = 4;
	std::string dimension = "3d";
	ThreeDGeometryDefinition geometry;
	MeshDefinition mesh;
};

namespace mesh_config_detail {

using config_detail::Find;
using config_detail::JsonValue;
using config_detail::RequireArray;
using config_detail::RequireBoolean;
using config_detail::RequireInteger;
using config_detail::RequireKnownKeys;
using config_detail::RequireNumber;
using config_detail::RequireObject;
using config_detail::RequireString;

inline const JsonValue& Required(
	const std::map<std::string, JsonValue>& object,
	const std::string& key,
	const std::string& context)
{
	const auto* value = Find(object, key);
	if (!value)
		throw std::runtime_error("simulation_config.json: " + context + " requires '" + key + "'");
	return *value;
}

inline double Positive(const JsonValue& value, const std::string& context)
{
	const double result = RequireNumber(value, context);
	if (!(result > 0.0))
		throw std::runtime_error("simulation_config.json: " + context + " must be positive");
	return result;
}

inline double UnitInterval(const JsonValue& value, const std::string& context)
{
	const double result = RequireNumber(value, context);
	if (result < 0.0 || result > 1.0)
		throw std::runtime_error("simulation_config.json: " + context + " must be in [0,1]");
	return result;
}

} // namespace mesh_config_detail

inline ThreeDGeometryDefinition ParseThreeDGeometryDefinition(
	const config_detail::JsonValue& value)
{
	using namespace mesh_config_detail;
	const auto& object = RequireObject(value, "geometry");
	RequireKnownKeys(object, {"kind", "file", "length_scale_to_m"}, "geometry");
	ThreeDGeometryDefinition result;
	result.kind = RequireString(Required(object, "kind", "geometry"), "geometry.kind");
	result.file = RequireString(Required(object, "file", "geometry"), "geometry.file");
	result.length_scale_to_m = Positive(
		Required(object, "length_scale_to_m", "geometry"), "geometry.length_scale_to_m");
	if (result.kind != "swc_network" && result.kind != "obj_network")
		throw std::runtime_error("simulation_config.json: geometry.kind must be 'swc_network' or 'obj_network'");
	const std::filesystem::path path(result.file);
	if (result.file.empty() || path.is_absolute() || path.has_parent_path())
		throw std::runtime_error("simulation_config.json: geometry.file must name a file in the case directory");
	const bool obj = path.extension() == ".obj";
	if ((result.kind == "obj_network") != obj
		|| (!obj && path.extension() != ".swc"))
		throw std::runtime_error("simulation_config.json: geometry.kind and geometry.file extension do not match");
	return result;
}

inline MeshDefinition ParseMeshDefinition(const config_detail::JsonValue& value)
{
	using namespace mesh_config_detail;
	const auto& root = RequireObject(value, "mesh");
	RequireKnownKeys(root, {"smoothing", "centerline", "junction", "quality"}, "mesh");
	MeshDefinition result;

	const auto& smoothing = RequireObject(
		Required(root, "smoothing", "mesh"), "mesh.smoothing");
	RequireKnownKeys(smoothing, {"iterations", "bifurcation_ratio", "noise_ratio"}, "mesh.smoothing");
	result.smoothing.iterations = RequireInteger(
		Required(smoothing, "iterations", "mesh.smoothing"), "mesh.smoothing.iterations");
	result.smoothing.bifurcation_ratio = UnitInterval(
		Required(smoothing, "bifurcation_ratio", "mesh.smoothing"),
		"mesh.smoothing.bifurcation_ratio");
	result.smoothing.noise_ratio = UnitInterval(
		Required(smoothing, "noise_ratio", "mesh.smoothing"), "mesh.smoothing.noise_ratio");

	const auto& centerline = RequireObject(
		Required(root, "centerline", "mesh"), "mesh.centerline");
	RequireKnownKeys(centerline, {"target_spacing", "max_spacing_over_diameter",
		"max_turn_degrees", "max_diameter_change_fraction",
		"maximum_curvature_radius_product"}, "mesh.centerline");
	result.centerline.target_spacing = Positive(
		Required(centerline, "target_spacing", "mesh.centerline"), "mesh.centerline.target_spacing");
	result.centerline.max_spacing_over_diameter = Positive(
		Required(centerline, "max_spacing_over_diameter", "mesh.centerline"),
		"mesh.centerline.max_spacing_over_diameter");
	result.centerline.max_turn_degrees = Positive(
		Required(centerline, "max_turn_degrees", "mesh.centerline"),
		"mesh.centerline.max_turn_degrees");
	if (result.centerline.max_turn_degrees > 90.0)
		throw std::runtime_error("simulation_config.json: mesh.centerline.max_turn_degrees must not exceed 90");
	result.centerline.max_diameter_change_fraction = Positive(
		Required(centerline, "max_diameter_change_fraction", "mesh.centerline"),
		"mesh.centerline.max_diameter_change_fraction");
	if (result.centerline.max_diameter_change_fraction > 1.0)
		throw std::runtime_error(
			"simulation_config.json: mesh.centerline.max_diameter_change_fraction must not exceed 1");
	result.centerline.maximum_curvature_radius_product = Positive(
		Required(centerline, "maximum_curvature_radius_product", "mesh.centerline"),
		"mesh.centerline.maximum_curvature_radius_product");
	if (result.centerline.maximum_curvature_radius_product >= 1.0)
		throw std::runtime_error(
			"simulation_config.json: mesh.centerline.maximum_curvature_radius_product must be below 1");

	const auto& junction = RequireObject(
		Required(root, "junction", "mesh"), "mesh.junction");
	RequireKnownKeys(junction, {"max_spacing_over_diameter",
		"upstream_clearance_over_diameter", "downstream_clearance_over_diameter",
		"minimum_angle_degrees", "maximum_radius_ratio", "optimization_iterations"},
		"mesh.junction");
	result.junction.max_spacing_over_diameter = Positive(
		Required(junction, "max_spacing_over_diameter", "mesh.junction"),
		"mesh.junction.max_spacing_over_diameter");
	result.junction.upstream_clearance_over_diameter = Positive(
		Required(junction, "upstream_clearance_over_diameter", "mesh.junction"),
		"mesh.junction.upstream_clearance_over_diameter");
	result.junction.downstream_clearance_over_diameter = Positive(
		Required(junction, "downstream_clearance_over_diameter", "mesh.junction"),
		"mesh.junction.downstream_clearance_over_diameter");
	result.junction.minimum_angle_degrees = Positive(
		Required(junction, "minimum_angle_degrees", "mesh.junction"),
		"mesh.junction.minimum_angle_degrees");
	if (result.junction.minimum_angle_degrees >= 90.0)
		throw std::runtime_error("simulation_config.json: mesh.junction.minimum_angle_degrees must be below 90");
	result.junction.maximum_radius_ratio = Positive(
		Required(junction, "maximum_radius_ratio", "mesh.junction"),
		"mesh.junction.maximum_radius_ratio");
	if (result.junction.maximum_radius_ratio < 1.0)
		throw std::runtime_error("simulation_config.json: mesh.junction.maximum_radius_ratio must be at least 1");
	result.junction.optimization_iterations = RequireInteger(
		Required(junction, "optimization_iterations", "mesh.junction"),
		"mesh.junction.optimization_iterations");
	if (result.junction.optimization_iterations > 1000)
		throw std::runtime_error(
			"simulation_config.json: mesh.junction.optimization_iterations must not exceed 1000");

	const auto& quality = RequireObject(
		Required(root, "quality", "mesh"), "mesh.quality");
	RequireKnownKeys(quality, {"minimum_scaled_jacobian", "check_self_intersection",
		"collision_safety_factor"}, "mesh.quality");
	result.quality.minimum_scaled_jacobian = Positive(
		Required(quality, "minimum_scaled_jacobian", "mesh.quality"),
		"mesh.quality.minimum_scaled_jacobian");
	if (result.quality.minimum_scaled_jacobian > 1.0)
		throw std::runtime_error("simulation_config.json: mesh.quality.minimum_scaled_jacobian must not exceed 1");
	result.quality.check_self_intersection = RequireBoolean(
		Required(quality, "check_self_intersection", "mesh.quality"),
		"mesh.quality.check_self_intersection");
	result.quality.collision_safety_factor = Positive(
		Required(quality, "collision_safety_factor", "mesh.quality"),
		"mesh.quality.collision_safety_factor");
	return result;
}

inline ThreeDMeshCaseConfiguration ParseThreeDMeshCaseConfiguration(const std::string& text)
{
	using namespace mesh_config_detail;
	const auto value = config_detail::JsonParser(text).Parse();
	const auto& root = RequireObject(value, "root");
	const int version = RequireInteger(Required(root, "schema_version", "root"), "schema_version");
	if (version != 4)
		throw std::runtime_error("simulation_config.json: 3d mesh configuration requires schema_version 4");
	const std::string dimension = RequireString(Required(root, "dimension", "root"), "dimension");
	if (dimension != "3d")
		throw std::runtime_error("simulation_config.json: 3d mesh configuration requires dimension '3d'");
	ThreeDMeshCaseConfiguration result;
	result.schema_version = version;
	result.dimension = dimension;
	result.geometry = ParseThreeDGeometryDefinition(Required(root, "geometry", "root"));
	result.mesh = ParseMeshDefinition(Required(root, "mesh", "root"));
	return result;
}

inline ThreeDMeshCaseConfiguration ReadThreeDMeshCaseConfiguration(
	const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open mesh case configuration: " + path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	return ParseThreeDMeshCaseConfiguration(contents.str());
}

} // namespace iga

#endif

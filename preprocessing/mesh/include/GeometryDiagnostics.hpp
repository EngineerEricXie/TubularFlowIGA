#pragma once

#include "SwcGraph.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace tubular {

struct SegmentDiagnostic
{
	int parent = -1;
	int child = -1;
	double length = 0.0;
	double length_over_diameter = 0.0;
	double diameter_change_fraction = 0.0;
};

struct JunctionDiagnostic
{
	int node = -1;
	double minimum_angle_degrees = 0.0;
	double radius_ratio = 0.0;
	double upstream_clearance_over_diameter = 0.0;
	double minimum_downstream_clearance_over_diameter = 0.0;
};

struct CollisionDiagnostic
{
	int first_parent = -1;
	int first_child = -1;
	int second_parent = -1;
	int second_child = -1;
	double distance = 0.0;
	double required_distance = 0.0;
};

struct GeometryDiagnostics
{
	std::vector<SegmentDiagnostic> segments;
	std::vector<JunctionDiagnostic> junctions;
	std::vector<CollisionDiagnostic> collisions;
	std::vector<double> node_min_length_over_diameter;
	std::vector<double> node_curvature_radius_product;
	std::vector<double> node_junction_min_angle_degrees;
	std::vector<int> node_risk;
	std::vector<std::string> warnings;
	std::vector<std::string> errors;

	bool valid() const { return errors.empty(); }
};

GeometryDiagnostics AnalyzeSkeletonGeometry(
	const SwcGraph& graph,
	const MeshParameters& parameters);

void WriteGeometryDiagnosticsJson(
	const GeometryDiagnostics& diagnostics,
	const SwcGraph& graph,
	const MeshParameters& parameters,
	const std::filesystem::path& path);

void WriteGeometryDiagnosticsVtp(
	const GeometryDiagnostics& diagnostics,
	const SwcGraph& graph,
	const std::filesystem::path& path);

void RequireValidGeometry(const GeometryDiagnostics& diagnostics);

} // namespace tubular

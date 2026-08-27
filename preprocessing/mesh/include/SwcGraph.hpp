#pragma once

#include "Geometry.hpp"

#include <filesystem>
#include <vector>

namespace tubular {

struct MeshParameters
{
	int noise_iterations = 0;
	double bifurcation_smoothing = 0.0;
	double noise_smoothing = 0.0;
	double segment_length = 0.0;
	double max_spacing_over_diameter = 1.0;
	double max_turn_degrees = 12.0;
	double max_diameter_change_fraction = 0.15;
	double maximum_curvature_radius_product = 0.8;
	double bifurcation_refinement = 0.25;
	double upstream_clearance_over_diameter = 1.0;
	double downstream_clearance_over_diameter = 1.5;
	double minimum_bifurcation_angle_degrees = 10.0;
	double maximum_junction_radius_ratio = 8.0;
	int junction_optimization_iterations = 4;
	double minimum_scaled_jacobian = 0.1;
	bool check_self_intersection = true;
	double collision_safety_factor = 1.0;

	static MeshParameters Read(const std::filesystem::path& path);
	void Validate() const;
};

struct SwcNode
{
	int id = -1;
	int type = 2;
	Vec3 position;
	double diameter = 0.0;
	int parent = -1;
	std::vector<int> children;
};

class SwcGraph
{
public:
	static SwcGraph Read(const std::filesystem::path& path);
	void Write(const std::filesystem::path& path) const;
	void WriteNormalized(const std::filesystem::path& path) const;
	void WriteVisualizationVtp(const std::filesystem::path& path) const;
	void RebuildChildren();
	void Validate() const;

	int root() const;
	bool is_branch(int node) const { return nodes.at(node).children.size() == 2; }
	bool is_terminal(int node) const { return nodes.at(node).children.empty(); }
	std::vector<std::vector<int>> Sections() const;

	std::vector<SwcNode> nodes;
};

SwcGraph SmoothSkeleton(const SwcGraph& input, const MeshParameters& parameters);

} // namespace tubular

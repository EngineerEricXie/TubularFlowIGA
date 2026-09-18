#pragma once

#include "Geometry.hpp"
#include <array>
#include <filesystem>
#include <vector>

namespace tubular {

using TemplateFace = std::array<int,4>;
struct CrossSectionTemplates {
    std::vector<Vec3> circle, merge;
    std::vector<TemplateFace> circle_faces, merge_faces;
    std::vector<TemplateFace> branch_bottom, branch_left, branch_right;
    std::vector<Vec3> branch_bottom_points, branch_left_points, branch_right_points;
    std::vector<int> boundary_circle, boundary_merge;
    std::vector<int> bottom_boundary, left_boundary, right_boundary;
    std::vector<std::array<int,4>> axis_projections;
    std::array<std::vector<int>,4> quarter_turn_nodes;
    bool generated = false;
    double target_size = 0.0;
    int core_divisions = 0;
    int outer_layers = 0;
};

CrossSectionTemplates GenerateCircularTemplates(double target_size);
CrossSectionTemplates ReadCrossSectionTemplates(const std::filesystem::path& directory);
void WriteCrossSectionTemplatePreviews(const CrossSectionTemplates& templates,
    const std::filesystem::path& directory);

} // namespace tubular

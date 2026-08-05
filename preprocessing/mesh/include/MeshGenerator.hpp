#pragma once

#include "HexMesh.hpp"
#include "SwcGraph.hpp"

#include <filesystem>
#include <vector>

namespace tubular {

struct ControlMesh
{
	std::vector<Vec3> points;
	std::vector<Hex> elements;
	std::vector<int> labels;
	std::vector<Vec3> velocity;
	QualityResult quality;
};

ControlMesh GenerateControlMesh(
	const SwcGraph& skeleton,
	const MeshParameters& parameters,
	const std::filesystem::path& template_directory,
	double minimum_scaled_jacobian = 1.0e-3);

void WriteControlMeshVtk(const ControlMesh& mesh, const std::filesystem::path& path);
void WriteVelocity(const ControlMesh& mesh, const std::filesystem::path& path);

} // namespace tubular

#include "MeshGenerator.hpp"
#include "SwcGraph.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::filesystem::path SkeletonPath(const std::filesystem::path& directory)
{
	const auto swc = directory/"skeleton_initial.swc";
	const auto obj = directory/"skeleton_initial.obj";
	const bool has_swc = std::filesystem::exists(swc);
	const bool has_obj = std::filesystem::exists(obj);
	if (has_swc == has_obj)
		throw std::runtime_error("case directory must contain exactly one skeleton_initial.swc or skeleton_initial.obj");
	return has_swc ? swc : obj;
}

void PrintMesh(const tubular::ControlMesh& mesh)
{
	std::cout << "points=" << mesh.points.size()
		<< " elements=" << mesh.elements.size()
		<< " min_detJ=" << mesh.quality.minimum_determinant
		<< " min_scaled_J=" << mesh.quality.minimum_scaled_jacobian
		<< " bad_elements=" << mesh.quality.bad_elements << '\n';
}

}

int main(int argc, char** argv)
{
	try {
		if (argc < 2) throw std::runtime_error(
			"usage: tubular_mesh smooth|generate|pipeline ...");
		const std::string command = argv[1];
		if (command == "smooth") {
			if (argc != 5) throw std::runtime_error(
				"usage: tubular_mesh smooth INPUT.swc|INPUT.obj mesh_parameter.txt OUTPUT.swc");
			const auto parameters = tubular::MeshParameters::Read(argv[3]);
			const auto output = tubular::SmoothSkeleton(tubular::SwcGraph::Read(argv[2]), parameters);
			output.Write(argv[4]);
			std::cout << "smoothed_nodes=" << output.nodes.size()
				<< " sections=" << output.Sections().size() << '\n';
			return 0;
		}
		if (command == "generate") {
			if (argc != 8) throw std::runtime_error(
				"usage: tubular_mesh generate SMOOTH.swc mesh_parameter.txt TEMPLATE_DIR controlmesh.vtk initial_velocityfield.txt MIN_SCALED_J");
			const auto parameters = tubular::MeshParameters::Read(argv[3]);
			const auto mesh = tubular::GenerateControlMesh(
				tubular::SwcGraph::Read(argv[2]), parameters, argv[4], std::stod(argv[7]));
			tubular::WriteControlMeshVtk(mesh, argv[5]);
			tubular::WriteVelocity(mesh, argv[6]);
			PrintMesh(mesh);
			return 0;
		}
		if (command == "pipeline") {
			if (argc < 3 || argc > 5) throw std::runtime_error(
				"usage: tubular_mesh pipeline CASE_DIR [TEMPLATE_DIR] [MIN_SCALED_J]");
			const std::filesystem::path directory = argv[2];
			const std::filesystem::path templates = argc >= 4 ? argv[3] : "meshgeneration/template";
			const double minimum_scaled = argc >= 5 ? std::stod(argv[4]) : 1.0e-3;
			const auto parameters = tubular::MeshParameters::Read(directory/"mesh_parameter.txt");
			const auto input = tubular::SwcGraph::Read(SkeletonPath(directory));
			const auto normalized_path = directory/"skeleton_normalized.swc";
			input.WriteNormalized(normalized_path);
			input.WriteVisualizationVtp(directory/"skeleton.vtp");
			const auto smooth = tubular::SmoothSkeleton(
				tubular::SwcGraph::Read(normalized_path), parameters);
			const auto smooth_path = directory/"skeleton_smooth.swc";
			smooth.Write(smooth_path);
			const auto quantized_smooth = tubular::SwcGraph::Read(smooth_path);
			const auto mesh = tubular::GenerateControlMesh(quantized_smooth, parameters, templates, minimum_scaled);
			tubular::WriteControlMeshVtk(mesh, directory/"controlmesh.vtk");
			tubular::WriteVelocity(mesh, directory/"initial_velocityfield.txt");
			std::cout << "smoothed_nodes=" << smooth.nodes.size() << ' ';
			PrintMesh(mesh);
			return 0;
		}
		throw std::runtime_error("unknown command: "+command);
	} catch (const std::exception& error) {
		std::cerr << "tubular_mesh: " << error.what() << '\n';
		return 1;
	}
}

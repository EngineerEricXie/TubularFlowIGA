#include "MeshGenerator.hpp"
#include "SwcGraph.hpp"
#include "GeometryDiagnostics.hpp"
#include "MeshConfig.hpp"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct PipelineInput
{
	std::filesystem::path skeleton;
	tubular::MeshParameters parameters;
	bool modern = false;
};

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

tubular::MeshParameters ConvertParameters(const iga::MeshDefinition& input)
{
	tubular::MeshParameters result;
	result.noise_iterations = input.smoothing.iterations;
	result.bifurcation_smoothing = input.smoothing.bifurcation_ratio;
	result.noise_smoothing = input.smoothing.noise_ratio;
	result.segment_length = input.centerline.target_spacing;
	result.max_spacing_over_diameter = input.centerline.max_spacing_over_diameter;
	result.max_turn_degrees = input.centerline.max_turn_degrees;
	result.max_diameter_change_fraction = input.centerline.max_diameter_change_fraction;
	result.maximum_curvature_radius_product = input.centerline.maximum_curvature_radius_product;
	result.bifurcation_refinement = input.junction.max_spacing_over_diameter;
	result.upstream_clearance_over_diameter = input.junction.upstream_clearance_over_diameter;
	result.downstream_clearance_over_diameter = input.junction.downstream_clearance_over_diameter;
	result.minimum_bifurcation_angle_degrees = input.junction.minimum_angle_degrees;
	result.maximum_junction_radius_ratio = input.junction.maximum_radius_ratio;
	result.junction_optimization_iterations = input.junction.optimization_iterations;
	result.minimum_scaled_jacobian = input.quality.minimum_scaled_jacobian;
	result.check_self_intersection = input.quality.check_self_intersection;
	result.collision_safety_factor = input.quality.collision_safety_factor;
	result.Validate();
	return result;
}

PipelineInput ReadPipelineInput(const std::filesystem::path& directory)
{
	const auto simulation = directory/"simulation_config.json";
	const auto legacy = directory/"mesh_parameter.txt";
	if (std::filesystem::exists(simulation)) {
		std::ifstream input(simulation);
		std::ostringstream contents;
		contents << input.rdbuf();
		const auto root = iga::config_detail::RequireObject(
			iga::config_detail::JsonParser(contents.str()).Parse(), "root");
		const auto* version = iga::config_detail::Find(root, "schema_version");
		if (version && iga::config_detail::RequireInteger(*version, "schema_version") == 4) {
			if (std::filesystem::exists(legacy))
				throw std::runtime_error(
					"schema-v4 simulation_config.json and mesh_parameter.txt cannot both be present");
			const auto configuration = iga::ParseThreeDMeshCaseConfiguration(contents.str());
			const auto skeleton = directory/configuration.geometry.file;
			if (!std::filesystem::is_regular_file(skeleton))
				throw std::runtime_error("configured geometry file does not exist: "+skeleton.string());
			return {skeleton, ConvertParameters(configuration.mesh), true};
		}
	}
	if (!std::filesystem::is_regular_file(legacy))
		throw std::runtime_error(
			"case requires schema-v4 simulation_config.json with geometry/mesh blocks or legacy mesh_parameter.txt");
	return {SkeletonPath(directory), tubular::MeshParameters::Read(legacy), false};
}

void PrintMesh(const tubular::ControlMesh& mesh)
{
	std::cout << "points=" << mesh.points.size()
		<< " elements=" << mesh.elements.size()
		<< " min_detJ=" << mesh.quality.minimum_determinant
		<< " min_scaled_J=" << mesh.quality.minimum_scaled_jacobian
		<< " bad_elements=" << mesh.quality.bad_elements
		<< " surface_intersections=" << mesh.surface_intersections.intersections << '\n';
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
			const double minimum_scaled = std::stod(argv[7]);
			const auto mesh = tubular::GenerateControlMesh(
				tubular::SwcGraph::Read(argv[2]), parameters, argv[4], minimum_scaled);
			tubular::WriteControlMeshVtk(mesh, argv[5]);
			tubular::WriteMeshQualityJson(
				mesh, minimum_scaled, std::filesystem::path(argv[5]).parent_path()/"mesh_quality.json");
			tubular::WriteVelocity(mesh, argv[6]);
			PrintMesh(mesh);
			return 0;
		}
		if (command == "pipeline") {
			const bool allow_preflight_failure = argc >= 3
				&& std::string(argv[argc-1]) == "--allow-preflight-failure";
			const int positional_argc = argc-(allow_preflight_failure ? 1 : 0);
			if (positional_argc < 3 || positional_argc > 5) throw std::runtime_error(
				"usage: tubular_mesh pipeline CASE_DIR [TEMPLATE_DIR] [MIN_SCALED_J] [--allow-preflight-failure]");
			const std::filesystem::path directory = argv[2];
			const std::filesystem::path templates = positional_argc >= 4 ? argv[3] : "meshgeneration/template";
			const double minimum_scaled = positional_argc >= 5 ? std::stod(argv[4]) : 1.0e-3;
			const auto pipeline = ReadPipelineInput(directory);
			auto parameters = pipeline.parameters;
			if (pipeline.modern && positional_argc >= 5)
				throw std::runtime_error("MIN_SCALED_J override is not accepted with schema-v4 mesh.quality");
			if (!pipeline.modern && positional_argc >= 5) parameters.minimum_scaled_jacobian = minimum_scaled;
			const auto input = tubular::SwcGraph::Read(pipeline.skeleton);
			const auto normalized_path = directory/"skeleton_normalized.swc";
			input.WriteNormalized(normalized_path);
			input.WriteVisualizationVtp(directory/"skeleton.vtp");
			const auto smooth = tubular::SmoothSkeleton(
				tubular::SwcGraph::Read(normalized_path), parameters);
			const auto smooth_path = directory/"skeleton_smooth.swc";
			smooth.Write(smooth_path);
			const auto quantized_smooth = tubular::SwcGraph::Read(smooth_path);
			const auto diagnostics = tubular::AnalyzeSkeletonGeometry(quantized_smooth, parameters);
			tubular::WriteGeometryDiagnosticsJson(
				diagnostics, quantized_smooth, parameters, directory/"mesh_diagnostics.json");
			tubular::WriteGeometryDiagnosticsVtp(
				diagnostics, quantized_smooth, directory/"skeleton_diagnostics.vtp");
			std::cout<<"geometry_preflight_errors="<<diagnostics.errors.size()
				<<" warnings="<<diagnostics.warnings.size()
				<<" candidate_collisions="<<diagnostics.collisions.size()<<'\n';
			for(const auto& warning:diagnostics.warnings)
				std::cerr<<"tubular_mesh warning: "<<warning<<'\n';
			if (allow_preflight_failure && !diagnostics.valid())
				std::cerr<<"tubular_mesh warning: continuing after geometry preflight failure because "
					"--allow-preflight-failure was specified\n";
			else
				tubular::RequireValidGeometry(diagnostics);
			const auto mesh = tubular::GenerateControlMesh(
				quantized_smooth, parameters, templates, parameters.minimum_scaled_jacobian);
			tubular::WriteControlMeshVtk(mesh, directory/"controlmesh.vtk");
			tubular::WriteMeshQualityJson(
				mesh, parameters.minimum_scaled_jacobian, directory/"mesh_quality.json");
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

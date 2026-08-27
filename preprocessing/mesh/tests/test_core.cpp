#include "BSpline.hpp"
#include "GeometryDiagnostics.hpp"
#include "HexMesh.hpp"
#include "MeshGenerator.hpp"
#include "SwcGraph.hpp"
#include "MeshConfig.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

template<class Function>
void RequireFailure(Function function, const std::string& expected)
{
	try { function(); }
	catch (const std::exception& error) {
		Require(std::string(error.what()).find(expected) != std::string::npos,
			"failure message did not contain expected text");
		return;
	}
	throw std::runtime_error("expected operation to fail");
}

}

int main(int argc, char** argv)
{
	using namespace tubular;
	if (argc != 2) throw std::runtime_error("usage: mesh_core_test TEMPLATE_DIR");
	Require(Norm(RotateSurface({1,0,0}, {0,0,1}, {0,1,0})-Vec3{1,0,0}) < 1.0e-12,
		"shortest rotation changed an orthogonal vector");

	const auto samples = SampleBranch({{0,0,0},{5,0,0},{10,0,0}}, {1,1,1}, 0.5, 4);
	Require(samples.size() >= 4, "B-spline produced too few samples");
	Require(Norm(samples.front().point-Vec3{0,0,0}) < 1.0e-12, "B-spline start mismatch");
	Require(Norm(samples.back().point-Vec3{10,0,0}) < 1.0e-12, "B-spline end mismatch");
	for (const auto& sample : samples)
		Require(sample.tangent.x > 0.0 && std::abs(sample.tangent.y) < 1.0e-12, "B-spline tangent mismatch");
	BranchSamplingOptions short_options;
	short_options.target_spacing=0.25;
	RequireFailure([&] { SampleBranch({{0,0,0},{1,0,0}}, {1,1}, short_options, 2, "short arm"); },
		"insufficient bifurcation clearance");
	BranchSamplingOptions adaptive_options;
	adaptive_options.target_spacing=10.0;
	adaptive_options.max_spacing_over_diameter=10.0;
	adaptive_options.max_turn_degrees=5.0;
	adaptive_options.max_diameter_change_fraction=0.05;
	const auto adaptive_samples=SampleBranch(
		{{0,0,0},{1,0,0},{1,1,0}}, {1,3,1}, adaptive_options, 4, "adaptive bend");
	Require(adaptive_samples.size()>10,"adaptive bend was not refined");
	for(std::size_t i=1;i<adaptive_samples.size();++i) {
		const auto& first=adaptive_samples[i-1];
		const auto& second=adaptive_samples[i];
		const double angle=std::acos(ClampUnit(Dot(
			Normalized(first.tangent,"test tangent"),Normalized(second.tangent,"test tangent"))))
			*180.0/std::acos(-1.0);
		Require(angle<=5.000001,"adaptive sampling exceeded the tangent-turn limit");
		Require(std::abs(second.diameter-first.diameter)
			/std::min(first.diameter,second.diameter)<=0.050001,
			"adaptive sampling exceeded the diameter-change limit");
	}
	const auto legacy_parameters_path=std::filesystem::temp_directory_path()
		/"tubularflowiga-mesh-parameters.txt";
	{
		std::ofstream output(legacy_parameters_path);
		output<<"n_noisesmooth 0\nratio_bifur_node 0.2\nratio_noisesmooth 0\n"
			<<"seg_length 0.5\nratio_refine 0.25\n";
	}
	Require(std::abs(MeshParameters::Read(legacy_parameters_path).segment_length-0.5)<1.0e-12,
		"valid legacy mesh parameters were rejected");
	{
		std::ofstream output(legacy_parameters_path);
		output<<"misspelled_name 0\nratio_bifur_node 0.2\nratio_noisesmooth 0\n"
			<<"seg_length 0.5\nratio_refine 0.25\n";
	}
	RequireFailure([&] { MeshParameters::Read(legacy_parameters_path); },
		"legacy mesh parameter names or order are invalid");
	std::filesystem::remove(legacy_parameters_path);

	const auto mesh_configuration=iga::ParseThreeDMeshCaseConfiguration(R"json({
		"schema_version":4,"dimension":"3d",
		"geometry":{"kind":"swc_network","file":"skeleton_initial.swc","length_scale_to_m":0.001},
		"mesh":{
			"smoothing":{"iterations":2,"bifurcation_ratio":0.1,"noise_ratio":0.05},
			"centerline":{"target_spacing":0.5,"max_spacing_over_diameter":1.0,
				"max_turn_degrees":12,"max_diameter_change_fraction":0.15,
				"maximum_curvature_radius_product":0.8},
			"junction":{"max_spacing_over_diameter":0.25,
				"upstream_clearance_over_diameter":1.0,"downstream_clearance_over_diameter":1.5,
				"minimum_angle_degrees":10,"maximum_radius_ratio":8,"optimization_iterations":4},
			"quality":{"minimum_scaled_jacobian":0.1,"check_self_intersection":true,
				"collision_safety_factor":1.0}
		}
	})json");
	Require(mesh_configuration.mesh.smoothing.iterations==2,"schema-v4 smoothing parse mismatch");
	Require(mesh_configuration.geometry.file=="skeleton_initial.swc","schema-v4 geometry parse mismatch");

	std::vector<Vec3> cube{{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
		{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
	std::vector<Hex> elements{{0,1,2,3,4,5,6,7}};
	const auto quality = EvaluateHexQuality(cube, elements);
	Require(quality.bad_elements == 0, "unit cube reported invalid");
	Require(std::abs(quality.minimum_determinant-1.0) < 1.0e-12, "unit cube determinant mismatch");
	Require(std::abs(quality.minimum_scaled_jacobian-1.0) < 1.0e-12, "unit cube scaled Jacobian mismatch");
	std::vector<Vec3> small_cube=cube;
	for(auto& point:small_cube) point*=1.0e-6;
	const auto small_quality=EvaluateHexQuality(small_cube,elements);
	Require(small_quality.bad_elements==0
		&&std::abs(small_quality.minimum_scaled_jacobian-1.0)<1.0e-12,
		"scaled Jacobian changed with coordinate units");
	std::vector<Vec3> overlapping=cube;
	auto shifted_cube=cube;
	for(auto& point:shifted_cube) point+=Vec3{1.0,0.25,0.5};
	overlapping.insert(overlapping.end(),shifted_cube.begin(),shifted_cube.end());
	std::vector<Hex> overlapping_elements{{0,1,2,3,4,5,6,7},{8,9,10,11,12,13,14,15}};
	Require(EvaluateBoundarySelfIntersections(overlapping,overlapping_elements).intersections==1,
		"overlapping disconnected cubes were not detected");
	for(auto& point:overlapping) point*=1.0e-12;
	Require(EvaluateBoundarySelfIntersections(overlapping,overlapping_elements).intersections==1,
		"scale-independent surface intersection check failed");

	const auto obj_path = std::filesystem::temp_directory_path()/"tubularflowiga-radius-skeleton.obj";
	{
		std::ofstream output(obj_path);
		output << "v 0 0 0 2 0 0\n"
			<< "v 1 0 0 1 0 0\n"
			<< "v 2 1 0 0.5 0 0\n"
			<< "v 2 -1 0 0.5 0 0\n"
			<< "l 1 2\n"
			<< "l 2 3\n"
			<< "l 2 4\n";
	}
	const auto obj_graph = SwcGraph::Read(obj_path);
	Require(obj_graph.nodes.size() == 4, "OBJ skeleton vertex count mismatch");
	Require(obj_graph.root() == 0, "OBJ skeleton inferred root mismatch");
	Require(obj_graph.nodes[0].diameter == 4.0, "OBJ skeleton radius conversion mismatch");
	Require(obj_graph.nodes[1].children.size() == 2, "OBJ skeleton orientation mismatch");
	const auto normalized_path = std::filesystem::temp_directory_path()/"tubularflowiga-normalized.swc";
	const auto vtp_path = std::filesystem::temp_directory_path()/"tubularflowiga-skeleton.vtp";
	obj_graph.WriteNormalized(normalized_path);
	obj_graph.WriteVisualizationVtp(vtp_path);
	Require(SwcGraph::Read(normalized_path).nodes.size() == 4,
		"normalized OBJ-to-SWC output is invalid");
	{
		std::ifstream input(vtp_path);
		const std::string contents((std::istreambuf_iterator<char>(input)),
			std::istreambuf_iterator<char>());
		Require(contents.find("Name=\"radius\"") != std::string::npos,
			"skeleton VTP is missing radius");
		Require(contents.find("Name=\"branch_id\"") != std::string::npos,
			"skeleton VTP is missing branch_id");
	}
	std::filesystem::remove(obj_path);
	std::filesystem::remove(normalized_path);
	std::filesystem::remove(vtp_path);

	const auto invalid_obj_path = std::filesystem::temp_directory_path()/"tubularflowiga-surface.obj";
	{
		std::ofstream output(invalid_obj_path);
		output << "v 0 0 0 1 0 0\n"
			<< "v 1 0 0 1 0 0\n"
			<< "f 1 2 1\n";
	}
	RequireFailure([&] { SwcGraph::Read(invalid_obj_path); }, "supports only v and l");
	std::filesystem::remove(invalid_obj_path);

	SwcGraph graph;
	graph.nodes.resize(4);
	graph.nodes[0].parent=-1;
	graph.nodes[1].parent=0;
	graph.nodes[2].parent=1;
	graph.nodes[3].parent=1;
	for (int i=0; i<4; ++i) { graph.nodes[i].position={double(i),0,0}; graph.nodes[i].diameter=1.0; }
	graph.nodes[2].position={2,1,0};
	graph.nodes[3].position={2,-1,0};
	graph.RebuildChildren();
	graph.Validate();
	Require(graph.Sections().size() == 3, "branch section extraction mismatch");
	SwcGraph multiway = graph;
	multiway.nodes.push_back({});
	multiway.nodes.back().parent = 1;
	multiway.nodes.back().position = {2,0,1};
	multiway.nodes.back().diameter = 1.0;
	multiway.RebuildChildren();
	RequireFailure([&] { multiway.Validate(); },
		"3d mesh generation supports at most two children; node 2 has 3");

	SwcGraph y;
	y.nodes.resize(5);
	y.nodes[0].parent=-1;y.nodes[0].position={0,0,0};y.nodes[0].diameter=1.0;
	y.nodes[1].parent=0;y.nodes[1].position={2,0,0};y.nodes[1].diameter=1.0;
	y.nodes[2].parent=1;y.nodes[2].position={4,0,0};y.nodes[2].diameter=1.0;
	y.nodes[3].parent=2;y.nodes[3].position={6,2,0};y.nodes[3].diameter=0.8;
	y.nodes[4].parent=2;y.nodes[4].position={6,-2,0};y.nodes[4].diameter=0.8;
	y.RebuildChildren();y.Validate();
	MeshParameters y_parameters;
	y_parameters.segment_length=0.5;
	y_parameters.bifurcation_refinement=0.25;
	const auto y_diagnostics=AnalyzeSkeletonGeometry(y,y_parameters);
	Require(y_diagnostics.valid(),"well-spaced Y bifurcation failed geometry preflight");
	SwcGraph short_y=y;
	short_y.nodes[3].position={4.25,0.1,0};
	short_y.nodes[4].position={4.25,-0.1,0};
	const auto short_diagnostics=AnalyzeSkeletonGeometry(short_y,y_parameters);
	Require(!short_diagnostics.valid(),"short thick Y bifurcation passed geometry preflight");
	Require(short_diagnostics.errors.front().find("bifurcation node")!=std::string::npos,
		"pathological Y diagnostic omitted the node context");
	SwcGraph narrow_y=y;
	narrow_y.nodes[3].position={6.0,0.10,0};
	narrow_y.nodes[4].position={6.0,-0.10,0};
	Require(!AnalyzeSkeletonGeometry(narrow_y,y_parameters).valid(),
		"near-collinear Y bifurcation passed geometry preflight");
	SwcGraph asymmetric_y=y;
	asymmetric_y.nodes[4].diameter=0.05;
	Require(!AnalyzeSkeletonGeometry(asymmetric_y,y_parameters).valid(),
		"extreme junction radius ratio passed geometry preflight");
	SwcGraph tight_bend;
	tight_bend.nodes.resize(3);
	tight_bend.nodes[0].parent=-1;tight_bend.nodes[0].position={0,0,0};
	tight_bend.nodes[1].parent=0;tight_bend.nodes[1].position={1,0,0};
	tight_bend.nodes[2].parent=1;tight_bend.nodes[2].position={1,0.1,0};
	for(auto& node:tight_bend.nodes) node.diameter=2.0;
	tight_bend.RebuildChildren();tight_bend.Validate();
	Require(!AnalyzeSkeletonGeometry(tight_bend,y_parameters).valid(),
		"curvature-radius singularity passed geometry preflight");
	SwcGraph crossing;
	crossing.nodes.resize(4);
	crossing.nodes[0].parent=-1;crossing.nodes[0].position={-2,0,0};
	crossing.nodes[1].parent=0;crossing.nodes[1].position={2,0,0};
	crossing.nodes[2].parent=1;crossing.nodes[2].position={2,2,0};
	crossing.nodes[3].parent=2;crossing.nodes[3].position={-2,-2,0};
	for(auto& node:crossing.nodes) node.diameter=0.5;
	crossing.RebuildChildren();crossing.Validate();
	const auto crossing_diagnostics=AnalyzeSkeletonGeometry(crossing,y_parameters);
	Require(!crossing_diagnostics.collisions.empty(),"crossing centerline tubes were not diagnosed");
	Require(!crossing_diagnostics.warnings.empty(),"crossing centerline omitted its warning");
	const auto y_mesh=GenerateControlMesh(y,y_parameters,std::filesystem::path(argv[1]));
	Require(y_mesh.quality.bad_elements==0,"Y integration generated invalid elements");
	Require(y_mesh.surface_intersections.intersections==0,"Y integration self-intersected");
	const auto quality_path=std::filesystem::temp_directory_path()/"tubularflowiga-mesh-quality.json";
	WriteMeshQualityJson(y_mesh,y_parameters.minimum_scaled_jacobian,quality_path);
	{
		std::ifstream input(quality_path);
		const std::string contents((std::istreambuf_iterator<char>(input)),
			std::istreambuf_iterator<char>());
		Require(contents.find("\"surface_intersections\": 0")!=std::string::npos,
			"mesh quality report omitted surface-intersection result");
	}
	std::filesystem::remove(quality_path);

	SwcGraph pipe;
	pipe.nodes.resize(3);
	pipe.nodes[0].parent=-1; pipe.nodes[0].position={0,0,0}; pipe.nodes[0].diameter=1.0;
	pipe.nodes[1].parent=0; pipe.nodes[1].position={5,0,0}; pipe.nodes[1].diameter=1.0;
	pipe.nodes[2].parent=1; pipe.nodes[2].position={10,0,0}; pipe.nodes[2].diameter=1.0;
	pipe.RebuildChildren();
	MeshParameters parameters;
	parameters.segment_length=0.5;
	const auto mesh=GenerateControlMesh(pipe,parameters,std::filesystem::path(argv[1]));
	Require(mesh.points.size()==2211, "pipe integration point count mismatch");
	Require(mesh.elements.size()==1800, "pipe integration element count mismatch");
	Require(mesh.quality.bad_elements==0, "pipe integration generated invalid elements");

	std::cout << "mesh_core_test: PASS\n";
	return 0;
}

#include "BSpline.hpp"
#include "HexMesh.hpp"
#include "MeshGenerator.hpp"
#include "SwcGraph.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
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

	std::vector<Vec3> cube{{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
		{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
	std::vector<Hex> elements{{0,1,2,3,4,5,6,7}};
	const auto quality = EvaluateHexQuality(cube, elements);
	Require(quality.bad_elements == 0, "unit cube reported invalid");
	Require(std::abs(quality.minimum_determinant-1.0) < 1.0e-12, "unit cube determinant mismatch");
	Require(std::abs(quality.minimum_scaled_jacobian-1.0) < 1.0e-12, "unit cube scaled Jacobian mismatch");

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

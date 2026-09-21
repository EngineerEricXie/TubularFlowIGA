#include "NativeTetMeshComponents.hpp"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

template <class Function>
void Reject(Function&& function, const char* message)
{
	bool rejected = false;
	try {
		function();
	} catch (const std::invalid_argument& error) {
		rejected = std::string(error.what()).find(message) != std::string::npos;
	}
	assert(rejected);
}

} // namespace

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points.resize(8);
	mesh.cells.push_back({1, {{0, 1, 2, 3}}});
	assert(iga::CountNativeTetFaceComponents(mesh) == 1);
	mesh.cells.push_back({2, {{0, 1, 2, 4}}});
	assert(iga::CountNativeTetFaceComponents(mesh) == 1);
	mesh.cells[1].nodes = {{0, 4, 5, 6}};
	assert(iga::CountNativeTetFaceComponents(mesh) == 2);
	Reject([&] { iga::RequireNativeTetSingleFaceComponent(mesh, "fluid"); },
		"2 disconnected face components");
	mesh.cells[1].nodes = {{4, 5, 6, 7}};
	assert(iga::CountNativeTetFaceComponents(mesh) == 2);
	mesh.cells[1].nodes = {{0, 0, 2, 4}};
	Reject([&] { iga::CountNativeTetFaceComponents(mesh); }, "repeated node");
	mesh.cells[1].nodes = {{0, 1, 2, 8}};
	Reject([&] { iga::CountNativeTetFaceComponents(mesh); }, "out of range");
	mesh.cells[1].nodes = {{0, 1, 2, 4}};
	mesh.cells[1].nodes = {{3, 2, 1, 0}};
	Reject([&] { iga::CountNativeTetFaceComponents(mesh); }, "geometry is duplicated");
	mesh.cells[1].nodes = {{0, 1, 2, 4}};
	mesh.cells[1].id = 1;
	Reject([&] { iga::CountNativeTetFaceComponents(mesh); }, "cell ID is duplicated");
	mesh.cells[1].id = 2;
	mesh.cells.push_back({3, {{0, 1, 2, 5}}});
	Reject([&] { iga::CountNativeTetFaceComponents(mesh); }, "more than two owners");
	std::cout << "native_tet_mesh_components_test: PASS\n";
}

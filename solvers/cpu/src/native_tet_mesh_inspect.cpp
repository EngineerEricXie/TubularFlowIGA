#include "NativeTetFem.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

int main(int argc, char** argv)
{
	try {
		if (argc != 2) throw std::runtime_error("usage: native_tet_mesh_inspect mesh.msh");
		std::ifstream input(argv[1]);
		if (!input) throw std::runtime_error("cannot open native tetra mesh");
		const auto mesh = iga::ReadNativeTetMeshGmsh41(input);
		const auto topology = iga::BuildNativeTaylorHoodTopology(mesh);
		for (const auto& cell : mesh.cells) iga::EvaluateNativeTetGeometry(mesh, cell);
		std::map<int, std::size_t> labels;
		for (const auto& face : mesh.boundary_triangles) ++labels[face.boundary_label];
		std::cout << "native_tet_mesh_inspect: PASS nodes=" << mesh.points.size()
			<< " tetrahedra=" << mesh.cells.size()
			<< " p2_edges=" << topology.edges.size()
			<< " boundary_triangles=" << mesh.boundary_triangles.size()
			<< " labels=";
		bool first = true;
		for (const auto& item : labels) {
			if (!first) std::cout << ',';
			std::cout << item.first << ':' << item.second;
			first = false;
		}
		std::cout << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "native_tet_mesh_inspect: ERROR: " << error.what() << '\n';
		return 2;
	}
}

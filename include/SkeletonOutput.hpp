#ifndef IGA_SKELETON_OUTPUT_HPP
#define IGA_SKELETON_OUTPUT_HPP

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace iga {

struct SkeletonOutputNode {
	int id = -1;
	int type = 2;
	std::array<double, 3> position{{0.0, 0.0, 0.0}};
	double radius = 0.0;
	int parent_id = -1;
};

struct SkeletonOutputTopology {
	std::vector<int> parent;
	std::vector<std::vector<int>> children;
	std::vector<int> incoming_branch;
	int root = -1;
	int branches = 0;
};

inline SkeletonOutputTopology AnalyzeSkeletonOutput(
	const std::vector<SkeletonOutputNode>& nodes)
{
	if (nodes.size() < 2) throw std::runtime_error("skeleton output requires at least two nodes");
	SkeletonOutputTopology topology;
	topology.parent.assign(nodes.size(), -1);
	topology.children.resize(nodes.size());
	topology.incoming_branch.assign(nodes.size(), -1);
	std::unordered_map<int, int> index;
	for (std::size_t i = 0; i < nodes.size(); ++i) {
		const auto& node = nodes[i];
		if (node.id < 1 || (node.parent_id != -1 && node.parent_id < 1)
			|| !std::isfinite(node.position[0]) || !std::isfinite(node.position[1])
			|| !std::isfinite(node.position[2]) || !std::isfinite(node.radius)
			|| !(node.radius > 0.0))
			throw std::runtime_error("skeleton output requires finite coordinates, positive node IDs and radii, and parent -1 or positive");
		if (!index.emplace(node.id, static_cast<int>(i)).second)
			throw std::runtime_error("skeleton output contains duplicate node id "+std::to_string(node.id));
	}
	for (std::size_t i = 0; i < nodes.size(); ++i) {
		if (nodes[i].parent_id < 0) {
			if (topology.root >= 0) throw std::runtime_error("skeleton output contains multiple roots");
			topology.root = static_cast<int>(i);
			continue;
		}
		const auto found = index.find(nodes[i].parent_id);
		if (found == index.end())
			throw std::runtime_error("skeleton output parent id is undefined: "
				+std::to_string(nodes[i].parent_id));
		topology.parent[i] = found->second;
		topology.children[static_cast<std::size_t>(found->second)].push_back(static_cast<int>(i));
	}
	if (topology.root < 0) throw std::runtime_error("skeleton output contains no root");
	std::vector<int> state(nodes.size(), 0);
	int next_branch = 0;
	std::size_t visited = 0;
	const auto visit = [&](const auto& self, int current) -> void {
		if (state[static_cast<std::size_t>(current)] == 1)
			throw std::runtime_error("skeleton output contains a cycle");
		if (state[static_cast<std::size_t>(current)] == 2) return;
		state[static_cast<std::size_t>(current)] = 1;
		++visited;
		for (const int child : topology.children[static_cast<std::size_t>(current)]) {
			if (current == topology.root
				|| topology.children[static_cast<std::size_t>(current)].size() != 1)
				topology.incoming_branch[static_cast<std::size_t>(child)] = next_branch++;
			else topology.incoming_branch[static_cast<std::size_t>(child)] =
				topology.incoming_branch[static_cast<std::size_t>(current)];
			self(self, child);
		}
		state[static_cast<std::size_t>(current)] = 2;
	};
	visit(visit, topology.root);
	if (visited != nodes.size()) throw std::runtime_error("skeleton output contains disconnected nodes");
	topology.branches = next_branch;
	return topology;
}

inline void WriteNormalizedSkeletonSwc(const std::filesystem::path& path,
	const std::vector<SkeletonOutputNode>& nodes)
{
	AnalyzeSkeletonOutput(nodes);
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create normalized SWC: "+path.string());
	output << "# TubularFlowIGA normalized skeleton\n"
		<< "# id type x y z radius parent\n" << std::setprecision(17);
	for (const auto& node : nodes)
		output << node.id << ' ' << node.type << ' ' << node.position[0] << ' '
			<< node.position[1] << ' ' << node.position[2] << ' ' << node.radius << ' '
			<< node.parent_id << '\n';
	if (!output) throw std::runtime_error("cannot write normalized SWC: "+path.string());
}

inline void WriteSkeletonVtp(const std::filesystem::path& path,
	const std::vector<SkeletonOutputNode>& nodes)
{
	const auto topology = AnalyzeSkeletonOutput(nodes);
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot create skeleton VTP: "+path.string());
	output << std::setprecision(17)
		<< "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<< "  <PolyData>\n"
		<< "    <Piece NumberOfPoints=\"" << nodes.size()
		<< "\" NumberOfVerts=\"0\" NumberOfLines=\"" << nodes.size()-1
		<< "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n"
		<< "      <PointData Scalars=\"radius\">\n";
	auto integer_array = [&](const char* name, const auto& value) {
		output << "        <DataArray type=\"Int32\" Name=\"" << name
			<< "\" format=\"ascii\">\n          ";
		for (std::size_t i = 0; i < nodes.size(); ++i) output << value(i) << ' ';
		output << "\n        </DataArray>\n";
	};
	auto real_array = [&](const char* name, const auto& value) {
		output << "        <DataArray type=\"Float64\" Name=\"" << name
			<< "\" format=\"ascii\">\n          ";
		for (std::size_t i = 0; i < nodes.size(); ++i) output << value(i) << ' ';
		output << "\n        </DataArray>\n";
	};
	real_array("radius", [&](std::size_t i) { return nodes[i].radius; });
	real_array("diameter", [&](std::size_t i) { return 2.0*nodes[i].radius; });
	integer_array("node_id", [&](std::size_t i) { return nodes[i].id; });
	integer_array("parent_id", [&](std::size_t i) { return nodes[i].parent_id; });
	integer_array("degree", [&](std::size_t i) {
		return static_cast<int>(topology.children[i].size())+(topology.parent[i] >= 0 ? 1 : 0);
	});
	integer_array("role", [&](std::size_t i) {
		if (static_cast<int>(i) == topology.root) return 1;
		if (topology.children[i].size() > 1) return 2;
		if (topology.children[i].empty()) return 3;
		return 0;
	});
	output << "      </PointData>\n"
		<< "      <CellData>\n"
		<< "        <DataArray type=\"Int32\" Name=\"segment_id\" format=\"ascii\">\n          ";
	int segment = 0;
	for (std::size_t i = 0; i < nodes.size(); ++i)
		if (topology.parent[i] >= 0) output << segment++ << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int32\" Name=\"branch_id\" format=\"ascii\">\n          ";
	for (std::size_t i = 0; i < nodes.size(); ++i)
		if (topology.parent[i] >= 0) output << topology.incoming_branch[i] << ' ';
	output << "\n        </DataArray>\n"
		<< "      </CellData>\n"
		<< "      <Points>\n"
		<< "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n          ";
	for (const auto& node : nodes)
		output << node.position[0] << ' ' << node.position[1] << ' ' << node.position[2] << ' ';
	output << "\n        </DataArray>\n"
		<< "      </Points>\n"
		<< "      <Lines>\n"
		<< "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          ";
	for (std::size_t i = 0; i < nodes.size(); ++i)
		if (topology.parent[i] >= 0) output << topology.parent[i] << ' ' << i << ' ';
	output << "\n        </DataArray>\n"
		<< "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
	for (std::size_t i = 1; i < nodes.size(); ++i) output << 2*i << ' ';
	output << "\n        </DataArray>\n"
		<< "      </Lines>\n"
		<< "    </Piece>\n"
		<< "  </PolyData>\n"
		<< "</VTKFile>\n";
	if (!output) throw std::runtime_error("cannot write skeleton VTP: "+path.string());
}

} // namespace iga

#endif

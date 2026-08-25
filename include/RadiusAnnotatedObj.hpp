#ifndef IGA_RADIUS_ANNOTATED_OBJ_HPP
#define IGA_RADIUS_ANNOTATED_OBJ_HPP

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct RadiusAnnotatedObjNode {
	int id = -1;
	std::array<double, 3> position{{0.0, 0.0, 0.0}};
	double radius = 0.0;
	int parent = -1;
	std::vector<int> children;
};

struct RadiusAnnotatedObjTree {
	std::vector<RadiusAnnotatedObjNode> nodes;
	int root = -1;
};

inline bool IsRadiusAnnotatedObjPath(const std::filesystem::path& path)
{
	auto extension = path.extension().string();
	std::transform(extension.begin(), extension.end(), extension.begin(),
		[](unsigned char character) { return static_cast<char>(std::tolower(character)); });
	return extension == ".obj";
}

inline int RadiusAnnotatedObjIndex(const std::string& token,
	const std::filesystem::path& path, int line_number)
{
	std::size_t used = 0;
	int index = 0;
	try { index = std::stoi(token, &used); }
	catch (const std::exception&) { used = 0; }
	if (used != token.size() || index < 1)
		throw std::runtime_error(path.string()+":"+std::to_string(line_number)
			+": line indices must be positive vertex-only OBJ indices");
	return index-1;
}

inline RadiusAnnotatedObjTree ReadRadiusAnnotatedObj(
	const std::filesystem::path& path, int requested_root_id = 0)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open radius-annotated OBJ skeleton: "+path.string());
	RadiusAnnotatedObjTree result;
	std::vector<std::pair<int, int>> edges;
	std::set<std::pair<int, int>> unique_edges;
	std::string line;
	int line_number = 0;
	while (std::getline(input, line)) {
		++line_number;
		const auto first = line.find_first_not_of(" \t\r");
		if (first == std::string::npos || line[first] == '#') continue;
		std::istringstream row(line);
		std::string tag;
		row >> tag;
		if (tag == "v") {
			RadiusAnnotatedObjNode node;
			double auxiliary0 = 0.0;
			double auxiliary1 = 0.0;
			if (!(row >> node.position[0] >> node.position[1] >> node.position[2]
				>> node.radius >> auxiliary0 >> auxiliary1))
				throw std::runtime_error(path.string()+":"+std::to_string(line_number)
					+": expected 'v x y z radius auxiliary auxiliary'");
			std::string extra;
			if (row >> extra)
				throw std::runtime_error(path.string()+":"+std::to_string(line_number)
					+": unexpected radius-annotated OBJ vertex value");
			if (!std::isfinite(node.position[0]) || !std::isfinite(node.position[1])
				|| !std::isfinite(node.position[2]) || !std::isfinite(node.radius)
				|| !std::isfinite(auxiliary0) || !std::isfinite(auxiliary1)
				|| !(node.radius > 0.0))
				throw std::runtime_error(path.string()+":"+std::to_string(line_number)
					+": OBJ coordinates and positive radius must be finite");
			node.id = static_cast<int>(result.nodes.size())+1;
			result.nodes.push_back(node);
		} else if (tag == "l") {
			std::vector<int> indices;
			std::string token;
			while (row >> token) indices.push_back(RadiusAnnotatedObjIndex(token, path, line_number));
			if (indices.size() < 2)
				throw std::runtime_error(path.string()+":"+std::to_string(line_number)
					+": an OBJ line requires at least two vertex indices");
			for (std::size_t i = 1; i < indices.size(); ++i) {
				int first_index = indices[i-1];
				int second_index = indices[i];
				if (first_index == second_index)
					throw std::runtime_error(path.string()+":"+std::to_string(line_number)
						+": OBJ skeleton contains a self edge");
				if (first_index > second_index) std::swap(first_index, second_index);
				if (!unique_edges.emplace(first_index, second_index).second)
					throw std::runtime_error(path.string()+":"+std::to_string(line_number)
						+": OBJ skeleton contains a duplicate edge");
				edges.emplace_back(first_index, second_index);
			}
		} else {
			throw std::runtime_error(path.string()+":"+std::to_string(line_number)
				+": radius-annotated OBJ skeleton supports only v and l records");
		}
	}
	if (result.nodes.size() < 2 || edges.empty())
		throw std::runtime_error("radius-annotated OBJ skeleton requires vertices and edges");
	std::vector<std::vector<int>> adjacency(result.nodes.size());
	for (const auto& edge : edges) {
		if (edge.first >= static_cast<int>(result.nodes.size())
			|| edge.second >= static_cast<int>(result.nodes.size()))
			throw std::runtime_error("radius-annotated OBJ edge references an undefined vertex");
		const auto& first_node = result.nodes[static_cast<std::size_t>(edge.first)];
		const auto& second_node = result.nodes[static_cast<std::size_t>(edge.second)];
		double distance_squared = 0.0;
		for (int axis = 0; axis < 3; ++axis) {
			const double delta = second_node.position[static_cast<std::size_t>(axis)]
				-first_node.position[static_cast<std::size_t>(axis)];
			distance_squared += delta*delta;
		}
		if (!(distance_squared > 0.0))
			throw std::runtime_error("radius-annotated OBJ skeleton contains a zero-length edge");
		adjacency[static_cast<std::size_t>(edge.first)].push_back(edge.second);
		adjacency[static_cast<std::size_t>(edge.second)].push_back(edge.first);
	}
	if (requested_root_id != 0) {
		if (requested_root_id < 1 || requested_root_id > static_cast<int>(result.nodes.size()))
			throw std::runtime_error("OBJ root_node_id is outside the vertex range");
		result.root = requested_root_id-1;
	} else {
		double largest_radius = -std::numeric_limits<double>::infinity();
		for (std::size_t i = 0; i < result.nodes.size(); ++i)
			if (adjacency[i].size() == 1 && result.nodes[i].radius > largest_radius) {
				largest_radius = result.nodes[i].radius;
				result.root = static_cast<int>(i);
			}
		if (result.root < 0)
			throw std::runtime_error("radius-annotated OBJ skeleton has no terminal vertex for root inference");
	}
	std::vector<int> visited(result.nodes.size(), 0);
	std::queue<int> pending;
	visited[static_cast<std::size_t>(result.root)] = 1;
	pending.push(result.root);
	std::size_t visited_count = 0;
	while (!pending.empty()) {
		const int node = pending.front();
		pending.pop();
		++visited_count;
		auto& neighbors = adjacency[static_cast<std::size_t>(node)];
		std::sort(neighbors.begin(), neighbors.end());
		for (const int neighbor : neighbors) {
			if (neighbor == result.nodes[static_cast<std::size_t>(node)].parent) continue;
			if (visited[static_cast<std::size_t>(neighbor)])
				throw std::runtime_error("radius-annotated OBJ skeleton contains a cycle");
			visited[static_cast<std::size_t>(neighbor)] = 1;
			result.nodes[static_cast<std::size_t>(neighbor)].parent = node;
			result.nodes[static_cast<std::size_t>(node)].children.push_back(neighbor);
			pending.push(neighbor);
		}
	}
	if (visited_count != result.nodes.size())
		throw std::runtime_error("radius-annotated OBJ skeleton contains disconnected vertices");
	return result;
}

} // namespace iga

#endif

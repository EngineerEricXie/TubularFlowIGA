#ifndef IGA_ONE_D_NETWORK_HPP
#define IGA_ONE_D_NETWORK_HPP

#include "OneDConfig.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace iga {

constexpr double OneDPi = 3.141592653589793238462643383279502884;

struct OneDNode {
	int id = -1;
	int type = 0;
	std::array<double, 3> position{{0.0, 0.0, 0.0}};
	double radius = 0.0;
	int parent_id = -1;
	int parent = -1;
	std::vector<int> children;
};

struct OneDSegment {
	int index = -1;
	int parent = -1;
	int child = -1;
	double length = 0.0;
	double baseline_radius0 = 0.0;
	double radius0 = 0.0;
	double area0 = 0.0;
	double resistance = 0.0;
	std::array<double, 3> tangent{{0.0, 0.0, 0.0}};
	int cell_offset = 0;
	int cells = 1;
};

struct OneDNetwork {
	std::vector<OneDNode> nodes;
	std::vector<OneDSegment> segments;
	std::unordered_map<int, int> node_index;
	std::vector<int> topological_nodes;
	std::vector<int> outlet_nodes;
	int root = -1;
	int cells = 0;
};

inline double OneDDistance(const std::array<double, 3>& first,
	const std::array<double, 3>& second)
{
	const double dx = second[0]-first[0];
	const double dy = second[1]-first[1];
	const double dz = second[2]-first[2];
	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

inline OneDNetwork ReadOneDNetwork(const std::filesystem::path& path,
	double length_scale_to_m, int cells_per_segment, double dynamic_viscosity)
{
	if (!(length_scale_to_m > 0.0) || cells_per_segment < 1 || !(dynamic_viscosity > 0.0))
		throw std::runtime_error("invalid 1d network construction parameters");
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open SWC network: " + path.string());
	OneDNetwork network;
	std::string line;
	int line_number = 0;
	while (std::getline(input, line)) {
		++line_number;
		const auto first = line.find_first_not_of(" \t\r");
		if (first == std::string::npos || line[first] == '#') continue;
		std::istringstream row(line);
		OneDNode node;
		if (!(row >> node.id >> node.type >> node.position[0] >> node.position[1]
			>> node.position[2] >> node.radius >> node.parent_id))
			throw std::runtime_error(path.string() + ":" + std::to_string(line_number)
				+ ": expected seven SWC columns");
		std::string extra;
		if (row >> extra) throw std::runtime_error(path.string() + ":" + std::to_string(line_number)
			+ ": unexpected SWC column");
		if (node.id < 0 || !std::isfinite(node.radius) || !(node.radius > 0.0))
			throw std::runtime_error(path.string() + ":" + std::to_string(line_number)
				+ ": node id and radius must be valid");
		for (double& coordinate : node.position) {
			if (!std::isfinite(coordinate)) throw std::runtime_error("SWC coordinates must be finite");
			coordinate *= length_scale_to_m;
		}
		node.radius *= length_scale_to_m;
		if (!network.node_index.emplace(node.id, static_cast<int>(network.nodes.size())).second)
			throw std::runtime_error("duplicate SWC node id " + std::to_string(node.id));
		network.nodes.push_back(std::move(node));
	}
	if (network.nodes.empty()) throw std::runtime_error("SWC network is empty");

	int roots = 0;
	for (std::size_t i = 0; i < network.nodes.size(); ++i) {
		auto& node = network.nodes[i];
		if (node.parent_id == -1) {
			network.root = static_cast<int>(i);
			++roots;
			continue;
		}
		const auto found = network.node_index.find(node.parent_id);
		if (found == network.node_index.end())
			throw std::runtime_error("SWC node " + std::to_string(node.id)
				+ " references missing parent " + std::to_string(node.parent_id));
		node.parent = found->second;
		network.nodes[static_cast<std::size_t>(node.parent)].children.push_back(static_cast<int>(i));
	}
	if (roots != 1) throw std::runtime_error("SWC network must contain exactly one root");

	std::vector<int> state(network.nodes.size(), 0);
	std::function<void(int)> visit = [&](int index) {
		if (state[static_cast<std::size_t>(index)] == 1) throw std::runtime_error("SWC network contains a cycle");
		if (state[static_cast<std::size_t>(index)] == 2) return;
		state[static_cast<std::size_t>(index)] = 1;
		network.topological_nodes.push_back(index);
		for (const int child : network.nodes[static_cast<std::size_t>(index)].children) visit(child);
		state[static_cast<std::size_t>(index)] = 2;
	};
	visit(network.root);
	if (network.topological_nodes.size() != network.nodes.size())
		throw std::runtime_error("SWC network contains nodes disconnected from the root");

	for (const int parent : network.topological_nodes) {
		const auto& parent_node = network.nodes[static_cast<std::size_t>(parent)];
		if (parent_node.children.empty()) network.outlet_nodes.push_back(parent);
		for (const int child : parent_node.children) {
			const auto& child_node = network.nodes[static_cast<std::size_t>(child)];
			OneDSegment segment;
			segment.index = static_cast<int>(network.segments.size());
			segment.parent = parent;
			segment.child = child;
			segment.length = OneDDistance(parent_node.position, child_node.position);
			if (!(segment.length > 0.0))
				throw std::runtime_error("SWC segment " + std::to_string(parent_node.id)
					+ "->" + std::to_string(child_node.id) + " has zero length");
			segment.radius0 = 0.5*(parent_node.radius+child_node.radius);
			segment.baseline_radius0 = segment.radius0;
			segment.area0 = OneDPi*segment.radius0*segment.radius0;
			segment.resistance = 8.0*dynamic_viscosity*segment.length
				/(OneDPi*std::pow(segment.radius0, 4.0));
			for (int axis = 0; axis < 3; ++axis)
				segment.tangent[static_cast<std::size_t>(axis)] =
					(child_node.position[static_cast<std::size_t>(axis)]
						-parent_node.position[static_cast<std::size_t>(axis)])/segment.length;
			segment.cell_offset = network.cells;
			segment.cells = cells_per_segment;
			network.cells += cells_per_segment;
			network.segments.push_back(segment);
		}
	}
	return network;
}

inline int OneDSegmentIntoNode(const OneDNetwork& network, int node)
{
	for (const auto& segment : network.segments)
		if (segment.child == node) return segment.index;
	return -1;
}

inline std::vector<int> OneDSegmentsOutOfNode(const OneDNetwork& network, int node)
{
	std::vector<int> result;
	for (const auto& segment : network.segments)
		if (segment.parent == node) result.push_back(segment.index);
	return result;
}

inline double OneDBranchAngleDegrees(const OneDSegment& parent,
	const OneDSegment& child)
{
	double cosine = 0.0;
	for (int axis = 0; axis < 3; ++axis)
		cosine += parent.tangent[static_cast<std::size_t>(axis)]
			*child.tangent[static_cast<std::size_t>(axis)];
	cosine = std::max(-1.0, std::min(1.0, cosine));
	return std::acos(cosine)*180.0/OneDPi;
}

inline void ValidateOneDTopologyReferences(const OneDConfiguration& configuration,
	const OneDNetwork& network)
{
	bool inlet = false;
	std::set<int> configured_outlets;
	for (const auto& boundary : configuration.boundaries) {
		for (const int id : boundary.node_ids) {
			const auto found = network.node_index.find(id);
			if (found == network.node_index.end())
				throw std::runtime_error("1d boundary references unknown SWC node " + std::to_string(id));
			if (boundary.role == "inlet" && found->second != network.root)
				throw std::runtime_error("1d inlet node must be the SWC root");
			if (boundary.role == "outlet"
				&& !network.nodes[static_cast<std::size_t>(found->second)].children.empty())
				throw std::runtime_error("1d outlet node must be a leaf");
			if (boundary.role == "outlet") configured_outlets.insert(found->second);
		}
		if (boundary.role == "inlet") {
			if (inlet) throw std::runtime_error("1d configuration requires exactly one inlet boundary");
			inlet = true;
			if (!boundary.node_ids.empty()
				&& (boundary.node_ids.size() != 1 || network.node_index.at(boundary.node_ids.front()) != network.root))
				throw std::runtime_error("1d inlet boundary must select only the root");
		}
	}
	if (!inlet) throw std::runtime_error("1d configuration requires an inlet boundary");
	for (const auto& flow : configuration.flow_systems)
		for (const auto& coefficient : flow.junctions.node_coefficients) {
			const auto found = network.node_index.find(coefficient.first);
			if (found == network.node_index.end()
				|| network.nodes[static_cast<std::size_t>(found->second)].children.size() < 2)
				throw std::runtime_error("junction coefficient must reference a branching SWC node");
		}
	(void)configured_outlets;
}

} // namespace iga

#endif

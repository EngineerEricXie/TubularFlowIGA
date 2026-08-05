#include "SwcGraph.hpp"
#include "BSpline.hpp"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace tubular {

MeshParameters MeshParameters::Read(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open mesh parameters: "+path.string());
	MeshParameters result;
	std::string name;
	if (!(input >> name >> result.noise_iterations)
		|| !(input >> name >> result.bifurcation_smoothing)
		|| !(input >> name >> result.noise_smoothing)
		|| !(input >> name >> result.segment_length)
		|| !(input >> name >> result.bifurcation_refinement))
		throw std::runtime_error("mesh parameter file must contain five name/value pairs");
	if (result.noise_iterations < 0)
		throw std::runtime_error("noise iterations must be non-negative");
	if (!std::isfinite(result.bifurcation_smoothing) || result.bifurcation_smoothing < 0.0 || result.bifurcation_smoothing > 1.0)
		throw std::runtime_error("bifurcation smoothing ratio must be in [0,1]");
	if (!std::isfinite(result.noise_smoothing) || result.noise_smoothing < 0.0 || result.noise_smoothing > 1.0)
		throw std::runtime_error("noise smoothing ratio must be in [0,1]");
	if (!std::isfinite(result.segment_length) || result.segment_length <= 0.0)
		throw std::runtime_error("segment length must be positive");
	if (!std::isfinite(result.bifurcation_refinement) || result.bifurcation_refinement < 0.0)
		throw std::runtime_error("bifurcation refinement ratio must be non-negative");
	return result;
}

SwcGraph SwcGraph::Read(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open SWC file: "+path.string());
	struct Raw { int id; int type; Vec3 p; double radius; int parent; };
	std::vector<Raw> raw;
	std::string line;
	while (std::getline(input, line)) {
		const auto first = line.find_first_not_of(" \t\r");
		if (first == std::string::npos || line[first] == '#') continue;
		std::istringstream row(line);
		Raw value;
		if (!(row >> value.id >> value.type >> value.p.x >> value.p.y >> value.p.z >> value.radius >> value.parent))
			throw std::runtime_error("invalid SWC row: "+line);
		if (value.id <= 0 || !IsFinite(value.p) || !std::isfinite(value.radius) || value.radius <= 0.0)
			throw std::runtime_error("invalid SWC node values at id "+std::to_string(value.id));
		raw.push_back(value);
	}
	if (raw.empty()) throw std::runtime_error("SWC file contains no nodes");
	std::unordered_map<int,int> index;
	for (std::size_t i=0; i<raw.size(); ++i) {
		if (!index.emplace(raw[i].id, static_cast<int>(i)).second)
			throw std::runtime_error("duplicate SWC id "+std::to_string(raw[i].id));
	}
	SwcGraph graph;
	graph.nodes.resize(raw.size());
	for (std::size_t i=0; i<raw.size(); ++i) {
		auto& node = graph.nodes[i];
		node.type = raw[i].type;
		node.position = raw[i].p;
		node.diameter = 2.0*raw[i].radius;
		if (raw[i].parent >= 0) {
			const auto found = index.find(raw[i].parent);
			if (found == index.end())
				throw std::runtime_error("SWC parent id not found: "+std::to_string(raw[i].parent));
			node.parent = found->second;
		}
	}
	graph.RebuildChildren();
	graph.Validate();
	return graph;
}

void SwcGraph::Write(const std::filesystem::path& path) const
{
	std::ofstream output(path);
	if (!output) throw std::runtime_error("cannot write SWC file: "+path.string());
	output << "# TubularFlowIGA smoothed skeleton\n"
		<< "# id type x y z radius parent\n";
	output << std::fixed << std::setprecision(8);
	for (std::size_t i=0; i<nodes.size(); ++i) {
		const auto& n = nodes[i];
		output << i+1 << ' ' << 2 << ' ' << n.position.x << ' ' << n.position.y << ' '
			<< n.position.z << ' ' << n.diameter/2.0 << ' '
			<< (n.parent < 0 ? -1 : n.parent+1) << '\n';
	}
}

void SwcGraph::RebuildChildren()
{
	for (auto& node : nodes) node.children.clear();
	for (std::size_t i=0; i<nodes.size(); ++i)
		if (nodes[i].parent >= 0) nodes.at(nodes[i].parent).children.push_back(static_cast<int>(i));
	for (auto& node : nodes) std::sort(node.children.begin(), node.children.end());
}

int SwcGraph::root() const
{
	int result = -1;
	for (std::size_t i=0; i<nodes.size(); ++i) {
		if (nodes[i].parent < 0) {
			if (result >= 0) throw std::runtime_error("SWC graph has multiple roots");
			result = static_cast<int>(i);
		}
	}
	if (result < 0) throw std::runtime_error("SWC graph has no root");
	return result;
}

void SwcGraph::Validate() const
{
	const int root_index = root();
	std::vector<int> state(nodes.size(), 0);
	std::function<void(int)> visit = [&](int i) {
		if (state[i] == 1) throw std::runtime_error("SWC graph contains a cycle");
		if (state[i] == 2) return;
		state[i] = 1;
		if (nodes[i].children.size() > 2)
			throw std::runtime_error("only binary branching is supported; node "+std::to_string(i+1)+" has more than two children");
		if (!IsFinite(nodes[i].position) || !std::isfinite(nodes[i].diameter) || nodes[i].diameter <= 0.0)
			throw std::runtime_error("invalid geometry at SWC node "+std::to_string(i+1));
		for (int child : nodes[i].children) visit(child);
		state[i] = 2;
	};
	visit(root_index);
	if (std::find(state.begin(), state.end(), 0) != state.end())
		throw std::runtime_error("SWC graph contains nodes disconnected from the root");
}

std::vector<std::vector<int>> SwcGraph::Sections() const
{
	std::vector<std::vector<int>> result;
	std::function<void(int)> visit = [&](int start) {
		for (int child : nodes[start].children) {
			std::vector<int> path{start, child};
			int current = child;
			while (nodes[current].children.size() == 1) {
				current = nodes[current].children.front();
				path.push_back(current);
			}
			result.push_back(path);
			if (!nodes[current].children.empty()) visit(current);
		}
	};
	visit(root());
	if (result.empty()) throw std::runtime_error("SWC graph contains no section");
	return result;
}

SwcGraph SmoothSkeleton(const SwcGraph& input, const MeshParameters& parameters)
{
	SwcGraph work = input;
	const auto sections = work.Sections();
	std::vector<int> branches;
	for (std::size_t i=0; i<work.nodes.size(); ++i)
		if (work.is_branch(static_cast<int>(i))) branches.push_back(static_cast<int>(i));

	for (int iteration=0; iteration<parameters.noise_iterations; ++iteration) {
		for (int b : branches) {
			auto& node = work.nodes[b];
			const auto& parent = work.nodes.at(node.parent);
			const auto& child0 = work.nodes.at(node.children[0]);
			const auto& child1 = work.nodes.at(node.children[1]);
			const Vec3 average = (parent.position+child0.position+child1.position)/3.0;
			const double average_d = (parent.diameter+child0.diameter+child1.diameter)/3.0;
			node.position = node.position*(1.0-parameters.bifurcation_smoothing)+average*parameters.bifurcation_smoothing;
			node.diameter = node.diameter*(1.0-parameters.bifurcation_smoothing)+average_d*parameters.bifurcation_smoothing;
		}
		for (const auto& section : sections) {
			for (std::size_t j=1; j+1<section.size(); ++j) {
				auto& current = work.nodes[section[j]];
				const auto& previous = work.nodes[section[j-1]];
				const auto& next = work.nodes[section[j+1]];
				current.position += ((previous.position+next.position)/2.0-current.position)*parameters.noise_smoothing;
				current.diameter += ((previous.diameter+next.diameter)/2.0-current.diameter)*parameters.noise_smoothing;
			}
		}
	}

	std::vector<int> critical_map(work.nodes.size(), -1);
	SwcGraph output;
	for (std::size_t i=0; i<work.nodes.size(); ++i) {
		if (work.nodes[i].parent < 0 || work.nodes[i].children.size() != 1) {
			critical_map[i] = static_cast<int>(output.nodes.size());
			SwcNode node = work.nodes[i];
			node.type = 2;
			node.parent = -1;
			node.children.clear();
			output.nodes.push_back(node);
		}
	}
	const bool no_bifurcations = branches.empty();
	for (const auto& section : sections) {
		std::vector<Vec3> points;
		std::vector<double> diameters;
		for (int index : section) {
			points.push_back(work.nodes[index].position);
			diameters.push_back(work.nodes[index].diameter);
		}
		int mode = 0;
		if (no_bifurcations) mode = 4;
		else if (work.is_branch(section.front()) && work.is_branch(section.back())) mode = 1;
		else if (work.is_terminal(section.back())) mode = 2;
		else if (section.front() == work.root()) mode = 3;
		else throw std::runtime_error("cannot classify skeleton section");
		const auto samples = SampleBranch(points, diameters, parameters.segment_length, mode);
		int parent = critical_map.at(section.front());
		if (parent < 0) throw std::runtime_error("section start is not a critical node");
		for (std::size_t j=1; j+1<samples.size(); ++j) {
			SwcNode node;
			node.type = 2;
			node.position = samples[j].point;
			node.diameter = samples[j].diameter;
			node.parent = parent;
			parent = static_cast<int>(output.nodes.size());
			output.nodes.push_back(node);
		}
		const int end = critical_map.at(section.back());
		if (end < 0) throw std::runtime_error("section end is not a critical node");
		if (output.nodes[end].parent >= 0)
			throw std::runtime_error("critical node received multiple parents");
		output.nodes[end].parent = parent;
	}
	output.RebuildChildren();
	output.Validate();
	return output;
}

} // namespace tubular

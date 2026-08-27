#include "SwcGraph.hpp"
#include "BSpline.hpp"
#include "RadiusAnnotatedObj.hpp"
#include "SkeletonOutput.hpp"

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
	std::string noise_name, bifurcation_name, smoothing_name, segment_name, refinement_name;
	if (!(input >> noise_name >> result.noise_iterations)
		|| !(input >> bifurcation_name >> result.bifurcation_smoothing)
		|| !(input >> smoothing_name >> result.noise_smoothing)
		|| !(input >> segment_name >> result.segment_length)
		|| !(input >> refinement_name >> result.bifurcation_refinement))
		throw std::runtime_error("mesh parameter file must contain five name/value pairs");
	if (noise_name != "n_noisesmooth" || bifurcation_name != "ratio_bifur_node"
		|| smoothing_name != "ratio_noisesmooth" || segment_name != "seg_length"
		|| refinement_name != "ratio_refine")
		throw std::runtime_error("legacy mesh parameter names or order are invalid");
	std::string extra;
	if (input >> extra) throw std::runtime_error("legacy mesh parameter file contains unexpected data");
	result.max_spacing_over_diameter = 1.1;
	result.minimum_scaled_jacobian = 1.0e-3;
	result.junction_optimization_iterations = 0;
	result.check_self_intersection = true;
	result.Validate();
	return result;
}

void MeshParameters::Validate() const
{
	if (noise_iterations < 0)
		throw std::runtime_error("noise iterations must be non-negative");
	if (!std::isfinite(bifurcation_smoothing) || bifurcation_smoothing < 0.0 || bifurcation_smoothing > 1.0)
		throw std::runtime_error("bifurcation smoothing ratio must be in [0,1]");
	if (!std::isfinite(noise_smoothing) || noise_smoothing < 0.0 || noise_smoothing > 1.0)
		throw std::runtime_error("noise smoothing ratio must be in [0,1]");
	if (!std::isfinite(segment_length) || segment_length <= 0.0)
		throw std::runtime_error("segment length must be positive");
	if (!std::isfinite(max_spacing_over_diameter) || max_spacing_over_diameter <= 0.0
		|| !std::isfinite(max_turn_degrees) || max_turn_degrees <= 0.0 || max_turn_degrees > 90.0
		|| !std::isfinite(max_diameter_change_fraction) || max_diameter_change_fraction <= 0.0
		|| max_diameter_change_fraction > 1.0
		|| !std::isfinite(maximum_curvature_radius_product)
		|| maximum_curvature_radius_product <= 0.0 || maximum_curvature_radius_product >= 1.0)
		throw std::runtime_error("adaptive centerline controls are invalid");
	if (!std::isfinite(bifurcation_refinement) || bifurcation_refinement < 0.0)
		throw std::runtime_error("bifurcation refinement ratio must be non-negative");
	if (!std::isfinite(upstream_clearance_over_diameter) || upstream_clearance_over_diameter <= 0.0
		|| !std::isfinite(downstream_clearance_over_diameter) || downstream_clearance_over_diameter <= 0.0
		|| !std::isfinite(minimum_bifurcation_angle_degrees)
		|| minimum_bifurcation_angle_degrees <= 0.0 || minimum_bifurcation_angle_degrees >= 90.0
		|| !std::isfinite(maximum_junction_radius_ratio) || maximum_junction_radius_ratio < 1.0
		|| junction_optimization_iterations < 0 || junction_optimization_iterations > 1000)
		throw std::runtime_error("junction controls are invalid");
	if (!std::isfinite(minimum_scaled_jacobian) || minimum_scaled_jacobian <= 0.0
		|| minimum_scaled_jacobian > 1.0 || !std::isfinite(collision_safety_factor)
		|| collision_safety_factor <= 0.0)
		throw std::runtime_error("mesh quality controls are invalid");
}

SwcGraph SwcGraph::Read(const std::filesystem::path& path)
{
	if (iga::IsRadiusAnnotatedObjPath(path)) {
		const auto input = iga::ReadRadiusAnnotatedObj(path);
		SwcGraph graph;
		graph.nodes.resize(input.nodes.size());
		for (std::size_t i = 0; i < input.nodes.size(); ++i) {
			const auto& source = input.nodes[i];
			auto& target = graph.nodes[i];
			target.id = source.id;
			target.type = 2;
			target.position = {source.position[0], source.position[1], source.position[2]};
			target.diameter = 2.0*source.radius;
			target.parent = source.parent;
		}
		graph.RebuildChildren();
		graph.Validate();
		return graph;
	}
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open SWC file: "+path.string());
	struct Raw { int id; int type; Vec3 p; double radius; int parent; };
	std::vector<Raw> raw;
	std::string line;
	int line_number = 0;
	while (std::getline(input, line)) {
		++line_number;
		const auto first = line.find_first_not_of(" \t\r");
		if (first == std::string::npos || line[first] == '#') continue;
		std::istringstream row(line);
		Raw value;
		if (!(row >> value.id >> value.type >> value.p.x >> value.p.y >> value.p.z >> value.radius >> value.parent))
			throw std::runtime_error(path.string()+":"+std::to_string(line_number)
				+": expected seven SWC columns");
		std::string extra;
		if (row >> extra) throw std::runtime_error(path.string()+":"+std::to_string(line_number)
			+": unexpected SWC column");
		if (value.id <= 0 || (value.parent != -1 && value.parent <= 0)
			|| !IsFinite(value.p) || !std::isfinite(value.radius) || value.radius <= 0.0)
			throw std::runtime_error(path.string()+":"+std::to_string(line_number)
				+": invalid SWC node values at id "+std::to_string(value.id));
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
		node.id = raw[i].id;
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

namespace {

std::vector<iga::SkeletonOutputNode> SkeletonOutputNodes(const SwcGraph& graph)
{
	std::vector<iga::SkeletonOutputNode> result(graph.nodes.size());
	for (std::size_t i = 0; i < graph.nodes.size(); ++i) {
		const auto& source = graph.nodes[i];
		auto& target = result[i];
		target.id = source.id > 0 ? source.id : static_cast<int>(i)+1;
		target.type = source.type;
		target.position = {{source.position.x, source.position.y, source.position.z}};
		target.radius = source.diameter/2.0;
		target.parent_id = source.parent < 0 ? -1
			: (graph.nodes[static_cast<std::size_t>(source.parent)].id > 0
				? graph.nodes[static_cast<std::size_t>(source.parent)].id : source.parent+1);
	}
	return result;
}

}

void SwcGraph::WriteNormalized(const std::filesystem::path& path) const
{
	Validate();
	iga::WriteNormalizedSkeletonSwc(path, SkeletonOutputNodes(*this));
}

void SwcGraph::WriteVisualizationVtp(const std::filesystem::path& path) const
{
	Validate();
	iga::WriteSkeletonVtp(path, SkeletonOutputNodes(*this));
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
			throw std::runtime_error("3d mesh generation supports at most two children; node "
				+std::to_string(nodes[i].id > 0 ? nodes[i].id : static_cast<int>(i)+1)
				+" has "+std::to_string(nodes[i].children.size()));
		if (!IsFinite(nodes[i].position) || !std::isfinite(nodes[i].diameter) || nodes[i].diameter <= 0.0)
			throw std::runtime_error("invalid geometry at SWC node "
				+std::to_string(nodes[i].id > 0 ? nodes[i].id : static_cast<int>(i)+1));
		if (nodes[i].parent >= 0
			&& !(Norm(nodes[i].position-nodes[static_cast<std::size_t>(nodes[i].parent)].position) > 0.0))
			throw std::runtime_error("zero-length segment at SWC node "
				+std::to_string(nodes[i].id > 0 ? nodes[i].id : static_cast<int>(i)+1));
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
	parameters.Validate();
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
	BranchSamplingOptions sampling;
	sampling.target_spacing = parameters.segment_length;
	sampling.max_spacing_over_diameter = parameters.max_spacing_over_diameter;
	sampling.max_turn_degrees = parameters.max_turn_degrees;
	sampling.max_diameter_change_fraction = parameters.max_diameter_change_fraction;
	sampling.upstream_clearance_over_diameter = parameters.upstream_clearance_over_diameter;
	sampling.downstream_clearance_over_diameter = parameters.downstream_clearance_over_diameter;
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
		const auto samples = SampleBranch(points, diameters, sampling, mode,
			"section "+std::to_string(work.nodes[section.front()].id)
			+"->"+std::to_string(work.nodes[section.back()].id));
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

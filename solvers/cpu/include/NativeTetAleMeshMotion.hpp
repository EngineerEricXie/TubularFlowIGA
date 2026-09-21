#ifndef IGA_NATIVE_TET_ALE_MESH_MOTION_HPP
#define IGA_NATIVE_TET_ALE_MESH_MOTION_HPP

#include "NativeTetFem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetAleHarmonicResult
{
	std::vector<std::array<double,3>> displacement_m;
	std::size_t constrained_nodes = 0;
	std::size_t free_nodes = 0;
	double maximum_free_residual = 0.0;
};

namespace native_tet_ale_motion_detail {

inline std::vector<double> SolveDense(std::vector<double> matrix,
	std::vector<double> right_hand_side)
{
	const std::size_t size = right_hand_side.size();
	if (matrix.size() != size*size)
		throw std::invalid_argument("native ALE harmonic matrix size is invalid");
	for (std::size_t column = 0; column < size; ++column) {
		std::size_t pivot = column;
		for (std::size_t row = column+1; row < size; ++row)
			if (std::abs(matrix[row*size+column]) > std::abs(matrix[pivot*size+column]))
				pivot = row;
		if (!(std::abs(matrix[pivot*size+column]) >
			100.0*std::numeric_limits<double>::epsilon()))
			throw std::runtime_error("native ALE harmonic system is singular");
		if (pivot != column) {
			for (std::size_t entry = column; entry < size; ++entry)
				std::swap(matrix[column*size+entry],matrix[pivot*size+entry]);
			std::swap(right_hand_side[column],right_hand_side[pivot]);
		}
		for (std::size_t row = column+1; row < size; ++row) {
			const double factor = matrix[row*size+column]/matrix[column*size+column];
			matrix[row*size+column] = 0.0;
			for (std::size_t entry = column+1; entry < size; ++entry)
				matrix[row*size+entry] -= factor*matrix[column*size+entry];
			right_hand_side[row] -= factor*right_hand_side[column];
		}
	}
	std::vector<double> solution(size,0.0);
	for (std::size_t reverse = size; reverse > 0; --reverse) {
		const std::size_t row = reverse-1;
		double value = right_hand_side[row];
		for (std::size_t column = row+1; column < size; ++column)
			value -= matrix[row*size+column]*solution[column];
		solution[row] = value/matrix[row*size+row];
	}
	return solution;
}

} // namespace native_tet_ale_motion_detail

// Componentwise P1 harmonic extension on the reference tetrahedral mesh.
// The caller supplies all moving/fixed/sliding Dirichlet values explicitly.
inline NativeTetAleHarmonicResult SolveNativeTetAleHarmonicMotion(
	const NativeTetMesh& mesh,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_displacement_m)
{
	if (mesh.points.empty() || mesh.cells.empty() || prescribed_displacement_m.empty())
		throw std::invalid_argument("native ALE harmonic motion requires a mesh and constraints");
	const std::size_t node_count = mesh.points.size();
	std::vector<double> stiffness(node_count*node_count,0.0);
	for (const auto& cell : mesh.cells) {
		for (const auto node : cell.nodes)
			if (node >= node_count)
				throw std::invalid_argument("native ALE harmonic cell index is out of range");
		const auto geometry = EvaluateNativeTetGeometry(mesh,cell);
		const double volume = geometry.determinant/6.0;
		for (std::size_t row = 0; row < 4; ++row)
			for (std::size_t column = 0; column < 4; ++column) {
				double value = 0.0;
				for (int component = 0; component < 3; ++component)
					value += geometry.barycentric_gradients[row][component]
						*geometry.barycentric_gradients[column][component];
				stiffness[cell.nodes[row]*node_count+cell.nodes[column]] += volume*value;
			}
	}
	std::vector<bool> constrained(node_count,false);
	NativeTetAleHarmonicResult result;
	result.displacement_m.assign(node_count,{{0,0,0}});
	for (const auto& prescribed : prescribed_displacement_m) {
		if (prescribed.first >= node_count)
			throw std::invalid_argument("native ALE harmonic constraint index is out of range");
		for (const double value : prescribed.second)
			if (!std::isfinite(value))
				throw std::invalid_argument("native ALE harmonic constraint is nonfinite");
		constrained[prescribed.first] = true;
		result.displacement_m[prescribed.first] = prescribed.second;
	}
	std::vector<std::size_t> free_nodes;
	for (std::size_t node = 0; node < node_count; ++node)
		if (!constrained[node]) free_nodes.push_back(node);
	result.constrained_nodes = prescribed_displacement_m.size();
	result.free_nodes = free_nodes.size();
	if (free_nodes.empty()) return result;
	const std::size_t free_count = free_nodes.size();
	std::vector<double> reduced(free_count*free_count,0.0);
	for (std::size_t row = 0; row < free_count; ++row)
		for (std::size_t column = 0; column < free_count; ++column)
			reduced[row*free_count+column] =
				stiffness[free_nodes[row]*node_count+free_nodes[column]];
	for (int component = 0; component < 3; ++component) {
		std::vector<double> right_hand_side(free_count,0.0);
		for (std::size_t row = 0; row < free_count; ++row)
			for (const auto& prescribed : prescribed_displacement_m)
				right_hand_side[row] -= stiffness[free_nodes[row]*node_count+prescribed.first]
					*prescribed.second[component];
		const auto solution = native_tet_ale_motion_detail::SolveDense(
			reduced,right_hand_side);
		for (std::size_t row = 0; row < free_count; ++row)
			result.displacement_m[free_nodes[row]][component] = solution[row];
	}
	for (const auto node : free_nodes)
		for (int component = 0; component < 3; ++component) {
			double residual = 0.0;
			for (std::size_t column = 0; column < node_count; ++column)
				residual += stiffness[node*node_count+column]
					*result.displacement_m[column][component];
			result.maximum_free_residual = std::max(result.maximum_free_residual,
				std::abs(residual));
		}
	return result;
}

} // namespace iga

#endif

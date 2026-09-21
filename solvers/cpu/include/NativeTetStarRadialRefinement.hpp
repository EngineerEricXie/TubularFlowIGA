#ifndef IGA_NATIVE_TET_STAR_RADIAL_REFINEMENT_HPP
#define IGA_NATIVE_TET_STAR_RADIAL_REFINEMENT_HPP

#include "NativeTetAleKinematics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iga {

// Refine a center-fan chamber without changing a single exterior triangle.
// Each radial triangular prism uses the same global vertex ordering, so the
// diagonals on shared radial quadrilaterals agree across neighboring cells.
inline NativeTetMesh RefineNativeTetStarRadialLayers(const NativeTetMesh& source,
	std::uint32_t center,const std::vector<double>& radial_fractions)
{
	if(center>=source.points.size()||source.cells.empty()
		||radial_fractions.empty())
		throw std::invalid_argument("native star radial refinement input is invalid");
	double previous_fraction=0.0;
	for(const double fraction:radial_fractions){
		if(!(fraction>previous_fraction&&fraction<1.0)||!std::isfinite(fraction))
			throw std::invalid_argument("native star radial fractions must increase in (0,1)");
		previous_fraction=fraction;
	}
	if(radial_fractions.size()>(std::numeric_limits<std::size_t>::max()-1)/3)
		throw std::invalid_argument("native star radial refinement has too many layers");
	const std::size_t surface_nodes=source.points.size()-1;
	if(surface_nodes>(std::numeric_limits<std::uint32_t>::max()-source.points.size())
		/radial_fractions.size())
		throw std::invalid_argument("native star radial refinement has too many nodes");
	NativeTetMesh result;
	result.points=source.points;
	result.boundary_triangles=source.boundary_triangles;
	std::vector<std::vector<std::uint32_t>> shells;
	shells.reserve(radial_fractions.size());
	for(const double fraction:radial_fractions){
		std::vector<std::uint32_t> shell(source.points.size(),
			std::numeric_limits<std::uint32_t>::max());
		for(std::size_t node=0;node<source.points.size();++node){
			if(node==center)continue;
			shell[node]=static_cast<std::uint32_t>(result.points.size());
			std::array<double,3> point{};
			for(int component=0;component<3;++component)
				point[component]=source.points[center][component]
					+fraction*(source.points[node][component]
						-source.points[center][component]);
			result.points.push_back(point);
		}
		shells.push_back(std::move(shell));
	}
	if(source.cells.size()>std::numeric_limits<std::size_t>::max()
		/(1+3*radial_fractions.size()))
		throw std::invalid_argument("native star radial refinement has too many cells");
	result.cells.reserve((1+3*radial_fractions.size())*source.cells.size());
	std::uint64_t id=1;
	const auto append=[&](std::array<std::uint32_t,4> nodes){
		const double determinant=NativeAleDeterminant(result.points[nodes[0]],
			result.points[nodes[1]],result.points[nodes[2]],result.points[nodes[3]]);
		if(!(std::abs(determinant)>0.0)||!std::isfinite(determinant))
			throw std::runtime_error("native star radial refinement has a degenerate cell");
		if(determinant<0.0)std::swap(nodes[2],nodes[3]);
		result.cells.push_back({id++,nodes});
	};
	for(const auto& cell:source.cells){
		std::array<std::uint32_t,3> outer{};std::size_t count=0;
		for(const auto node:cell.nodes){
			if(node>=source.points.size())
				throw std::invalid_argument("native star cell node is out of range");
			if(node!=center){
				if(count>=3)throw std::invalid_argument("native star cell omits center");
				outer[count++]=node;
			}
		}
		if(count!=3||outer[0]==outer[1]||outer[0]==outer[2]||outer[1]==outer[2])
			throw std::invalid_argument("native star cell is not center plus a triangle");
		std::sort(outer.begin(),outer.end());
		const auto a=outer[0],b=outer[1],c=outer[2];
		append({{center,shells.front()[a],shells.front()[b],shells.front()[c]}});
		for(std::size_t layer=0;layer<shells.size();++layer){
			const auto& inner=shells[layer];
			const auto outer_a=layer+1<shells.size()?shells[layer+1][a]:a;
			const auto outer_b=layer+1<shells.size()?shells[layer+1][b]:b;
			const auto outer_c=layer+1<shells.size()?shells[layer+1][c]:c;
			append({{inner[a],inner[b],inner[c],outer_a}});
			append({{outer_a,inner[b],inner[c],outer_b}});
			append({{outer_a,outer_b,inner[c],outer_c}});
		}
	}
	return result;
}

inline NativeTetMesh RefineNativeTetStarRadially(const NativeTetMesh& source,
	std::uint32_t center,double outer_fraction)
{
	return RefineNativeTetStarRadialLayers(source,center,{outer_fraction});
}

} // namespace iga

#endif

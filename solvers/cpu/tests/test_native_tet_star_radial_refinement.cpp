#include "IdealizedLeftVentricleFixture.hpp"
#include "NativeTetStarRadialRefinement.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

int main()
{
	const auto surface=iga::IdealizedLeftVentricleFixture::Motion().Frames().front().surface;
	iga::NativeTetMesh star;
	star.points=surface.vertices;
	const auto center=static_cast<std::uint32_t>(star.points.size());
	star.points.push_back({{0.,0.,.045}});
	std::uint64_t id=1;
	for(const auto& triangle:surface.triangles){
		auto nodes=std::array<std::uint32_t,4>{{center,
			static_cast<std::uint32_t>(triangle.indices[0]),
			static_cast<std::uint32_t>(triangle.indices[1]),
			static_cast<std::uint32_t>(triangle.indices[2])}};
		if(iga::NativeAleDeterminant(star.points[nodes[0]],star.points[nodes[1]],
			star.points[nodes[2]],star.points[nodes[3]])<0.0)std::swap(nodes[2],nodes[3]);
		star.cells.push_back({id,nodes});
		star.boundary_triangles.push_back({id,{{
			static_cast<std::uint32_t>(triangle.indices[0]),
			static_cast<std::uint32_t>(triangle.indices[1]),
			static_cast<std::uint32_t>(triangle.indices[2])}},
			static_cast<int>(triangle.boundary_id)});++id;
	}
	const auto refined=iga::RefineNativeTetStarRadially(star,center,0.5);
	const auto two_shells=iga::RefineNativeTetStarRadialLayers(star,center,{1./3.,2./3.});
	assert(refined.points.size()==2*star.points.size()-1);
	assert(refined.cells.size()==4*star.cells.size());
	assert(two_shells.points.size()==3*star.points.size()-2);
	assert(two_shells.cells.size()==7*star.cells.size());
	double original_volume=0.0;
	for(const auto& cell:star.cells)
		original_volume+=iga::EvaluateNativeTetGeometry(star,cell).determinant/6.0;
	const auto verify=[&](const iga::NativeTetMesh& mesh){
		assert(mesh.boundary_triangles.size()==star.boundary_triangles.size());
		for(std::size_t i=0;i<star.boundary_triangles.size();++i){
			assert(mesh.boundary_triangles[i].nodes==star.boundary_triangles[i].nodes);
			assert(mesh.boundary_triangles[i].boundary_label==star.boundary_triangles[i].boundary_label);
		}
		double volume=0.0;
		for(const auto& cell:mesh.cells)
			volume+=iga::EvaluateNativeTetGeometry(mesh,cell).determinant/6.0;
		assert(std::abs(original_volume-volume)<1e-16);
		using Face=std::array<std::uint32_t,3>;
		std::map<Face,unsigned> uses;
		for(const auto& cell:mesh.cells)
			for(std::size_t opposite=0;opposite<4;++opposite){
				Face face{};std::size_t entry=0;
				for(std::size_t local=0;local<4;++local)
					if(local!=opposite)face[entry++]=cell.nodes[local];
				std::sort(face.begin(),face.end());++uses[face];
			}
		std::set<Face> boundary;
		for(const auto& item:uses){
			assert(item.second==1||item.second==2);
			if(item.second==1)boundary.insert(item.first);
		}
		std::set<Face> labelled;
		for(const auto& triangle:mesh.boundary_triangles){
			auto face=triangle.nodes;std::sort(face.begin(),face.end());
			labelled.insert(face);
		}
		assert(boundary==labelled);
		return volume;
	};
	const double refined_volume=verify(refined);
	verify(two_shells);
	bool rejected=false;
	try{iga::RefineNativeTetStarRadialLayers(star,center,{0.6,0.4});}
	catch(const std::invalid_argument&){rejected=true;}
	assert(rejected);
	std::cout<<"native LV radial tetra refinement passed cells="<<refined.cells.size()
		<<" two_shell_cells="<<two_shells.cells.size()
		<<" volume_m3="<<refined_volume<<'\n';
}

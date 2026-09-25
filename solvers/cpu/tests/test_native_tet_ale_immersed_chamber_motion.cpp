#include "IdealizedLeftVentricleFixture.hpp"
#include "NativeTetAleConservation.hpp"
#include "NativeTetAleKinematics.hpp"
#include "NativeTetAleMeshMotion.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

namespace {

iga::NativeTetMesh StarMesh(const iga::RawSurfaceSoup& surface)
{
	iga::NativeTetMesh mesh;mesh.points=surface.vertices;
	const auto center=static_cast<std::uint32_t>(mesh.points.size());
	mesh.points.push_back({{0.,0.,.045}});std::uint64_t id=1;
	for(const auto& triangle:surface.triangles){
		auto nodes=std::array<std::uint32_t,4>{{center,static_cast<std::uint32_t>(triangle.indices[0]),
			static_cast<std::uint32_t>(triangle.indices[1]),static_cast<std::uint32_t>(triangle.indices[2])}};
		if(iga::NativeAleDeterminant(mesh.points[nodes[0]],mesh.points[nodes[1]],
			mesh.points[nodes[2]],mesh.points[nodes[3]])<0.)std::swap(nodes[2],nodes[3]);
		mesh.cells.push_back({id,nodes});mesh.boundary_triangles.push_back({id,
			{{static_cast<std::uint32_t>(triangle.indices[0]),static_cast<std::uint32_t>(triangle.indices[1]),
			static_cast<std::uint32_t>(triangle.indices[2])}},static_cast<int>(triangle.boundary_id)});++id;
	}
	for(const auto& cell:mesh.cells)assert(iga::EvaluateNativeTetGeometry(mesh,cell).determinant>0.);
	return mesh;
}

double Volume(const iga::NativeTetMesh& mesh)
{
	double result=0.;
	for(const auto& cell:mesh.cells)
		result+=iga::EvaluateNativeTetGeometry(mesh,cell).determinant/6.;
	return result;
}

}

int main()
{
	const auto motion=iga::IdealizedLeftVentricleFixture::Motion();
	const auto reference=StarMesh(motion.Frames().front().surface);
	assert(reference.cells.size()==motion.Frames().front().surface.triangles.size());
	const double reference_surface_volume=motion.Evaluate(0.,0.,.05).Surface().Diagnostics().volume_m3;
	assert(std::abs(Volume(reference)-reference_surface_volume)<1e-16);
	iga::NativeTetAleKinematics kinematics(reference.points,reference.cells);auto previous=reference;
	constexpr int steps=16;const double dt=iga::IdealizedLeftVentricleFixture::PeriodS/steps;
	double maximum_volume_relative_error=0.,maximum_position_error=0.,maximum_velocity_error=0.;
	double maximum_gcl_relative_error=0.;
	for(int step=1;step<=steps;++step){
		const double start=(step-1)*dt,end=step*dt;
		const auto evaluated=motion.Evaluate(end,start,end);
		std::map<std::uint32_t,std::array<double,3>> boundary_displacement;
		for(std::uint32_t node=0;node<motion.Frames().front().surface.vertices.size();++node){
			for(int component=0;component<3;++component)
				boundary_displacement[node][component]=evaluated.SourceVerticesM()[node][component]
					-reference.points[node][component];
		}
		const auto harmonic=iga::SolveNativeTetAleHarmonicMotion(reference,boundary_displacement);
		assert(harmonic.maximum_free_residual<1e-12);
		const auto& trial=kinematics.BeginTrial(harmonic.displacement_m,end);
		const auto current=iga::BuildNativeTetAleCurrentMesh(reference,trial);
		for(std::size_t node=0;node<evaluated.SourceVerticesM().size();++node)
			for(int component=0;component<3;++component){
				maximum_position_error=std::max(maximum_position_error,
					std::abs(current.points[node][component]-evaluated.SourceVerticesM()[node][component]));
				maximum_velocity_error=std::max(maximum_velocity_error,
					std::abs(trial.mesh_velocity_m_s[node][component]
						-evaluated.SourceVertexVelocitiesMPerS()[node][component]));
			}
		const double immersed_volume=evaluated.Surface().Diagnostics().volume_m3,ale_volume=Volume(current);
		maximum_volume_relative_error=std::max(maximum_volume_relative_error,
			std::abs(ale_volume-immersed_volume)/immersed_volume);
		const auto conservation=iga::EvaluateNativeTetAleConservation(previous,current,
			trial.mesh_velocity_m_s,trial.mesh_velocity_m_s,trial.mesh_velocity_m_s,dt);
		maximum_gcl_relative_error=std::max(maximum_gcl_relative_error,
			std::abs(conservation.gcl_residual_m3_s)/std::max(1e-30,std::abs(conservation.volume_rate_m3_s)));
		kinematics.CommitTrial();previous=current;
	}
	assert(maximum_position_error<1e-15&&maximum_velocity_error<1e-14);
	assert(maximum_volume_relative_error<2e-14&&maximum_gcl_relative_error<2e-12);
	assert(std::abs(Volume(previous)-reference_surface_volume)<1e-16);
	std::cout<<"native ALE/immersed shared chamber motion passed volume_relative_error="
		<<maximum_volume_relative_error<<" position_error_m="<<maximum_position_error
		<<" velocity_error_m_s="<<maximum_velocity_error<<" gcl_relative_error="
		<<maximum_gcl_relative_error<<'\n';return 0;
}

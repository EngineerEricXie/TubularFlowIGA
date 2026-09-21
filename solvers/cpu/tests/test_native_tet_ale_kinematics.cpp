#include "NativeTetAleKinematics.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
	const std::vector<std::array<double,3>> points{{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}};
	iga::NativeTetCell cell;cell.nodes={{0,1,2,3}};
	iga::NativeTetAleKinematics state(points,{cell});
	std::vector<std::array<double,3>> displacement(4,{{0,0,0}});
	const auto& zero=state.BeginTrial(displacement,0.1);
	assert(std::abs(zero.quality.minimum_determinant_ratio-1.0)<1e-14);
	for(const auto& velocity:zero.mesh_velocity_m_s)
		for(const double value:velocity) assert(value==0.0);
	state.CommitTrial();assert(state.CommittedTime()==0.1);
	for(auto& value:displacement) value={{0.2,-0.1,0.3}};
	const auto& translation=state.BeginTrial(displacement,0.2);
	assert(std::abs(translation.quality.minimum_determinant_ratio-1.0)<1e-14);
	for(const auto& velocity:translation.mesh_velocity_m_s) {
		assert(std::abs(velocity[0]-2.0)<1e-14);
		assert(std::abs(velocity[1]+1.0)<1e-14);
		assert(std::abs(velocity[2]-3.0)<1e-14);
	}
	state.RejectTrial();assert(state.CommittedTime()==0.1&&!state.HasTrial());
	const double angle=0.7,c=std::cos(angle),s=std::sin(angle);
	for(std::size_t i=0;i<points.size();++i) displacement[i]={{
		c*points[i][0]-s*points[i][1]-points[i][0],
		s*points[i][0]+c*points[i][1]-points[i][1],0}};
	const auto& rotation=state.BeginTrial(displacement,0.2);
	assert(std::abs(rotation.quality.minimum_determinant_ratio-1.0)<1e-13);
	state.RejectTrial();
	for(std::size_t i=0;i<points.size();++i)
		for(int component=0;component<3;++component)
			displacement[i][component]=0.1*points[i][component];
	const auto& expansion=state.BeginTrial(displacement,0.2);
	assert(std::abs(expansion.quality.minimum_determinant_ratio-1.331)<1e-13);
	iga::NativeTetMesh reference_mesh;
	reference_mesh.points=points;reference_mesh.cells={cell};
	const auto current_mesh=iga::BuildNativeTetAleCurrentMesh(reference_mesh,expansion);
	const auto current_geometry=iga::EvaluateNativeTetGeometry(current_mesh,current_mesh.cells[0]);
	assert(std::abs(current_geometry.determinant-1.331)<1e-13);
	assert(current_mesh.cells[0].nodes==reference_mesh.cells[0].nodes);
	state.RejectTrial();
	displacement.assign(4,{{0,0,0}});displacement[1]={{-2,0,0}};
	try { state.BeginTrial(displacement,0.2);assert(false); }
	catch(const std::runtime_error&) {}
	assert(!state.HasTrial()&&state.CommittedTime()==0.1);
	for(const auto& value:state.CommittedDisplacement())
		for(const double component:value) assert(component==0.0);
	std::cout<<"native tetrahedral ALE kinematics tests passed\n";
	return 0;
}

#include "NativeTetAleBackflow.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={{1,{{0,1,2,3}}}};
	mesh.boundary_triangles={{1,{{1,2,3}},2}};
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const auto owners=iga::NativeTetAleBoundaryFaceOwners(mesh);
	assert(owners.size()==1&&owners[0]==0);
	std::array<double,30> velocity{};
	for(std::size_t node=0;node<10;++node)
		for(int component=0;component<3;++component)
			velocity[3*node+component]=-0.2+0.003*node+0.001*component;
	const std::array<std::array<double,3>,4> grid{};
	const double density=1050.0,beta=0.5;
	const auto system=iga::BuildNativeTetAleBackflowFaceSystem(mesh,topology,0,
		mesh.boundary_triangles[0],velocity,grid,density,beta);
	assert(system.inward_relative_flux_m3_s>0.0);
	double energy=0.0;
	for(std::size_t dof=0;dof<30;++dof)energy+=system.residual[dof]*velocity[dof];
	assert(energy>0.0);
	const double epsilon=1e-7;
	double error_squared=0.0,reference_squared=0.0;
	for(std::size_t column=0;column<30;++column){
		auto plus=velocity,minus=velocity;
		plus[column]+=epsilon;minus[column]-=epsilon;
		const auto positive=iga::BuildNativeTetAleBackflowFaceSystem(mesh,topology,0,
			mesh.boundary_triangles[0],plus,grid,density,beta);
		const auto negative=iga::BuildNativeTetAleBackflowFaceSystem(mesh,topology,0,
			mesh.boundary_triangles[0],minus,grid,density,beta);
		for(std::size_t row=0;row<30;++row){
			const double finite_difference=(positive.residual[row]-negative.residual[row])
				/(2.0*epsilon);
			const double error=system.jacobian[30*row+column]-finite_difference;
			error_squared+=error*error;
			reference_squared+=finite_difference*finite_difference;
		}
	}
	assert(std::sqrt(error_squared/reference_squared)<1e-8);
	for(double& value:velocity)value=-value;
	const auto outward=iga::BuildNativeTetAleBackflowFaceSystem(mesh,topology,0,
		mesh.boundary_triangles[0],velocity,grid,density,beta);
	assert(outward.inward_relative_flux_m3_s==0.0);
	for(double value:outward.residual)assert(value==0.0);
	for(double value:outward.jacobian)assert(value==0.0);
	std::cout<<"native ALE optional backflow energy/Jacobian tests passed\n";
}

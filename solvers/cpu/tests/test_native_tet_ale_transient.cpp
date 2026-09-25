#include "NativeTetAleTransient.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0.1,-0.05,0.02}},{{1.1,-0.05,0.02}},
		{{0.1,0.95,0.02}},{{0.1,-0.05,1.02}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};
	using System=iga::NativeTaylorHoodElementSystem;
	std::array<double,System::dofs> current{},committed{};
	for(std::size_t index=0;index<System::dofs;++index) {
		current[index]=0.01*std::sin(static_cast<double>(index+1));
		committed[index]=0.008*std::cos(static_cast<double>(2*index+1));
	}
	const std::array<std::array<double,3>,4> grid{{{{0.2,-0.1,0.04}},
		{{0.2,-0.1,0.04}},{{0.2,-0.1,0.04}},{{0.2,-0.1,0.04}}}};
	const iga::NativeNavierStokesParameters parameters{997.0,0.0035};
	const double dt=0.025;
	const auto system=iga::BuildNativeTaylorHoodAleTransientElement(mesh,mesh.cells[0],
		current,committed,grid,parameters,dt);
	const double epsilon=1e-7;
	double difference_squared=0.0,reference_squared=0.0;
	for(std::size_t column=0;column<System::dofs;++column) {
		auto plus=current,minus=current;plus[column]+=epsilon;minus[column]-=epsilon;
		const auto plus_system=iga::BuildNativeTaylorHoodAleTransientElement(
			mesh,mesh.cells[0],plus,committed,grid,parameters,dt);
		const auto minus_system=iga::BuildNativeTaylorHoodAleTransientElement(
			mesh,mesh.cells[0],minus,committed,grid,parameters,dt);
		for(std::size_t row=0;row<System::dofs;++row) {
			const double finite_difference=(plus_system.residual[row]
				-minus_system.residual[row])/(2.0*epsilon);
			const double analytic=system.jacobian[row*System::dofs+column];
			difference_squared+=(analytic-finite_difference)*(analytic-finite_difference);
			reference_squared+=finite_difference*finite_difference;
		}
	}
	assert(std::sqrt(difference_squared/reference_squared)<2e-9);

	std::array<double,System::dofs> uniform{};
	for(std::size_t node=0;node<10;++node) {
		uniform[3*node]=0.2;uniform[3*node+1]=-0.1;uniform[3*node+2]=0.04;
	}
	const auto preserved=iga::BuildNativeTaylorHoodAleTransientElement(mesh,mesh.cells[0],
		uniform,uniform,grid,parameters,dt);
	for(const auto value:preserved.residual) assert(std::abs(value)<1e-11);
	std::cout<<"native tetrahedral ALE transient element tests passed\n";
	return 0;
}

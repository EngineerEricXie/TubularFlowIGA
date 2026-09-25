#include "NativeTetMixedHyperelasticSolid.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
iga::NativeTetMesh Mesh(){iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};return mesh;}
template<class Function>void Reject(Function&& function){bool rejected=false;try{function();}
	catch(const std::exception&){rejected=true;}assert(rejected);}
}

int main()
{
	const auto mesh=Mesh();const auto material=iga::NativeTetMixedSolidMaterialFromYoungPoisson(1200.,.4999,1000.);
	std::array<double,12> displacement{};std::array<double,4> pressure{};
	for(std::size_t i=0;i<displacement.size();++i)displacement[i]=.003*std::sin(double(i+1));
	for(std::size_t i=0;i<pressure.size();++i)pressure[i]=.2*std::cos(double(i+1));
	const auto system=iga::BuildNativeTetMixedHyperelasticSolidElement(
		mesh,mesh.cells[0],displacement,pressure,material);
	const double epsilon=1e-7;double difference2=0.,reference2=0.;
	for(std::size_t column=0;column<16;++column){
		auto up=displacement,um=displacement;auto pp=pressure,pm=pressure;
		if(column<12){up[column]+=epsilon;um[column]-=epsilon;}
		else{pp[column-12]+=epsilon;pm[column-12]-=epsilon;}
		const auto plus=iga::BuildNativeTetMixedHyperelasticSolidElement(mesh,mesh.cells[0],up,pp,material);
		const auto minus=iga::BuildNativeTetMixedHyperelasticSolidElement(mesh,mesh.cells[0],um,pm,material);
		for(std::size_t row=0;row<16;++row){
			const double finite=(plus.residual[row]-minus.residual[row])/(2.*epsilon);
			const double analytic=system.tangent[row*16+column];difference2+=(analytic-finite)*(analytic-finite);
			reference2+=finite*finite;
		}
	}
	const double tangent_error=std::sqrt(difference2/reference2);assert(tangent_error<3e-8);
	for(std::size_t row=0;row<16;++row)for(std::size_t column=0;column<16;++column)
		assert(std::abs(system.tangent[row*16+column]-system.tangent[column*16+row])<2e-8);
	std::array<double,12> zero_u{};std::array<double,4> zero_p{};
	const auto zero=iga::BuildNativeTetMixedHyperelasticSolidElement(mesh,mesh.cells[0],zero_u,zero_p,material);
	for(double value:zero.residual)assert(std::abs(value)<1e-12);
	assert(std::abs(zero.mixed_potential_j)<1e-12);
	Reject([&]{auto invalid=material;invalid.bulk_modulus_pa=0.;
		(void)iga::BuildNativeTetMixedHyperelasticSolidElement(mesh,mesh.cells[0],zero_u,zero_p,invalid);});
	std::cout<<"native mixed tetrahedral hyperelastic element passed tangent_error="<<tangent_error<<'\n';
	return 0;
}

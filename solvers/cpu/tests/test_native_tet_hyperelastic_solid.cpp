#include "NativeTetHyperelasticSolid.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

iga::NativeTetMesh Mesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};
	return mesh;
}

bool Close(double a,double b,double tolerance=1e-11)
{
	return std::abs(a-b)<=tolerance*std::max({1.0,std::abs(a),std::abs(b)});
}

template<class Function> void Reject(Function&& function)
{
	bool rejected=false;try{function();}catch(const std::exception&){rejected=true;}assert(rejected);
}

std::array<double,12> AffineDisplacement(const iga::NativeTetMesh& mesh,
	const iga::NativeTetSolidMatrix3& f,const std::array<double,3>& translation)
{
	std::array<double,12> result{};
	for(std::size_t node=0;node<4;++node) for(int i=0;i<3;++i) {
		result[3*node+i]=translation[i]-mesh.points[node][i];
		for(int j=0;j<3;++j) result[3*node+i]+=f[i][j]*mesh.points[node][j];
	}
	return result;
}

} // namespace

int main()
{
	const auto mesh=Mesh();const auto cell=mesh.cells[0];
	const iga::NativeTetSolidMaterial material{1200.0,0.3,1000.0};
	std::array<double,12> zero{};
	const auto reference=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,zero,material);
	assert(Close(reference.reference_volume_m3,1.0/6.0));
	assert(Close(reference.deformation_jacobian,1.0));
	assert(Close(reference.strain_energy_j,0.0));
	for(double value:reference.residual_n) assert(Close(value,0.0));

	iga::NativeTetSolidMatrix3 identity{};for(int i=0;i<3;++i)identity[i][i]=1.0;
	const auto translated=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,
		AffineDisplacement(mesh,identity,{{2.0,-3.0,4.0}}),material);
	assert(Close(translated.strain_energy_j,0.0));
	for(double value:translated.residual_n) assert(Close(value,0.0));

	const double angle=0.73,c=std::cos(angle),s=std::sin(angle);
	iga::NativeTetSolidMatrix3 rotation{{{{c,-s,0}},{{s,c,0}},{{0,0,1}}}};
	const auto rotated=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,
		AffineDisplacement(mesh,rotation,{{0.2,-0.4,0.7}}),material);
	assert(std::abs(rotated.strain_energy_j)<1e-10);
	for(double value:rotated.residual_n) assert(std::abs(value)<1e-9);

	iga::NativeTetSolidMatrix3 stretch{{{{1.12,0.08,0.0}},{{0.02,0.94,0.03}},{{0.0,0.01,1.05}}}};
	const auto displacement=AffineDisplacement(mesh,stretch,{{0.1,-0.2,0.3}});
	const auto deformed=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,displacement,material);
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)assert(Close(deformed.deformation_gradient[i][j],stretch[i][j]));
	assert(deformed.deformation_jacobian>0.0&&deformed.strain_energy_j>0.0);
	double force[3]={0,0,0},moment[3]={0,0,0};
	for(std::size_t node=0;node<4;++node){
		std::array<double,3> current=mesh.points[node];
		for(int i=0;i<3;++i)current[i]+=displacement[3*node+i];
		for(int i=0;i<3;++i)force[i]+=deformed.residual_n[3*node+i];
		moment[0]+=current[1]*deformed.residual_n[3*node+2]-current[2]*deformed.residual_n[3*node+1];
		moment[1]+=current[2]*deformed.residual_n[3*node]-current[0]*deformed.residual_n[3*node+2];
		moment[2]+=current[0]*deformed.residual_n[3*node+1]-current[1]*deformed.residual_n[3*node];
	}
	for(double value:force)assert(std::abs(value)<1e-10);
	for(double value:moment)assert(std::abs(value)<1e-10);

	const double epsilon=1e-7;
	double difference2=0.0,reference2=0.0;
	for(std::size_t column=0;column<12;++column){
		auto plus=displacement,minus=displacement;plus[column]+=epsilon;minus[column]-=epsilon;
		const auto rp=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,plus,material);
		const auto rm=iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,minus,material);
		for(std::size_t row=0;row<12;++row){
			const double finite=(rp.residual_n[row]-rm.residual_n[row])/(2*epsilon);
			const double analytic=deformed.tangent_n_m[row*12+column];
			difference2+=(analytic-finite)*(analytic-finite);reference2+=finite*finite;
		}
	}
	assert(std::sqrt(difference2/reference2)<2e-9);
	for(std::size_t row=0;row<12;++row)for(std::size_t column=0;column<12;++column)
		assert(Close(deformed.tangent_n_m[row*12+column],deformed.tangent_n_m[column*12+row],1e-10));
	for(int component=0;component<3;++component){
		double total=0.0;for(std::size_t a=0;a<4;++a)for(std::size_t b=0;b<4;++b)
			total+=deformed.mass_kg[(3*a+component)*12+3*b+component];
		assert(Close(total,material.density_kg_m3/6.0));
	}

	Reject([&]{auto invalid=material;invalid.poisson_ratio=0.49;
		(void)iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,zero,invalid);});
	iga::NativeTetSolidMatrix3 inverted{{{{-1,0,0}},{{0,1,0}},{{0,0,1}}}};
	Reject([&]{(void)iga::BuildNativeTetHyperelasticSolidElement(mesh,cell,
		AffineDisplacement(mesh,inverted,{{0,0,0}}),material);});
	std::cout<<"native tetrahedral hyperelastic solid element tests passed tangent_error="
		<<std::sqrt(difference2/reference2)<<'\n';
	return 0;
}

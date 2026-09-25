#include "NativeTetSolidPrestress.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {
iga::NativeTetMesh Mesh(){iga::NativeTetMesh mesh;mesh.points={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	mesh.cells={iga::NativeTetCell{1,{{0,1,2,3}}}};return mesh;}
template<class Function>void Reject(Function&& function){bool rejected=false;
	try{function();}catch(const std::exception&){rejected=true;}assert(rejected);}
}

int main()
{
	const auto reference=Mesh();const iga::NativeTetSolidMaterial material{1500.,.3,980.};
	const iga::NativeTetSolidMatrix3 target_f{{{{1.04,0,0}},{{0,.98,0}},{{0,0,1.02}}}};
	std::array<double,12> target{};
	for(std::size_t node=0;node<4;++node)for(int i=0;i<3;++i){
		target[3*node+i]=-reference.points[node][i];
		for(int j=0;j<3;++j)target[3*node+i]+=target_f[i][j]*reference.points[node][j];
	}
	const auto target_element=iga::BuildNativeTetHyperelasticSolidElement(
		reference,reference.cells[0],target,material);
	const std::map<std::size_t,double> supports{{0,0},{1,0},{2,0},{4,0},{5,0},{6,0},{8,0},{9,0},{10,0}};
	std::vector<double> external(12,0.0);
	for(std::size_t dof=0;dof<12;++dof)if(!supports.count(dof))external[dof]=target_element.residual_n[dof];
	auto loaded=reference;
	for(std::size_t node=0;node<4;++node)for(int component=0;component<3;++component)
		loaded.points[node][component]+=target[3*node+component];
	iga::NativeTetSolidPrestressOptions options;options.reference_update_relaxation=.7;
	options.absolute_reconstruction_tolerance_m=2e-10;options.forward_options.relative_tolerance=1e-11;
	options.forward_options.absolute_tolerance_n=1e-12;
	const auto result=iga::EstimateNativeTetSolidUnloadedReference(
		loaded,material,external,supports,options);
	assert(result.converged&&result.inverse_iterations>0&&result.maximum_reconstruction_error_m<=2e-10);
	double reference_error=0.;for(std::size_t node=0;node<4;++node)for(int component=0;component<3;++component)
		reference_error=std::max(reference_error,std::abs(result.unloaded_reference_mesh.points[node][component]
			-reference.points[node][component]));
	assert(reference_error<2e-8&&result.minimum_loaded_deformation_jacobian>0.0);
	assert(result.loaded_first_piola_pa.size()==1);
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)
		assert(std::abs(result.loaded_first_piola_pa[0][i][j]-target_element.first_piola_pa[i][j])<2e-5);
	for(std::size_t dof=0;dof<12;++dof)if(!supports.count(dof))
		assert(std::abs(result.equilibrating_internal_force_n[dof]-external[dof])<1e-8);
	const auto unloaded_control=iga::EstimateNativeTetSolidUnloadedReference(
		reference,material,std::vector<double>(12,0.0),supports,options);
	assert(unloaded_control.converged&&unloaded_control.inverse_iterations==0
		&&unloaded_control.maximum_reconstruction_error_m==0.0);
	for(const auto& stress:unloaded_control.loaded_first_piola_pa)
		for(const auto& row:stress)for(double value:row)assert(std::abs(value)<1e-12);
	Reject([&]{auto invalid=supports;invalid[0]=1e-3;
		(void)iga::EstimateNativeTetSolidUnloadedReference(loaded,material,external,invalid,options);});
	std::cout<<"native tetrahedral solid inverse prestress passed iterations="<<result.inverse_iterations
		<<" reconstruction_error_m="<<result.maximum_reconstruction_error_m
		<<" reference_error_m="<<reference_error<<'\n';
	return 0;
}

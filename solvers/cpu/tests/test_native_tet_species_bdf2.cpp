#include "NativeTetMovingSpeciesTransport.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh TwoCellMesh()
{
	iga::NativeTetMesh mesh;
	mesh.points={{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{0,0,-1}}};
	mesh.cells={{1,{{0,1,2,3}}},{2,{{0,2,1,4}}}};
	using Face=std::array<std::uint32_t,3>;
	std::map<Face,unsigned> uses;
	for(const auto& cell:mesh.cells)
		for(std::size_t opposite=0;opposite<4;++opposite){
			Face face{};std::size_t entry=0;
			for(std::size_t local=0;local<4;++local)
				if(local!=opposite)face[entry++]=cell.nodes[local];
			std::sort(face.begin(),face.end());++uses[face];
		}
	std::uint64_t id=1;
	for(const auto& item:uses)
		if(item.second==1)mesh.boundary_triangles.push_back({id++,item.first,1});
	return mesh;
}

std::vector<Vector> ZeroP2(const iga::NativeTetMesh& mesh)
{
	return std::vector<Vector>(mesh.points.size()
		+iga::BuildNativeTaylorHoodTopology(mesh).edges.size(),Vector{{0,0,0}});
}

// Closed mesh, no flow or diffusion: a uniform concentration decays exactly as
// the scalar ODE dc/dt = -k c, so the error isolates the time integrator.
double DecayError(iga::NativeTetSpeciesTimeIntegration scheme,double dt)
{
	const auto mesh=TwoCellMesh();
	const std::vector<Vector> grid(mesh.points.size(),Vector{{0,0,0}});
	const double rate=1.,final_time=1.;
	const int steps=static_cast<int>(std::lround(final_time/dt));
	std::vector<double> older,previous(mesh.points.size(),1.);
	for(int step=0;step<steps;++step){
		iga::NativeTetMovingSpeciesStepOptions options;
		if(scheme==iga::NativeTetSpeciesTimeIntegration::Bdf2&&step>0){
			options.time_integration=scheme;
			options.older_concentration_mol_m3=older;
		}
		const auto result=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),
			grid,previous,{},0.,0.,dt,false,rate,{},{},{},options);
		assert(std::abs(result.balance_defect_mol_s)<=1e-12);
		older=previous;
		previous=result.concentration_mol_m3;
	}
	return std::abs(previous.front()-std::exp(-rate*final_time));
}

bool Rejected(const std::function<void()>& call)
{
	try{call();}catch(const std::invalid_argument&){return true;}
	return false;
}

}

int main()
{
	using Scheme=iga::NativeTetSpeciesTimeIntegration;
	const std::array<double,3> dts{{0.1,0.05,0.025}};
	std::array<double,3> backward{},bdf2{};
	for(std::size_t level=0;level<dts.size();++level){
		backward[level]=DecayError(Scheme::BackwardEuler,dts[level]);
		bdf2[level]=DecayError(Scheme::Bdf2,dts[level]);
	}
	for(std::size_t level=1;level<dts.size();++level){
		const double order_be=std::log2(backward[level-1]/backward[level]);
		const double order_bdf2=std::log2(bdf2[level-1]/bdf2[level]);
		std::cout<<"dt="<<dts[level]<<" order_be="<<order_be<<" order_bdf2="<<order_bdf2<<'\n';
		assert(order_be>0.9&&order_be<1.1);
		assert(order_bdf2>1.85&&order_bdf2<2.15);
	}
	assert(bdf2.back()<0.1*backward.back());

	// Budget closes for BDF2 with a localized source and diffusion.
	const auto mesh=TwoCellMesh();
	const std::vector<Vector> grid(mesh.points.size(),Vector{{0,0,0}});
	const std::vector<double> c0={0.,0.2,0.4,0.1,0.3},c1={0.1,0.25,0.35,0.15,0.2};
	iga::NativeTetMovingSpeciesStepOptions options;
	options.time_integration=Scheme::Bdf2;
	options.older_concentration_mol_m3=c0;
	const auto step=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),grid,
		c1,{},0.7,0.,0.05,false,0.3,{},{},{1.,0.},options);
	assert(std::abs(step.balance_defect_mol_s)<=1e-12);
	assert(std::abs(step.inventory_rate_mol_s-(1.5*step.current_inventory_mol
		-2.*step.previous_inventory_mol+0.5*step.older_inventory_mol)/0.05)<=1e-12);

	// Invalid uses are rejected.
	iga::NativeTetMovingSpeciesStepOptions missing;
	missing.time_integration=Scheme::Bdf2;
	assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),grid,
		c1,{},0.,0.,0.1,false,0.,{},{},{},missing);}));
	iga::NativeTetMovingSpeciesStepOptions stray;
	stray.older_concentration_mol_m3=c0;
	assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),grid,
		c1,{},0.,0.,0.1,false,0.,{},{},{},stray);}));
	auto moved=mesh;
	for(auto& point:moved.points)point[0]+=0.01;
	const std::vector<Vector> translation(mesh.points.size(),Vector{{0.1,0,0}});
	assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,moved,ZeroP2(moved),
		translation,c1,{},0.,0.,0.1,false,0.,{},{},{},options);}));
	const std::vector<double> negative={0.,-1e-3,0.4,0.1,0.3};
	assert(Rejected([&]{iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),grid,
		negative,{},0.,0.,0.1);}));
	iga::NativeTetMovingSpeciesStepOptions signed_history;
	signed_history.allow_negative_concentration=true;
	const auto accepted=iga::SolveNativeTetMovingSpeciesDenseStep(mesh,mesh,ZeroP2(mesh),grid,
		negative,{},0.,0.,0.1,false,0.,{},{},{},signed_history);
	assert(std::abs(accepted.balance_defect_mol_s)<=1e-12);
	std::cout<<"native Tet species BDF2: PASS\n";
	return 0;
}

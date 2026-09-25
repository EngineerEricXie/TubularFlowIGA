#include "NativeTetMovingSpeciesCheckpoint.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Vector=std::array<double,3>;

iga::NativeTetMesh SmallMesh()
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

std::vector<Vector> P2VertexField(const iga::NativeTetMesh& mesh,
	const std::vector<Vector>& vertices)
{
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	std::vector<Vector> result=vertices;
	for(const auto& edge:topology.edges){
		Vector value{};
		for(int component=0;component<3;++component)
			value[component]=0.5*(vertices[edge[0]][component]
				+vertices[edge[1]][component]);
		result.push_back(value);
	}
	return result;
}

void Close(double actual,double expected,double tolerance=1e-10)
{
	assert(std::isfinite(actual));
	assert(std::abs(actual-expected)<=tolerance*std::max(1.,std::abs(expected)));
}

template<class Function> void Reject(Function action)
{
	bool rejected=false;
	try{action();}catch(const std::exception&){rejected=true;}
	assert(rejected);
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank=0,ranks=1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const auto reference=SmallMesh();
	auto state=iga::InitializeNativeTetMovingSpeciesState(reference,
		std::vector<double>(reference.points.size(),2.),0.);
	auto translated=reference;
	for(auto& point:translated.points)point[0]+=0.1;
	const std::vector<Vector> translation(reference.points.size(),Vector{{1,0,0}});
	const auto first=iga::AdvanceNativeTetMovingSpeciesState(state,translated,
		P2VertexField(translated,translation),translation,{},0.,0.1);
	Close(first.step.balance_defect_mol_s,0.);
	assert(state.accepted_steps==1);Close(state.time_s,0.1);
	const auto checkpoint=iga::CaptureNativeTetMovingSpeciesCheckpoint(state);
	const auto bytes=iga::SerializeNativeTetMovingSpeciesCheckpoint(checkpoint);
	const auto parsed=iga::ParseNativeTetMovingSpeciesCheckpoint(bytes);
	auto restarted=iga::RestoreNativeTetMovingSpeciesCheckpoint(parsed,reference,0.);
	assert(restarted.model_identity_sha256==state.model_identity_sha256);
	assert(restarted.accepted_steps==1);
	assert(restarted.current_mesh.points==translated.points);
	assert(restarted.concentration_mol_m3==state.concentration_mol_m3);
	const auto before_failed_step=iga::SerializeNativeTetMovingSpeciesCheckpoint(
		iga::CaptureNativeTetMovingSpeciesCheckpoint(restarted));
	std::vector<Vector> extra_flow(translation.size(),Vector{{2,0,0}});
	Reject([&]{iga::AdvanceNativeTetMovingSpeciesState(restarted,translated,
		P2VertexField(translated,extra_flow),
		std::vector<Vector>(translation.size(),Vector{{0,0,0}}),{},0.,0.1);});
	assert(before_failed_step==iga::SerializeNativeTetMovingSpeciesCheckpoint(
		iga::CaptureNativeTetMovingSpeciesCheckpoint(restarted)));
	auto wrong_trial=translated;wrong_trial.cells[0].id+=10;
	Reject([&]{iga::AdvanceNativeTetMovingSpeciesState(restarted,wrong_trial,
		P2VertexField(wrong_trial,translation),
		std::vector<Vector>(translation.size(),Vector{{0,0,0}}),{},0.,0.1);});
	assert(before_failed_step==iga::SerializeNativeTetMovingSpeciesCheckpoint(
		iga::CaptureNativeTetMovingSpeciesCheckpoint(restarted)));
	auto corrupted=bytes;corrupted[corrupted.size()/2]^=1;
	Reject([&]{iga::ParseNativeTetMovingSpeciesCheckpoint(corrupted);});
	Reject([&]{iga::ParseNativeTetMovingSpeciesCheckpoint(
		std::string_view(bytes.data(),bytes.size()-1));});
	Reject([&]{iga::RestoreNativeTetMovingSpeciesCheckpoint(parsed,reference,0.01);});
	auto wrong_reference=reference;wrong_reference.points[0][0]=0.01;
	Reject([&]{iga::RestoreNativeTetMovingSpeciesCheckpoint(parsed,wrong_reference,0.);});
	auto forged=parsed;forged.concentration_mol_m3[0]+=1.;
	Reject([&]{iga::RestoreNativeTetMovingSpeciesCheckpoint(forged,reference,0.);});
	auto expanded=translated;
	for(auto& point:expanded.points)point[2]*=1.2;
	std::vector<Vector> expansion(expanded.points.size());
	for(std::size_t node=0;node<expansion.size();++node)
		expansion[node][2]=(expanded.points[node][2]-translated.points[node][2])/0.1;
	const auto flow=P2VertexField(expanded,expansion);
	const auto uninterrupted=iga::AdvanceNativeTetMovingSpeciesState(state,expanded,
		flow,expansion,{},0.,0.1);
	const auto resumed=iga::AdvanceNativeTetMovingSpeciesState(restarted,expanded,
		flow,expansion,{},0.,0.1);
	assert(state.accepted_steps==2&&restarted.accepted_steps==2);
	Close(state.time_s,0.2);Close(restarted.time_s,0.2);
	for(std::size_t node=0;node<reference.points.size();++node){
		Close(state.concentration_mol_m3[node],5./3.);
		Close(state.concentration_mol_m3[node],restarted.concentration_mol_m3[node]);
	}
	Close(uninterrupted.step.current_inventory_mol,2./3.);
	Close(uninterrupted.step.current_inventory_mol,resumed.step.current_inventory_mol);
	assert(iga::SerializeNativeTetMovingSpeciesCheckpoint(
		iga::CaptureNativeTetMovingSpeciesCheckpoint(state))
		==iga::SerializeNativeTetMovingSpeciesCheckpoint(
			iga::CaptureNativeTetMovingSpeciesCheckpoint(restarted)));
	if(rank==0)std::cout<<"native moving species checkpoint passed ranks="<<ranks
		<<" bytes="<<bytes.size()<<" inventory_mol="
		<<resumed.step.current_inventory_mol<<'\n';
	PetscFinalize();
}

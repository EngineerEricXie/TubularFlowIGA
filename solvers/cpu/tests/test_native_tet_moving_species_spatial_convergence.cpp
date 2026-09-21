#include "NativeTetMovingSpeciesPetscRuntime.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

using Face=std::array<std::uint32_t,3>;

iga::NativeTetMesh UnitCube(int divisions)
{
	if(divisions<1)throw std::invalid_argument("native species cube division is invalid");
	iga::NativeTetMesh mesh;
	const int side=divisions+1;
	const auto index=[side](int i,int j,int k){
		return static_cast<std::uint32_t>(i+side*(j+side*k));
	};
	for(int k=0;k<side;++k)
		for(int j=0;j<side;++j)
			for(int i=0;i<side;++i)
				mesh.points.push_back({{double(i)/divisions,double(j)/divisions,
					double(k)/divisions}});
	const std::array<std::array<int,4>,6> local{{
		{{0,1,3,7}},{{0,3,2,7}},{{0,2,6,7}},
		{{0,6,4,7}},{{0,4,5,7}},{{0,5,1,7}}}};
	std::map<Face,unsigned> face_uses;
	for(int k=0;k<divisions;++k)
		for(int j=0;j<divisions;++j)
			for(int i=0;i<divisions;++i){
				const std::array<std::uint32_t,8> vertices{{
					index(i,j,k),index(i+1,j,k),index(i,j+1,k),index(i+1,j+1,k),
					index(i,j,k+1),index(i+1,j,k+1),index(i,j+1,k+1),index(i+1,j+1,k+1)}};
				for(const auto& tet:local){
					iga::NativeTetCell cell;
					cell.id=mesh.cells.size()+1;
					for(std::size_t node=0;node<4;++node)
						cell.nodes[node]=vertices[tet[node]];
					iga::EvaluateNativeTetGeometry(mesh,cell);
					mesh.cells.push_back(cell);
					for(std::size_t opposite=0;opposite<4;++opposite){
						Face face{};std::size_t entry=0;
						for(std::size_t node=0;node<4;++node)
							if(node!=opposite)face[entry++]=cell.nodes[node];
						std::sort(face.begin(),face.end());
						++face_uses[face];
					}
				}
			}
	for(const auto& face:face_uses)
		if(face.second==1)
			mesh.boundary_triangles.push_back({
				mesh.boundary_triangles.size()+1,face.first,1});
		else if(face.second!=2)
			throw std::runtime_error("native species cube has nonmanifold faces");
	return mesh;
}

template<class Exact>
double RelativeL2Error(const iga::NativeTetMesh& mesh,
	const std::vector<double>& concentration,Exact exact,double reference_norm_squared)
{
	double squared_error=0.;
	const auto quadrature=iga::NativeTetDegreeFiveQuadrature();
	for(const auto& cell:mesh.cells){
		const double determinant=iga::EvaluateNativeTetGeometry(mesh,cell).determinant;
		for(const auto& point:quadrature){
			const std::array<double,4> basis{{1.-point.reference[0]
				-point.reference[1]-point.reference[2],point.reference[0],
				point.reference[1],point.reference[2]}};
			double x=0.,numerical=0.;
			for(std::size_t node=0;node<4;++node){
				x+=basis[node]*mesh.points[cell.nodes[node]][0];
				numerical+=basis[node]*concentration[cell.nodes[node]];
			}
			const double difference=numerical-exact(x);
			squared_error+=determinant*point.weight*difference*difference;
		}
	}
	return std::sqrt(squared_error/reference_norm_squared);
}

double SolveLevel(int divisions,bool monotone)
{
	const auto mesh=UnitCube(divisions);
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::vector<std::array<double,3>> fluid_velocity(
		mesh.points.size()+topology.edges.size(),{{0.,0.,0.}});
	const std::vector<std::array<double,3>> mesh_velocity(mesh.points.size(),{{0.,0.,0.}});
	const double pi=std::acos(-1.),diffusivity=0.1,dt=0.1;
	std::vector<double> initial;
	initial.reserve(mesh.points.size());
	for(const auto& point:mesh.points)initial.push_back(2.+std::cos(pi*point[0]));
	const auto result=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,mesh,
		fluid_velocity,mesh_velocity,initial,{},diffusivity,0.,dt,monotone);
	if(result.converged_reason<=0
		||std::abs(result.step.current_inventory_mol-2.)>1e-9
		||std::abs(result.step.balance_defect_mol_s)>1e-9)
		throw std::runtime_error("native species diffusion spatial balance failed");
	const double amplitude=1./(1.+diffusivity*pi*pi*dt);
	return RelativeL2Error(mesh,result.step.concentration_mol_m3,
		[=](double x){return 2.+amplitude*std::cos(pi*x);},
		amplitude*amplitude/2.);
}

struct FrontResult
{
	double relative_l2=0.;
	double minimum_mol_m3=0.;
	double inventory_mol=0.;
	double outward_flux_mol_s=0.;
	std::vector<double> concentration_mol_m3;
};

FrontResult SolveAdvectionFrontLevel(int divisions,bool monotone,bool moving=false)
{
	const auto mesh=UnitCube(divisions);
	auto current_mesh=mesh;
	const double dt=0.2,grid_speed=moving?0.1:0.;
	if(moving)
		for(auto& point:current_mesh.points)point[0]+=grid_speed*dt;
	const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::vector<std::array<double,3>> fluid_velocity(
		mesh.points.size()+topology.edges.size(),{{1.+grid_speed,0.,0.}});
	const std::vector<std::array<double,3>> mesh_velocity(mesh.points.size(),
		{{grid_speed,0.,0.}});
	const std::vector<double> initial(mesh.points.size(),0.);
	const double front_length=dt;
	const auto result=iga::SolveNativeTetMovingSpeciesPetscStep(mesh,current_mesh,
		fluid_velocity,mesh_velocity,initial,{{1,1.}},0.,0.,dt,monotone);
	if(result.converged_reason<=0
		||!(result.step.current_inventory_mol>0.)
		||std::abs(result.step.balance_defect_mol_s)>1e-9)
		throw std::runtime_error("native species pure-advection front balance failed");
	FrontResult front;
	front.minimum_mol_m3=*std::min_element(result.step.concentration_mol_m3.begin(),
		result.step.concentration_mol_m3.end());
	front.inventory_mol=result.step.current_inventory_mol;
	front.outward_flux_mol_s=result.step.outward_advective_flux_mol_s;
	front.concentration_mol_m3=result.step.concentration_mol_m3;
	const double norm_squared=front_length/2.
		*(1.-std::exp(-2./front_length));
	front.relative_l2=RelativeL2Error(current_mesh,result.step.concentration_mol_m3,
		[=](double x){return std::exp(-(x-grid_speed*dt)/front_length);},
		norm_squared);
	return front;
}

double SolveAffineExpansionLevel(int divisions,bool monotone)
{
	const auto previous_mesh=UnitCube(divisions);
	auto current_mesh=previous_mesh;
	const double dt=0.2,stretch=1.2,initial_concentration=2.;
	for(auto& point:current_mesh.points)point[0]*=stretch;
	const auto topology=iga::BuildNativeTaylorHoodTopology(previous_mesh);
	std::vector<std::array<double,3>> fluid_velocity(
		previous_mesh.points.size()+topology.edges.size(),{{0.,0.,0.}});
	std::vector<std::array<double,3>> mesh_velocity(previous_mesh.points.size(),
		{{0.,0.,0.}});
	for(std::size_t node=0;node<previous_mesh.points.size();++node){
		const double velocity=(stretch-1.)*previous_mesh.points[node][0]/dt;
		fluid_velocity[node][0]=velocity;
		mesh_velocity[node][0]=velocity;
	}
	for(std::size_t edge=0;edge<topology.edges.size();++edge){
		const auto& endpoints=topology.edges[edge];
		fluid_velocity[previous_mesh.points.size()+edge][0]=0.5*
			(mesh_velocity[endpoints[0]][0]+mesh_velocity[endpoints[1]][0]);
	}
	const std::vector<double> initial(previous_mesh.points.size(),initial_concentration);
	const auto result=iga::SolveNativeTetMovingSpeciesPetscStep(previous_mesh,
		current_mesh,fluid_velocity,mesh_velocity,initial,{},0.,0.,dt,monotone);
	const double expected=initial_concentration/stretch;
	double maximum_error=0.;
	for(const double concentration:result.step.concentration_mol_m3)
		maximum_error=std::max(maximum_error,std::abs(concentration-expected));
	if(result.converged_reason<=0
		||std::abs(result.step.previous_inventory_mol-initial_concentration)>1e-9
		||std::abs(result.step.current_inventory_mol-initial_concentration)>1e-9
		||std::abs(result.step.outward_advective_flux_mol_s)>1e-9
		||std::abs(result.step.balance_defect_mol_s)>1e-9)
		throw std::runtime_error("native species affine ALE expansion balance failed");
	return maximum_error;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int result_code=0;
	try{
		std::array<double,3> galerkin{},monotone{};
		std::array<FrontResult,3> advective_galerkin{},advective_monotone{};
		double maximum_ale_field_difference=0.,maximum_ale_inventory_difference=0.;
		double maximum_ale_flux_difference=0.,maximum_ale_error_difference=0.;
		double maximum_expansion_concentration_error=0.;
		const std::array<int,3> divisions{{2,4,8}};
		for(std::size_t level=0;level<divisions.size();++level){
			galerkin[level]=SolveLevel(divisions[level],false);
			monotone[level]=SolveLevel(divisions[level],true);
			advective_galerkin[level]=SolveAdvectionFrontLevel(divisions[level],false);
			advective_monotone[level]=SolveAdvectionFrontLevel(divisions[level],true);
			for(const bool low_order:{false,true})
				maximum_expansion_concentration_error=std::max(
					maximum_expansion_concentration_error,
					SolveAffineExpansionLevel(divisions[level],low_order));
			for(const bool low_order:{false,true}){
				const auto moving=SolveAdvectionFrontLevel(divisions[level],low_order,true);
				const auto& fixed=low_order?advective_monotone[level]
					:advective_galerkin[level];
				if(moving.concentration_mol_m3.size()!=fixed.concentration_mol_m3.size())
					throw std::runtime_error("native species ALE front field size differs");
				for(std::size_t node=0;node<fixed.concentration_mol_m3.size();++node)
					maximum_ale_field_difference=std::max(maximum_ale_field_difference,
						std::abs(moving.concentration_mol_m3[node]
							-fixed.concentration_mol_m3[node]));
				maximum_ale_inventory_difference=std::max(maximum_ale_inventory_difference,
					std::abs(moving.inventory_mol-fixed.inventory_mol));
				maximum_ale_flux_difference=std::max(maximum_ale_flux_difference,
					std::abs(moving.outward_flux_mol_s-fixed.outward_flux_mol_s));
				maximum_ale_error_difference=std::max(maximum_ale_error_difference,
					std::abs(moving.relative_l2-fixed.relative_l2));
			}
			if(!(monotone[level]>galerkin[level]))
				throw std::runtime_error("native species monotone diffusion accuracy tradeoff changed");
		}
		const double first_order=std::log2(galerkin[0]/galerkin[1]);
		const double second_order=std::log2(galerkin[1]/galerkin[2]);
		const double monotone_first_order=std::log2(monotone[0]/monotone[1]);
		const double monotone_second_order=std::log2(monotone[1]/monotone[2]);
		const double front_galerkin_first_order=std::log2(
			advective_galerkin[0].relative_l2/advective_galerkin[1].relative_l2);
		const double front_galerkin_second_order=std::log2(
			advective_galerkin[1].relative_l2/advective_galerkin[2].relative_l2);
		const double front_monotone_first_order=std::log2(
			advective_monotone[0].relative_l2/advective_monotone[1].relative_l2);
		const double front_monotone_second_order=std::log2(
			advective_monotone[1].relative_l2/advective_monotone[2].relative_l2);
		if(maximum_ale_field_difference>1e-10
			||maximum_expansion_concentration_error>1e-10
			||maximum_ale_inventory_difference>1e-10
			||maximum_ale_flux_difference>1e-10
			||maximum_ale_error_difference>1e-10
			||!(first_order>1.5&&second_order>1.5)
			||!(monotone_first_order>1.2&&monotone_second_order>1.2)
			||!(front_galerkin_first_order>1.4&&front_galerkin_second_order>1.4)
			||!(front_monotone_first_order>0.4&&front_monotone_second_order>0.4)
			||!(advective_galerkin[0].minimum_mol_m3<-0.01))
			throw std::runtime_error("native species spatial errors do not decrease");
		for(std::size_t level=0;level<divisions.size();++level)
			if(advective_monotone[level].minimum_mol_m3<-1e-12
				||!(advective_monotone[level].relative_l2
					>advective_galerkin[level].relative_l2))
				throw std::runtime_error("native species advective front loses monotone tradeoff");
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cout<<std::setprecision(17)
			<<"native moving species spatial convergence passed"
			<<" galerkin_relative_l2="<<galerkin[0]<<','<<galerkin[1]<<','<<galerkin[2]
			<<" galerkin_orders="<<first_order<<','<<second_order
			<<" monotone_relative_l2="<<monotone[0]<<','<<monotone[1]<<','<<monotone[2]
			<<" monotone_orders="<<monotone_first_order<<','<<monotone_second_order
			<<" advective_galerkin_relative_l2="
			<<advective_galerkin[0].relative_l2<<','
			<<advective_galerkin[1].relative_l2<<','
			<<advective_galerkin[2].relative_l2
			<<" advective_galerkin_orders="
			<<front_galerkin_first_order<<','<<front_galerkin_second_order
			<<" advective_galerkin_minimum_mol_m3="
			<<advective_galerkin[0].minimum_mol_m3<<','
			<<advective_galerkin[1].minimum_mol_m3<<','
			<<advective_galerkin[2].minimum_mol_m3
			<<" advective_monotone_relative_l2="
			<<advective_monotone[0].relative_l2<<','
			<<advective_monotone[1].relative_l2<<','
			<<advective_monotone[2].relative_l2
			<<" advective_monotone_orders="
			<<front_monotone_first_order<<','<<front_monotone_second_order
			<<" advective_monotone_minimum_mol_m3="
			<<advective_monotone[0].minimum_mol_m3<<','
			<<advective_monotone[1].minimum_mol_m3<<','
			<<advective_monotone[2].minimum_mol_m3
			<<" maximum_ale_field_difference_mol_m3="<<maximum_ale_field_difference
			<<" maximum_ale_inventory_difference_mol="<<maximum_ale_inventory_difference
			<<" maximum_ale_flux_difference_mol_s="<<maximum_ale_flux_difference
			<<" maximum_ale_relative_l2_difference="<<maximum_ale_error_difference
			<<" maximum_expansion_concentration_error_mol_m3="
			<<maximum_expansion_concentration_error
			<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<error.what()<<'\n';
		result_code=1;
	}
	PetscFinalize();return result_code;
}

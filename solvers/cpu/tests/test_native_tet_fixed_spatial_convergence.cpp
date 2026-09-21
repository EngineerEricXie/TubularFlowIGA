#include "NativeTetAleDenseRuntime.hpp"
#include "NativeTetAleKinematics.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace {

std::uint32_t Node(int n,int i,int j,int k)
{return static_cast<std::uint32_t>((i*(n+1)+j)*(n+1)+k);}

iga::NativeTetMesh Mesh(int n)
{
	iga::NativeTetMesh mesh;
	for(int i=0;i<=n;++i)for(int j=0;j<=n;++j)for(int k=0;k<=n;++k)
		mesh.points.push_back({{double(i)/n,double(j)/n,double(k)/n}});
	std::uint64_t id=1;
	for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){
		const auto a=Node(n,i,j,k),b=Node(n,i+1,j,k),c=Node(n,i,j+1,k),d=Node(n,i+1,j+1,k);
		const auto e=Node(n,i,j,k+1),f=Node(n,i+1,j,k+1),g=Node(n,i,j+1,k+1),h=Node(n,i+1,j+1,k+1);
		const std::array<std::array<std::uint32_t,4>,6> cells{{{{a,b,d,h}},{{a,d,c,h}},{{a,c,g,h}},
			{{a,g,e,h}},{{a,e,f,h}},{{a,f,b,h}}}};
		for(auto nodes:cells){
			if(iga::NativeAleDeterminant(mesh.points[nodes[0]],mesh.points[nodes[1]],
				mesh.points[nodes[2]],mesh.points[nodes[3]])<0.0)std::swap(nodes[1],nodes[2]);
			mesh.cells.push_back({id++,nodes});
		}
	}
	using Face=std::array<std::uint32_t,3>;std::map<Face,std::pair<int,Face>> faces;
	for(const auto& cell:mesh.cells)for(std::size_t opposite=0;opposite<4;++opposite){
		Face oriented{};std::size_t entry=0;
		for(std::size_t local=0;local<4;++local)if(local!=opposite)oriented[entry++]=cell.nodes[local];
		auto key=oriented;std::sort(key.begin(),key.end());auto& record=faces[key];
		++record.first;record.second=oriented;
	}
	id=1;for(const auto& face:faces)if(face.second.first==1)
		mesh.boundary_triangles.push_back({id++,face.second.second,0});
	return mesh;
}

std::array<double,3> Coordinate(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology,std::uint32_t node)
{
	if(node<mesh.points.size())return mesh.points[node];
	const auto& edge=topology.edges.at(node-mesh.points.size());std::array<double,3> point{};
	for(int component=0;component<3;++component)
		point[component]=.5*(mesh.points[edge[0]][component]+mesh.points[edge[1]][component]);
	return point;
}

constexpr double amplitude=1.e-4;
double Exact(double y){return amplitude*y*y*y;}

double Run(int n)
{
	const auto mesh=Mesh(n);const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	std::map<std::uint32_t,std::array<double,3>> boundary;
	for(auto node:topology.boundary_velocity_nodes.at(0))
		boundary[node]={{Exact(Coordinate(mesh,topology,node)[1]),0,0}};
	iga::NativeNavierStokesParameters parameters;parameters.density=1.;parameters.dynamic_viscosity=.004;
	parameters.body_acceleration_m_s2[0][2]=-6.*parameters.dynamic_viscosity*amplitude/parameters.density;
	const auto solved=iga::SolveNativeTetDenseSteady(mesh,
		std::vector<double>(3*velocity_nodes+mesh.points.size(),0.0),boundary,0,parameters,1e-9,12);
	assert(solved.final_free_residual_l2<=1e-9&&solved.newton_iterations>0);
	double error_squared=0.,exact_squared=0.;
	for(std::size_t cell=0;cell<mesh.cells.size();++cell){
		const auto geometry=iga::EvaluateNativeTetGeometry(mesh,mesh.cells[cell]);
		const auto dofs=iga::native_tet_ale_dense_detail::ElementDofs(mesh,topology,cell);
		for(const auto& quadrature:iga::NativeTetDegreeFiveQuadrature()){
			const auto basis=iga::EvaluateNativeTaylorHoodPhysicalBasis(geometry,quadrature.reference);
			double numerical=0.,y=0.;
			for(std::size_t local=0;local<10;++local)numerical+=basis.velocity[local]*solved.state[dofs[3*local]];
			for(std::size_t local=0;local<4;++local)y+=basis.pressure[local]*mesh.points[mesh.cells[cell].nodes[local]][1];
			const double exact=Exact(y),weight=quadrature.weight*geometry.determinant;
			error_squared+=weight*(numerical-exact)*(numerical-exact);exact_squared+=weight*exact*exact;
		}
	}
	return std::sqrt(error_squared/exact_squared);
}

}

int main()
{
	const std::array<int,3> subdivisions{{2,3,4}};std::array<double,3> error{};
	for(std::size_t level=0;level<subdivisions.size();++level){
		std::cout<<"manufactured spatial level n="<<subdivisions[level]<<std::endl;
		error[level]=Run(subdivisions[level]);
	}
	const double order12=std::log(error[0]/error[1])/std::log(1.5);
	const double order23=std::log(error[1]/error[2])/std::log(4.0/3.0);
	assert(error[2]<error[1]&&error[1]<error[0]);
	assert(order12>2.2&&order23>2.2);
	std::cout<<"native fixed-domain manufactured spatial convergence passed errors="
		<<error[0]<<','<<error[1]<<','<<error[2]<<" orders="<<order12<<','<<order23<<'\n';
	return 0;
}

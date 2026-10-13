// The distributed ALE runtime reproduces the replicated runtime: same Newton
// history and the same state to round-off, on any rank count, for a channel
// with inflow, no-slip walls, and a pressure outlet (no gauge).
#include "NativeTetAleDistributedRuntime.hpp"
#include "NativeTetAlePetscRuntime.hpp"

#include <petscksp.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

iga::NativeTetMesh ChannelMesh(int nx,int ny,int nz,double length,double height,double width)
{
	iga::NativeTetMesh mesh;
	const auto node=[&](int i,int j,int k){return static_cast<std::uint32_t>((i*(ny+1)+j)*(nz+1)+k);};
	for(int i=0;i<=nx;++i)for(int j=0;j<=ny;++j)for(int k=0;k<=nz;++k)
		mesh.points.push_back({{length*i/nx,height*j/ny,width*k/nz}});
	const auto determinant=[&](std::array<std::uint32_t,4> n){
		const auto& a=mesh.points[n[0]];const auto& b=mesh.points[n[1]];const auto& c=mesh.points[n[2]];const auto& d=mesh.points[n[3]];
		const double x[3]{b[0]-a[0],b[1]-a[1],b[2]-a[2]},y[3]{c[0]-a[0],c[1]-a[1],c[2]-a[2]},z[3]{d[0]-a[0],d[1]-a[1],d[2]-a[2]};
		return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0]);};
	std::uint64_t id=1;
	for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<nz;++k){
		const auto a=node(i,j,k),b=node(i+1,j,k),c=node(i,j+1,k),d=node(i+1,j+1,k);
		const auto e=node(i,j,k+1),f=node(i+1,j,k+1),g=node(i,j+1,k+1),h=node(i+1,j+1,k+1);
		for(auto cell:{std::array<std::uint32_t,4>{{a,b,d,h}},{{a,d,c,h}},{{a,c,g,h}},{{a,g,e,h}},{{a,e,f,h}},{{a,f,b,h}}}){
			if(determinant(cell)<0)std::swap(cell[1],cell[2]);
			mesh.cells.push_back({id++,cell});
		}
	}
	std::map<std::array<std::uint32_t,3>,std::pair<int,std::array<std::uint32_t,3>>> faces;
	for(const auto& cell:mesh.cells)for(int opposite=0;opposite<4;++opposite){
		std::array<std::uint32_t,3> face{};int entry=0;
		for(int vertex=0;vertex<4;++vertex)if(vertex!=opposite)face[entry++]=cell.nodes[vertex];
		auto key=face;std::sort(key.begin(),key.end());
		auto& use=faces.try_emplace(key,std::make_pair(0,face)).first->second;++use.first;
	}
	std::uint64_t face_id=1;
	for(const auto& face:faces){
		if(face.second.first!=1)continue;
		const auto& nodes=face.second.second;
		int label=3;
		if(mesh.points[nodes[0]][0]==0&&mesh.points[nodes[1]][0]==0&&mesh.points[nodes[2]][0]==0)label=1;
		else if(mesh.points[nodes[0]][0]==length&&mesh.points[nodes[1]][0]==length&&mesh.points[nodes[2]][0]==length)label=2;
		iga::NativeTetTriangle triangle;triangle.id=face_id++;triangle.boundary_label=label;triangle.nodes=nodes;
		mesh.boundary_triangles.push_back(triangle);
	}
	return mesh;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try{
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		const auto mesh=ChannelMesh(6,2,3,0.02,0.004,0.004);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> committed(3*velocity_nodes+mesh.points.size(),0.0);
		auto trial=committed;
		std::map<std::uint32_t,std::array<double,3>> prescribed;
		for(const auto node:topology.boundary_velocity_nodes.at(1))prescribed[node]={{0.2,0.0,0.0}};
		for(const auto node:topology.boundary_velocity_nodes.at(3))prescribed.emplace(node,std::array<double,3>{{0.0,0.0,0.0}});
		for(const auto& item:prescribed)for(int component=0;component<3;++component)trial[3*item.first+component]=item.second[component];
		iga::NativeTetAleBoundaryConditions conditions;conditions.prescribed_pressure_pa[2]=0.0;
		const iga::NativeNavierStokesParameters parameters{1000.0,0.004};
		const std::vector<std::array<double,3>> mesh_velocity(mesh.points.size(),{{0.0,0.0,0.0}});
		const double dt=0.01,tolerance=1e-10;
		const auto reference=iga::SolveNativeTetAlePetscTransient(mesh,mesh_velocity,committed,trial,prescribed,
			std::numeric_limits<std::uint32_t>::max(),parameters,dt,tolerance,12,conditions);
		const iga::NativeTetAleDistributedMesh distributed(mesh,PETSC_COMM_WORLD);
		const auto result=iga::SolveNativeTetAleDistributed(distributed,topology,mesh.points,mesh_velocity,committed,trial,
			prescribed,std::numeric_limits<std::uint32_t>::max(),parameters,dt,false,tolerance,12,conditions);
		double difference=0.0,scale=0.0;
		for(std::size_t dof=0;dof<reference.replicated_state.size();++dof){
			const double d=reference.replicated_state[dof]-result.replicated_state[dof];
			difference+=d*d;scale+=reference.replicated_state[dof]*reference.replicated_state[dof];
		}
		const double relative=std::sqrt(difference/scale);
		std::uint32_t owned_cells=0;
		for(const auto owned:distributed.partition.local_cell_owned)owned_cells+=owned;
		std::uint32_t total_owned=0;
		MPI_Allreduce(&owned_cells,&total_owned,1,MPI_UINT32_T,MPI_SUM,PETSC_COMM_WORLD);
		if(total_owned!=mesh.cells.size())throw std::runtime_error("partition does not cover every cell once");
		if(result.global_rows!=static_cast<PetscInt>(reference.replicated_state.size()))
			throw std::runtime_error("distributed row count differs from the replicated layout");
		if(result.newton_iterations!=reference.newton_iterations)
			throw std::runtime_error("distributed Newton iteration count differs");
		if(!(relative<1e-10))throw std::runtime_error("distributed state differs from the replicated runtime: "+std::to_string(relative));
		if(rank==0)std::cout<<"native ALE distributed runtime passed ranks="<<ranks<<" cells="<<mesh.cells.size()
			<<" local_cells="<<distributed.partition.local.cells.size()<<" owned_rows="<<result.owned_rows
			<<" newton="<<result.newton_iterations<<" relative_difference="<<relative<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<"native ALE distributed runtime test: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}

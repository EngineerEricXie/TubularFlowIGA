// The partitioned mesh motion, interface traction, and solid solve reproduce
// their replicated counterparts to round-off on any rank count, including a
// solid mesh carrying a vertex no cell uses.
#include "NativeTetAleDistributedRuntime.hpp"
#include "NativeTetSolidDistributedSolver.hpp"

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
		else if(mesh.points[nodes[0]][1]==height&&mesh.points[nodes[1]][1]==height&&mesh.points[nodes[2]][1]==height)label=7;
		iga::NativeTetTriangle triangle;triangle.id=face_id++;triangle.boundary_label=label;triangle.nodes=nodes;
		mesh.boundary_triangles.push_back(triangle);
	}
	return mesh;
}

double RelativeDifference(const std::vector<double>& a,const std::vector<double>& b)
{
	if(a.size()!=b.size())throw std::runtime_error("compared vectors differ in size");
	double difference=0.0,scale=0.0;
	for(std::size_t i=0;i<a.size();++i){difference+=(a[i]-b[i])*(a[i]-b[i]);scale+=a[i]*a[i];}
	return std::sqrt(difference/std::max(scale,1e-300));
}

std::vector<double> Flatten(const std::vector<std::array<double,3>>& values)
{
	std::vector<double> flat;for(const auto& value:values)for(double v:value)flat.push_back(v);return flat;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);int status=0;
	try{
		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		const double length=0.02,height=0.004,width=0.004;
		const auto mesh=ChannelMesh(10,3,3,length,height,width);
		const auto topology=iga::BuildNativeTaylorHoodTopology(mesh);
		iga::NativeTetAleDistributedMesh distributed(mesh,PETSC_COMM_WORLD);
		distributed.BindReplicatedTopology(topology);
		// Mesh motion: the top wall (label 7) bulges, every other boundary vertex stays.
		std::map<std::uint32_t,std::array<double,3>> prescribed;
		for(const auto& triangle:mesh.boundary_triangles)for(const auto node:triangle.nodes){
			const auto& p=mesh.points[node];
			const double bulge=triangle.boundary_label==7?2e-4*std::sin(M_PI*p[0]/length)*std::sin(M_PI*p[2]/width):0.0;
			prescribed[node]={{0.0,bulge,0.0}};
		}
		iga::NativeTetPetscLinearSolveDefaults motion;motion.ksp_type=KSPCG;motion.pc_type=PCJACOBI;motion.factor_package.clear();
		const auto replicated_motion=iga::SolveNativeTetAleHarmonicMotion(mesh,prescribed,iga::MakeNativeTetPetscLinearSolve("fsi_mesh_",motion));
		const auto distributed_motion=iga::SolveNativeTetAleHarmonicMotionDistributed(distributed,mesh,prescribed);
		const double motion_difference=RelativeDifference(Flatten(replicated_motion.displacement_m),Flatten(distributed_motion.displacement_m));
		if(replicated_motion.free_nodes!=distributed_motion.free_nodes||!(motion_difference<1e-9))
			throw std::runtime_error("distributed mesh motion differs: "+std::to_string(motion_difference));
		// Traction of a smooth state on the moved mesh.
		std::vector<std::array<double,3>> current(mesh.points.size());
		for(std::size_t node=0;node<current.size();++node)for(int c=0;c<3;++c)
			current[node][c]=mesh.points[node][c]+distributed_motion.displacement_m[node][c];
		iga::NativeTetMesh current_mesh=mesh;current_mesh.points=current;
		const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
		std::vector<double> state(3*velocity_nodes+mesh.points.size());
		const auto velocity_position=[&](std::size_t node){
			if(node<mesh.points.size())return current[node];
			const auto& edge=topology.edges[node-mesh.points.size()];
			std::array<double,3> mid{};for(int c=0;c<3;++c)mid[c]=0.5*(current[edge[0]][c]+current[edge[1]][c]);return mid;};
		for(std::size_t node=0;node<velocity_nodes;++node){
			const auto p=velocity_position(node);
			state[3*node]=p[1]*(height-p[1])*(1.0+p[2]);state[3*node+1]=0.1*p[0]*p[2];state[3*node+2]=0.05*p[1]*p[1];
		}
		for(std::size_t node=0;node<mesh.points.size();++node)state[3*velocity_nodes+node]=100.0*(length-current[node][0])+3.0*current[node][1];
		const auto replicated_traction=iga::EvaluateNativeTetFluidTractionOnMatchingStructure(current_mesh,topology,state,0.004,7);
		const auto distributed_traction=iga::EvaluateNativeTetFluidTractionDistributed(distributed,mesh,current,state,0.004,7);
		if(replicated_traction.size()!=distributed_traction.size())throw std::runtime_error("distributed traction count differs");
		std::vector<double> a,b;
		for(std::size_t i=0;i<replicated_traction.size();++i){
			if(replicated_traction[i].triangle_id!=distributed_traction[i].triangle_id)throw std::runtime_error("distributed traction order differs");
			for(int c=0;c<3;++c){a.push_back(replicated_traction[i].traction_on_structure_pa[c]);b.push_back(distributed_traction[i].traction_on_structure_pa[c]);}
		}
		const double traction_difference=RelativeDifference(a,b);
		if(!(traction_difference<1e-12))throw std::runtime_error("distributed traction differs: "+std::to_string(traction_difference));
		// Solid: the same box as a wall with one unused vertex, fixed at x=0, loaded on the top.
		auto solid=mesh;solid.points.push_back({{-1.0,-1.0,-1.0}});
		const iga::NativeTetMeshPartition solid_partition=iga::PartitionNativeTetMesh(solid,PETSC_COMM_WORLD);
		std::map<std::size_t,double> fixed;
		for(const auto& triangle:solid.boundary_triangles)if(triangle.boundary_label==1)
			for(const auto node:triangle.nodes)for(int c=0;c<3;++c)fixed.emplace(3*node+c,0.0);
		for(int c=0;c<3;++c)fixed.emplace(3*(solid.points.size()-1)+c,0.0);
		std::vector<double> force(3*solid.points.size(),0.0);
		for(const auto& triangle:solid.boundary_triangles)if(triangle.boundary_label==7)
			for(const auto node:triangle.nodes)force[3*node+1]-=1e-6;
		const iga::NativeTetSolidMaterial material{2.0e5,0.3,1000.0};
		iga::NativeTetSolidStaticOptions options;options.relative_tolerance=1e-9;options.absolute_tolerance_n=1e-14;
		options.linear_solve=iga::MakeNativeTetPetscLinearSolve("fsi_solid_");
		iga::NativeTetSolidStaticResult replicated_solid,distributed_solid;
		try{replicated_solid=iga::SolveNativeTetSolidStatic(solid,material,force,fixed,{},options);}
		catch(const std::exception& error){throw std::runtime_error(std::string("replicated solid: ")+error.what());}
		try{distributed_solid=iga::SolveNativeTetSolidStaticDistributed(solid_partition,solid,material,force,fixed,{},options,PETSC_COMM_WORLD);}
		catch(const std::exception& error){throw std::runtime_error(std::string("distributed solid: ")+error.what());}
		const double solid_difference=RelativeDifference(replicated_solid.displacement_m,distributed_solid.displacement_m);
		if(replicated_solid.iterations!=distributed_solid.iterations||!(solid_difference<1e-9)
			||std::abs(replicated_solid.strain_energy_j-distributed_solid.strain_energy_j)>1e-9*std::abs(replicated_solid.strain_energy_j))
			throw std::runtime_error("distributed solid differs: "+std::to_string(solid_difference)+" iterations "
				+std::to_string(replicated_solid.iterations)+"/"+std::to_string(distributed_solid.iterations));
		if(rank==0)std::cout<<"native FSI distributed parts passed ranks="<<ranks<<" motion="<<motion_difference
			<<" traction="<<traction_difference<<" solid="<<solid_difference<<" solid_iterations="<<distributed_solid.iterations<<'\n';
	}catch(const std::exception& error){
		int rank=0;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
		if(rank==0)std::cerr<<"native FSI distributed parts test: "<<error.what()<<'\n';
		status=1;
	}
	PetscFinalize();return status;
}

#ifndef IGA_CUDA_NATIVE_TET_SPECIES_KERNELS_CUH
#define IGA_CUDA_NATIVE_TET_SPECIES_KERNELS_CUH

#include <cuda_runtime.h>
#include <cstddef>

namespace iga::cuda {

struct NativeTetSpeciesCell
{
	int nodes[4];
	int velocity_nodes[10];
	double gradient[4][3];
	double position[4][3];
	double volume;
};

struct NativeTetSpeciesFace
{
	int nodes[3];
	int velocity_nodes[10];
	int local_nodes[3];
	int cell_index;
	int opposite;
	double area_vector[3];
	double area;
	double inlet_concentration;
	double transfer;
	double external_concentration;
};

__device__ inline void FaceBarycentric(int sample,double* barycentric)
{
	if(sample==0){barycentric[0]=1./3.;barycentric[1]=1./3.;
		barycentric[2]=1./3.;return;}
	if(sample<=3){
		for(int vertex=0;vertex<3;++vertex)
			barycentric[vertex]=vertex==sample-1?0.059715871789770
				:0.470142064105115;
		return;
	}
	for(int vertex=0;vertex<3;++vertex)
		barycentric[vertex]=vertex==sample-4?0.797426985353087
			:0.101286507323456;
}

__device__ inline double FaceWeight(int sample)
{
	return sample==0?0.225:(sample<=3?0.132394152788506:0.125939180544827);
}

__device__ inline void TaylorHoodFaceBasis(const double* lambda,double* shape)
{
	for(int node=0;node<4;++node)shape[node]=lambda[node]*(2.*lambda[node]-1.);
	shape[4]=4.*lambda[0]*lambda[1];
	shape[5]=4.*lambda[0]*lambda[2];
	shape[6]=4.*lambda[0]*lambda[3];
	shape[7]=4.*lambda[1]*lambda[2];
	shape[8]=4.*lambda[1]*lambda[3];
	shape[9]=4.*lambda[2]*lambda[3];
}

__global__ void AssembleNativeTetSpeciesCells(
	const NativeTetSpeciesCell* cells,int cell_count,int dofs,
	double dt,double diffusivity,double source,double decay,
	const double* velocity,const double* velocity_moments,
	const double* rt0_face_flow,
	const double* previous,double* matrix,double* rhs)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=cell_count*16)return;
	const auto& cell=cells[index/16];
	const int row=(index%16)/4,column=index%4;
	const double mass=cell.volume*(row==column?2.:1.)/20.;
	const double diffusion=diffusivity*cell.volume*(
		cell.gradient[row][0]*cell.gradient[column][0]
		+cell.gradient[row][1]*cell.gradient[column][1]
		+cell.gradient[row][2]*cell.gradient[column][2]);
	double advection=0.;
	if(rt0_face_flow){
		for(int vertex=0;vertex<4;++vertex){
			const double mass_factor=(vertex==column?2.:1.)/20.;
			for(int opposite=0;opposite<4;++opposite){
				const double flow=rt0_face_flow[4*(index/16)+opposite];
				for(int axis=0;axis<3;++axis)
					advection-=mass_factor*flow/3.*cell.gradient[row][axis]
						*(cell.position[vertex][axis]-cell.position[opposite][axis]);
			}
		}
	}else for(int local=0;local<10;++local){
		const double moment=velocity_moments[10*column+local];
		for(int axis=0;axis<3;++axis)
			advection-=cell.volume*cell.gradient[row][axis]*moment
				*velocity[3*cell.velocity_nodes[local]+axis];
	}
	atomicAdd(matrix+static_cast<std::size_t>(cell.nodes[column])*dofs
		+cell.nodes[row],mass/dt+decay*mass+diffusion+advection);
	atomicAdd(rhs+cell.nodes[row],mass*previous[cell.nodes[column]]/dt);
	if(column==0)atomicAdd(rhs+cell.nodes[row],source*cell.volume/4.);
}

__global__ void AssembleNativeTetSpeciesFaces(
	const NativeTetSpeciesFace* faces,int face_count,int dofs,
	const double* velocity,const double* rt0_face_flow,
	double* matrix,double* rhs)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=face_count*9)return;
	const auto& face=faces[index/9];
	const int row=(index%9)/3,column=index%3;
	const double mass=face.area*(row==column?2.:1.)/12.;
	double outflow=0.,inflow=0.;
	for(int sample=0;sample<7;++sample){
		double face_lambda[3],lambda[4]={0.,0.,0.,0.},shape[10];
		FaceBarycentric(sample,face_lambda);
		for(int node=0;node<3;++node)
			lambda[face.local_nodes[node]]=face_lambda[node];
		TaylorHoodFaceBasis(lambda,shape);
		double normal_flux=0.;
		if(rt0_face_flow)
			normal_flux=rt0_face_flow[4*face.cell_index+face.opposite];
		else for(int local=0;local<10;++local)
			for(int axis=0;axis<3;++axis)
				normal_flux+=shape[local]
					*velocity[3*face.velocity_nodes[local]+axis]
					*face.area_vector[axis];
		normal_flux*=FaceWeight(sample);
		outflow+=face_lambda[row]*face_lambda[column]*fmax(normal_flux,0.);
		if(column==0)
			inflow-=face_lambda[row]*fmin(normal_flux,0.)
				*face.inlet_concentration;
	}
	atomicAdd(matrix+static_cast<std::size_t>(face.nodes[column])*dofs
		+face.nodes[row],outflow+face.transfer*mass);
	if(column==0){
		atomicAdd(rhs+face.nodes[row],inflow);
		atomicAdd(rhs+face.nodes[row],
			face.transfer*face.external_concentration*face.area/3.);
	}
}

} // namespace iga::cuda

#endif

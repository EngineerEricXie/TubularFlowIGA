#ifndef IGA_CUDA_NATIVE_TET_DARCY_KERNELS_CUH
#define IGA_CUDA_NATIVE_TET_DARCY_KERNELS_CUH

#include <cuda_runtime.h>
#include <cstddef>

namespace iga::cuda {

struct NativeTetDarcyCell
{
	int nodes[4];
	double gradient[4][3];
	double volume;
	double mobility;
	double source;
};

__global__ void AssembleNativeTetDarcy(
	const NativeTetDarcyCell* cells,int cell_count,int dofs,
	double* matrix,double* right_hand_side)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=cell_count*20)return;
	const int cell=index/20,local=index%20;
	const auto& element=cells[cell];
	if(local<16){
		const int row=local/4,column=local%4;
		double stiffness=0.;
		for(int axis=0;axis<3;++axis)
			stiffness+=element.gradient[row][axis]*element.gradient[column][axis];
		stiffness*=element.mobility*element.volume;
		atomicAdd(matrix+static_cast<std::size_t>(element.nodes[column])*dofs
			+element.nodes[row],stiffness);
	}else{
		const int row=local-16;
		atomicAdd(right_hand_side+element.nodes[row],
			element.source*element.volume/4.);
	}
}

__global__ void SetNativeTetDarcyBoundaryValues(
	const int* fixed,const double* pressure,int dofs,double* right_hand_side)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index<dofs&&fixed[index])right_hand_side[index]=pressure[index];
}

__global__ void ApplyNativeTetDarcyBoundaryRows(
	const int* fixed,int dofs,double* matrix)
{
	const std::size_t index=static_cast<std::size_t>(blockIdx.x)*blockDim.x
		+threadIdx.x;
	if(index>=static_cast<std::size_t>(dofs)*dofs)return;
	const int row=static_cast<int>(index%dofs);
	const int column=static_cast<int>(index/dofs);
	if(fixed[row])matrix[index]=row==column?1.:0.;
}

__global__ void EvaluateNativeTetDarcyCellFlux(
	const NativeTetDarcyCell* cells,int cell_count,
	const double* pressure,double* cell_flux)
{
	const int cell=blockIdx.x*blockDim.x+threadIdx.x;
	if(cell>=cell_count)return;
	const auto& element=cells[cell];
	for(int axis=0;axis<3;++axis){
		double gradient=0.;
		for(int node=0;node<4;++node)
			gradient+=pressure[element.nodes[node]]*element.gradient[node][axis];
		cell_flux[3*cell+axis]=-element.mobility*gradient;
	}
}

} // namespace iga::cuda

#endif

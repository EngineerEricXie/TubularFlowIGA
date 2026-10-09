#ifndef IGA_CUDA_NATIVE_TET_DARCY_KERNELS_CUH
#define IGA_CUDA_NATIVE_TET_DARCY_KERNELS_CUH

#include "SparseKernels.cuh"

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
	const NativeTetDarcyCell* cells,int cell_count,DevicePatternView pattern,
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
		atomicAdd(matrix+FindBlock(pattern,element.nodes[row],element.nodes[column]),
			stiffness);
	}else{
		const int row=local-16;
		atomicAdd(right_hand_side+element.nodes[row],
			element.source*element.volume/4.);
	}
}

// Impose pressure on fixed nodes symmetrically, keeping the matrix SPD for CG:
// a fixed row becomes an identity row with the prescribed value, and a free
// row moves its fixed-column entries to the right-hand side. Each thread owns
// one row, so no atomics are needed.
__global__ void EliminateNativeTetDarcyPressure(
	DevicePatternView pattern,const int* fixed,const double* pressure,
	double* matrix,double* right_hand_side)
{
	const int row=blockIdx.x*blockDim.x+threadIdx.x;
	if(row>=pattern.nodes)return;
	if(fixed[row]){
		for(int entry=pattern.row_offsets[row];entry<pattern.row_offsets[row+1];++entry)
			matrix[entry]=pattern.columns[entry]==row?1.:0.;
		right_hand_side[row]=pressure[row];
		return;
	}
	for(int entry=pattern.row_offsets[row];entry<pattern.row_offsets[row+1];++entry){
		const int column=pattern.columns[entry];
		if(!fixed[column])continue;
		right_hand_side[row]-=matrix[entry]*pressure[column];
		matrix[entry]=0.;
	}
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

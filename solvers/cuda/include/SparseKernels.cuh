#ifndef IGA_CUDA_SPARSE_KERNELS_CUH
#define IGA_CUDA_SPARSE_KERNELS_CUH

#include "BlockCsr.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>

namespace iga::cuda {

__device__ inline int FindBlock(DevicePatternView pattern, int row, int column)
{
	int first = pattern.row_offsets[row];
	int last = pattern.row_offsets[row+1];
	while (first < last) {
		const int middle = first + (last-first)/2;
		if (pattern.columns[middle] < column) first = middle+1;
		else last = middle;
	}
	return first;
}

template <int Fields>
__global__ void BlockSpmvKernel(DevicePatternView pattern, const double* values,
	const double* x, double* y)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= pattern.nodes) return;
	double result[Fields]{};
	for (int block = pattern.row_offsets[row]; block < pattern.row_offsets[row+1]; ++block) {
		const int column = pattern.columns[block];
		const double* local = values+static_cast<std::size_t>(block)*Fields*Fields;
		for (int i = 0; i < Fields; ++i)
			for (int j = 0; j < Fields; ++j)
				result[i] += local[i*Fields+j]*x[column*Fields+j];
	}
	for (int i = 0; i < Fields; ++i) y[row*Fields+i] = result[i];
}

template <int Fields>
__global__ void BuildBlockInverseKernel(DevicePatternView pattern, const double* values,
	double* inverse, unsigned int* singular)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= pattern.nodes) return;
	const double* source = values+static_cast<std::size_t>(pattern.diagonal[node])*Fields*Fields;
	double augmented[Fields][2*Fields];
	for (int i = 0; i < Fields; ++i)
		for (int j = 0; j < 2*Fields; ++j)
			augmented[i][j] = j < Fields ? source[i*Fields+j] : (j-Fields == i ? 1.0 : 0.0);
	bool failed = false;
	for (int pivot = 0; pivot < Fields; ++pivot) {
		int best = pivot;
		for (int row = pivot+1; row < Fields; ++row)
			if (fabs(augmented[row][pivot]) > fabs(augmented[best][pivot])) best = row;
		if (!isfinite(augmented[best][pivot]) || fabs(augmented[best][pivot]) < 1e-14) {
			failed = true;
			break;
		}
		if (best != pivot)
			for (int column = 0; column < 2*Fields; ++column) {
				const double temporary = augmented[pivot][column];
				augmented[pivot][column] = augmented[best][column];
				augmented[best][column] = temporary;
			}
		const double scale = 1.0/augmented[pivot][pivot];
		for (int column = 0; column < 2*Fields; ++column) augmented[pivot][column] *= scale;
		for (int row = 0; row < Fields; ++row) {
			if (row == pivot) continue;
			const double factor = augmented[row][pivot];
			for (int column = 0; column < 2*Fields; ++column)
				augmented[row][column] -= factor*augmented[pivot][column];
		}
	}
	double* destination = inverse+static_cast<std::size_t>(node)*Fields*Fields;
	if (!failed) {
		for (int i = 0; i < Fields; ++i)
			for (int j = 0; j < Fields; ++j)
				destination[i*Fields+j] = augmented[i][Fields+j];
	} else {
		atomicAdd(singular, 1u);
		for (int i = 0; i < Fields; ++i)
			for (int j = 0; j < Fields; ++j)
				destination[i*Fields+j] = i == j && fabs(source[i*Fields+i]) > 1e-14
					? 1.0/source[i*Fields+i] : (i == j ? 1.0 : 0.0);
	}
}

template <int Fields>
__global__ void ApplyBlockInverseKernel(int nodes, const double* inverse,
	const double* x, double* y)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= nodes) return;
	for (int i = 0; i < Fields; ++i) {
		double result = 0.0;
		for (int j = 0; j < Fields; ++j)
			result += inverse[(static_cast<std::size_t>(node)*Fields+i)*Fields+j]
				*x[node*Fields+j];
		y[node*Fields+i] = result;
	}
}

// Scalar (Fields=1) Jacobi inverse without an absolute pivot threshold: P1
// diagonals scale with element volume and can be far below 1e-14 on fine
// meshes. Only zero or non-finite diagonals count as singular.
__global__ void BuildScalarJacobiKernel(DevicePatternView pattern, const double* values,
	double* inverse, unsigned int* singular)
{
	const int row = blockIdx.x*blockDim.x+threadIdx.x;
	if (row >= pattern.nodes) return;
	const double diagonal = values[pattern.diagonal[row]];
	if (diagonal == 0.0 || !isfinite(diagonal)) {
		atomicAdd(singular, 1u);
		inverse[row] = 1.0;
	} else inverse[row] = 1.0/diagonal;
}

} // namespace iga::cuda

#endif

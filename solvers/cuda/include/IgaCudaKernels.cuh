#ifndef IGA_CUDA_KERNELS_CUH
#define IGA_CUDA_KERNELS_CUH

#include "BlockCsr.hpp"
#include "CaseInput.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdint>

namespace iga::cuda {

constexpr int kMaximumBasis = 80;
constexpr int kQuadraturePoints = 64;
constexpr int kPairTile = 256;

struct ReferenceView {
	const double* basis = nullptr;
	const double* gradient = nullptr;
	const double* hessian = nullptr;
	const double* weight = nullptr;
	const double* surface_basis = nullptr;
	const double* surface_gradient = nullptr;
	const double* surface_weight = nullptr;
};

struct GeometryView {
	const double* determinant = nullptr;
	const double* inverse = nullptr;
	const double* inverse_second = nullptr;
};

struct ElementTilesView {
	int count = 0;
	const int* element = nullptr;
	const int* first_pair = nullptr;
};

__device__ inline double Determinant3(const double a[9])
{
	return a[0]*(a[4]*a[8]-a[5]*a[7])
		- a[1]*(a[3]*a[8]-a[5]*a[6])
		+ a[2]*(a[3]*a[7]-a[4]*a[6]);
}

__device__ inline void Inverse3(const double a[9], double inverse[9], double determinant)
{
	const double scale = 1.0 / determinant;
	inverse[0] = (a[4]*a[8]-a[5]*a[7])*scale;
	inverse[1] = (a[2]*a[7]-a[1]*a[8])*scale;
	inverse[2] = (a[1]*a[5]-a[2]*a[4])*scale;
	inverse[3] = (a[5]*a[6]-a[3]*a[8])*scale;
	inverse[4] = (a[0]*a[8]-a[2]*a[6])*scale;
	inverse[5] = (a[2]*a[3]-a[0]*a[5])*scale;
	inverse[6] = (a[3]*a[7]-a[4]*a[6])*scale;
	inverse[7] = (a[1]*a[6]-a[0]*a[7])*scale;
	inverse[8] = (a[0]*a[4]-a[1]*a[3])*scale;
}

__global__ void BuildGeometryKernel(DeviceMeshView mesh, ReferenceView reference,
	double* determinant, double* inverse, double* inverse_second, int with_hessian)
{
	const int index = blockIdx.x * blockDim.x + threadIdx.x;
	const int total = mesh.elements * kQuadraturePoints;
	if (index >= total) return;
	const int element = index / kQuadraturePoints;
	const int q = index - element * kQuadraturePoints;
	const double* points = mesh.bezier_points + static_cast<std::size_t>(element) * 64 * 3;
	double jacobian[9]{};
	for (int p = 0; p < 64; ++p)
		for (int physical = 0; physical < 3; ++physical)
			for (int parameter = 0; parameter < 3; ++parameter)
				jacobian[physical*3+parameter] += points[p*3+physical]
					* reference.gradient[(q*64+p)*3+parameter];
	const double raw_determinant = Determinant3(jacobian);
	determinant[index] = 0.125 * raw_determinant;
	double* inverse_out = inverse + static_cast<std::size_t>(index) * 9;
	if (!isfinite(raw_determinant) || fabs(raw_determinant) < 1e-14) {
		for (int i = 0; i < 9; ++i) inverse_out[i] = 0.0;
		if (with_hessian)
			for (int i = 0; i < 27; ++i)
				inverse_second[static_cast<std::size_t>(index)*27+i] = 0.0;
		return;
	}
	double inverse_local[9];
	Inverse3(jacobian, inverse_local, raw_determinant);
	for (int i = 0; i < 9; ++i) inverse_out[i] = inverse_local[i];
	if (!with_hessian) return;

	double geometry_second[27]{};
	for (int physical = 0; physical < 3; ++physical)
		for (int a = 0; a < 3; ++a)
			for (int b = 0; b < 3; ++b)
				for (int p = 0; p < 64; ++p)
					geometry_second[physical*9+a*3+b] += points[p*3+physical]
						* reference.hessian[(q*64+p)*9+a*3+b];
	double* second_out = inverse_second + static_cast<std::size_t>(index) * 27;
	for (int c = 0; c < 3; ++c)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j) {
				double value = 0.0;
				for (int physical = 0; physical < 3; ++physical)
					for (int a = 0; a < 3; ++a)
						for (int b = 0; b < 3; ++b)
							value -= geometry_second[physical*9+a*3+b]
								* inverse_local[a*3+i] * inverse_local[b*3+j]
								* inverse_local[c*3+physical];
				second_out[c*9+i*3+j] = value;
			}
}

__device__ inline void EvaluateBasisDevice(DeviceMeshView mesh, ReferenceView reference,
	GeometryView geometry, int element, int q, int local_basis, double& value,
	double gradient[3], double* hessian)
{
	const int basis_index = mesh.element_offsets[element] + local_basis;
	const int extraction_begin = mesh.extraction_offsets[basis_index];
	const int extraction_end = mesh.extraction_offsets[basis_index+1];
	double parametric_gradient[3]{};
	double parametric_hessian[9]{};
	value = 0.0;
	for (int entry = extraction_begin; entry < extraction_end; ++entry) {
		const int p = mesh.extraction_columns[entry];
		const double coefficient = mesh.extraction_values[entry];
		value += coefficient * reference.basis[q*64+p];
		for (int a = 0; a < 3; ++a)
			parametric_gradient[a] += coefficient * reference.gradient[(q*64+p)*3+a];
		if (hessian)
			for (int a = 0; a < 9; ++a)
				parametric_hessian[a] += coefficient * reference.hessian[(q*64+p)*9+a];
	}
	const int geometry_index = element*kQuadraturePoints + q;
	const double* inverse = geometry.inverse + static_cast<std::size_t>(geometry_index)*9;
	for (int physical = 0; physical < 3; ++physical) {
		gradient[physical] = 0.0;
		for (int parameter = 0; parameter < 3; ++parameter)
			gradient[physical] += parametric_gradient[parameter] * inverse[parameter*3+physical];
	}
	if (!hessian) return;
	const double* inverse_second = geometry.inverse_second
		+ static_cast<std::size_t>(geometry_index)*27;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			double result = 0.0;
			for (int a = 0; a < 3; ++a) {
				for (int b = 0; b < 3; ++b)
					result += parametric_hessian[a*3+b] * inverse[a*3+i] * inverse[b*3+j];
				result += parametric_gradient[a] * inverse_second[a*9+i*3+j];
			}
			hessian[i*3+j] = result;
		}
}

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

__global__ void AssembleTransportKernel(DeviceMeshView mesh, ReferenceView reference,
	GeometryView geometry, ElementTilesView tiles, DevicePatternView pattern,
	const double* velocity, iga::TransportParameters parameters,
	double* left, double* previous)
{
	const int tile = blockIdx.x;
	if (tile >= tiles.count) return;
	const int element = tiles.element[tile];
	const int begin = mesh.element_offsets[element];
	const int nen = mesh.element_offsets[element+1] - begin;
	const int pair = tiles.first_pair[tile] + threadIdx.x;
	const bool active = pair < nen*nen;
	const int a = active ? pair/nen : 0;
	const int b = active ? pair-a*nen : 0;
	int block = 0;
	if (active) {
		const int row = mesh.connectivity[begin+a];
		const int column = mesh.connectivity[begin+b];
		block = FindBlock(pattern, row, column);
	}
	__shared__ double basis[kMaximumBasis];
	__shared__ double gradient[kMaximumBasis*3];
	__shared__ double velocity_q[3];
	__shared__ double tau;
	double l00 = 0.0, l0p = 0.0, lp0 = 0.0, lpp = 0.0;
	double p00 = 0.0, ppp = 0.0;
	for (int q = 0; q < kQuadraturePoints; ++q) {
		if (threadIdx.x < nen) {
			double local_gradient[3];
			EvaluateBasisDevice(mesh, reference, geometry, element, q, threadIdx.x,
				basis[threadIdx.x], local_gradient, nullptr);
			for (int d = 0; d < 3; ++d) gradient[threadIdx.x*3+d] = local_gradient[d];
		}
		__syncthreads();
		if (threadIdx.x == 0) {
			for (int d = 0; d < 3; ++d) velocity_q[d] = 0.0;
			for (int local = 0; local < nen; ++local) {
				const int node = mesh.connectivity[begin+local];
				for (int d = 0; d < 3; ++d)
					velocity_q[d] += basis[local] * velocity[node*3+d];
			}
			double inverse_length = 0.0;
			for (int local = 0; local < nen; ++local)
				inverse_length += fabs(velocity_q[0]*gradient[local*3]
					+ velocity_q[1]*gradient[local*3+1]
					+ velocity_q[2]*gradient[local*3+2]);
			const double tau_space = inverse_length > 0.0 ? 1.0/inverse_length : 0.0;
			const double tau_time = parameters.dt/2.0;
			tau = tau_space > 0.0 ? 1.0/sqrt(1.0/(tau_space*tau_space)
				+ 1.0/(tau_time*tau_time)) : 0.0;
		}
		__syncthreads();
		if (active) {
			const double* ga = gradient+a*3;
			const double* gb = gradient+b*3;
			const double adv_a = velocity_q[0]*ga[0]+velocity_q[1]*ga[1]+velocity_q[2]*ga[2];
			const double adv_b = velocity_q[0]*gb[0]+velocity_q[1]*gb[1]+velocity_q[2]*gb[2];
			const double test_plus = basis[a] + tau*adv_a;
			const double gradient_dot = ga[0]*gb[0]+ga[1]*gb[1]+ga[2]*gb[2];
			const double measure = reference.weight[q]
				* geometry.determinant[element*kQuadraturePoints+q];
			const double mass = basis[a]*basis[b]*measure;
			l00 += ((1.0+parameters.dt*(parameters.kplus+parameters.kminus))*basis[a]*basis[b]
				+ parameters.dt*parameters.diffusion*gradient_dot)*measure;
			l0p += -parameters.detach_plus*parameters.dt*mass;
			lp0 += -parameters.kplus*parameters.dt*test_plus*basis[b]*measure;
			lpp += ((1.0+parameters.dt*parameters.detach_plus)*test_plus*basis[b]
				+ parameters.dt*test_plus*adv_b
				+ parameters.dt*parameters.vplus*parameters.artificial_diffusion
					*test_plus*gradient_dot)*measure;
			p00 += mass;
			ppp += test_plus*basis[b]*measure;
		}
		__syncthreads();
	}
	if (active) {
		double* l = left + static_cast<std::size_t>(block)*4;
		double* p = previous + static_cast<std::size_t>(block)*4;
		atomicAdd(l, l00); atomicAdd(l+1, l0p);
		atomicAdd(l+2, lp0); atomicAdd(l+3, lpp);
		atomicAdd(p, p00); atomicAdd(p+3, ppp);
	}
}

__device__ inline void StabilizationDevice(const double* inverse,
	const double state[4], double kinematic_viscosity, double dt,
	double& tau_m, double& tau_c)
{
	double metric[9]{};
	double direction[3]{};
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k)
				metric[i*3+j] += inverse[k*3+i]*inverse[k*3+j];
			direction[i] += inverse[j*3+i];
		}
	double metric_norm = 0.0, direction_norm = 0.0, velocity_metric = 0.0;
	for (int i = 0; i < 3; ++i) {
		direction_norm += direction[i]*direction[i];
		for (int j = 0; j < 3; ++j) {
			metric_norm += metric[i*3+j]*metric[i*3+j];
			velocity_metric += state[i]*metric[i*3+j]*state[j];
		}
	}
	const double temporal_scale = dt > 0.0 ? 4.0/(dt*dt) : 0.0;
	tau_m = 1.0/sqrt(temporal_scale+velocity_metric
		+(1.0/12.0)*kinematic_viscosity*kinematic_viscosity*metric_norm);
	tau_c = 1.0/(tau_m*direction_norm);
}

__global__ void AssembleNavierStokesKernel(DeviceMeshView mesh,
	ReferenceView reference, GeometryView geometry, ElementTilesView tiles,
	DevicePatternView pattern, const double* nodal_state, const double* previous_nodal_state,
	double density, double viscosity, double dt,
	double* jacobian, double* negative_residual)
{
	const int tile = blockIdx.x;
	if (tile >= tiles.count) return;
	const int element = tiles.element[tile];
	const int begin = mesh.element_offsets[element];
	const int nen = mesh.element_offsets[element+1] - begin;
	const int pair = tiles.first_pair[tile] + threadIdx.x;
	const bool active = pair < nen*nen;
	const int a = active ? pair/nen : 0;
	const int b = active ? pair-a*nen : 0;
	int block = 0;
	if (active) {
		const int row = mesh.connectivity[begin+a];
		const int column = mesh.connectivity[begin+b];
		block = FindBlock(pattern, row, column);
	}
	__shared__ double basis[kMaximumBasis];
	__shared__ double basis_gradient[kMaximumBasis*3];
	__shared__ double basis_hessian[kMaximumBasis*9];
	__shared__ double state[4];
	__shared__ double state_gradient[12];
	__shared__ double state_hessian[36];
	__shared__ double previous_velocity[3];
	__shared__ double time_derivative[3];
	__shared__ double fine_velocity[3];
	__shared__ double fine_pressure;
	__shared__ double tau_m;
	__shared__ double tau_c;
	double tangent_sum[16]{};
	double residual_sum[4]{};
	for (int q = 0; q < kQuadraturePoints; ++q) {
		if (threadIdx.x < nen) {
			double local_gradient[3], local_hessian[9];
			EvaluateBasisDevice(mesh, reference, geometry, element, q, threadIdx.x,
				basis[threadIdx.x], local_gradient, local_hessian);
			for (int d = 0; d < 3; ++d)
				basis_gradient[threadIdx.x*3+d] = local_gradient[d];
			for (int d = 0; d < 9; ++d)
				basis_hessian[threadIdx.x*9+d] = local_hessian[d];
		}
		__syncthreads();
		if (threadIdx.x == 0) {
			for (int field = 0; field < 4; ++field) state[field] = 0.0;
			for (int component = 0; component < 3; ++component) previous_velocity[component] = 0.0;
			for (int i = 0; i < 12; ++i) state_gradient[i] = 0.0;
			for (int i = 0; i < 36; ++i) state_hessian[i] = 0.0;
			for (int local = 0; local < nen; ++local) {
				const int node = mesh.connectivity[begin+local];
				for (int field = 0; field < 4; ++field) {
					const double coefficient = nodal_state[node*4+field];
					state[field] += coefficient*basis[local];
					for (int i = 0; i < 3; ++i) {
						state_gradient[field*3+i] += coefficient*basis_gradient[local*3+i];
						for (int j = 0; j < 3; ++j)
							state_hessian[field*9+i*3+j] += coefficient
								*basis_hessian[local*9+i*3+j];
					}
				}
				if (dt > 0.0)
					for (int component = 0; component < 3; ++component)
						previous_velocity[component] += previous_nodal_state[node*4+component]
							*basis[local];
			}
			for (int component = 0; component < 3; ++component)
				time_derivative[component] = dt > 0.0
					? (state[component]-previous_velocity[component])/dt : 0.0;
			const double* inverse = geometry.inverse
				+ static_cast<std::size_t>(element*kQuadraturePoints+q)*9;
			const double kinematic_viscosity = viscosity/density;
			StabilizationDevice(inverse, state, kinematic_viscosity, dt, tau_m, tau_c);
			for (int component = 0; component < 3; ++component) {
				const double convection = state[0]*state_gradient[component*3]
					+ state[1]*state_gradient[component*3+1]
					+ state[2]*state_gradient[component*3+2];
				const double pressure_gradient = state_gradient[9+component];
				const double laplacian = state_hessian[component*9]
					+ state_hessian[component*9+4]+state_hessian[component*9+8];
				fine_velocity[component] = -tau_m*(time_derivative[component]+convection
					+pressure_gradient/density-kinematic_viscosity*laplacian);
			}
			fine_pressure = -density*tau_c
				*(state_gradient[0]+state_gradient[4]+state_gradient[8]);
		}
		__syncthreads();
		if (active) {
			const double na = basis[a], nb = basis[b];
			const double* ga = basis_gradient+a*3;
			const double* gb = basis_gradient+b*3;
			const double measure = reference.weight[q]
				* geometry.determinant[element*kQuadraturePoints+q];
			if (b == 0) {
				double residual[4]{};
				for (int component = 0; component < 3; ++component) {
					residual[component] = density*na*time_derivative[component]
						-ga[component]*state[3]-ga[component]*fine_pressure;
					for (int direction = 0; direction < 3; ++direction) {
						residual[component] += viscosity*ga[direction]
							*(state_gradient[component*3+direction]
								+state_gradient[direction*3+component]);
						residual[component] += density*na*(state[direction]+fine_velocity[direction])
							*state_gradient[component*3+direction];
						residual[component] -= density*ga[direction]*fine_velocity[component]
							*(state[direction]+fine_velocity[direction]);
					}
				}
				residual[3] = na*(state_gradient[0]+state_gradient[4]+state_gradient[8])
					-ga[0]*fine_velocity[0]-ga[1]*fine_velocity[1]-ga[2]*fine_velocity[2];
				for (int field = 0; field < 4; ++field)
					residual_sum[field] -= residual[field]*measure;
			}
			const double convection_b = state[0]*gb[0]+state[1]*gb[1]+state[2]*gb[2];
			const double streamline_a = state[0]*ga[0]+state[1]*ga[1]+state[2]*ga[2];
			const double gradient_dot = ga[0]*gb[0]+ga[1]*gb[1]+ga[2]*gb[2];
			double tangent[16]{};
			const double mass_b = dt > 0.0 ? nb/dt : 0.0;
			const double diagonal = density*na*(mass_b+convection_b)+viscosity*gradient_dot
				+density*tau_m*streamline_a*(mass_b+convection_b);
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					tangent[i*4+j] = viscosity*ga[j]*gb[i]+density*tau_c*ga[i]*gb[j];
			for (int i = 0; i < 3; ++i) tangent[i*4+i] += diagonal;
			for (int i = 0; i < 3; ++i) {
				tangent[i*4+3] = -ga[i]*nb+tau_m*streamline_a*gb[i];
				tangent[12+i] = na*gb[i]+tau_m*ga[i]*(mass_b+convection_b);
			}
			tangent[15] = (tau_m/density)*gradient_dot;
			for (int i = 0; i < 16; ++i) tangent_sum[i] += tangent[i]*measure;
		}
		__syncthreads();
	}
	if (active) {
		double* values = jacobian+static_cast<std::size_t>(block)*16;
		for (int i = 0; i < 16; ++i) atomicAdd(values+i, tangent_sum[i]);
		if (b == 0) {
			const int row = mesh.connectivity[begin+a];
			for (int field = 0; field < 4; ++field)
				atomicAdd(negative_residual+row*4+field, residual_sum[field]);
		}
	}
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

__global__ void ApplyTransportBoundaryKernel(DevicePatternView pattern,
	const int* constrained, double* left, double* previous)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= pattern.nodes || !constrained[node]) return;
	for (int block = pattern.row_offsets[node]; block < pattern.row_offsets[node+1]; ++block)
		for (int i = 0; i < 4; ++i) {
			left[static_cast<std::size_t>(block)*4+i] = 0.0;
			previous[static_cast<std::size_t>(block)*4+i] = 0.0;
		}
	double* diagonal = left+static_cast<std::size_t>(pattern.diagonal[node])*4;
	diagonal[0] = 1.0;
	diagonal[3] = 1.0;
}

__global__ void SetTransportBoundaryVectorKernel(int nodes, const int* constrained,
	const double* n0, const double* nplus, double* vector)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node < nodes && constrained[node]) {
		vector[node*2] = n0[node];
		vector[node*2+1] = nplus[node];
	}
}

__global__ void ApplyNavierStokesBoundaryKernel(DevicePatternView pattern,
	const int* velocity_constrained, const int* pressure_constrained, double* jacobian)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= pattern.nodes) return;
	for (int field = 0; field < 4; ++field) {
		const bool constrained = field < 3 ? velocity_constrained[node] : pressure_constrained[node];
		if (!constrained) continue;
		for (int block = pattern.row_offsets[node]; block < pattern.row_offsets[node+1]; ++block)
			for (int column_field = 0; column_field < 4; ++column_field)
				jacobian[static_cast<std::size_t>(block)*16+field*4+column_field] = 0.0;
		jacobian[static_cast<std::size_t>(pattern.diagonal[node])*16+field*4+field] = 1.0;
	}
}

__global__ void SetNavierStokesBoundaryRhsKernel(int nodes,
	const int* velocity_constrained, const int* pressure_constrained,
	const double* velocity, const double* pressure, const double* state, double* rhs)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= nodes) return;
	for (int field = 0; field < 4; ++field) {
		const bool constrained = field < 3 ? velocity_constrained[node] : pressure_constrained[node];
		if (!constrained) continue;
		const double target = field < 3 ? velocity[node*3+field] : pressure[node];
		rhs[node*4+field] = target-state[node*4+field];
	}
}

} // namespace iga::cuda

#endif

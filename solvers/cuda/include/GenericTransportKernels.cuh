#ifndef IGA_CUDA_GENERIC_TRANSPORT_KERNELS_CUH
#define IGA_CUDA_GENERIC_TRANSPORT_KERNELS_CUH

#include "IgaCudaKernels.cuh"

namespace iga::cuda {

struct DeviceWeakTerm {
	int kind;
	int equation;
	int trial;
	double coefficient;
};

struct DeviceSurfaceView {
	int count = 0;
	const int* element = nullptr;
	const int* face = nullptr;
	const int* first_pair = nullptr;
	const double* flux = nullptr;
	const double* robin = nullptr;
	const double* robin_source = nullptr;
};

template <int Fields>
__global__ void AssembleGenericTransportKernel(DeviceMeshView mesh, ReferenceView reference,
	GeometryView geometry, ElementTilesView tiles, DevicePatternView pattern,
	const double* velocity, const DeviceWeakTerm* terms, int term_count,
	const int* supg, double dt, double* left, double* previous, double* source)
{
	const int tile = blockIdx.x;
	if (tile >= tiles.count) return;
	const int element = tiles.element[tile];
	const int begin = mesh.element_offsets[element];
	const int nen = mesh.element_offsets[element+1]-begin;
	const int pair = tiles.first_pair[tile]+threadIdx.x;
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
	double local_left[Fields*Fields]{};
	double local_previous[Fields*Fields]{};
	double local_source[Fields]{};
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
				for (int d = 0; d < 3; ++d) velocity_q[d] += basis[local]*velocity[node*3+d];
			}
			double inverse_length = 0.0;
			for (int local = 0; local < nen; ++local)
				inverse_length += fabs(velocity_q[0]*gradient[local*3]
					+velocity_q[1]*gradient[local*3+1]+velocity_q[2]*gradient[local*3+2]);
			const double tau_space = inverse_length > 0.0 ? 1.0/inverse_length : 0.0;
			const double tau_time = dt/2.0;
			tau = tau_space > 0.0 ? 1.0/sqrt(1.0/(tau_space*tau_space)
				+1.0/(tau_time*tau_time)) : 0.0;
		}
		__syncthreads();
		if (active) {
			const double* ga = gradient+a*3;
			const double* gb = gradient+b*3;
			const double adv_a = velocity_q[0]*ga[0]+velocity_q[1]*ga[1]+velocity_q[2]*ga[2];
			const double adv_b = velocity_q[0]*gb[0]+velocity_q[1]*gb[1]+velocity_q[2]*gb[2];
			const double gradient_dot = ga[0]*gb[0]+ga[1]*gb[1]+ga[2]*gb[2];
			const double measure = reference.weight[q]*geometry.determinant[element*kQuadraturePoints+q];
			for (int t = 0; t < term_count; ++t) {
				const DeviceWeakTerm term = terms[t];
				const double test = basis[a]+(supg[term.equation] ? tau*adv_a : 0.0);
				if (term.kind == static_cast<int>(TermKind::VolumeSource)) {
					if (b == 0) local_source[term.equation] += dt*term.coefficient*test*measure;
					continue;
				}
				const int entry = term.equation*Fields+term.trial;
				if (term.kind == static_cast<int>(TermKind::TimeDerivative)) {
					const double mass = term.coefficient*test*basis[b]*measure;
					local_left[entry] += mass;
					local_previous[entry] += mass;
				} else if (term.kind == static_cast<int>(TermKind::Diffusion)) {
					const double weak_gradient = supg[term.equation] ? test*gradient_dot : gradient_dot;
					local_left[entry] += dt*term.coefficient*weak_gradient*measure;
				} else if (term.kind == static_cast<int>(TermKind::Advection)) {
					local_left[entry] += dt*term.coefficient*test*adv_b*measure;
				} else if (term.kind == static_cast<int>(TermKind::LinearCoupling)) {
					local_left[entry] += dt*term.coefficient*test*basis[b]*measure;
				}
			}
		}
		__syncthreads();
	}
	if (active) {
		double* l = left+static_cast<std::size_t>(block)*Fields*Fields;
		double* p = previous+static_cast<std::size_t>(block)*Fields*Fields;
		for (int i = 0; i < Fields*Fields; ++i) {
			atomicAdd(l+i, local_left[i]);
			atomicAdd(p+i, local_previous[i]);
		}
		if (b == 0) {
			const int row = mesh.connectivity[begin+a];
			for (int field = 0; field < Fields; ++field)
				atomicAdd(source+row*Fields+field, local_source[field]);
		}
	}
}

template <int Fields>
__global__ void AssembleGenericSurfaceKernel(DeviceMeshView mesh, ReferenceView reference,
	DeviceSurfaceView surfaces, DevicePatternView pattern, double dt,
	double* left, double* source)
{
	const int tile = blockIdx.x;
	if (tile >= surfaces.count) return;
	const int element = surfaces.element[tile];
	const int face = surfaces.face[tile];
	const int begin = mesh.element_offsets[element];
	const int nen = mesh.element_offsets[element+1]-begin;
	const int pair = surfaces.first_pair[tile]+threadIdx.x;
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
	__shared__ double measure;
	double local_left[Fields]{};
	double local_source[Fields]{};
	for (int qface = 0; qface < 16; ++qface) {
		const int q = face*16+qface;
		if (threadIdx.x < nen) {
			const int basis_index = begin+threadIdx.x;
			double value = 0.0;
			for (int entry = mesh.extraction_offsets[basis_index];
				entry < mesh.extraction_offsets[basis_index+1]; ++entry) {
				const int p = mesh.extraction_columns[entry];
				value += mesh.extraction_values[entry]*reference.surface_basis[q*64+p];
			}
			basis[threadIdx.x] = value;
		}
		if (threadIdx.x == 0) {
			constexpr int varying_axes[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 2}, {1, 2}, {0, 1}};
			double tangent[2][3]{};
			const double* points = mesh.bezier_points+static_cast<std::size_t>(element)*64*3;
			for (int p = 0; p < 64; ++p)
				for (int d = 0; d < 3; ++d)
					for (int direction = 0; direction < 2; ++direction)
						tangent[direction][d] += points[p*3+d]
							* reference.surface_gradient[(q*64+p)*3+varying_axes[face][direction]];
			const double cross[3] = {
				tangent[0][1]*tangent[1][2]-tangent[0][2]*tangent[1][1],
				tangent[0][2]*tangent[1][0]-tangent[0][0]*tangent[1][2],
				tangent[0][0]*tangent[1][1]-tangent[0][1]*tangent[1][0]};
			measure = reference.surface_weight[q]
				* sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
		}
		__syncthreads();
		if (active) {
			const std::size_t condition = static_cast<std::size_t>(tile)*Fields;
			for (int field = 0; field < Fields; ++field) {
				local_left[field] += dt*surfaces.robin[condition+field]*basis[a]*basis[b]*measure;
				if (b == 0)
					local_source[field] += dt*(surfaces.flux[condition+field]
						+surfaces.robin_source[condition+field])*basis[a]*measure;
			}
		}
		__syncthreads();
	}
	if (!active) return;
	double* matrix = left+static_cast<std::size_t>(block)*Fields*Fields;
	for (int field = 0; field < Fields; ++field)
		atomicAdd(matrix+field*Fields+field, local_left[field]);
	if (b == 0) {
		const int row = mesh.connectivity[begin+a];
		for (int field = 0; field < Fields; ++field)
			atomicAdd(source+row*Fields+field, local_source[field]);
	}
}

template <int Fields>
__global__ void ApplyGenericBoundaryKernel(DevicePatternView pattern,
	const int* constrained, double* left, double* previous)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= pattern.nodes) return;
	for (int field = 0; field < Fields; ++field) {
		if (!constrained[node*Fields+field]) continue;
		for (int block = pattern.row_offsets[node]; block < pattern.row_offsets[node+1]; ++block)
			for (int column = 0; column < Fields; ++column) {
				left[static_cast<std::size_t>(block)*Fields*Fields+field*Fields+column] = 0.0;
				previous[static_cast<std::size_t>(block)*Fields*Fields+field*Fields+column] = 0.0;
			}
		left[static_cast<std::size_t>(pattern.diagonal[node])*Fields*Fields+field*Fields+field] = 1.0;
	}
}

template <int Fields>
__global__ void SetGenericBoundaryVectorKernel(int nodes, const int* constrained,
	const double* boundary_value, double* vector)
{
	const int node = blockIdx.x*blockDim.x+threadIdx.x;
	if (node >= nodes) return;
	for (int field = 0; field < Fields; ++field)
		if (constrained[node*Fields+field]) vector[node*Fields+field] = boundary_value[node*Fields+field];
}

__global__ void AddVectorKernel(int size, const double* addend, double* destination)
{
	const int index = blockIdx.x*blockDim.x+threadIdx.x;
	if (index < size) destination[index] += addend[index];
}

} // namespace iga::cuda

#endif

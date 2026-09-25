#ifndef IGA_CUDA_NATIVE_TET_FLOW_KERNELS_CUH
#define IGA_CUDA_NATIVE_TET_FLOW_KERNELS_CUH

#include <cuda_runtime.h>
#include <cstddef>

namespace iga::cuda {

constexpr int kTetFlowDofs=34;
constexpr int kTetFlowQuadrature=64;
constexpr int kTetFlowEntries=kTetFlowDofs+kTetFlowDofs*kTetFlowDofs;

struct NativeTetFlowQuadrature
{
	double weight;
	double velocity[10];
	double gradient[10][3];
	double pressure[4];
	double position[3];
};

struct NativeTetFlowField
{
	double velocity[3];
	double previous_velocity[3];
	double grid_velocity[3];
	double gradient[3][3];
	double pressure;
};

__global__ void EvaluateNativeTetFlowFields(
	const NativeTetFlowQuadrature* quadrature,
	const int* rows,const double* state,const double* previous_state,
	const double* mesh_velocity,int cells,NativeTetFlowField* fields)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=cells*kTetFlowQuadrature)return;
	const int cell=index/kTetFlowQuadrature;
	const auto& basis=quadrature[index];
	auto& field=fields[index];
	const int* local=rows+cell*kTetFlowDofs;
	for(int component=0;component<3;++component){
		field.velocity[component]=0.;
		field.previous_velocity[component]=0.;
		field.grid_velocity[component]=0.;
		for(int axis=0;axis<3;++axis)field.gradient[component][axis]=0.;
	}
	field.pressure=0.;
	for(int node=0;node<10;++node)
		for(int component=0;component<3;++component){
			const double value=state[local[3*node+component]];
			field.velocity[component]+=basis.velocity[node]*value;
			field.previous_velocity[component]+=
				basis.velocity[node]*previous_state[local[3*node+component]];
			for(int axis=0;axis<3;++axis)
				field.gradient[component][axis]+=basis.gradient[node][axis]*value;
		}
	for(int node=0;node<4;++node){
		field.pressure+=basis.pressure[node]*state[local[30+node]];
		for(int component=0;component<3;++component)
			field.grid_velocity[component]+=
				basis.pressure[node]*mesh_velocity[cell*12+3*node+component];
	}
}

__global__ void AssembleNativeTetFlowElements(
	const NativeTetFlowQuadrature* quadrature,
	const NativeTetFlowField* fields,const double* body_acceleration,
	double density,double viscosity,double dt_s,int cells,
	double* element_residual,double* element_jacobian)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=cells*kTetFlowEntries)return;
	const int cell=index/kTetFlowEntries;
	const int slot=index-cell*kTetFlowEntries;
	const bool residual=slot<kTetFlowDofs;
	const int entry=residual?slot:slot-kTetFlowDofs;
	const int row=residual?entry:entry/kTetFlowDofs;
	const int column=residual?-1:entry%kTetFlowDofs;
	double total=0.;
	for(int q=0;q<kTetFlowQuadrature;++q){
		const auto& basis=quadrature[cell*kTetFlowQuadrature+q];
		const auto& field=fields[cell*kTetFlowQuadrature+q];
		const double weight=basis.weight;
		if(row<30){
			const int node=row/3,component=row%3;
			if(residual){
				double viscous=0.,convection=0.;
				for(int axis=0;axis<3;++axis){
					viscous+=viscosity*(field.gradient[component][axis]
						+field.gradient[axis][component])*basis.gradient[node][axis];
					convection+=(field.velocity[axis]-field.grid_velocity[axis])
						*field.gradient[component][axis];
				}
				double acceleration=body_acceleration[4*component];
				for(int axis=0;axis<3;++axis)
					acceleration+=body_acceleration[4*component+axis+1]
						*basis.position[axis];
				total+=weight*(viscous+density*basis.velocity[node]*convection
					-field.pressure*basis.gradient[node][component]
					-density*basis.velocity[node]*acceleration);
				if(dt_s>0.)total+=weight*density*basis.velocity[node]
					*(field.velocity[component]-field.previous_velocity[component])/dt_s;
			}else if(column<30){
				const int other=column/3,unknown=column%3;
				double tangent=viscosity*(
					(component==unknown?
						basis.gradient[other][0]*basis.gradient[node][0]
						+basis.gradient[other][1]*basis.gradient[node][1]
						+basis.gradient[other][2]*basis.gradient[node][2]:0.)
					+basis.gradient[other][component]*basis.gradient[node][unknown]);
				tangent+=density*basis.velocity[node]*(
					basis.velocity[other]*field.gradient[component][unknown]
					+(component==unknown?
						(field.velocity[0]-field.grid_velocity[0])*basis.gradient[other][0]
						+(field.velocity[1]-field.grid_velocity[1])*basis.gradient[other][1]
						+(field.velocity[2]-field.grid_velocity[2])*basis.gradient[other][2]
						:0.));
				if(dt_s>0.&&component==unknown)
					tangent+=density*basis.velocity[node]*basis.velocity[other]/dt_s;
				total+=weight*tangent;
			}else total-=weight*basis.pressure[column-30]
				*basis.gradient[node][component];
		}else{
			const int node=row-30;
			if(residual){
				const double divergence=field.gradient[0][0]
					+field.gradient[1][1]+field.gradient[2][2];
				total+=weight*basis.pressure[node]*divergence;
			}else if(column<30)
				total+=weight*basis.pressure[node]
					*basis.gradient[column/3][column%3];
		}
	}
	if(residual)element_residual[cell*kTetFlowDofs+row]=total;
	else element_jacobian[cell*kTetFlowDofs*kTetFlowDofs+entry]=total;
}

__global__ void ScatterNativeTetFlowElements(
	const int* rows,const double* element_residual,
	const double* element_jacobian,int cells,int global_dofs,
	double* residual,double* jacobian)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=cells*kTetFlowEntries)return;
	const int cell=index/kTetFlowEntries;
	const int slot=index-cell*kTetFlowEntries;
	const int* local=rows+cell*kTetFlowDofs;
	if(slot<kTetFlowDofs){
		atomicAdd(residual+local[slot],element_residual[cell*kTetFlowDofs+slot]);
	}else{
		const int entry=slot-kTetFlowDofs;
		const int row=entry/kTetFlowDofs,column=entry%kTetFlowDofs;
		atomicAdd(jacobian+static_cast<std::size_t>(local[column])*global_dofs
			+local[row],element_jacobian[cell*kTetFlowDofs*kTetFlowDofs+entry]);
	}
}

__global__ void ApplyNativeTetFlowBoundaryRows(
	const int* constrained,int dofs,double* jacobian)
{
	const std::size_t index=static_cast<std::size_t>(blockIdx.x)*blockDim.x
		+threadIdx.x;
	if(index>=static_cast<std::size_t>(dofs)*dofs)return;
	const int row=static_cast<int>(index%dofs);
	const int column=static_cast<int>(index/dofs);
	if(constrained[row])jacobian[index]=row==column?1.:0.;
}

__global__ void NativeTetFlowNewtonRightHandSide(
	const double* residual,const double* boundary_load,
	const double* state,const int* constrained,const double* prescribed,
	int dofs,double* right_hand_side)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index>=dofs)return;
	right_hand_side[index]=constrained[index]
		?prescribed[index]-state[index]
		:-residual[index]-boundary_load[index];
}

__global__ void UpdateNativeTetFlowState(double* state,const double* increment,int dofs)
{
	const int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index<dofs)state[index]+=increment[index];
}

} // namespace iga::cuda

#endif

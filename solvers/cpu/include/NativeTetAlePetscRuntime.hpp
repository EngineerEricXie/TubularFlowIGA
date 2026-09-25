#ifndef IGA_NATIVE_TET_ALE_PETSC_RUNTIME_HPP
#define IGA_NATIVE_TET_ALE_PETSC_RUNTIME_HPP

#include "NativeTetAleBoundaryControl.hpp"
#include "NativeTetAleBackflow.hpp"
#include "NativeTetAlePetscAssembly.hpp"
#include "NativeTetAleTransient.hpp"

#include <petscksp.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct NativeTetAlePetscSolveResult
{
	std::vector<double> replicated_state;
	std::size_t newton_iterations = 0;
	double final_residual_l2 = 0.0;
	PetscInt total_linear_iterations = 0;
	PetscInt final_linear_converged_reason = 0;
	std::string linear_solver_type;
	std::string preconditioner_type;
	bool system_scaling_enabled = false;
	bool user_schur_enabled = false;
	std::vector<double> flow_controller_multipliers_pa;
};

inline NativeTetAlePetscSolveResult SolveNativeTetAlePetsc(
	const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& vertex_mesh_velocity_m_s,
	const std::vector<double>& committed_state,const std::vector<double>& initial_trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double dt_s,bool steady,double absolute_tolerance=1e-10,
	std::size_t maximum_iterations=10,
	const NativeTetAleBoundaryConditions& boundary_conditions={})
{
	int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	const auto topology=BuildNativeTaylorHoodTopology(current_mesh);
	const auto boundary_operators=BuildNativeTetAleBoundaryFluxOperators(
		current_mesh,topology);
	ValidateNativeTetAleBoundaryConditions(boundary_conditions,boundary_operators);
	const auto backflow_owners=boundary_conditions.backflow_stabilization_beta.empty()
		?std::vector<std::size_t>{}:NativeTetAleBoundaryFaceOwners(current_mesh);
	const std::size_t velocity_nodes=current_mesh.points.size()+topology.edges.size();
	const std::size_t physical_dofs_size=3*velocity_nodes+current_mesh.points.size();
	const std::size_t controller_dofs=boundary_conditions.flow_rate_controls.size();
	const std::size_t dofs_size=physical_dofs_size+controller_dofs;
	const bool pressure_gauge_enabled=
		pressure_gauge_node!=std::numeric_limits<std::uint32_t>::max();
	if(pressure_gauge_enabled&&!boundary_conditions.prescribed_pressure_pa.empty())
		throw std::invalid_argument("native ALE pressure outlet and point gauge conflict");
	if(!pressure_gauge_enabled&&boundary_conditions.prescribed_pressure_pa.empty())
		throw std::invalid_argument("native ALE fluid pressure needs a gauge or outlet");
	if(dofs_size>static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())
		||vertex_mesh_velocity_m_s.size()!=current_mesh.points.size()
		||committed_state.size()!=physical_dofs_size
		||initial_trial_state.size()!=physical_dofs_size
		||(pressure_gauge_enabled&&pressure_gauge_node>=current_mesh.points.size())
		||(!steady&&!(dt_s>0.0))
		||!(absolute_tolerance>0.0)||maximum_iterations==0)
		throw std::invalid_argument("native ALE PETSc runtime input is invalid");
	const PetscInt dofs=static_cast<PetscInt>(dofs_size);
	Mat matrix=nullptr,schur_preconditioner=nullptr;
	Vec state=nullptr,residual=nullptr,right=nullptr,update=nullptr,base_state=nullptr;
	Vec row_scale=nullptr,column_scale=nullptr;
	KSP solver=nullptr;VecScatter scatter=nullptr;Vec replicated=nullptr;
	NativeTetAlePetscCheck(MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
		dofs,dofs,80,nullptr,80,nullptr,&matrix),"MatCreateAIJ runtime");
	NativeTetAlePetscCheck(MatSetOption(matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),
		"MatSetOption runtime");
	NativeTetAlePetscCheck(MatCreateVecs(matrix,&state,&residual),"MatCreateVecs runtime");
	NativeTetAlePetscCheck(VecDuplicate(residual,&right),"VecDuplicate right");
	NativeTetAlePetscCheck(VecDuplicate(state,&update),"VecDuplicate update");
	NativeTetAlePetscCheck(VecDuplicate(state,&base_state),"VecDuplicate Newton base state");
	NativeTetAlePetscCheck(VecDuplicate(state,&row_scale),"VecDuplicate row scale runtime");
	NativeTetAlePetscCheck(VecDuplicate(state,&column_scale),"VecDuplicate column scale runtime");
	PetscInt begin=0,end=0;VecGetOwnershipRange(state,&begin,&end);
	PetscScalar* local_state=nullptr;VecGetArray(state,&local_state);
	for(PetscInt row=begin;row<end;++row)
		local_state[row-begin]=static_cast<std::size_t>(row)<physical_dofs_size
			?initial_trial_state[static_cast<std::size_t>(row)]:0.0;
	VecRestoreArray(state,&local_state);VecAssemblyBegin(state);VecAssemblyEnd(state);
	PetscBool scale_system=PETSC_FALSE;
	PetscOptionsGetBool(nullptr,nullptr,"-native_ale_scale",&scale_system,nullptr);
	if(scale_system&&controller_dofs)
		throw std::invalid_argument("native ALE flow controllers do not yet support system scaling");
	if(scale_system) {
		double characteristic_velocity=0.0;
		for(const auto& prescribed:prescribed_velocity) {
			const auto& value=prescribed.second;
			characteristic_velocity=std::max(characteristic_velocity,
				std::sqrt(value[0]*value[0]+value[1]*value[1]+value[2]*value[2]));
		}
		if(!(characteristic_velocity>0.0)) characteristic_velocity=1.0;
		std::array<double,3> minimum=current_mesh.points.front();
		std::array<double,3> maximum=minimum;
		for(const auto& point:current_mesh.points)
			for(int component=0;component<3;++component) {
				minimum[component]=std::min(minimum[component],point[component]);
				maximum[component]=std::max(maximum[component],point[component]);
			}
		double characteristic_length=std::numeric_limits<double>::max();
		for(int component=0;component<3;++component) {
			const double extent=maximum[component]-minimum[component];
			if(extent>0.0) characteristic_length=std::min(characteristic_length,extent);
		}
		if(!std::isfinite(characteristic_length)||!(characteristic_length>0.0))
			throw std::runtime_error("native ALE scaling cannot determine a length scale");
		const double pressure_scale=(steady?0.0:parameters.density*characteristic_velocity
			*characteristic_length/dt_s)+parameters.dynamic_viscosity
			*characteristic_velocity/characteristic_length+parameters.density
			*characteristic_velocity*characteristic_velocity;
		PetscInt scale_begin=0,scale_end=0;
		VecGetOwnershipRange(row_scale,&scale_begin,&scale_end);
		PetscScalar* row_values=nullptr;PetscScalar* column_values=nullptr;
		VecGetArray(row_scale,&row_values);VecGetArray(column_scale,&column_values);
		const PetscInt velocity_dofs=static_cast<PetscInt>(3*velocity_nodes);
		for(PetscInt row=scale_begin;row<scale_end;++row) {
			const bool velocity=row<velocity_dofs;
			row_values[row-scale_begin]=velocity
				?1.0/(pressure_scale*characteristic_length)
				:1.0/(characteristic_velocity*characteristic_length);
			column_values[row-scale_begin]=velocity?characteristic_velocity:pressure_scale;
		}
		VecRestoreArray(row_scale,&row_values);VecRestoreArray(column_scale,&column_values);
	}
	std::map<PetscInt,double> constraints;
	for(const auto& prescribed:prescribed_velocity) {
		if(prescribed.first>=velocity_nodes)
			throw std::invalid_argument("native ALE PETSc boundary node is out of range");
		for(int component=0;component<3;++component)
			constraints.emplace(static_cast<PetscInt>(3*prescribed.first+component),
				prescribed.second[component]);
	}
	if(pressure_gauge_enabled)
		constraints.emplace(static_cast<PetscInt>(3*velocity_nodes+pressure_gauge_node),0.0);
	std::vector<PetscInt> local_constraint_rows;
	std::vector<PetscScalar> local_constraint_values;
	for(const auto& constraint:constraints)
		if(constraint.first>=begin&&constraint.first<end) {
			local_constraint_rows.push_back(constraint.first);
			local_constraint_values.push_back(constraint.second);
		}
	NativeTetAlePetscCheck(VecSetValues(state,static_cast<PetscInt>(local_constraint_rows.size()),
		local_constraint_rows.data(),local_constraint_values.data(),INSERT_VALUES),
		"VecSetValues initial constraints");
	VecAssemblyBegin(state);VecAssemblyEnd(state);
	NativeTetAlePetscCheck(VecScatterCreateToAll(state,&scatter,&replicated),
		"VecScatterCreateToAll runtime");
	NativeTetAlePetscCheck(KSPCreate(PETSC_COMM_WORLD,&solver),"KSPCreate runtime");
	KSPSetType(solver,KSPPREONLY);PC pc=nullptr;KSPGetPC(solver,&pc);PCSetType(pc,PCLU);
	PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);KSPSetFromOptions(solver);
	PCType pc_type=nullptr;
	NativeTetAlePetscCheck(PCGetType(pc,&pc_type),"PCGetType runtime");
	PetscBool use_user_schur=PETSC_FALSE;
	if(pc_type&&std::string(pc_type)==PCFIELDSPLIT) {
		IS velocity_is=nullptr,pressure_is=nullptr;
		PetscInt matrix_begin=0,matrix_end=0;
		NativeTetAlePetscCheck(MatGetOwnershipRange(matrix,&matrix_begin,&matrix_end),
			"MatGetOwnershipRange fieldsplit runtime");
		const PetscInt velocity_end=static_cast<PetscInt>(3*velocity_nodes);
		std::vector<PetscInt> velocity_indices,pressure_indices;
		velocity_indices.reserve(static_cast<std::size_t>(matrix_end-matrix_begin));
		pressure_indices.reserve(static_cast<std::size_t>(matrix_end-matrix_begin));
		for(PetscInt row=matrix_begin;row<matrix_end;++row)
			(row<velocity_end?velocity_indices:pressure_indices).push_back(row);
		NativeTetAlePetscCheck(ISCreateGeneral(PETSC_COMM_WORLD,
			static_cast<PetscInt>(velocity_indices.size()),velocity_indices.data(),
			PETSC_COPY_VALUES,&velocity_is),"ISCreateGeneral velocity runtime");
		NativeTetAlePetscCheck(ISCreateGeneral(PETSC_COMM_WORLD,
			static_cast<PetscInt>(pressure_indices.size()),pressure_indices.data(),
			PETSC_COPY_VALUES,&pressure_is),"ISCreateGeneral pressure runtime");
		NativeTetAlePetscCheck(PCFieldSplitSetIS(pc,"velocity",velocity_is),
			"PCFieldSplitSetIS velocity runtime");
		NativeTetAlePetscCheck(PCFieldSplitSetIS(pc,"pressure",pressure_is),
			"PCFieldSplitSetIS pressure runtime");
		NativeTetAlePetscCheck(PCSetUseAmat(pc,PETSC_TRUE),"PCSetUseAmat runtime");
		use_user_schur=PETSC_TRUE;
		PetscOptionsGetBool(nullptr,nullptr,"-native_ale_user_schur",
			&use_user_schur,nullptr);
		if(use_user_schur) {
			const PetscInt pressure_dofs=static_cast<PetscInt>(
				current_mesh.points.size()+controller_dofs);
			const PetscInt local_pressure_dofs=static_cast<PetscInt>(pressure_indices.size());
			NativeTetAlePetscCheck(MatCreateAIJ(PETSC_COMM_WORLD,local_pressure_dofs,
				local_pressure_dofs,pressure_dofs,pressure_dofs,24,nullptr,24,nullptr,
				&schur_preconditioner),"MatCreateAIJ Schur runtime");
			NativeTetAlePetscCheck(MatSetOption(schur_preconditioner,
				MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE),"MatSetOption Schur runtime");
			for(std::size_t cell=static_cast<std::size_t>(rank);cell<current_mesh.cells.size();
					cell+=static_cast<std::size_t>(ranks)) {
				const auto geometry=EvaluateNativeTetGeometry(current_mesh,current_mesh.cells[cell]);
				std::array<PetscInt,4> rows{};std::array<PetscScalar,16> values{};
				for(std::size_t i=0;i<4;++i) {
					rows[i]=static_cast<PetscInt>(current_mesh.cells[cell].nodes[i]);
					for(std::size_t j=0;j<4;++j)
						values[4*i+j]=-geometry.determinant/(i==j?60.0:120.0);
				}
				NativeTetAlePetscCheck(MatSetValues(schur_preconditioner,4,rows.data(),
					4,rows.data(),values.data(),ADD_VALUES),"MatSetValues Schur runtime");
			}
			NativeTetAlePetscCheck(MatAssemblyBegin(schur_preconditioner,MAT_FINAL_ASSEMBLY),
				"MatAssemblyBegin Schur runtime");
			NativeTetAlePetscCheck(MatAssemblyEnd(schur_preconditioner,MAT_FINAL_ASSEMBLY),
				"MatAssemblyEnd Schur runtime");
			for(std::size_t control=0;control<controller_dofs;++control){
				const PetscInt row=static_cast<PetscInt>(current_mesh.points.size()+control);
				NativeTetAlePetscCheck(MatSetValue(schur_preconditioner,row,row,-1.0,
					INSERT_VALUES),"MatSetValue controller Schur runtime");
			}
			NativeTetAlePetscCheck(MatAssemblyBegin(schur_preconditioner,MAT_FINAL_ASSEMBLY),
				"MatAssemblyBegin controller Schur runtime");
			NativeTetAlePetscCheck(MatAssemblyEnd(schur_preconditioner,MAT_FINAL_ASSEMBLY),
				"MatAssemblyEnd controller Schur runtime");
			NativeTetAlePetscCheck(PCFieldSplitSetSchurPre(pc,PC_FIELDSPLIT_SCHUR_PRE_USER,
				schur_preconditioner),"PCFieldSplitSetSchurPre runtime");
		}
		ISDestroy(&velocity_is);ISDestroy(&pressure_is);
	}
	NativeTetAlePetscSolveResult result;
	KSPType ksp_type=nullptr;KSPGetType(solver,&ksp_type);
	result.linear_solver_type=ksp_type?ksp_type:"unknown";
	result.preconditioner_type=pc_type?pc_type:"unknown";
	result.system_scaling_enabled=scale_system==PETSC_TRUE;
	result.user_schur_enabled=use_user_schur==PETSC_TRUE;
	bool line_search_pending=false;
	double base_residual_norm=0.0,step_length=1.0;
	for(std::size_t iteration=0;iteration<=maximum_iterations;++iteration) {
		VecScatterBegin(scatter,state,replicated,INSERT_VALUES,SCATTER_FORWARD);
		VecScatterEnd(scatter,state,replicated,INSERT_VALUES,SCATTER_FORWARD);
		const PetscScalar* global=nullptr;VecGetArrayRead(replicated,&global);
		auto builder=[&](std::size_t cell,const auto& rows) {
			std::array<double,NativeTaylorHoodElementSystem::dofs> trial{},committed{};
			for(std::size_t entry=0;entry<rows.size();++entry) {
				trial[entry]=PetscRealPart(global[rows[entry]]);
				committed[entry]=committed_state[static_cast<std::size_t>(rows[entry])];
			}
			std::array<std::array<double,3>,4> grid{};
			for(std::size_t node=0;node<4;++node)
				grid[node]=vertex_mesh_velocity_m_s[current_mesh.cells[cell].nodes[node]];
			if(steady) return BuildNativeTaylorHoodAleNavierStokesElement(current_mesh,
				current_mesh.cells[cell],trial,grid,parameters);
			return BuildNativeTaylorHoodAleTransientElement(current_mesh,
				current_mesh.cells[cell],trial,committed,grid,parameters,dt_s);
		};
		AssembleNativeTetAlePetscElements(current_mesh,topology,rank,ranks,builder,
			matrix,residual);
		for(std::size_t face=static_cast<std::size_t>(rank);
			face<current_mesh.boundary_triangles.size();face+=static_cast<std::size_t>(ranks)){
			const auto& triangle=current_mesh.boundary_triangles[face];
			const auto stabilization=boundary_conditions.backflow_stabilization_beta.find(
				triangle.boundary_label);
			if(stabilization==boundary_conditions.backflow_stabilization_beta.end()
				||stabilization->second==0.0)continue;
			const std::size_t cell=backflow_owners[face];
			std::array<double,30> local_velocity{};
			std::array<std::array<double,3>,4> local_grid{};
			std::array<PetscInt,30> rows{};
			for(std::size_t node=0;node<10;++node)
				for(int component=0;component<3;++component){
					const std::size_t local=3*node+component;
					rows[local]=static_cast<PetscInt>(3*topology.cell_velocity_nodes[cell][node]
						+component);
					local_velocity[local]=PetscRealPart(global[rows[local]]);
				}
			for(std::size_t node=0;node<4;++node)
				local_grid[node]=vertex_mesh_velocity_m_s[current_mesh.cells[cell].nodes[node]];
			const auto contribution=BuildNativeTetAleBackflowFaceSystem(current_mesh,
				topology,cell,triangle,local_velocity,local_grid,parameters.density,
				stabilization->second);
			NativeTetAlePetscCheck(VecSetValues(residual,30,rows.data(),
				contribution.residual.data(),ADD_VALUES),"VecSetValues backflow runtime");
			NativeTetAlePetscCheck(MatSetValues(matrix,30,rows.data(),30,rows.data(),
				contribution.jacobian.data(),ADD_VALUES),"MatSetValues backflow runtime");
		}
		if(rank==0){
			for(const auto& prescribed:boundary_conditions.prescribed_pressure_pa){
				const auto& coefficients=boundary_operators.at(prescribed.first)
					.velocity_coefficients;
				for(std::size_t dof=0;dof<coefficients.size();++dof)
					if(coefficients[dof]!=0.0)
						NativeTetAlePetscCheck(VecSetValue(residual,static_cast<PetscInt>(dof),
							prescribed.second*coefficients[dof],ADD_VALUES),
							"VecSetValue prescribed pressure runtime");
			}
			for(std::size_t control=0;control<controller_dofs;++control){
				const auto& specification=boundary_conditions.flow_rate_controls[control];
				const auto& coefficients=boundary_operators.at(specification.boundary_label)
					.velocity_coefficients;
				const PetscInt controller_row=static_cast<PetscInt>(physical_dofs_size+control);
				const double multiplier=PetscRealPart(global[controller_row]);
				double flow=-specification.target_outward_flow_m3_s;
				for(std::size_t dof=0;dof<coefficients.size();++dof)
					if(coefficients[dof]!=0.0){
						const PetscInt velocity_row=static_cast<PetscInt>(dof);
						flow+=coefficients[dof]*PetscRealPart(global[velocity_row]);
						NativeTetAlePetscCheck(VecSetValue(residual,velocity_row,
							multiplier*coefficients[dof],ADD_VALUES),
							"VecSetValue flow controller traction runtime");
						NativeTetAlePetscCheck(MatSetValue(matrix,velocity_row,controller_row,
							coefficients[dof],ADD_VALUES),
							"MatSetValue flow controller column runtime");
						NativeTetAlePetscCheck(MatSetValue(matrix,controller_row,velocity_row,
							coefficients[dof],ADD_VALUES),
							"MatSetValue flow controller row runtime");
					}
				NativeTetAlePetscCheck(VecSetValue(residual,controller_row,flow,ADD_VALUES),
					"VecSetValue flow controller residual runtime");
			}
		}
		VecRestoreArrayRead(replicated,&global);
		NativeTetAlePetscCheck(MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY),
			"MatAssemblyBegin boundary controls runtime");
		NativeTetAlePetscCheck(MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY),
			"MatAssemblyEnd boundary controls runtime");
		NativeTetAlePetscCheck(VecAssemblyBegin(residual),
			"VecAssemblyBegin boundary controls runtime");
		NativeTetAlePetscCheck(VecAssemblyEnd(residual),
			"VecAssemblyEnd boundary controls runtime");
		NativeTetAlePetscCheck(MatZeroRowsColumns(matrix,
			static_cast<PetscInt>(local_constraint_rows.size()),local_constraint_rows.data(),
			1.0,nullptr,nullptr),"MatZeroRowsColumns runtime");
		std::vector<PetscScalar> values(local_constraint_rows.size());
		VecGetValues(state,static_cast<PetscInt>(local_constraint_rows.size()),
			local_constraint_rows.data(),values.data());
		for(std::size_t entry=0;entry<values.size();++entry)
			values[entry]-=local_constraint_values[entry];
		VecSetValues(residual,static_cast<PetscInt>(local_constraint_rows.size()),
			local_constraint_rows.data(),values.data(),INSERT_VALUES);
		VecAssemblyBegin(residual);VecAssemblyEnd(residual);
		PetscReal norm=0.0;VecNorm(residual,NORM_2,&norm);
		if(line_search_pending){
			if(std::isfinite(static_cast<double>(norm))&&norm<base_residual_norm)
				line_search_pending=false;
			else {
				step_length*=0.5;
				if(step_length<1.0/16384.0){
					PetscInt residual_begin=0,residual_end=0;
					VecGetOwnershipRange(residual,&residual_begin,&residual_end);
					const PetscScalar* residual_values=nullptr;
					VecGetArrayRead(residual,&residual_values);
					double local_squares[3]{},global_squares[3]{};
					for(PetscInt row=residual_begin;row<residual_end;++row){
						const auto value=static_cast<double>(PetscRealPart(
							residual_values[row-residual_begin]));
						const int group=row<static_cast<PetscInt>(3*velocity_nodes)?0
							:row<static_cast<PetscInt>(physical_dofs_size)?1:2;
						local_squares[group]+=value*value;
					}
					VecRestoreArrayRead(residual,&residual_values);
					MPI_Allreduce(local_squares,global_squares,3,MPI_DOUBLE,MPI_SUM,
						PETSC_COMM_WORLD);
					std::ostringstream diagnostic;
					diagnostic<<std::setprecision(17)
						<<"native ALE PETSc Newton line search failed: iteration="
						<<iteration<<" base_residual="<<base_residual_norm
						<<" trial_residual="<<norm<<" rejected_step="<<step_length
						<<" velocity_residual="<<std::sqrt(global_squares[0])
						<<" continuity_residual="<<std::sqrt(global_squares[1])
						<<" controller_residual="<<std::sqrt(global_squares[2])
						<<" flow_controllers="<<controller_dofs;
					throw std::runtime_error(diagnostic.str());
				}
				NativeTetAlePetscCheck(VecCopy(base_state,state),
					"VecCopy Newton rollback state");
				NativeTetAlePetscCheck(VecAXPY(state,step_length,update),
					"VecAXPY Newton backtrack");
				--iteration;
				continue;
			}
		}
		result.final_residual_l2=norm;result.newton_iterations=iteration;
		if(norm<=absolute_tolerance) break;
		if(iteration==maximum_iterations)
			throw std::runtime_error("native ALE PETSc Newton did not converge: residual="
				+std::to_string(static_cast<double>(norm)));
		if(scale_system) {
			NativeTetAlePetscCheck(MatDiagonalScale(matrix,row_scale,column_scale),
				"MatDiagonalScale runtime");
			NativeTetAlePetscCheck(VecPointwiseMult(right,residual,row_scale),
				"VecPointwiseMult residual scale runtime");
		} else VecCopy(residual,right);
		VecScale(right,-1.0);KSPSetOperators(solver,matrix,matrix);
		NativeTetAlePetscCheck(KSPSolve(solver,right,update),"KSPSolve runtime");
		KSPConvergedReason reason;KSPGetConvergedReason(solver,&reason);
		PetscInt iterations=0;KSPGetIterationNumber(solver,&iterations);
		result.final_linear_converged_reason=static_cast<PetscInt>(reason);
		if(reason<=0) throw std::runtime_error("native ALE PETSc linear solve did not converge: "
			"reason="+std::to_string(static_cast<int>(reason))
			+" iterations="+std::to_string(static_cast<long long>(iterations)));
		result.total_linear_iterations+=iterations;
		if(scale_system) NativeTetAlePetscCheck(VecPointwiseMult(update,update,column_scale),
			"VecPointwiseMult update scale runtime");
		NativeTetAlePetscCheck(VecCopy(state,base_state),"VecCopy Newton base state");
		base_residual_norm=norm;step_length=1.0;line_search_pending=true;
		NativeTetAlePetscCheck(VecAXPY(state,step_length,update),"VecAXPY Newton update");
	}
	VecScatterBegin(scatter,state,replicated,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd(scatter,state,replicated,INSERT_VALUES,SCATTER_FORWARD);
	const PetscScalar* final_values=nullptr;VecGetArrayRead(replicated,&final_values);
	result.replicated_state.resize(physical_dofs_size);
	for(std::size_t dof=0;dof<physical_dofs_size;++dof)
		result.replicated_state[dof]=PetscRealPart(final_values[dof]);
	result.flow_controller_multipliers_pa.resize(controller_dofs);
	for(std::size_t control=0;control<controller_dofs;++control)
		result.flow_controller_multipliers_pa[control]=
			PetscRealPart(final_values[physical_dofs_size+control]);
	VecRestoreArrayRead(replicated,&final_values);
	VecScatterDestroy(&scatter);VecDestroy(&replicated);KSPDestroy(&solver);
	MatDestroy(&schur_preconditioner);
	VecDestroy(&column_scale);VecDestroy(&row_scale);VecDestroy(&base_state);
	VecDestroy(&update);VecDestroy(&right);
	VecDestroy(&residual);VecDestroy(&state);
	MatDestroy(&matrix);return result;
}

inline NativeTetAlePetscSolveResult SolveNativeTetAlePetscTransient(
	const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& vertex_mesh_velocity_m_s,
	const std::vector<double>& committed_state,const std::vector<double>& initial_trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double dt_s,double absolute_tolerance=1e-10,std::size_t maximum_iterations=10,
	const NativeTetAleBoundaryConditions& boundary_conditions={})
{
	return SolveNativeTetAlePetsc(current_mesh,vertex_mesh_velocity_m_s,committed_state,
		initial_trial_state,prescribed_velocity,pressure_gauge_node,parameters,dt_s,false,
		absolute_tolerance,maximum_iterations,boundary_conditions);
}

inline NativeTetAlePetscSolveResult SolveNativeTetAlePetscSteady(
	const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& vertex_mesh_velocity_m_s,
	const std::vector<double>& initial_trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double absolute_tolerance=1e-10,std::size_t maximum_iterations=10,
	const NativeTetAleBoundaryConditions& boundary_conditions={})
{
	return SolveNativeTetAlePetsc(current_mesh,vertex_mesh_velocity_m_s,initial_trial_state,
		initial_trial_state,prescribed_velocity,pressure_gauge_node,parameters,1.0,true,
		absolute_tolerance,maximum_iterations,boundary_conditions);
}

} // namespace iga

#endif

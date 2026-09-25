#ifndef IGA_NATIVE_TET_ALE_DENSE_RUNTIME_HPP
#define IGA_NATIVE_TET_ALE_DENSE_RUNTIME_HPP

#include "NativeTetAleTransient.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

struct NativeTetAleDenseSolveResult
{
	std::vector<double> state;
	std::size_t newton_iterations = 0;
	double final_free_residual_l2 = 0.0;
};

namespace native_tet_ale_dense_detail {

inline std::vector<double> Solve(std::vector<double> matrix,std::vector<double> right)
{
	const std::size_t size=right.size();
	if(matrix.size()!=size*size)
		throw std::invalid_argument("native ALE dense matrix size is invalid");
	for(std::size_t column=0;column<size;++column) {
		std::size_t pivot=column;
		for(std::size_t row=column+1;row<size;++row)
			if(std::abs(matrix[row*size+column])>std::abs(matrix[pivot*size+column]))
				pivot=row;
		if(!(std::abs(matrix[pivot*size+column])>
			100.0*std::numeric_limits<double>::epsilon()))
			throw std::runtime_error("native ALE dense Jacobian is singular");
		if(pivot!=column) {
			for(std::size_t entry=column;entry<size;++entry)
				std::swap(matrix[column*size+entry],matrix[pivot*size+entry]);
			std::swap(right[column],right[pivot]);
		}
		for(std::size_t row=column+1;row<size;++row) {
			const double factor=matrix[row*size+column]/matrix[column*size+column];
			matrix[row*size+column]=0.0;
			for(std::size_t entry=column+1;entry<size;++entry)
				matrix[row*size+entry]-=factor*matrix[column*size+entry];
			right[row]-=factor*right[column];
		}
	}
	std::vector<double> solution(size,0.0);
	for(std::size_t reverse=size;reverse>0;--reverse) {
		const std::size_t row=reverse-1;double value=right[row];
		for(std::size_t column=row+1;column<size;++column)
			value-=matrix[row*size+column]*solution[column];
		solution[row]=value/matrix[row*size+row];
	}
	return solution;
}

inline std::array<std::size_t,NativeTaylorHoodElementSystem::dofs> ElementDofs(
	const NativeTetMesh& mesh,const NativeTaylorHoodTopology& topology,std::size_t cell)
{
	std::array<std::size_t,NativeTaylorHoodElementSystem::dofs> result{};
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	for(std::size_t local=0;local<10;++local)
		for(int component=0;component<3;++component)
			result[3*local+component]=3*topology.cell_velocity_nodes[cell][local]+component;
	for(std::size_t local=0;local<4;++local)
		result[30+local]=3*velocity_nodes+mesh.cells[cell].nodes[local];
	return result;
}

} // namespace native_tet_ale_dense_detail

// Small-case reference runtime. It exercises global native assembly and
// nonlinear lifecycle; production-size cases should replace only the dense
// algebra with the PETSc/MPI path.
inline NativeTetAleDenseSolveResult SolveNativeTetAleDenseTransient(
	const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& vertex_mesh_velocity_m_s,
	const std::vector<double>& committed_state,std::vector<double> trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double dt_s,double absolute_tolerance=1.0e-10,std::size_t maximum_iterations=12)
{
	const auto topology=BuildNativeTaylorHoodTopology(current_mesh);
	const std::size_t velocity_nodes=current_mesh.points.size()+topology.edges.size();
	const std::size_t dofs=3*velocity_nodes+current_mesh.points.size();
	if(vertex_mesh_velocity_m_s.size()!=current_mesh.points.size()
		||committed_state.size()!=dofs||trial_state.size()!=dofs
		||pressure_gauge_node>=current_mesh.points.size()
		||!(absolute_tolerance>0.0)||maximum_iterations==0)
		throw std::invalid_argument("native ALE dense runtime input is invalid");
	std::map<std::size_t,double> constraints;
	for(const auto& prescribed:prescribed_velocity) {
		if(prescribed.first>=velocity_nodes)
			throw std::invalid_argument("native ALE prescribed velocity node is out of range");
		for(int component=0;component<3;++component) {
			if(!std::isfinite(prescribed.second[component]))
				throw std::invalid_argument("native ALE prescribed velocity is nonfinite");
			constraints.emplace(3*prescribed.first+component,prescribed.second[component]);
		}
	}
	constraints.emplace(3*velocity_nodes+pressure_gauge_node,0.0);
	for(const auto& constraint:constraints) trial_state[constraint.first]=constraint.second;
	NativeTetAleDenseSolveResult result;
	for(std::size_t iteration=0;iteration<=maximum_iterations;++iteration) {
		std::vector<double> residual(dofs,0.0),jacobian(dofs*dofs,0.0);
		for(std::size_t cell=0;cell<current_mesh.cells.size();++cell) {
			const auto element_dofs=native_tet_ale_dense_detail::ElementDofs(
				current_mesh,topology,cell);
			std::array<double,NativeTaylorHoodElementSystem::dofs> local_trial{},local_committed{};
			for(std::size_t local=0;local<element_dofs.size();++local) {
				local_trial[local]=trial_state[element_dofs[local]];
				local_committed[local]=committed_state[element_dofs[local]];
			}
			std::array<std::array<double,3>,4> local_mesh_velocity{};
			for(std::size_t local=0;local<4;++local)
				local_mesh_velocity[local]=vertex_mesh_velocity_m_s[current_mesh.cells[cell].nodes[local]];
			const auto element=BuildNativeTaylorHoodAleTransientElement(current_mesh,
				current_mesh.cells[cell],local_trial,local_committed,local_mesh_velocity,
				parameters,dt_s);
			for(std::size_t row=0;row<element_dofs.size();++row) {
				residual[element_dofs[row]]+=element.residual[row];
				for(std::size_t column=0;column<element_dofs.size();++column)
					jacobian[element_dofs[row]*dofs+element_dofs[column]]
						+=element.jacobian[row*element_dofs.size()+column];
			}
		}
		for(const auto& constraint:constraints) {
			const std::size_t dof=constraint.first;
			for(std::size_t other=0;other<dofs;++other) {
				jacobian[dof*dofs+other]=0.0;jacobian[other*dofs+dof]=0.0;
			}
			jacobian[dof*dofs+dof]=1.0;
			residual[dof]=trial_state[dof]-constraint.second;
		}
		double norm_squared=0.0;
		for(const auto value:residual) norm_squared+=value*value;
		result.final_free_residual_l2=std::sqrt(norm_squared);
		result.newton_iterations=iteration;
		if(result.final_free_residual_l2<=absolute_tolerance) {
			result.state=std::move(trial_state);return result;
		}
		if(iteration==maximum_iterations)
			throw std::runtime_error("native ALE dense Newton did not converge");
		for(auto& value:residual) value=-value;
		const auto update=native_tet_ale_dense_detail::Solve(std::move(jacobian),
			std::move(residual));
		for(std::size_t dof=0;dof<dofs;++dof) trial_state[dof]+=update[dof];
	}
	throw std::logic_error("native ALE dense runtime reached unreachable state");
}

// Small-case steady reference solve used for operator verification, including
// manufactured body-force cases. Production meshes use the PETSc/MPI algebra.
inline NativeTetAleDenseSolveResult SolveNativeTetDenseSteady(
	const NativeTetMesh& mesh,std::vector<double> trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double absolute_tolerance=1.0e-10,std::size_t maximum_iterations=12)
{
	const auto topology=BuildNativeTaylorHoodTopology(mesh);
	const std::size_t velocity_nodes=mesh.points.size()+topology.edges.size();
	const std::size_t dofs=3*velocity_nodes+mesh.points.size();
	if(trial_state.size()!=dofs||pressure_gauge_node>=mesh.points.size()
		||!(absolute_tolerance>0.0)||maximum_iterations==0)
		throw std::invalid_argument("native steady dense runtime input is invalid");
	std::map<std::size_t,double> constraints;
	for(const auto& prescribed:prescribed_velocity){
		if(prescribed.first>=velocity_nodes)
			throw std::invalid_argument("native steady prescribed velocity node is out of range");
		for(int component=0;component<3;++component){
			if(!std::isfinite(prescribed.second[component]))
				throw std::invalid_argument("native steady prescribed velocity is nonfinite");
			constraints.emplace(3*prescribed.first+component,prescribed.second[component]);
		}
	}
	constraints.emplace(3*velocity_nodes+pressure_gauge_node,0.0);
	for(const auto& constraint:constraints)trial_state[constraint.first]=constraint.second;
	NativeTetAleDenseSolveResult result;
	for(std::size_t iteration=0;iteration<=maximum_iterations;++iteration){
		std::vector<double> residual(dofs,0.0),jacobian(dofs*dofs,0.0);
		for(std::size_t cell=0;cell<mesh.cells.size();++cell){
			const auto element_dofs=native_tet_ale_dense_detail::ElementDofs(mesh,topology,cell);
			std::array<double,NativeTaylorHoodElementSystem::dofs> local{};
			for(std::size_t i=0;i<element_dofs.size();++i)local[i]=trial_state[element_dofs[i]];
			const auto element=BuildNativeTaylorHoodNavierStokesElement(mesh,mesh.cells[cell],local,parameters);
			for(std::size_t row=0;row<element_dofs.size();++row){
				residual[element_dofs[row]]+=element.residual[row];
				for(std::size_t column=0;column<element_dofs.size();++column)
					jacobian[element_dofs[row]*dofs+element_dofs[column]]
						+=element.jacobian[row*element_dofs.size()+column];
			}
		}
		for(const auto& constraint:constraints){
			const std::size_t dof=constraint.first;
			for(std::size_t other=0;other<dofs;++other){
				jacobian[dof*dofs+other]=0.0;jacobian[other*dofs+dof]=0.0;
			}
			jacobian[dof*dofs+dof]=1.0;residual[dof]=trial_state[dof]-constraint.second;
		}
		double norm_squared=0.0;for(double value:residual)norm_squared+=value*value;
		result.final_free_residual_l2=std::sqrt(norm_squared);result.newton_iterations=iteration;
		if(result.final_free_residual_l2<=absolute_tolerance){result.state=std::move(trial_state);return result;}
		if(iteration==maximum_iterations)
			throw std::runtime_error("native steady dense Newton did not converge; residual="
				+std::to_string(result.final_free_residual_l2));
		for(double& value:residual)value=-value;
		const auto update=native_tet_ale_dense_detail::Solve(std::move(jacobian),std::move(residual));
		for(std::size_t dof=0;dof<dofs;++dof)trial_state[dof]+=update[dof];
	}
	throw std::logic_error("native steady dense runtime reached unreachable state");
}

} // namespace iga

#endif

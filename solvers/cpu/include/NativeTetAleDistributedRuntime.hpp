#ifndef IGA_NATIVE_TET_ALE_DISTRIBUTED_RUNTIME_HPP
#define IGA_NATIVE_TET_ALE_DISTRIBUTED_RUNTIME_HPP

// P2/P1 ALE Navier--Stokes Newton solve on a partitioned tetrahedral mesh.
// Each rank assembles its owned cells and owned boundary faces into a PETSc
// system laid out by NativeTetDistributedDofs, keeps a ghosted state vector,
// and never needs another rank's cells. Inputs and the returned state use the
// replicated layout of the serial runtime (velocity of every vertex and edge
// node, then pressure of every vertex), so callers can switch between the two
// runtimes; the replicated result costs one Allreduce of that layout.
// Flow-rate controllers and system scaling are not supported here.

#include "NativeTetAleBoundaryControl.hpp"
#include "NativeTetAleBackflow.hpp"
#include "NativeTetAleTransient.hpp"
#include "NativeTetAlePetscAssembly.hpp"
#include "NativeTetDistributedDofs.hpp"
#include "NativeTetMeshPartition.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct NativeTetAleDistributedResult
{
	std::vector<double> replicated_state;
	std::size_t newton_iterations=0;
	double final_residual_l2=0.0;
	PetscInt total_linear_iterations=0;
	std::string linear_solver_type,preconditioner_type;
	PetscInt owned_rows=0,global_rows=0;
};

// Partition and row layout of one reference mesh, reused across time steps.
struct NativeTetAleDistributedMesh
{
	NativeTetMeshPartition partition;
	NativeTetDistributedDofs dofs;
	MPI_Comm communicator=MPI_COMM_NULL;
	std::vector<std::uint32_t> global_to_local_vertex;        // global vertex -> local or max
	std::map<std::array<std::uint32_t,2>,std::uint32_t> local_edge_of;  // (min,max) global -> local edge

	NativeTetAleDistributedMesh(const NativeTetMesh& reference,MPI_Comm comm)
		:partition(PartitionNativeTetMesh(reference,comm)),dofs(BuildNativeTetDistributedDofs(partition,comm)),
		 communicator(comm),global_to_local_vertex(reference.points.size(),std::numeric_limits<std::uint32_t>::max())
	{
		for(std::size_t local=0;local<partition.local_to_global_node.size();++local)
			global_to_local_vertex[partition.local_to_global_node[local]]=static_cast<std::uint32_t>(local);
		for(std::size_t edge=0;edge<dofs.edge_global.size();++edge)
			local_edge_of.emplace(dofs.edge_global[edge],static_cast<std::uint32_t>(edge));
	}
};

namespace native_tet_ale_distributed_detail {

struct Objects
{
	Mat matrix=nullptr;
	Vec state=nullptr,residual=nullptr,right=nullptr,update=nullptr,base_state=nullptr;
	KSP solver=nullptr;
	~Objects()
	{
		if(solver)KSPDestroy(&solver);
		if(base_state)VecDestroy(&base_state);
		if(update)VecDestroy(&update);
		if(right)VecDestroy(&right);
		if(residual)VecDestroy(&residual);
		if(state)VecDestroy(&state);
		if(matrix)MatDestroy(&matrix);
	}
};

} // namespace native_tet_ale_distributed_detail

inline NativeTetAleDistributedResult SolveNativeTetAleDistributed(
	const NativeTetAleDistributedMesh& distributed,
	const NativeTaylorHoodTopology& replicated_topology,
	const std::vector<std::array<double,3>>& current_points_m,
	const std::vector<std::array<double,3>>& vertex_mesh_velocity_m_s,
	const std::vector<double>& committed_state,const std::vector<double>& initial_trial_state,
	const std::map<std::uint32_t,std::array<double,3>>& prescribed_velocity,
	std::uint32_t pressure_gauge_node,const NativeNavierStokesParameters& parameters,
	double dt_s,bool steady,double absolute_tolerance=1e-10,std::size_t maximum_iterations=10,
	const NativeTetAleBoundaryConditions& boundary_conditions={})
{
	using namespace native_tet_ale_distributed_detail;
	const auto& partition=distributed.partition;
	const auto& dofs=distributed.dofs;
	const MPI_Comm comm=distributed.communicator;
	const int rank=partition.rank;
	const std::size_t global_vertices=partition.node_owner.size();
	const std::size_t replicated_velocity_nodes=global_vertices+replicated_topology.edges.size();
	const std::size_t replicated_dofs=3*replicated_velocity_nodes+global_vertices;
	const bool gauge_enabled=pressure_gauge_node!=std::numeric_limits<std::uint32_t>::max();
	if(!boundary_conditions.flow_rate_controls.empty())
		throw std::invalid_argument("native ALE distributed runtime does not support flow-rate controllers");
	// Every prescribed-pressure or backflow label must exist on the global
	// boundary: collect the labels of owned faces and reduce.
	{
		std::vector<int> labels;
		for(const auto& item:boundary_conditions.prescribed_pressure_pa)labels.push_back(item.first);
		for(const auto& item:boundary_conditions.backflow_stabilization_beta)labels.push_back(item.first);
		std::vector<int> present(labels.size(),0),anywhere(labels.size(),0);
		for(std::size_t triangle=0;triangle<partition.local.boundary_triangles.size();++triangle){
			if(!partition.local_triangle_owned[triangle])continue;
			for(std::size_t index=0;index<labels.size();++index)
				if(partition.local.boundary_triangles[triangle].boundary_label==labels[index])present[index]=1;
		}
		if(!labels.empty()&&MPI_Allreduce(present.data(),anywhere.data(),static_cast<int>(labels.size()),MPI_INT,MPI_MAX,comm)!=MPI_SUCCESS)
			throw std::runtime_error("native ALE distributed runtime: boundary label reduction failed");
		for(std::size_t index=0;index<labels.size();++index)if(!anywhere[index])
			throw std::invalid_argument("native ALE distributed runtime: boundary label "+std::to_string(labels[index])
				+" has no faces");
		for(const auto& item:boundary_conditions.prescribed_pressure_pa)if(!std::isfinite(item.second))
			throw std::invalid_argument("native ALE distributed runtime: prescribed pressure is nonfinite");
		for(const auto& item:boundary_conditions.backflow_stabilization_beta)if(!std::isfinite(item.second)||item.second<0.0)
			throw std::invalid_argument("native ALE distributed runtime: backflow beta is invalid");
	}
	if(gauge_enabled&&!boundary_conditions.prescribed_pressure_pa.empty())
		throw std::invalid_argument("native ALE pressure outlet and point gauge conflict");
	if(!gauge_enabled&&boundary_conditions.prescribed_pressure_pa.empty())
		throw std::invalid_argument("native ALE fluid pressure needs a gauge or outlet");
	if(current_points_m.size()!=global_vertices||vertex_mesh_velocity_m_s.size()!=global_vertices
		||committed_state.size()!=replicated_dofs||initial_trial_state.size()!=replicated_dofs
		||(gauge_enabled&&pressure_gauge_node>=global_vertices)||(!steady&&!(dt_s>0.0))
		||!(absolute_tolerance>0.0)||maximum_iterations==0)
		throw std::invalid_argument("native ALE distributed runtime input is invalid");
	PetscBool scale_system=PETSC_FALSE;
	PetscOptionsGetBool(nullptr,nullptr,"-native_ale_scale",&scale_system,nullptr);
	if(scale_system)throw std::invalid_argument("native ALE distributed runtime does not support -native_ale_scale");
	// Current local mesh and the replicated-layout index of every local unknown.
	NativeTetMesh current=partition.local;
	for(std::size_t local=0;local<current.points.size();++local)
		current.points[local]=current_points_m[partition.local_to_global_node[local]];
	std::map<std::array<std::uint32_t,2>,std::uint32_t> replicated_edge_index;
	for(std::size_t edge=0;edge<replicated_topology.edges.size();++edge){
		const auto& ends=replicated_topology.edges[edge];
		replicated_edge_index.emplace(std::array<std::uint32_t,2>{{std::min(ends[0],ends[1]),std::max(ends[0],ends[1])}},
			static_cast<std::uint32_t>(edge));
	}
	const std::size_t local_vertices=current.points.size(),local_edges=dofs.topology.edges.size();
	std::vector<std::size_t> vertex_replicated_velocity(local_vertices),vertex_replicated_pressure(local_vertices),
		edge_replicated_velocity(local_edges);
	for(std::size_t local=0;local<local_vertices;++local){
		const auto global=partition.local_to_global_node[local];
		vertex_replicated_velocity[local]=3*static_cast<std::size_t>(global);
		vertex_replicated_pressure[local]=3*replicated_velocity_nodes+global;
	}
	for(std::size_t edge=0;edge<local_edges;++edge){
		const auto found=replicated_edge_index.find(dofs.edge_global[edge]);
		if(found==replicated_edge_index.end())
			throw std::invalid_argument("native ALE distributed runtime: local edge is absent from the replicated topology");
		edge_replicated_velocity[edge]=3*(global_vertices+found->second);
	}
	const auto cell_replicated=[&](std::size_t cell){
		std::array<std::size_t,34> rows{};std::size_t entry=0;
		for(const auto node:dofs.topology.cell_velocity_nodes[cell]){
			const std::size_t base=node<local_vertices?vertex_replicated_velocity[node]:edge_replicated_velocity[node-local_vertices];
			for(int component=0;component<3;++component)rows[entry++]=base+component;
		}
		for(const auto node:current.cells[cell].nodes)rows[entry++]=vertex_replicated_pressure[node];
		return rows;
	};
	// Matrix with exact preallocation of the owned rows from every local cell.
	const PetscInt begin=dofs.row_begin,end=dofs.row_end,local_rows=end-begin,global_rows=dofs.global_rows;
	Objects petsc;
	{
		std::vector<std::vector<PetscInt>> columns(static_cast<std::size_t>(local_rows));
		for(PetscInt row=begin;row<end;++row)columns[row-begin].push_back(row);
		for(std::size_t cell=0;cell<current.cells.size();++cell){
			const auto rows=dofs.CellRows(cell);
			for(const auto row:rows)if(row>=begin&&row<end){
				auto& entries=columns[row-begin];entries.insert(entries.end(),rows.begin(),rows.end());}
		}
		std::vector<PetscInt> diagonal(local_rows),off_diagonal(local_rows);
		for(std::size_t row=0;row<columns.size();++row){
			auto& entries=columns[row];std::sort(entries.begin(),entries.end());
			entries.erase(std::unique(entries.begin(),entries.end()),entries.end());
			for(const auto column:entries)++(column>=begin&&column<end?diagonal[row]:off_diagonal[row]);
		}
		NativeTetAlePetscCheck(MatCreateAIJ(comm,local_rows,local_rows,global_rows,global_rows,0,diagonal.data(),
			0,off_diagonal.data(),&petsc.matrix),"MatCreateAIJ distributed");
		NativeTetAlePetscCheck(MatSetOption(petsc.matrix,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE),"MatSetOption distributed");
	}
	NativeTetAlePetscCheck(VecCreateGhost(comm,local_rows,global_rows,static_cast<PetscInt>(dofs.ghost_rows.size()),
		dofs.ghost_rows.data(),&petsc.state),"VecCreateGhost distributed");
	NativeTetAlePetscCheck(VecCreateMPI(comm,local_rows,global_rows,&petsc.residual),"VecCreateMPI residual distributed");
	NativeTetAlePetscCheck(VecDuplicate(petsc.residual,&petsc.right),"VecDuplicate right distributed");
	NativeTetAlePetscCheck(VecDuplicate(petsc.residual,&petsc.update),"VecDuplicate update distributed");
	NativeTetAlePetscCheck(VecDuplicate(petsc.residual,&petsc.base_state),"VecDuplicate base distributed");
	// Initial owned state from the replicated trial state.
	{
		PetscScalar* owned=nullptr;NativeTetAlePetscCheck(VecGetArray(petsc.state,&owned),"VecGetArray state distributed");
		for(std::size_t vertex=0;vertex<dofs.owned_vertices;++vertex){
			const PetscInt velocity=dofs.vertex_velocity_local[vertex],pressure=dofs.vertex_pressure_local[vertex];
			for(int component=0;component<3;++component)
				owned[velocity+component]=initial_trial_state[vertex_replicated_velocity[vertex]+component];
			owned[pressure]=initial_trial_state[vertex_replicated_pressure[vertex]];
		}
		for(std::size_t edge=0;edge<local_edges;++edge)if(dofs.edge_owned[edge])
			for(int component=0;component<3;++component)
				owned[dofs.edge_velocity_local[edge]+component]=initial_trial_state[edge_replicated_velocity[edge]+component];
		NativeTetAlePetscCheck(VecRestoreArray(petsc.state,&owned),"VecRestoreArray state distributed");
	}
	// Constraints on owned rows: prescribed velocity (replicated node ids) and the gauge.
	std::vector<PetscInt> constraint_rows;std::vector<PetscScalar> constraint_values;
	for(const auto& prescribed:prescribed_velocity){
		PetscInt base=-1;
		if(prescribed.first<global_vertices){
			if(partition.node_owner[prescribed.first]!=rank)continue;
			base=dofs.vertex_velocity_row[distributed.global_to_local_vertex[prescribed.first]];
		}else{
			const std::size_t replicated_edge=prescribed.first-global_vertices;
			if(replicated_edge>=replicated_topology.edges.size())
				throw std::invalid_argument("native ALE distributed runtime: prescribed node is out of range");
			const auto& ends=replicated_topology.edges[replicated_edge];
			const std::array<std::uint32_t,2> key{{std::min(ends[0],ends[1]),std::max(ends[0],ends[1])}};
			const auto local=distributed.local_edge_of.find(key);
			if(local==distributed.local_edge_of.end()||!dofs.edge_owned[local->second])continue;
			base=dofs.edge_velocity_row[local->second];
		}
		for(int component=0;component<3;++component){
			constraint_rows.push_back(base+component);constraint_values.push_back(prescribed.second[component]);
		}
	}
	if(gauge_enabled&&partition.node_owner[pressure_gauge_node]==rank){
		constraint_rows.push_back(dofs.vertex_pressure_row[distributed.global_to_local_vertex[pressure_gauge_node]]);
		constraint_values.push_back(0.0);
	}
	// Owned boundary faces with their adjacent local cell and opposite vertex.
	struct OwnedFace{std::size_t triangle,cell;int opposite;};
	std::vector<OwnedFace> owned_faces;
	{
		std::map<std::array<std::uint32_t,3>,std::pair<std::size_t,int>> face_cell;
		for(std::size_t cell=0;cell<current.cells.size();++cell)
			for(int opposite=0;opposite<4;++opposite){
				std::array<std::uint32_t,3> face{};int entry=0;
				for(int vertex=0;vertex<4;++vertex)if(vertex!=opposite)face[entry++]=current.cells[cell].nodes[vertex];
				std::sort(face.begin(),face.end());face_cell.emplace(face,std::make_pair(cell,opposite));
			}
		for(std::size_t triangle=0;triangle<current.boundary_triangles.size();++triangle){
			if(!partition.local_triangle_owned[triangle])continue;
			auto face=current.boundary_triangles[triangle].nodes;std::sort(face.begin(),face.end());
			const auto found=face_cell.find(face);
			if(found==face_cell.end())throw std::runtime_error("native ALE distributed runtime: boundary face has no local cell");
			owned_faces.push_back({triangle,found->second.first,found->second.second});
		}
	}
	const std::array<std::array<double,3>,3> face_quadrature{{{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},
		{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
	// Solver: direct by default; fieldsplit through options with owned-row index sets.
	NativeTetAlePetscCheck(KSPCreate(comm,&petsc.solver),"KSPCreate distributed");
	KSPSetType(petsc.solver,KSPPREONLY);PC pc=nullptr;KSPGetPC(petsc.solver,&pc);PCSetType(pc,PCLU);
	PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);KSPSetFromOptions(petsc.solver);
	PCType pc_type=nullptr;NativeTetAlePetscCheck(PCGetType(pc,&pc_type),"PCGetType distributed");
	if(pc_type&&std::string(pc_type)==PCFIELDSPLIT){
		PetscBool user_schur=PETSC_FALSE;
		PetscOptionsGetBool(nullptr,nullptr,"-native_ale_user_schur",&user_schur,nullptr);
		if(user_schur)throw std::invalid_argument("native ALE distributed runtime needs -native_ale_user_schur false");
		IS velocity_is=nullptr,pressure_is=nullptr;
		NativeTetAlePetscCheck(ISCreateStride(comm,dofs.velocity_rows_end-begin,begin,1,&velocity_is),"ISCreateStride velocity");
		NativeTetAlePetscCheck(ISCreateStride(comm,end-dofs.velocity_rows_end,dofs.velocity_rows_end,1,&pressure_is),"ISCreateStride pressure");
		NativeTetAlePetscCheck(PCFieldSplitSetIS(pc,"velocity",velocity_is),"PCFieldSplitSetIS velocity distributed");
		NativeTetAlePetscCheck(PCFieldSplitSetIS(pc,"pressure",pressure_is),"PCFieldSplitSetIS pressure distributed");
		NativeTetAlePetscCheck(PCSetUseAmat(pc,PETSC_TRUE),"PCSetUseAmat distributed");
		ISDestroy(&velocity_is);ISDestroy(&pressure_is);
	}
	NativeTetAleDistributedResult result;
	KSPType ksp_type=nullptr;KSPGetType(petsc.solver,&ksp_type);
	result.linear_solver_type=ksp_type?ksp_type:"unknown";result.preconditioner_type=pc_type?pc_type:"unknown";
	result.owned_rows=local_rows;result.global_rows=global_rows;
	bool line_search_pending=false;double base_residual_norm=0.0,step_length=1.0;
	for(std::size_t iteration=0;iteration<=maximum_iterations;++iteration){
		NativeTetAlePetscCheck(VecGhostUpdateBegin(petsc.state,INSERT_VALUES,SCATTER_FORWARD),"VecGhostUpdateBegin");
		NativeTetAlePetscCheck(VecGhostUpdateEnd(petsc.state,INSERT_VALUES,SCATTER_FORWARD),"VecGhostUpdateEnd");
		Vec local_form=nullptr;NativeTetAlePetscCheck(VecGhostGetLocalForm(petsc.state,&local_form),"VecGhostGetLocalForm");
		const PetscScalar* local_state=nullptr;NativeTetAlePetscCheck(VecGetArrayRead(local_form,&local_state),"VecGetArrayRead local form");
		NativeTetAlePetscCheck(MatZeroEntries(petsc.matrix),"MatZeroEntries distributed");
		NativeTetAlePetscCheck(VecSet(petsc.residual,0.0),"VecSet residual distributed");
		for(std::size_t cell=0;cell<current.cells.size();++cell){
			if(!partition.local_cell_owned[cell])continue;
			const auto rows=dofs.CellRows(cell);const auto local=dofs.CellLocal(cell);const auto replicated=cell_replicated(cell);
			std::array<double,NativeTaylorHoodElementSystem::dofs> trial{},committed{};
			for(std::size_t entry=0;entry<rows.size();++entry){
				trial[entry]=PetscRealPart(local_state[local[entry]]);committed[entry]=committed_state[replicated[entry]];}
			std::array<std::array<double,3>,4> grid{};
			for(std::size_t node=0;node<4;++node)
				grid[node]=vertex_mesh_velocity_m_s[partition.local_to_global_node[current.cells[cell].nodes[node]]];
			const auto system=steady?BuildNativeTaylorHoodAleNavierStokesElement(current,current.cells[cell],trial,grid,parameters)
				:BuildNativeTaylorHoodAleTransientElement(current,current.cells[cell],trial,committed,grid,parameters,dt_s);
			NativeTetAlePetscCheck(MatSetValues(petsc.matrix,34,rows.data(),34,rows.data(),system.jacobian.data(),ADD_VALUES),
				"MatSetValues element distributed");
			NativeTetAlePetscCheck(VecSetValues(petsc.residual,34,rows.data(),system.residual.data(),ADD_VALUES),
				"VecSetValues element distributed");
		}
		for(const auto& face:owned_faces){
			const auto& triangle=current.boundary_triangles[face.triangle];
			const auto stabilization=boundary_conditions.backflow_stabilization_beta.find(triangle.boundary_label);
			const auto pressure=boundary_conditions.prescribed_pressure_pa.find(triangle.boundary_label);
			const bool backflow=stabilization!=boundary_conditions.backflow_stabilization_beta.end()&&stabilization->second!=0.0;
			if(!backflow&&pressure==boundary_conditions.prescribed_pressure_pa.end())continue;
			const auto rows=dofs.CellRows(face.cell);const auto local=dofs.CellLocal(face.cell);
			std::array<PetscInt,30> velocity_rows{};std::array<double,30> local_velocity{};
			for(std::size_t entry=0;entry<30;++entry){velocity_rows[entry]=rows[entry];local_velocity[entry]=PetscRealPart(local_state[local[entry]]);}
			if(backflow){
				std::array<std::array<double,3>,4> local_grid{};
				for(std::size_t node=0;node<4;++node)
					local_grid[node]=vertex_mesh_velocity_m_s[partition.local_to_global_node[current.cells[face.cell].nodes[node]]];
				const auto contribution=BuildNativeTetAleBackflowFaceSystem(current,dofs.topology,face.cell,triangle,
					local_velocity,local_grid,parameters.density,stabilization->second);
				NativeTetAlePetscCheck(VecSetValues(petsc.residual,30,velocity_rows.data(),contribution.residual.data(),ADD_VALUES),
					"VecSetValues backflow distributed");
				NativeTetAlePetscCheck(MatSetValues(petsc.matrix,30,velocity_rows.data(),30,velocity_rows.data(),
					contribution.jacobian.data(),ADD_VALUES),"MatSetValues backflow distributed");
			}
			if(pressure!=boundary_conditions.prescribed_pressure_pa.end()){
				// Outward area vector and the P2 face quadrature of the serial operator.
				const auto& p0=current.points[triangle.nodes[0]];const auto& p1=current.points[triangle.nodes[1]];
				const auto& p2=current.points[triangle.nodes[2]];
				std::array<double,3> first{},second{},area_vector{};
				for(int component=0;component<3;++component){first[component]=p1[component]-p0[component];second[component]=p2[component]-p0[component];}
				area_vector={{0.5*(first[1]*second[2]-first[2]*second[1]),0.5*(first[2]*second[0]-first[0]*second[2]),
					0.5*(first[0]*second[1]-first[1]*second[0])}};
				const auto& interior=current.points[current.cells[face.cell].nodes[face.opposite]];
				double inward=0.0;
				for(int component=0;component<3;++component)inward+=area_vector[component]*(interior[component]-p0[component]);
				if(inward>0.0)for(double& component:area_vector)component=-component;
				std::array<std::size_t,3> local_vertex{};
				for(std::size_t vertex=0;vertex<3;++vertex){
					local_vertex[vertex]=4;
					for(std::size_t node=0;node<4;++node)if(current.cells[face.cell].nodes[node]==triangle.nodes[vertex])local_vertex[vertex]=node;
					if(local_vertex[vertex]>=4)throw std::runtime_error("native ALE distributed runtime: outlet vertex is absent from its cell");
				}
				std::array<double,30> load{};
				for(const auto& face_lambda:face_quadrature){
					std::array<double,4> lambda{{0,0,0,0}};
					for(std::size_t vertex=0;vertex<3;++vertex)lambda[local_vertex[vertex]]=face_lambda[vertex];
					const auto basis=EvaluateNativeTaylorHoodBasis(lambda[1],lambda[2],lambda[3]);
					for(std::size_t node=0;node<10;++node)for(int component=0;component<3;++component)
						load[3*node+component]+=pressure->second*basis.velocity[node]*area_vector[component]/3.0;
				}
				NativeTetAlePetscCheck(VecSetValues(petsc.residual,30,velocity_rows.data(),load.data(),ADD_VALUES),
					"VecSetValues outlet pressure distributed");
			}
		}
		NativeTetAlePetscCheck(VecRestoreArrayRead(local_form,&local_state),"VecRestoreArrayRead local form");
		NativeTetAlePetscCheck(VecGhostRestoreLocalForm(petsc.state,&local_form),"VecGhostRestoreLocalForm");
		NativeTetAlePetscCheck(MatAssemblyBegin(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyBegin distributed");
		NativeTetAlePetscCheck(MatAssemblyEnd(petsc.matrix,MAT_FINAL_ASSEMBLY),"MatAssemblyEnd distributed");
		NativeTetAlePetscCheck(VecAssemblyBegin(petsc.residual),"VecAssemblyBegin distributed");
		NativeTetAlePetscCheck(VecAssemblyEnd(petsc.residual),"VecAssemblyEnd distributed");
		NativeTetAlePetscCheck(MatZeroRowsColumns(petsc.matrix,static_cast<PetscInt>(constraint_rows.size()),
			constraint_rows.data(),1.0,nullptr,nullptr),"MatZeroRowsColumns distributed");
		std::vector<PetscScalar> values(constraint_rows.size());
		NativeTetAlePetscCheck(VecGetValues(petsc.state,static_cast<PetscInt>(constraint_rows.size()),constraint_rows.data(),values.data()),
			"VecGetValues constraints distributed");
		for(std::size_t entry=0;entry<values.size();++entry)values[entry]-=constraint_values[entry];
		NativeTetAlePetscCheck(VecSetValues(petsc.residual,static_cast<PetscInt>(constraint_rows.size()),constraint_rows.data(),
			values.data(),INSERT_VALUES),"VecSetValues constraints distributed");
		NativeTetAlePetscCheck(VecAssemblyBegin(petsc.residual),"VecAssemblyBegin constraints");
		NativeTetAlePetscCheck(VecAssemblyEnd(petsc.residual),"VecAssemblyEnd constraints");
		PetscReal norm=0.0;NativeTetAlePetscCheck(VecNorm(petsc.residual,NORM_2,&norm),"VecNorm distributed");
		if(line_search_pending){
			if(std::isfinite(static_cast<double>(norm))&&norm<base_residual_norm)line_search_pending=false;
			else{
				step_length*=0.5;
				if(step_length<1.0/16384.0){
					std::ostringstream diagnostic;
					diagnostic<<std::setprecision(17)<<"native ALE distributed Newton line search failed: iteration="
						<<iteration<<" base_residual="<<base_residual_norm<<" trial_residual="<<norm
						<<" rejected_step="<<step_length;
					throw std::runtime_error(diagnostic.str());
				}
				NativeTetAlePetscCheck(VecCopy(petsc.base_state,petsc.state),"VecCopy rollback distributed");
				NativeTetAlePetscCheck(VecAXPY(petsc.state,step_length,petsc.update),"VecAXPY backtrack distributed");
				--iteration;continue;
			}
		}
		result.final_residual_l2=norm;result.newton_iterations=iteration;
		if(norm<=absolute_tolerance)break;
		if(iteration==maximum_iterations)
			throw std::runtime_error("native ALE distributed Newton did not converge: residual="+std::to_string(static_cast<double>(norm)));
		NativeTetAlePetscCheck(VecCopy(petsc.residual,petsc.right),"VecCopy right distributed");
		NativeTetAlePetscCheck(VecScale(petsc.right,-1.0),"VecScale right distributed");
		NativeTetAlePetscCheck(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators distributed");
		NativeTetAlePetscCheck(KSPSolve(petsc.solver,petsc.right,petsc.update),"KSPSolve distributed");
		KSPConvergedReason reason=KSP_CONVERGED_ITERATING;KSPGetConvergedReason(petsc.solver,&reason);
		PetscInt iterations=0;KSPGetIterationNumber(petsc.solver,&iterations);result.total_linear_iterations+=iterations;
		if(reason<=0)throw std::runtime_error("native ALE distributed linear solve did not converge: reason="
			+std::to_string(static_cast<int>(reason))+" iterations="+std::to_string(static_cast<long long>(iterations)));
		NativeTetAlePetscCheck(VecCopy(petsc.state,petsc.base_state),"VecCopy base distributed");
		base_residual_norm=norm;step_length=1.0;line_search_pending=true;
		NativeTetAlePetscCheck(VecAXPY(petsc.state,step_length,petsc.update),"VecAXPY update distributed");
	}
	// Replicated result: owned values scattered into the serial layout and summed.
	std::vector<double> local_layout(replicated_dofs,0.0);
	{
		const PetscScalar* owned=nullptr;NativeTetAlePetscCheck(VecGetArrayRead(petsc.state,&owned),"VecGetArrayRead result distributed");
		for(std::size_t vertex=0;vertex<dofs.owned_vertices;++vertex){
			for(int component=0;component<3;++component)
				local_layout[vertex_replicated_velocity[vertex]+component]=PetscRealPart(owned[dofs.vertex_velocity_local[vertex]+component]);
			local_layout[vertex_replicated_pressure[vertex]]=PetscRealPart(owned[dofs.vertex_pressure_local[vertex]]);
		}
		for(std::size_t edge=0;edge<local_edges;++edge)if(dofs.edge_owned[edge])
			for(int component=0;component<3;++component)
				local_layout[edge_replicated_velocity[edge]+component]=PetscRealPart(owned[dofs.edge_velocity_local[edge]+component]);
		NativeTetAlePetscCheck(VecRestoreArrayRead(petsc.state,&owned),"VecRestoreArrayRead result distributed");
	}
	result.replicated_state.assign(replicated_dofs,0.0);
	if(MPI_Allreduce(local_layout.data(),result.replicated_state.data(),static_cast<int>(replicated_dofs),MPI_DOUBLE,MPI_SUM,comm)
		!=MPI_SUCCESS)throw std::runtime_error("native ALE distributed runtime: result gather failed");
	return result;
}

} // namespace iga

#endif

#ifndef IGA_PERFUSED_ORGAN_HEMODYNAMICS_HPP
#define IGA_PERFUSED_ORGAN_HEMODYNAMICS_HPP

#include "NativeTetDarcyPetsc.hpp"
#include "OneDNetwork.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

// Steady organ hemodynamics: rigid Poiseuille vascular trees exchange with one
// P1 Darcy tissue compartment only at their terminals.
//
// Each terminal k owns a territory of tissue cells i with weights w_ki summing
// to one. The unresolved microvasculature between terminal k and cell i carries
//   q_ki = G_k w_ki e_i (p_k - pbar_i)    [m^3/s, vessel -> tissue positive]
// where pbar_i is the mean of the cell's four nodal pressures and
// G_k = g_tree V_k scales the tree's specific conductance g_tree [1/(Pa s)]
// with the territory volume V_k (sum of its cell volumes). With volume
// weights w_ki = V_i/V_k each cell exchanges g_tree e_i V_i (p_k - pbar_i);
// e_i > 0 is an optional per-cell exchange multiplier (default one). Segment flows are pi r^4 dp/(8 mu L). The
// tree nodes, Darcy nodes and exchange terms form one symmetric positive
// definite system whose only Dirichlet data are the tree root pressures; the
// tissue boundary carries natural zero flux. Inflow and outflow trees obey the
// same equations; their role only labels results.
struct PerfusedOrganTree
{
	std::string id;
	bool inflow=true;
	OneDNetwork network;
	// Terminal node index -> (tissue cell index, weight), weights summing to one.
	std::map<int,std::vector<std::pair<std::size_t,double>>> territories;
};

struct PerfusedOrganFlowInput
{
	double viscosity_pa_s=0.;
	std::vector<double> root_pressure_pa;              // per tree
	std::vector<double> specific_conductance_inv_pa_s; // per tree
	std::vector<double> mobility_by_cell_m2_pa_s;      // per tissue cell
	// Optional per-cell multiplier of every terminal exchange conductance into
	// that cell (local microvascular state, e.g. sinusoidal compression).
	std::vector<double> exchange_scale_by_cell;
};

struct PerfusedOrganExchange
{
	std::size_t tree=0;
	int terminal=-1;
	std::size_t cell=0;
	double flow_m3_s=0.; // vessel -> tissue positive
};

struct PerfusedOrganFlowResult
{
	std::vector<double> tissue_pressure_pa;               // per tissue node
	std::vector<double> cell_mean_pressure_pa;            // per tissue cell
	std::vector<std::vector<double>> tree_pressure_pa;    // per tree node
	std::vector<std::vector<double>> segment_flow_m3_s;   // parent -> child positive
	std::vector<double> root_inflow_m3_s;                 // per tree, into the organ positive
	std::vector<std::map<int,double>> terminal_flow_m3_s; // per tree, vessel -> tissue positive
	std::vector<PerfusedOrganExchange> exchanges;
	std::vector<double> cell_exchange_in_m3_s,cell_exchange_out_m3_s; // both nonnegative
	NativeTetDarcyResult darcy;
	double exchange_imbalance_m3_s=0.; // net tissue source removed before flux recovery
	PetscInt linear_iterations=0;
	double residual_norm=0.;
};

inline double PerfusedOrganCellVolume(const NativeTetMesh& mesh,std::size_t cell)
{
	return EvaluateNativeTetGeometry(mesh,mesh.cells[cell]).determinant/6.;
}

inline std::vector<double> PerfusedOrganTerritoryVolumes(const NativeTetMesh& mesh,
	const PerfusedOrganTree& tree)
{
	std::vector<double> result;
	for(const auto& territory:tree.territories){
		double volume=0.;
		for(const auto& entry:territory.second)volume+=PerfusedOrganCellVolume(mesh,entry.first);
		result.push_back(volume);
	}
	return result;
}

inline void ValidatePerfusedOrganTree(const NativeTetMesh& mesh,const PerfusedOrganTree& tree)
{
	const auto& network=tree.network;
	if(network.root<0||network.nodes.empty())
		throw std::invalid_argument("perfused organ tree "+tree.id+" has no root");
	if(network.has_cycles)
		throw std::invalid_argument("perfused organ tree "+tree.id+" must be acyclic");
	if(tree.territories.empty())
		throw std::invalid_argument("perfused organ tree "+tree.id+" has no territories");
	for(const auto& territory:tree.territories){
		const int node=territory.first;
		if(node<0||node>=static_cast<int>(network.nodes.size())||node==network.root)
			throw std::invalid_argument("perfused organ tree "+tree.id+" territory node is invalid");
		if(std::find(network.outlet_nodes.begin(),network.outlet_nodes.end(),node)
			==network.outlet_nodes.end())
			throw std::invalid_argument("perfused organ tree "+tree.id
				+" territory node is not a terminal");
		if(territory.second.empty())
			throw std::invalid_argument("perfused organ tree "+tree.id+" territory is empty");
		double sum=0.;
		for(const auto& entry:territory.second){
			if(entry.first>=mesh.cells.size()||!(entry.second>0.)||!std::isfinite(entry.second))
				throw std::invalid_argument("perfused organ tree "+tree.id
					+" territory cell or weight is invalid");
			sum+=entry.second;
		}
		if(std::abs(sum-1.)>1e-9)
			throw std::invalid_argument("perfused organ tree "+tree.id
				+" territory weights do not sum to one");
	}
}

namespace perfused_organ_detail {

inline void Check(PetscErrorCode error,const char* context)
{
	if(error)throw std::runtime_error(std::string("perfused organ PETSc failure: ")+context);
}

struct Triplets
{
	std::vector<PetscInt> row,column;
	std::vector<PetscScalar> value;
	void Add(PetscInt i,PetscInt j,double v){row.push_back(i);column.push_back(j);value.push_back(v);}
};

struct FlowObjects
{
	Mat matrix=nullptr;
	Vec rhs=nullptr,solution=nullptr;
	KSP solver=nullptr;
	~FlowObjects()
	{
		if(solver)KSPDestroy(&solver);
		if(solution)VecDestroy(&solution);
		if(rhs)VecDestroy(&rhs);
		if(matrix)MatDestroy(&matrix);
	}
};

} // namespace perfused_organ_detail

// Solves on PETSC_COMM_SELF. The default solver is a sparse Cholesky
// factorization (MUMPS when available); `-organ_flow_ksp_type` and
// `-organ_flow_pc_type` override it.
inline PerfusedOrganFlowResult SolvePerfusedOrganFlow(const NativeTetMesh& mesh,
	const std::vector<PerfusedOrganTree>& trees,const PerfusedOrganFlowInput& input)
{
	using namespace perfused_organ_detail;
	const std::size_t tissue_nodes=mesh.points.size(),cells=mesh.cells.size();
	if(!(input.viscosity_pa_s>0.)||!std::isfinite(input.viscosity_pa_s))
		throw std::invalid_argument("perfused organ viscosity is invalid");
	if(trees.empty()||input.root_pressure_pa.size()!=trees.size()
		||input.specific_conductance_inv_pa_s.size()!=trees.size())
		throw std::invalid_argument("perfused organ tree inputs do not match the trees");
	if(input.mobility_by_cell_m2_pa_s.size()!=cells)
		throw std::invalid_argument("perfused organ mobility does not match the cells");
	for(const double value:input.mobility_by_cell_m2_pa_s)
		if(!(value>0.)||!std::isfinite(value))
			throw std::invalid_argument("perfused organ mobility must be positive");
	if(!input.exchange_scale_by_cell.empty()&&input.exchange_scale_by_cell.size()!=cells)
		throw std::invalid_argument("perfused organ exchange scale does not match the cells");
	for(const double value:input.exchange_scale_by_cell)
		if(!(value>0.)||!std::isfinite(value))
			throw std::invalid_argument("perfused organ exchange scale must be positive");
	const auto exchange_scale=[&](std::size_t cell){
		return input.exchange_scale_by_cell.empty()?1.:input.exchange_scale_by_cell[cell];};
	for(std::size_t t=0;t<trees.size();++t){
		if(!std::isfinite(input.root_pressure_pa[t]))
			throw std::invalid_argument("perfused organ root pressure is invalid");
		if(!(input.specific_conductance_inv_pa_s[t]>0.)
			||!std::isfinite(input.specific_conductance_inv_pa_s[t]))
			throw std::invalid_argument("perfused organ terminal conductance must be positive");
		ValidatePerfusedOrganTree(mesh,trees[t]);
	}

	std::vector<NativeTetGeometry> geometry(cells);
	std::vector<double> volume(cells);
	for(std::size_t i=0;i<cells;++i){
		geometry[i]=EvaluateNativeTetGeometry(mesh,mesh.cells[i]);
		volume[i]=geometry[i].determinant/6.;
	}
	// Unknown numbering: tissue nodes first, then every non-root tree node.
	std::vector<std::vector<PetscInt>> tree_unknown(trees.size());
	PetscInt unknowns=static_cast<PetscInt>(tissue_nodes);
	for(std::size_t t=0;t<trees.size();++t){
		tree_unknown[t].assign(trees[t].network.nodes.size(),-1);
		for(std::size_t node=0;node<trees[t].network.nodes.size();++node)
			if(static_cast<int>(node)!=trees[t].network.root)tree_unknown[t][node]=unknowns++;
	}
	std::vector<double> rhs(static_cast<std::size_t>(unknowns),0.);
	Triplets entries;
	for(std::size_t i=0;i<cells;++i){
		const auto& cell=mesh.cells[i];
		const double scale=input.mobility_by_cell_m2_pa_s[i]*volume[i];
		for(std::size_t a=0;a<4;++a)
			for(std::size_t b=0;b<4;++b){
				double dot=0.;
				for(int axis=0;axis<3;++axis)
					dot+=geometry[i].barycentric_gradients[a][axis]
						*geometry[i].barycentric_gradients[b][axis];
				entries.Add(cell.nodes[a],cell.nodes[b],scale*dot);
			}
	}
	// Adds c (x_u - x_v)^2 / 2 energy with x_v possibly fixed at a root.
	const auto add_link=[&](PetscInt u,PetscInt v,double fixed_v,double c){
		entries.Add(u,u,c);
		if(v<0){rhs[static_cast<std::size_t>(u)]+=c*fixed_v;return;}
		entries.Add(v,v,c);
		entries.Add(u,v,-c);
		entries.Add(v,u,-c);
	};
	std::vector<std::vector<double>> segment_conductance(trees.size());
	std::vector<std::vector<double>> terminal_conductance(trees.size());
	for(std::size_t t=0;t<trees.size();++t){
		const auto& network=trees[t].network;
		for(const auto& segment:network.segments){
			const double conductance=OneDPi*std::pow(segment.radius0,4.)
				/(8.*input.viscosity_pa_s*segment.length);
			segment_conductance[t].push_back(conductance);
			const PetscInt parent=tree_unknown[t][static_cast<std::size_t>(segment.parent)];
			const PetscInt child=tree_unknown[t][static_cast<std::size_t>(segment.child)];
			if(child<0)throw std::invalid_argument("perfused organ tree root is a child");
			add_link(child,parent,input.root_pressure_pa[t],conductance);
		}
		const auto territory_volume=PerfusedOrganTerritoryVolumes(mesh,trees[t]);
		std::size_t index=0;
		for(const auto& territory:trees[t].territories){
			const double g=input.specific_conductance_inv_pa_s[t]*territory_volume[index++];
			terminal_conductance[t].push_back(g);
			const PetscInt terminal=tree_unknown[t][static_cast<std::size_t>(territory.first)];
			for(const auto& entry:territory.second){
				// Energy c (p_k - mean_n p_n)^2 / 2 over the cell's four nodes.
				const double c=g*entry.second*exchange_scale(entry.first);
				const auto& nodes=mesh.cells[entry.first].nodes;
				entries.Add(terminal,terminal,c);
				for(std::size_t a=0;a<4;++a){
					entries.Add(terminal,nodes[a],-0.25*c);
					entries.Add(nodes[a],terminal,-0.25*c);
					for(std::size_t b=0;b<4;++b)entries.Add(nodes[a],nodes[b],0.0625*c);
				}
			}
		}
	}

	// Triplets -> CSR with duplicates summed.
	std::vector<std::size_t> order(entries.row.size());
	std::iota(order.begin(),order.end(),0);
	std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){
		return entries.row[a]<entries.row[b]
			||(entries.row[a]==entries.row[b]&&entries.column[a]<entries.column[b]);});
	std::vector<PetscInt> row_start(static_cast<std::size_t>(unknowns)+1,0),columns;
	std::vector<PetscScalar> values;
	columns.reserve(order.size()/2);
	values.reserve(order.size()/2);
	for(std::size_t k=0;k<order.size();){
		const PetscInt r=entries.row[order[k]],c=entries.column[order[k]];
		double sum=0.;
		while(k<order.size()&&entries.row[order[k]]==r&&entries.column[order[k]]==c)
			sum+=PetscRealPart(entries.value[order[k++]]);
		columns.push_back(c);
		values.push_back(sum);
		++row_start[static_cast<std::size_t>(r)+1];
	}
	for(std::size_t r=0;r<static_cast<std::size_t>(unknowns);++r)
		row_start[r+1]+=row_start[r];
	for(std::size_t r=0;r<static_cast<std::size_t>(unknowns);++r)
		if(row_start[r+1]==row_start[r])
			throw std::invalid_argument("perfused organ system has an unconnected unknown");

	FlowObjects petsc;
	Check(MatCreate(PETSC_COMM_SELF,&petsc.matrix),"MatCreate");
	Check(MatSetSizes(petsc.matrix,unknowns,unknowns,unknowns,unknowns),"MatSetSizes");
	Check(MatSetType(petsc.matrix,MATSEQAIJ),"MatSetType");
	Check(MatSeqAIJSetPreallocationCSR(petsc.matrix,row_start.data(),columns.data(),
		values.data()),"MatSeqAIJSetPreallocationCSR");
	Check(MatSetOption(petsc.matrix,MAT_SYMMETRIC,PETSC_TRUE),"MatSetOption symmetric");
	Check(MatSetOption(petsc.matrix,MAT_SPD,PETSC_TRUE),"MatSetOption spd");
	Check(VecCreateSeq(PETSC_COMM_SELF,unknowns,&petsc.rhs),"VecCreateSeq");
	Check(VecDuplicate(petsc.rhs,&petsc.solution),"VecDuplicate");
	{
		PetscScalar* data=nullptr;
		Check(VecGetArray(petsc.rhs,&data),"VecGetArray rhs");
		for(std::size_t r=0;r<rhs.size();++r)data[r]=rhs[r];
		Check(VecRestoreArray(petsc.rhs,&data),"VecRestoreArray rhs");
	}
	Check(KSPCreate(PETSC_COMM_SELF,&petsc.solver),"KSPCreate");
	Check(KSPSetOptionsPrefix(petsc.solver,"organ_flow_"),"KSPSetOptionsPrefix");
	Check(KSPSetOperators(petsc.solver,petsc.matrix,petsc.matrix),"KSPSetOperators");
	PC preconditioner=nullptr;
	Check(KSPGetPC(petsc.solver,&preconditioner),"KSPGetPC");
	Check(KSPSetType(petsc.solver,KSPPREONLY),"KSPSetType");
	Check(PCSetType(preconditioner,PCCHOLESKY),"PCSetType");
#if defined(PETSC_HAVE_MUMPS)
	Check(PCFactorSetMatSolverType(preconditioner,MATSOLVERMUMPS),"PCFactorSetMatSolverType");
#endif
	Check(KSPSetTolerances(petsc.solver,1e-12,1e-30,PETSC_DEFAULT,2000),"KSPSetTolerances");
	Check(KSPSetFromOptions(petsc.solver),"KSPSetFromOptions");
	Check(KSPSolve(petsc.solver,petsc.rhs,petsc.solution),"KSPSolve");
	PerfusedOrganFlowResult result;
	KSPConvergedReason reason=KSP_CONVERGED_ITERATING;
	Check(KSPGetConvergedReason(petsc.solver,&reason),"KSPGetConvergedReason");
	Check(KSPGetIterationNumber(petsc.solver,&result.linear_iterations),"KSPGetIterationNumber");
	if(reason<0)throw std::runtime_error("perfused organ flow solve did not converge");
	std::vector<double> solution(static_cast<std::size_t>(unknowns));
	{
		const PetscScalar* data=nullptr;
		Check(VecGetArrayRead(petsc.solution,&data),"VecGetArrayRead");
		for(std::size_t r=0;r<solution.size();++r)solution[r]=PetscRealPart(data[r]);
		Check(VecRestoreArrayRead(petsc.solution,&data),"VecRestoreArrayRead");
	}
	{
		Vec residual=nullptr;
		Check(VecDuplicate(petsc.rhs,&residual),"VecDuplicate residual");
		Check(MatMult(petsc.matrix,petsc.solution,residual),"MatMult");
		Check(VecAXPY(residual,-1.,petsc.rhs),"VecAXPY");
		PetscReal norm=0.;
		Check(VecNorm(residual,NORM_2,&norm),"VecNorm");
		result.residual_norm=static_cast<double>(norm);
		VecDestroy(&residual);
	}
	for(const double value:solution)
		if(!std::isfinite(value))throw std::runtime_error("perfused organ flow solution is not finite");

	result.tissue_pressure_pa.assign(solution.begin(),
		solution.begin()+static_cast<std::ptrdiff_t>(tissue_nodes));
	result.cell_mean_pressure_pa.resize(cells);
	for(std::size_t i=0;i<cells;++i){
		double sum=0.;
		for(const auto node:mesh.cells[i].nodes)sum+=result.tissue_pressure_pa[node];
		result.cell_mean_pressure_pa[i]=0.25*sum;
	}
	result.cell_exchange_in_m3_s.assign(cells,0.);
	result.cell_exchange_out_m3_s.assign(cells,0.);
	std::vector<double> source(cells,0.);
	for(std::size_t t=0;t<trees.size();++t){
		const auto& network=trees[t].network;
		std::vector<double> pressure(network.nodes.size());
		for(std::size_t node=0;node<pressure.size();++node){
			const PetscInt unknown=tree_unknown[t][node];
			pressure[node]=unknown<0?input.root_pressure_pa[t]
				:solution[static_cast<std::size_t>(unknown)];
		}
		result.tree_pressure_pa.push_back(pressure);
		std::vector<double> flows;
		double root_flow=0.;
		for(std::size_t s=0;s<network.segments.size();++s){
			const auto& segment=network.segments[s];
			const double flow=segment_conductance[t][s]
				*(pressure[static_cast<std::size_t>(segment.parent)]
					-pressure[static_cast<std::size_t>(segment.child)]);
			flows.push_back(flow);
			if(segment.parent==network.root)root_flow+=flow;
		}
		result.segment_flow_m3_s.push_back(flows);
		result.root_inflow_m3_s.push_back(root_flow);
		std::map<int,double> terminals;
		std::size_t index=0;
		for(const auto& territory:trees[t].territories){
			const double g=terminal_conductance[t][index++];
			double total=0.;
			for(const auto& entry:territory.second){
				const double flow=g*entry.second*exchange_scale(entry.first)
					*(pressure[static_cast<std::size_t>(territory.first)]
						-result.cell_mean_pressure_pa[entry.first]);
				total+=flow;
				source[entry.first]+=flow;
				if(flow>=0.)result.cell_exchange_in_m3_s[entry.first]+=flow;
				else result.cell_exchange_out_m3_s[entry.first]-=flow;
				result.exchanges.push_back({t,territory.first,entry.first,flow});
			}
			terminals[territory.first]=total;
		}
		result.terminal_flow_m3_s.push_back(terminals);
	}
	// Discrete conservation holds to the linear-solve residual; remove that
	// tiny net source in proportion to |source| so the closed tissue's face-flow
	// recovery is compatible, and report the amount removed.
	double net=0.,magnitude=0.;
	for(const double value:source){net+=value;magnitude+=std::abs(value);}
	result.exchange_imbalance_m3_s=net;
	std::vector<double> source_s_inv(cells,0.);
	for(std::size_t i=0;i<cells;++i){
		const double corrected=magnitude>0.?source[i]-net*std::abs(source[i])/magnitude:0.;
		source_s_inv[i]=corrected/volume[i];
	}
	result.darcy.pressure_pa=result.tissue_pressure_pa;
	CompleteNativeTetDarcyFlux(mesh,input.mobility_by_cell_m2_pa_s,source_s_inv,{},{},result.darcy);
	return result;
}

} // namespace iga

#endif

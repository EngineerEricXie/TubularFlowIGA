#ifndef IGA_NATIVE_TET_WALL_RESERVOIR_EXCHANGE_HPP
#define IGA_NATIVE_TET_WALL_RESERVOIR_EXCHANGE_HPP

#include "NativeTetMovingSpeciesPetscRuntime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

// A fixed-volume, well-mixed species store. No hydraulic volume flow or
// pressure is implied by diffusive wall exchange.
struct NativeTetWallReservoirModel
{
	int boundary_label=-1;
	double transfer_coefficient_m_s=0.;
	double volume_m3=0.;
};

struct NativeTetWallReservoirState
{
	double amount_mol=0.;
};

struct NativeTetWallReservoirResult
{
	NativeTetMovingSpeciesPetscResult vessel;
	NativeTetWallReservoirState reservoir;
	double previous_reservoir_concentration_mol_m3=0.;
	double current_reservoir_concentration_mol_m3=0.;
	double outward_wall_flux_mol_s=0.;
	double reservoir_balance_defect_mol_s=0.;
	double combined_balance_defect_mol_s=0.;
};

struct NativeTetWallReservoirRegion
{
	NativeTetWallReservoirModel model;
	NativeTetWallReservoirState committed;
};

struct NativeTetWallReservoirsResult
{
	NativeTetMovingSpeciesPetscResult vessel;
	std::map<int,NativeTetWallReservoirState> reservoirs;
	std::map<int,double> outward_wall_flux_mol_s;
	std::map<int,double> reservoir_balance_defect_mol_s;
	double combined_balance_defect_mol_s=0.;
};

// Each distinct label owns exactly one fixed-volume store. The native FEM
// supplies the full cross-label flux response; the small dense system is the
// backward-Euler 0D Schur complement, not a replacement FEM discretization.
inline NativeTetWallReservoirsResult SolveNativeTetWallReservoirsStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	const std::map<int,NativeTetWallReservoirRegion>& regions,
	bool monotone=false,double first_order_decay_rate_s_inv=0.)
{
	if(regions.empty()||regions.size()>32||!(dt_s>0.)||!std::isfinite(dt_s))
		throw std::invalid_argument("native wall reservoir region count or time is invalid");
	std::vector<int> labels;
	std::map<int,NativeTetWallExchange> walls;
	for(const auto& item:regions){
		const auto& model=item.second.model;
		if(item.first!=model.boundary_label||item.first<0
			||!(model.transfer_coefficient_m_s>0.)
			||!std::isfinite(model.transfer_coefficient_m_s)
			||!(model.volume_m3>0.)||!std::isfinite(model.volume_m3)
			||!(item.second.committed.amount_mol>=0.)
			||!std::isfinite(item.second.committed.amount_mol))
			throw std::invalid_argument("native wall reservoir region is invalid");
		labels.push_back(item.first);
		walls.emplace(item.first,NativeTetWallExchange{model.transfer_coefficient_m_s,0.});
	}
	const auto solve=[&](const std::vector<double>& exterior){
		for(std::size_t index=0;index<labels.size();++index)
			walls.at(labels[index]).external_concentration_mol_m3=exterior[index];
		return SolveNativeTetMovingSpeciesPetscStep(previous_mesh,current_mesh,
			fluid_velocity_nodes_m_s,mesh_velocity_nodes_m_s,
			previous_concentration_mol_m3,inflow_concentration_mol_m3,
			diffusivity_m2_s,source_mol_m3_s,dt_s,monotone,
			first_order_decay_rate_s_inv,walls);
	};
	const std::size_t count=labels.size();
	std::vector<double> exterior(count,0.);
	const auto zero=solve(exterior);
	std::vector<double> baseline(count,0.),matrix(count*count,0.),rhs(count,0.);
	for(std::size_t row=0;row<count;++row)
		baseline[row]=zero.step.outward_wall_exchange_by_label_mol_s.at(labels[row]);
	for(std::size_t column=0;column<count;++column){
		exterior[column]=1.;
		const auto unit=solve(exterior);
		exterior[column]=0.;
		for(std::size_t row=0;row<count;++row)
			matrix[row*count+column]=-dt_s*(
				unit.step.outward_wall_exchange_by_label_mol_s.at(labels[row])
				-baseline[row]);
	}
	for(std::size_t row=0;row<count;++row){
		const auto& region=regions.at(labels[row]);
		matrix[row*count+row]+=region.model.volume_m3;
		rhs[row]=region.committed.amount_mol+dt_s*baseline[row];
	}
	exterior=native_tet_moving_species_detail::SolveDense(matrix,rhs);
	for(const double value:exterior)
		if(!(value>=0.)||!std::isfinite(value))
			throw std::runtime_error("native wall reservoir coupled concentration is invalid");
	NativeTetWallReservoirsResult result;
	result.vessel=solve(exterior);
	double previous_amount=0.,current_amount=0.;
	for(std::size_t row=0;row<count;++row){
		const int label=labels[row];
		const auto& region=regions.at(label);
		const double amount=region.model.volume_m3*exterior[row];
		const double flux=result.vessel.step.outward_wall_exchange_by_label_mol_s.at(label);
		const double defect=(amount-region.committed.amount_mol)/dt_s-flux;
		const double scale=std::max({1e-30,
			std::abs(amount/dt_s),std::abs(region.committed.amount_mol/dt_s),
			std::abs(flux)});
		if(!std::isfinite(amount)||!std::isfinite(defect)
			||std::abs(defect)>1e-10*scale+1e-14)
			throw std::runtime_error("native wall reservoir region balance failed");
		result.reservoirs.emplace(label,NativeTetWallReservoirState{amount});
		result.outward_wall_flux_mol_s.emplace(label,flux);
		result.reservoir_balance_defect_mol_s.emplace(label,defect);
		previous_amount+=region.committed.amount_mol;
		current_amount+=amount;
	}
	const auto& step=result.vessel.step;
	result.combined_balance_defect_mol_s=(step.current_inventory_mol
		-step.previous_inventory_mol+current_amount-previous_amount)/dt_s
		+step.outward_advective_flux_mol_s+step.reaction_sink_mol_s-step.source_mol_s;
	const double scale=std::max({1e-30,std::abs(step.previous_inventory_mol/dt_s),
		std::abs(previous_amount/dt_s),std::abs(step.source_mol_s),
		std::abs(step.outward_wall_exchange_mol_s)});
	if(!std::isfinite(result.combined_balance_defect_mol_s)
		||std::abs(result.combined_balance_defect_mol_s)>1e-10*scale+1e-14)
		throw std::runtime_error("native wall reservoir combined balance failed");
	return result;
}

// The FEM matrix is independent of the prescribed exterior concentration.
// Two native FEM/PETSc solves determine F(C)=F(0)+C*(F(1)-F(0)); the scalar
// backward-Euler reservoir balance then gives the monolithic linear solution.
// A final native solve evaluates the accepted field and both budgets. No
// committed state is modified until the caller accepts this return value.
inline NativeTetWallReservoirResult SolveNativeTetWallReservoirStep(
	const NativeTetMesh& previous_mesh,const NativeTetMesh& current_mesh,
	const std::vector<std::array<double,3>>& fluid_velocity_nodes_m_s,
	const std::vector<std::array<double,3>>& mesh_velocity_nodes_m_s,
	const std::vector<double>& previous_concentration_mol_m3,
	const std::map<int,double>& inflow_concentration_mol_m3,
	double diffusivity_m2_s,double source_mol_m3_s,double dt_s,
	const NativeTetWallReservoirModel& model,
	const NativeTetWallReservoirState& committed,bool monotone=false,
	double first_order_decay_rate_s_inv=0.)
{
	if(model.boundary_label<0||!(model.transfer_coefficient_m_s>0.)
		||!std::isfinite(model.transfer_coefficient_m_s)
		||!(model.volume_m3>0.)||!std::isfinite(model.volume_m3)
		||!(committed.amount_mol>=0.)||!std::isfinite(committed.amount_mol)
		||!(dt_s>0.)||!std::isfinite(dt_s))
		throw std::invalid_argument("native wall reservoir model or state is invalid");
	const auto solve=[&](double exterior){
		return SolveNativeTetMovingSpeciesPetscStep(previous_mesh,current_mesh,
			fluid_velocity_nodes_m_s,mesh_velocity_nodes_m_s,
			previous_concentration_mol_m3,inflow_concentration_mol_m3,
			diffusivity_m2_s,source_mol_m3_s,dt_s,monotone,
			first_order_decay_rate_s_inv,
			{{model.boundary_label,{model.transfer_coefficient_m_s,exterior}}});
	};
	const auto zero=solve(0.);
	const auto unit=solve(1.);
	const double flux_zero=zero.step.outward_wall_exchange_by_label_mol_s.at(
		model.boundary_label);
	const double flux_unit=unit.step.outward_wall_exchange_by_label_mol_s.at(
		model.boundary_label);
	const double slope=flux_unit-flux_zero;
	const double scale=std::max({1.,std::abs(flux_zero),std::abs(flux_unit)});
	if(!std::isfinite(slope)||slope>1e-10*scale)
		throw std::runtime_error("native wall reservoir exchange response is not passive");
	const double denominator=model.volume_m3-dt_s*slope;
	if(!(denominator>0.)||!std::isfinite(denominator))
		throw std::runtime_error("native wall reservoir coupled denominator is invalid");
	const double exterior=(committed.amount_mol+dt_s*flux_zero)/denominator;
	if(!(exterior>=0.)||!std::isfinite(exterior))
		throw std::runtime_error("native wall reservoir coupled concentration is invalid");
	NativeTetWallReservoirResult result;
	result.vessel=solve(exterior);
	result.previous_reservoir_concentration_mol_m3=
		committed.amount_mol/model.volume_m3;
	result.current_reservoir_concentration_mol_m3=exterior;
	result.reservoir.amount_mol=model.volume_m3*exterior;
	result.outward_wall_flux_mol_s=
		result.vessel.step.outward_wall_exchange_by_label_mol_s.at(
			model.boundary_label);
	result.reservoir_balance_defect_mol_s=(result.reservoir.amount_mol
		-committed.amount_mol)/dt_s-result.outward_wall_flux_mol_s;
	result.combined_balance_defect_mol_s=(result.vessel.step.current_inventory_mol
		-result.vessel.step.previous_inventory_mol
		+result.reservoir.amount_mol-committed.amount_mol)/dt_s
		+result.vessel.step.outward_advective_flux_mol_s
		+result.vessel.step.reaction_sink_mol_s-result.vessel.step.source_mol_s;
	const double budget_scale=std::max({1e-30,
		std::abs(result.vessel.step.previous_inventory_mol/dt_s),
		std::abs(committed.amount_mol/dt_s),
		std::abs(result.outward_wall_flux_mol_s),
		std::abs(result.vessel.step.source_mol_s)});
	if(!std::isfinite(result.reservoir.amount_mol)
		||!std::isfinite(result.reservoir_balance_defect_mol_s)
		||!std::isfinite(result.combined_balance_defect_mol_s)
		||std::abs(result.reservoir_balance_defect_mol_s)>1e-10*budget_scale+1e-14
		||std::abs(result.combined_balance_defect_mol_s)>1e-10*budget_scale+1e-14)
		throw std::runtime_error("native wall reservoir paired species balance failed");
	return result;
}

} // namespace iga

#endif

#ifndef IGA_ZERO_D_SPECIES_RESERVOIR_HPP
#define IGA_ZERO_D_SPECIES_RESERVOIR_HPP

#include "CoupledDomainRuntime.hpp"

#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace iga {

// A well-mixed, variable-volume 0D species compartment. Hydraulic pressure
// and flow are supplied by a separate model; this model never infers them.
struct ZeroDSpeciesReservoirModel
{
	std::set<std::string> port_ids;
	std::set<std::string> species_ids;
	std::map<std::string,double> source_rate_mol_s;
};

struct ZeroDSpeciesReservoirState
{
	double volume_m3=0.0;
	std::map<std::string,double> amount_mol;
};

struct ZeroDSpeciesReservoirPortInput
{
	double outward_flow_m3_s=0.0;
	// Required for every species only when flow enters this compartment.
	std::map<std::string,double> incoming_concentration_mol_m3;
};

struct ZeroDSpeciesReservoirBalance
{
	double initial_volume_m3=0.0;
	double final_volume_m3=0.0;
	double outward_port_volume_m3=0.0;
	double volume_residual_m3=0.0;
	std::map<std::string,SpeciesStepAccounting> species;
};

struct ZeroDSpeciesReservoirTrial
{
	ZeroDSpeciesReservoirState state;
	std::map<std::string,PortState> ports;
	ZeroDSpeciesReservoirBalance balance;
};

inline void ValidateZeroDSpeciesReservoirModel(
	const ZeroDSpeciesReservoirModel& model)
{
	if(model.port_ids.empty()||model.species_ids.empty())
		throw std::invalid_argument("0D species reservoir requires ports and species");
	for(const auto& id:model.port_ids)
		if(id.empty())throw std::invalid_argument("0D species reservoir port id is empty");
	for(const auto& id:model.species_ids)
		if(id.empty())throw std::invalid_argument("0D species reservoir species id is empty");
	for(const auto& source:model.source_rate_mol_s)
		if(!model.species_ids.count(source.first)||!std::isfinite(source.second))
			throw std::invalid_argument("0D species reservoir source is invalid");
}

inline void ValidateZeroDSpeciesReservoirState(
	const ZeroDSpeciesReservoirModel& model,
	const ZeroDSpeciesReservoirState& state)
{
	ValidateZeroDSpeciesReservoirModel(model);
	if(!(state.volume_m3>0.0)||!std::isfinite(state.volume_m3))
		throw std::invalid_argument("0D species reservoir volume must be positive and finite");
	if(state.amount_mol.size()!=model.species_ids.size())
		throw std::invalid_argument("0D species reservoir amount set differs from model");
	for(const auto& species:model.species_ids){
		const auto found=state.amount_mol.find(species);
		if(found==state.amount_mol.end()||!(found->second>=0.0)
			||!std::isfinite(found->second))
			throw std::invalid_argument("0D species reservoir amount is invalid");
	}
}

inline ZeroDSpeciesReservoirTrial EvaluateZeroDSpeciesReservoirTrial(
	const ZeroDSpeciesReservoirModel& model,
	const ZeroDSpeciesReservoirState& committed,
	const std::map<std::string,ZeroDSpeciesReservoirPortInput>& inputs,
	double dt_s,double end_time_s)
{
	ValidateZeroDSpeciesReservoirState(model,committed);
	if(!(dt_s>0.0)||!std::isfinite(dt_s)||!std::isfinite(end_time_s)
		||inputs.size()!=model.port_ids.size())
		throw std::invalid_argument("0D species reservoir step or port set is invalid");
	double outward_flow=0.0,outgoing_flow=0.0;
	for(const auto& port:model.port_ids){
		const auto found=inputs.find(port);
		if(found==inputs.end()||!std::isfinite(found->second.outward_flow_m3_s))
			throw std::invalid_argument("0D species reservoir port flow is missing or invalid");
		const auto& input=found->second;
		outward_flow+=input.outward_flow_m3_s;
		if(input.outward_flow_m3_s>0.0)
			outgoing_flow+=input.outward_flow_m3_s;
		if(input.outward_flow_m3_s<0.0){
			if(input.incoming_concentration_mol_m3.size()!=model.species_ids.size())
				throw std::invalid_argument("0D species reservoir inflow needs every concentration");
			for(const auto& species:model.species_ids){
				const auto concentration=input.incoming_concentration_mol_m3.find(species);
				if(concentration==input.incoming_concentration_mol_m3.end()
					||!(concentration->second>=0.0)||!std::isfinite(concentration->second))
					throw std::invalid_argument("0D species reservoir inflow concentration is invalid");
			}
		}else if(!input.incoming_concentration_mol_m3.empty())
			throw std::invalid_argument("0D species reservoir non-inflow has donor concentration");
	}
	ZeroDSpeciesReservoirTrial result;
	result.state.volume_m3=committed.volume_m3-dt_s*outward_flow;
	if(!(result.state.volume_m3>0.0)||!std::isfinite(result.state.volume_m3))
		throw std::runtime_error("0D species reservoir trial volume is nonpositive or nonfinite");
	result.balance.initial_volume_m3=committed.volume_m3;
	result.balance.final_volume_m3=result.state.volume_m3;
	result.balance.outward_port_volume_m3=dt_s*outward_flow;
	result.balance.volume_residual_m3=result.state.volume_m3
		-committed.volume_m3+result.balance.outward_port_volume_m3;
	for(const auto& species:model.species_ids){
		const auto source=model.source_rate_mol_s.find(species);
		const double source_rate=source==model.source_rate_mol_s.end()?0.0:source->second;
		double incoming_flux=0.0;
		for(const auto& port:model.port_ids){
			const auto& input=inputs.at(port);
			if(input.outward_flow_m3_s<0.0)
				incoming_flux+=input.outward_flow_m3_s
					*input.incoming_concentration_mol_m3.at(species);
		}
		const double amount=(committed.amount_mol.at(species)
			+dt_s*(source_rate-incoming_flux))
			/(1.0+dt_s*outgoing_flow/result.state.volume_m3);
		if(!(amount>=0.0)||!std::isfinite(amount))
			throw std::runtime_error("0D species reservoir trial amount is invalid");
		result.state.amount_mol.emplace(species,amount);
		SpeciesStepAccounting accounting;
		accounting.initial_mass=committed.amount_mol.at(species);
		accounting.final_mass=amount;
		accounting.source_amount=dt_s*source_rate;
		result.balance.species.emplace(species,std::move(accounting));
	}
	for(const auto& port:model.port_ids){
		const auto& input=inputs.at(port);
		PortState state;
		state.time_s=end_time_s;
		state.outward_flow_m3_s=input.outward_flow_m3_s;
		for(const auto& species:model.species_ids){
			const double mixed=result.state.amount_mol.at(species)/result.state.volume_m3;
			const double donor=input.outward_flow_m3_s<0.0
				?input.incoming_concentration_mol_m3.at(species):mixed;
			const double flux=input.outward_flow_m3_s*donor;
			state.concentration.emplace(species,mixed);
			state.outward_species_flux.emplace(species,flux);
			auto& accounting=result.balance.species.at(species);
			accounting.outward_port_amount.emplace(port,dt_s*flux);
		}
		result.ports.emplace(port,std::move(state));
	}
	for(auto& species:result.balance.species){
		auto& accounting=species.second;
		accounting.residual=accounting.final_mass-accounting.initial_mass
			-accounting.source_amount;
		for(const auto& port:accounting.outward_port_amount)
			accounting.residual+=port.second;
		if(!std::isfinite(accounting.residual))
			throw std::runtime_error("0D species reservoir amount residual is invalid");
	}
	return result;
}

} // namespace iga

#endif

#include "ZeroDSpeciesReservoir.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace {

bool Near(double first,double second)
{
	return std::abs(first-second)<=1e-12
		*std::max(1.0,std::max(std::abs(first),std::abs(second)));
}

void RequireRejected(const std::function<void()>& action)
{
	bool rejected=false;
	try{action();}catch(const std::exception&){rejected=true;}
	assert(rejected);
}

} // namespace

int main()
{
	iga::ZeroDSpeciesReservoirModel model;
	model.port_ids={"inlet","outlet"};
	model.species_ids={"oxygen","tracer"};
	model.source_rate_mol_s={{"tracer",0.12}};
	iga::ZeroDSpeciesReservoirState committed;
	committed.volume_m3=1.0;
	committed.amount_mol={{"oxygen",2.0},{"tracer",1.0}};
	std::map<std::string,iga::ZeroDSpeciesReservoirPortInput> inputs;
	inputs["inlet"]={-0.2,{{"oxygen",5.0},{"tracer",1.0}}};
	inputs["outlet"]={0.2,{}};
	const auto trial=iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);
	assert(Near(trial.state.volume_m3,1.0));
	assert(Near(trial.state.amount_mol.at("oxygen"),2.5));
	assert(Near(trial.state.amount_mol.at("tracer"),1.1));
	assert(Near(*trial.ports.at("inlet").outward_flow_m3_s,-0.2));
	assert(Near(trial.ports.at("inlet").outward_species_flux.at("oxygen"),-1.0));
	assert(Near(trial.ports.at("outlet").outward_species_flux.at("oxygen"),0.5));
	assert(Near(trial.ports.at("outlet").concentration.at("oxygen"),2.5));
	assert(Near(trial.ports.at("outlet").outward_species_flux.at("tracer"),0.22));
	assert(Near(trial.balance.volume_residual_m3,0.0));
	for(const auto& species:trial.balance.species)
		assert(Near(species.second.residual,0.0));
	assert(Near(trial.balance.species.at("tracer").source_amount,0.12));

	// Reverse the hydraulic direction. The new inlet needs a donor concentration;
	// the old inlet now publishes the compartment's mixed concentration.
	inputs["inlet"]={0.1,{}};
	inputs["outlet"]={-0.1,{{"oxygen",3.0},{"tracer",0.0}}};
	const auto reverse=iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,0.5,0.5);
	assert(Near(reverse.state.volume_m3,1.0));
	assert(Near(reverse.state.amount_mol.at("oxygen"),2.15/1.05));
	assert(Near(reverse.ports.at("outlet").outward_species_flux.at("oxygen"),-0.3));
	assert(Near(reverse.ports.at("inlet").outward_species_flux.at("oxygen"),
		0.1*reverse.state.amount_mol.at("oxygen")));
	for(const auto& species:reverse.balance.species)
		assert(Near(species.second.residual,0.0));

	// Pure influx expands the 0D volume without inventing or losing amount.
	inputs["inlet"]={-0.5,{{"oxygen",0.0},{"tracer",0.0}}};
	inputs["outlet"]={0.0,{}};
	const auto expanded=iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);
	assert(Near(expanded.state.volume_m3,1.5));
	assert(Near(expanded.state.amount_mol.at("oxygen"),2.0));
	assert(Near(expanded.ports.at("outlet").concentration.at("oxygen"),4.0/3.0));
	assert(Near(expanded.balance.volume_residual_m3,0.0));
	inputs["inlet"]={0.0,{}};
	inputs["outlet"]={0.0,{}};
	const auto still=iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);
	assert(Near(still.state.volume_m3,1.0));
	assert(Near(still.state.amount_mol.at("oxygen"),2.0));
	assert(Near(still.state.amount_mol.at("tracer"),1.12));
	assert(Near(still.ports.at("inlet").outward_species_flux.at("oxygen"),0.0));
	for(const auto& species:still.balance.species)
		assert(Near(species.second.residual,0.0));

	inputs["inlet"]={-0.1,{{"oxygen",1.0}}};
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);}); // Missing tracer donor.
	inputs["inlet"]={0.1,{{"oxygen",1.0},{"tracer",1.0}}};
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);}); // Outgoing port cannot silently ignore donor data.
	inputs["inlet"]={2.0,{}};
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);}); // Drained volume.
	inputs["inlet"]={std::numeric_limits<double>::quiet_NaN(),{}};
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);});
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,0.0,1.0);});
	inputs.erase("outlet");
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);});
	auto invalid_state=committed;
	invalid_state.amount_mol["oxygen"]=-1.0;
	RequireRejected([&]{iga::ValidateZeroDSpeciesReservoirState(model,invalid_state);});
	model.source_rate_mol_s["oxygen"]=-4.0;
	inputs["outlet"]={0.0,{}};
	inputs["inlet"]={0.0,{}};
	RequireRejected([&]{iga::EvaluateZeroDSpeciesReservoirTrial(model,committed,
		inputs,1.0,1.0);});
	model.source_rate_mol_s["unknown"]=1.0;
	RequireRejected([&]{iga::ValidateZeroDSpeciesReservoirModel(model);});
	return 0;
}

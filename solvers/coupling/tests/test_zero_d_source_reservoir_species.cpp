#include "ZeroDSourceReservoirSpecies.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
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
	iga::ZeroDFlowModel hydraulic;
	hydraulic.role=iga::ZeroDFlowRole::SourceReservoir;
	hydraulic.source={0.5,2.0,0.2};
	const iga::ZeroDFlowState pressure{1.0};
	iga::ZeroDSpeciesReservoirModel species;
	species.port_ids={"graph","pump"};
	species.species_ids={"a","b"};
	const iga::ZeroDSpeciesReservoirState amount{1.0,{{"a",2.0},{"b",3.0}}};
	const std::map<std::string,double> donor{{"a",5.0},{"b",7.0}};
	const auto forward=iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,0.0,{},donor,0.1,0.0);
	assert(*forward.hydraulic.port.outward_flow_m3_s>0.0);
	assert(Near(forward.transport.state.volume_m3-amount.volume_m3,
		forward.hydraulic.storage.stored_volume_change_m3));
	assert(Near(forward.transport.ports.at("pump").outward_species_flux.at("a"),-1.0));
	assert(forward.transport.ports.at("graph").outward_species_flux.at("a")>0.0);
	assert(Near(forward.transport.balance.species.at("a").residual,0.0));
	assert(Near(forward.transport.balance.volume_residual_m3,0.0));
	RequireRejected([&]{iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,0.0,{},{},0.1,0.0);});

	const auto reverse=iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,10.0,donor,donor,0.1,0.0);
	assert(*reverse.hydraulic.port.outward_flow_m3_s<0.0);
	assert(reverse.transport.ports.at("graph").outward_species_flux.at("a")<0.0);
	assert(Near(reverse.transport.balance.species.at("a").residual,0.0));
	RequireRejected([&]{iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,10.0,{},donor,0.1,0.0);});

	hydraulic.source.prescribed_flow_m3_s=-0.2;
	const auto pump_reverse=iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,0.0,{}, {},0.1,0.0);
	assert(pump_reverse.transport.ports.at("pump").outward_species_flux.at("a")>0.0);
	assert(Near(pump_reverse.transport.balance.species.at("b").residual,0.0));
	RequireRejected([&]{iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,0.0,{},donor,0.1,0.0);});
	species.port_ids={"graph"};
	RequireRejected([&]{iga::EvaluateZeroDSourceReservoirSpeciesTrial(hydraulic,
		pressure,species,amount,0.0,{}, {},0.1,0.0);});
	return 0;
}

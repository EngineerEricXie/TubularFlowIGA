#include "ZeroDTerminalRcrSpecies.hpp"

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
	hydraulic.role=iga::ZeroDFlowRole::TerminalRcr;
	hydraulic.terminal={1.5,5.0,2.0,7.0};
	const iga::ZeroDFlowState pressure{9.0};
	iga::ZeroDSpeciesReservoirModel species;
	species.port_ids={"graph","distal"};
	species.species_ids={"tracer"};
	const iga::ZeroDSpeciesReservoirState amount{1.0,{{"tracer",2.0}}};
	const auto forward=iga::EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic,
		pressure,species,amount,-2.0,{{"tracer",5.0}},{},0.1,0.0);
	const double distal_flow=forward.hydraulic.storage.distal_sink_amount_m3/0.1;
	assert(distal_flow>0.0);
	assert(Near(forward.transport.state.volume_m3-amount.volume_m3,
		forward.hydraulic.storage.stored_volume_change_m3));
	assert(Near(forward.transport.ports.at("graph").outward_species_flux.at("tracer"),
		-10.0));
	assert(Near(forward.transport.ports.at("distal").outward_species_flux.at("tracer"),
		distal_flow*forward.transport.state.amount_mol.at("tracer")
			/forward.transport.state.volume_m3));
	assert(Near(forward.transport.balance.species.at("tracer").residual,0.0));
	assert(Near(forward.transport.balance.volume_residual_m3,0.0));
	RequireRejected([&]{iga::EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic,
		pressure,species,amount,-2.0,{},{},0.1,0.0);});

	// Reverse graph flow also reverses the distal RCR branch for this load.
	hydraulic.terminal.distal_pressure_pa=20.0;
	const auto reverse=iga::EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic,
		pressure,species,amount,0.75,{},{{"tracer",3.0}},0.1,0.0);
	assert(reverse.hydraulic.storage.distal_sink_amount_m3<0.0);
	assert(reverse.transport.ports.at("graph").outward_species_flux.at("tracer")>0.0);
	assert(reverse.transport.ports.at("distal").outward_species_flux.at("tracer")<0.0);
	assert(Near(reverse.transport.state.volume_m3-amount.volume_m3,
		reverse.hydraulic.storage.stored_volume_change_m3));
	assert(Near(reverse.transport.balance.species.at("tracer").residual,0.0));
	RequireRejected([&]{iga::EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic,
		pressure,species,amount,0.75,{},{},0.1,0.0);});
	species.port_ids={"graph"};
	RequireRejected([&]{iga::EvaluateZeroDTerminalRcrSpeciesTrial(hydraulic,
		pressure,species,amount,0.75,{},{},0.1,0.0);});
	return 0;
}

#include "ZeroDSourceReservoirSpeciesDomainRuntime.hpp"

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

iga::PortBoundaryData Pressure(double time,double value)
{
	iga::PortBoundaryData result;
	result.time_s=time;result.mean_pressure_pa=value;
	return result;
}

} // namespace

int main()
{
	iga::ZeroDFlowModel flow_model;
	flow_model.role=iga::ZeroDFlowRole::SourceReservoir;
	flow_model.source={0.5,2.0,0.2};
	iga::ZeroDSpeciesReservoirModel species_model;
	species_model.port_ids={"graph","pump"};
	species_model.species_ids={"tracer"};
	RequireRejected([&]{iga::ZeroDSourceReservoirSpeciesDomainRuntime missing(
		"missing",flow_model,{1.0},species_model,{1.0,{{"tracer",2.0}}});});
	iga::ZeroDSourceReservoirSpeciesDomainRuntime runtime("source",flow_model,{1.0},
		species_model,{1.0,{{"tracer",2.0}}},{{"tracer",5.0}});
	assert(runtime.Ports().size()==1);
	assert(runtime.Ports().front().species==std::set<std::string>{"tracer"});
	runtime.BeginStep({0,0.0,0.1});
	runtime.SetPortInput("port",Pressure(0.1,0.0));
	runtime.SolveHydraulicTrial();
	assert(*runtime.GetHydraulicPortState("port").outward_flow_m3_s>0.0);
	const auto before=runtime.GetTransportPortState("port");
	assert(before.outward_species_flux.at("tracer")>0.0);
	RequireRejected([&]{runtime.SetTransportConcentration("port",0.1,{{"tracer",3.0}});});
	runtime.SolveTransportTrial();
	const auto accounting=runtime.GetSpeciesStepAccounting().at("tracer");
	assert(Near(accounting.residual,0.0));
	assert(accounting.outward_port_amount.at("pump")<0.0);
	assert(accounting.outward_port_amount.at("port")>0.0);
	runtime.RollbackTransportTrial();
	RequireRejected([&]{runtime.GetSpeciesStepAccounting();});
	runtime.RollbackHydraulicTrial();
	runtime.SetPortInput("port",Pressure(0.1,0.0));
	runtime.SolveTrial();
	runtime.PrepareCommitStep();
	assert(Near(runtime.CommittedSpeciesState().amount_mol.at("tracer"),2.0));
	runtime.FinalizeCommitStep();
	assert(runtime.CommittedStepIndex()==0);
	const double committed_amount=runtime.CommittedSpeciesState().amount_mol.at("tracer");
	const double committed_pressure=runtime.CommittedHydraulicState().stored_pressure_pa;
	runtime.BeginStep({1,0.1,0.1});
	runtime.SetPortInput("port",Pressure(0.2,10.0));
	runtime.SolveHydraulicTrial();
	assert(*runtime.GetHydraulicPortState("port").outward_flow_m3_s<0.0);
	RequireRejected([&]{runtime.GetTransportPortState("port");});
	RequireRejected([&]{runtime.SolveTransportTrial();});
	RequireRejected([&]{runtime.SetTransportConcentration("port",0.2,{{"wrong",3.0}});});
	runtime.SetTransportConcentration("port",0.2,{{"tracer",3.0}});
	runtime.SolveTransportTrial();
	assert(runtime.GetTransportPortState("port").outward_species_flux.at("tracer")<0.0);
	runtime.AbortStep();
	assert(runtime.CommittedStepIndex()==0);
	assert(Near(runtime.CommittedSpeciesState().amount_mol.at("tracer"),committed_amount));
	assert(Near(runtime.CommittedHydraulicState().stored_pressure_pa,committed_pressure));
	return 0;
}

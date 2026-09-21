#include "ZeroDTerminalRcrSpeciesDomainRuntime.hpp"

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

iga::PortBoundaryData Flow(double time,double value)
{
	iga::PortBoundaryData result;
	result.time_s=time;result.outward_flow_m3_s=value;
	return result;
}

} // namespace

int main()
{
	iga::ZeroDFlowModel flow_model;
	flow_model.role=iga::ZeroDFlowRole::TerminalRcr;
	flow_model.terminal={1.5,5.0,2.0,7.0};
	iga::ZeroDSpeciesReservoirModel species_model;
	species_model.port_ids={"graph","distal"};
	species_model.species_ids={"tracer"};
	iga::ZeroDTerminalRcrSpeciesDomainRuntime runtime("terminal",flow_model,{9.0},
		species_model,{1.0,{{"tracer",2.0}}},{{"tracer",3.0}});
	assert(runtime.Ports().size()==1);
	assert(runtime.Ports().front().species==std::set<std::string>{"tracer"});
	RequireRejected([&]{runtime.SolveHydraulicTrial();});
	runtime.BeginStep({0,0.0,0.1});
	runtime.SetPortInput("port",Flow(0.1,-2.0));
	runtime.SolveHydraulicTrial();
	assert(Near(*runtime.GetHydraulicPortState("port").outward_flow_m3_s,-2.0));
	RequireRejected([&]{runtime.GetTransportPortState("port");});
	RequireRejected([&]{runtime.SolveTransportTrial();});
	RequireRejected([&]{runtime.SetTransportConcentration("port",0.1,{{"wrong",5.0}});});
	runtime.SetTransportConcentration("port",0.1,{{"tracer",5.0}});
	runtime.SolveTransportTrial();
	const auto first=runtime.GetTransportPortState("port");
	assert(Near(first.outward_species_flux.at("tracer"),-10.0));
	assert(first.mean_pressure_pa.has_value());
	assert(Near(runtime.GetSpeciesStepAccounting().at("tracer").residual,0.0));
	assert(runtime.GetSpeciesStepAccounting().at("tracer").outward_port_amount.count("port"));
	runtime.RollbackTransportTrial();
	RequireRejected([&]{runtime.GetSpeciesStepAccounting();});
	runtime.SetTransportConcentration("port",0.1,{{"tracer",5.0}});
	runtime.SolveTransportTrial();
	runtime.PrepareCommitStep();
	assert(Near(runtime.CommittedSpeciesState().amount_mol.at("tracer"),2.0));
	runtime.FinalizeCommitStep();
	assert(runtime.CommittedStepIndex()==0);
	assert(runtime.CommittedSpeciesState().amount_mol.at("tracer")>2.0);
	const double committed_amount=runtime.CommittedSpeciesState().amount_mol.at("tracer");
	const double committed_pressure=runtime.CommittedHydraulicState().stored_pressure_pa;

	runtime.BeginStep({1,0.1,0.1});
	runtime.SetPortInput("port",Flow(0.2,0.75));
	runtime.SolveHydraulicTrial();
	const auto donor=runtime.GetTransportPortState("port");
	assert(donor.outward_species_flux.at("tracer")>0.0);
	RequireRejected([&]{runtime.SetTransportConcentration("port",0.2,{{"tracer",5.0}});});
	runtime.RollbackHydraulicTrial();
	assert(Near(runtime.CommittedSpeciesState().amount_mol.at("tracer"),committed_amount));
	runtime.SetPortInput("port",Flow(0.2,0.75));
	runtime.SolveHydraulicTrial();
	runtime.SolveTransportTrial();
	assert(Near(runtime.GetSpeciesStepAccounting().at("tracer").residual,0.0));
	runtime.AbortStep();
	assert(runtime.CommittedStepIndex()==0);
	assert(Near(runtime.CommittedSpeciesState().amount_mol.at("tracer"),committed_amount));
	assert(Near(runtime.CommittedHydraulicState().stored_pressure_pa,committed_pressure));

	flow_model.terminal.distal_pressure_pa=20.0;
	iga::ZeroDTerminalRcrSpeciesDomainRuntime reverse("reverse",flow_model,{9.0},
		species_model,{1.0,{{"tracer",2.0}}},{{"tracer",3.0}});
	reverse.BeginStep({0,0.0,0.1});
	reverse.SetPortInput("port",Flow(0.1,0.75));
	reverse.SolveHydraulicTrial();
	assert(reverse.GetTransportPortState("port").outward_species_flux.at("tracer")>0.0);
	reverse.SolveTransportTrial();
	assert(reverse.GetSpeciesStepAccounting().at("tracer").outward_port_amount.at("distal")<0.0);
	reverse.PrepareCommitStep();reverse.FinalizeCommitStep();
	iga::ZeroDTerminalRcrSpeciesDomainRuntime missing_distal("missing",flow_model,{9.0},
		species_model,{1.0,{{"tracer",2.0}}});
	missing_distal.BeginStep({0,0.0,0.1});
	missing_distal.SetPortInput("port",Flow(0.1,0.75));
	missing_distal.SolveHydraulicTrial();
	RequireRejected([&]{missing_distal.SolveTransportTrial();});
	missing_distal.AbortStep();
	assert(missing_distal.CommittedStepIndex()==-1);
	assert(Near(missing_distal.CommittedSpeciesState().amount_mol.at("tracer"),2.0));
	return 0;
}

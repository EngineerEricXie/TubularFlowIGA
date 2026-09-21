#include "SpeciesPressureFlowComponentExecutor.hpp"
#include "ZeroDTerminalRcrSpeciesDomainRuntime.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>

namespace {

class ControlledSource final : public iga::CoupledDomainRuntime,
	public iga::StagedFlowTransportDomainRuntime
{
public:
	ControlledSource()
	{
		iga::CouplingPort port;
		port.id="port";port.subsystem_id="source";
		port.locator_kind="controlled_source";port.locator="port";
		port.provides={iga::PortQuantity::FlowRate,iga::PortQuantity::MeanPressure,
			iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
		port.requires={iga::PortQuantity::MeanPressure,
			iga::PortQuantity::SpeciesConcentration,iga::PortQuantity::SpeciesFlux};
		port.species={"tracer"};ports_.push_back(port);
		model_.port_ids={"port"};model_.species_ids={"tracer"};
		committed_.volume_m3=2.0;committed_.amount_mol={{"tracer",10.0}};
	}
	const std::string& DomainId()const noexcept override{return id_;}
	iga::DomainKind Kind()const noexcept override{return iga::DomainKind::OneDFlow;}
	const std::vector<iga::CouplingPort>& Ports()const noexcept override{return ports_;}
	double CommittedAmount()const{return committed_.amount_mol.at("tracer");}
	void SetFlow(double outward_flow_m3_s)
	{
		if(phase_!=Phase::Committed||!std::isfinite(outward_flow_m3_s))
			throw std::runtime_error("controlled source flow change needs committed phase");
		flow_m3_s_=outward_flow_m3_s;
	}
	void BeginStep(const iga::DomainStepContext& step)override
	{
		if(phase_!=Phase::Committed)throw std::runtime_error("source begin phase");
		step.Validate();step_=step;pressure_.reset();donor_.reset();trial_.reset();
		phase_=Phase::Ready;
	}
	void SetPortInput(const std::string& port,const iga::PortBoundaryData& data)override
	{
		if(phase_!=Phase::Ready||port!="port"||data.time_s!=step_.EndTime()
			||!data.mean_pressure_pa||data.outward_flow_m3_s)
			throw std::runtime_error("source pressure input");
		pressure_=*data.mean_pressure_pa;
	}
	void SolveHydraulicTrial()override
	{
		if(phase_!=Phase::Ready||!pressure_)throw std::runtime_error("source hydraulic phase");
		phase_=Phase::Hydraulic;
	}
	iga::PortState GetHydraulicPortState(const std::string& port)const override
	{
		if(port!="port"||(phase_!=Phase::Hydraulic&&phase_!=Phase::Transport
			&&phase_!=Phase::Prepared))throw std::runtime_error("source hydraulic state");
		iga::PortState result;result.time_s=step_.EndTime();
		result.mean_pressure_pa=*pressure_;result.outward_flow_m3_s=flow_m3_s_;
		return result;
	}
	void RollbackHydraulicTrial()override
	{
		if(phase_!=Phase::Hydraulic)throw std::runtime_error("source hydraulic rollback");
		pressure_.reset();phase_=Phase::Ready;
	}
	void SetTransportConcentration(const std::string& port,double time_s,
		const std::map<std::string,double>& concentration)override
	{
		if(phase_!=Phase::Hydraulic||port!="port"||time_s!=step_.EndTime()
			||flow_m3_s_>=0.0||donor_||concentration.size()!=1
			||!concentration.count("tracer")
			||!(concentration.at("tracer")>=0.0)
			||!std::isfinite(concentration.at("tracer")))
			throw std::runtime_error("controlled source donor is invalid");
		donor_=concentration;
	}
	void SolveTransportTrial()override
	{
		if(phase_!=Phase::Hydraulic)throw std::runtime_error("source transport phase");
		if(flow_m3_s_<0.0&&!donor_)throw std::runtime_error("source inflow needs donor");
		const std::map<std::string,double> empty;
		trial_=iga::EvaluateZeroDSpeciesReservoirTrial(model_,committed_,
			{{"port",{flow_m3_s_,donor_?*donor_:empty}}},step_.dt_s,step_.EndTime());
		phase_=Phase::Transport;
	}
	iga::PortState GetTransportPortState(const std::string& port)const override
	{
		if(port!="port"||(phase_!=Phase::Transport&&phase_!=Phase::Prepared))
			throw std::runtime_error("source transport state");
		auto result=trial_->ports.at("port");
		result.mean_pressure_pa=*pressure_;return result;
	}
	void RollbackTransportTrial()override
	{
		if(phase_!=Phase::Transport)throw std::runtime_error("source transport rollback");
		trial_.reset();donor_.reset();phase_=Phase::Hydraulic;
	}
	std::map<std::string,iga::SpeciesStepAccounting> GetSpeciesStepAccounting()const override
	{
		if(phase_!=Phase::Transport&&phase_!=Phase::Prepared)
			throw std::runtime_error("source accounting phase");
		return trial_->balance.species;
	}
	void SolveTrial()override{SolveHydraulicTrial();SolveTransportTrial();}
	iga::PortState GetPortState(const std::string& port)const override
	{return GetTransportPortState(port);}
	void RollbackTrial()override{RollbackTransportTrial();RollbackHydraulicTrial();}
	void AbortStep()override
	{
		if(phase_==Phase::Committed)throw std::runtime_error("source abort phase");
		pressure_.reset();donor_.reset();trial_.reset();prepared_.reset();phase_=Phase::Committed;
	}
	void PrepareCommitStep()override
	{
		if(phase_!=Phase::Transport)throw std::runtime_error("source prepare phase");
		prepared_=trial_->state;phase_=Phase::Prepared;
	}
	void FinalizeCommitStep()noexcept override
	{
		if(phase_!=Phase::Prepared||!prepared_)std::terminate();
		committed_.volume_m3=prepared_->volume_m3;
		committed_.amount_mol.swap(prepared_->amount_mol);
		pressure_.reset();donor_.reset();trial_.reset();prepared_.reset();phase_=Phase::Committed;
	}
private:
	enum class Phase{Committed,Ready,Hydraulic,Transport,Prepared};
	std::string id_="source";
	std::vector<iga::CouplingPort> ports_;
	iga::ZeroDSpeciesReservoirModel model_;
	iga::ZeroDSpeciesReservoirState committed_;
	iga::DomainStepContext step_;
	std::optional<double> pressure_;
	std::optional<std::map<std::string,double>> donor_;
	std::optional<iga::ZeroDSpeciesReservoirTrial> trial_;
	std::optional<iga::ZeroDSpeciesReservoirState> prepared_;
	Phase phase_=Phase::Committed;
	double flow_m3_s_=2.0;
};

bool Near(double first,double second)
{
	return std::abs(first-second)<=1e-12
		*std::max(1.0,std::max(std::abs(first),std::abs(second)));
}

} // namespace

int main()
{
	iga::ZeroDFlowModel flow;
	flow.role=iga::ZeroDFlowRole::TerminalRcr;
	flow.terminal={1.5,5.0,2.0,7.0};
	iga::ZeroDSpeciesReservoirModel species;
	species.port_ids={"graph","distal"};species.species_ids={"tracer"};
	auto source=std::make_unique<ControlledSource>();
	auto terminal=std::make_unique<iga::ZeroDTerminalRcrSpeciesDomainRuntime>(
		"terminal",flow,iga::ZeroDFlowState{9.0},species,
		iga::ZeroDSpeciesReservoirState{1.0,{{"tracer",2.0}}});
	auto* source_ptr=source.get();auto* terminal_ptr=terminal.get();
	iga::SimulationGraph graph({
		{"source",iga::DomainKind::OneDFlow,source->Ports(),{{"tracer","source_tracer"}}},
		{"terminal",iga::DomainKind::ZeroDFlow,terminal->Ports(),
			{{"tracer","terminal_tracer"}}}},
		{{"edge",{"source","port"},{"terminal","port"},
			iga::CouplingLaw::PressureFlow,{"tracer"}}},{{"tracer","mol/m3"}});
	std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> owners;
	owners.push_back(std::move(source));owners.push_back(std::move(terminal));
	iga::DomainRuntimeRegistry registry(graph,std::move(owners));
	iga::SpeciesPressureFlowExecutionControls controls;
	controls.amount_tolerances.emplace("tracer",iga::SpeciesAmountTolerance{});
	iga::SpeciesPressureFlowComponentExecutor executor(registry,"source",controls);
	bool rejected=false;
	try{
		executor.Advance({0,0.0,0.1},{{"edge",9.0}},
			[](const auto&){throw std::runtime_error("injected precommit failure");});
	}catch(const std::runtime_error&){rejected=true;}
	assert(rejected);
	assert(Near(source_ptr->CommittedAmount(),10.0));
	assert(Near(terminal_ptr->CommittedSpeciesState().amount_mol.at("tracer"),2.0));
	assert(Near(terminal_ptr->CommittedHydraulicState().stored_pressure_pa,9.0));
	const auto result=executor.Advance({0,0.0,0.1},{{"edge",9.0}});
	assert(result.edge_amounts.size()==1);
	assert(Near(result.edge_amounts.front().first_outward_amount,1.0));
	assert(Near(result.edge_amounts.front().second_outward_amount,-1.0));
	assert(Near(result.global_balances.at("tracer").residual,0.0));
	assert(terminal_ptr->CommittedStepIndex()==0);
	assert(source_ptr->CommittedAmount()<10.0);
	assert(terminal_ptr->CommittedSpeciesState().amount_mol.at("tracer")>2.0);
	source_ptr->SetFlow(-0.75);
	const auto reverse=executor.Advance({1,0.1,0.1},{{"edge",9.0}});
	assert(reverse.edge_amounts.size()==1);
	assert(reverse.edge_amounts.front().first_outward_amount<0.0);
	assert(reverse.edge_amounts.front().second_outward_amount>0.0);
	assert(Near(reverse.edge_amounts.front().residual,0.0));
	assert(Near(reverse.global_balances.at("tracer").residual,0.0));
	assert(terminal_ptr->CommittedStepIndex()==1);
	assert(source_ptr->CommittedAmount()>9.0);
	return 0;
}

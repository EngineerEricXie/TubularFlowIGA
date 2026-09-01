#ifndef IGA_ONE_D_FLOW_DOMAIN_ADAPTER_HPP
#define IGA_ONE_D_FLOW_DOMAIN_ADAPTER_HPP

#include "CoupledDomainRuntime.hpp"
#include "OneDRuntime.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class OneDInletPolicy { ConfiguredOpenLoop, CoupledRoot };

class OneDFlowDomainAdapter : public CoupledDomainRuntime {
public:
	OneDFlowDomainAdapter(std::string domain_id, OneDFlowRuntime& runtime,
		std::vector<CouplingPort> ports, OneDInletPolicy inlet_policy)
		: domain_id_(std::move(domain_id)), runtime_(runtime), ports_(std::move(ports)),
		  inlet_policy_(inlet_policy)
	{
		if (domain_id_.empty()) throw std::runtime_error("1D domain adapter id must be nonempty");
		ValidateCouplingPorts(ports_);
		int root_flow_receivers = 0;
		for (const auto& port : ports_) {
			if (port.subsystem_id != domain_id_ || port.locator_kind != "runtime_port")
				throw std::runtime_error("1D domain adapter requires matching runtime_port metadata");
			if (port.orientation.native_to_outward_sign != 1)
				throw std::runtime_error(
					"1D runtime_port metadata must use canonical outward orientation +1");
			if (port.locator == "root") {
				if (!port.requires.empty()
					&& port.requires != std::set<PortQuantity>{PortQuantity::FlowRate})
					throw std::runtime_error("1D root adapter port must receive flow_rate");
				if (port.requires.count(PortQuantity::FlowRate)) ++root_flow_receivers;
			} else if (port.locator.compare(0, 7, "outlet:") == 0) {
				if (!port.requires.empty()
					&& port.requires != std::set<PortQuantity>{PortQuantity::MeanPressure})
					throw std::runtime_error("1D outlet adapter port must receive mean_pressure");
			} else {
				throw std::runtime_error("1D domain adapter runtime port must be root or outlet:<node-id>");
			}
		}
		if (inlet_policy_ == OneDInletPolicy::ConfiguredOpenLoop && root_flow_receivers != 0)
			throw std::runtime_error("configured-open-loop 1D domain cannot expose a coupled root flow receiver");
		if (inlet_policy_ == OneDInletPolicy::CoupledRoot && root_flow_receivers != 1)
			throw std::runtime_error("coupled-root 1D domain requires exactly one root flow receiver");
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::OneDFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		step.Validate();
		runtime_.BeginStep(step.start_time_s, step.dt_s);
		if (inlet_policy_ == OneDInletPolicy::ConfiguredOpenLoop)
			runtime_.SetConfiguredOpenLoopInlet();
	}

	void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) override
	{
		runtime_.SetPortInput(Port(port_id).locator, input);
	}

	void SolveTrial() override { runtime_.SolveTrial(); }

	PortState GetPortState(const std::string& port_id) const override
	{
		return runtime_.GetPortState(Port(port_id).locator);
	}

	void RollbackTrial() override { runtime_.RollbackTrial(); }
	void AbortStep() override { runtime_.AbortStep(); }
	void PrepareCommitStep() override { runtime_.PrepareCommitStep(); }
	void FinalizeCommitStep() noexcept override { runtime_.FinalizeCommitStep(); }

private:
	const CouplingPort& Port(const std::string& port_id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(),
			[&](const CouplingPort& port) { return port.id == port_id; });
		if (found == ports_.end())
			throw std::runtime_error("1D domain adapter has no port '"+port_id+"'");
		return *found;
	}

	std::string domain_id_;
	OneDFlowRuntime& runtime_;
	std::vector<CouplingPort> ports_;
	OneDInletPolicy inlet_policy_;
};

} // namespace iga

#endif

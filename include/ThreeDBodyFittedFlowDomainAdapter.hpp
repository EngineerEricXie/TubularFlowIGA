#ifndef IGA_THREE_D_BODY_FITTED_FLOW_DOMAIN_ADAPTER_HPP
#define IGA_THREE_D_BODY_FITTED_FLOW_DOMAIN_ADAPTER_HPP

#include "CoupledDomainRuntime.hpp"
#include "TemporalFunction.hpp"
#include "ThreeDFlowCoupling.hpp"
#include "TransientFlowRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ThreeDFlowDomainControls {
	int maximum_newton = 12;
	double nonlinear_relative_tolerance = 1.0e-8;
	double nonlinear_absolute_tolerance = 1.0e-8;
	double mass_relative_tolerance = 3.0;

	void Validate() const
	{
		if (maximum_newton < 1 || !(nonlinear_relative_tolerance > 0.0)
			|| !(nonlinear_absolute_tolerance > 0.0) || !(mass_relative_tolerance > 0.0)
			|| !std::isfinite(nonlinear_relative_tolerance)
			|| !std::isfinite(nonlinear_absolute_tolerance)
			|| !std::isfinite(mass_relative_tolerance))
			throw std::runtime_error("3D flow domain controls require positive finite tolerances");
	}
};

class ThreeDBodyFittedFlowDomainAdapter : public CoupledDomainRuntime {
public:
	ThreeDBodyFittedFlowDomainAdapter(std::string domain_id, TransientFlowRuntime& runtime,
		std::vector<CouplingPort> ports, SimulationConfiguration base_configuration,
		std::filesystem::path case_directory,
		std::map<std::string, double> reference_outward_flow_m3_s,
		ThreeDFlowDomainControls controls = {})
		: domain_id_(std::move(domain_id)), runtime_(runtime), ports_(std::move(ports)),
		  base_configuration_(std::move(base_configuration)),
		  case_directory_(std::move(case_directory)),
		  reference_outward_flow_m3_s_(std::move(reference_outward_flow_m3_s)),
		  controls_(controls)
	{
		ValidateThreeDBodyFittedFlowDomainMetadata(domain_id_, ports_);
		controls_.Validate();
		if (base_configuration_.equation_systems.size() != 1
			|| base_configuration_.equation_systems.front().kind != EquationKind::NavierStokes)
			throw std::runtime_error("3D domain adapter requires one Navier-Stokes equation system");
		for (const auto& port : ports_) {
			if (port.requires.count(PortQuantity::FlowRate)) {
				const auto reference = reference_outward_flow_m3_s_.find(port.id);
				if (reference == reference_outward_flow_m3_s_.end()
					|| reference->second == 0.0 || !std::isfinite(reference->second))
					throw std::runtime_error("3D flow receiver requires a finite nonzero reference flow");
			}
		}
		for (const auto& reference : reference_outward_flow_m3_s_) {
			const auto& port = Port(reference.first);
			if (!port.requires.count(PortQuantity::FlowRate) || reference.second == 0.0
				|| !std::isfinite(reference.second))
				throw std::runtime_error("3D domain adapter reference flow requires a finite nonzero flow receiver");
		}
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ThreeDBodyFittedFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		step.Validate();
		step_ = step;
		inputs_.clear();
		runtime_.BeginStep(step.step_index, step.EndTime(), controls_.maximum_newton,
			controls_.nonlinear_relative_tolerance, controls_.nonlinear_absolute_tolerance,
			controls_.mass_relative_tolerance);
	}

	void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) override
	{
		if (runtime_.Phase() != FlowStepPhase::TrialReady)
			throw std::runtime_error("3D domain adapter input requires an active trial-ready step");
		const auto& port = Port(port_id);
		ValidatePortBoundaryData(input);
		const auto time_scale = std::max({1.0, std::abs(step_.EndTime()), std::abs(input.time_s)});
		if (std::abs(input.time_s-step_.EndTime()) > 1.0e-12*time_scale)
			throw std::runtime_error("3D domain adapter input time does not match the active step");
		const int supplied = static_cast<int>(input.outward_flow_m3_s.has_value())
			+static_cast<int>(input.mean_pressure_pa.has_value())
			+static_cast<int>(input.mean_normal_traction_pa.has_value())
			+static_cast<int>(input.total_pressure_pa.has_value());
		if (supplied != 1 || !input.concentration.empty() || !input.outward_species_flux.empty())
			throw std::runtime_error("3D flow domain input requires exactly one flow or pressure quantity");
		if ((input.outward_flow_m3_s && !port.requires.count(PortQuantity::FlowRate))
			|| (input.mean_pressure_pa && !port.requires.count(PortQuantity::MeanPressure))
			|| (input.mean_normal_traction_pa
				&& !port.requires.count(PortQuantity::MeanNormalTraction))
			|| input.total_pressure_pa)
			throw std::runtime_error("3D flow domain input quantity is not required by its port");
		inputs_[port_id] = input;
	}

	void SolveTrial() override
	{
		if (runtime_.Phase() != FlowStepPhase::TrialReady)
			throw std::runtime_error("3D domain adapter solve requires an active trial-ready step");
		for (const auto& port : ports_)
			for (const auto quantity : port.requires)
				if ((quantity == PortQuantity::FlowRate || quantity == PortQuantity::MeanPressure
					|| quantity == PortQuantity::MeanNormalTraction) && !inputs_.count(port.id))
					throw std::runtime_error("3D flow domain is missing trial input for port '"
						+port.id+"'");
		auto trial_configuration = MaterializeBoundaryWaveforms(base_configuration_,
			case_directory_.string(), step_.EndTime());
		for (const auto& input : inputs_) {
			const auto& port = Port(input.first);
			if (input.second.outward_flow_m3_s)
				ApplyThreeDReferenceProfileInput(trial_configuration,
					trial_configuration.equation_systems.front(), port, input.second,
					reference_outward_flow_m3_s_.at(port.id));
		}
		runtime_.SetTrialBoundaryConfiguration(trial_configuration);
		for (const auto& input : inputs_)
			if (!input.second.outward_flow_m3_s)
				runtime_.SetPortInput(Port(input.first), input.second);
		runtime_.SolveTrial();
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		return runtime_.GetPortState(Port(port_id));
	}

	void RollbackTrial() override { runtime_.RollbackTrial(); }
	void AbortStep() override
	{
		inputs_.clear();
		runtime_.AbortStep();
	}
	void PrepareCommitStep() override { runtime_.PrepareCommitStep(); }
	void FinalizeCommitStep() noexcept override
	{
		inputs_.clear();
		runtime_.FinalizeCommitStep();
	}

private:
	const CouplingPort& Port(const std::string& port_id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(),
			[&](const CouplingPort& port) { return port.id == port_id; });
		if (found == ports_.end())
			throw std::runtime_error("3D domain adapter has no port '"+port_id+"'");
		return *found;
	}

	std::string domain_id_;
	TransientFlowRuntime& runtime_;
	std::vector<CouplingPort> ports_;
	SimulationConfiguration base_configuration_;
	std::filesystem::path case_directory_;
	std::map<std::string, double> reference_outward_flow_m3_s_;
	ThreeDFlowDomainControls controls_;
	DomainStepContext step_;
	std::map<std::string, PortBoundaryData> inputs_;
};

} // namespace iga

#endif

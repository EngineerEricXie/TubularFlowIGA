#ifndef IGA_THREE_D_IMMERSED_FLOW_DOMAIN_HPP
#define IGA_THREE_D_IMMERSED_FLOW_DOMAIN_HPP

// The immersed backend is intentionally quasi-static in Phase 6: the domain
// step clock is coupling metadata while the assembled Navier--Stokes operator
// retains its construction-time dt=0 contract.
#include "CoupledDomainRuntime.hpp"
#include "FlowDomainPortMetadata.hpp"
#include "../solvers/cpu/include/ImmersedStaticFlowRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ThreeDImmersedFlowDomainDiagnostics {
	double committed_time_s = 0.0;
	int committed_step_index = -1;
	std::size_t committed_steps = 0;
	std::size_t rollback_count = 0;
	std::size_t abort_count = 0;
	std::size_t prepared_count = 0;
};

class ThreeDImmersedFlowDomain final : public CoupledDomainRuntime {
public:
	ThreeDImmersedFlowDomain(std::string domain_id, ImmersedStaticFlowRuntime& runtime,
		std::vector<CouplingPort> ports, double initial_time_s = 0.0)
		: domain_id_(std::move(domain_id)), runtime_(runtime), ports_(std::move(ports))
	{
		if (!std::isfinite(initial_time_s))
			throw std::invalid_argument("immersed flow domain initial time must be finite");
		ValidateThreeDImmersedFlowDomainMetadata(domain_id_, ports_);
		ValidateRuntimeBinding();
		diagnostics_.committed_time_s = initial_time_s;
		committed_controls_ = ControlValues();
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ThreeDImmersedFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }
	const ThreeDImmersedFlowDomainDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const std::map<std::string, PortState>& CommittedPortStates() const noexcept
	{ return committed_measurements_; }
	std::vector<PetscScalar> CommittedBackendState() const { return runtime_.CommittedState(); }

	void BeginStep(const DomainStepContext& step) override
	{
		RequirePhase(Phase::Committed, "begin"); step.Validate();
		if (step.start_time_s != diagnostics_.committed_time_s)
			throw std::runtime_error("immersed flow domain step start time is not the committed time");
		if (step.step_index != diagnostics_.committed_step_index+1)
			throw std::runtime_error("immersed flow domain step index is not the next committed index");
		step_ = step; inputs_.clear(); trial_measurements_.clear();
		trial_controls_.clear(); phase_ = Phase::TrialReady;
	}

	void SetPortInput(const std::string& port_id, const PortBoundaryData& input) override
	{
		RequirePhase(Phase::TrialReady, "set input");
		const auto& port = Port(port_id); ValidatePortBoundaryData(input);
		if (input.time_s != step_.EndTime())
			throw std::runtime_error("immersed flow domain input time does not match the active step end");
		if (inputs_.count(port_id)) throw std::runtime_error("immersed flow domain port input was supplied twice");
		const int fields = static_cast<int>(input.outward_flow_m3_s.has_value())
			+static_cast<int>(input.mean_pressure_pa.has_value())
			+static_cast<int>(input.mean_normal_traction_pa.has_value())
			+static_cast<int>(input.total_pressure_pa.has_value());
		if (fields != 1 || input.total_pressure_pa || !input.concentration.empty()
			|| !input.outward_species_flux.empty())
			throw std::runtime_error("immersed flow domain input requires exactly one hydraulic field");
		const PortQuantity quantity = input.outward_flow_m3_s ? PortQuantity::FlowRate
			: input.mean_pressure_pa ? PortQuantity::MeanPressure : PortQuantity::MeanNormalTraction;
		if (port.requires.size() != 1 || !port.requires.count(quantity))
			throw std::runtime_error("immersed flow domain input quantity is not declared by its port");
		inputs_.emplace(port_id, input);
	}

	void SolveTrial() override
	{
		RequirePhase(Phase::TrialReady, "solve");
		for (const auto& port : ports_)
			if (!inputs_.count(port.id))
				throw std::runtime_error("immersed flow domain is missing input for controlled port '"+port.id+"'");
		const auto saved = committed_controls_;
		try {
			trial_controls_.clear();
			for (const auto& port : ports_) {
				const double value = InputValue(inputs_.at(port.id));
				runtime_.SetPortControlValue(port.id, value);
				trial_controls_.emplace(port.id, value);
			}
			if (!runtime_.SolveTrial()) throw std::runtime_error("immersed flow trial did not converge");
			std::map<std::string, PortState> measured;
			for (const auto& port : ports_) measured.emplace(port.id, MakePortState(port));
			trial_measurements_.swap(measured);
			phase_ = Phase::TrialSolved;
		} catch (...) {
			// SolveTrial internally rolls its PETSc image back on failure.  Restore
			// the mutable control values as well, leaving this trial retryable.
			if (runtime_.Diagnostics().trial_active) runtime_.Rollback();
			RestoreControls(saved); trial_controls_.clear(); trial_measurements_.clear();
			phase_ = Phase::TrialReady;
			throw;
		}
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::Prepared)
			throw std::runtime_error("immersed flow domain has no solved trial port state");
		const auto found = trial_measurements_.find(port_id);
		if (found == trial_measurements_.end()) throw std::runtime_error("immersed flow domain has no port '"+port_id+"'");
		return found->second;
	}

	void RollbackTrial() override
	{
		if (phase_ != Phase::TrialSolved)
			throw std::runtime_error("immersed flow domain rollback requires a solved trial");
		if (runtime_.Diagnostics().trial_active) runtime_.Rollback();
		RestoreControls(committed_controls_); inputs_.clear(); trial_controls_.clear();
		trial_measurements_.clear(); phase_ = Phase::TrialReady; ++diagnostics_.rollback_count;
	}

	void AbortStep() override
	{
		if (phase_ != Phase::TrialReady && phase_ != Phase::TrialSolved && phase_ != Phase::Prepared)
			throw std::runtime_error("immersed flow domain abort requires an active trial");
		if (phase_ == Phase::Prepared) runtime_.AbortPrepared();
		if (runtime_.Diagnostics().trial_active) runtime_.Rollback();
		RestoreControls(committed_controls_); inputs_.clear(); trial_controls_.clear();
		trial_measurements_.clear(); prepared_measurements_.clear(); phase_ = Phase::Committed;
		++diagnostics_.abort_count;
	}

	void PrepareCommitStep() override
	{
		RequirePhase(Phase::TrialSolved, "prepare commit");
		// Allocation/validation precedes backend preparation.  Thus a throwing
		// port-state copy cannot publish either half of the transaction.
		std::map<std::string, PortState> prepared;
		for (const auto& state : trial_measurements_) { ValidatePortState(state.second); prepared.emplace(state); }
		runtime_.PrepareCommit();
		prepared_measurements_.swap(prepared); phase_ = Phase::Prepared;
		++diagnostics_.prepared_count;
	}

	void FinalizeCommitStep() noexcept override
	{
		if (phase_ != Phase::Prepared) return;
		// No PETSc call or allocation is permitted below this point.
		runtime_.FinalizeCommit();
		committed_measurements_.swap(prepared_measurements_);
		committed_controls_.swap(trial_controls_);
		diagnostics_.committed_time_s = step_.EndTime();
		diagnostics_.committed_step_index = step_.step_index;
		++diagnostics_.committed_steps;
		inputs_.clear(); trial_measurements_.clear(); phase_ = Phase::Committed;
	}

private:
	enum class Phase { Committed, TrialReady, TrialSolved, Prepared };

	void RequirePhase(Phase expected, const char* operation) const
	{
		if (phase_ != expected)
			throw std::runtime_error(std::string("immersed flow domain cannot ")+operation+" in its current lifecycle phase");
	}
	const CouplingPort& Port(const std::string& id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(), [&id](const CouplingPort& port) { return port.id == id; });
		if (found == ports_.end()) throw std::runtime_error("immersed flow domain has no port '"+id+"'");
		return *found;
	}
	static PortQuantity QuantityFor(ImmersedFlowPortControlMode mode)
	{
		if (mode == ImmersedFlowPortControlMode::FlowRate) return PortQuantity::FlowRate;
		if (mode == ImmersedFlowPortControlMode::Pressure) return PortQuantity::MeanPressure;
		if (mode == ImmersedFlowPortControlMode::MeanNormalTraction) return PortQuantity::MeanNormalTraction;
		throw std::runtime_error("immersed total-pressure control is unsupported");
	}
	void ValidateRuntimeBinding() const
	{
		const auto& definitions = runtime_.PortDefinitions();
		if (definitions.size() != ports_.size()) throw std::runtime_error("immersed runtime and coupling-port counts differ");
		for (const auto& definition : definitions) {
			const auto& port = Port(definition.id);
			if (ParseThreeDImmersedFlowBoundaryLabel(port) != definition.boundary_label
				|| port.requires != std::set<PortQuantity>{QuantityFor(definition.control_mode)})
				throw std::runtime_error("immersed runtime port definition does not exactly match coupling metadata");
		}
	}
	std::map<std::string, double> ControlValues() const
	{
		std::map<std::string, double> values;
		for (const auto& definition : runtime_.PortDefinitions()) values.emplace(definition.id, definition.value);
		return values;
	}
	void RestoreControls(const std::map<std::string, double>& controls)
	{
		for (const auto& control : controls) runtime_.SetPortControlValue(control.first, control.second);
	}
	static double InputValue(const PortBoundaryData& input)
	{
		const double result = input.outward_flow_m3_s ? *input.outward_flow_m3_s
			: input.mean_pressure_pa ? *input.mean_pressure_pa : *input.mean_normal_traction_pa;
		if (!std::isfinite(result)) throw std::runtime_error("immersed flow input is not finite");
		return result;
	}
	PortState MakePortState(const CouplingPort& port) const
	{
		const auto found = std::find_if(runtime_.Diagnostics().ports.begin(), runtime_.Diagnostics().ports.end(), [&port](const auto& item) { return item.id == port.id; });
		if (found == runtime_.Diagnostics().ports.end()) throw std::runtime_error("immersed runtime did not measure a configured port");
		PortState result; result.time_s = step_.EndTime();
		if (port.provides.count(PortQuantity::Area)) result.area_m2 = found->measurement.area_m2;
		if (port.provides.count(PortQuantity::FlowRate)) result.outward_flow_m3_s = found->measurement.outward_flow_m3_s;
		if (port.provides.count(PortQuantity::MeanPressure)) result.mean_pressure_pa = found->measurement.mean_pressure_pa;
		if (port.provides.count(PortQuantity::MeanNormalTraction)) result.mean_normal_traction_pa = found->measurement.mean_normal_traction_pa;
		ValidatePortState(result); return result;
	}

	std::string domain_id_;
	ImmersedStaticFlowRuntime& runtime_;
	std::vector<CouplingPort> ports_;
	Phase phase_ = Phase::Committed;
	DomainStepContext step_;
	ThreeDImmersedFlowDomainDiagnostics diagnostics_;
	std::map<std::string, PortBoundaryData> inputs_;
	std::map<std::string, double> committed_controls_, trial_controls_;
	std::map<std::string, PortState> trial_measurements_, prepared_measurements_, committed_measurements_;
};

} // namespace iga

#endif

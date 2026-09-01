#ifndef IGA_THREE_D_BODY_FITTED_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP
#define IGA_THREE_D_BODY_FITTED_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP

#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "TransientTransportRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ThreeDFlowTransportDomainControls {
	ThreeDFlowDomainControls flow;
	double species_flow_epsilon_m3_s = 0.0;

	void Validate() const
	{
		flow.Validate();
		if (!(species_flow_epsilon_m3_s >= 0.0)
			|| !std::isfinite(species_flow_epsilon_m3_s))
			throw std::runtime_error(
				"3D flow/transport domain requires a finite nonnegative species flow epsilon");
	}
};

inline void SetThreeDPortSpeciesBoundary(SimulationConfiguration& configuration,
	const CouplingPort& port, const std::string& native_field,
	const std::optional<double>& concentration)
{
	const int label = ParseThreeDFlowBoundaryLabel(port);
	auto boundary = std::find_if(configuration.boundaries.begin(),
		configuration.boundaries.end(),
		[label](const NamedBoundaryDefinition& value) { return value.label == label; });
	if (boundary == configuration.boundaries.end())
		throw std::runtime_error("3D coupled species port boundary label is not configured");
	auto condition = boundary->conditions.end();
	for (auto candidate = boundary->conditions.begin();
		candidate != boundary->conditions.end(); ++candidate)
		if (candidate->field == native_field) {
			if (condition != boundary->conditions.end())
				throw std::runtime_error(
					"3D coupled species boundary has duplicate native field conditions");
			condition = candidate;
		}
	if (condition == boundary->conditions.end()) {
		boundary->conditions.push_back({});
		condition = std::prev(boundary->conditions.end());
		condition->field = native_field;
	} else if (condition->kind != FieldBoundaryKind::Dirichlet
		&& condition->kind != FieldBoundaryKind::NoFlux
		&& condition->kind != FieldBoundaryKind::AdvectiveOutflow)
		throw std::runtime_error(
			"3D coupled species boundary may switch only among Dirichlet, no_flux, and advective_outflow");
	condition->profile.clear();
	condition->waveform.clear();
	condition->scale = 1.0;
	if (concentration) {
		condition->kind = FieldBoundaryKind::Dirichlet;
		condition->value = {*concentration};
	} else {
		condition->kind = FieldBoundaryKind::AdvectiveOutflow;
		condition->value.clear();
	}
}

class ThreeDBodyFittedFlowTransportDomainAdapter : public CoupledDomainRuntime {
public:
	ThreeDBodyFittedFlowTransportDomainAdapter(std::string domain_id,
		TransientFlowRuntime& flow_runtime, TransientTransportRuntime& transport_runtime,
		std::vector<CouplingPort> ports, SimulationConfiguration base_configuration,
		std::filesystem::path case_directory,
		std::map<std::string, std::string> species_bindings,
		std::map<std::string, double> reference_outward_flow_m3_s,
		ThreeDFlowTransportDomainControls controls = {})
		: domain_id_(std::move(domain_id)), flow_runtime_(flow_runtime),
		  transport_runtime_(transport_runtime), ports_(std::move(ports)),
		  base_configuration_(std::move(base_configuration)),
		  case_directory_(std::move(case_directory)),
		  species_bindings_(std::move(species_bindings)), controls_(controls),
		  flow_adapter_(domain_id_, flow_runtime_, ports_, base_configuration_,
			case_directory_, std::move(reference_outward_flow_m3_s), controls_.flow)
	{
		controls_.Validate();
		if (flow_runtime_.RequiredNodes() != transport_runtime_.RequiredNodes())
			throw std::runtime_error(
				"3D flow and transport runtimes require identical ordered node sets");
		std::set<std::string> native_fields;
		for (const auto& binding : species_bindings_) {
			if (binding.first.empty() || binding.second.empty()
				|| !transport_runtime_.System().field_index.count(binding.second))
				throw std::runtime_error(
					"3D species binding must name a transported native scalar field");
			if (!native_fields.insert(binding.second).second)
				throw std::runtime_error(
					"3D species bindings require unique transported native fields");
		}
		for (const auto& port : ports_)
			for (const auto& species : port.species)
				if (!species_bindings_.count(species))
					throw std::runtime_error("3D coupled port species '"+species
						+"' has no native transport binding");
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ThreeDBodyFittedFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }

	void BeginStep(const DomainStepContext& step) override
	{
		step.Validate();
		flow_adapter_.BeginStep(step);
		try {
			transport_runtime_.BeginStep();
		} catch (...) {
			flow_adapter_.AbortStep();
			throw;
		}
		step_ = step;
		inputs_.clear();
		trial_solve_succeeded_ = false;
	}

	void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) override
	{
		if (flow_runtime_.Phase() != FlowStepPhase::TrialReady
			|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen)
			throw std::runtime_error(
				"3D flow/transport input requires an active trial-open step");
		const auto& port = Port(port_id);
		ValidatePortBoundaryData(input);
		const auto time_scale = std::max({1.0, std::abs(step_.EndTime()),
			std::abs(input.time_s)});
		if (std::abs(input.time_s-step_.EndTime()) > 1.0e-12*time_scale)
			throw std::runtime_error(
				"3D flow/transport input time does not match the active step");
		ValidateInputCapabilities(port, input);
		const auto found = inputs_.find(port_id);
		if (found == inputs_.end()) inputs_.emplace(port_id, input);
		else MergeInput(found->second, input);
	}

	void SolveTrial() override
	{
		if (flow_runtime_.Phase() != FlowStepPhase::TrialReady
			|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen)
			throw std::runtime_error(
				"3D flow/transport solve requires trial-ready runtimes");
		trial_solve_succeeded_ = false;
		for (const auto& port : ports_) {
			const auto found = inputs_.find(port.id);
			const int required_hydraulic = HydraulicRequirementCount(port);
			const int supplied_hydraulic = found == inputs_.end()
				? 0 : HydraulicValueCount(found->second);
			if (required_hydraulic != supplied_hydraulic)
				throw std::runtime_error("3D flow/transport domain is missing or has an extra hydraulic input for port '"
					+port.id+"'");
			if (supplied_hydraulic == 1)
				flow_adapter_.SetPortInput(port.id, HydraulicInput(found->second));
		}
		flow_adapter_.SolveTrial();

		auto trial_configuration = MaterializeBoundaryWaveforms(base_configuration_,
			case_directory_, step_.EndTime());
		for (const auto& port : ports_) {
			if (port.species.empty()) continue;
			const auto input = inputs_.find(port.id);
			const std::map<std::string, double> empty;
			const auto& concentration = input == inputs_.end()
				? empty : input->second.concentration;
			const bool has_concentration = !concentration.empty();
			if (has_concentration && SpeciesKeys(concentration) != port.species)
				throw std::runtime_error("3D coupled species concentration set is incomplete for port '"
					+port.id+"'");
			const auto flow = flow_adapter_.GetPortState(port.id);
			if (!flow.outward_flow_m3_s)
				throw std::runtime_error("3D coupled species port omitted outward flow");
			if (*flow.outward_flow_m3_s < -controls_.species_flow_epsilon_m3_s
				&& !has_concentration)
				throw std::runtime_error("3D inward species port '"+port.id
					+"' requires every coupled concentration");
			if (*flow.outward_flow_m3_s > controls_.species_flow_epsilon_m3_s
				&& has_concentration)
				throw std::runtime_error("3D outward species port '"+port.id
					+"' must not impose a remote concentration");
			for (const auto& species : port.species) {
				std::optional<double> value;
				if (has_concentration) value = concentration.at(species);
				SetThreeDPortSpeciesBoundary(trial_configuration, port,
					species_bindings_.at(species), value);
			}
		}
		transport_runtime_.SolveTrial(trial_configuration,
			flow_runtime_.RequiredNodes(), flow_runtime_.GatherRequiredVelocity());
		trial_solve_succeeded_ = true;
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		const auto& port = Port(port_id);
		if (port.species.empty()) return flow_adapter_.GetPortState(port_id);
		if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
			|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved
			|| !trial_solve_succeeded_)
			throw std::runtime_error(
				"3D flow/transport port state requires successful flow and transport trials");
		auto state = flow_runtime_.MeasurePorts({port}, step_.EndTime(),
			transport_runtime_.System().fields, transport_runtime_.GatherRequiredState(),
			&transport_runtime_.System()).at(port.id);
		std::map<std::string, double> logical_concentration;
		std::map<std::string, double> logical_flux;
		for (const auto& species : port.species) {
			const auto& native = species_bindings_.at(species);
			logical_concentration.emplace(species, state.concentration.at(native));
			logical_flux.emplace(species, state.outward_species_flux.at(native));
		}
		state.concentration = std::move(logical_concentration);
		state.outward_species_flux = std::move(logical_flux);
		ValidatePortState(state);
		return state;
	}

	void RollbackTrial() override
	{
		if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
			|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved)
			throw std::runtime_error(
				"3D flow/transport rollback requires two solved trial phases");
		transport_runtime_.RollbackTrial();
		flow_adapter_.RollbackTrial();
		inputs_.clear();
		trial_solve_succeeded_ = false;
	}

	void AbortStep() override
	{
		std::exception_ptr failure;
		try { transport_runtime_.AbortStep(); }
		catch (...) { failure = std::current_exception(); }
		try { flow_adapter_.AbortStep(); }
		catch (...) { if (!failure) failure = std::current_exception(); }
		inputs_.clear();
		trial_solve_succeeded_ = false;
		if (failure) std::rethrow_exception(failure);
	}

	void PrepareCommitStep() override
	{
		if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
			|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved
			|| !trial_solve_succeeded_)
			throw std::runtime_error(
				"3D flow/transport prepare requires two successful trial solves");
		transport_runtime_.PrepareCommitStep();
		flow_adapter_.PrepareCommitStep();
	}

	void FinalizeCommitStep() noexcept override
	{
		transport_runtime_.FinalizeCommitStep();
		flow_adapter_.FinalizeCommitStep();
		inputs_.clear();
		trial_solve_succeeded_ = false;
	}

private:
	const CouplingPort& Port(const std::string& port_id) const
	{
		const auto found = std::find_if(ports_.begin(), ports_.end(),
			[&](const CouplingPort& port) { return port.id == port_id; });
		if (found == ports_.end())
			throw std::runtime_error("3D flow/transport domain has no port '"+port_id+"'");
		return *found;
	}

	static int HydraulicRequirementCount(const CouplingPort& port)
	{
		return static_cast<int>(port.requires.count(PortQuantity::FlowRate))
			+static_cast<int>(port.requires.count(PortQuantity::MeanPressure))
			+static_cast<int>(port.requires.count(PortQuantity::MeanNormalTraction));
	}

	static int HydraulicValueCount(const PortBoundaryData& input)
	{
		return static_cast<int>(input.outward_flow_m3_s.has_value())
			+static_cast<int>(input.mean_pressure_pa.has_value())
			+static_cast<int>(input.mean_normal_traction_pa.has_value())
			+static_cast<int>(input.total_pressure_pa.has_value());
	}

	static PortBoundaryData HydraulicInput(const PortBoundaryData& input)
	{
		PortBoundaryData result;
		result.time_s = input.time_s;
		result.outward_flow_m3_s = input.outward_flow_m3_s;
		result.mean_pressure_pa = input.mean_pressure_pa;
		result.mean_normal_traction_pa = input.mean_normal_traction_pa;
		result.total_pressure_pa = input.total_pressure_pa;
		return result;
	}

	static std::set<std::string> SpeciesKeys(
		const std::map<std::string, double>& values)
	{
		std::set<std::string> result;
		for (const auto& value : values) result.insert(value.first);
		return result;
	}

	void ValidateInputCapabilities(const CouplingPort& port,
		const PortBoundaryData& input) const
	{
		if ((input.outward_flow_m3_s && !port.requires.count(PortQuantity::FlowRate))
			|| (input.mean_pressure_pa && !port.requires.count(PortQuantity::MeanPressure))
			|| (input.mean_normal_traction_pa
				&& !port.requires.count(PortQuantity::MeanNormalTraction))
			|| input.total_pressure_pa)
			throw std::runtime_error(
				"3D flow/transport input contains an unsupported hydraulic quantity");
		if (!input.concentration.empty()
			&& !port.requires.count(PortQuantity::SpeciesConcentration))
			throw std::runtime_error(
				"3D flow/transport port does not accept species concentration");
		// Total species flux is a measured port output.  The graph executor owns
		// the two-sided conservative residual and must check it before prepare;
		// imposing concentration and peer total flux here would overconstrain the PDE.
		if (!input.outward_species_flux.empty())
			throw std::runtime_error(
				"3D peer species flux is executor-owned residual data, not a boundary input");
		for (const auto& value : input.concentration)
			if (!port.species.count(value.first) || !species_bindings_.count(value.first))
				throw std::runtime_error("3D flow/transport input species '"
					+value.first+"' is not declared on its port");
	}

	static void MergeInput(PortBoundaryData& destination,
		const PortBoundaryData& source)
	{
		const auto scale = std::max({1.0, std::abs(destination.time_s),
			std::abs(source.time_s)});
		if (std::abs(destination.time_s-source.time_s) > 1.0e-12*scale)
			throw std::runtime_error("3D flow/transport input updates have inconsistent times");
		destination.time_s = source.time_s;
		auto merge_optional = [](auto& target, const auto& value) {
			if (!value) return;
			if (target) throw std::runtime_error(
				"3D flow/transport input quantity was supplied more than once");
			target = value;
		};
		merge_optional(destination.outward_flow_m3_s, source.outward_flow_m3_s);
		merge_optional(destination.mean_pressure_pa, source.mean_pressure_pa);
		merge_optional(destination.mean_normal_traction_pa, source.mean_normal_traction_pa);
		merge_optional(destination.total_pressure_pa, source.total_pressure_pa);
		for (const auto& value : source.concentration)
			if (!destination.concentration.emplace(value).second)
				throw std::runtime_error(
					"3D flow/transport concentration was supplied more than once");
		for (const auto& value : source.outward_species_flux)
			if (!destination.outward_species_flux.emplace(value).second)
				throw std::runtime_error(
					"3D flow/transport species flux was supplied more than once");
	}

	std::string domain_id_;
	TransientFlowRuntime& flow_runtime_;
	TransientTransportRuntime& transport_runtime_;
	std::vector<CouplingPort> ports_;
	SimulationConfiguration base_configuration_;
	std::filesystem::path case_directory_;
	std::map<std::string, std::string> species_bindings_;
	ThreeDFlowTransportDomainControls controls_;
	ThreeDBodyFittedFlowDomainAdapter flow_adapter_;
	DomainStepContext step_;
	std::map<std::string, PortBoundaryData> inputs_;
	bool trial_solve_succeeded_ = false;
};

} // namespace iga

#endif

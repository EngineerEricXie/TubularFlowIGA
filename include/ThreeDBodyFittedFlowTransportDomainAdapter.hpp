#ifndef IGA_THREE_D_BODY_FITTED_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP
#define IGA_THREE_D_BODY_FITTED_FLOW_TRANSPORT_DOMAIN_ADAPTER_HPP

#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "TransientTransportRuntime.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iterator>
#include <iomanip>
#include <sstream>
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

class ThreeDBodyFittedFlowTransportDomainAdapter : public CoupledDomainRuntime,
	public StagedFlowTransportDomainRuntime {
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
		std::string plan;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged begin preparation", [&] {
			step.Validate();
			if (transport_runtime_.Phase() != TransportStepPhase::Committed)
				throw std::runtime_error("3D staged BeginStep requires committed transport");
			plan = PlanSignature(step);
			const auto dt_scale = std::max({1.0, std::abs(step.dt_s),
				std::abs(transport_runtime_.System().dt)});
			if (std::abs(step.dt_s-transport_runtime_.System().dt) > 1.0e-12*dt_scale)
				throw std::runtime_error(
					"3D staged macro timestep must match the compiled transport timestep");
		});
		RequireCollectiveSameText(flow_runtime_.Communicator(), "3d staged plan agreement", plan);
		std::map<std::string, double> initial_mass;
		flow_adapter_.BeginStep(step);
		try {
			transport_runtime_.BeginStep();
			initial_mass = LogicalMass(transport_runtime_.TotalMass());
		} catch (...) {
			AbortStep();
			throw;
		}
		step_ = step;
		hydraulic_inputs_.clear();
		concentration_inputs_.clear();
		initial_species_mass_.swap(initial_mass);
		hydraulic_trial_succeeded_ = false;
		transport_trial_succeeded_ = false;
	}

	void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) override
	{
		decltype(hydraulic_inputs_) hydraulic_candidate;
		decltype(concentration_inputs_) concentration_candidate;
		std::string signature;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged input", [&] {
			hydraulic_candidate = hydraulic_inputs_;
			concentration_candidate = concentration_inputs_;
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
			const auto hydraulic = HydraulicInput(input);
			if (HydraulicValueCount(hydraulic) > 0) {
				auto found = hydraulic_candidate.find(port_id);
				if (found == hydraulic_candidate.end()) hydraulic_candidate.emplace(port_id, hydraulic);
				else MergeInput(found->second, hydraulic);
			}
			if (!input.concentration.empty()) {
				if (concentration_candidate.count(port_id))
					throw std::runtime_error(
						"3D flow/transport concentration was supplied more than once");
				concentration_candidate.emplace(port_id, input.concentration);
			}
			signature = InputSignature(port_id, input);
		});
		RequireCollectiveSameText(flow_runtime_.Communicator(), "3d staged input agreement", signature);
		hydraulic_inputs_.swap(hydraulic_candidate);
		concentration_inputs_.swap(concentration_candidate);
	}

	void SolveTrial() override
	{
		SolveHydraulicTrial();
		SolveTransportTrial();
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		return GetTransportPortState(port_id);
	}

	void RollbackTrial() override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged rollback preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved)
				throw std::runtime_error(
					"3D flow/transport rollback requires two solved trial phases");
		});
		RollbackTransportTrial();
		RollbackHydraulicTrial();
	}

	void AbortStep() override
	{
		std::exception_ptr failure;
		try { transport_runtime_.AbortStep(); }
		catch (...) { failure = std::current_exception(); }
		try { flow_adapter_.AbortStep(); }
		catch (...) { if (!failure) failure = std::current_exception(); }
		hydraulic_inputs_.clear();
		concentration_inputs_.clear();
		initial_species_mass_.clear();
		hydraulic_trial_succeeded_ = false;
		transport_trial_succeeded_ = false;
		if (failure) std::rethrow_exception(failure);
	}

	void PrepareCommitStep() override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged prepare commit", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved
				|| !hydraulic_trial_succeeded_ || !transport_trial_succeeded_)
				throw std::runtime_error(
					"3D flow/transport prepare requires two successful trial solves");
		});
		transport_runtime_.PrepareCommitStep();
		flow_adapter_.PrepareCommitStep();
	}

	void FinalizeCommitStep() noexcept override
	{
		transport_runtime_.FinalizeCommitStep();
		flow_adapter_.FinalizeCommitStep();
		hydraulic_inputs_.clear();
		concentration_inputs_.clear();
		initial_species_mass_.clear();
		hydraulic_trial_succeeded_ = false;
		transport_trial_succeeded_ = false;
	}

	void SolveHydraulicTrial() override
	{
		std::vector<std::pair<std::string, PortBoundaryData>> pending;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged hydraulic preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialReady
				|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen)
				throw std::runtime_error(
					"3D staged hydraulic solve requires flow trial-ready and transport trial-open");
			for (const auto& port : ports_) {
				const auto found = hydraulic_inputs_.find(port.id);
				const int supplied = found == hydraulic_inputs_.end()
					? 0 : HydraulicValueCount(found->second);
				if (HydraulicRequirementCount(port) != supplied)
					throw std::runtime_error("3D flow/transport domain is missing or has an extra hydraulic input for port '"
						+port.id+"'");
				if (supplied == 1) pending.emplace_back(port.id, found->second);
			}
		});
		hydraulic_trial_succeeded_ = false;
		for (const auto& input : pending) flow_adapter_.SetPortInput(input.first, input.second);
		flow_adapter_.SolveTrial();
		hydraulic_trial_succeeded_ = true;
	}

	PortState GetHydraulicPortState(const std::string& port_id) const override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged hydraulic port preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| !hydraulic_trial_succeeded_)
				throw std::runtime_error(
					"3D staged hydraulic port state requires a successful hydraulic trial");
		});
		return flow_adapter_.GetPortState(port_id);
	}

	void RollbackHydraulicTrial() override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged hydraulic rollback preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen)
				throw std::runtime_error(
					"3D staged hydraulic rollback requires solved flow and open transport");
		});
		// The flow-only adapter deliberately retains its legacy rollback behavior.
		// Reopen only the flow transaction so its compatibility input cache is
		// fresh, while the independently open transport transaction is untouched.
		flow_adapter_.AbortStep();
		flow_adapter_.BeginStep(step_);
		hydraulic_inputs_.clear();
		hydraulic_trial_succeeded_ = false;
	}

	void SetTransportConcentration(const std::string& port_id, double time_s,
		const std::map<std::string, double>& concentration) override
	{
		decltype(concentration_inputs_) candidate;
		std::string signature;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged concentration", [&] {
			candidate = concentration_inputs_;
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen
				|| !hydraulic_trial_succeeded_)
				throw std::runtime_error(
					"3D staged concentration input requires accepted hydraulic trial and open transport");
			const auto& port = Port(port_id);
			if (!port.requires.count(PortQuantity::SpeciesConcentration)
				|| port.species.empty() || SpeciesKeys(concentration) != port.species)
				throw std::runtime_error(
					"3D staged concentration input requires the complete declared logical species map");
			for (const auto& value : concentration) RequireFinitePortValue("3D staged concentration", value.second);
			ValidateTransportInputTime(time_s);
			if (candidate.count(port_id))
				throw std::runtime_error("3D staged concentration was supplied more than once");
			candidate.emplace(port_id, concentration);
			PortBoundaryData input;
			input.time_s = time_s; input.concentration = concentration;
			signature = InputSignature(port_id, input);
		});
		RequireCollectiveSameText(flow_runtime_.Communicator(), "3d staged concentration agreement", signature);
		concentration_inputs_.swap(candidate);
	}

	void SolveTransportTrial() override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged transport solve preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialOpen
				|| !hydraulic_trial_succeeded_)
				throw std::runtime_error(
					"3D staged transport solve requires a successful hydraulic trial and open transport");
		});
		transport_trial_succeeded_ = false;
		SimulationConfiguration trial_configuration;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged transport configuration", [&] {
			trial_configuration = MaterializeBoundaryWaveforms(base_configuration_, case_directory_, step_.EndTime());
		});
		for (const auto& port : ports_) ConfigureTransportPort(trial_configuration, port);
		transport_runtime_.SolveTrial(trial_configuration,
			flow_runtime_.RequiredNodes(), flow_runtime_.GatherRequiredVelocity());
		transport_trial_succeeded_ = true;
	}

	PortState GetTransportPortState(const std::string& port_id) const override
	{
		const CouplingPort* port = nullptr;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged transport port preparation", [&] {
			port = &Port(port_id);
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved
				|| !hydraulic_trial_succeeded_ || !transport_trial_succeeded_)
				throw std::runtime_error(
					"3D staged transport port state requires successful hydraulic and transport trials");
		});
		return LogicalTransportPortState(*port);
	}

	void RollbackTransportTrial() override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged transport rollback preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved)
				throw std::runtime_error(
					"3D staged transport rollback requires accepted flow and solved transport");
		});
		transport_runtime_.RollbackTrial();
		concentration_inputs_.clear();
		transport_trial_succeeded_ = false;
	}

	std::map<std::string, SpeciesStepAccounting>
	GetSpeciesStepAccounting() const override
	{
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d accounting preparation", [&] {
			if (flow_runtime_.Phase() != FlowStepPhase::TrialSolved
				|| transport_runtime_.Phase() != TransportStepPhase::TrialSolved
				|| !hydraulic_trial_succeeded_ || !transport_trial_succeeded_)
				throw std::runtime_error(
					"3D species accounting requires successful hydraulic and transport trials");
		});
		const auto final_mass = LogicalMass(transport_runtime_.TotalMass());
		const auto source_rate = LogicalSourceRate();
		std::map<std::string, PortState> states;
		for (const auto& port : ports_) {
			if (port.species.empty()) continue;
			auto state = GetTransportPortState(port.id);
			CollectiveLocalStage(flow_runtime_.Communicator(), "3d accounting port cache", [&] {
				states.emplace(port.id, std::move(state));
			});
		}
		std::map<std::string, SpeciesStepAccounting> result;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d accounting result", [&] {
			for (const auto& binding : species_bindings_) {
				SpeciesStepAccounting accounting;
				accounting.initial_mass = initial_species_mass_.at(binding.first);
				accounting.final_mass = final_mass.at(binding.first);
				// Backward Euler reports physical end-step compiled flux/source rates.
				// The amount is dt times that rate; residual is diagnostic only because
				// stabilized SUPG algebra is not claimed to be exactly conservative.
				accounting.source_amount = step_.dt_s*source_rate.at(binding.first);
				for (const auto& port : ports_) if (port.species.count(binding.first))
					accounting.outward_port_amount.emplace(port.id,
						step_.dt_s*states.at(port.id).outward_species_flux.at(binding.first));
				for (const auto& amount : accounting.outward_port_amount)
					accounting.residual += amount.second;
				accounting.residual += accounting.final_mass-accounting.initial_mass-accounting.source_amount;
				result.emplace(binding.first, std::move(accounting));
			}
		});
		return result;
	}

private:
	static std::string InputSignature(const std::string& port_id, const PortBoundaryData& input)
	{
		std::ostringstream text;
		text.exceptions(std::ios::badbit | std::ios::failbit);
		text << std::quoted(port_id) << ' ' << std::hexfloat << input.time_s;
		for (const auto& value : {input.outward_flow_m3_s, input.mean_pressure_pa,
			input.mean_normal_traction_pa, input.total_pressure_pa}) {
			text << ' ' << value.has_value();
			if (value) text << ' ' << *value;
		}
		for (const auto* values : {&input.concentration, &input.outward_species_flux}) {
			text << ' ' << values->size();
			for (const auto& value : *values) text << ' ' << std::quoted(value.first) << ' ' << value.second;
		}
		return text.str();
	}

	std::string PlanSignature(const DomainStepContext& step) const
	{
		std::ostringstream text;
		text.exceptions(std::ios::badbit | std::ios::failbit);
		text << step.step_index << ' ' << std::hexfloat << step.start_time_s << ' ' << step.dt_s
			<< ' ' << controls_.species_flow_epsilon_m3_s << ' ' << ports_.size();
		for (const auto& port : ports_) {
			text << ' ' << std::quoted(port.id) << ' ' << ParseThreeDFlowBoundaryLabel(port)
				<< ' ' << port.orientation.native_to_outward_sign;
			for (const auto* quantities : {&port.requires, &port.provides}) {
				text << ' ' << quantities->size();
				for (const auto value : *quantities) text << ' ' << static_cast<int>(value);
			}
			text << ' ' << port.species.size();
			for (const auto& field : port.species) text << ' ' << std::quoted(field);
		}
		text << ' ' << species_bindings_.size();
		for (const auto& binding : species_bindings_)
			text << ' ' << std::quoted(binding.first) << ' ' << std::quoted(binding.second);
		return text.str();
	}

	void ValidateTransportInputTime(double time_s) const
	{
		if (!std::isfinite(time_s))
			throw std::runtime_error("3D staged concentration input time must be finite");
		const auto scale = std::max({1.0, std::abs(step_.EndTime()), std::abs(time_s)});
		if (std::abs(time_s-step_.EndTime()) > 1.0e-12*scale)
			throw std::runtime_error(
				"3D staged concentration input time does not match the active step");
	}

	void ConfigureTransportPort(SimulationConfiguration& configuration,
		const CouplingPort& port) const
	{
		if (port.species.empty()
			|| !port.requires.count(PortQuantity::SpeciesConcentration)) return;
		const auto found = concentration_inputs_.find(port.id);
		const bool has_concentration = found != concentration_inputs_.end();
		const auto flow = GetHydraulicPortState(port.id);
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d staged species boundary", [&] {
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
				if (has_concentration) value = found->second.at(species);
				SetThreeDPortSpeciesBoundary(configuration, port,
					species_bindings_.at(species), value);
			}
		});
	}

	PortState LogicalTransportPortState(const CouplingPort& port) const
	{
		std::vector<CouplingPort> ports;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d logical port preparation", [&] { ports.push_back(port); });
		const auto values = transport_runtime_.GatherRequiredState();
		const auto measured = flow_runtime_.MeasurePorts(ports, step_.EndTime(),
			transport_runtime_.System().fields, values, &transport_runtime_.System());
		PortState state;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d logical port result", [&] {
			state = measured.at(port.id);
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
		});
		return state;
	}

	std::map<std::string, double> LogicalMass(
		const std::map<std::string, double>& native_mass) const
	{
		std::map<std::string, double> result;
		CollectiveLocalStage(flow_runtime_.Communicator(), "3d logical mass", [&] {
			for (const auto& binding : species_bindings_)
				result.emplace(binding.first, native_mass.at(binding.second));
		});
		return result;
	}

	std::map<std::string, double> LogicalSourceRate() const
	{
		return LogicalMass(transport_runtime_.SourceIntegrals());
	}

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
	std::map<std::string, PortBoundaryData> hydraulic_inputs_;
	std::map<std::string, std::map<std::string, double>> concentration_inputs_;
	std::map<std::string, double> initial_species_mass_;
	bool hydraulic_trial_succeeded_ = false;
	bool transport_trial_succeeded_ = false;
};

} // namespace iga

#endif

#ifndef IGA_ZERO_D_FLOW_DOMAIN_HPP
#define IGA_ZERO_D_FLOW_DOMAIN_HPP

// Dependency-free 0D hydraulic contracts.  Flow is always positive outward
// from this 0D subsystem, and all quantities use SI units.
#include "CoupledDomainRuntime.hpp"
#include "Sha256.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

enum class ZeroDFlowRole : std::uint8_t {
	SourceReservoir,
	TerminalRcr
};

inline const char* ZeroDFlowRoleName(ZeroDFlowRole role)
{
	if (role == ZeroDFlowRole::SourceReservoir) return "source_reservoir";
	if (role == ZeroDFlowRole::TerminalRcr) return "terminal_rcr";
	return "unknown";
}

inline bool IsKnownZeroDFlowRole(ZeroDFlowRole role)
{
	return role == ZeroDFlowRole::SourceReservoir || role == ZeroDFlowRole::TerminalRcr;
}

struct ZeroDSourceReservoirModel {
	double capacitance_m3_pa = 0.0;
	double resistance_pa_s_m3 = 0.0;
	double prescribed_flow_m3_s = 0.0;
};

struct ZeroDTerminalRcrModel {
	double proximal_resistance_pa_s_m3 = 0.0;
	double distal_resistance_pa_s_m3 = 0.0;
	double capacitance_m3_pa = 0.0;
	double distal_pressure_pa = 0.0;
};

struct ZeroDFlowModel {
	ZeroDFlowRole role = ZeroDFlowRole::SourceReservoir;
	ZeroDSourceReservoirModel source;
	ZeroDTerminalRcrModel terminal;
};

struct ZeroDFlowState {
	double stored_pressure_pa = 0.0;
};

struct ZeroDFlowModelConfiguration {
	ZeroDFlowModel model;
	ZeroDFlowState initial_state;
};

struct ZeroDFlowStorageBalance {
	double stored_volume_change_m3 = 0.0;
	double prescribed_source_amount_m3 = 0.0;
	double distal_sink_amount_m3 = 0.0;
	double outward_graph_port_amount_m3 = 0.0;
	double residual_m3 = 0.0;
};

// Per-step volume accounting is intentionally separate from the hydraulic
// state.  It gives callers a compact conservation check without adding a
// second numerical update path to either 0D model.
struct ZeroDFlowStepAccounting {
	double initial_stored_volume_m3 = 0.0;
	double final_stored_volume_m3 = 0.0;
	double prescribed_source_amount_m3 = 0.0;
	double distal_sink_amount_m3 = 0.0;
	double outward_graph_port_amount_m3 = 0.0;
	double residual_m3 = 0.0;
};

struct ZeroDFlowTrial {
	ZeroDFlowState state;
	PortState port;
	ZeroDFlowStorageBalance storage;
};

// Complete accepted publication. Initial (pre-step) state and active trials
// deliberately have no checkpoint representation in the v1 graph contract.
struct ZeroDFlowCheckpointState {
	std::string domain_id;
	std::string model_identity_sha256;
	DomainStepContext accepted_step;
	ZeroDFlowState state;
	PortState port;
	ZeroDFlowStepAccounting accounting;
};

namespace zero_d_flow_detail {

inline void AppendString(Sha256& hash, const std::string& value)
{
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(value.size()));
	hash.Append(value.data(), value.size());
}

inline void RequirePositiveFinite(const std::string& name, double value)
{
	if (!(value > 0.0) || !std::isfinite(value))
		throw std::runtime_error(name+" must be positive and finite");
}

inline void RequireNonnegativeFinite(const std::string& name, double value)
{
	if (!(value >= 0.0) || !std::isfinite(value))
		throw std::runtime_error(name+" must be nonnegative and finite");
}

} // namespace zero_d_flow_detail

inline void ValidateZeroDFlowModel(const ZeroDFlowModel& model)
{
	if (!IsKnownZeroDFlowRole(model.role))
		throw std::runtime_error("0D flow model has an unknown role");
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		zero_d_flow_detail::RequirePositiveFinite("0D source capacitance_m3_pa",
			model.source.capacitance_m3_pa);
		zero_d_flow_detail::RequirePositiveFinite("0D source resistance_pa_s_m3",
			model.source.resistance_pa_s_m3);
		RequireFinitePortValue("0D source prescribed_flow_m3_s",
			model.source.prescribed_flow_m3_s);
	} else {
		zero_d_flow_detail::RequireNonnegativeFinite("0D terminal proximal_resistance_pa_s_m3",
			model.terminal.proximal_resistance_pa_s_m3);
		zero_d_flow_detail::RequirePositiveFinite("0D terminal distal_resistance_pa_s_m3",
			model.terminal.distal_resistance_pa_s_m3);
		zero_d_flow_detail::RequirePositiveFinite("0D terminal capacitance_m3_pa",
			model.terminal.capacitance_m3_pa);
		RequireFinitePortValue("0D terminal distal_pressure_pa", model.terminal.distal_pressure_pa);
	}
}

inline void ValidateZeroDFlowState(const ZeroDFlowState& state)
{
	RequireFinitePortValue("0D flow stored_pressure_pa", state.stored_pressure_pa);
}

inline std::string BuildZeroDFlowModelIdentitySha256(const ZeroDFlowModel& model)
{
	ValidateZeroDFlowModel(model);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash, "ZeroDFlowModel/v1");
	hash.AppendLittleEndian32(static_cast<std::uint32_t>(model.role));
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		hash.AppendNormalizedDouble(model.source.capacitance_m3_pa);
		hash.AppendNormalizedDouble(model.source.resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.source.prescribed_flow_m3_s);
	} else {
		hash.AppendNormalizedDouble(model.terminal.proximal_resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.terminal.distal_resistance_pa_s_m3);
		hash.AppendNormalizedDouble(model.terminal.capacitance_m3_pa);
		hash.AppendNormalizedDouble(model.terminal.distal_pressure_pa);
	}
	return hash.Hex();
}

inline std::string BuildZeroDFlowStateIdentitySha256(const ZeroDFlowModel& model,
	const ZeroDFlowState& state)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(state);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash, "ZeroDFlowState/v1");
	zero_d_flow_detail::AppendString(hash, BuildZeroDFlowModelIdentitySha256(model));
	hash.AppendNormalizedDouble(state.stored_pressure_pa);
	return hash.Hex();
}

inline void ValidateZeroDFlowStepAccounting(const ZeroDFlowStepAccounting& accounting)
{
	RequireFinitePortValue("0D flow initial stored volume", accounting.initial_stored_volume_m3);
	RequireFinitePortValue("0D flow final stored volume", accounting.final_stored_volume_m3);
	RequireFinitePortValue("0D flow prescribed source amount", accounting.prescribed_source_amount_m3);
	RequireFinitePortValue("0D flow distal sink amount", accounting.distal_sink_amount_m3);
	RequireFinitePortValue("0D flow outward graph-port amount", accounting.outward_graph_port_amount_m3);
	RequireFinitePortValue("0D flow storage residual", accounting.residual_m3);
}

inline std::string BuildZeroDFlowStepAccountingIdentitySha256(
	const ZeroDFlowModel& model, const ZeroDFlowStepAccounting& accounting)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowStepAccounting(accounting);
	Sha256 hash;
	zero_d_flow_detail::AppendString(hash, "ZeroDFlowStepAccounting/v1");
	zero_d_flow_detail::AppendString(hash, BuildZeroDFlowModelIdentitySha256(model));
	hash.AppendNormalizedDouble(accounting.initial_stored_volume_m3);
	hash.AppendNormalizedDouble(accounting.final_stored_volume_m3);
	hash.AppendNormalizedDouble(accounting.prescribed_source_amount_m3);
	hash.AppendNormalizedDouble(accounting.distal_sink_amount_m3);
	hash.AppendNormalizedDouble(accounting.outward_graph_port_amount_m3);
	hash.AppendNormalizedDouble(accounting.residual_m3);
	return hash.Hex();
}

inline CouplingPort MakeZeroDFlowPort(const std::string& domain_id, ZeroDFlowRole role)
{
	if (domain_id.empty()) throw std::runtime_error("0D flow domain id must be nonempty");
	if (!IsKnownZeroDFlowRole(role)) throw std::runtime_error("0D flow role is unknown");
	CouplingPort port;
	port.id = "port";
	port.subsystem_id = domain_id;
	port.locator_kind = "zero_d_port";
	port.locator = "port";
	port.provides = {PortQuantity::MeanPressure, PortQuantity::FlowRate};
	port.requires = role == ZeroDFlowRole::SourceReservoir
		? std::set<PortQuantity>{PortQuantity::MeanPressure}
		: std::set<PortQuantity>{PortQuantity::FlowRate};
	return port;
}

inline void ValidateZeroDFlowDomainMetadata(const std::string& domain_id,
	const std::vector<CouplingPort>& ports, ZeroDFlowRole role)
{
	if (ports.size() != 1) throw std::runtime_error("0D flow domain requires exactly one port");
	const auto expected = MakeZeroDFlowPort(domain_id, role);
	const auto& port = ports.front();
	ValidateCouplingPort(port);
	if (port.id != expected.id || port.subsystem_id != expected.subsystem_id
		|| port.locator_kind != expected.locator_kind || port.locator != expected.locator
		|| port.orientation.native_to_outward_sign != 1 || port.provides != expected.provides
		|| port.requires != expected.requires || !port.species.empty())
		throw std::runtime_error("0D flow domain port does not match its model role");
}

inline ZeroDFlowTrial EvaluateZeroDFlowTrial(const ZeroDFlowModel& model,
	const ZeroDFlowState& committed, double port_input, double dt_s, double time_s = 0.0)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(committed);
	zero_d_flow_detail::RequirePositiveFinite("0D flow dt_s", dt_s);
	RequireFinitePortValue("0D flow port input", port_input);
	RequireFinitePortValue("0D flow time_s", time_s);
	ZeroDFlowTrial result;
	result.port.time_s = time_s+dt_s;
	if (model.role == ZeroDFlowRole::SourceReservoir) {
		const auto& source = model.source;
		const double pressure = ((source.capacitance_m3_pa/dt_s)*committed.stored_pressure_pa
			+source.prescribed_flow_m3_s+port_input/source.resistance_pa_s_m3)
			/((source.capacitance_m3_pa/dt_s)+1.0/source.resistance_pa_s_m3);
		const double outward_flow = (pressure-port_input)/source.resistance_pa_s_m3;
		result.state.stored_pressure_pa = pressure;
		// The source consumes the interface pressure; it is reported unchanged
		// because this one-port contract has no fabricated internal boundary area.
		result.port.mean_pressure_pa = port_input;
		result.port.outward_flow_m3_s = outward_flow;
		result.storage.stored_volume_change_m3 = source.capacitance_m3_pa
			*(pressure-committed.stored_pressure_pa);
		result.storage.prescribed_source_amount_m3 = source.prescribed_flow_m3_s*dt_s;
		result.storage.outward_graph_port_amount_m3 = outward_flow*dt_s;
		result.storage.residual_m3 = result.storage.stored_volume_change_m3
			-result.storage.prescribed_source_amount_m3
			+result.storage.outward_graph_port_amount_m3;
	} else {
		const auto& terminal = model.terminal;
		// Input uses the 0D port convention: Q_port=-Qin.
		const double inlet_flow = -port_input;
		const double pressure = ((terminal.capacitance_m3_pa/dt_s)
			*committed.stored_pressure_pa+inlet_flow
			+terminal.distal_pressure_pa/terminal.distal_resistance_pa_s_m3)
			/((terminal.capacitance_m3_pa/dt_s)+1.0/terminal.distal_resistance_pa_s_m3);
		const double distal_flow = (pressure-terminal.distal_pressure_pa)
			/terminal.distal_resistance_pa_s_m3;
		result.state.stored_pressure_pa = pressure;
		result.port.mean_pressure_pa = pressure
			+terminal.proximal_resistance_pa_s_m3*inlet_flow;
		result.port.outward_flow_m3_s = port_input;
		result.storage.stored_volume_change_m3 = terminal.capacitance_m3_pa
			*(pressure-committed.stored_pressure_pa);
		result.storage.distal_sink_amount_m3 = distal_flow*dt_s;
		result.storage.outward_graph_port_amount_m3 = port_input*dt_s;
		result.storage.residual_m3 = result.storage.stored_volume_change_m3
			+result.storage.distal_sink_amount_m3
			+result.storage.outward_graph_port_amount_m3;
	}
	ValidateZeroDFlowState(result.state);
	ValidatePortState(result.port);
	return result;
}

inline ZeroDFlowStepAccounting MakeZeroDFlowStepAccounting(const ZeroDFlowModel& model,
	const ZeroDFlowState& initial, const ZeroDFlowTrial& trial)
{
	ValidateZeroDFlowModel(model);
	ValidateZeroDFlowState(initial);
	ValidateZeroDFlowState(trial.state);
	const double capacitance = model.role == ZeroDFlowRole::SourceReservoir
		? model.source.capacitance_m3_pa : model.terminal.capacitance_m3_pa;
	ZeroDFlowStepAccounting result;
	result.initial_stored_volume_m3 = capacitance*initial.stored_pressure_pa;
	result.final_stored_volume_m3 = capacitance*trial.state.stored_pressure_pa;
	result.prescribed_source_amount_m3 = trial.storage.prescribed_source_amount_m3;
	result.distal_sink_amount_m3 = trial.storage.distal_sink_amount_m3;
	result.outward_graph_port_amount_m3 = trial.storage.outward_graph_port_amount_m3;
	result.residual_m3 = trial.storage.residual_m3;
	ValidateZeroDFlowStepAccounting(result);
	return result;
}

// Transactional production owner for the dependency-free 0D kernels.  The
// only mutable numerical state is committed pressure; every SolveTrial call
// evaluates from the step's committed t_n snapshot, never from a prior trial.
class ZeroDFlowDomainRuntime final : public CoupledDomainRuntime {
public:
	ZeroDFlowDomainRuntime(std::string domain_id, ZeroDFlowModel model,
		ZeroDFlowState initial_state, std::vector<CouplingPort> ports,
		double initial_time_s = 0.0)
		: domain_id_(std::move(domain_id)), model_(std::move(model)),
		  ports_(std::move(ports)), committed_state_(initial_state)
	{
		if (domain_id_.empty()) throw std::runtime_error("0D flow runtime domain id must be nonempty");
		if (!std::isfinite(initial_time_s))
			throw std::runtime_error("0D flow runtime initial time must be finite");
		ValidateZeroDFlowModel(model_);
		ValidateZeroDFlowState(committed_state_);
		ValidateZeroDFlowDomainMetadata(domain_id_, ports_, model_.role);
		committed_time_s_ = initial_time_s;
	}

	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ZeroDFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }
	const ZeroDFlowModel& Model() const noexcept { return model_; }
	const ZeroDFlowState& CommittedState() const noexcept { return committed_state_; }
	double CommittedTime() const noexcept { return committed_time_s_; }
	int CommittedStepIndex() const noexcept { return committed_step_index_; }
	std::size_t CommittedStepCount() const noexcept { return committed_step_count_; }
	const std::optional<PortState>& CommittedPortState() const noexcept { return committed_port_; }
	const std::optional<ZeroDFlowStepAccounting>& CommittedStepAccounting() const noexcept
	{ return committed_accounting_; }
	std::string ModelIdentitySha256() const { return BuildZeroDFlowModelIdentitySha256(model_); }
	std::string CommittedStateIdentitySha256() const
	{ return BuildZeroDFlowStateIdentitySha256(model_, committed_state_); }

	ZeroDFlowCheckpointState CaptureCheckpointState() const
	{
		RequirePhase(Phase::Committed, "capture checkpoint");
		if (!committed_step_count_ || !committed_port_ || !committed_accounting_)
			throw std::runtime_error("0D checkpoint requires an accepted step");
		ZeroDFlowCheckpointState result{domain_id_, ModelIdentitySha256(), committed_step_,
			committed_state_, *committed_port_, *committed_accounting_};
		ValidateCheckpointState(result);
		return result;
	}

	// Local candidate initialization. The graph coordinator must restore all
	// domains into unpublished owners, agree on success, then publish together.
	// Validation/copies finish before the first accepted datum is changed.
	void RestoreCheckpointState(ZeroDFlowCheckpointState value)
	{
		RequirePhase(Phase::Committed, "restore checkpoint");
		if (committed_step_count_)
			throw std::runtime_error("0D checkpoint restore requires a fresh runtime");
		ValidateCheckpointState(value);
		Staged restored{value.state, std::move(value.port), value.accounting};
		static_assert(noexcept(committed_port_.swap(restored.port)) && noexcept(committed_accounting_.swap(restored.accounting)),
			"validated 0D checkpoint publication must not throw");
		committed_port_.swap(restored.port); committed_accounting_.swap(restored.accounting);
		committed_state_ = restored.state;
		committed_step_ = value.accepted_step;
		committed_time_s_ = committed_step_.EndTime();
		committed_step_index_ = committed_step_.step_index;
		committed_step_count_ = static_cast<std::size_t>(committed_step_index_)+1;
	}

	void BeginStep(const DomainStepContext& step) override
	{
		RequirePhase(Phase::Committed, "begin step");
		step.Validate();
		if (step.start_time_s != committed_time_s_)
			throw std::runtime_error("0D flow runtime step start time is not the committed time");
		if (committed_step_index_ == std::numeric_limits<int>::max()
			|| step.step_index != committed_step_index_+1)
			throw std::runtime_error("0D flow runtime step index is not the next committed index");
		step_ = step;
		base_state_ = committed_state_;
		input_.reset(); trial_.reset(); trial_accounting_.reset(); staged_.reset();
		phase_ = Phase::TrialReady;
	}

	void SetPortInput(const std::string& port_id, const PortBoundaryData& input) override
	{
		RequirePhase(Phase::TrialReady, "set port input");
		if (port_id != ports_.front().id)
			throw std::runtime_error("0D flow runtime has no port '"+port_id+"'");
		ValidatePortBoundaryData(input);
		if (input.time_s != step_.EndTime())
			throw std::runtime_error("0D flow runtime input time does not match the active step end");
		if (input_)
			throw std::runtime_error("0D flow runtime port input was supplied twice");
		const bool source = model_.role == ZeroDFlowRole::SourceReservoir;
		if (source) {
			if (!input.mean_pressure_pa || input.outward_flow_m3_s || input.mean_normal_traction_pa
				|| input.total_pressure_pa || !input.concentration.empty() || !input.outward_species_flux.empty())
				throw std::runtime_error("0D source runtime input requires exactly mean_pressure_pa");
		} else if (!input.outward_flow_m3_s || input.mean_pressure_pa || input.mean_normal_traction_pa
			|| input.total_pressure_pa || !input.concentration.empty() || !input.outward_species_flux.empty()) {
			throw std::runtime_error("0D terminal runtime input requires exactly outward_flow_m3_s");
		}
		input_ = input;
	}

	void SolveTrial() override
	{
		RequirePhase(Phase::TrialReady, "solve trial");
		if (!input_) throw std::runtime_error("0D flow runtime is missing its required port input");
		const double value = model_.role == ZeroDFlowRole::SourceReservoir
			? *input_->mean_pressure_pa : *input_->outward_flow_m3_s;
		// Evaluate from base_state_ every time; that snapshot is made exactly once
		// at BeginStep and is untouched by rollback/retry activity.
		auto trial = EvaluateZeroDFlowTrial(model_, base_state_, value, step_.dt_s,
			step_.start_time_s);
		auto accounting = MakeZeroDFlowStepAccounting(model_, base_state_, trial);
		ValidatePortState(trial.port);
		ValidateZeroDFlowStepAccounting(accounting);
		trial_ = std::move(trial);
		trial_accounting_ = std::move(accounting);
		phase_ = Phase::TrialSolved;
	}

	PortState GetPortState(const std::string& port_id) const override
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::Prepared)
			throw std::runtime_error("0D flow runtime has no solved trial port state");
		if (port_id != ports_.front().id)
			throw std::runtime_error("0D flow runtime has no port '"+port_id+"'");
		return trial_->port;
	}

	// Return a snapshot so callers cannot retain a reference into trial storage
	// that is invalidated by rollback, abort, or final promotion.
	ZeroDFlowStepAccounting TrialStepAccounting() const
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::Prepared)
			throw std::runtime_error("0D flow runtime has no solved trial accounting");
		return *trial_accounting_;
	}

	void RollbackTrial() override
	{
		RequirePhase(Phase::TrialSolved, "rollback trial");
		input_.reset(); trial_.reset(); trial_accounting_.reset();
		phase_ = Phase::TrialReady;
	}

	void AbortStep() override
	{
		if (phase_ != Phase::TrialReady && phase_ != Phase::TrialSolved && phase_ != Phase::Prepared)
			throw std::runtime_error("0D flow runtime abort requires an active step");
		input_.reset(); trial_.reset(); trial_accounting_.reset(); staged_.reset();
		phase_ = Phase::Committed;
	}

	void PrepareCommitStep() override
	{
		RequirePhase(Phase::TrialSolved, "prepare commit");
		// Do all potentially throwing checks and copies before the no-throw
		// promotion path.  No committed datum is published here.
		ValidateZeroDFlowState(trial_->state);
		ValidatePortState(trial_->port);
		ValidateZeroDFlowStepAccounting(*trial_accounting_);
		Staged committed{trial_->state, trial_->port, *trial_accounting_};
		staged_ = std::move(committed);
		phase_ = Phase::Prepared;
	}

	void FinalizeCommitStep() noexcept override
	{
		if (phase_ != Phase::Prepared) return;
		// Only scalar/optional swaps and moves remain: this promotion cannot
		// allocate or call the numerical kernel.
		committed_state_ = staged_->state;
		committed_port_.swap(staged_->port);
		committed_accounting_.swap(staged_->accounting);
		committed_time_s_ = step_.EndTime();
		committed_step_ = step_;
		committed_step_index_ = step_.step_index;
		++committed_step_count_;
		input_.reset(); trial_.reset(); trial_accounting_.reset(); staged_.reset();
		phase_ = Phase::Committed;
	}

private:
	enum class Phase { Committed, TrialReady, TrialSolved, Prepared };

	struct Staged {
		ZeroDFlowState state;
		std::optional<PortState> port;
		std::optional<ZeroDFlowStepAccounting> accounting;
	};

	void ValidateCheckpointState(const ZeroDFlowCheckpointState& value) const
	{
		if (value.domain_id != domain_id_ || value.model_identity_sha256 != ModelIdentitySha256())
			throw std::runtime_error("0D checkpoint domain or model differs");
		value.accepted_step.Validate(); ValidateZeroDFlowState(value.state);
		ValidatePortState(value.port); ValidateZeroDFlowStepAccounting(value.accounting);
		const auto& port = value.port; const auto& accounting = value.accounting;
		if (port.time_s != value.accepted_step.EndTime() || !port.outward_flow_m3_s || !port.mean_pressure_pa
			|| port.area_m2 || port.mean_normal_traction_pa || port.total_pressure_pa
			|| !port.concentration.empty() || !port.outward_species_flux.empty())
			throw std::runtime_error("0D checkpoint accepted port is inconsistent");
		const double dt = value.accepted_step.dt_s, pressure = value.state.stored_pressure_pa;
		const double flow = *port.outward_flow_m3_s;
		const bool source = model_.role == ZeroDFlowRole::SourceReservoir;
		const double capacitance = source ? model_.source.capacitance_m3_pa : model_.terminal.capacitance_m3_pa;
		if (accounting.final_stored_volume_m3 != capacitance*pressure
			|| accounting.outward_graph_port_amount_m3 != flow*dt)
			throw std::runtime_error("0D checkpoint state and accounting differ");
		if (source) {
			if (flow != (pressure-*port.mean_pressure_pa)/model_.source.resistance_pa_s_m3
				|| accounting.prescribed_source_amount_m3 != model_.source.prescribed_flow_m3_s*dt
				|| accounting.distal_sink_amount_m3 != 0.0)
				throw std::runtime_error("0D source checkpoint is inconsistent with its model");
		} else {
			if (*port.mean_pressure_pa != pressure+model_.terminal.proximal_resistance_pa_s_m3*(-flow)
				|| accounting.distal_sink_amount_m3 != ((pressure-model_.terminal.distal_pressure_pa)
					/model_.terminal.distal_resistance_pa_s_m3)*dt
				|| accounting.prescribed_source_amount_m3 != 0.0)
				throw std::runtime_error("0D terminal checkpoint is inconsistent with its model");
		}
		// The kernel computes C*(p_new-p_old), whereas the stored accounting
		// exposes C*p_new and C*p_old. Bound their different rounding orders.
		const long double balance = static_cast<long double>(accounting.final_stored_volume_m3)
			-accounting.initial_stored_volume_m3-accounting.prescribed_source_amount_m3
			+accounting.distal_sink_amount_m3+accounting.outward_graph_port_amount_m3;
		long double scale = 0;
		for (double amount : {accounting.initial_stored_volume_m3, accounting.final_stored_volume_m3,
			accounting.prescribed_source_amount_m3, accounting.distal_sink_amount_m3,
			accounting.outward_graph_port_amount_m3}) scale += std::abs(static_cast<long double>(amount));
		const long double roundoff = 128*std::numeric_limits<double>::epsilon()*scale
			+16*static_cast<long double>(std::numeric_limits<double>::denorm_min());
		if (std::abs(balance-accounting.residual_m3) > roundoff || std::abs(balance) > roundoff)
			throw std::runtime_error("0D checkpoint storage balance is inconsistent");
	}

	void RequirePhase(Phase expected, const char* operation) const
	{
		if (phase_ != expected)
			throw std::runtime_error(std::string("0D flow runtime cannot ")+operation
				+" in its current lifecycle phase");
	}

	std::string domain_id_;
	const ZeroDFlowModel model_;
	const std::vector<CouplingPort> ports_;
	Phase phase_ = Phase::Committed;
	DomainStepContext step_;
	DomainStepContext committed_step_;
	ZeroDFlowState committed_state_;
	ZeroDFlowState base_state_;
	double committed_time_s_ = 0.0;
	int committed_step_index_ = -1;
	std::size_t committed_step_count_ = 0;
	std::optional<PortBoundaryData> input_;
	std::optional<ZeroDFlowTrial> trial_;
	std::optional<ZeroDFlowStepAccounting> trial_accounting_;
	std::optional<Staged> staged_;
	std::optional<PortState> committed_port_;
	std::optional<ZeroDFlowStepAccounting> committed_accounting_;
};

} // namespace iga

#endif

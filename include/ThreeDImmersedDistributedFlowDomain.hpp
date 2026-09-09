#ifndef IGA_THREE_D_IMMERSED_DISTRIBUTED_FLOW_DOMAIN_HPP
#define IGA_THREE_D_IMMERSED_DISTRIBUTED_FLOW_DOMAIN_HPP

#include "ThreeDImmersedFlowDomain.hpp"
#include "../solvers/cpu/include/ImmersedStaticDistributedRuntime.hpp"

namespace iga {

// Steady graph transactions on a distributed backend. Mutating operations are
// collective on runtime.Communicator(); port-state getters are local. Geometry
// remains fixed and coupling time never changes the steady dt=0 contract.
class ThreeDImmersedDistributedFlowDomain final : public CoupledDomainRuntime {
public:
	ThreeDImmersedDistributedFlowDomain(const std::string& domain_id,ImmersedStaticDistributedRuntime& runtime,
		const std::vector<CouplingPort>& ports,double initial_time_s = 0.0)
		: runtime_(runtime),communicator_(runtime.Communicator())
	{
		std::string signature;
		Local("distributed immersed domain metadata",[&] {
			if (!std::isfinite(initial_time_s)) throw std::invalid_argument("immersed initial time must be finite");
			ValidateThreeDImmersedFlowDomainMetadata(domain_id,ports);
			domain_id_ = domain_id; ports_ = ports;
			if (runtime_.PortDefinitions().size() != ports_.size()) throw std::invalid_argument("immersed runtime and graph port counts differ");
			for (const auto& definition : runtime_.PortDefinitions()) {
				const auto& port = Port(definition.id);
				if (ParseThreeDImmersedFlowBoundaryLabel(port) != definition.boundary_label
					|| port.requires != std::set<PortQuantity>{Quantity(definition.control_mode)})
					throw std::invalid_argument("immersed runtime and graph port binding differ");
				committed_controls_.emplace(definition.id,definition.value);
			}
			diagnostics_.committed_time_s = initial_time_s;
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10);
			Append(text,domain_id_); text << initial_time_s << ':' << ports_.size() << ':';
			for (const auto& port : ports_) {
				Append(text,port.id); Append(text,port.subsystem_id); Append(text,port.locator_kind); Append(text,port.locator);
				text << port.orientation.native_to_outward_sign << ':' << port.provides.size() << ':';
				for (auto q : port.provides) text << static_cast<int>(q) << ':';
				text << port.requires.size() << ':';
				for (auto q : port.requires) text << static_cast<int>(q) << ':';
			}
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_,"distributed immersed domain agreement",signature);
	}
	const std::string& DomainId() const noexcept override { return domain_id_; }
	DomainKind Kind() const noexcept override { return DomainKind::ThreeDImmersedFlow; }
	const std::vector<CouplingPort>& Ports() const noexcept override { return ports_; }
	const ThreeDImmersedFlowDomainDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	const std::map<std::string,PortState>& CommittedPortStates() const noexcept { return committed_measurements_; }
	PetscInt RowBegin() const noexcept { return runtime_.RowBegin(); }
	PetscInt RowEnd() const noexcept { return runtime_.RowEnd(); }
	std::vector<PetscScalar> CommittedOwnedBackendState() const
	{
		PetscReadArray view; view.Acquire(runtime_.CommittedState());
		std::vector<PetscScalar> result(static_cast<std::size_t>(RowEnd()-RowBegin()));
		std::copy_n(view.Data(),result.size(),result.data()); view.Restore(); return result;
	}
	void BeginStep(const DomainStepContext& step) override
	{
		std::string signature;
		Local("distributed immersed begin preflight",[&] {
			Require(Phase::Committed,"begin"); step.Validate();
			if (step.start_time_s != diagnostics_.committed_time_s
				|| diagnostics_.committed_step_index == std::numeric_limits<int>::max()
				|| step.step_index != diagnostics_.committed_step_index+1)
				throw std::invalid_argument("immersed step is not the next committed time/index");
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			text << std::setprecision(std::numeric_limits<double>::max_digits10) << step.start_time_s << ':' << step.dt_s << ':' << step.step_index;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_,"distributed immersed step agreement",signature);
		step_ = step; inputs_.clear(); trial_controls_.clear(); trial_measurements_.clear(); phase_ = Phase::TrialReady;
	}
	void SetPortInput(const std::string& id,const PortBoundaryData& input) override
	{
		std::map<std::string,double> candidate;
		std::string signature;
		Local("distributed immersed input preflight",[&] {
			Require(Phase::TrialReady,"set input"); ValidatePortBoundaryData(input);
			const auto& port = Port(id);
			const int count = static_cast<int>(input.outward_flow_m3_s.has_value())+static_cast<int>(input.mean_pressure_pa.has_value())+static_cast<int>(input.mean_normal_traction_pa.has_value());
			if (input.time_s != step_.EndTime() || inputs_.count(id) || count != 1 || input.total_pressure_pa
				|| !input.concentration.empty() || !input.outward_species_flux.empty()) throw std::invalid_argument("invalid immersed hydraulic input");
			const auto quantity = input.outward_flow_m3_s ? PortQuantity::FlowRate : input.mean_pressure_pa ? PortQuantity::MeanPressure : PortQuantity::MeanNormalTraction;
			if (port.requires != std::set<PortQuantity>{quantity}) throw std::invalid_argument("immersed input quantity differs from port metadata");
			const double value = input.outward_flow_m3_s ? *input.outward_flow_m3_s : input.mean_pressure_pa ? *input.mean_pressure_pa : *input.mean_normal_traction_pa;
			candidate = inputs_; candidate.emplace(id,value);
			std::ostringstream text; text.exceptions(std::ios::badbit | std::ios::failbit);
			Append(text,id); text << std::setprecision(std::numeric_limits<double>::max_digits10) << input.time_s << ':' << static_cast<int>(quantity) << ':' << value;
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_,"distributed immersed input agreement",signature);
		inputs_.swap(candidate);
	}
	void SolveTrial() override
	{
		std::map<std::string,double> controls;
		Local("distributed immersed solve preflight",[&] {
			Require(Phase::TrialReady,"solve");
			for (const auto& port : ports_) if (!inputs_.count(port.id)) throw std::invalid_argument("immersed trial is missing controlled input");
			controls = inputs_;
		});
		try {
			for (const auto& control : controls) runtime_.SetPortControlValue(control.first,control.second);
			if (!runtime_.SolveTrial()) throw std::runtime_error("immersed trial did not converge");
			std::map<std::string,PortState> measured;
			Local("distributed immersed trial measurements",[&] {
				for (const auto& port : ports_) measured.emplace(port.id,Measure(port));
			});
			trial_controls_.swap(controls); trial_measurements_.swap(measured); phase_ = Phase::TrialSolved;
		} catch (...) {
			if (runtime_.Diagnostics().trial_active) runtime_.Rollback();
			RestoreControls(); trial_controls_.clear(); trial_measurements_.clear(); phase_ = Phase::TrialReady; throw;
		}
	}
	PortState GetPortState(const std::string& id) const override
	{
		if (phase_ != Phase::TrialSolved && phase_ != Phase::Prepared) throw std::logic_error("immersed domain has no solved trial port state");
		return trial_measurements_.at(id);
	}
	void RollbackTrial() override
	{
		Local("distributed immersed rollback preflight",[&] { Require(Phase::TrialSolved,"rollback"); });
		runtime_.Rollback(); RestoreControls();
		inputs_.clear(); trial_controls_.clear(); trial_measurements_.clear(); phase_ = Phase::TrialReady; ++diagnostics_.rollback_count;
	}
	void AbortStep() override
	{
		Local("distributed immersed abort preflight",[&] { if (phase_ == Phase::Committed) throw std::logic_error("immersed abort requires an active step"); });
		if (phase_ == Phase::Prepared) runtime_.AbortPrepared();
		if (runtime_.Diagnostics().trial_active) runtime_.Rollback();
		RestoreControls(); inputs_.clear(); trial_controls_.clear(); trial_measurements_.clear(); prepared_measurements_.clear();
		phase_ = Phase::Committed; ++diagnostics_.abort_count;
	}
	void PrepareCommitStep() override
	{
		std::map<std::string,PortState> candidate;
		Local("distributed immersed prepare preflight",[&] {
			Require(Phase::TrialSolved,"prepare");
			for (const auto& entry : trial_measurements_) { ValidatePortState(entry.second); candidate.emplace(entry); }
		});
		runtime_.PrepareCommit(); prepared_measurements_.swap(candidate); phase_ = Phase::Prepared; ++diagnostics_.prepared_count;
	}
	void FinalizeCommitStep() noexcept override
	{
		if (phase_ != Phase::Prepared || !runtime_.Diagnostics().prepared) return;
		runtime_.FinalizeCommit(); committed_measurements_.swap(prepared_measurements_); committed_controls_.swap(trial_controls_);
		diagnostics_.committed_time_s = step_.EndTime(); diagnostics_.committed_step_index = step_.step_index; ++diagnostics_.committed_steps;
		inputs_.clear(); trial_measurements_.clear(); phase_ = Phase::Committed;
	}
private:
	enum class Phase { Committed,TrialReady,TrialSolved,Prepared };
	template<class Function> void Local(const char* stage,Function&& function) const { CollectiveLocalStage(communicator_,stage,std::forward<Function>(function)); }
	void Require(Phase phase,const char* action) const { if (phase_ != phase) throw std::logic_error(std::string("immersed domain cannot ")+action+" in this phase"); }
	static void Append(std::ostream& text,const std::string& value) { text << value.size() << ':' << value << ':'; }
	const CouplingPort& Port(const std::string& id) const
	{
		for (const auto& port : ports_) if (port.id == id) return port;
		throw std::out_of_range("immersed graph port is absent");
	}
	static PortQuantity Quantity(ImmersedFlowPortControlMode mode)
	{
		if (mode == ImmersedFlowPortControlMode::FlowRate) return PortQuantity::FlowRate;
		if (mode == ImmersedFlowPortControlMode::Pressure) return PortQuantity::MeanPressure;
		if (mode == ImmersedFlowPortControlMode::MeanNormalTraction) return PortQuantity::MeanNormalTraction;
		throw std::invalid_argument("unsupported immersed hydraulic control");
	}
	void RestoreControls() { for (const auto& control : committed_controls_) runtime_.SetPortControlValue(control.first,control.second); }
	PortState Measure(const CouplingPort& port) const
	{
		for (const auto& diagnostic : runtime_.Diagnostics().ports) if (diagnostic.id == port.id) {
			PortState result; result.time_s = step_.EndTime();
			if (port.provides.count(PortQuantity::Area)) result.area_m2 = diagnostic.measurement.area_m2;
			if (port.provides.count(PortQuantity::FlowRate)) result.outward_flow_m3_s = diagnostic.measurement.outward_flow_m3_s;
			if (port.provides.count(PortQuantity::MeanPressure)) result.mean_pressure_pa = diagnostic.measurement.mean_pressure_pa;
			if (port.provides.count(PortQuantity::MeanNormalTraction)) result.mean_normal_traction_pa = diagnostic.measurement.mean_normal_traction_pa;
			ValidatePortState(result); return result;
		}
		throw std::logic_error("immersed runtime did not measure a configured port");
	}
	ImmersedStaticDistributedRuntime& runtime_;
	MPI_Comm communicator_;
	std::string domain_id_;
	std::vector<CouplingPort> ports_;
	Phase phase_ = Phase::Committed;
	DomainStepContext step_;
	ThreeDImmersedFlowDomainDiagnostics diagnostics_;
	std::map<std::string,double> inputs_,committed_controls_,trial_controls_;
	std::map<std::string,PortState> trial_measurements_,prepared_measurements_,committed_measurements_;
};

} // namespace iga
#endif

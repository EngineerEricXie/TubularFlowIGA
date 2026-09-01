#ifndef IGA_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_TRANSIENT_FLOW_RUNTIME_HPP

#include "BoundaryFlow.hpp"
#include "CouplingPort.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "NavierStokesElement.hpp"
#include "OutletModel.hpp"
#include "OwnedRowAssembler.hpp"
#include "PressureTraction.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

#if defined(__GNUC__) || defined(__clang__)
#define IGA_FLOW_NOINLINE __attribute__((noinline))
#else
#define IGA_FLOW_NOINLINE
#endif

struct FlowPortMeasurements {
	std::map<int, double> flows;
	std::map<int, double> pressures;
	std::map<int, std::map<std::string, double>> species_fluxes;
	std::map<int, std::map<std::string, double>> species_concentrations;
};

struct FlowRuntimeSummary {
	PetscInt linear_iterations = 0;
	PetscReal state_l2 = 0.0;
	double velocity_l2 = 0.0;
	double pressure_l2 = 0.0;
};

struct FlowConvergenceMetrics {
	double continuity_l2 = 0.0;
	double continuity_sum = 0.0;
	double net_boundary_flow = 0.0;
	double absolute_boundary_flow = 0.0;
	double relative_mass_imbalance = 0.0;
};

enum class FlowStepPhase {
	Committed,
	TrialReady,
	TrialSolved,
	CommitPrepared
};

class TransientFlowRuntime {
public:
	TransientFlowRuntime(Database& database, MPI_Comm communicator, bool configured,
		bool transient, NavierStokesParameters parameters,
		ResolvedBoundaryConditions initial_boundaries,
		std::vector<int> labels,
		std::vector<std::array<double, 3>> boundary_velocity,
		std::set<std::int32_t> wall_trace_basis,
		std::vector<OutletModelState> outlet_models)
		: database_(database), communicator_(communicator), configured_(configured),
			transient_(transient), parameters_(parameters), boundaries_(std::move(initial_boundaries)),
			labels_(std::move(labels)), boundary_velocity_(std::move(boundary_velocity)),
			wall_trace_basis_(std::move(wall_trace_basis)),
			outlet_models_(std::move(outlet_models)), assembler_(database, communicator, 4)
	{
		if (labels_.size() != database_.header().nodes
			|| boundary_velocity_.size() != database_.header().nodes)
			throw std::runtime_error("flow boundary data do not match database nodes");
		constraint_velocity_mask_ = boundaries_.velocity_constrained;
		constraint_pressure_mask_ = boundaries_.pressure_constrained;
		MPI_Comm_rank(communicator_, &rank_);
		owned_elements_ = database_.LoadOwned(rank_);
		BuildBoundaryLabelIndex();
		jacobian_ = assembler_.CreateMatrix(true);
		state_ = assembler_.CreateVector();
		previous_ = assembler_.CreateVector();
		committed_state_ = assembler_.CreateVector();
		update_ = assembler_.CreateVector();
		rhs_ = assembler_.CreateVector();
		VecSet(state_, 0.0);
		VecSet(previous_, 0.0);
		VecSet(committed_state_, 0.0);
		VecSet(update_, 0.0);
		VecSet(rhs_, 0.0);
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
			const auto index = static_cast<std::size_t>(node);
			if (wall_trace_basis_.count(static_cast<std::int32_t>(node))
				|| boundaries_.velocity_constrained[index])
				for (int field = 0; field < 3; ++field)
					boundary_rows_.push_back(static_cast<PetscInt>(4*node+field));
			if (boundaries_.pressure_constrained[index])
				boundary_rows_.push_back(static_cast<PetscInt>(4*node+3));
		}
		BuildGhostScatter();
		KSPCreate(communicator_, &solver_);
		KSPSetType(solver_, KSPFGMRES);
		KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
		PC preconditioner = nullptr;
		KSPGetPC(solver_, &preconditioner);
		const bool scalable_default = database_.header().nodes
			>= kScalablePreconditionerNodeThreshold;
		if (scalable_default) {
			PCSetType(preconditioner, PCFIELDSPLIT);
			PCFieldSplitSetBlockSize(preconditioner, 4);
			const PetscInt velocity_fields[] = {0, 1, 2};
			const PetscInt pressure_fields[] = {3};
			PCFieldSplitSetFields(preconditioner, "0", 3, velocity_fields, velocity_fields);
			PCFieldSplitSetFields(preconditioner, "1", 1, pressure_fields, pressure_fields);
			PCFieldSplitSetType(preconditioner, PC_COMPOSITE_SCHUR);
			PCFieldSplitSetSchurFactType(preconditioner, PC_FIELDSPLIT_SCHUR_FACT_FULL);
			PCFieldSplitSetSchurPre(preconditioner, PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr);
			char requested_preconditioner[64]{};
			PetscBool has_preconditioner_override = PETSC_FALSE;
			PetscOptionsGetString(nullptr, nullptr, "-pc_type", requested_preconditioner,
				sizeof(requested_preconditioner), &has_preconditioner_override);
			if (!has_preconditioner_override
				|| std::string(requested_preconditioner) == PCFIELDSPLIT) {
				SetDefaultPetscOption("-fieldsplit_0_ksp_type", "preonly");
				SetDefaultPetscOption("-fieldsplit_0_pc_type", "gamg");
				SetDefaultPetscOption("-fieldsplit_1_ksp_type", "preonly");
				SetDefaultPetscOption("-fieldsplit_1_pc_type", "gamg");
			}
		} else PCSetType(preconditioner, PCBJACOBI);
		KSPSetFromOptions(solver_);
	}

	~TransientFlowRuntime()
	{
		KSPDestroy(&solver_);
		VecScatterDestroy(&scatter_);
		ISDestroy(&destination_rows_);
		VecDestroy(&ghost_previous_);
		VecDestroy(&ghost_state_);
		ISDestroy(&source_rows_);
		VecDestroy(&rhs_);
		VecDestroy(&update_);
		VecDestroy(&previous_);
		VecDestroy(&committed_state_);
		VecDestroy(&state_);
		MatDestroy(&jacobian_);
	}

	void InitializeState()
	{
		if (!transient_) return;
		std::vector<PetscScalar> values;
		values.reserve(boundary_rows_.size());
		for (const auto row : boundary_rows_) values.push_back(BoundaryValue(row));
		VecSetValues(state_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
			values.data(), INSERT_VALUES);
		OwnedRowAssembler::Assemble(state_);
		VecCopy(state_, previous_);
	}

	void InitializeState(const SimulationConfiguration& initial_configuration)
	{
		if (configured_) UpdateConfiguredBoundaries(initial_configuration);
		InitializeState();
	}

	void CopyStateToPrevious()
	{
		RequirePhase(FlowStepPhase::Committed, "CopyStateToPrevious");
		VecCopy(state_, previous_);
	}

	void BeginStep(int step, double physical_time, int maximum_newton,
		double nonlinear_relative_tolerance, double nonlinear_absolute_tolerance,
		double mass_relative_tolerance)
	{
		RequirePhase(FlowStepPhase::Committed, "BeginStep");
		ValidateSolveControls(step, physical_time, maximum_newton,
			nonlinear_relative_tolerance, nonlinear_absolute_tolerance,
			mass_relative_tolerance);
		VecCopy(state_, committed_state_);
		VecCopy(committed_state_, previous_);
		committed_boundaries_ = boundaries_;
		committed_pressure_tractions_ = pressure_tractions_;
		committed_outlet_models_ = outlet_models_;
		trial_step_ = step;
		trial_time_ = physical_time;
		trial_maximum_newton_ = maximum_newton;
		trial_nonlinear_relative_tolerance_ = nonlinear_relative_tolerance;
		trial_nonlinear_absolute_tolerance_ = nonlinear_absolute_tolerance;
		trial_mass_relative_tolerance_ = mass_relative_tolerance;
		trial_linear_iterations_ = 0;
		trial_solve_succeeded_ = false;
		has_trial_configuration_ = false;
		trial_pressure_overrides_.clear();
		phase_ = FlowStepPhase::TrialReady;
	}

	void SetTrialBoundaryConfiguration(const SimulationConfiguration& step_configuration)
	{
		RequirePhase(FlowStepPhase::TrialReady, "SetTrialBoundaryConfiguration");
		if (!configured_)
			throw std::runtime_error(
				"SetTrialBoundaryConfiguration requires a configured flow runtime");
		UpdateConfiguredBoundaries(step_configuration);
		trial_configuration_ = step_configuration;
		has_trial_configuration_ = true;
	}

	void SetPortInput(const CouplingPort& port, const PortBoundaryData& input)
	{
		RequirePhase(FlowStepPhase::TrialReady, "SetPortInput");
		if (configured_ && !has_trial_configuration_)
			throw std::runtime_error(
				"configured 3D flow trial requires SetTrialBoundaryConfiguration before SetPortInput");
		ValidateCouplingPort(port);
		ValidatePortBoundaryData(input);
		const auto time_scale = std::max({1.0, std::abs(trial_time_), std::abs(input.time_s)});
		if (std::abs(input.time_s-trial_time_) > 1e-12*time_scale)
			throw std::runtime_error("3D trial port input time does not match the active step");
		if (input.outward_flow_m3_s || input.total_pressure_pa
			|| !input.concentration.empty() || !input.outward_species_flux.empty())
			throw std::runtime_error(
				"3D flow trial ports currently support only mean pressure or mean normal traction");
		const int supplied = static_cast<int>(input.mean_pressure_pa.has_value())
			+static_cast<int>(input.mean_normal_traction_pa.has_value());
		if (supplied != 1)
			throw std::runtime_error(
				"3D flow trial port input requires exactly one of mean pressure or mean normal traction");
		const auto quantity = input.mean_pressure_pa ? PortQuantity::MeanPressure
			: PortQuantity::MeanNormalTraction;
		if (!port.requires.count(quantity))
			throw std::runtime_error(std::string("3D flow trial port does not require ")
				+PortQuantityName(quantity));
		const auto label = ParseBoundaryLabelLocator(port);
		if (!boundary_label_index_.count(label))
			throw std::runtime_error("3D flow trial port boundary label is absent from the database");
		for (std::size_t node = 0; node < labels_.size(); ++node)
			if (labels_[node] == label && boundaries_.pressure_constrained[node])
				throw std::runtime_error(
					"3D flow trial pressure cannot override a pressure Dirichlet boundary");
		for (const auto& model : committed_outlet_models_)
			if (model.label == label)
				throw std::runtime_error(
					"3D flow trial port pressure cannot override an active outlet model");
		const auto pressure = input.mean_pressure_pa
			? *input.mean_pressure_pa : -*input.mean_normal_traction_pa;
		trial_pressure_overrides_[label] = pressure;
		pressure_tractions_[label] = pressure;
	}

	IGA_FLOW_NOINLINE void SolveTrial()
	{
		RequirePhase(FlowStepPhase::TrialReady, "SolveTrial");
		if (configured_ && !has_trial_configuration_)
			throw std::runtime_error(
				"configured 3D flow trial requires SetTrialBoundaryConfiguration before SolveTrial");
		VecCopy(committed_state_, state_);
		VecCopy(committed_state_, previous_);
		outlet_models_ = committed_outlet_models_;
		trial_linear_iterations_ = 0;
		trial_solve_succeeded_ = false;
		try {
			SolveCurrentTrial();
			trial_solve_succeeded_ = true;
			phase_ = FlowStepPhase::TrialSolved;
		} catch (...) {
			phase_ = FlowStepPhase::TrialSolved;
			throw;
		}
	}

	void RollbackTrial()
	{
		RequirePhase(FlowStepPhase::TrialSolved, "RollbackTrial");
		RestoreCommittedSnapshot();
		phase_ = FlowStepPhase::TrialReady;
	}

	void AbortStep()
	{
		if (phase_ == FlowStepPhase::Committed) return;
		if (phase_ != FlowStepPhase::TrialReady && phase_ != FlowStepPhase::TrialSolved
			&& phase_ != FlowStepPhase::CommitPrepared)
			throw std::runtime_error("AbortStep is invalid in the current 3D flow lifecycle phase");
		RestoreCommittedSnapshot();
		trial_configuration_ = {};
		trial_step_ = -1;
		trial_time_ = 0.0;
		trial_maximum_newton_ = 0;
		trial_nonlinear_relative_tolerance_ = 0.0;
		trial_nonlinear_absolute_tolerance_ = 0.0;
		trial_mass_relative_tolerance_ = 0.0;
		phase_ = FlowStepPhase::Committed;
	}

	void PrepareCommitStep()
	{
		RequirePhase(FlowStepPhase::TrialSolved, "PrepareCommitStep");
		if (!trial_solve_succeeded_)
			throw std::runtime_error("PrepareCommitStep requires a successful 3D flow trial solve");
		phase_ = FlowStepPhase::CommitPrepared;
	}

	void FinalizeCommitStep() noexcept
	{
		if (phase_ != FlowStepPhase::CommitPrepared) std::terminate();
		total_linear_iterations_ += trial_linear_iterations_;
		trial_linear_iterations_ = 0;
		trial_solve_succeeded_ = false;
		has_trial_configuration_ = false;
		trial_pressure_overrides_.clear();
		phase_ = FlowStepPhase::Committed;
	}

	void CommitStep()
	{
		PrepareCommitStep();
		FinalizeCommitStep();
	}

	IGA_FLOW_NOINLINE void Advance(const SimulationConfiguration& step_configuration, int step,
		double physical_time, int maximum_newton, double nonlinear_relative_tolerance,
		double nonlinear_absolute_tolerance, double mass_relative_tolerance)
	{
		BeginStep(step, physical_time, maximum_newton, nonlinear_relative_tolerance,
			nonlinear_absolute_tolerance, mass_relative_tolerance);
		try {
			if (configured_) SetTrialBoundaryConfiguration(step_configuration);
			SolveTrial();
			CommitStep();
		} catch (...) {
			if (phase_ != FlowStepPhase::Committed) {
				RestoreCommittedSnapshot();
				phase_ = FlowStepPhase::Committed;
			}
			throw;
		}
	}

	IGA_FLOW_NOINLINE PortState GetPortState(const CouplingPort& port) const
	{
		RequirePhase(FlowStepPhase::TrialSolved, "GetPortState");
		if (!trial_solve_succeeded_)
			throw std::runtime_error("GetPortState requires a successful 3D flow trial solve");
		return MeasurePorts({port}, trial_time_, {}, {}).at(port.id);
	}

	FlowStepPhase Phase() const { return phase_; }
	PetscInt TrialLinearIterations() const { return trial_linear_iterations_; }
	std::optional<double> PressureTractionValue(int label) const
	{
		const auto found = pressure_tractions_.find(label);
		return found == pressure_tractions_.end()
			? std::optional<double>{} : std::optional<double>{found->second};
	}

private:
	IGA_FLOW_NOINLINE void SolveCurrentTrial()
	{
		std::vector<double> previous_capacitor_pressure(outlet_models_.size());
		for (std::size_t i = 0; i < outlet_models_.size(); ++i)
			previous_capacitor_pressure[i] = outlet_models_[i].capacitor_pressure;
		bool outlet_converged = false;
		const int maximum_outlet_iterations = outlet_models_.empty() ? 1 : 12;
		for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
			if (configured_) UpdateConfiguredBoundaries(trial_configuration_);
			const auto converged = SolveNonlinearStep(trial_step_, trial_time_, trial_maximum_newton_,
				trial_nonlinear_relative_tolerance_, trial_nonlinear_absolute_tolerance_,
				trial_mass_relative_tolerance_);
			if (!converged)
				throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON at physical step "
					+std::to_string(trial_step_+1));
			if (outlet_models_.empty()) {
				outlet_converged = true;
				break;
			}
			const auto flows = MeasureOutletFlows();
			const auto evaluated = EvaluateOutletCoupling(outlet_models_,
				previous_capacitor_pressure, flows, parameters_.dt);
			const auto tolerance = OutletCouplingTolerance(evaluated);
			if (rank_ == 0) std::cout << "step=" << trial_step_+1 << " time=" << trial_time_
				<< " outlet_iteration=" << coupling
				<< " pressure_change=" << evaluated.maximum_pressure_change
				<< " tolerance=" << tolerance << '\n';
			if (evaluated.maximum_pressure_change <= tolerance) {
				CommitOutletCoupling(outlet_models_, evaluated);
				for (const auto& model : outlet_models_)
					if (rank_ == 0) std::cout << "outlet label=" << model.label
						<< " flow=" << model.flow << " pressure=" << model.pressure
						<< " capacitor_pressure=" << model.capacitor_pressure << '\n';
				outlet_converged = true;
				break;
			}
			RelaxOutletCoupling(outlet_models_, evaluated);
		}
		if (!outlet_converged)
			throw std::runtime_error("outlet fixed-point iteration did not converge at physical step "
				+std::to_string(trial_step_+1));
	}

public:

	IGA_FLOW_NOINLINE double ReferenceBoundaryFlow(int label) const
	{
		double local = 0.0;
		for (const auto& element : owned_elements_) {
			std::vector<std::array<double, 4>> nodal(element.connectivity.size());
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				const auto node = static_cast<std::size_t>(element.connectivity[a]);
				for (int field = 0; field < 3; ++field)
					nodal[a][field] = boundary_velocity_[node][static_cast<std::size_t>(field)];
			}
			BodyFittedSurface4x4QuadratureProvider quadrature(element);
			local += IntegrateBoundaryFlow(element, nodal, quadrature.Rule(), label);
		}
		double global = 0.0;
		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		return global;
	}

	IGA_FLOW_NOINLINE std::map<std::string, PortState> MeasurePorts(
		const std::vector<CouplingPort>& ports, double physical_time,
		const std::vector<std::string>& species_fields,
		const std::vector<double>& species_state,
		const CompiledLinearSystem* transport_system = nullptr) const
	{
		RequireFinitePortValue("3D port measurement time_s", physical_time);
		if (ports.empty()) throw std::runtime_error("3D port measurement requires at least one port");
		ValidateCouplingPorts(ports);
		if (!species_fields.empty()
			&& species_state.size() != ghost_nodes_.size()*species_fields.size())
			throw std::runtime_error("3D port species state size does not match flow-required nodes");
		if (transport_system && transport_system->fields != species_fields)
			throw std::runtime_error(
				"3D port transport system fields do not match the supplied species fields");
		std::vector<std::vector<double>> advection(species_fields.size(),
			std::vector<double>(species_fields.size(), 0.0));
		std::vector<std::vector<double>> diffusion(species_fields.size(),
			std::vector<double>(species_fields.size(), 0.0));
		if (transport_system) {
			for (const auto& term : transport_system->terms) {
				if (term.kind == TermKind::Advection)
					advection.at(term.equation).at(term.trial) += term.coefficient;
				else if (term.kind == TermKind::Diffusion)
					diffusion.at(term.equation).at(term.trial) += term.coefficient;
			}
		} else
			for (std::size_t field = 0; field < species_fields.size(); ++field)
				advection[field][field] = 1.0;
		std::map<int, std::size_t> index;
		std::set<std::string> ids;
		for (std::size_t i = 0; i < ports.size(); ++i) {
			if (!ids.insert(ports[i].id).second)
				throw std::runtime_error("3D port measurement ids must be unique");
			const int label = ParseBoundaryLabelLocator(ports[i]);
			if (!index.emplace(label, i).second)
				throw std::runtime_error("3D port measurement boundary_label locators must be unique");
		}
		std::vector<double> local_flow(index.size(), 0.0), global_flow(index.size(), 0.0);
		std::vector<double> local_pressure(index.size(), 0.0), global_pressure(index.size(), 0.0);
		std::vector<double> local_area(index.size(), 0.0), global_area(index.size(), 0.0);
		std::vector<double> local_species(index.size()*species_fields.size(), 0.0);
		std::vector<double> global_species(index.size()*species_fields.size(), 0.0);
		std::vector<double> local_concentration(index.size()*species_fields.size(), 0.0);
		std::vector<double> global_concentration(index.size()*species_fields.size(), 0.0);
		ScatterState();
		const PetscScalar* values = nullptr;
		VecGetArrayRead(ghost_state_, &values);
		for (const auto& element : owned_elements_) {
			std::vector<std::array<double, 4>> nodal(element.connectivity.size());
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				const auto position = ghost_position_.at(element.connectivity[a]);
				for (int field = 0; field < 4; ++field)
					nodal[a][field] = PetscRealPart(values[4*position+field]);
			}
			BodyFittedSurface4x4QuadratureProvider quadrature(element);
			for (const auto& entry : index) {
				const auto port = entry.second;
				local_flow[port] += IntegrateBoundaryFlow(element, nodal, quadrature.Rule(), entry.first);
				const auto pressure = IntegrateBoundaryScalarAndArea(element, nodal,
					quadrature.Rule(), entry.first);
				local_pressure[port] += pressure[0];
				local_area[port] += pressure[1];
			}
			for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
				const auto found = index.find(element.boundary_labels[face]);
				if (found == index.end()) continue;
				const auto port = found->second;
				std::vector<std::vector<double>> species(element.connectivity.size(),
					std::vector<double>(species_fields.size()));
				for (std::size_t a = 0; a < element.connectivity.size(); ++a)
					for (std::size_t field = 0; field < species_fields.size(); ++field)
						species[a][field] = species_state[
							ghost_position_.at(element.connectivity[a])*species_fields.size()+field];
				for (std::size_t field = 0; field < species_fields.size(); ++field) {
					const auto measured = IntegrateBoundaryTransportFlux(element, face,
						nodal, species, field, advection[field], diffusion[field]);
					local_concentration[port*species_fields.size()+field]
						+= measured.concentration_integral;
					local_species[port*species_fields.size()+field]
						+= measured.total_outward_flux;
				}
			}
		}
		VecRestoreArrayRead(ghost_state_, &values);
		MPI_Allreduce(local_flow.data(), global_flow.data(), static_cast<int>(global_flow.size()),
			MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(local_pressure.data(), global_pressure.data(), static_cast<int>(global_pressure.size()),
			MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(local_area.data(), global_area.data(), static_cast<int>(global_area.size()),
			MPI_DOUBLE, MPI_SUM, communicator_);
		if (!species_fields.empty())
			MPI_Allreduce(local_species.data(), global_species.data(),
				static_cast<int>(global_species.size()), MPI_DOUBLE, MPI_SUM, communicator_);
		if (!species_fields.empty())
			MPI_Allreduce(local_concentration.data(), global_concentration.data(),
				static_cast<int>(global_concentration.size()), MPI_DOUBLE, MPI_SUM, communicator_);
		std::map<std::string, PortState> result;
		for (std::size_t i = 0; i < ports.size(); ++i) {
			if (!(global_area[i] > 0.0))
				throw std::runtime_error("3D port measurement has zero boundary area");
			PortState state;
			state.time_s = physical_time;
			state.area_m2 = global_area[i];
			state.outward_flow_m3_s = ports[i].orientation.ToOutward(global_flow[i]);
			state.mean_pressure_pa = global_pressure[i]/global_area[i];
			for (std::size_t field = 0; field < species_fields.size(); ++field) {
				state.concentration.emplace(species_fields[field],
					global_concentration[i*species_fields.size()+field]/global_area[i]);
				state.outward_species_flux.emplace(species_fields[field],
					ports[i].orientation.ToOutward(global_species[i*species_fields.size()+field]));
			}
			ValidatePortState(state);
			result.emplace(ports[i].id, std::move(state));
		}
		return result;
	}

	IGA_FLOW_NOINLINE FlowPortMeasurements MeasurePorts(const ThreeDVascularPortDefinition& ports,
		const std::vector<std::string>& species_fields,
		const std::vector<double>& species_state,
		const CompiledLinearSystem* transport_system = nullptr) const
	{
		std::vector<CouplingPort> generic_ports;
		std::vector<int> labels;
		if (ports.inlet_label >= 0) labels.push_back(ports.inlet_label);
		labels.insert(labels.end(), ports.outlet_labels.begin(), ports.outlet_labels.end());
		generic_ports.reserve(labels.size());
		for (const auto label : labels) {
			CouplingPort port;
			port.id = VcaPortId(label);
			port.subsystem_id = "three_d_vca";
			port.locator_kind = "boundary_label";
			port.locator = std::to_string(label);
			port.provides = {PortQuantity::Area, PortQuantity::FlowRate,
				PortQuantity::MeanPressure};
			if (!species_fields.empty()) port.provides.insert(PortQuantity::SpeciesFlux);
			generic_ports.push_back(std::move(port));
		}
		const auto measured = MeasurePorts(generic_ports, 0.0, species_fields,
			species_state, transport_system);
		FlowPortMeasurements result;
		for (const auto label : labels) {
			const auto found = measured.find(VcaPortId(label));
			if (found == measured.end() || !found->second.outward_flow_m3_s
				|| !found->second.mean_pressure_pa)
				throw std::runtime_error("VCA port measurement is incomplete");
			result.flows.emplace(label, *found->second.outward_flow_m3_s);
			result.pressures.emplace(label, *found->second.mean_pressure_pa);
			if (!species_fields.empty())
				result.species_fluxes.emplace(label, found->second.outward_species_flux);
			if (!species_fields.empty())
				result.species_concentrations.emplace(label, found->second.concentration);
		}
		return result;
	}

	IGA_FLOW_NOINLINE std::vector<std::array<double, 3>> GatherRequiredVelocity() const
	{
		ScatterState();
		const PetscScalar* values = nullptr;
		VecGetArrayRead(ghost_state_, &values);
		std::vector<std::array<double, 3>> velocity(ghost_nodes_.size());
		for (std::size_t node = 0; node < ghost_nodes_.size(); ++node)
			for (int component = 0; component < 3; ++component)
				velocity[static_cast<std::size_t>(node)][static_cast<std::size_t>(component)]
					= PetscRealPart(values[4*node+component]);
		VecRestoreArrayRead(ghost_state_, &values);
		return velocity;
	}

	IGA_FLOW_NOINLINE FlowRuntimeSummary Summary() const
	{
		FlowRuntimeSummary result;
		result.linear_iterations = total_linear_iterations_;
		VecNorm(state_, NORM_2, &result.state_l2);
		const PetscScalar* values = nullptr;
		VecGetArrayRead(state_, &values);
		double local_velocity_squared = 0.0;
		double local_pressure_squared = 0.0;
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
			const auto local = static_cast<std::size_t>(node-assembler_.node_begin())*4;
			for (int field = 0; field < 3; ++field) {
				const auto value = PetscRealPart(values[local+field]);
				local_velocity_squared += value*value;
			}
			const auto pressure = PetscRealPart(values[local+3]);
			local_pressure_squared += pressure*pressure;
		}
		VecRestoreArrayRead(state_, &values);
		double velocity_squared = 0.0, pressure_squared = 0.0;
		MPI_Allreduce(&local_velocity_squared, &velocity_squared, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(&local_pressure_squared, &pressure_squared, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		result.velocity_l2 = std::sqrt(velocity_squared);
		result.pressure_l2 = std::sqrt(pressure_squared);
		return result;
	}

	Vec State() const { return state_; }
	const OwnedRowAssembler& Assembler() const { return assembler_; }
	const std::vector<Element>& Elements() const { return assembler_.elements(); }
	const std::vector<Element>& OwnedElements() const { return owned_elements_; }
	const std::vector<std::int32_t>& RequiredNodes() const { return ghost_nodes_; }
	std::vector<OutletModelState>& OutletModels()
	{
		RequirePhase(FlowStepPhase::Committed, "mutable OutletModels");
		return outlet_models_;
	}
	const std::vector<OutletModelState>& OutletModels() const { return outlet_models_; }

private:
	static constexpr std::uint64_t kScalablePreconditionerNodeThreshold = 1000;

	void RestoreCommittedSnapshot()
	{
		VecCopy(committed_state_, state_);
		VecCopy(committed_state_, previous_);
		boundaries_ = committed_boundaries_;
		pressure_tractions_ = committed_pressure_tractions_;
		outlet_models_ = committed_outlet_models_;
		trial_linear_iterations_ = 0;
		trial_solve_succeeded_ = false;
		has_trial_configuration_ = false;
		trial_pressure_overrides_.clear();
	}

	void RequirePhase(FlowStepPhase required, const char* operation) const
	{
		if (phase_ != required)
			throw std::runtime_error(std::string(operation)
				+" is invalid in the current 3D flow lifecycle phase");
	}

	static void ValidateSolveControls(int step, double physical_time, int maximum_newton,
		double nonlinear_relative_tolerance, double nonlinear_absolute_tolerance,
		double mass_relative_tolerance)
	{
		if (step < 0) throw std::invalid_argument("physical step must be nonnegative");
		if (!std::isfinite(physical_time) || physical_time < 0.0)
			throw std::invalid_argument("physical time must be finite and nonnegative");
		if (maximum_newton <= 0)
			throw std::invalid_argument("maximum Newton iterations must be positive");
		if (!std::isfinite(nonlinear_relative_tolerance)
			|| nonlinear_relative_tolerance <= 0.0)
			throw std::invalid_argument(
				"nonlinear relative tolerance must be finite and positive");
		if (!std::isfinite(nonlinear_absolute_tolerance)
			|| nonlinear_absolute_tolerance <= 0.0)
			throw std::invalid_argument(
				"nonlinear absolute tolerance must be finite and positive");
		if (!std::isfinite(mass_relative_tolerance) || mass_relative_tolerance <= 0.0)
			throw std::invalid_argument("mass relative tolerance must be finite and positive");
	}

	static int ParseBoundaryLabelLocator(const CouplingPort& port)
	{
		if (port.locator_kind != "boundary_label")
			throw std::runtime_error("3D port measurement supports only boundary_label locators");
		int label = -1;
		const auto parsed = std::from_chars(port.locator.data(),
			port.locator.data()+port.locator.size(), label);
		if (parsed.ec != std::errc{} || parsed.ptr != port.locator.data()+port.locator.size()
			|| label < 0)
			throw std::runtime_error("3D boundary_label locator must be a nonnegative integer");
		return label;
	}

	static std::string VcaPortId(int label)
	{
		return "vca_boundary_label_"+std::to_string(label);
	}

	static void SetDefaultPetscOption(const char* name, const char* value)
	{
		PetscBool present = PETSC_FALSE;
		PetscOptionsHasName(nullptr, nullptr, name, &present);
		if (!present) PetscOptionsSetValue(nullptr, name, value);
	}

	void BuildBoundaryLabelIndex()
	{
		std::set<int> local;
		for (const auto& element : owned_elements_)
			for (const auto label : element.boundary_labels)
				if (label >= 0) local.insert(label);
		const int local_count = static_cast<int>(local.size());
		int size = 0;
		MPI_Comm_size(communicator_, &size);
		std::vector<int> counts(static_cast<std::size_t>(size));
		MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator_);
		std::vector<int> offsets(static_cast<std::size_t>(size), 0);
		for (int i = 1; i < size; ++i)
			offsets[static_cast<std::size_t>(i)] = offsets[static_cast<std::size_t>(i-1)]
				+counts[static_cast<std::size_t>(i-1)];
		const auto total = offsets.back()+counts.back();
		std::vector<int> local_labels(local.begin(), local.end());
		std::vector<int> gathered(static_cast<std::size_t>(total));
		MPI_Allgatherv(local_labels.data(), local_count, MPI_INT, gathered.data(), counts.data(),
			offsets.data(), MPI_INT, communicator_);
		std::sort(gathered.begin(), gathered.end());
		gathered.erase(std::unique(gathered.begin(), gathered.end()), gathered.end());
		boundary_labels_ = std::move(gathered);
		for (std::size_t i = 0; i < boundary_labels_.size(); ++i)
			boundary_label_index_.emplace(boundary_labels_[i], i);
	}

	FlowConvergenceMetrics MeasureConvergence(const PetscScalar* ghost_values) const
	{
		FlowConvergenceMetrics result;
		const PetscScalar* residual_values = nullptr;
		VecGetArrayRead(rhs_, &residual_values);
		double local_continuity_squared = 0.0;
		double local_continuity_sum = 0.0;
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
			const auto local = static_cast<std::size_t>(node-assembler_.node_begin())*4+3;
			const auto value = PetscRealPart(residual_values[local]);
			local_continuity_squared += value*value;
			local_continuity_sum += value;
		}
		VecRestoreArrayRead(rhs_, &residual_values);
		double global_continuity_squared = 0.0;
		MPI_Allreduce(&local_continuity_squared, &global_continuity_squared, 1,
			MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(&local_continuity_sum, &result.continuity_sum, 1,
			MPI_DOUBLE, MPI_SUM, communicator_);
		result.continuity_l2 = std::sqrt(global_continuity_squared);

		std::vector<double> local_flow(boundary_labels_.size(), 0.0);
		std::vector<double> global_flow(boundary_labels_.size(), 0.0);
		for (const auto& element : owned_elements_) {
			std::vector<std::array<double, 4>> nodal(element.connectivity.size());
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				const auto position = ghost_position_.at(element.connectivity[a]);
				for (int field = 0; field < 4; ++field)
					nodal[a][field] = PetscRealPart(ghost_values[4*position+field]);
			}
			BodyFittedSurface4x4QuadratureProvider quadrature(element);
			for (const auto& entry : boundary_label_index_)
				local_flow[entry.second] += IntegrateBoundaryFlow(element, nodal,
					quadrature.Rule(), entry.first);
		}
		if (!local_flow.empty())
			MPI_Allreduce(local_flow.data(), global_flow.data(), static_cast<int>(local_flow.size()),
				MPI_DOUBLE, MPI_SUM, communicator_);
		for (const auto flow : global_flow) {
			result.net_boundary_flow += flow;
			result.absolute_boundary_flow += std::abs(flow);
		}
		if (result.absolute_boundary_flow > 0.0)
			result.relative_mass_imbalance = 2.0*std::abs(result.net_boundary_flow)
				/result.absolute_boundary_flow;
		return result;
	}

	static std::string PreciseNumber(double value)
	{
		std::ostringstream stream;
		stream << std::scientific << std::setprecision(16) << value;
		return stream.str();
	}

	void BuildGhostScatter()
	{
		std::vector<std::int32_t> nodes;
		for (const auto& element : assembler_.elements())
			nodes.insert(nodes.end(), element.connectivity.begin(), element.connectivity.end());
		for (const auto& element : owned_elements_)
			nodes.insert(nodes.end(), element.connectivity.begin(), element.connectivity.end());
		std::sort(nodes.begin(), nodes.end());
		nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
		ghost_nodes_ = std::move(nodes);
		std::vector<PetscInt> rows;
		rows.reserve(4*ghost_nodes_.size());
		for (std::size_t i = 0; i < ghost_nodes_.size(); ++i) {
			ghost_position_.emplace(ghost_nodes_[i], i);
			for (int field = 0; field < 4; ++field) rows.push_back(4*ghost_nodes_[i]+field);
		}
		ISCreateGeneral(communicator_, static_cast<PetscInt>(rows.size()), rows.data(),
			PETSC_COPY_VALUES, &source_rows_);
		VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), &ghost_state_);
		VecDuplicate(ghost_state_, &ghost_previous_);
		ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), 0, 1,
			&destination_rows_);
		VecScatterCreate(state_, source_rows_, ghost_state_, destination_rows_, &scatter_);
	}

	void ScatterState() const
	{
		VecScatterBegin(scatter_, state_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter_, state_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD);
	}

	double BoundaryValue(PetscInt row) const
	{
		const auto node = static_cast<std::uint64_t>(row/4);
		const auto field = row%4;
		const auto index = static_cast<std::size_t>(node);
		return field < 3 ? (wall_trace_basis_.count(static_cast<std::int32_t>(node)) ? 0.0
			: boundaries_.velocity[index][static_cast<std::size_t>(field)])
			: boundaries_.pressure[index];
	}

	void UpdateConfiguredBoundaries(const SimulationConfiguration& step_configuration)
	{
		const auto configuration = MaterializeOutletPressures(step_configuration, outlet_models_);
		const auto& system = FirstNavierStokesSystem(configuration);
		auto pressure_tractions = ExtractPressureTractions(configuration, system);
		auto boundaries = ResolveFlowBoundaries(configuration, system, labels_, boundary_velocity_);
		ValidateConstraintTopology(boundaries);
		for (const auto& override : trial_pressure_overrides_)
			pressure_tractions[override.first] = override.second;
		pressure_tractions_ = std::move(pressure_tractions);
		boundaries_ = std::move(boundaries);
	}

	void ValidateConstraintTopology(const ResolvedBoundaryConditions& boundaries) const
	{
		if (boundaries.velocity_constrained != constraint_velocity_mask_
			|| boundaries.pressure_constrained != constraint_pressure_mask_)
			throw std::runtime_error(
				"configured 3D flow boundary constraint topology changed; PETSc boundary rows are fixed");
	}

	IGA_FLOW_NOINLINE bool SolveNonlinearStep(int step, double physical_time, int maximum_newton,
		double nonlinear_relative_tolerance, double nonlinear_absolute_tolerance,
		double mass_relative_tolerance)
	{
		PetscReal initial_residual = -1.0;
		for (int nonlinear = 0; nonlinear < maximum_newton; ++nonlinear) {
			const auto iteration_start = std::chrono::steady_clock::now();
			MatZeroEntries(jacobian_);
			VecSet(rhs_, 0.0);
			ScatterState();
			if (transient_) {
				VecScatterBegin(scatter_, previous_, ghost_previous_, INSERT_VALUES, SCATTER_FORWARD);
				VecScatterEnd(scatter_, previous_, ghost_previous_, INSERT_VALUES, SCATTER_FORWARD);
			}
			const PetscScalar* values = nullptr;
			VecGetArrayRead(ghost_state_, &values);
			const PetscScalar* previous_values = nullptr;
			if (transient_) VecGetArrayRead(ghost_previous_, &previous_values);
			for (const auto& element : assembler_.elements()) {
				std::vector<std::array<double, 4>> nodal(element.connectivity.size());
				std::vector<std::array<double, 4>> previous_nodal;
				if (transient_) previous_nodal.resize(element.connectivity.size());
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto position = ghost_position_.at(element.connectivity[a]);
					for (int field = 0; field < 4; ++field) {
						nodal[a][field] = PetscRealPart(values[4*position+field]);
						if (transient_)
							previous_nodal[a][field] = PetscRealPart(previous_values[4*position+field]);
					}
				}
				FullCell4x4x4VolumeQuadratureProvider volume_quadrature(element);
				auto local = BuildNavierStokesElement(element, nodal, previous_nodal,
					parameters_, volume_quadrature.Rule());
				BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
				for (const auto& traction : pressure_tractions_) {
					const auto surface = IntegrateBoundaryPressureTraction(element,
						surface_quadrature.Rule(), traction.first, traction.second);
					for (std::size_t row = 0; row < surface.size(); ++row)
						local.negative_residual[row] += surface[row];
				}
				assembler_.AddElementMatrix(jacobian_, element, local.jacobian);
				assembler_.AddElementVector(rhs_, element, local.negative_residual);
			}
			if (transient_) VecRestoreArrayRead(ghost_previous_, &previous_values);
			OwnedRowAssembler::Assemble(jacobian_);
			OwnedRowAssembler::Assemble(rhs_);
			const auto convergence = MeasureConvergence(values);
			VecRestoreArrayRead(ghost_state_, &values);
			const PetscScalar* owned = nullptr;
			VecGetArrayRead(state_, &owned);
			std::vector<PetscScalar> boundary_update;
			boundary_update.reserve(boundary_rows_.size());
			for (const auto row : boundary_rows_) {
				const auto local = static_cast<std::size_t>(row-4*assembler_.node_begin());
				boundary_update.push_back(BoundaryValue(row)-PetscRealPart(owned[local]));
			}
			VecRestoreArrayRead(state_, &owned);
			MatZeroRows(jacobian_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
				1.0, nullptr, nullptr);
			VecSetValues(rhs_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
				boundary_update.data(), INSERT_VALUES);
			OwnedRowAssembler::Assemble(rhs_);
			PetscReal residual = 0.0;
			VecNorm(rhs_, NORM_2, &residual);
			if (initial_residual < 0.0) initial_residual = residual;
			const auto residual_scale = std::sqrt(
				static_cast<PetscReal>(assembler_.global_rows()));
			const auto residual_rms = residual/residual_scale;
			const auto tolerance = std::max<PetscReal>(
				nonlinear_absolute_tolerance*residual_scale,
				nonlinear_relative_tolerance*initial_residual);
			if (residual <= tolerance
				&& convergence.relative_mass_imbalance <= mass_relative_tolerance) {
				if (rank_ == 0) std::cout << "step=" << step+1 << " time=" << physical_time
					<< " converged newton=" << nonlinear << " residual_l2=" << residual
					<< " residual_rms=" << residual_rms << " tolerance=" << tolerance
					<< " absolute_rms_tolerance=" << nonlinear_absolute_tolerance
					<< " continuity_l2=" << convergence.continuity_l2
					<< " continuity_sum=" << convergence.continuity_sum
					<< " net_boundary_flow=" << convergence.net_boundary_flow
					<< " relative_mass_imbalance=" << convergence.relative_mass_imbalance
					<< " mass_tolerance=" << mass_relative_tolerance << " assembly_s="
					<< std::chrono::duration<double>(std::chrono::steady_clock::now()-iteration_start).count() << '\n';
				return true;
			}
			const auto linear_start = std::chrono::steady_clock::now();
			KSPSetOperators(solver_, jacobian_, jacobian_);
			KSPSolve(solver_, rhs_, update_);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver_, &reason);
			if (reason <= 0) {
				PetscInt failed_iterations = 0;
				PetscReal failed_residual = 0.0;
				KSPGetIterationNumber(solver_, &failed_iterations);
				KSPGetResidualNorm(solver_, &failed_residual);
				throw std::runtime_error("Navier-Stokes linear solve failed at nonlinear iteration "
					+std::to_string(nonlinear) + " (KSP reason="
					+std::to_string(static_cast<int>(reason)) + ", iterations="
					+std::to_string(static_cast<long long>(failed_iterations)) + ", residual="
					+PreciseNumber(static_cast<double>(failed_residual)) + ", nonlinear_residual="
					+PreciseNumber(static_cast<double>(residual)) + ", continuity_l2="
					+PreciseNumber(convergence.continuity_l2) + ", net_boundary_flow="
					+PreciseNumber(convergence.net_boundary_flow) + ", relative_mass_imbalance="
					+PreciseNumber(convergence.relative_mass_imbalance) + ")");
			}
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver_, &iterations);
			trial_linear_iterations_ += iterations;
			PetscReal linear_residual = 0.0, update_norm = 0.0;
			KSPGetResidualNorm(solver_, &linear_residual);
			VecNorm(update_, NORM_2, &update_norm);
			VecAXPY(state_, 1.0, update_);
			if (rank_ == 0) std::cout << "step=" << step+1 << " time=" << physical_time
				<< " newton=" << nonlinear << " residual_l2=" << residual
				<< " residual_rms=" << residual_rms
				<< " continuity_l2=" << convergence.continuity_l2
				<< " continuity_sum=" << convergence.continuity_sum
				<< " net_boundary_flow=" << convergence.net_boundary_flow
				<< " relative_mass_imbalance=" << convergence.relative_mass_imbalance
				<< " update_l2=" << update_norm << " linear_iterations=" << iterations
				<< " linear_residual=" << linear_residual << " assembly_s="
				<< std::chrono::duration<double>(linear_start-iteration_start).count()
				<< " linear_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-linear_start).count() << '\n';
		}
		return false;
	}

	IGA_FLOW_NOINLINE std::vector<double> MeasureOutletFlows() const
	{
		std::vector<double> local(outlet_models_.size(), 0.0), global(outlet_models_.size(), 0.0);
		std::unordered_map<int, std::size_t> indices;
		for (std::size_t i = 0; i < outlet_models_.size(); ++i)
			indices.emplace(outlet_models_[i].label, i);
		ScatterState();
		const PetscScalar* values = nullptr;
		VecGetArrayRead(ghost_state_, &values);
		for (const auto& element : owned_elements_) {
			std::vector<std::array<double, 4>> nodal(element.connectivity.size());
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				const auto position = ghost_position_.at(element.connectivity[a]);
				for (int field = 0; field < 4; ++field)
					nodal[a][field] = PetscRealPart(values[4*position+field]);
			}
			BodyFittedSurface4x4QuadratureProvider quadrature(element);
			for (const auto& entry : indices)
				local[entry.second] += IntegrateBoundaryFlow(element, nodal,
					quadrature.Rule(), entry.first);
		}
		VecRestoreArrayRead(ghost_state_, &values);
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()),
			MPI_DOUBLE, MPI_SUM, communicator_);
		return global;
	}

	Database& database_;
	MPI_Comm communicator_;
	bool configured_ = false;
	bool transient_ = false;
	NavierStokesParameters parameters_;
	ResolvedBoundaryConditions boundaries_;
	std::vector<int> constraint_velocity_mask_;
	std::vector<int> constraint_pressure_mask_;
	std::vector<int> labels_;
	std::vector<std::array<double, 3>> boundary_velocity_;
	std::set<std::int32_t> wall_trace_basis_;
	std::vector<OutletModelState> outlet_models_;
	std::vector<Element> owned_elements_;
	std::vector<int> boundary_labels_;
	std::unordered_map<int, std::size_t> boundary_label_index_;
	OwnedRowAssembler assembler_;
	std::vector<PetscInt> boundary_rows_;
	std::vector<std::int32_t> ghost_nodes_;
	std::unordered_map<std::int32_t, std::size_t> ghost_position_;
	std::map<int, double> pressure_tractions_;
	Mat jacobian_ = nullptr;
	Vec state_ = nullptr, previous_ = nullptr, committed_state_ = nullptr;
	Vec update_ = nullptr, rhs_ = nullptr;
	IS source_rows_ = nullptr, destination_rows_ = nullptr;
	Vec ghost_state_ = nullptr, ghost_previous_ = nullptr;
	VecScatter scatter_ = nullptr;
	KSP solver_ = nullptr;
	PetscInt total_linear_iterations_ = 0;
	PetscInt trial_linear_iterations_ = 0;
	FlowStepPhase phase_ = FlowStepPhase::Committed;
	ResolvedBoundaryConditions committed_boundaries_;
	std::map<int, double> committed_pressure_tractions_;
	std::vector<OutletModelState> committed_outlet_models_;
	SimulationConfiguration trial_configuration_;
	std::map<int, double> trial_pressure_overrides_;
	int trial_step_ = -1;
	double trial_time_ = 0.0;
	int trial_maximum_newton_ = 0;
	double trial_nonlinear_relative_tolerance_ = 0.0;
	double trial_nonlinear_absolute_tolerance_ = 0.0;
	double trial_mass_relative_tolerance_ = 0.0;
	bool has_trial_configuration_ = false;
	bool trial_solve_succeeded_ = false;
	int rank_ = 0;
};

} // namespace iga

#undef IGA_FLOW_NOINLINE

#endif

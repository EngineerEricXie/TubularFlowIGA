#ifndef IGA_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_TRANSIENT_FLOW_RUNTIME_HPP

#include "BoundaryFlow.hpp"
#include "CouplingPort.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "NavierStokesElement.hpp"
#include "OutletModel.hpp"
#include "OwnedRowAssembler.hpp"
#include "ElementAssemblyExecution.hpp"
#include "CollectivePetscOptions.hpp"
#include <optional>
#include "PetscReadArray.hpp"
#include "PressureTraction.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <charconv>
#include "PetscPhaseProfile.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
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
	// Borrows communicator for the runtime lifetime. The owner must destroy
	// the runtime before freeing the communicator and before PetscFinalize.
	TransientFlowRuntime(Database& database, MPI_Comm communicator, bool configured,
		bool transient, NavierStokesParameters parameters,
		const ResolvedBoundaryConditions& initial_boundaries,
		const std::vector<int>& labels,
		const std::vector<std::array<double, 3>>& boundary_velocity,
		const std::set<std::int32_t>& wall_trace_basis,
		const std::vector<OutletModelState>& outlet_models,
		const std::set<std::string>& application_options = {})
		: database_(database), communicator_(communicator), configured_(configured),
			transient_(transient), parameters_(parameters), assembler_(database, communicator, 4),
			trial_configuration_(PrepareRuntimeConstructionInput<SimulationConfiguration>(communicator,
				"flow trial configuration preparation"))
	{
		MPI_Comm_rank(communicator_, &rank_);
		std::array<const char*, 4> added_defaults{};
		std::size_t added_default_count = 0;
		try {
			CollectiveLocalStage(communicator_, "flow assembly resources", [&] { assembly_execution_.emplace(); });
			RuntimeConstructionStage(communicator_, "flow runtime input", [&] {
				boundaries_ = initial_boundaries;
				labels_ = labels;
				ProbeRuntimeConstruction("flow boundary input copied");
				boundary_velocity_ = boundary_velocity;
				wall_trace_basis_ = wall_trace_basis;
				outlet_models_ = outlet_models;
				if (labels_.size() != database_.header().nodes
					|| boundary_velocity_.size() != database_.header().nodes
					|| boundaries_.velocity_constrained.size() != database_.header().nodes
					|| boundaries_.pressure_constrained.size() != database_.header().nodes
					|| boundaries_.velocity.size() != database_.header().nodes
					|| boundaries_.pressure.size() != database_.header().nodes)
					throw std::runtime_error("flow boundary data do not match database nodes");
				constraint_velocity_mask_ = boundaries_.velocity_constrained;
				constraint_pressure_mask_ = boundaries_.pressure_constrained;
				owned_elements_ = database_.LoadOwned(rank_);
			});
			BuildBoundaryLabelIndex();
			jacobian_ = assembler_.CreateMatrix(true);
			state_ = assembler_.CreateVector();
			previous_ = assembler_.CreateVector();
			committed_state_ = assembler_.CreateVector();
			update_ = assembler_.CreateVector();
			rhs_ = assembler_.CreateVector();
			RequireCollectivePetscSuccess(communicator_, "flow initial state_ clear", VecSet(state_, 0.0));
			RequireCollectivePetscSuccess(communicator_, "flow initial previous_ clear", VecSet(previous_, 0.0));
			RequireCollectivePetscSuccess(communicator_, "flow initial committed_state_ clear", VecSet(committed_state_, 0.0));
			RequireCollectivePetscSuccess(communicator_, "flow initial update_ clear", VecSet(update_, 0.0));
			RequireCollectivePetscSuccess(communicator_, "flow initial rhs_ clear", VecSet(rhs_, 0.0));
			RuntimeConstructionStage(communicator_, "flow boundary row preparation", [&] {
				for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
					const auto index = static_cast<std::size_t>(node);
					if (wall_trace_basis_.count(static_cast<std::int32_t>(node))
						|| boundaries_.velocity_constrained[index])
						for (int field = 0; field < 3; ++field)
							boundary_rows_.push_back(static_cast<PetscInt>(4*node+field));
					if (boundaries_.pressure_constrained[index])
						boundary_rows_.push_back(static_cast<PetscInt>(4*node+3));
				}
			});
			BuildGhostScatter();
			// The caller has separately validated these application controls;
			// all remaining (including prefixed/unused) solver options agree.
			RequireCollectivePetscOptions(communicator_, nullptr, application_options);
			RequireCollectivePetscSuccess(communicator_, "flow solver creation", KSPCreate(communicator_, &solver_));
			ObserveConstructedObject(communicator_, "flow solver created", reinterpret_cast<PetscObject>(solver_));
			RequireCollectivePetscSuccess(communicator_, "flow solver type", KSPSetType(solver_, KSPFGMRES));
			RequireCollectivePetscSuccess(communicator_, "flow solver tolerances", KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 5000));
			PC preconditioner = nullptr;
			RequireCollectivePetscSuccess(communicator_, "flow preconditioner lookup", KSPGetPC(solver_, &preconditioner));
			const bool scalable_default = database_.header().nodes
				>= kScalablePreconditionerNodeThreshold;
			if (scalable_default) {
				RequireCollectivePetscSuccess(communicator_, "flow fieldsplit type", PCSetType(preconditioner, PCFIELDSPLIT));
				RequireCollectivePetscSuccess(communicator_, "flow fieldsplit block size", PCFieldSplitSetBlockSize(preconditioner, 4));
				const PetscInt velocity_fields[] = {0, 1, 2};
				const PetscInt pressure_fields[] = {3};
				RequireCollectivePetscSuccess(communicator_, "flow velocity fields", PCFieldSplitSetFields(preconditioner, "0", 3, velocity_fields, velocity_fields));
				RequireCollectivePetscSuccess(communicator_, "flow pressure fields", PCFieldSplitSetFields(preconditioner, "1", 1, pressure_fields, pressure_fields));
				RequireCollectivePetscSuccess(communicator_, "flow Schur composition", PCFieldSplitSetType(preconditioner, PC_COMPOSITE_SCHUR));
				RequireCollectivePetscSuccess(communicator_, "flow Schur factorization", PCFieldSplitSetSchurFactType(preconditioner, PC_FIELDSPLIT_SCHUR_FACT_FULL));
				RequireCollectivePetscSuccess(communicator_, "flow Schur preconditioner", PCFieldSplitSetSchurPre(preconditioner, PC_FIELDSPLIT_SCHUR_PRE_A11, nullptr));
				char requested_preconditioner[64]{};
				PetscBool has_preconditioner_override = PETSC_FALSE;
				RequireCollectivePetscSuccess(communicator_, "flow preconditioner override",
					PetscOptionsGetString(nullptr, nullptr, "-pc_type", requested_preconditioner,
						sizeof(requested_preconditioner), &has_preconditioner_override));
				if (!has_preconditioner_override
					|| std::strcmp(requested_preconditioner, PCFIELDSPLIT) == 0) {
					RuntimeConstructionStage(communicator_, "flow default solver options", [&] {
						SetDefaultPetscOption("-fieldsplit_0_ksp_type", "preonly", added_defaults, added_default_count);
						SetDefaultPetscOption("-fieldsplit_0_pc_type", "gamg", added_defaults, added_default_count);
						SetDefaultPetscOption("-fieldsplit_1_ksp_type", "preonly", added_defaults, added_default_count);
						SetDefaultPetscOption("-fieldsplit_1_pc_type", "gamg", added_defaults, added_default_count);
					});
				}
			} else RequireCollectivePetscSuccess(communicator_, "flow block Jacobi type", PCSetType(preconditioner, PCBJACOBI));
			RequireCollectivePetscSuccess(communicator_, "flow solver options", KSPSetFromOptions(solver_));
			RuntimeConstructionStage(communicator_, "flow runtime ready", [] {});
		} catch (...) {
			// Only remove defaults absent on entry. Each rank may have reached
			// a different insertion before a local failure was coordinated.
			PetscErrorCode rollback_error = 0;
			while (added_default_count) {
				const auto code = PetscOptionsClearValue(nullptr, added_defaults[--added_default_count]);
				if (!rollback_error) rollback_error = code;
			}
			DestroyPetsc();
			RequireCollectivePetscSuccess(communicator_, "flow construction options rollback", rollback_error);
			throw;
		}
	}

	~TransientFlowRuntime() { DestroyPetsc(); }
	TransientFlowRuntime(const TransientFlowRuntime&) = delete;
	TransientFlowRuntime& operator=(const TransientFlowRuntime&) = delete;

	void InitializeState()
	{
		InitializeStateImpl(nullptr);
	}

	void InitializeState(const SimulationConfiguration& initial_configuration)
	{
		InitializeStateImpl(&initial_configuration);
	}

private:
	void DestroyPetsc() noexcept
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

	void InitializeStateImpl(const SimulationConfiguration* initial_configuration)
	{
		ResolvedBoundaryConditions candidate;
		std::map<int, double> tractions;
		std::vector<PetscScalar> values;
		std::string signature;
		CollectiveLocalStage(communicator_, "flow initialization preparation", [&] {
			RequirePhase(FlowStepPhase::Committed, "InitializeState");
			if (configured_ && initial_configuration)
				ResolveConfiguredBoundaries(*initial_configuration, candidate, tractions);
			else {
				candidate = boundaries_;
				tractions = pressure_tractions_;
			}
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << configured_ << ' ' << transient_ << ' ' << (initial_configuration != nullptr)
				<< ' ' << labels_.size() << std::hexfloat;
			for (std::size_t node = 0; node < labels_.size(); ++node) {
				const bool wall = wall_trace_basis_.count(static_cast<std::int32_t>(node));
				const bool velocity = candidate.velocity_constrained.at(node) || wall;
				const bool pressure = candidate.pressure_constrained.at(node);
				text << ' ' << velocity << ' ' << pressure;
				for (int field = 0; field < 4; ++field) if (field < 3 ? velocity : pressure) {
					const auto value = BoundaryValue(static_cast<PetscInt>(4*node+field), candidate);
					if (!std::isfinite(value)) throw std::runtime_error("nonfinite initial flow boundary value");
					text << ' ' << value;
				}
			}
			text << ' ' << tractions.size();
			for (const auto& traction : tractions) {
				if (!std::isfinite(traction.second)) throw std::runtime_error("nonfinite initial pressure traction");
				text << ' ' << traction.first << ' ' << traction.second;
			}
			signature = text.str();
			if (transient_) {
				values.reserve(boundary_rows_.size());
				for (const auto row : boundary_rows_) values.push_back(BoundaryValue(row, candidate));
			}
		});
		RequireCollectiveSameText(communicator_, "flow initialization agreement", signature);
		if (transient_) {
			// Existing scratch vectors have identical ownership/layout. Construct
			// both candidate fields before publishing into the stable Vec handles.
			RequireCollectivePetscSuccess(communicator_, "flow initialization copy candidate", VecCopy(state_, update_));
			CollectiveLocalStage(communicator_, "flow initialization insertion", [&] {
				if (VecSetValues(update_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
					values.data(), INSERT_VALUES)) throw std::runtime_error("initial flow VecSetValues failed");
			});
			OwnedRowAssembler::Assemble(update_, communicator_);
			RequireCollectivePetscSuccess(communicator_, "flow initialization history candidate", VecCopy(update_, rhs_));
			RequireCollectivePetscSuccess(communicator_, "flow initialization publish state", VecSwap(state_, update_));
			RequireCollectivePetscSuccess(communicator_, "flow initialization publish history", VecSwap(previous_, rhs_));
		}
		std::swap(boundaries_, candidate);
		pressure_tractions_.swap(tractions);
	}

public:

	void CopyStateToPrevious()
	{
		CollectiveLocalStage(communicator_, "flow copy history preparation", [&] {
			RequirePhase(FlowStepPhase::Committed, "CopyStateToPrevious");
		});
		RequireCollectivePetscSuccess(communicator_, "flow copy history", VecCopy(state_, previous_));
	}

	void BeginStep(int step, double physical_time, int maximum_newton,
		double nonlinear_relative_tolerance, double nonlinear_absolute_tolerance,
		double mass_relative_tolerance)
	{
		ResolvedBoundaryConditions boundaries;
		std::map<int, double> tractions;
		std::vector<OutletModelState> outlets;
		std::string controls;
		CollectiveLocalStage(communicator_, "flow begin step preparation", [&] {
			RequirePhase(FlowStepPhase::Committed, "BeginStep");
			ValidateSolveControls(step, physical_time, maximum_newton,
				nonlinear_relative_tolerance, nonlinear_absolute_tolerance, mass_relative_tolerance);
			boundaries = boundaries_; tractions = pressure_tractions_; outlets = outlet_models_;
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << step << ' ' << maximum_newton << std::hexfloat << ' ' << physical_time
				<< ' ' << nonlinear_relative_tolerance << ' ' << nonlinear_absolute_tolerance
				<< ' ' << mass_relative_tolerance;
			controls = text.str();
		});
		RequireCollectiveSameText(communicator_, "flow step controls agreement", controls);
		RequireCollectivePetscSuccess(communicator_, "flow save state", VecCopy(state_, committed_state_));
		RequireCollectivePetscSuccess(communicator_, "flow save history", VecCopy(committed_state_, previous_));
		std::swap(committed_boundaries_, boundaries);
		committed_pressure_tractions_.swap(tractions);
		committed_outlet_models_.swap(outlets);
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
		assembly_execution_->RequireCaller();
		std::vector<OutletModelState> outlets;
		CollectiveLocalStage(communicator_, "flow solve preparation", [&] {
			RequirePhase(FlowStepPhase::TrialReady, "SolveTrial");
			if (configured_ && !has_trial_configuration_)
				throw std::runtime_error(
					"configured 3D flow trial requires SetTrialBoundaryConfiguration before SolveTrial");
			outlets = committed_outlet_models_;
		});
		RequireCollectivePetscSuccess(communicator_, "flow trial state", VecCopy(committed_state_, state_));
		RequireCollectivePetscSuccess(communicator_, "flow trial history", VecCopy(committed_state_, previous_));
		outlet_models_.swap(outlets);
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
		CollectiveLocalStage(communicator_, "flow rollback preparation", [&] {
			RequirePhase(FlowStepPhase::TrialSolved, "RollbackTrial");
		});
		RestoreCommittedSnapshot();
		phase_ = FlowStepPhase::TrialReady;
	}

	void AbortStep()
	{
		RequireCollectiveSameInt(communicator_, "flow abort phase agreement", static_cast<int>(phase_));
		CollectiveLocalStage(communicator_, "flow abort preparation", [&] {
			if (phase_ != FlowStepPhase::Committed && phase_ != FlowStepPhase::TrialReady
				&& phase_ != FlowStepPhase::TrialSolved && phase_ != FlowStepPhase::CommitPrepared)
				throw std::runtime_error("AbortStep is invalid in the current 3D flow lifecycle phase");
		});
		if (phase_ == FlowStepPhase::Committed) return;
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
		CollectiveLocalStage(communicator_, "flow prepare commit", [&] {
			RequirePhase(FlowStepPhase::TrialSolved, "PrepareCommitStep");
			if (!trial_solve_succeeded_)
				throw std::runtime_error("PrepareCommitStep requires a successful 3D flow trial solve");
		});
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
			CollectiveLocalStage(communicator_, "flow advance boundaries", [&] {
				if (configured_) SetTrialBoundaryConfiguration(step_configuration);
			});
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
		std::vector<CouplingPort> ports;
		CollectiveLocalStage(communicator_, "flow port state preparation", [&] {
			RequirePhase(FlowStepPhase::TrialSolved, "GetPortState");
			if (!trial_solve_succeeded_)
				throw std::runtime_error("GetPortState requires a successful 3D flow trial solve");
			ports.push_back(port);
		});
		const auto measured = MeasurePorts(ports, trial_time_, {}, {});
		PortState result;
		CollectiveLocalStage(communicator_, "flow port state result", [&] { result = measured.at(port.id); });
		return result;
	}

	MPI_Comm Communicator() const noexcept { return communicator_; }

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
		std::vector<double> previous_capacitor_pressure;
		CollectiveLocalStage(communicator_, "flow outlet preparation", [&] {
			if (outlet_models_.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::runtime_error("too many flow outlet models");
			previous_capacitor_pressure.resize(outlet_models_.size());
			for (std::size_t i = 0; i < outlet_models_.size(); ++i)
				previous_capacitor_pressure[i] = outlet_models_[i].capacitor_pressure;
		});
		RequireCollectiveSameInt(communicator_, "flow outlet count agreement", static_cast<int>(outlet_models_.size()));
		bool outlet_converged = false;
		const int maximum_outlet_iterations = outlet_models_.empty() ? 1 : 12;
		for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
			CollectiveLocalStage(communicator_, "flow trial boundaries", [&] {
				if (configured_) UpdateConfiguredBoundaries(trial_configuration_);
			});
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
			std::vector<OutletModelState> candidate;
			CollectiveLocalStage(communicator_, "flow outlet evaluation", [&] {
				candidate = outlet_models_;
				const auto evaluated = EvaluateOutletCoupling(candidate,
					previous_capacitor_pressure, flows, parameters_.dt);
				const auto tolerance = OutletCouplingTolerance(evaluated);
				if (rank_ == 0) std::cout << "step=" << trial_step_+1 << " time=" << trial_time_
					<< " outlet_iteration=" << coupling
					<< " pressure_change=" << evaluated.maximum_pressure_change
					<< " tolerance=" << tolerance << '\n';
				outlet_converged = evaluated.maximum_pressure_change <= tolerance;
				if (outlet_converged) {
					CommitOutletCoupling(candidate, evaluated);
					for (const auto& model : candidate)
						if (rank_ == 0) std::cout << "outlet label=" << model.label
							<< " flow=" << model.flow << " pressure=" << model.pressure
							<< " capacitor_pressure=" << model.capacitor_pressure << '\n';
				} else RelaxOutletCoupling(candidate, evaluated);
			});
			RequireCollectiveSameInt(communicator_, "flow outlet convergence agreement", outlet_converged ? 1 : 0);
			outlet_models_.swap(candidate);
			if (outlet_converged) break;
		}
		if (!outlet_converged)
			throw std::runtime_error("outlet fixed-point iteration did not converge at physical step "
				+std::to_string(trial_step_+1));
	}

public:

	IGA_FLOW_NOINLINE double ReferenceBoundaryFlow(int label) const
	{
		RequireCollectiveSameInt(communicator_, "flow reference label agreement", label);
		double local = 0.0;
		CollectiveLocalStage(communicator_, "flow reference integration", [&] {
			for (const auto& element : owned_elements_) {
				std::vector<std::array<double, 4>> nodal(element.connectivity.size());
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto node = static_cast<std::size_t>(element.connectivity[a]);
					for (int field = 0; field < 3; ++field)
						nodal[a][field] = boundary_velocity_.at(node)[static_cast<std::size_t>(field)];
				}
				BodyFittedSurface4x4QuadratureProvider quadrature(element);
				local += IntegrateBoundaryFlow(element, nodal, quadrature.Rule(), label);
			}
		});
		double global = 0.0;
		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		CollectiveLocalStage(communicator_, "flow reference result", [&] {
			if (!std::isfinite(global)) throw std::runtime_error("nonfinite reference boundary flow");
		});
		return global;
	}

	IGA_FLOW_NOINLINE std::map<std::string, PortState> MeasurePorts(
		const std::vector<CouplingPort>& ports, double physical_time,
		const std::vector<std::string>& species_fields,
		const std::vector<double>& species_state,
		const CompiledLinearSystem* transport_system = nullptr) const
	{
		std::vector<std::vector<double>> advection, diffusion;
		std::map<int, std::size_t> index;
		std::vector<double> local_flow, global_flow, local_pressure, global_pressure, local_area, global_area;
		std::vector<double> local_species, global_species, local_concentration, global_concentration;
		std::string signature;
		CollectiveLocalStage(communicator_, "flow port measurement preparation", [&] {
			RequireFinitePortValue("3D port measurement time_s", physical_time);
			if (ports.empty()) throw std::runtime_error("3D port measurement requires at least one port");
			ValidateCouplingPorts(ports);
			if (!species_fields.empty()
				&& species_state.size() != ghost_nodes_.size()*species_fields.size())
				throw std::runtime_error("3D port species state size does not match flow-required nodes");
			if (transport_system && transport_system->fields != species_fields)
				throw std::runtime_error(
					"3D port transport system fields do not match the supplied species fields");
			advection.assign(species_fields.size(), std::vector<double>(species_fields.size(), 0.0));
			diffusion.assign(species_fields.size(), std::vector<double>(species_fields.size(), 0.0));
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
			std::set<std::string> ids;
			for (std::size_t i = 0; i < ports.size(); ++i) {
				if (!ids.insert(ports[i].id).second)
					throw std::runtime_error("3D port measurement ids must be unique");
				const int label = ParseBoundaryLabelLocator(ports[i]);
				if (!index.emplace(label, i).second)
					throw std::runtime_error("3D port measurement boundary_label locators must be unique");
			}
			const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
			if (index.size() > limit || species_fields.size() > limit/index.size())
				throw std::runtime_error("3D port measurement exceeds MPI count capacity");
			const auto count = index.size()*species_fields.size();
			local_flow.assign(index.size(), 0.0); global_flow.assign(index.size(), 0.0);
			local_pressure.assign(index.size(), 0.0); global_pressure.assign(index.size(), 0.0);
			local_area.assign(index.size(), 0.0); global_area.assign(index.size(), 0.0);
			local_species.assign(count, 0.0); global_species.assign(count, 0.0);
			local_concentration.assign(count, 0.0); global_concentration.assign(count, 0.0);
			// Frame string lengths so user-provided IDs cannot alias the signature.
			const auto append = [&](const std::string& value) {
				signature += std::to_string(value.size())+":"+value;
			};
			append(std::to_string(ports.size()));
			for (const auto& port : ports) {
				append(port.id); append(std::to_string(ParseBoundaryLabelLocator(port)));
				append(std::to_string(port.orientation.native_to_outward_sign));
			}
			append(std::to_string(species_fields.size()));
			for (const auto& field : species_fields) append(field);
			std::ostringstream numbers;
			numbers.exceptions(std::ios::badbit | std::ios::failbit);
			numbers << std::hexfloat << physical_time;
			for (const auto& row : advection) for (double value : row) numbers << ' ' << value;
			for (const auto& row : diffusion) for (double value : row) numbers << ' ' << value;
			append(numbers.str());
		});
		RequireCollectiveSameText(communicator_, "flow port measurement agreement", signature);
		ScatterState();
		CollectiveLocalStage(communicator_, "flow port measurement integration", [&] {
			PetscReadArray view;
			view.Acquire(ghost_state_);
			const auto* values = view.Data();
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
				if (!species_fields.empty()) {
					std::vector<std::vector<double>> species(element.connectivity.size(),
						std::vector<double>(species_fields.size()));
					for (std::size_t a = 0; a < element.connectivity.size(); ++a)
						for (std::size_t field = 0; field < species_fields.size(); ++field)
							species[a][field] = species_state[
								ghost_position_.at(element.connectivity[a])*species_fields.size()+field];
					for (const auto& entry : index) {
						const auto port = entry.second;
						for (std::size_t field = 0; field < species_fields.size(); ++field) {
							const auto measured = IntegrateBoundaryTransportFlux(element, nodal,
								species, field, advection[field], diffusion[field], quadrature.Rule(),
								entry.first);
							local_concentration[port*species_fields.size()+field]
								+= measured.concentration_integral;
							local_species[port*species_fields.size()+field]
								+= measured.total_outward_flux;
						}
					}
				}
			}
			view.Restore();
		});
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
		CollectiveLocalStage(communicator_, "flow port measurement result", [&] {
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
		});
		return result;
	}

	IGA_FLOW_NOINLINE FlowPortMeasurements MeasurePorts(const ThreeDVascularPortDefinition& ports,
		const std::vector<std::string>& species_fields,
		const std::vector<double>& species_state,
		const CompiledLinearSystem* transport_system = nullptr) const
	{
		std::vector<CouplingPort> generic_ports;
		std::vector<int> labels;
		CollectiveLocalStage(communicator_, "flow VCA port preparation", [&] {
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
		});
		const auto measured = MeasurePorts(generic_ports, 0.0, species_fields,
			species_state, transport_system);
		FlowPortMeasurements result;
		CollectiveLocalStage(communicator_, "flow VCA port result", [&] {
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
		});
		return result;
	}

	IGA_FLOW_NOINLINE std::vector<std::array<double, 3>> GatherRequiredVelocity() const
	{
		ScatterState();
		std::vector<std::array<double, 3>> velocity;
		CollectiveLocalStage(communicator_, "flow required velocity", [&] {
			PetscReadArray view;
			view.Acquire(ghost_state_);
			velocity.resize(ghost_nodes_.size());
			for (std::size_t node = 0; node < ghost_nodes_.size(); ++node)
				for (int component = 0; component < 3; ++component)
					velocity[node][static_cast<std::size_t>(component)]
						= PetscRealPart(view.Data()[4*node+component]);
			view.Restore();
		});
		return velocity;
	}

	IGA_FLOW_NOINLINE FlowRuntimeSummary Summary() const
	{
		FlowRuntimeSummary result;
		result.linear_iterations = total_linear_iterations_;
		double local_velocity_squared = 0.0;
		double local_pressure_squared = 0.0;
		CollectiveLocalStage(communicator_, "flow summary local state", [&] {
			PetscInt count = 0;
			if (VecGetLocalSize(state_, &count)) throw std::runtime_error("VecGetLocalSize failed");
			if (count < 0 || static_cast<std::uint64_t>(count)
				!= 4*(assembler_.node_end()-assembler_.node_begin()))
				throw std::runtime_error("flow summary state ownership mismatch");
			PetscReadArray view;
			view.Acquire(state_);
			const auto* values = view.Data();
			for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
				const auto local = static_cast<std::size_t>(node-assembler_.node_begin())*4;
				for (int field = 0; field < 4; ++field)
					if (!std::isfinite(PetscRealPart(values[local+field])))
						throw std::runtime_error("nonfinite flow summary state");
				for (int field = 0; field < 3; ++field) {
					const auto value = PetscRealPart(values[local+field]);
					local_velocity_squared += value*value;
				}
				const auto pressure = PetscRealPart(values[local+3]);
				local_pressure_squared += pressure*pressure;
			}
			view.Restore();
		});
		RequireCollectivePetscSuccess(communicator_, "flow summary state norm", VecNorm(state_, NORM_2, &result.state_l2));
		double velocity_squared = 0.0, pressure_squared = 0.0;
		MPI_Allreduce(&local_velocity_squared, &velocity_squared, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(&local_pressure_squared, &pressure_squared, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		CollectiveLocalStage(communicator_, "flow summary result", [&] {
			result.velocity_l2 = std::sqrt(velocity_squared);
			result.pressure_l2 = std::sqrt(pressure_squared);
			if (!std::isfinite(result.state_l2) || !std::isfinite(result.velocity_l2)
				|| !std::isfinite(result.pressure_l2))
				throw std::runtime_error("nonfinite flow summary norm");
		});
		return result;
	}

	Vec State() const { return state_; }
	const ElementBatchStatistics& LastVolumeBatchStatistics() const noexcept { return last_volume_batch_; }
#ifdef IGA_FLOW_RUNTIME_TESTING
	void SetVolumeProbeForTesting(std::function<void(std::size_t)> probe) { volume_probe_for_testing_=std::move(probe); }
	std::vector<std::array<double, 3>>& ReferenceVelocityForTesting() { return boundary_velocity_; }
#endif
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
		ResolvedBoundaryConditions boundaries;
		std::map<int, double> tractions;
		std::vector<OutletModelState> outlets;
		CollectiveLocalStage(communicator_, "flow restore preparation", [&] {
			boundaries = committed_boundaries_;
			tractions = committed_pressure_tractions_;
			outlets = committed_outlet_models_;
		});
		RequireCollectivePetscSuccess(communicator_, "flow restore state", VecCopy(committed_state_, state_));
		RequireCollectivePetscSuccess(communicator_, "flow restore history", VecCopy(committed_state_, previous_));
		std::swap(boundaries_, boundaries);
		pressure_tractions_.swap(tractions);
		outlet_models_.swap(outlets);
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

	static void SetDefaultPetscOption(const char* name, const char* value,
		std::array<const char*, 4>& added, std::size_t& count)
	{
		PetscBool present = PETSC_FALSE;
		if (PetscOptionsHasName(nullptr, nullptr, name, &present))
			throw std::runtime_error("cannot query default PETSc option");
		if (!present) {
			if (count == added.size()) throw std::logic_error("too many default PETSc options");
			// Record before the call, so a partially failed insertion is also
			// rolled back. Names are static strings, requiring no allocation.
			added[count++] = name;
			if (PetscOptionsSetValue(nullptr, name, value))
				throw std::runtime_error("cannot set default PETSc option");
			ProbeRuntimeConstruction(name);
		}
	}

	void BuildBoundaryLabelIndex()
	{
		int size = 0;
		MPI_Comm_size(communicator_, &size);
		std::vector<int> counts, offsets, local_labels, gathered;
		RuntimeConstructionStage(communicator_, "flow boundary catalog preparation", [&] {
			std::set<int> local;
			for (const auto& element : owned_elements_)
				for (const auto label : element.boundary_labels)
					if (label >= 0) local.insert(label);
			if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::overflow_error("flow boundary label count exceeds MPI limits");
			local_labels.assign(local.begin(), local.end());
			counts.resize(static_cast<std::size_t>(size));
			offsets.resize(static_cast<std::size_t>(size));
		});
		const int local_count = static_cast<int>(local_labels.size());
		MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator_);
		RuntimeConstructionStage(communicator_, "flow boundary catalog layout", [&] {
			int total = 0;
			for (int i = 0; i < size; ++i) {
				if (counts[i] < 0 || counts[i] > std::numeric_limits<int>::max() - total)
					throw std::overflow_error("flow boundary label gather exceeds MPI limits");
				offsets[i] = total;
				total += counts[i];
			}
			gathered.resize(static_cast<std::size_t>(total));
		});
		MPI_Allgatherv(local_labels.data(), local_count, MPI_INT, gathered.data(), counts.data(),
			offsets.data(), MPI_INT, communicator_);
		RuntimeConstructionStage(communicator_, "flow boundary catalog publication", [&] {
			std::sort(gathered.begin(), gathered.end());
			gathered.erase(std::unique(gathered.begin(), gathered.end()), gathered.end());
			boundary_labels_ = std::move(gathered);
			for (std::size_t i = 0; i < boundary_labels_.size(); ++i)
				boundary_label_index_.emplace(boundary_labels_[i], i);
		});
	}

	FlowConvergenceMetrics MeasureConvergence(const PetscScalar* ghost_values) const
	{
		FlowConvergenceMetrics result;
		double local_continuity_squared = 0.0, local_continuity_sum = 0.0;
		CollectiveLocalStage(communicator_, "flow continuity diagnostics", [&] {
			PetscReadArray residual_view;
			residual_view.Acquire(rhs_);
			const auto* residual_values = residual_view.Data();
			for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node) {
				const auto local = static_cast<std::size_t>(node-assembler_.node_begin())*4+3;
				const auto value = PetscRealPart(residual_values[local]);
				local_continuity_squared += value*value;
				local_continuity_sum += value;
			}
			residual_view.Restore();
		});
		double global_continuity_squared = 0.0;
		MPI_Allreduce(&local_continuity_squared, &global_continuity_squared, 1,
			MPI_DOUBLE, MPI_SUM, communicator_);
		MPI_Allreduce(&local_continuity_sum, &result.continuity_sum, 1,
			MPI_DOUBLE, MPI_SUM, communicator_);
		result.continuity_l2 = std::sqrt(global_continuity_squared);

		std::vector<double> local_flow, global_flow;
		CollectiveLocalStage(communicator_, "flow boundary integration", [&] {
			local_flow.assign(boundary_labels_.size(), 0.0);
			global_flow.assign(boundary_labels_.size(), 0.0);
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
		});
		if (!local_flow.empty())
			MPI_Allreduce(local_flow.data(), global_flow.data(), static_cast<int>(local_flow.size()),
				MPI_DOUBLE, MPI_SUM, communicator_);
		CollectiveLocalStage(communicator_, "flow conservation result", [&] {
			for (const auto flow : global_flow) {
				result.net_boundary_flow += flow;
				result.absolute_boundary_flow += std::abs(flow);
			}
			if (result.absolute_boundary_flow > 0.0)
				result.relative_mass_imbalance = 2.0*std::abs(result.net_boundary_flow)
					/result.absolute_boundary_flow;
			for (const auto value : {result.continuity_l2, result.continuity_sum,
				result.net_boundary_flow, result.absolute_boundary_flow, result.relative_mass_imbalance})
				if (!std::isfinite(value)) throw std::runtime_error("nonfinite flow conservation diagnostics");
		});
		return result;
	}

	static std::string PreciseNumber(double value)
	{
		std::ostringstream stream;
		stream.exceptions(std::ios::badbit | std::ios::failbit);
		stream << std::scientific << std::setprecision(16) << value;
		return stream.str();
	}

	void BuildGhostScatter()
	{
		std::vector<PetscInt> rows;
		RuntimeConstructionStage(communicator_, "flow halo preparation", [&] {
			std::vector<std::int32_t> nodes;
			for (const auto& element : assembler_.elements())
				nodes.insert(nodes.end(), element.connectivity.begin(), element.connectivity.end());
			for (const auto& element : owned_elements_)
				nodes.insert(nodes.end(), element.connectivity.begin(), element.connectivity.end());
			std::sort(nodes.begin(), nodes.end());
			nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
			ghost_nodes_ = std::move(nodes);
			rows = assembler_.RequiredRows(ghost_nodes_);
			for (std::size_t i = 0; i < ghost_nodes_.size(); ++i)
				ghost_position_.emplace(ghost_nodes_[i], i);
		});
		RequireCollectivePetscSuccess(communicator_, "flow halo source index creation",
			ISCreateGeneral(communicator_, static_cast<PetscInt>(rows.size()), rows.data(), PETSC_COPY_VALUES, &source_rows_));
		ObserveConstructedObject(communicator_, "flow halo source index created", reinterpret_cast<PetscObject>(source_rows_));
		RequireCollectivePetscSuccess(communicator_, "flow halo state creation",
			VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), &ghost_state_));
		ObserveConstructedObject(communicator_, "flow halo state created", reinterpret_cast<PetscObject>(ghost_state_));
		RequireCollectivePetscSuccess(communicator_, "flow halo previous creation", VecDuplicate(ghost_state_, &ghost_previous_));
		ObserveConstructedObject(communicator_, "flow halo previous created", reinterpret_cast<PetscObject>(ghost_previous_));
		RequireCollectivePetscSuccess(communicator_, "flow halo destination index creation",
			ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), 0, 1, &destination_rows_));
		ObserveConstructedObject(communicator_, "flow halo destination index created", reinterpret_cast<PetscObject>(destination_rows_));
		RequireCollectivePetscSuccess(communicator_, "flow halo scatter creation",
			VecScatterCreate(state_, source_rows_, ghost_state_, destination_rows_, &scatter_));
		ObserveConstructedObject(communicator_, "flow halo scatter created", reinterpret_cast<PetscObject>(scatter_));
	}

	void ScatterState() const
	{
		PhaseScope communication_phase(ProfilePhase::Communication);
		RequireCollectivePetscSuccess(communicator_, "flow state scatter begin", VecScatterBegin(scatter_, state_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD));
		RequireCollectivePetscSuccess(communicator_, "flow state scatter end", VecScatterEnd(scatter_, state_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD));
	}

	double BoundaryValue(PetscInt row) const
	{
		return BoundaryValue(row, boundaries_);
	}

	double BoundaryValue(PetscInt row, const ResolvedBoundaryConditions& boundaries) const
	{
		if (row < 0) throw std::out_of_range("negative flow boundary row");
		const auto node = static_cast<std::uint64_t>(row/4);
		const auto field = row%4;
		const auto index = static_cast<std::size_t>(node);
		return field < 3 ? (wall_trace_basis_.count(static_cast<std::int32_t>(node)) ? 0.0
			: boundaries.velocity.at(index)[static_cast<std::size_t>(field)])
			: boundaries.pressure.at(index);
	}

	void ResolveConfiguredBoundaries(const SimulationConfiguration& step_configuration,
		ResolvedBoundaryConditions& boundaries, std::map<int, double>& pressure_tractions) const
	{
		const auto configuration = MaterializeOutletPressures(step_configuration, outlet_models_);
		const auto& system = FirstNavierStokesSystem(configuration);
		pressure_tractions = ExtractPressureTractions(configuration, system);
		boundaries = ResolveFlowBoundaries(configuration, system, labels_, boundary_velocity_);
		ValidateConstraintTopology(boundaries);
		for (const auto& override : trial_pressure_overrides_)
			pressure_tractions[override.first] = override.second;
	}

	void UpdateConfiguredBoundaries(const SimulationConfiguration& step_configuration)
	{
		ResolvedBoundaryConditions boundaries;
		std::map<int, double> pressure_tractions;
		ResolveConfiguredBoundaries(step_configuration, boundaries, pressure_tractions);
		pressure_tractions_.swap(pressure_tractions);
		std::swap(boundaries_, boundaries);
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
		RequireCollectiveSameInt(communicator_, "flow transient mode agreement", transient_ ? 1 : 0);
		PetscReal initial_residual = -1.0;
		for (int nonlinear = 0; nonlinear < maximum_newton; ++nonlinear) {
			PhaseScope assembly_phase(ProfilePhase::Assembly);
			const auto iteration_start = std::chrono::steady_clock::now();
			// All ranks set this flag before entering the local assembly stage,
			// and clear it only after collectively completing both objects.
			// A coordinated worker/insertion failure therefore leaves every
			// rank on the same recovery branch, including empty ranks.
			if(assembly_pending_) {
				// FLUSH clears insertion mode without compressing preallocated
				// entries that a first, failed assembly has not filled yet.
				RequireCollectivePetscSuccess(communicator_, "flow recovery flush begin",
					MatAssemblyBegin(jacobian_, MAT_FLUSH_ASSEMBLY));
				RequireCollectivePetscSuccess(communicator_, "flow recovery flush end",
					MatAssemblyEnd(jacobian_, MAT_FLUSH_ASSEMBLY));
				OwnedRowAssembler::Assemble(rhs_, communicator_);
			}
			assembly_pending_=true;
			RequireCollectivePetscSuccess(communicator_, "flow clear Jacobian", MatZeroEntries(jacobian_));
			RequireCollectivePetscSuccess(communicator_, "flow clear residual", VecSet(rhs_, 0.0));
			ScatterState();
			if (transient_) {
				PhaseScope communication_phase(ProfilePhase::Communication);
				RequireCollectivePetscSuccess(communicator_, "flow history scatter begin", VecScatterBegin(scatter_, previous_, ghost_previous_, INSERT_VALUES, SCATTER_FORWARD));
				RequireCollectivePetscSuccess(communicator_, "flow history scatter end", VecScatterEnd(scatter_, previous_, ghost_previous_, INSERT_VALUES, SCATTER_FORWARD));
			}
			PetscReadArray current_view, previous_view;
			CollectiveLocalStage(communicator_, "flow element assembly", [&] {
				current_view.Acquire(ghost_state_);
				const auto* values = current_view.Data();
				if (transient_) previous_view.Acquire(ghost_previous_);
				const auto* previous_values = previous_view.Data();
				struct PreparedVolume {
					const Element* element;
					std::vector<std::array<double,4>> nodal, previous_nodal;
				};
				const auto batch=ForEachElementBatch(assembler_.elements().size(),assembly_execution_->Options(),
					[&](std::size_t index) {
					const auto& element=assembler_.elements()[index];
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
					return PreparedVolume{&element,std::move(nodal),std::move(previous_nodal)};
					},
					[this](const PreparedVolume& input,std::size_t index) {
					FullCell4x4x4VolumeQuadratureProvider volume_quadrature(*input.element);
					auto local=BuildNavierStokesElement(*input.element,input.nodal,input.previous_nodal,
						parameters_,volume_quadrature.Rule());
#ifdef IGA_FLOW_RUNTIME_TESTING
					if(volume_probe_for_testing_) volume_probe_for_testing_(index);
#else
					(void)index;
#endif
					return local;
					},
					[&](const PreparedVolume& input,NavierStokesSystem& local,std::size_t) {
					const auto& element=*input.element;
					BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
					for (const auto& traction : pressure_tractions_) {
						const auto surface = IntegrateBoundaryPressureTraction(element,
							surface_quadrature.Rule(), traction.first, traction.second);
						for (std::size_t row = 0; row < surface.size(); ++row)
							local.negative_residual[row] += surface[row];
					}
					assembler_.AddElementMatrix(jacobian_, element, local.jacobian);
					assembler_.AddElementVector(rhs_, element, local.negative_residual);
					});
				last_volume_batch_=batch;
				if(CurrentPhaseProfile().Enabled())
					std::cout << "body_fitted_element_assembly rank=" << rank_
						<< " threads_requested=" << assembly_execution_->Options().threads
						<< " team_size=" << batch.maximum_team_size << " elements=" << batch.items
						<< " batches=" << batch.batches << " maximum_resident_items=" << batch.maximum_resident_items << '\n';
				previous_view.Restore();
			});
			{
				PhaseScope communication_phase(ProfilePhase::Communication);
				OwnedRowAssembler::Assemble(jacobian_, communicator_);
				OwnedRowAssembler::Assemble(rhs_, communicator_);
			}
			assembly_pending_=false;
			PhaseScope diagnostics_phase(ProfilePhase::Diagnostics);
			const auto convergence = MeasureConvergence(current_view.Data());
			diagnostics_phase.Stop();
			std::vector<PetscScalar> boundary_update;
			CollectiveLocalStage(communicator_, "flow boundary values", [&] {
				current_view.Restore();
				PetscReadArray owned_view;
				owned_view.Acquire(state_);
				const auto* owned = owned_view.Data();
				boundary_update.reserve(boundary_rows_.size());
				for (const auto row : boundary_rows_) {
					const auto local = static_cast<std::size_t>(row-4*assembler_.node_begin());
					boundary_update.push_back(BoundaryValue(row)-PetscRealPart(owned[local]));
				}
				owned_view.Restore();
			});
			RequireCollectivePetscSuccess(communicator_, "flow boundary rows",
				MatZeroRows(jacobian_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
					1.0, nullptr, nullptr));
			CollectiveLocalStage(communicator_, "flow boundary insertion", [&] {
				if (VecSetValues(rhs_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
					boundary_update.data(), INSERT_VALUES))
					throw std::runtime_error("flow boundary VecSetValues failed");
			});
			OwnedRowAssembler::Assemble(rhs_, communicator_);
			PetscReal residual = 0.0;
			RequireCollectivePetscSuccess(communicator_, "flow residual norm", VecNorm(rhs_, NORM_2, &residual));
			if (initial_residual < 0.0) initial_residual = residual;
			const auto residual_scale = std::sqrt(
				static_cast<PetscReal>(assembler_.global_rows()));
			const auto residual_rms = residual/residual_scale;
			const auto tolerance = std::max<PetscReal>(
				nonlinear_absolute_tolerance*residual_scale,
				nonlinear_relative_tolerance*initial_residual);
			bool converged = false;
			CollectiveLocalStage(communicator_, "flow nonlinear convergence", [&] {
				if (!std::isfinite(residual) || !std::isfinite(tolerance) || !std::isfinite(residual_rms))
					throw std::runtime_error("nonfinite Navier-Stokes nonlinear residual/tolerance");
				converged = residual <= tolerance
					&& convergence.relative_mass_imbalance <= mass_relative_tolerance;
			});
			RequireCollectiveSameInt(communicator_, "flow nonlinear convergence agreement", converged ? 1 : 0);
			if (converged) {
				CollectiveLocalStage(communicator_, "flow convergence logging", [&] {
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
				});
				return true;
			}
			const auto linear_start = std::chrono::steady_clock::now();
			assembly_phase.Stop();
			RequireKspFactorBackend(solver_, jacobian_, communicator_);
			RequireCollectivePetscSuccess(communicator_, "flow solver operators", KSPSetOperators(solver_, jacobian_, jacobian_));
			{
				PhaseScope setup_phase(ProfilePhase::SolverSetup);
				RequireCollectivePetscSuccess(communicator_, "flow solver setup", KSPSetUp(solver_));
				RequireCollectivePetscSuccess(communicator_, "flow block solver setup", KSPSetUpOnBlocks(solver_));
			}
			{
				PhaseScope solve_phase(ProfilePhase::LinearSolve);
				RequireCollectivePetscSuccess(communicator_, "flow linear solve", KSPSolve(solver_, rhs_, update_));
			}
			PetscInt iterations = 0;
			PetscReal linear_residual = 0.0, update_norm = 0.0;
			CollectiveLocalStage(communicator_, "flow linear convergence", [&] {
				KSPConvergedReason reason;
				if (KSPGetConvergedReason(solver_, &reason)
					|| KSPGetIterationNumber(solver_, &iterations)
					|| KSPGetResidualNorm(solver_, &linear_residual))
					throw std::runtime_error("cannot query Navier-Stokes linear convergence");
				if (reason <= 0)
					throw std::runtime_error("Navier-Stokes linear solve failed at nonlinear iteration "
						+std::to_string(nonlinear) + " (KSP reason="
						+std::to_string(static_cast<int>(reason)) + ", iterations="
						+std::to_string(static_cast<long long>(iterations)) + ", residual="
						+PreciseNumber(static_cast<double>(linear_residual)) + ", nonlinear_residual="
						+PreciseNumber(static_cast<double>(residual)) + ", continuity_l2="
						+PreciseNumber(convergence.continuity_l2) + ", net_boundary_flow="
						+PreciseNumber(convergence.net_boundary_flow) + ", relative_mass_imbalance="
						+PreciseNumber(convergence.relative_mass_imbalance) + ")");
				if (!std::isfinite(linear_residual) || iterations < 0
					|| iterations > std::numeric_limits<PetscInt>::max()-trial_linear_iterations_)
					throw std::runtime_error("invalid Navier-Stokes linear diagnostics");
			});
			RequireCollectivePetscSuccess(communicator_, "flow update norm", VecNorm(update_, NORM_2, &update_norm));
			CollectiveLocalStage(communicator_, "flow update validation", [&] {
				if (!std::isfinite(update_norm)) throw std::runtime_error("nonfinite Navier-Stokes update");
			});
			RequireCollectivePetscSuccess(communicator_, "flow update state", VecAXPY(state_, 1.0, update_));
			trial_linear_iterations_ += iterations;
			CollectiveLocalStage(communicator_, "flow iteration logging", [&] {
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
			});
		}
		return false;
	}

	IGA_FLOW_NOINLINE std::vector<double> MeasureOutletFlows() const
	{
		std::vector<double> local, global;
		std::map<int, std::size_t> indices;
		std::string labels;
		CollectiveLocalStage(communicator_, "flow outlet measurement preparation", [&] {
			if (outlet_models_.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::runtime_error("too many flow outlet models");
			local.assign(outlet_models_.size(), 0.0); global.resize(local.size());
			for (std::size_t i = 0; i < outlet_models_.size(); ++i) {
				indices.emplace(outlet_models_[i].label, i);
				labels += std::to_string(outlet_models_[i].label)+" ";
			}
		});
		RequireCollectiveSameText(communicator_, "flow outlet labels agreement", labels);
		ScatterState();
		CollectiveLocalStage(communicator_, "flow outlet integration", [&] {
			PetscReadArray view;
			view.Acquire(ghost_state_);
			const auto* values = view.Data();
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
			view.Restore();
		});
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()),
			MPI_DOUBLE, MPI_SUM, communicator_);
		return global;
	}

	std::optional<ElementAssemblyExecution> assembly_execution_;
	ElementBatchStatistics last_volume_batch_;
	bool assembly_pending_=false;
#ifdef IGA_FLOW_RUNTIME_TESTING
	std::function<void(std::size_t)> volume_probe_for_testing_;
#endif
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

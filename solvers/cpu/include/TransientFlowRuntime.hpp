#ifndef IGA_TRANSIENT_FLOW_RUNTIME_HPP
#define IGA_TRANSIENT_FLOW_RUNTIME_HPP

#include "BoundaryFlow.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "NavierStokesElement.hpp"
#include "OutletModel.hpp"
#include "OwnedRowAssembler.hpp"
#include "PressureTraction.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
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
		MPI_Comm_rank(communicator_, &rank_);
		owned_elements_ = database_.LoadOwned(rank_);
		BuildBoundaryLabelIndex();
		jacobian_ = assembler_.CreateMatrix(true);
		state_ = assembler_.CreateVector();
		previous_ = assembler_.CreateVector();
		update_ = assembler_.CreateVector();
		rhs_ = assembler_.CreateVector();
		VecSet(state_, 0.0);
		VecSet(previous_, 0.0);
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
		VecCopy(state_, previous_);
	}

	IGA_FLOW_NOINLINE void Advance(const SimulationConfiguration& step_configuration, int step,
		double physical_time, int maximum_newton, double nonlinear_relative_tolerance,
		double nonlinear_absolute_tolerance, double mass_relative_tolerance)
	{
		if (maximum_newton <= 0) throw std::invalid_argument("maximum Newton iterations must be positive");
		if (!std::isfinite(nonlinear_relative_tolerance) || nonlinear_relative_tolerance <= 0.0)
			throw std::invalid_argument("nonlinear relative tolerance must be finite and positive");
		if (!std::isfinite(nonlinear_absolute_tolerance) || nonlinear_absolute_tolerance <= 0.0)
			throw std::invalid_argument("nonlinear absolute tolerance must be finite and positive");
		if (!std::isfinite(mass_relative_tolerance) || mass_relative_tolerance <= 0.0)
			throw std::invalid_argument("mass relative tolerance must be finite and positive");
		if (step > 0) VecCopy(state_, previous_);
		std::vector<double> previous_capacitor_pressure(outlet_models_.size());
		for (std::size_t i = 0; i < outlet_models_.size(); ++i)
			previous_capacitor_pressure[i] = outlet_models_[i].capacitor_pressure;
		bool outlet_converged = false;
		const int maximum_outlet_iterations = outlet_models_.empty() ? 1 : 12;
		for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
			if (configured_) UpdateConfiguredBoundaries(step_configuration);
			const auto converged = SolveNonlinearStep(step, physical_time, maximum_newton,
				nonlinear_relative_tolerance, nonlinear_absolute_tolerance,
				mass_relative_tolerance);
			if (!converged)
				throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON at physical step "
					+std::to_string(step+1));
			if (outlet_models_.empty()) {
				outlet_converged = true;
				break;
			}
			const auto flows = MeasureOutletFlows();
			const auto evaluated = EvaluateOutletCoupling(outlet_models_,
				previous_capacitor_pressure, flows, parameters_.dt);
			const auto tolerance = OutletCouplingTolerance(evaluated);
			if (rank_ == 0) std::cout << "step=" << step+1 << " time=" << physical_time
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
				+std::to_string(step+1));
	}

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
			for (std::size_t face = 0; face < element.boundary_labels.size(); ++face)
				if (element.boundary_labels[face] == label)
					local += IntegrateBoundaryFlow(element, face, nodal);
		}
		double global = 0.0;
		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, communicator_);
		return global;
	}

	IGA_FLOW_NOINLINE FlowPortMeasurements MeasurePorts(const ThreeDVascularPortDefinition& ports,
		const std::vector<std::string>& species_fields,
		const std::vector<double>& species_state) const
	{
		if (!species_fields.empty()
			&& species_state.size() != ghost_nodes_.size()*species_fields.size())
			throw std::runtime_error("VCA species state size does not match flow-required nodes");
		std::map<int, std::size_t> index;
		for (std::size_t i = 0; i < ports.outlet_labels.size(); ++i)
			index.emplace(ports.outlet_labels[i], i);
		std::vector<double> local_flow(index.size(), 0.0), global_flow(index.size(), 0.0);
		std::vector<double> local_pressure(index.size(), 0.0), global_pressure(index.size(), 0.0);
		std::vector<double> local_area(index.size(), 0.0), global_area(index.size(), 0.0);
		std::vector<double> local_species(index.size()*species_fields.size(), 0.0);
		std::vector<double> global_species(index.size()*species_fields.size(), 0.0);
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
			for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
				const auto found = index.find(element.boundary_labels[face]);
				if (found == index.end()) continue;
				const auto port = found->second;
				local_flow[port] += IntegrateBoundaryFlow(element, face, nodal);
				const auto pressure = IntegrateBoundaryScalarAndArea(element, face, nodal);
				local_pressure[port] += pressure[0];
				local_area[port] += pressure[1];
				for (std::size_t field = 0; field < species_fields.size(); ++field) {
					std::vector<double> species(element.connectivity.size());
					for (std::size_t a = 0; a < element.connectivity.size(); ++a)
						species[a] = species_state[ghost_position_.at(element.connectivity[a])
							*species_fields.size()+field];
					local_species[port*species_fields.size()+field]
						+= IntegrateBoundarySpeciesFlux(element, face, nodal, species);
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
		FlowPortMeasurements result;
		for (std::size_t i = 0; i < ports.outlet_labels.size(); ++i) {
			if (!(global_area[i] > 0.0)) throw std::runtime_error("VCA outlet has zero boundary area");
			result.flows.emplace(ports.outlet_labels[i], global_flow[i]);
			result.pressures.emplace(ports.outlet_labels[i], global_pressure[i]/global_area[i]);
			for (std::size_t field = 0; field < species_fields.size(); ++field)
				result.species_fluxes[ports.outlet_labels[i]].emplace(species_fields[field],
					global_species[i*species_fields.size()+field]);
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
	std::vector<OutletModelState>& OutletModels() { return outlet_models_; }
	const std::vector<OutletModelState>& OutletModels() const { return outlet_models_; }

private:
	static constexpr std::uint64_t kScalablePreconditionerNodeThreshold = 1000;

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
			for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
				const auto found = boundary_label_index_.find(element.boundary_labels[face]);
				if (found != boundary_label_index_.end())
					local_flow[found->second] += IntegrateBoundaryFlow(element, face, nodal);
			}
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
		pressure_tractions_ = ExtractPressureTractions(configuration, FirstNavierStokesSystem(configuration));
		boundaries_ = ResolveFlowBoundaries(configuration, FirstNavierStokesSystem(configuration),
			labels_, boundary_velocity_);
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
				auto local = BuildNavierStokesElement(element, nodal, previous_nodal, parameters_);
				for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
					const auto traction = pressure_tractions_.find(element.boundary_labels[face]);
					if (traction == pressure_tractions_.end()) continue;
					const auto surface = IntegrateBoundaryPressureTraction(element, face, traction->second);
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
			total_linear_iterations_ += iterations;
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
			for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
				const auto found = indices.find(element.boundary_labels[face]);
				if (found != indices.end())
					local[found->second] += IntegrateBoundaryFlow(element, face, nodal);
			}
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
	Vec state_ = nullptr, previous_ = nullptr, update_ = nullptr, rhs_ = nullptr;
	IS source_rows_ = nullptr, destination_rows_ = nullptr;
	Vec ghost_state_ = nullptr, ghost_previous_ = nullptr;
	VecScatter scatter_ = nullptr;
	KSP solver_ = nullptr;
	PetscInt total_linear_iterations_ = 0;
	int rank_ = 0;
};

} // namespace iga

#undef IGA_FLOW_NOINLINE

#endif

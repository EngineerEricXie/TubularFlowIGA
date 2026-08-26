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
#include <iostream>
#include <map>
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

class TransientFlowRuntime {
public:
	TransientFlowRuntime(Database& database, MPI_Comm communicator, bool configured,
		bool transient, NavierStokesParameters parameters,
		ResolvedBoundaryConditions initial_boundaries,
		std::vector<int> labels,
		std::vector<std::array<double, 3>> boundary_velocity,
		std::vector<OutletModelState> outlet_models)
		: database_(database), communicator_(communicator), configured_(configured),
			transient_(transient), parameters_(parameters), boundaries_(std::move(initial_boundaries)),
			labels_(std::move(labels)), boundary_velocity_(std::move(boundary_velocity)),
			outlet_models_(std::move(outlet_models)), assembler_(database, communicator, 4)
	{
		if (labels_.size() != database_.header().nodes
			|| boundary_velocity_.size() != database_.header().nodes)
			throw std::runtime_error("flow boundary data do not match database nodes");
		MPI_Comm_rank(communicator_, &rank_);
		owned_elements_ = database_.LoadOwned(rank_);
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
			if (boundaries_.velocity_constrained[index])
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
		PCSetType(preconditioner, PCBJACOBI);
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
		double physical_time, int maximum_newton)
	{
		if (maximum_newton <= 0) throw std::invalid_argument("maximum Newton iterations must be positive");
		if (step > 0) VecCopy(state_, previous_);
		std::vector<double> previous_capacitor_pressure(outlet_models_.size());
		for (std::size_t i = 0; i < outlet_models_.size(); ++i)
			previous_capacitor_pressure[i] = outlet_models_[i].capacitor_pressure;
		bool outlet_converged = false;
		const int maximum_outlet_iterations = outlet_models_.empty() ? 1 : 12;
		for (int coupling = 0; coupling < maximum_outlet_iterations; ++coupling) {
			if (configured_) UpdateConfiguredBoundaries(step_configuration);
			const auto converged = SolveNonlinearStep(step, physical_time, maximum_newton);
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
	const std::vector<std::int32_t>& RequiredNodes() const { return ghost_nodes_; }
	std::vector<OutletModelState>& OutletModels() { return outlet_models_; }
	const std::vector<OutletModelState>& OutletModels() const { return outlet_models_; }

private:
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
		return field < 3 ? boundaries_.velocity[index][static_cast<std::size_t>(field)]
			: boundaries_.pressure[index];
	}

	void UpdateConfiguredBoundaries(const SimulationConfiguration& step_configuration)
	{
		const auto configuration = MaterializeOutletPressures(step_configuration, outlet_models_);
		pressure_tractions_ = ExtractPressureTractions(configuration, FirstNavierStokesSystem(configuration));
		boundaries_ = ResolveFlowBoundaries(configuration, FirstNavierStokesSystem(configuration),
			labels_, boundary_velocity_);
	}

	IGA_FLOW_NOINLINE bool SolveNonlinearStep(int step, double physical_time, int maximum_newton)
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
			VecRestoreArrayRead(ghost_state_, &values);
			OwnedRowAssembler::Assemble(jacobian_);
			OwnedRowAssembler::Assemble(rhs_);
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
			const auto tolerance = std::max<PetscReal>(1e-10, 1e-5*initial_residual);
			if (residual <= tolerance) {
				if (rank_ == 0) std::cout << "step=" << step+1 << " time=" << physical_time
					<< " converged newton=" << nonlinear << " residual_l2=" << residual
					<< " tolerance=" << tolerance << " assembly_s="
					<< std::chrono::duration<double>(std::chrono::steady_clock::now()-iteration_start).count() << '\n';
				return true;
			}
			const auto linear_start = std::chrono::steady_clock::now();
			KSPSetOperators(solver_, jacobian_, jacobian_);
			KSPSolve(solver_, rhs_, update_);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver_, &reason);
			if (reason <= 0)
				throw std::runtime_error("Navier-Stokes linear solve failed at nonlinear iteration "
					+std::to_string(nonlinear));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver_, &iterations);
			total_linear_iterations_ += iterations;
			PetscReal linear_residual = 0.0, update_norm = 0.0;
			KSPGetResidualNorm(solver_, &linear_residual);
			VecNorm(update_, NORM_2, &update_norm);
			VecAXPY(state_, 1.0, update_);
			if (rank_ == 0) std::cout << "step=" << step+1 << " time=" << physical_time
				<< " newton=" << nonlinear << " residual_l2=" << residual
				<< " update_l2=" << update_norm << " linear_iterations=" << iterations
				<< " linear_residual=" << linear_residual << " assembly_s="
				<< std::chrono::duration<double>(linear_start-iteration_start).count()
				<< " linear_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-linear_start).count() << '\n';
			if (update_norm < 1e-10) return true;
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
	std::vector<OutletModelState> outlet_models_;
	std::vector<Element> owned_elements_;
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

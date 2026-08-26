#ifndef IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP
#define IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP

#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "OwnedRowAssembler.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

class TransientTransportRuntime {
public:
	TransientTransportRuntime(Database& database, MPI_Comm communicator,
		const SimulationConfiguration& configuration,
		CompiledLinearSystem system, const std::vector<int>& labels)
		: communicator_(communicator), configuration_(configuration),
			system_(std::move(system)), assembler_(database, communicator, system_.fields.size()),
			labels_(labels)
	{
		MPI_Comm_rank(communicator_, &rank_);
		if (system_.velocity_source != "prescribed")
			throw std::runtime_error("in-process VCA transport requires velocity_source prescribed");
		if (labels_.size() != database.header().nodes)
			throw std::runtime_error("transport boundary labels do not match database nodes");
		const auto boundaries = ResolveScalarBoundaries(configuration_, system_, labels_);
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node)
			for (std::size_t field = 0; field < system_.fields.size(); ++field)
				if (boundaries.constrained[static_cast<std::size_t>(node)*system_.fields.size()+field])
					boundary_rows_.push_back(static_cast<PetscInt>(node*system_.fields.size()+field));
		left_ = assembler_.CreateMatrix();
		previous_ = assembler_.CreateMatrix();
		forcing_ = assembler_.CreateVector();
		current_ = assembler_.CreateVector();
		next_ = assembler_.CreateVector();
		rhs_ = assembler_.CreateVector();
		MatSetOption(left_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
		MatSetOption(previous_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
		MatSetOption(previous_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
		const auto initial = InitialScalarValues(configuration_, system_);
		PetscScalar* values = nullptr;
		VecGetArray(current_, &values);
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node)
			for (std::size_t field = 0; field < system_.fields.size(); ++field) {
				const auto local = static_cast<std::size_t>(node-assembler_.node_begin())
					*system_.fields.size()+field;
				const auto global = static_cast<std::size_t>(node)*system_.fields.size()+field;
				values[local] = boundaries.constrained[global] ? boundaries.value[global] : initial[field];
			}
		VecRestoreArray(current_, &values);
		BuildGhostScatter();
		KSPCreate(communicator_, &solver_);
		KSPSetType(solver_, KSPGMRES);
		KSPGMRESSetRestart(solver_, 50);
		KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
		PC preconditioner = nullptr;
		KSPGetPC(solver_, &preconditioner);
		PCSetType(preconditioner, PCBJACOBI);
		KSPSetFromOptions(solver_);
	}

	~TransientTransportRuntime()
	{
		KSPDestroy(&solver_);
		VecScatterDestroy(&scatter_);
		ISDestroy(&destination_rows_);
		VecDestroy(&ghost_state_);
		ISDestroy(&source_rows_);
		VecDestroy(&rhs_); VecDestroy(&next_); VecDestroy(&current_); VecDestroy(&forcing_);
		MatDestroy(&previous_); MatDestroy(&left_);
	}

	void Advance(const SimulationConfiguration& step_configuration,
		const std::vector<std::int32_t>& velocity_nodes,
		const std::vector<std::array<double, 3>>& velocity)
	{
		if (velocity_nodes != ghost_nodes_ || velocity.size() != ghost_nodes_.size())
			throw std::runtime_error("VCA transport velocity nodes do not match required transport nodes");
		const auto boundaries = ResolveScalarBoundaries(step_configuration, system_, labels_);
		MatZeroEntries(left_);
		MatZeroEntries(previous_);
		VecSet(forcing_, 0.0);
		for (const auto& element : assembler_.elements()) {
			const auto matrices = BuildGenericTransportElementWithVelocity(element,
				[&velocity, this](std::int32_t node) -> const std::array<double, 3>& {
					return velocity.at(ghost_position_.at(node));
				}, system_, step_configuration);
			assembler_.AddElementMatrix(left_, element, matrices.left);
			assembler_.AddElementMatrix(previous_, element, matrices.previous);
			assembler_.AddElementVector(forcing_, element, matrices.source);
		}
		OwnedRowAssembler::Assemble(left_);
		OwnedRowAssembler::Assemble(previous_);
		OwnedRowAssembler::Assemble(forcing_);
		MatZeroRows(left_, static_cast<PetscInt>(boundary_rows_.size()),
			boundary_rows_.data(), 1.0, nullptr, nullptr);
		MatZeroRows(previous_, static_cast<PetscInt>(boundary_rows_.size()),
			boundary_rows_.data(), 0.0, nullptr, nullptr);
		MatMult(previous_, current_, rhs_);
		VecAXPY(rhs_, 1.0, forcing_);
		std::vector<PetscScalar> boundary_values;
		boundary_values.reserve(boundary_rows_.size());
		for (const auto row : boundary_rows_)
			boundary_values.push_back(boundaries.value[static_cast<std::size_t>(row)]);
		VecSetValues(rhs_, static_cast<PetscInt>(boundary_rows_.size()), boundary_rows_.data(),
			boundary_values.data(), INSERT_VALUES);
		OwnedRowAssembler::Assemble(rhs_);
		if (steps_ > 0) VecCopy(current_, next_);
		KSPSetOperators(solver_, left_, left_);
		KSPSetInitialGuessNonzero(solver_, steps_ > 0 ? PETSC_TRUE : PETSC_FALSE);
		KSPSetUp(solver_);
		KSPSolve(solver_, rhs_, next_);
		KSPConvergedReason reason;
		KSPGetConvergedReason(solver_, &reason);
		if (reason <= 0) throw std::runtime_error("VCA transport linear solve did not converge");
		VecSwap(current_, next_);
		++steps_;
	}

	std::vector<double> GatherState() const
	{
		Vec all = nullptr;
		VecScatter scatter = nullptr;
		VecScatterCreateToAll(current_, &scatter, &all);
		VecScatterBegin(scatter, current_, all, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter, current_, all, INSERT_VALUES, SCATTER_FORWARD);
		const PetscScalar* values = nullptr;
		VecGetArrayRead(all, &values);
		std::vector<double> result(static_cast<std::size_t>(labels_.size())*system_.fields.size());
		for (std::size_t i = 0; i < result.size(); ++i) result[i] = PetscRealPart(values[i]);
		VecRestoreArrayRead(all, &values);
		VecScatterDestroy(&scatter);
		VecDestroy(&all);
		return result;
	}

	std::vector<double> GatherRequiredState() const
	{
		ScatterState();
		const PetscScalar* values = nullptr;
		VecGetArrayRead(ghost_state_, &values);
		std::vector<double> result(ghost_nodes_.size()*system_.fields.size());
		for (std::size_t i = 0; i < result.size(); ++i) result[i] = PetscRealPart(values[i]);
		VecRestoreArrayRead(ghost_state_, &values);
		return result;
	}

	const CompiledLinearSystem& System() const { return system_; }
	const std::vector<std::int32_t>& RequiredNodes() const { return ghost_nodes_; }
	void WriteState(const std::filesystem::path& path) const
	{
		PetscViewer viewer = nullptr;
		PetscViewerBinaryOpen(communicator_, path.string().c_str(), FILE_MODE_WRITE, &viewer);
		VecView(current_, viewer);
		PetscViewerDestroy(&viewer);
	}

	void ReadState(const std::filesystem::path& path)
	{
		PetscViewer viewer = nullptr;
		PetscViewerBinaryOpen(communicator_, path.string().c_str(), FILE_MODE_READ, &viewer);
		VecLoad(current_, viewer);
		PetscViewerDestroy(&viewer);
		steps_ = 1;
	}

	std::map<std::string, double> TotalMass(const std::vector<double>& state) const
	{
		if (state.size() != labels_.size()*system_.fields.size())
			throw std::runtime_error("VCA transport state size is invalid");
		constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
			0.6699905217924281, 0.9305681557970262}};
		constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
			0.6521451548625461, 0.3478548451374539}};
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& element : assembler_.elements()) {
			if (element.owner != rank_) continue;
			for (std::size_t qz = 0; qz < 4; ++qz)
				for (std::size_t qy = 0; qy < 4; ++qy)
					for (std::size_t qx = 0; qx < 4; ++qx) {
						const auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz]);
						const double measure = weights[qx]*weights[qy]*weights[qz]*basis.determinant;
						for (std::size_t field = 0; field < system_.fields.size(); ++field)
							for (std::size_t a = 0; a < element.connectivity.size(); ++a)
								local[field] += measure*basis.value[a]*state[
									static_cast<std::size_t>(element.connectivity[a])*system_.fields.size()+field];
					}
		}
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_DOUBLE,
			MPI_SUM, communicator_);
		std::map<std::string, double> result;
		for (std::size_t field = 0; field < system_.fields.size(); ++field)
			result.emplace(system_.fields[field], global[field]);
		return result;
	}

	std::map<std::string, double> TotalMass() const
	{
		const auto state = GatherRequiredState();
		constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
			0.6699905217924281, 0.9305681557970262}};
		constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
			0.6521451548625461, 0.3478548451374539}};
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& element : assembler_.elements()) {
			if (element.owner != rank_) continue;
			for (std::size_t qz = 0; qz < 4; ++qz)
				for (std::size_t qy = 0; qy < 4; ++qy)
					for (std::size_t qx = 0; qx < 4; ++qx) {
						const auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz]);
						const double measure = weights[qx]*weights[qy]*weights[qz]*basis.determinant;
						for (std::size_t field = 0; field < system_.fields.size(); ++field)
							for (std::size_t a = 0; a < element.connectivity.size(); ++a)
								local[field] += measure*basis.value[a]*state[
									ghost_position_.at(element.connectivity[a])*system_.fields.size()+field];
					}
		}
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_DOUBLE,
			MPI_SUM, communicator_);
		std::map<std::string, double> result;
		for (std::size_t field = 0; field < system_.fields.size(); ++field)
			result.emplace(system_.fields[field], global[field]);
		return result;
	}

	std::map<std::string, double> SourceIntegrals() const
	{
		constexpr std::array<double, 4> points{{0.06943184420297371, 0.33000947820757187,
			0.6699905217924281, 0.9305681557970262}};
		constexpr std::array<double, 4> weights{{0.3478548451374539, 0.6521451548625461,
			0.6521451548625461, 0.3478548451374539}};
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& term : system_.terms)
			if (term.kind == TermKind::VolumeSource)
				for (const auto& element : assembler_.elements()) {
					if (element.owner != rank_) continue;
					for (std::size_t qz = 0; qz < 4; ++qz)
						for (std::size_t qy = 0; qy < 4; ++qy)
							for (std::size_t qx = 0; qx < 4; ++qx) {
								const auto basis = EvaluateBasis(element, points[qx], points[qy], points[qz]);
								local[term.equation] += term.coefficient*weights[qx]*weights[qy]
									*weights[qz]*basis.determinant;
							}
				}
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_DOUBLE,
			MPI_SUM, communicator_);
		std::map<std::string, double> result;
		for (std::size_t field = 0; field < system_.fields.size(); ++field)
			result.emplace(system_.fields[field], global[field]);
		return result;
	}

private:
	void BuildGhostScatter()
	{
		for (const auto& element : assembler_.elements())
			ghost_nodes_.insert(ghost_nodes_.end(), element.connectivity.begin(), element.connectivity.end());
		std::sort(ghost_nodes_.begin(), ghost_nodes_.end());
		ghost_nodes_.erase(std::unique(ghost_nodes_.begin(), ghost_nodes_.end()), ghost_nodes_.end());
		std::vector<PetscInt> rows;
		rows.reserve(ghost_nodes_.size()*system_.fields.size());
		for (std::size_t i = 0; i < ghost_nodes_.size(); ++i) {
			ghost_position_.emplace(ghost_nodes_[i], i);
			for (std::size_t field = 0; field < system_.fields.size(); ++field)
				rows.push_back(static_cast<PetscInt>(ghost_nodes_[i]*system_.fields.size()+field));
		}
		ISCreateGeneral(communicator_, static_cast<PetscInt>(rows.size()), rows.data(),
			PETSC_COPY_VALUES, &source_rows_);
		VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), &ghost_state_);
		ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), 0, 1,
			&destination_rows_);
		VecScatterCreate(current_, source_rows_, ghost_state_, destination_rows_, &scatter_);
	}

	void ScatterState() const
	{
		VecScatterBegin(scatter_, current_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(scatter_, current_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD);
	}

	MPI_Comm communicator_;
	SimulationConfiguration configuration_;
	CompiledLinearSystem system_;
	OwnedRowAssembler assembler_;
	std::vector<int> labels_;
	std::vector<PetscInt> boundary_rows_;
	std::vector<std::int32_t> ghost_nodes_;
	std::unordered_map<std::int32_t, std::size_t> ghost_position_;
	Mat left_ = nullptr, previous_ = nullptr;
	Vec forcing_ = nullptr, current_ = nullptr, next_ = nullptr, rhs_ = nullptr;
	IS source_rows_ = nullptr, destination_rows_ = nullptr;
	Vec ghost_state_ = nullptr;
	VecScatter scatter_ = nullptr;
	KSP solver_ = nullptr;
	int steps_ = 0;
	int rank_ = 0;
};

} // namespace iga

#endif

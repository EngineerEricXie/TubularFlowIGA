#ifndef IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP
#define IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP

#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "OwnedRowAssembler.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace iga {

enum class TransportStepPhase { Committed, TrialOpen, TrialSolved, CommitPrepared };

template <class FieldValueAt>
inline void IntegrateTransportFieldMass(const Element& element, std::size_t fields,
	FieldValueAt&& value_at, const VolumeQuadratureRule& quadrature,
	std::vector<double>& result)
{
	if (result.size() != fields)
		throw std::invalid_argument("transport field-mass result has the wrong field count");
	ValidateVolumeQuadratureRule(element, quadrature);
	for (const auto& point : quadrature.Points()) {
		const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1],
			point.parametric[2]);
		const auto measure = point.weight*basis.raw_determinant;
		for (std::size_t field = 0; field < fields; ++field)
			for (std::size_t a = 0; a < element.connectivity.size(); ++a)
				result[field] += measure*basis.value[a]*value_at(element.connectivity[a], field);
	}
}

inline double IntegrateTransportPhysicalVolume(const Element& element,
	const VolumeQuadratureRule& quadrature)
{
	ValidateVolumeQuadratureRule(element, quadrature);
	double result = 0.0;
	for (const auto& point : quadrature.Points()) {
		const auto basis = EvaluateBasis(element, point.parametric[0], point.parametric[1],
			point.parametric[2]);
		result += point.weight*basis.raw_determinant;
	}
	return result;
}

class TransientTransportRuntime {
public:
	TransientTransportRuntime(Database& database, MPI_Comm communicator,
		const SimulationConfiguration& configuration,
		CompiledLinearSystem system, const std::vector<int>& labels,
		std::map<std::uint64_t, VolumeQuadratureRule> volume_rules = {})
		: communicator_(communicator), configuration_(configuration),
			system_(std::move(system)),
			coupling_patterns_(BuildTransportCouplingPatterns(system_, configuration_)),
			assembler_(database, communicator, system_.fields.size()),
			element_matrices_(coupling_patterns_),
			labels_(labels)
	{
		if (system_.velocity_source != "prescribed")
			throw std::runtime_error("in-process VCA transport requires velocity_source prescribed");
		if (labels_.size() != database.header().nodes)
			throw std::runtime_error("transport boundary labels do not match database nodes");
		BindVolumeRules(std::move(volume_rules));
		const auto boundaries = ResolveScalarBoundaries(configuration_, system_, labels_);
		left_ = assembler_.CreateMatrix(coupling_patterns_.left);
		previous_ = assembler_.CreateMatrix(coupling_patterns_.previous);
		forcing_ = assembler_.CreateVector();
		current_ = assembler_.CreateVector();
		committed_ = assembler_.CreateVector();
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
		VecCopy(current_, committed_);
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
		VecDestroy(&rhs_); VecDestroy(&next_); VecDestroy(&committed_);
		VecDestroy(&current_); VecDestroy(&forcing_);
		MatDestroy(&previous_); MatDestroy(&left_);
	}

	void Advance(const SimulationConfiguration& step_configuration,
		const std::vector<std::int32_t>& velocity_nodes,
		const std::vector<std::array<double, 3>>& velocity)
	{
		BeginStep();
		try {
			SolveTrial(step_configuration, velocity_nodes, velocity);
			PrepareCommitStep();
			FinalizeCommitStep();
		} catch (...) {
			AbortStep();
			throw;
		}
	}

	void BeginStep()
	{
		RequirePhase(TransportStepPhase::Committed, "BeginStep");
		VecCopy(current_, committed_);
		committed_steps_ = steps_;
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::TrialOpen;
	}

	void SolveTrial(const SimulationConfiguration& step_configuration,
		const std::vector<std::int32_t>& velocity_nodes,
		const std::vector<std::array<double, 3>>& velocity)
	{
		RequirePhase(TransportStepPhase::TrialOpen, "SolveTrial");
		if (velocity_nodes != ghost_nodes_ || velocity.size() != ghost_nodes_.size())
			throw std::runtime_error("VCA transport velocity nodes do not match required transport nodes");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		try {
		const auto boundaries = ResolveScalarBoundaries(step_configuration, system_, labels_);
		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node)
			for (std::size_t field = 0; field < system_.fields.size(); ++field)
				if (boundaries.constrained[static_cast<std::size_t>(node)*system_.fields.size()+field])
					boundary_rows.push_back(static_cast<PetscInt>(node*system_.fields.size()+field));
		MatZeroEntries(left_);
		MatZeroEntries(previous_);
		VecSet(forcing_, 0.0);
		for (const auto& element : assembler_.elements()) {
			BodyFittedSurface4x4QuadratureProvider surface_quadrature(element);
			BuildGenericTransportElementWithVelocity(element,
				[&velocity, this](std::int32_t node) -> const std::array<double, 3>& {
					return velocity.at(ghost_position_.at(node));
				}, system_, step_configuration, element_matrices_, VolumeRule(element),
				surface_quadrature.Rule());
			assembler_.AddElementMatrix(left_, element, element_matrices_.left);
			assembler_.AddElementMatrix(previous_, element, element_matrices_.previous);
			assembler_.AddElementVector(forcing_, element, element_matrices_.source);
		}
		OwnedRowAssembler::Assemble(left_);
		OwnedRowAssembler::Assemble(previous_);
		OwnedRowAssembler::Assemble(forcing_);
		MatZeroRows(left_, static_cast<PetscInt>(boundary_rows.size()),
			boundary_rows.data(), 1.0, nullptr, nullptr);
		MatZeroRows(previous_, static_cast<PetscInt>(boundary_rows.size()),
			boundary_rows.data(), 0.0, nullptr, nullptr);
		MatMult(previous_, current_, rhs_);
		VecAXPY(rhs_, 1.0, forcing_);
		std::vector<PetscScalar> boundary_values;
		boundary_values.reserve(boundary_rows.size());
		for (const auto row : boundary_rows)
			boundary_values.push_back(boundaries.value[static_cast<std::size_t>(row)]);
		VecSetValues(rhs_, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(),
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
		steps_ = committed_steps_+1;
		trial_solve_succeeded_ = true;
		phase_ = TransportStepPhase::TrialSolved;
		} catch (...) {
			phase_ = TransportStepPhase::TrialSolved;
			throw;
		}
	}

	void RollbackTrial()
	{
		RequirePhase(TransportStepPhase::TrialSolved, "RollbackTrial");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::TrialOpen;
	}

	void AbortStep()
	{
		if (phase_ == TransportStepPhase::Committed) return;
		if (phase_ != TransportStepPhase::TrialOpen
			&& phase_ != TransportStepPhase::TrialSolved
			&& phase_ != TransportStepPhase::CommitPrepared)
			throw std::runtime_error("illegal transport runtime transition: AbortStep");
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::Committed;
	}

	void PrepareCommitStep()
	{
		RequirePhase(TransportStepPhase::TrialSolved, "PrepareCommitStep");
		if (!trial_solve_succeeded_)
			throw std::runtime_error(
				"transport PrepareCommitStep requires a successful trial solve");
		phase_ = TransportStepPhase::CommitPrepared;
	}

	void FinalizeCommitStep() noexcept
	{
		if (phase_ != TransportStepPhase::CommitPrepared) std::terminate();
		phase_ = TransportStepPhase::Committed;
	}

	void CommitStep()
	{
		PrepareCommitStep();
		FinalizeCommitStep();
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
		RequirePhase(TransportStepPhase::Committed, "ReadState");
		PetscViewer viewer = nullptr;
		PetscViewerBinaryOpen(communicator_, path.string().c_str(), FILE_MODE_READ, &viewer);
		VecLoad(current_, viewer);
		PetscViewerDestroy(&viewer);
		steps_ = 1;
		VecCopy(current_, committed_);
		committed_steps_ = steps_;
	}

	TransportStepPhase Phase() const noexcept { return phase_; }
	int Steps() const noexcept { return steps_; }

	std::map<std::string, double> TotalMass(const std::vector<double>& state) const
	{
		if (state.size() != labels_.size()*system_.fields.size())
			throw std::runtime_error("VCA transport state size is invalid");
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& element : assembler_.elements()) {
			if (!assembler_.OwnsElementByMinimumNode(element)) continue;
			IntegrateTransportFieldMass(element, system_.fields.size(),
				[&state, this](std::int32_t node, std::size_t field) {
					return state[static_cast<std::size_t>(node)*system_.fields.size()+field];
				}, VolumeRule(element), local);
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
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& element : assembler_.elements()) {
			if (!assembler_.OwnsElementByMinimumNode(element)) continue;
			IntegrateTransportFieldMass(element, system_.fields.size(),
				[&state, this](std::int32_t node, std::size_t field) {
					return state[ghost_position_.at(node)*system_.fields.size()+field];
				}, VolumeRule(element), local);
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
		std::vector<double> local(system_.fields.size(), 0.0), global(system_.fields.size(), 0.0);
		for (const auto& term : system_.terms)
			if (term.kind == TermKind::VolumeSource)
				for (const auto& element : assembler_.elements()) {
					if (!assembler_.OwnsElementByMinimumNode(element)) continue;
					local[term.equation] += term.coefficient
						*IntegrateTransportPhysicalVolume(element, VolumeRule(element));
				}
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_DOUBLE,
			MPI_SUM, communicator_);
		std::map<std::string, double> result;
		for (std::size_t field = 0; field < system_.fields.size(); ++field)
			result.emplace(system_.fields[field], global[field]);
		return result;
	}

private:
	void BindVolumeRules(std::map<std::uint64_t, VolumeQuadratureRule> volume_rules)
	{
		if (volume_rules.empty())
			for (const auto& element : assembler_.elements()) {
				FullCell4x4x4VolumeQuadratureProvider provider(element);
				volume_rules.emplace(element.id, provider.Rule());
			}
		if (volume_rules.size() != assembler_.elements().size())
			throw std::runtime_error("transport volume quadrature catalog does not cover each local element exactly once");
		for (const auto& element : assembler_.elements()) {
			const auto found = volume_rules.find(element.id);
			if (found == volume_rules.end())
				throw std::runtime_error("transport volume quadrature catalog is missing an element id");
			ValidateVolumeQuadratureRule(element, found->second);
		}
		volume_rules_ = std::move(volume_rules);
	}

	const VolumeQuadratureRule& VolumeRule(const Element& element) const
	{
		return volume_rules_.at(element.id);
	}

	void RequirePhase(TransportStepPhase required, const char* operation) const
	{
		if (phase_ != required)
			throw std::runtime_error(std::string("transport ")+operation
				+" called in an invalid lifecycle phase");
	}

	void RestoreCommitted()
	{
		VecCopy(committed_, current_);
		steps_ = committed_steps_;
	}

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
	TransportCouplingPatterns coupling_patterns_;
	OwnedRowAssembler assembler_;
	GenericTransportMatrices element_matrices_;
	std::map<std::uint64_t, VolumeQuadratureRule> volume_rules_;
	std::vector<int> labels_;
	std::vector<std::int32_t> ghost_nodes_;
	std::unordered_map<std::int32_t, std::size_t> ghost_position_;
	Mat left_ = nullptr, previous_ = nullptr;
	Vec forcing_ = nullptr, current_ = nullptr, committed_ = nullptr;
	Vec next_ = nullptr, rhs_ = nullptr;
	IS source_rows_ = nullptr, destination_rows_ = nullptr;
	Vec ghost_state_ = nullptr;
	VecScatter scatter_ = nullptr;
	KSP solver_ = nullptr;
	int steps_ = 0;
	int committed_steps_ = 0;
	bool trial_solve_succeeded_ = false;
	TransportStepPhase phase_ = TransportStepPhase::Committed;
};

} // namespace iga

#endif

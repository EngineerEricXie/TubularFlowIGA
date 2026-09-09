#ifndef IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP
#define IGA_TRANSIENT_TRANSPORT_RUNTIME_HPP

#include "RuntimeCleanup.hpp"
#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "OwnedRowAssembler.hpp"
#include "CollectivePetscOptions.hpp"
#include "PetscSolverOptions.hpp"
#include "PetscReadArray.hpp"
#include "OwnedCheckpointVector.hpp"

#include <petscksp.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <map>
#include "PetscGather.hpp"
#include "PetscCheckpointRead.hpp"
#include "PetscCheckpointWrite.hpp"
#include <limits>
#include <iomanip>
#include <sstream>
#include <set>
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

struct TransportAcceptedCheckpointState {
	std::string configuration_identity_sha256;
	int accepted_steps = 0;
	OwnedCheckpointVector field;
};

class TransientTransportRuntime {
public:
	// Borrows communicator for all solves, reductions and checkpoint viewers.
	// Destroy this runtime before the owner frees it and before PetscFinalize.
	TransientTransportRuntime(Database& database, MPI_Comm communicator,
		const SimulationConfiguration& configuration,
		const CompiledLinearSystem& system, const std::vector<int>& labels,
		const std::map<std::uint64_t, VolumeQuadratureRule>& volume_rules = {},
		const std::set<std::string>& application_options = {},
		const std::string& checkpoint_identity_sha256 = {},
		const std::string& solver_options_prefix = {})
		: communicator_(communicator),
			configuration_(PrepareRuntimeConstructionInput<SimulationConfiguration>(communicator,
				"transport configuration preparation", configuration)),
			system_(PrepareRuntimeConstructionInput<CompiledLinearSystem>(communicator,
			"transport compiled input preparation", system)),
			assembler_(database, communicator, system_.fields.size())
	{
		try {
			ResolvedScalarBoundaries boundaries;
			std::vector<double> initial;
			RuntimeConstructionStage(communicator_, "transport runtime input", [&] {
				checkpoint_identity_sha256_ = checkpoint_identity_sha256;
				if (!checkpoint_identity_sha256_.empty()) ValidateCheckpointConfigurationIdentity(checkpoint_identity_sha256_);
				labels_ = labels;
				if (system_.velocity_source != "prescribed")
					throw std::runtime_error("in-process VCA transport requires velocity_source prescribed");
				if (labels_.size() != database.header().nodes)
					throw std::runtime_error("transport boundary labels do not match database nodes");
				coupling_patterns_ = BuildTransportCouplingPatterns(system_, configuration_);
				element_matrices_ = GenericTransportMatrices(coupling_patterns_);
				BindVolumeRules(volume_rules);
				boundaries = ResolveScalarBoundaries(configuration_, system_, labels_);
				initial = InitialScalarValues(configuration_, system_);
			});
			left_ = assembler_.CreateMatrix(coupling_patterns_.left);
			previous_ = assembler_.CreateMatrix(coupling_patterns_.previous);
			forcing_ = assembler_.CreateVector();
			current_ = assembler_.CreateVector();
			committed_ = assembler_.CreateVector();
			next_ = assembler_.CreateVector();
			rhs_ = assembler_.CreateVector();
			RequireCollectivePetscSuccess(communicator_, "transport left pattern policy",
				MatSetOption(left_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
			RequireCollectivePetscSuccess(communicator_, "transport previous pattern policy",
				MatSetOption(previous_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
			RequireCollectivePetscSuccess(communicator_, "transport previous zero policy",
				MatSetOption(previous_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
			RuntimeConstructionStage(communicator_, "transport initial values", [&] {
				// VecCreateMPI fixes the CPU MPI vector type; array access and
				// restoration here are local, including for an empty rank.
				PetscScalar* values = nullptr;
				if (VecGetArray(current_, &values)) throw std::runtime_error("cannot access initial transport array");
				try {
					for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node)
						for (std::size_t field = 0; field < system_.fields.size(); ++field) {
							const auto local = static_cast<std::size_t>(node - assembler_.node_begin()) * system_.fields.size() + field;
							const auto global = static_cast<std::size_t>(node) * system_.fields.size() + field;
							values[local] = boundaries.constrained[global] ? boundaries.value[global] : initial[field];
						}
					ProbeRuntimeConstruction("transport initial array active");
				} catch (...) { VecRestoreArray(current_, &values); throw; }
				if (VecRestoreArray(current_, &values)) throw std::runtime_error("cannot restore initial transport array");
			});
			RequireCollectivePetscSuccess(communicator_, "transport initial state copy", VecCopy(current_, committed_));
			BuildGhostScatter();
			RequireCollectivePetscOptions(communicator_, nullptr, application_options);
			RequireCollectivePetscSuccess(communicator_, "transport solver creation", KSPCreate(communicator_, &solver_));
			ObserveConstructedObject(communicator_, "transport solver created", reinterpret_cast<PetscObject>(solver_));
			const PetscOptionEntries solver_defaults;
			solver_options_ = AllocateCollectiveRuntime<PetscSolverOptions>(communicator_, communicator_, solver_options_prefix,
				nullptr, solver_defaults, application_options);
			solver_options_->Attach(solver_);
			RequireCollectivePetscSuccess(communicator_, "transport solver type", KSPSetType(solver_, KSPGMRES));
			RequireCollectivePetscSuccess(communicator_, "transport solver restart", KSPGMRESSetRestart(solver_, 50));
			RequireCollectivePetscSuccess(communicator_, "transport solver tolerances",
				KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000));
			PC preconditioner = nullptr;
			RequireCollectivePetscSuccess(communicator_, "transport preconditioner lookup", KSPGetPC(solver_, &preconditioner));
			RequireCollectivePetscSuccess(communicator_, "transport preconditioner type", PCSetType(preconditioner, PCBJACOBI));
			solver_options_->Call("transport solver options", [&] { return KSPSetFromOptions(solver_); });
			solver_options_->RecordUsed();
			RuntimeConstructionStage(communicator_, "transport runtime ready", [] {});
		} catch (...) { DestroyPetsc(); throw; }
	}

	~TransientTransportRuntime() { DestroyPetsc(); }

	// Terminal collective operation. Call before reporting success and before
	// releasing the borrowed communicator. Only Close/destruction may follow.
	void Close()
	{
		DestroyPetsc().Check(communicator_, "transport runtime cleanup");
	}

	PetscKspConfiguration SolverConfiguration() const
	{
		if (cleanup_started_) throw std::logic_error("solver configuration requested after runtime close");
		return CaptureKspConfiguration(solver_);
	}

	TransientTransportRuntime(const TransientTransportRuntime&) = delete;
	TransientTransportRuntime& operator=(const TransientTransportRuntime&) = delete;

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
		CollectiveLocalStage(communicator_, "transport begin step preparation", [&] {
			RequirePhase(TransportStepPhase::Committed, "BeginStep");
			if (steps_ == std::numeric_limits<int>::max()) throw std::runtime_error("transport step counter overflows");
		});
		RequireCollectiveSameInt(communicator_, "transport step agreement", steps_);
		RequireCollectivePetscSuccess(communicator_, "transport save state", VecCopy(current_, committed_));
		committed_steps_ = steps_;
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::TrialOpen;
	}

	void SolveTrial(const SimulationConfiguration& step_configuration,
		const std::vector<std::int32_t>& velocity_nodes,
		const std::vector<std::array<double, 3>>& velocity)
	{
		CollectiveLocalStage(communicator_, "transport trial input", [&] {
			RequirePhase(TransportStepPhase::TrialOpen, "SolveTrial");
			if (velocity_nodes != ghost_nodes_ || velocity.size() != ghost_nodes_.size())
				throw std::runtime_error("VCA transport velocity nodes do not match required transport nodes");
		});
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		try {
		ResolvedScalarBoundaries boundaries;
		std::vector<PetscInt> boundary_rows;
		CollectiveLocalStage(communicator_, "transport trial boundaries", [&] {
			boundaries = ResolveScalarBoundaries(step_configuration, system_, labels_);
			for (std::uint64_t node = assembler_.node_begin(); node < assembler_.node_end(); ++node)
				for (std::size_t field = 0; field < system_.fields.size(); ++field)
					if (boundaries.constrained[static_cast<std::size_t>(node)*system_.fields.size()+field])
						boundary_rows.push_back(static_cast<PetscInt>(node*system_.fields.size()+field));
		});
		RequireCollectivePetscSuccess(communicator_, "transport clear left", MatZeroEntries(left_));
		RequireCollectivePetscSuccess(communicator_, "transport clear previous", MatZeroEntries(previous_));
		RequireCollectivePetscSuccess(communicator_, "transport clear forcing", VecSet(forcing_, 0.0));
		CollectiveLocalStage(communicator_, "transport element assembly", [&] {
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
		});
		OwnedRowAssembler::Assemble(left_, communicator_);
		OwnedRowAssembler::Assemble(previous_, communicator_);
		OwnedRowAssembler::Assemble(forcing_, communicator_);
		RequireCollectivePetscSuccess(communicator_, "transport left boundary rows",
			MatZeroRows(left_, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 1.0, nullptr, nullptr));
		RequireCollectivePetscSuccess(communicator_, "transport previous boundary rows",
			MatZeroRows(previous_, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), 0.0, nullptr, nullptr));
		RequireCollectivePetscSuccess(communicator_, "transport previous product", MatMult(previous_, current_, rhs_));
		RequireCollectivePetscSuccess(communicator_, "transport forcing", VecAXPY(rhs_, 1.0, forcing_));
		CollectiveLocalStage(communicator_, "transport boundary insertion", [&] {
			std::vector<PetscScalar> boundary_values;
			boundary_values.reserve(boundary_rows.size());
			for (const auto row : boundary_rows)
				boundary_values.push_back(boundaries.value[static_cast<std::size_t>(row)]);
			if (VecSetValues(rhs_, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(),
				boundary_values.data(), INSERT_VALUES))
				throw std::runtime_error("transport boundary VecSetValues failed");
		});
		OwnedRowAssembler::Assemble(rhs_, communicator_);
		// PREONLY applies the preconditioner directly and rejects a nonzero
		// initial guess. Preserve warm starts for iterative solvers only.
		PetscBool preonly = PETSC_FALSE;
		CollectiveLocalStage(communicator_, "transport solver type", [&] {
			if (PetscObjectTypeCompare(reinterpret_cast<PetscObject>(solver_), KSPPREONLY, &preonly))
				throw std::runtime_error("cannot query transport solver type");
		});
		const bool warm_start = steps_ > 0 && !preonly;
		RequireCollectiveSameInt(communicator_, "transport warm start agreement", warm_start ? 1 : 0);
		if (warm_start)
			RequireCollectivePetscSuccess(communicator_, "transport warm start", VecCopy(current_, next_));
		RequireKspFactorBackend(solver_, left_, communicator_);
		RequireCollectivePetscSuccess(communicator_, "transport solver operators", KSPSetOperators(solver_, left_, left_));
		RequireCollectivePetscSuccess(communicator_, "transport solver initial guess", KSPSetInitialGuessNonzero(solver_, warm_start ? PETSC_TRUE : PETSC_FALSE));
		solver_options_->Call("transport solver setup", [&] { return KSPSetUp(solver_); });
		solver_options_->Call("transport solve", [&] { return KSPSolve(solver_, rhs_, next_); });
		solver_options_->RecordUsed();
		CollectiveLocalStage(communicator_, "transport convergence", [&] {
			KSPConvergedReason reason;
			if (KSPGetConvergedReason(solver_, &reason))
				throw std::runtime_error("cannot query transport convergence reason");
			if (reason <= 0) throw std::runtime_error("VCA transport linear solve did not converge");
		});
		RequireCollectivePetscSuccess(communicator_, "transport publish state", VecSwap(current_, next_));
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
		CollectiveLocalStage(communicator_, "transport rollback preparation", [&] {
			RequirePhase(TransportStepPhase::TrialSolved, "RollbackTrial");
		});
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::TrialOpen;
	}

	void AbortStep()
	{
		RequireCollectiveSameInt(communicator_, "transport abort phase agreement", static_cast<int>(phase_));
		CollectiveLocalStage(communicator_, "transport abort preparation", [&] {
			if (phase_ != TransportStepPhase::Committed && phase_ != TransportStepPhase::TrialOpen
				&& phase_ != TransportStepPhase::TrialSolved && phase_ != TransportStepPhase::CommitPrepared)
				throw std::runtime_error("illegal transport runtime transition: AbortStep");
		});
		if (phase_ == TransportStepPhase::Committed) return;
		RestoreCommitted();
		trial_solve_succeeded_ = false;
		phase_ = TransportStepPhase::Committed;
	}

	void PrepareCommitStep()
	{
		CollectiveLocalStage(communicator_, "transport prepare commit", [&] {
			RequirePhase(TransportStepPhase::TrialSolved, "PrepareCommitStep");
			if (!trial_solve_succeeded_)
				throw std::runtime_error("transport PrepareCommitStep requires a successful trial solve");
		});
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
		std::uint64_t expected = 0;
		std::string signature;
		CollectiveLocalStage(communicator_, "transport global state preparation", [&] {
			if (system_.fields.empty() || labels_.size() > std::numeric_limits<std::uint64_t>::max()/system_.fields.size())
				throw std::runtime_error("invalid transport global state dimensions");
			expected = static_cast<std::uint64_t>(labels_.size())*system_.fields.size();
			if (expected != static_cast<std::uint64_t>(assembler_.global_rows()))
				throw std::runtime_error("transport global state layout differs from its assembled fields");
			std::set<std::string> names;
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << expected << ' ' << system_.fields.size();
			for (const auto& field : system_.fields) {
				if (field.empty() || !names.insert(field).second)
					throw std::runtime_error("global transport fields must be nonempty and unique");
				text << ' ' << std::quoted(field);
			}
			signature = text.str();
		});
		RequireCollectiveSameText(communicator_, "transport global state agreement", signature);
		return GatherAllPetscReal(current_, communicator_, expected);
	}

	std::vector<double> GatherRequiredState() const
	{
		ScatterState();
		std::vector<double> result;
		CollectiveLocalStage(communicator_, "transport required state", [&] {
			PetscReadArray view;
			view.Acquire(ghost_state_);
			result.resize(ghost_nodes_.size()*system_.fields.size());
			for (std::size_t i = 0; i < result.size(); ++i) result[i] = PetscRealPart(view.Data()[i]);
			view.Restore();
		});
		return result;
	}

	const CompiledLinearSystem& System() const { return system_; }
#ifdef IGA_TRANSPORT_RUNTIME_TESTING
	VolumeQuadratureRule& VolumeRuleForTesting(std::uint64_t element) { return volume_rules_.at(element); }
#endif
	const std::vector<std::int32_t>& RequiredNodes() const { return ghost_nodes_; }
	void WriteState(const std::filesystem::path& path) const
	{
		CollectiveLocalStage(communicator_, "transport checkpoint write preparation", [&] {
			RequirePhase(TransportStepPhase::Committed, "WriteState");
		});
		WritePetscCheckpointVector(current_, communicator_, path);
	}

	void ReadState(const std::filesystem::path& path)
	{
		CollectiveLocalStage(communicator_, "transport checkpoint preparation", [&] {
			RequirePhase(TransportStepPhase::Committed, "ReadState");
		});
		ReadPetscCheckpointVector(next_, communicator_, path);
		RequireCollectivePetscSuccess(communicator_, "transport checkpoint snapshot", VecCopy(next_, committed_));
		RequireCollectivePetscSuccess(communicator_, "transport checkpoint publish", VecSwap(current_, next_));
		steps_ = 1;
		committed_steps_ = steps_;
	}

	TransportAcceptedCheckpointState CreateCheckpointRestoreCandidate() const
	{
		TransportAcceptedCheckpointState value;
		CollectiveLocalStage(communicator_, "transport checkpoint candidate layout", [&] {
			RequirePhase(TransportStepPhase::Committed, "CreateCheckpointRestoreCandidate");
			if (cleanup_started_ || steps_ != 0) throw std::runtime_error("transport decoder requires a fresh open target");
			ValidateCheckpointConfigurationIdentity(checkpoint_identity_sha256_);
			value.configuration_identity_sha256 = checkpoint_identity_sha256_;
			value.field = CaptureOwnedCheckpointVector(current_);
		});
		return value;
	}

	TransportAcceptedCheckpointState CaptureCheckpointState() const
	{
		TransportAcceptedCheckpointState value;
		CollectiveLocalStage(communicator_, "transport accepted checkpoint capture", [&] {
			RequirePhase(TransportStepPhase::Committed, "CaptureCheckpointState");
			if (cleanup_started_) throw std::runtime_error("checkpoint runtime is closed");
			ValidateCheckpointConfigurationIdentity(checkpoint_identity_sha256_);
			if (steps_ <= 0) throw std::runtime_error("transport checkpoint requires an accepted step");
			value.configuration_identity_sha256 = checkpoint_identity_sha256_;
			value.accepted_steps = steps_; value.field = CaptureOwnedCheckpointVector(current_);
		});
		RequireCollectiveSameText(communicator_, "transport checkpoint identity agreement", value.configuration_identity_sha256);
		RequireCollectiveSameInt(communicator_, "transport checkpoint clock agreement", value.accepted_steps);
		return value;
	}

	// Restore only an unpublished fresh runtime. All local validation and
	// staging precede publication. A PETSc failure discards this whole candidate.
	void RestoreCheckpointState(const TransportAcceptedCheckpointState& value)
	{
		CollectiveLocalStage(communicator_, "transport accepted checkpoint validation", [&] {
			RequirePhase(TransportStepPhase::Committed, "RestoreCheckpointState");
			if (cleanup_started_) throw std::runtime_error("checkpoint runtime is closed");
			ValidateCheckpointConfigurationIdentity(checkpoint_identity_sha256_);
			if (steps_ != 0 || value.accepted_steps <= 0
				|| value.configuration_identity_sha256 != checkpoint_identity_sha256_)
				throw std::runtime_error("transport checkpoint requires a compatible fresh target and accepted count");
			ValidateOwnedCheckpointVector(current_, value.field);
		});
		RequireCollectiveSameText(communicator_, "transport checkpoint identity agreement", value.configuration_identity_sha256);
		RequireCollectiveSameInt(communicator_, "transport checkpoint clock agreement", value.accepted_steps);
		CollectiveLocalStage(communicator_, "transport accepted checkpoint staging", [&] { StageOwnedCheckpointVector(next_, value.field); });
		RequireCollectivePetscSuccess(communicator_, "transport accepted checkpoint snapshot", VecCopy(next_, committed_));
		RequireCollectivePetscSuccess(communicator_, "transport accepted checkpoint publication", VecSwap(current_, next_));
		steps_ = value.accepted_steps; committed_steps_ = steps_;
	}

	TransportStepPhase Phase() const noexcept { return phase_; }
	int Steps() const noexcept { return steps_; }

	std::map<std::string, double> TotalMass(const std::vector<double>& state) const
	{
		return IntegrateFields("transport mass preparation", "transport mass integration",
			[&](std::ostringstream&) {
				if (state.size() != labels_.size()*system_.fields.size())
					throw std::runtime_error("VCA transport state size is invalid");
			}, [&](std::vector<double>& local) {
				for (const auto& element : assembler_.elements()) {
					if (!assembler_.OwnsElementByMinimumNode(element)) continue;
					IntegrateTransportFieldMass(element, system_.fields.size(),
						[&state, this](std::int32_t node, std::size_t field) {
							return state.at(static_cast<std::size_t>(node)*system_.fields.size()+field);
						}, VolumeRule(element), local);
				}
			});
	}

	std::map<std::string, double> TotalMass() const
	{
		const auto state = GatherRequiredState();
		return IntegrateFields("transport mass preparation", "transport mass integration",
			[](std::ostringstream&) {}, [&](std::vector<double>& local) {
				for (const auto& element : assembler_.elements()) {
					if (!assembler_.OwnsElementByMinimumNode(element)) continue;
					IntegrateTransportFieldMass(element, system_.fields.size(),
						[&state, this](std::int32_t node, std::size_t field) {
							return state.at(ghost_position_.at(node)*system_.fields.size()+field);
						}, VolumeRule(element), local);
				}
			});
	}

	std::map<std::string, double> SourceIntegrals() const
	{
		return IntegrateFields("transport source preparation", "transport source integration",
			[&](std::ostringstream& text) {
				for (const auto& term : system_.terms) if (term.kind == TermKind::VolumeSource) {
					if (term.equation >= system_.fields.size() || !std::isfinite(term.coefficient))
						throw std::runtime_error("invalid transport volume source");
					text << ' ' << term.equation << ' ' << std::hexfloat << term.coefficient;
				}
			}, [&](std::vector<double>& local) {
				for (const auto& term : system_.terms)
					if (term.kind == TermKind::VolumeSource)
						for (const auto& element : assembler_.elements()) {
							if (!assembler_.OwnsElementByMinimumNode(element)) continue;
							local.at(term.equation) += term.coefficient
								*IntegrateTransportPhysicalVolume(element, VolumeRule(element));
						}
			});
	}

private:
	bool cleanup_started_ = false;
	RuntimeCleanupResult cleanup_result_;

	RuntimeCleanupResult DestroyPetsc() noexcept
	{
		if (cleanup_started_) return cleanup_result_;
		cleanup_started_ = true;
		cleanup_result_.Observe("transport solver", KSPDestroy(&solver_));
		cleanup_result_.Observe("transport scatter", VecScatterDestroy(&scatter_));
		cleanup_result_.Observe("transport destination_rows", ISDestroy(&destination_rows_));
		cleanup_result_.Observe("transport ghost_state", VecDestroy(&ghost_state_));
		cleanup_result_.Observe("transport source_rows", ISDestroy(&source_rows_));
		cleanup_result_.Observe("transport rhs", VecDestroy(&rhs_));
		cleanup_result_.Observe("transport next", VecDestroy(&next_));
		cleanup_result_.Observe("transport committed", VecDestroy(&committed_));
		cleanup_result_.Observe("transport current", VecDestroy(&current_));
		cleanup_result_.Observe("transport forcing", VecDestroy(&forcing_));
		cleanup_result_.Observe("transport previous", MatDestroy(&previous_));
		cleanup_result_.Observe("transport left", MatDestroy(&left_));
		return cleanup_result_;
	}

	// Callbacks are local work only. Both the reduction buffers and their
	// ordered physical meaning are agreed before integration/reduction.
	template <class Prepare, class Integrate>
	std::map<std::string, double> IntegrateFields(const char* preparation_stage,
		const char* integration_stage, Prepare&& prepare, Integrate&& integrate) const
	{
		std::vector<double> local, global;
		std::string signature;
		CollectiveLocalStage(communicator_, preparation_stage, [&] {
			if (system_.fields.empty()
				|| system_.fields.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
				throw std::runtime_error("invalid transport integral field count");
			std::set<std::string> names;
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			text << system_.fields.size();
			for (const auto& field : system_.fields) {
				if (field.empty() || !names.insert(field).second)
					throw std::runtime_error("transport integral fields must be nonempty and unique");
				text << ' ' << std::quoted(field);
			}
			prepare(text);
			signature = text.str();
			local.assign(system_.fields.size(), 0.0); global.resize(local.size());
		});
		RequireCollectiveSameText(communicator_, "transport integral agreement", signature);
		CollectiveLocalStage(communicator_, integration_stage, [&] { integrate(local); });
		MPI_Allreduce(local.data(), global.data(), static_cast<int>(global.size()), MPI_DOUBLE,
			MPI_SUM, communicator_);
		std::map<std::string, double> result;
		CollectiveLocalStage(communicator_, "transport integral result", [&] {
			for (std::size_t field = 0; field < system_.fields.size(); ++field) {
				if (!std::isfinite(global[field])) throw std::runtime_error("nonfinite transport integral");
				result.emplace(system_.fields[field], global[field]);
			}
		});
		return result;
	}

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
		RequireCollectivePetscSuccess(communicator_, "transport restore state", VecCopy(committed_, current_));
		steps_ = committed_steps_;
	}

	void BuildGhostScatter()
	{
		std::vector<PetscInt> rows;
		RuntimeConstructionStage(communicator_, "transport halo preparation", [&] {
			for (const auto& element : assembler_.elements())
				ghost_nodes_.insert(ghost_nodes_.end(), element.connectivity.begin(), element.connectivity.end());
			std::sort(ghost_nodes_.begin(), ghost_nodes_.end());
			ghost_nodes_.erase(std::unique(ghost_nodes_.begin(), ghost_nodes_.end()), ghost_nodes_.end());
			rows = assembler_.RequiredRows(ghost_nodes_);
			for (std::size_t i = 0; i < ghost_nodes_.size(); ++i)
				ghost_position_.emplace(ghost_nodes_[i], i);
		});
		RequireCollectivePetscSuccess(communicator_, "transport halo source index creation",
			ISCreateGeneral(communicator_, static_cast<PetscInt>(rows.size()), rows.data(), PETSC_COPY_VALUES, &source_rows_));
		ObserveConstructedObject(communicator_, "transport halo source index created", reinterpret_cast<PetscObject>(source_rows_));
		RequireCollectivePetscSuccess(communicator_, "transport halo state creation",
			VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), &ghost_state_));
		ObserveConstructedObject(communicator_, "transport halo state created", reinterpret_cast<PetscObject>(ghost_state_));
		RequireCollectivePetscSuccess(communicator_, "transport halo destination index creation",
			ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(rows.size()), 0, 1, &destination_rows_));
		ObserveConstructedObject(communicator_, "transport halo destination index created", reinterpret_cast<PetscObject>(destination_rows_));
		RequireCollectivePetscSuccess(communicator_, "transport halo scatter creation",
			VecScatterCreate(current_, source_rows_, ghost_state_, destination_rows_, &scatter_));
		ObserveConstructedObject(communicator_, "transport halo scatter created", reinterpret_cast<PetscObject>(scatter_));
	}

	void ScatterState() const
	{
		RequireCollectivePetscSuccess(communicator_, "transport state scatter begin", VecScatterBegin(scatter_, current_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD));
		RequireCollectivePetscSuccess(communicator_, "transport state scatter end", VecScatterEnd(scatter_, current_, ghost_state_, INSERT_VALUES, SCATTER_FORWARD));
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
	std::unique_ptr<PetscSolverOptions> solver_options_;
	KSP solver_ = nullptr;
	std::string checkpoint_identity_sha256_;
	int steps_ = 0;
	int committed_steps_ = 0;
	bool trial_solve_succeeded_ = false;
	TransportStepPhase phase_ = TransportStepPhase::Committed;
};

} // namespace iga

#endif

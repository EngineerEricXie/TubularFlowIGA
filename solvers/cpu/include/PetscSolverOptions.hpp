#ifndef IGA_PETSC_SOLVER_OPTIONS_HPP
#define IGA_PETSC_SOLVER_OPTIONS_HPP

#include "CollectivePetscOptions.hpp"
#include "Sha256.hpp"
#include <petscksp.h>

namespace iga {

inline std::string PetscDomainOptionsPrefix(const std::string& domain, const std::string& role)
{
	if (domain.empty() || role.empty() || role.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") != std::string::npos)
		throw std::invalid_argument("invalid solver domain or role");
	// Lowercase ordinary IDs remain readable. Separate namespaces avoid
	// collisions with encoded IDs, including case-insensitive PETSc lookups.
	if (domain.size() <= 64 && domain.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_") == std::string::npos)
		return "domain_"+domain+"_"+role+"_";
	Sha256 hash; hash.Append(domain.data(), domain.size());
	return "domainh_"+hash.Hex()+"_"+role+"_";
}

struct PetscKspConfiguration {
	std::string prefix, ksp, pc, factor_backend;
	PetscReal relative_tolerance = 0.0, absolute_tolerance = 0.0, divergence_tolerance = 0.0;
	PetscInt maximum_iterations = 0, last_iterations = 0;
	KSPConvergedReason last_reason = KSP_CONVERGED_ITERATING;
};

// Local, read-only description of the effective top-level solver. A factor
// backend named "auto" has not yet been selected; this is not a capability claim.
inline PetscKspConfiguration CaptureKspConfiguration(KSP solver)
{
	using petsc_options_detail::Check;
	PetscKspConfiguration result; const char* name = nullptr; PC pc = nullptr;
	Check(KSPGetOptionsPrefix(solver, &name), "KSPGetOptionsPrefix"); result.prefix = name ? name : "";
	Check(KSPGetType(solver, &name), "KSPGetType"); result.ksp = name ? name : "unset";
	Check(KSPGetPC(solver, &pc), "KSPGetPC"); Check(PCGetType(pc, &name), "PCGetType"); result.pc = name ? name : "unset";
	if (result.pc == PCLU || result.pc == PCCHOLESKY || result.pc == PCILU || result.pc == PCICC) {
		Check(PCFactorGetMatSolverType(pc, &name), "PCFactorGetMatSolverType"); result.factor_backend = name ? name : "auto";
	} else result.factor_backend = "none";
	Check(KSPGetTolerances(solver, &result.relative_tolerance, &result.absolute_tolerance, &result.divergence_tolerance, &result.maximum_iterations), "KSPGetTolerances");
	Check(KSPGetIterationNumber(solver, &result.last_iterations), "KSPGetIterationNumber");
	Check(KSPGetConvergedReason(solver, &result.last_reason), "KSPGetConvergedReason");
	return result;
}

// An immutable per-runtime snapshot. Explicit prefixed options override inherited
// unprefixed options, which override the runtime's supplied option defaults.
// Keep this owner alive until every attached PETSc object has been destroyed.
class PetscSolverOptions {
public:
	PetscSolverOptions(MPI_Comm communicator, const std::string& prefix,
		PetscOptions source = nullptr, const PetscOptionEntries& defaults = {},
		const std::set<std::string>& excluded = {}) : communicator_(communicator), source_(source)
	{
		std::string agreement;
		try {
			CollectiveLocalStage(communicator_, "solver options snapshot", [&] {
				if (prefix.size() > 128 || (!prefix.empty() && (prefix.back() != '_'
					|| prefix.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_") != std::string::npos)))
					throw std::invalid_argument("solver options prefix must contain only ASCII letters, digits or underscores and end in underscore");
				prefix_ = prefix;
				const auto inputs = CapturePetscOptionEntries(source, excluded);
				auto entries = defaults;
				for (const auto& entry : inputs) entries[entry.first] = entry.second;
				PetscOptionEntries effective;
				for (const auto& entry : entries) {
					if (entry.first.size() < 2 || entry.first.front() != '-') throw std::invalid_argument("invalid default solver option name");
					const auto name = "-"+prefix_+entry.first.substr(1);
					if (name.size() >= PETSC_MAX_OPTION_NAME) {
						if (defaults.count(entry.first)) throw std::invalid_argument("prefixed default solver option exceeds PETSc name capacity");
						// Unrelated long source keys remain available unchanged;
						// PETSc could not request this overlong prefixed alias.
						continue;
					}
					effective[name] = entry.second;
					if (inputs.count(entry.first)) origins_[name] = entry.first;
				}
				// Retain the original keys too: PETSc can query non-prefixed global
				// controls, and explicit keys must win over inherited aliases.
				for (const auto& entry : inputs) { effective[entry.first] = entry.second; origins_[entry.first] = entry.first; }
				agreement = std::to_string(prefix_.size())+":"+prefix_+SerializePetscOptionEntries(effective);
				petsc_options_detail::Check(PetscOptionsCreate(&options_), "PetscOptionsCreate solver snapshot");
				for (const auto& entry : effective)
					petsc_options_detail::Check(PetscOptionsSetValue(options_, entry.first.c_str(), entry.second ? entry.second->c_str() : nullptr), "PetscOptionsSetValue solver snapshot");
			});
			RequireCollectiveSameText(communicator_, "solver options snapshot agreement", agreement);
		} catch (...) { PetscOptionsDestroy(&options_); throw; }
	}
	~PetscSolverOptions() { PetscOptionsDestroy(&options_); }
	PetscSolverOptions(const PetscSolverOptions&) = delete;
	PetscSolverOptions& operator=(const PetscSolverOptions&) = delete;
	const std::string& Prefix() const noexcept { return prefix_; }
	PetscOptions Database() const noexcept { return options_; }

	// Attach before creating fieldsplit children and before SetFromOptions.
	// An already-created, unconfigured PC is updated as well.
	void Attach(KSP solver) const
	{
		PC pc = nullptr;
		RequireCollectivePetscSuccess(communicator_, "solver options KSP database", PetscObjectSetOptions(reinterpret_cast<PetscObject>(solver), options_));
		RequireCollectivePetscSuccess(communicator_, "solver options PC lookup", KSPGetPC(solver, &pc));
		RequireCollectivePetscSuccess(communicator_, "solver options PC database", PetscObjectSetOptions(reinterpret_cast<PetscObject>(pc), options_));
		RequireCollectivePetscSuccess(communicator_, "solver options prefix", KSPSetOptionsPrefix(solver, prefix_.c_str()));
	}

	// PETSc 3.15 propagates the prefix but not the database to some child KSPs.
	// Activate the private database for each SetFromOptions/SetUp/Solve call,
	// A returning error handler prevents the default handler from aborting peer
	// ranks on invalid configuration. Restore both stacks before error agreement.
	// Like PETSc itself, this scope belongs to the MPI initialization thread;
	// simultaneous solver calls from different host threads are unsupported.
	template<class Operation> void Call(const char* stage, Operation&& operation)
	{
		bool pushed = false, handler_pushed = false;
		try {
			CollectiveLocalStage(communicator_, "solver options activation", [&] {
				if (active_) throw std::logic_error("solver option database is already active");
				petsc_options_detail::Check(PetscOptionsPush(options_), "PetscOptionsPush solver snapshot");
				active_ = pushed = true;
				petsc_options_detail::Check(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr), "PetscPushErrorHandler solver scope");
				handler_pushed = true;
			});
		} catch (...) {
			if (handler_pushed) PetscPopErrorHandler();
			if (pushed) { PetscOptionsPop(); active_ = false; }
			throw;
		}
		PetscErrorCode code = 0;
		std::exception_ptr failure;
		try { code = operation(); } catch (...) { failure = std::current_exception(); }
		const auto handler_restored = PetscPopErrorHandler();
		const auto restored = PetscOptionsPop(); active_ = false;
		CollectiveLocalStage(communicator_, stage, [&] {
			if (failure) std::rethrow_exception(failure);
			if (code) throw std::runtime_error("PETSc returned error "+std::to_string(code));
			petsc_options_detail::Check(handler_restored, "PetscPopErrorHandler solver scope");
			petsc_options_detail::Check(restored, "PetscOptionsPop solver snapshot");
		});
	}

	// Copy only consumed flags back to the source, never option values. Call
	// after setup/solve while the source database still exists, so options_left
	// retains its meaning even for options first read by a nested solver.
	void RecordUsed() const
	{
		CollectiveLocalStage(communicator_, "solver option usage", [&] {
			for (const auto& origin : origins_) {
				PetscBool used = PETSC_FALSE;
				petsc_options_detail::Check(PetscOptionsUsed(options_, origin.first.c_str()+1, &used), "PetscOptionsUsed solver snapshot");
				if (!used) continue;
				const char* value = nullptr; PetscBool present = PETSC_FALSE;
				petsc_options_detail::Check(PetscOptionsFindPair(source_, nullptr, origin.second.c_str(), &value, &present), "PetscOptionsFindPair solver source");
			}
		});
	}

private:
	MPI_Comm communicator_;
	PetscOptions source_ = nullptr, options_ = nullptr;
	std::string prefix_;
	std::map<std::string, std::string> origins_;
	bool active_ = false;
};

} // namespace iga
#endif

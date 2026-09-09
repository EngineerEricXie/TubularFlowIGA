#ifndef IGA_PHASE_PROFILE_HPP
#define IGA_PHASE_PROFILE_HPP

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <ostream>

namespace iga {

enum class ProfilePhase : std::size_t {
	Input, Geometry, Assembly, SolverSetup, LinearSolve, Coupling,
	Communication, Output, Diagnostics, Count
};

class PhaseScope;

// Main-thread, rank-local wall-clock accounting. Nested scopes subtract their
// whole duration from the parent. Library-internal MPI stays in its caller's
// phase; Communication denotes only explicitly instrumented application calls.
class PhaseProfile {
public:
	using Now = double (*)();
	struct Counter {
		double inclusive = 0.0;
		double exclusive = 0.0;
		std::size_t calls = 0;
	};

	explicit PhaseProfile(bool enabled = false, Now now = SteadyNow)
		: enabled_(enabled), now_(now), start_(enabled ? now() : 0.0) {}
	PhaseProfile(const PhaseProfile&) = delete;
	PhaseProfile& operator=(const PhaseProfile&) = delete;
	bool Enabled() const noexcept { return enabled_; }

	void EnableFromEnvironment()
	{
		const char* value = std::getenv("IGA_PROFILE");
		enabled_ = value && value[0] == '1' && value[1] == '\0';
		start_ = enabled_ ? now_() : 0.0;
	}

	const Counter& Get(ProfilePhase phase) const
	{
		return counters_.at(static_cast<std::size_t>(phase));
	}

	void Write(std::ostream& out, int rank, int ranks, int status) const
	{
		if (!enabled_) return;
		if (top_) std::terminate();
		constexpr const char* names[] = {"input", "geometry", "assembly",
			"solver_setup", "linear_solve", "coupling", "communication",
			"output", "diagnostics"};
		const double elapsed = now_()-start_;
		double accounted = 0.0;
		const auto precision = out.precision();
		out << std::setprecision(17)
			<< "hpc_profile {\"schema_version\":1,\"rank\":" << rank
			<< ",\"ranks\":" << ranks << ",\"status\":" << status
			<< ",\"elapsed_s\":" << elapsed << ",\"phases\":{";
		for (std::size_t index = 0; index < counters_.size(); ++index) {
			if (index) out << ',';
			const auto& counter = counters_[index];
			accounted += counter.exclusive;
			out << '"' << names[index] << "\":{\"inclusive_s\":" << counter.inclusive
				<< ",\"exclusive_s\":" << counter.exclusive
				<< ",\"calls\":" << counter.calls << '}';
		}
		out << "},\"unscoped_s\":" << elapsed-accounted << "}\n";
		out.precision(precision);
	}

private:
	friend class PhaseScope;
	static double SteadyNow()
	{
		return std::chrono::duration<double>(
			std::chrono::steady_clock::now().time_since_epoch()).count();
	}
	bool enabled_;
	Now now_;
	double start_;
	std::array<Counter, static_cast<std::size_t>(ProfilePhase::Count)> counters_{};
	PhaseScope* top_ = nullptr;
};

inline PhaseProfile& CurrentPhaseProfile()
{
	// Workers start disabled. Threaded assembly should be timed by its enclosing
	// main-thread scope, not by adding overlapping worker wall times.
	static thread_local PhaseProfile profile;
	return profile;
}

class PhaseScope {
public:
	explicit PhaseScope(ProfilePhase phase) : PhaseScope(CurrentPhaseProfile(), phase) {}
	PhaseScope(PhaseProfile& profile, ProfilePhase phase)
		: profile_(profile.enabled_ ? &profile : nullptr), phase_(phase)
	{
		if (!profile_) return;
		parent_ = profile_->top_;
		profile_->top_ = this;
		start_ = profile_->now_();
	}
	PhaseScope(const PhaseScope&) = delete;
	PhaseScope& operator=(const PhaseScope&) = delete;
	~PhaseScope() { Stop(); }

	void Stop() noexcept
	{
		if (!profile_) return;
		if (profile_->top_ != this) std::terminate();
		const double elapsed = profile_->now_()-start_;
		auto& counter = profile_->counters_[static_cast<std::size_t>(phase_)];
		counter.inclusive += elapsed;
		counter.exclusive += elapsed-children_;
		++counter.calls;
		if (parent_) parent_->children_ += elapsed;
		profile_->top_ = parent_;
		profile_ = nullptr;
	}

private:
	PhaseProfile* profile_;
	ProfilePhase phase_;
	PhaseScope* parent_ = nullptr;
	double start_ = 0.0;
	double children_ = 0.0;
};

} // namespace iga

#endif

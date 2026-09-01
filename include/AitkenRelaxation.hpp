#ifndef IGA_AITKEN_RELAXATION_HPP
#define IGA_AITKEN_RELAXATION_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace iga {

enum class AitkenRelaxationStatus {
	Initial,
	Dynamic,
	ClampedMinimum,
	ClampedMaximum,
	TinyDifferenceFallback,
	NonfiniteCandidateFallback
};

inline int AitkenRelaxationStatusCode(AitkenRelaxationStatus status)
{
	return static_cast<int>(status);
}

struct AitkenRelaxationControls {
	double initial_relaxation = 0.5;
	double minimum_relaxation = 0.05;
	double maximum_relaxation = 1.0;
	double scaled_difference_threshold = 64.0*std::numeric_limits<double>::epsilon();
};

inline void ValidateAitkenRelaxationControls(const AitkenRelaxationControls& controls)
{
	if (!std::isfinite(controls.initial_relaxation) || !std::isfinite(controls.minimum_relaxation)
		|| !std::isfinite(controls.maximum_relaxation) || !std::isfinite(controls.scaled_difference_threshold)
		|| !(controls.minimum_relaxation > 0.0) || controls.minimum_relaxation > controls.initial_relaxation
		|| controls.initial_relaxation > controls.maximum_relaxation || controls.maximum_relaxation > 1.0
		|| !(controls.scaled_difference_threshold > 0.0))
		throw std::runtime_error("Aitken controls require 0 < min <= initial <= max <= 1 and finite positive threshold");
}

template <std::size_t N>
struct AitkenRelaxationProposal {
	double relaxation_factor = 0.0;
	double unclamped_relaxation_factor = 0.0;
	double scaled_numerator = 0.0;
	double scaled_denominator = 0.0;
	AitkenRelaxationStatus status = AitkenRelaxationStatus::Initial;
	bool has_previous_residual = false;
	std::array<double, N> next{};
};

template <std::size_t N>
class AitkenRelaxation {
public:
	explicit AitkenRelaxation(AitkenRelaxationControls controls) : controls_(controls)
	{
		ValidateAitkenRelaxationControls(controls_);
	}

	void Reset()
	{
		has_previous_ = false;
		previous_relaxation_ = controls_.initial_relaxation;
	}

	AitkenRelaxationProposal<N> Propose(const std::array<double, N>& x,
		const std::array<double, N>& residual, double reference) const
	{
		if (!(reference > 0.0) || !std::isfinite(reference))
			throw std::runtime_error("Aitken proposal requires positive finite reference pressure");
		for (std::size_t i = 0; i < N; ++i)
			if (!std::isfinite(x[i]) || !std::isfinite(residual[i]))
				throw std::runtime_error("Aitken proposal requires finite iterate and residual");
		AitkenRelaxationProposal<N> proposal;
		proposal.has_previous_residual = has_previous_;
		proposal.relaxation_factor = previous_relaxation_;
		proposal.unclamped_relaxation_factor = previous_relaxation_;
		if (!has_previous_) proposal.status = AitkenRelaxationStatus::Initial;
		else {
			double scale = reference;
			for (std::size_t i = 0; i < N; ++i) scale = std::max(scale, std::max(std::abs(residual[i]), std::abs(previous_residual_[i])));
			double numerator = 0.0, denominator = 0.0;
			for (std::size_t i = 0; i < N; ++i) {
				const double previous_scaled = previous_residual_[i]/scale;
				const double current_scaled = residual[i]/scale;
				const double difference = current_scaled-previous_scaled;
				numerator += previous_scaled*difference;
				denominator += difference*difference;
			}
			proposal.scaled_numerator = numerator;
			proposal.scaled_denominator = denominator;
			if (denominator <= controls_.scaled_difference_threshold) proposal.status = AitkenRelaxationStatus::TinyDifferenceFallback;
			else {
				const double candidate = -previous_relaxation_*numerator/denominator;
				proposal.unclamped_relaxation_factor = candidate;
				if (!std::isfinite(candidate)) { proposal.unclamped_relaxation_factor = previous_relaxation_; proposal.status = AitkenRelaxationStatus::NonfiniteCandidateFallback; }
				else if (candidate < controls_.minimum_relaxation) { proposal.relaxation_factor = controls_.minimum_relaxation; proposal.status = AitkenRelaxationStatus::ClampedMinimum; }
				else if (candidate > controls_.maximum_relaxation) { proposal.relaxation_factor = controls_.maximum_relaxation; proposal.status = AitkenRelaxationStatus::ClampedMaximum; }
				else { proposal.relaxation_factor = std::max(controls_.minimum_relaxation, std::min(controls_.maximum_relaxation, candidate)); proposal.status = AitkenRelaxationStatus::Dynamic; }
			}
		}
		for (std::size_t i = 0; i < N; ++i) {
			proposal.next[i] = x[i]+proposal.relaxation_factor*residual[i];
			if (!std::isfinite(proposal.next[i])) throw std::runtime_error("Aitken proposal produced nonfinite next iterate");
		}
		return proposal;
	}

	void AcceptApplied(const std::array<double, N>& residual, double relaxation)
	{
		for (const auto value : residual) if (!std::isfinite(value)) throw std::runtime_error("Aitken accepted residual must be finite");
		if (!(relaxation >= controls_.minimum_relaxation) || !(relaxation <= controls_.maximum_relaxation) || !std::isfinite(relaxation))
			throw std::runtime_error("Aitken accepted relaxation is outside configured bounds");
		previous_residual_ = residual;
		previous_relaxation_ = relaxation;
		has_previous_ = true;
	}

private:
	AitkenRelaxationControls controls_;
	std::array<double, N> previous_residual_{};
	double previous_relaxation_ = controls_.initial_relaxation;
	bool has_previous_ = false;
};

} // namespace iga

#endif

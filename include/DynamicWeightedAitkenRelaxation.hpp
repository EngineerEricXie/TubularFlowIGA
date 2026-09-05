#ifndef IGA_DYNAMIC_WEIGHTED_AITKEN_RELAXATION_HPP
#define IGA_DYNAMIC_WEIGHTED_AITKEN_RELAXATION_HPP

// Dependency-free distributed dynamic Aitken contract. Each owner keeps its
// unnormalized owned reference-area weights; the caller supplies reductions.
#include "AitkenRelaxation.hpp"
#include "DistributedSurfaceInterface.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

inline void ValidateDynamicAitkenWeights(const std::vector<double>& weights)
{
	if (weights.empty()) throw std::runtime_error("dynamic Aitken weights must be nonempty");
	for (const double weight : weights)
		if (!std::isfinite(weight) || !(weight > 0.0))
			throw std::runtime_error("dynamic Aitken weights must be finite and positive");
}

inline std::string BuildDynamicAitkenWeightIdentitySha256(const std::string& partition_identity_sha256,
	const std::vector<double>& local_weights, double global_weight_total)
{
	if (!IsLowercaseSha256(partition_identity_sha256))
		throw std::runtime_error("dynamic Aitken partition identity must be a lowercase SHA-256 hash");
	ValidateDynamicAitkenWeights(local_weights);
	if (!std::isfinite(global_weight_total) || !(global_weight_total > 0.0))
		throw std::runtime_error("dynamic Aitken global weight total must be finite and positive");
	Sha256 hash;
	distributed_surface_detail::AppendString(hash, "DynamicWeightedAitkenWeights/v2");
	distributed_surface_detail::AppendString(hash, partition_identity_sha256);
	hash.AppendNormalizedDouble(global_weight_total);
	hash.AppendLittleEndian64(static_cast<std::uint64_t>(local_weights.size()));
	for (const double weight : local_weights) hash.AppendNormalizedDouble(weight);
	return hash.Hex();
}

inline double DynamicWeightedAitkenInnerProduct(const std::vector<double>& first,
	const std::vector<double>& second, const std::vector<double>& weights)
{
	if (first.size() != second.size() || first.size() != weights.size())
		throw std::runtime_error("dynamic Aitken inner product vectors must have equal sizes");
	ValidateDynamicAitkenWeights(weights);
	double result = 0.0;
	for (std::size_t i = 0; i < first.size(); ++i) {
		if (!std::isfinite(first[i]) || !std::isfinite(second[i]))
			throw std::runtime_error("dynamic Aitken inner product values must be finite");
		result += weights[i]*first[i]*second[i];
	}
	if (!std::isfinite(result)) throw std::runtime_error("dynamic Aitken inner product is nonfinite");
	return result;
}

struct DynamicWeightedAitkenLocalContributions {
	// Sum these across ranks. The numerator uses accepted r_(k-1), not r_k.
	double numerator = 0.0;
	double denominator = 0.0;
};

struct DynamicWeightedAitkenProposal {
	double relaxation_factor = 0.0;
	double unclamped_relaxation_factor = 0.0;
	double global_numerator = 0.0;
	double global_denominator = 0.0;
	double residual_scale = 0.0;
	AitkenRelaxationStatus status = AitkenRelaxationStatus::Initial;
	bool has_previous_residual = false;
	// This identity covers only globally comparable control state. It must not
	// include partition-local residual history or ownership weights.
	std::string control_state_identity_sha256;
	// Globally identical proposal facts. This is safe to fold into the pending
	// collective control-state identity.
	std::string proposal_identity_sha256;
	// Locally exact transaction facts, including the bound partition and local
	// vectors. This is intentionally not collectively comparable.
	std::string local_proposal_identity_sha256;
	std::string partition_identity_sha256;
	std::string weight_identity_sha256;
	std::uint64_t control_generation = 0;
	std::uint64_t accepted_iteration_count = 0;
	std::vector<double> iterate;
	std::vector<double> residual;
	std::vector<double> update;
	std::vector<double> next;
};

class DynamicWeightedAitkenRelaxation {
public:
	DynamicWeightedAitkenRelaxation(std::string partition_identity_sha256,
		std::vector<double> local_weights, double global_weight_total,
		AitkenRelaxationControls controls = {})
		: partition_identity_sha256_(std::move(partition_identity_sha256)),
		  weight_identity_sha256_(BuildDynamicAitkenWeightIdentitySha256(
			  partition_identity_sha256_, local_weights, global_weight_total)),
		  controls_(controls), local_weights_(std::move(local_weights)),
		  global_weight_total_(global_weight_total)
	{
		ValidateAitkenRelaxationControls(controls_);
		for (const double weight : local_weights_) local_weight_total_ += weight;
		if (!std::isfinite(local_weight_total_) || !(local_weight_total_ > 0.0)
			|| local_weight_total_ > global_weight_total_)
			throw std::runtime_error("dynamic Aitken local weights are inconsistent with global total");
		previous_residual_.resize(local_weights_.size());
	}

	const std::string& PartitionIdentitySha256() const noexcept { return partition_identity_sha256_; }
	const std::string& WeightIdentitySha256() const noexcept { return weight_identity_sha256_; }
	const std::vector<double>& LocalWeights() const noexcept { return local_weights_; }
	double GlobalWeightTotal() const noexcept { return global_weight_total_; }
	double LocalWeightTotal() const noexcept { return local_weight_total_; }

	// Every rank in one coupling operation must have this same identity before
	// it participates in a global proposal. It intentionally excludes the
	// partition-local previous residual vector.
	std::string ControlStateIdentitySha256() const
	{
		Sha256 hash;
		distributed_surface_detail::AppendString(hash, "DynamicWeightedAitkenControlState/v2");
		hash.AppendLittleEndian64(control_generation_);
		hash.AppendLittleEndian64(accepted_iteration_count_);
		hash.AppendLittleEndian32(has_previous_ ? 1U : 0U);
		hash.AppendNormalizedDouble(previous_relaxation_);
		hash.AppendNormalizedDouble(controls_.initial_relaxation);
		hash.AppendNormalizedDouble(controls_.minimum_relaxation);
		hash.AppendNormalizedDouble(controls_.maximum_relaxation);
		hash.AppendNormalizedDouble(controls_.scaled_difference_threshold);
		hash.AppendNormalizedDouble(global_weight_total_);
		hash.AppendLittleEndian32(pending_proposal_ ? 1U : 0U);
		if (pending_proposal_)
			distributed_surface_detail::AppendString(hash, pending_proposal_->proposal_identity_sha256);
		return hash.Hex();
	}

	void Reset() noexcept
	{
		has_previous_ = false;
		previous_relaxation_ = controls_.initial_relaxation;
		accepted_iteration_count_ = 0;
		++control_generation_;
		ClearPendingProposal();
	}

	// Globally reduce max(reference_scale, LocalRequiredResidualScale(...)).
	double LocalRequiredResidualScale(const std::vector<double>& residual,
		double reference_scale, const std::string& partition_identity_sha256) const
	{
		ValidateResidualBinding(residual, partition_identity_sha256);
		ValidateReferenceScale(reference_scale);
		double scale = reference_scale;
		for (const double value : residual) scale = std::max(scale, std::abs(value));
		if (has_previous_)
			for (const double value : previous_residual_) scale = std::max(scale, std::abs(value));
		return scale;
	}

	DynamicWeightedAitkenLocalContributions LocalContributions(const std::vector<double>& residual,
		const std::string& partition_identity_sha256) const
	{
		ValidateResidualBinding(residual, partition_identity_sha256);
		DynamicWeightedAitkenLocalContributions result;
		if (!has_previous_) return result;
		for (std::size_t i = 0; i < residual.size(); ++i) {
			const double difference = residual[i]-previous_residual_[i];
			result.numerator += local_weights_[i]*previous_residual_[i]*difference;
			result.denominator += local_weights_[i]*difference*difference;
		}
		if (!std::isfinite(result.numerator) || !std::isfinite(result.denominator))
			throw std::runtime_error("dynamic Aitken local reduction contributions are nonfinite");
		return result;
	}

	// Every rank supplies the same globally reduced scalar terms and scale.
	DynamicWeightedAitkenProposal ProposeWithGlobalReduction(const std::vector<double>& iterate,
		const std::vector<double>& residual, double reference_scale, double global_residual_scale,
		double global_numerator, double global_denominator,
		const std::string& expected_control_state_identity_sha256,
		const std::string& partition_identity_sha256)
	{
		ValidateBinding(iterate, residual, partition_identity_sha256);
		if (pending_proposal_)
			throw std::runtime_error("dynamic Aitken has an unapplied pending proposal");
		if (!IsLowercaseSha256(expected_control_state_identity_sha256))
			throw std::runtime_error("dynamic Aitken expected control-state identity must be a lowercase SHA-256 hash");
		const std::string control_identity = ControlStateIdentitySha256();
		if (expected_control_state_identity_sha256 != control_identity)
			throw std::runtime_error("dynamic Aitken control state is desynchronized across ranks");
		ValidateReferenceScale(reference_scale);
		if (!std::isfinite(global_residual_scale) || global_residual_scale < reference_scale)
			throw std::runtime_error("dynamic Aitken global residual scale is invalid");
		if (!std::isfinite(global_numerator) || !std::isfinite(global_denominator)
			|| global_denominator < 0.0)
			throw std::runtime_error("dynamic Aitken global reduction values are invalid");
		if (global_residual_scale < LocalRequiredResidualScale(residual, reference_scale, partition_identity_sha256))
			throw std::runtime_error("dynamic Aitken global residual scale omits local residual data");

		DynamicWeightedAitkenProposal proposal;
		proposal.partition_identity_sha256 = partition_identity_sha256_;
		proposal.weight_identity_sha256 = weight_identity_sha256_;
		proposal.iterate = iterate;
		proposal.residual = residual;
		proposal.update.resize(iterate.size());
		proposal.next.resize(iterate.size());
		proposal.has_previous_residual = has_previous_;
		proposal.relaxation_factor = previous_relaxation_;
		proposal.unclamped_relaxation_factor = previous_relaxation_;
		proposal.global_numerator = global_numerator;
		proposal.global_denominator = global_denominator;
		proposal.residual_scale = global_residual_scale;
		proposal.control_state_identity_sha256 = control_identity;
		proposal.control_generation = control_generation_;
		proposal.accepted_iteration_count = accepted_iteration_count_;
		if (!has_previous_) proposal.status = AitkenRelaxationStatus::Initial;
		else {
			const double normalized_difference = (global_denominator/global_weight_total_)
				/(global_residual_scale*global_residual_scale);
			if (!std::isfinite(normalized_difference))
				proposal.status = AitkenRelaxationStatus::NonfiniteCandidateFallback;
			else if (normalized_difference <= controls_.scaled_difference_threshold)
				proposal.status = AitkenRelaxationStatus::TinyDifferenceFallback;
			else {
				const double candidate = -previous_relaxation_*global_numerator/global_denominator;
				proposal.unclamped_relaxation_factor = candidate;
				if (!std::isfinite(candidate)) {
					proposal.unclamped_relaxation_factor = previous_relaxation_;
					proposal.status = AitkenRelaxationStatus::NonfiniteCandidateFallback;
				} else if (candidate < controls_.minimum_relaxation) {
					proposal.relaxation_factor = controls_.minimum_relaxation;
					proposal.status = AitkenRelaxationStatus::ClampedMinimum;
				} else if (candidate > controls_.maximum_relaxation) {
					proposal.relaxation_factor = controls_.maximum_relaxation;
					proposal.status = AitkenRelaxationStatus::ClampedMaximum;
				} else {
					proposal.relaxation_factor = candidate;
					proposal.status = AitkenRelaxationStatus::Dynamic;
				}
			}
		}
		for (std::size_t i = 0; i < iterate.size(); ++i) {
			proposal.update[i] = proposal.relaxation_factor*residual[i];
			proposal.next[i] = iterate[i]+proposal.update[i];
			if (!std::isfinite(proposal.next[i]))
				throw std::runtime_error("dynamic Aitken proposal produced nonfinite next iterate");
		}
		proposal.proposal_identity_sha256 = BuildGlobalProposalIdentitySha256(proposal);
		proposal.local_proposal_identity_sha256 = BuildLocalProposalIdentitySha256(proposal);
		pending_proposal_ = proposal;
		return proposal;
	}

	// Exact single-rank path: it delegates to the global-reduction method.
	DynamicWeightedAitkenProposal Propose(const std::vector<double>& iterate,
		const std::vector<double>& residual, double reference_scale,
		const std::string& partition_identity_sha256)
	{
		if (local_weight_total_ != global_weight_total_)
			throw std::runtime_error("distributed dynamic Aitken requires global reduction");
		const double scale = LocalRequiredResidualScale(residual, reference_scale, partition_identity_sha256);
		const auto terms = LocalContributions(residual, partition_identity_sha256);
		return ProposeWithGlobalReduction(iterate, residual, reference_scale, scale,
			terms.numerator, terms.denominator, ControlStateIdentitySha256(), partition_identity_sha256);
	}

	void AcceptApplied(const DynamicWeightedAitkenProposal& proposal,
		const std::vector<double>& residual, double relaxation,
		const std::string& partition_identity_sha256)
	{
		if (partition_identity_sha256 != partition_identity_sha256_)
			throw std::runtime_error("dynamic Aitken cannot accept a changed partition identity");
		if (residual.size() != previous_residual_.size())
			throw std::runtime_error("dynamic Aitken cannot accept a changed residual size");
		if (!pending_proposal_
			|| proposal.proposal_identity_sha256 != BuildGlobalProposalIdentitySha256(proposal)
			|| proposal.local_proposal_identity_sha256 != BuildLocalProposalIdentitySha256(proposal)
			|| !ProposalExactlyEquals(proposal, *pending_proposal_)
			|| !DoubleVectorExactlyEquals(residual, proposal.residual))
			throw std::runtime_error("dynamic Aitken cannot accept a stale or foreign proposal");
		if (!DoubleExactlyEquals(relaxation, proposal.relaxation_factor))
			throw std::runtime_error("dynamic Aitken accepted relaxation must exactly match its pending proposal");
		if (!(relaxation >= controls_.minimum_relaxation) || !(relaxation <= controls_.maximum_relaxation)
			|| !std::isfinite(relaxation))
			throw std::runtime_error("dynamic Aitken accepted relaxation is outside configured bounds");
		for (const double value : residual)
			if (!std::isfinite(value)) throw std::runtime_error("dynamic Aitken accepted residual must be finite");
		std::copy(residual.begin(), residual.end(), previous_residual_.begin());
		previous_relaxation_ = relaxation;
		has_previous_ = true;
		++accepted_iteration_count_;
		ClearPendingProposal();
	}

private:
	static std::uint64_t DoubleBits(double value) noexcept
	{
		std::uint64_t bits = 0;
		static_assert(sizeof(bits) == sizeof(value), "unsupported double representation");
		std::memcpy(&bits, &value, sizeof(bits));
		return bits;
	}

	static void AppendExactDouble(Sha256& hash, double value)
	{
		hash.AppendLittleEndian64(DoubleBits(value));
	}

	static void AppendExactDoubleVector(Sha256& hash, const std::vector<double>& values)
	{
		hash.AppendLittleEndian64(static_cast<std::uint64_t>(values.size()));
		for (const double value : values) AppendExactDouble(hash, value);
	}

	static bool DoubleExactlyEquals(double first, double second) noexcept
	{
		return DoubleBits(first) == DoubleBits(second);
	}

	static bool DoubleVectorExactlyEquals(const std::vector<double>& first,
		const std::vector<double>& second) noexcept
	{
		if (first.size() != second.size()) return false;
		for (std::size_t i = 0; i < first.size(); ++i)
			if (!DoubleExactlyEquals(first[i], second[i])) return false;
		return true;
	}

	static std::string BuildGlobalProposalIdentitySha256(const DynamicWeightedAitkenProposal& proposal)
	{
		Sha256 hash;
		distributed_surface_detail::AppendString(hash, "DynamicWeightedAitkenGlobalProposal/v2");
		distributed_surface_detail::AppendString(hash, proposal.control_state_identity_sha256);
		hash.AppendLittleEndian64(proposal.control_generation);
		hash.AppendLittleEndian64(proposal.accepted_iteration_count);
		AppendExactDouble(hash, proposal.relaxation_factor);
		AppendExactDouble(hash, proposal.unclamped_relaxation_factor);
		AppendExactDouble(hash, proposal.global_numerator);
		AppendExactDouble(hash, proposal.global_denominator);
		AppendExactDouble(hash, proposal.residual_scale);
		hash.AppendLittleEndian32(static_cast<std::uint32_t>(proposal.status));
		hash.AppendLittleEndian32(proposal.has_previous_residual ? 1U : 0U);
		return hash.Hex();
	}

	static std::string BuildLocalProposalIdentitySha256(const DynamicWeightedAitkenProposal& proposal)
	{
		Sha256 hash;
		distributed_surface_detail::AppendString(hash, "DynamicWeightedAitkenLocalProposal/v2");
		distributed_surface_detail::AppendString(hash, proposal.proposal_identity_sha256);
		distributed_surface_detail::AppendString(hash, proposal.partition_identity_sha256);
		distributed_surface_detail::AppendString(hash, proposal.weight_identity_sha256);
		AppendExactDoubleVector(hash, proposal.iterate);
		AppendExactDoubleVector(hash, proposal.residual);
		AppendExactDoubleVector(hash, proposal.update);
		AppendExactDoubleVector(hash, proposal.next);
		return hash.Hex();
	}

	static bool ProposalExactlyEquals(const DynamicWeightedAitkenProposal& first,
		const DynamicWeightedAitkenProposal& second) noexcept
	{
		return DoubleExactlyEquals(first.relaxation_factor, second.relaxation_factor)
			&& DoubleExactlyEquals(first.unclamped_relaxation_factor, second.unclamped_relaxation_factor)
			&& DoubleExactlyEquals(first.global_numerator, second.global_numerator)
			&& DoubleExactlyEquals(first.global_denominator, second.global_denominator)
			&& DoubleExactlyEquals(first.residual_scale, second.residual_scale)
			&& first.status == second.status
			&& first.has_previous_residual == second.has_previous_residual
			&& first.control_state_identity_sha256 == second.control_state_identity_sha256
			&& first.proposal_identity_sha256 == second.proposal_identity_sha256
			&& first.local_proposal_identity_sha256 == second.local_proposal_identity_sha256
			&& first.partition_identity_sha256 == second.partition_identity_sha256
			&& first.weight_identity_sha256 == second.weight_identity_sha256
			&& first.control_generation == second.control_generation
			&& first.accepted_iteration_count == second.accepted_iteration_count
			&& DoubleVectorExactlyEquals(first.iterate, second.iterate)
			&& DoubleVectorExactlyEquals(first.residual, second.residual)
			&& DoubleVectorExactlyEquals(first.update, second.update)
			&& DoubleVectorExactlyEquals(first.next, second.next);
	}

	void ClearPendingProposal() noexcept
	{
		pending_proposal_.reset();
	}

	static void ValidateReferenceScale(double reference_scale)
	{
		if (!(reference_scale > 0.0) || !std::isfinite(reference_scale))
			throw std::runtime_error("dynamic Aitken reference scale must be finite and positive");
	}

	void ValidateResidualBinding(const std::vector<double>& residual,
		const std::string& partition_identity_sha256) const
	{
		if (partition_identity_sha256 != partition_identity_sha256_)
			throw std::runtime_error("dynamic Aitken requires its bound partition identity");
		if (residual.size() != local_weights_.size())
			throw std::runtime_error("dynamic Aitken residual size must match fixed local weights");
		for (const double value : residual)
			if (!std::isfinite(value)) throw std::runtime_error("dynamic Aitken residual must be finite");
	}

	void ValidateBinding(const std::vector<double>& iterate, const std::vector<double>& residual,
		const std::string& partition_identity_sha256) const
	{
		ValidateResidualBinding(residual, partition_identity_sha256);
		if (iterate.size() != local_weights_.size())
			throw std::runtime_error("dynamic Aitken iterate size must match fixed local weights");
		for (const double value : iterate)
			if (!std::isfinite(value)) throw std::runtime_error("dynamic Aitken iterate must be finite");
	}

	std::string partition_identity_sha256_;
	std::string weight_identity_sha256_;
	AitkenRelaxationControls controls_;
	std::vector<double> local_weights_;
	std::vector<double> previous_residual_;
	double global_weight_total_ = 0.0;
	double local_weight_total_ = 0.0;
	double previous_relaxation_ = controls_.initial_relaxation;
	bool has_previous_ = false;
	std::uint64_t control_generation_ = 0;
	std::uint64_t accepted_iteration_count_ = 0;
	std::optional<DynamicWeightedAitkenProposal> pending_proposal_;
};

} // namespace iga

#endif

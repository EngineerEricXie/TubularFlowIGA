#ifndef IGA_DISTRIBUTED_WEIGHTED_AITKEN_HPP
#define IGA_DISTRIBUTED_WEIGHTED_AITKEN_HPP

#include "DynamicWeightedAitkenRelaxation.hpp"
#include "CollectiveFailure.hpp"
#include <memory>

namespace iga {
// Collective state owner. Proposals and acceptances are prepared on a copy;
// only a successful collective validation swaps it into the authoritative
// state. The coordinator separately owns application of the proposed fields.
class DistributedWeightedAitken {
public:
	DistributedWeightedAitken(MPI_Comm comm,const std::string& partition_identity,
		std::vector<double> local_weights,AitkenRelaxationControls controls={}) : comm_(comm)
	{
		long double local=0.,global=0.;
		CollectiveLocalStage(comm_,"distributed Aitken weights",[&] {
			ValidateDynamicAitkenWeights(local_weights);
			for(double weight:local_weights)local+=static_cast<long double>(weight);
			if(!std::isfinite(local))throw std::overflow_error("Aitken weight sum overflow");
		});
		MPI_Allreduce(&local,&global,1,MPI_LONG_DOUBLE,MPI_SUM,comm_);
		std::string control;
		CollectiveLocalStage(comm_,"distributed Aitken state",[&] {
			if(!std::isfinite(global)||!(global>0.)||global>std::numeric_limits<double>::max())
				throw std::invalid_argument("invalid global Aitken weight sum");
			state_=std::make_unique<DynamicWeightedAitkenRelaxation>(partition_identity,std::move(local_weights),static_cast<double>(global),controls);
			control=state_->ControlStateIdentitySha256();
		});
		RequireCollectiveSameText(comm_,"distributed Aitken initial controls",control);
	}
	std::string ControlIdentity() const { return state_->ControlStateIdentitySha256(); }
	double GlobalWeightTotal() const noexcept { return state_->GlobalWeightTotal(); }

	DynamicWeightedAitkenProposal Propose(const std::vector<double>& iterate,
		const std::vector<double>& residual,double reference_scale)
	{
		std::string control,reference_identity;
		double local[2]{},global[2]{},scale=0.,global_scale=0.;
		CollectiveLocalStage(comm_,"distributed Aitken local terms",[&] {
			control=state_->ControlStateIdentitySha256();
			Sha256 hash;hash.AppendNormalizedDouble(reference_scale);reference_identity=hash.Hex();
			const auto terms=state_->LocalContributions(residual,state_->PartitionIdentitySha256());
			local[0]=terms.numerator;local[1]=terms.denominator;
			scale=state_->LocalRequiredResidualScale(residual,reference_scale,state_->PartitionIdentitySha256());
		});
		RequireCollectiveSameText(comm_,"distributed Aitken control agreement",control);
		RequireCollectiveSameText(comm_,"distributed Aitken reference scale",reference_identity);
		MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,comm_);
		MPI_Allreduce(&scale,&global_scale,1,MPI_DOUBLE,MPI_MAX,comm_);
		std::unique_ptr<DynamicWeightedAitkenRelaxation> candidate;
		DynamicWeightedAitkenProposal proposal;std::string pending;
		CollectiveLocalStage(comm_,"distributed Aitken proposal",[&] {
			candidate=std::make_unique<DynamicWeightedAitkenRelaxation>(*state_);
			proposal=candidate->ProposeWithGlobalReduction(iterate,residual,reference_scale,global_scale,global[0],global[1],control,state_->PartitionIdentitySha256());
			pending=candidate->ControlStateIdentitySha256();
		});
		RequireCollectiveSameText(comm_,"distributed Aitken pending agreement",pending);
		state_.swap(candidate);return proposal;
	}
	void AcceptApplied(const DynamicWeightedAitkenProposal& proposal,
		const std::vector<double>& applied_residual,double applied_relaxation)
	{
		std::unique_ptr<DynamicWeightedAitkenRelaxation> candidate;std::string accepted;
		CollectiveLocalStage(comm_,"distributed Aitken acceptance",[&] {
			candidate=std::make_unique<DynamicWeightedAitkenRelaxation>(*state_);
			candidate->AcceptApplied(proposal,applied_residual,applied_relaxation,state_->PartitionIdentitySha256());
			accepted=candidate->ControlStateIdentitySha256();
		});
		RequireCollectiveSameText(comm_,"distributed Aitken accepted agreement",accepted);state_.swap(candidate);
	}
	void Reset()
	{
		std::unique_ptr<DynamicWeightedAitkenRelaxation> candidate;std::string reset;
		CollectiveLocalStage(comm_,"distributed Aitken reset",[&] {
			candidate=std::make_unique<DynamicWeightedAitkenRelaxation>(*state_);candidate->Reset();reset=candidate->ControlStateIdentitySha256();
		});
		RequireCollectiveSameText(comm_,"distributed Aitken reset agreement",reset);state_.swap(candidate);
	}
private:
	MPI_Comm comm_;std::unique_ptr<DynamicWeightedAitkenRelaxation> state_;
};
} // namespace iga
#endif

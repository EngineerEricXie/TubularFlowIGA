#ifndef IGA_DISTRIBUTED_FSI_CONVERGENCE_HPP
#define IGA_DISTRIBUTED_FSI_CONVERGENCE_HPP

#include "CollectiveFailure.hpp"
#include "Sha256.hpp"
#include <cmath>
#include <limits>
#include <vector>

namespace iga {

struct DistributedFsiConvergence
{
	double area_weighted_rms_residual_m = 0.;
	double max_residual_m = 0.;
	double displacement_scale_m = 0.;
	double convergence_threshold_m = 0.;
	bool converged = false;
};

// Inputs contain uniquely owned scalar normal displacements only. The caller
// must have validated surface ownership and a common trial epoch beforehand.
// Empty partitions are legal; an entirely empty interface is rejected.
inline DistributedFsiConvergence EvaluateDistributedFsiConvergence(MPI_Comm comm,
	const std::vector<double>& owned_areas, const std::vector<double>& residual,
	const std::vector<double>& raw_displacement, const std::vector<double>& current_displacement,
	double reference_scale, double absolute_tolerance, double relative_tolerance)
{
	double local_max[3]{}, global_max[3]{};
	std::string controls;
	CollectiveLocalStage(comm, "FSI convergence inputs", [&] {
		if (owned_areas.size()!=residual.size() || residual.size()!=raw_displacement.size()
			|| residual.size()!=current_displacement.size())
			throw std::invalid_argument("FSI convergence tuple count mismatch");
		if (!std::isfinite(reference_scale) || reference_scale<=0.
			|| !std::isfinite(absolute_tolerance) || absolute_tolerance<0.
			|| !std::isfinite(relative_tolerance) || relative_tolerance<0.)
			throw std::invalid_argument("invalid FSI convergence controls");
		Sha256 hash;
		hash.AppendNormalizedDouble(reference_scale);
		hash.AppendNormalizedDouble(absolute_tolerance);
		hash.AppendNormalizedDouble(relative_tolerance);
		controls=hash.Hex();
		for (std::size_t i=0;i<residual.size();++i) {
			if (!std::isfinite(owned_areas[i]) || owned_areas[i]<=0.
				|| !std::isfinite(residual[i]) || !std::isfinite(raw_displacement[i])
				|| !std::isfinite(current_displacement[i]))
				throw std::invalid_argument("invalid owned FSI convergence value");
			local_max[0]=std::max(local_max[0],owned_areas[i]);
			local_max[1]=std::max(local_max[1],std::abs(residual[i]));
			local_max[2]=std::max(local_max[2],std::max(std::abs(raw_displacement[i]),std::abs(current_displacement[i])));
		}
	});
	RequireCollectiveSameText(comm,"FSI convergence control agreement",controls);
	MPI_Allreduce(local_max,global_max,3,MPI_DOUBLE,MPI_MAX,comm);
	long double local_sum[2]{}, global_sum[2]{};
	CollectiveLocalStage(comm,"FSI convergence scaled sums",[&] {
		if (!(global_max[0]>0.)) throw std::invalid_argument("empty global FSI interface");
		for (std::size_t i=0;i<residual.size();++i) {
			const long double weight=static_cast<long double>(owned_areas[i])/global_max[0];
			const long double scaled=global_max[1]>0. ? static_cast<long double>(residual[i])/global_max[1] : 0.;
			local_sum[0]+=weight;
			local_sum[1]+=weight*scaled*scaled;
		}
		if (!std::isfinite(local_sum[0]) || !std::isfinite(local_sum[1]))
			throw std::overflow_error("FSI convergence local sum overflow");
	});
	MPI_Allreduce(local_sum,global_sum,2,MPI_LONG_DOUBLE,MPI_SUM,comm);
	DistributedFsiConvergence result;
	CollectiveLocalStage(comm,"FSI convergence decision",[&] {
		if (!std::isfinite(global_sum[0]) || !(global_sum[0]>0.) || !std::isfinite(global_sum[1]))
			throw std::overflow_error("invalid global FSI convergence sums");
		// Normalization avoids forming area * residual squared in double.
		const long double ratio=std::min(1.L,global_sum[1]/global_sum[0]);
		result.area_weighted_rms_residual_m=static_cast<double>(global_max[1]*std::sqrt(ratio));
		result.max_residual_m=global_max[1];
		result.displacement_scale_m=std::max(reference_scale,global_max[2]);
		const long double threshold=static_cast<long double>(absolute_tolerance)
			+static_cast<long double>(relative_tolerance)*result.displacement_scale_m;
		if (!std::isfinite(threshold) || threshold>std::numeric_limits<double>::max())
			throw std::overflow_error("FSI convergence threshold overflow");
		result.convergence_threshold_m=static_cast<double>(threshold);
		// Preserve the serial coordinator's RMS stopping rule; max is diagnostic.
		result.converged=result.area_weighted_rms_residual_m<=result.convergence_threshold_m;
	});
	return result;
}

} // namespace iga
#endif

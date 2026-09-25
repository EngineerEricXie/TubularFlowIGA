#ifndef IGA_NATIVE_TET_ALE_STEP_CONTROL_HPP
#define IGA_NATIVE_TET_ALE_STEP_CONTROL_HPP

#include "NativeTetAleKinematics.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace iga {

struct NativeTetAleRetryPolicy
{
	std::size_t maximum_attempts = 5;
	double reduction_factor = 0.5;
	double minimum_dt_s = 1.0e-8;
	double minimum_determinant_ratio = 0.05;
	double minimum_scaled_jacobian = 1.0e-3;
};

struct NativeTetAleRetryResult
{
	std::size_t attempts = 0;
	std::size_t rejected_attempts = 0;
	double accepted_dt_s = 0.0;
};

// motion_at_time(t) must return absolute displacement relative to X. This
// begins a valid trial but never commits it; the coupled/fluid solve owns
// acceptance and calls CommitTrial only after all of its gates pass.
template<class MotionAtTime>
NativeTetAleRetryResult BeginNativeTetAleQualityControlledTrial(
	NativeTetAleKinematics& kinematics,double requested_next_time_s,
	MotionAtTime&& motion_at_time,const NativeTetAleRetryPolicy& policy = {})
{
	if (policy.maximum_attempts == 0 || !(policy.reduction_factor > 0.0)
		|| !(policy.reduction_factor < 1.0) || !(policy.minimum_dt_s > 0.0)
		|| !std::isfinite(policy.minimum_dt_s))
		throw std::invalid_argument("native ALE retry policy is invalid");
	const double committed_time = kinematics.CommittedTime();
	double dt = requested_next_time_s-committed_time;
	if (!(dt > 0.0) || !std::isfinite(dt))
		throw std::invalid_argument("native ALE requested time does not advance history");
	NativeTetAleRetryResult result;
	std::string last_failure;
	for (std::size_t attempt = 1; attempt <= policy.maximum_attempts; ++attempt) {
		if (dt < policy.minimum_dt_s) break;
		result.attempts = attempt;
		try {
			kinematics.BeginTrial(motion_at_time(committed_time+dt),committed_time+dt,
				policy.minimum_determinant_ratio,policy.minimum_scaled_jacobian);
			result.accepted_dt_s = dt;
			return result;
		} catch (const std::runtime_error& error) {
			last_failure = error.what();
			++result.rejected_attempts;
			kinematics.RejectTrial();
			dt *= policy.reduction_factor;
		}
	}
	throw std::runtime_error("native ALE quality retries exhausted: "+last_failure);
}

} // namespace iga

#endif

#ifndef IGA_COUPLED_DOMAIN_RUNTIME_HPP
#define IGA_COUPLED_DOMAIN_RUNTIME_HPP

#include "CouplingEdge.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct DomainStepContext {
	int step_index = 0;
	double start_time_s = 0.0;
	double dt_s = 0.0;

	double EndTime() const { return start_time_s+dt_s; }

	void Validate() const
	{
		if (step_index < 0 || !std::isfinite(start_time_s)
			|| !(dt_s > 0.0) || !std::isfinite(dt_s) || !std::isfinite(EndTime()))
			throw std::runtime_error("domain step requires a nonnegative index, finite start time, and positive finite dt");
	}
};

class CoupledDomainRuntime {
public:
	virtual ~CoupledDomainRuntime() = default;

	virtual const std::string& DomainId() const noexcept = 0;
	virtual DomainKind Kind() const noexcept = 0;
	virtual const std::vector<CouplingPort>& Ports() const noexcept = 0;

	virtual void BeginStep(const DomainStepContext& step) = 0;
	virtual void SetPortInput(const std::string& port_id,
		const PortBoundaryData& input) = 0;
	virtual void SolveTrial() = 0;
	virtual PortState GetPortState(const std::string& port_id) const = 0;
	virtual void RollbackTrial() = 0;
	virtual void AbortStep() = 0;
	virtual void PrepareCommitStep() = 0;
	virtual void FinalizeCommitStep() noexcept = 0;
};

} // namespace iga

#endif

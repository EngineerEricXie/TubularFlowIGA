#ifndef IGA_COUPLED_DOMAIN_RUNTIME_HPP
#define IGA_COUPLED_DOMAIN_RUNTIME_HPP

#include "CouplingEdge.hpp"

#include <cmath>
#include <map>
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

// This deliberately sits beside CoupledDomainRuntime rather than extending it:
// pressure/flow-only executors retain their established one-shot trial contract.
// A flow/transport executor can instead hold an accepted hydraulic trial while
// retrying scalar boundary data against its fixed trial velocity.
struct SpeciesStepAccounting {
	double initial_mass = 0.0;
	double final_mass = 0.0;
	double source_amount = 0.0;
	std::map<std::string, double> outward_port_amount;
	double residual = 0.0;
};

class StagedFlowTransportDomainRuntime {
public:
	virtual ~StagedFlowTransportDomainRuntime() = default;

	virtual void SolveHydraulicTrial() = 0;
	virtual PortState GetHydraulicPortState(const std::string& port_id) const = 0;
	virtual void RollbackHydraulicTrial() = 0;

	virtual void SetTransportConcentration(const std::string& port_id,
		double time_s, const std::map<std::string, double>& concentration) = 0;
	virtual void SolveTransportTrial() = 0;
	virtual PortState GetTransportPortState(const std::string& port_id) const = 0;
	virtual void RollbackTransportTrial() = 0;
	virtual std::map<std::string, SpeciesStepAccounting>
	GetSpeciesStepAccounting() const = 0;
};

} // namespace iga

#endif

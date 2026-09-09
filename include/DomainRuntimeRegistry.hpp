#ifndef IGA_DOMAIN_RUNTIME_REGISTRY_HPP
#define IGA_DOMAIN_RUNTIME_REGISTRY_HPP

#include "CoupledDomainRuntime.hpp"
#include "SimulationGraph.hpp"

#include <algorithm>
#include <exception>
#include <functional>
#include <type_traits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

inline bool SameCouplingPort(const CouplingPort& first, const CouplingPort& second)
{
	return first.id == second.id && first.subsystem_id == second.subsystem_id
		&& first.locator_kind == second.locator_kind && first.locator == second.locator
		&& first.orientation.native_to_outward_sign == second.orientation.native_to_outward_sign
		&& first.provides == second.provides && first.requires == second.requires
		&& first.species == second.species;
}

class DomainRuntimeRegistry {
public:
	using ConstructionOutcome = std::function<void(std::exception_ptr)>;

	// Validate while the caller still owns every runtime. Distributed callers
	// synchronize the local outcome before the noexcept ownership transfer.
	// On failure runtimes remains untouched, including when indexing allocates.
	DomainRuntimeRegistry(const SimulationGraph& graph,
		std::vector<std::unique_ptr<CoupledDomainRuntime>>&& runtimes,
		const ConstructionOutcome& synchronize = {})
		: graph_(graph)
	{
		static_assert(std::is_nothrow_move_assignable<decltype(owned_)>::value,
			"registry ownership transfer must not throw");
		std::exception_ptr error;
		try { ValidateAndIndex(runtimes); }
		catch (...) { error = std::current_exception(); }
		if (synchronize) synchronize(error);
		if (error) std::rethrow_exception(error);
		owned_ = std::move(runtimes);
	}

	CoupledDomainRuntime& Runtime(const std::string& domain_id) const
	{
		const auto found = by_id_.find(domain_id);
		if (found == by_id_.end())
			throw std::runtime_error("runtime registry has no domain '"+domain_id+"'");
		return *found->second;
	}

	const SimulationGraph& Graph() const noexcept { return graph_; }

private:
	void ValidateAndIndex(const std::vector<std::unique_ptr<CoupledDomainRuntime>>& runtimes)
	{
		if (runtimes.size() != graph_.Domains().size())
			throw std::runtime_error("runtime registry requires exactly one runtime per graph domain");
		for (const auto& runtime : runtimes) {
			if (!runtime) throw std::runtime_error("runtime registry cannot own a null runtime");
			const auto& node = graph_.Domain(runtime->DomainId());
			if (node.kind != runtime->Kind())
				throw std::runtime_error("runtime kind does not match graph domain '"+node.id+"'");
			if (node.ports.size() != runtime->Ports().size())
				throw std::runtime_error("runtime ports do not match graph domain '"+node.id+"'");
			for (const auto& port : node.ports) {
				const auto found = std::find_if(runtime->Ports().begin(), runtime->Ports().end(),
					[&](const CouplingPort& candidate) { return candidate.id == port.id; });
				if (found == runtime->Ports().end() || !SameCouplingPort(port, *found))
					throw std::runtime_error("runtime port metadata does not match graph port '"
						+node.id+"."+port.id+"'");
			}
			if (!by_id_.emplace(node.id, runtime.get()).second)
				throw std::runtime_error("runtime registry domain ids must be unique");
		}
	}

	const SimulationGraph& graph_;
	std::vector<std::unique_ptr<CoupledDomainRuntime>> owned_;
	std::map<std::string, CoupledDomainRuntime*> by_id_;
};

} // namespace iga

#endif

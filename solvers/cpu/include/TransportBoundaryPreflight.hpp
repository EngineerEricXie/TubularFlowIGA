#ifndef IGA_TRANSPORT_BOUNDARY_PREFLIGHT_HPP
#define IGA_TRANSPORT_BOUNDARY_PREFLIGHT_HPP

#include "OwnedRowAssembler.hpp"
#include "SimulationConfig.hpp"

#include <map>
#include <stdexcept>
#include <string>

namespace iga {

inline std::map<int, long long> CountConfiguredScalarSurfaceFaces(
	const SimulationConfiguration& configuration, const CompiledLinearSystem& system,
	const OwnedRowAssembler& assembler, MPI_Comm communicator)
{
	std::map<int, long long> counts;
	for (const auto& boundary : configuration.boundaries)
		for (const auto& condition : boundary.conditions)
			if (system.field_index.count(condition.field)
				&& (condition.kind == FieldBoundaryKind::Flux
					|| condition.kind == FieldBoundaryKind::Robin))
				counts.emplace(boundary.label, 0);
	for (const auto& element : assembler.elements()) {
		if (!assembler.OwnsElementByMinimumNode(element)) continue;
		for (const auto label : element.boundary_labels) {
			auto found = counts.find(label);
			if (found != counts.end()) ++found->second;
		}
	}
	for (auto& entry : counts) {
		long long global = 0;
		MPI_Allreduce(&entry.second, &global, 1, MPI_LONG_LONG, MPI_SUM, communicator);
		entry.second = global;
	}
	return counts;
}

inline void RequireConfiguredScalarSurfaceFaces(
	const SimulationConfiguration& configuration, const CompiledLinearSystem& system,
	const OwnedRowAssembler& assembler, MPI_Comm communicator)
{
	for (const auto& entry : CountConfiguredScalarSurfaceFaces(
		configuration, system, assembler, communicator))
		if (entry.second == 0)
			throw std::runtime_error("configured scalar surface boundary label "
				+std::to_string(entry.first)
				+" has no boundary faces in the .ntiga database; repack with iga_pack");
}

} // namespace iga

#endif

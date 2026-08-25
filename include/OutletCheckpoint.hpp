#ifndef IGA_OUTLET_CHECKPOINT_HPP
#define IGA_OUTLET_CHECKPOINT_HPP

#include "FlowCheckpoint.hpp"
#include "OutletModel.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

inline const char* OutletKindName(FieldBoundaryKind kind)
{
	if (kind == FieldBoundaryKind::Resistance) return "resistance";
	if (kind == FieldBoundaryKind::WindkesselRC) return "windkessel_rc";
	if (kind == FieldBoundaryKind::WindkesselRCR) return "windkessel_rcr";
	throw std::runtime_error("unsupported outlet model kind");
}

inline void RestoreOutletCheckpoint(const FlowCheckpointMetadata& metadata,
	std::vector<OutletModelState>& models)
{
	if (metadata.outlets.size() != models.size())
		throw std::runtime_error("flow checkpoint outlet count does not match configuration");
	for (auto& model : models) {
		const auto found = std::find_if(metadata.outlets.begin(), metadata.outlets.end(),
			[&](const FlowCheckpointOutlet& outlet) { return outlet.label == model.label; });
		if (found == metadata.outlets.end() || found->kind != OutletKindName(model.kind))
			throw std::runtime_error(
				"flow checkpoint outlet labels or kinds do not match configuration");
		model.flow = found->flow;
		model.pressure = found->pressure;
		model.capacitor_pressure = found->capacitor_pressure;
	}
}

inline void AppendOutletCheckpoint(const std::vector<OutletModelState>& models,
	FlowCheckpointMetadata& metadata)
{
	metadata.outlets.clear();
	for (const auto& model : models)
		metadata.outlets.push_back({model.label, OutletKindName(model.kind),
			model.flow, model.pressure, model.capacitor_pressure});
}

} // namespace iga

#endif

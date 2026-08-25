#ifndef IGA_OUTLET_FLOW_HPP
#define IGA_OUTLET_FLOW_HPP

#include "BoundaryFlow.hpp"
#include "OutletModel.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace iga {

inline std::vector<double> IntegrateOutletModelFlows(
	const std::vector<OutletModelState>& models,
	const std::vector<Element>& elements, const std::vector<double>& state)
{
	if (state.size()%4 != 0)
		throw std::runtime_error("flow state must contain four values per node");
	std::vector<double> flows(models.size(), 0.0);
	std::unordered_map<int, std::size_t> model_index;
	for (std::size_t i = 0; i < models.size(); ++i)
		if (!model_index.emplace(models[i].label, i).second)
			throw std::runtime_error("outlet model labels must be unique");
	for (const auto& element : elements) {
		std::vector<std::array<double, 4>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
			if (element.connectivity[a] < 0)
				throw std::runtime_error("outlet element connectivity is negative");
			const auto node = static_cast<std::size_t>(element.connectivity[a]);
			if (node >= state.size()/4)
				throw std::runtime_error("outlet element connectivity exceeds flow state");
			for (int field = 0; field < 4; ++field)
				nodal[a][field] = state[4*node+static_cast<std::size_t>(field)];
		}
		for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
			const auto found = model_index.find(element.boundary_labels[face]);
			if (found != model_index.end())
				flows[found->second] += IntegrateBoundaryFlow(element, face, nodal);
		}
	}
	return flows;
}

} // namespace iga

#endif

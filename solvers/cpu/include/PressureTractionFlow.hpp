#ifndef IGA_PRESSURE_TRACTION_FLOW_HPP
#define IGA_PRESSURE_TRACTION_FLOW_HPP

#include "BoundaryFlow.hpp"
#include "IgaDatabase.hpp"
#include "PressureTraction.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace iga {

inline std::vector<Element> LoadPressureTractionElements(
	Database& database, const std::map<int, double>& tractions)
{
	if (tractions.empty()) return {};
	std::map<int, std::size_t> face_counts;
	for (const auto& traction : tractions) face_counts.emplace(traction.first, 0);
	std::vector<Element> result;
	for (std::uint64_t index = 0; index < database.header().elements; ++index) {
		auto element = database.Load(index);
		bool retained = false;
		for (const auto label : element.boundary_labels) {
			const auto found = face_counts.find(label);
			if (found == face_counts.end()) continue;
			++found->second;
			retained = true;
		}
		if (retained) result.push_back(std::move(element));
	}
	for (const auto& count : face_counts)
		if (count.second == 0)
			throw std::runtime_error("pressure traction label "+std::to_string(count.first)
				+" has no boundary faces in the .ntiga database; repack with iga_pack");
	return result;
}

inline std::vector<double> IntegratePressureTractionForces(
	const std::map<int, double>& tractions, const std::vector<Element>& elements,
	std::size_t nodes)
{
	std::vector<double> result(4*nodes, 0.0);
	for (const auto& element : elements)
		for (std::size_t face = 0; face < element.boundary_labels.size(); ++face) {
			const auto traction = tractions.find(element.boundary_labels[face]);
			if (traction == tractions.end()) continue;
			const auto local = IntegrateBoundaryPressureTraction(
				element, face, traction->second);
			for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
				if (element.connectivity[a] < 0
					|| static_cast<std::size_t>(element.connectivity[a]) >= nodes)
					throw std::runtime_error(
						"pressure traction connectivity exceeds node count");
				const auto node = static_cast<std::size_t>(element.connectivity[a]);
				for (int component = 0; component < 3; ++component)
					result[4*node+static_cast<std::size_t>(component)]
						+= local[4*a+static_cast<std::size_t>(component)];
			}
		}
	return result;
}

} // namespace iga

#endif

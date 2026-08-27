#include "IgaDatabase.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>

int main(int argc, char** argv)
{
	try {
		if (argc != 2) {
			std::cerr << "usage: iga_inspect DATABASE.ntiga\n";
			return 2;
		}
		iga::Database db(argv[1]);
		std::vector<std::uint64_t> counts(db.header().ranks, 0);
		for (auto owner : db.owners()) ++counts.at(static_cast<std::size_t>(owner));
		std::uint64_t total_basis = 0;
		std::size_t max_basis = 0;
		std::uint64_t total_required = 0;
		std::uint64_t min_required = std::numeric_limits<std::uint64_t>::max();
		std::uint64_t max_required = 0;
		std::map<std::int32_t, std::uint64_t> boundary_faces;
		for (std::uint32_t rank = 0; rank < db.header().ranks; ++rank) {
			const auto n = db.RequiredElementIndices(static_cast<std::int32_t>(rank)).size();
			total_required += n;
			min_required = std::min(min_required, static_cast<std::uint64_t>(n));
			max_required = std::max(max_required, static_cast<std::uint64_t>(n));
		}
		for (std::uint64_t i = 0; i < db.header().elements; ++i) {
			auto element = db.Load(i);
			total_basis += element.connectivity.size();
			max_basis = std::max(max_basis, element.connectivity.size());
			for (const auto label : element.boundary_labels)
				if (label >= 0) ++boundary_faces[label];
		}
		auto [minimum, maximum] = std::minmax_element(counts.begin(), counts.end());
		std::cout << "elements: " << db.header().elements << '\n'
			<< "nodes: " << db.header().nodes << '\n'
			<< "ranks: " << db.header().ranks << '\n'
			<< "source origin: [" << db.header().geometry_transform.source_origin[0] << ", "
			<< db.header().geometry_transform.source_origin[1] << ", "
			<< db.header().geometry_transform.source_origin[2] << "]\n"
			<< "source units/normalized unit: "
			<< db.header().geometry_transform.source_units_per_normalized_unit << '\n'
			<< "source length scale to m: "
			<< db.header().geometry_transform.source_length_scale_to_m << '\n'
			<< "basis/element mean: " << static_cast<double>(total_basis) / db.header().elements << '\n'
			<< "basis/element max: " << max_basis << '\n'
			<< "elements/rank min: " << *minimum << '\n'
			<< "elements/rank max: " << *maximum << '\n';
		std::cout << "row-touching elements/rank min: " << min_required << '\n'
			<< "row-touching elements/rank max: " << max_required << '\n'
			<< "row-touching duplication: " << static_cast<double>(total_required) / db.header().elements << "x\n";
		for (const auto& entry : boundary_faces)
			std::cout << "boundary_faces[" << entry.first << "]: " << entry.second << '\n';
		return 0;
	} catch (const std::exception& e) {
		std::cerr << "iga_inspect: " << e.what() << '\n';
		return 1;
	}
}

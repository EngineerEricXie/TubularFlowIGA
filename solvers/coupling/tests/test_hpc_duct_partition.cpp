#include "CouplingFixture.hpp"
#include <iostream>
#include <iterator>

int main(int argc, char** argv)
{
	try {
		if (argc != 2) throw std::runtime_error("expected a temporary output directory");
		const std::filesystem::path root = argv[1];
		std::filesystem::create_directories(root);
		for (const std::uint32_t ranks : {1u, 4u, 8u, 64u}) {
			const auto full = root/("full-"+std::to_string(ranks)+".ntiga");
			const auto compact = root/("compact-"+std::to_string(ranks)+".ntiga");
			iga::test::WriteC2SquareDuctDatabase(full, 4, ranks, 8);
			iga::test::WriteC2SquareDuctDatabase(compact, 4, ranks, 8, 1.0, true);
			iga::Database database(compact.string());
			const auto prefix = database.header().rank_index_offset;
			std::vector<char> before(prefix), after(prefix);
			std::ifstream old_input(full, std::ios::binary), new_input(compact, std::ios::binary);
			old_input.read(before.data(), before.size());
			new_input.read(after.data(), after.size());
			if (!old_input || !new_input || before != after)
				throw std::runtime_error("compact index changed physical records or partition owners");
			std::uint64_t total = 0;
			for (std::uint32_t rank = 0; rank < ranks; ++rank) {
				const auto range = database.NodeRange(rank);
				std::vector<std::uint64_t> expected;
				for (std::uint64_t id = 0; id < database.header().elements; ++id) {
					const auto element = database.Load(id);
					if (std::any_of(element.connectivity.begin(), element.connectivity.end(),
						[&](auto node) { return static_cast<std::uint64_t>(node) >= range.first
							&& static_cast<std::uint64_t>(node) < range.second; })) expected.push_back(id);
				}
				const auto actual = database.RequiredElementIndices(rank);
				if (actual != expected) throw std::runtime_error("required index differs from owned-row oracle");
				total += actual.size();
			}
			if (ranks > 1 && total >= ranks*database.header().elements)
				throw std::runtime_error("compact index still replicates all elements");
			std::cout << "ranks=" << ranks << " elements=" << database.header().elements
				<< " required_total=" << total << " physical_records=identical coverage=passed\n";
		}
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

int main(int argc, char** argv)
{
	try {
		if (argc != 4) throw std::runtime_error("usage: native_graph_checkpoint_fixture_test ROOT zero|flow|species|multi RANKS");
		const fs::path root = argv[1]; const std::string mode = argv[2];
		const auto ranks = std::stoul(argv[3]);
		if (ranks < 1 || ranks > 128 || !fs::create_directories(root)) throw std::runtime_error("invalid ranks or existing fixture root");
		if (mode == "zero") {
			for (const auto* name : {"zero_source", "zero_terminal", "zero_junction"}) fs::create_directory(root/name);
			WriteZeroDModel(root/"zero_source", true); WriteZeroDModel(root/"zero_terminal", false);
			WriteZeroDThreeDCase(root/"zero_junction"); WriteDatabase(root/"zero_junction.ntiga", ranks); WriteZeroDClockGraph(root);
		} else if (mode == "flow") {
			for (const auto* name : {"source", "branch_a", "branch_b", "three_d"}) fs::create_directory(root/name);
			WriteOneDCase(root/"source", 1e-3); WriteOneDCase(root/"branch_a", 6e-4); WriteOneDCase(root/"branch_b", 4e-4);
			WriteThreeDCase(root/"three_d"); WriteDatabase(root/"graph.ntiga", ranks); WriteGraph(root, "graph.ntiga", "explicit");
		} else if (mode == "species") {
			for (const auto* name : {"species_source", "species_leaf_a", "species_leaf_b", "junction"}) fs::create_directory(root/name);
			WriteOneDTransportCase(root/"species_source", 1e-3, "source_red", "source_blue");
			WriteOneDTransportCase(root/"species_leaf_a", 6e-4, "leaf_a_red", "leaf_a_blue");
			WriteOneDTransportCase(root/"species_leaf_b", 4e-4, "leaf_b_red", "leaf_b_blue");
			WriteThreeDTransportCase(root/"junction"); WriteDatabase(root/"group.ntiga", ranks); WriteSpeciesGraph(root, "group.ntiga", true);
		} else if (mode == "multi") {
			for (const auto* name : {"three_a", "three_b", "source", "bridge", "leaf_a", "leaf_b", "leaf_c"}) fs::create_directory(root/name);
			WriteThreeDCase(root/"three_a"); WriteThreeDCase(root/"three_b"); WriteOneDCase(root/"source", 1e-3);
			WriteOneDCase(root/"bridge", 5.5e-4); WriteOneDCase(root/"leaf_a", 4.5e-4); WriteOneDCase(root/"leaf_b", 3e-4); WriteOneDCase(root/"leaf_c", 2.5e-4);
			WriteDatabase(root/"a.ntiga", ranks); WriteDatabase(root/"b.ntiga", ranks); WriteMultiIslandGraph(root, "explicit");
		} else throw std::runtime_error("unknown fixture mode");
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

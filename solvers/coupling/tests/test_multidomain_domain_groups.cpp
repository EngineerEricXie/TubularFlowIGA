#include "MultidomainRunner.hpp"

#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

namespace {

void Prepare(const fs::path& root, std::uint32_t database_ranks, bool grouped)
{
	for (const auto& domain : {"three_a", "three_b", "source", "bridge",
		"leaf_a", "leaf_b", "leaf_c"}) fs::create_directories(root/domain);
	WriteThreeDCase(root/"three_a");
	WriteThreeDCase(root/"three_b");
	WriteOneDCase(root/"source", 1.0e-3);
	WriteOneDCase(root/"bridge", 5.5e-4);
	WriteOneDCase(root/"leaf_a", 4.5e-4);
	WriteOneDCase(root/"leaf_b", 3.0e-4);
	WriteOneDCase(root/"leaf_c", 2.5e-4);
	WriteDatabase(root/"a.ntiga", database_ranks);
	WriteDatabase(root/"b.ntiga", database_ranks);
	WriteMultiIslandGraph(root, "explicit");
	if (!grouped) return;
	auto config = ReadFile(root/"simulation_config.json");
	const auto position = config.find("  \"domains\":[");
	if (position == std::string::npos) throw std::runtime_error("domain catalog marker missing");
	config.insert(position,
		"  \"resources\":{\"mode\":\"domain_groups\",\"groups\":["
		"{\"id\":\"small\",\"ranks\":1,\"domains\":[\"source\",\"bridge\",\"leaf_a\",\"leaf_b\",\"leaf_c\"]},"
		"{\"id\":\"island_a\",\"ranks\":2,\"domains\":[\"island_a\"]},"
		"{\"id\":\"island_b\",\"ranks\":2,\"domains\":[\"island_b\"]}]},\n");
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << config;
}

void PrepareSpecies(const fs::path& root, std::uint32_t database_ranks, bool grouped)
{
	for (const auto& domain : {"junction", "species_source", "species_leaf_a",
		"species_leaf_b"}) fs::create_directories(root/domain);
	WriteThreeDTransportCase(root/"junction");
	WriteOneDTransportCase(root/"species_source", 1.0e-3, "source_red", "source_blue");
	WriteOneDTransportCase(root/"species_leaf_a", 6.0e-4, "leaf_a_red", "leaf_a_blue");
	WriteOneDTransportCase(root/"species_leaf_b", 4.0e-4, "leaf_b_red", "leaf_b_blue");
	WriteDatabase(root/"junction.ntiga", database_ranks);
	WriteSpeciesGraph(root, "junction.ntiga");
	if (!grouped) return;
	auto config = ReadFile(root/"simulation_config.json");
	const auto position = config.find("  \"domains\":[");
	if (position == std::string::npos) throw std::runtime_error("species domain catalog marker missing");
	config.insert(position,
		"  \"resources\":{\"mode\":\"domain_groups\",\"groups\":["
		"{\"id\":\"small\",\"ranks\":1,\"domains\":[\"source\",\"leaf_a\",\"leaf_b\"]},"
		"{\"id\":\"junction\",\"ranks\":2,\"domains\":[\"junction\"]}]},\n");
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << config;
}

void PrepareZeroD(const fs::path& root, std::uint32_t database_ranks, bool grouped)
{
	for (const auto& domain : {"zero_source", "zero_terminal", "zero_junction"})
		fs::create_directories(root/domain);
	WriteZeroDModel(root/"zero_source", true);
	WriteZeroDModel(root/"zero_terminal", false);
	WriteZeroDThreeDCase(root/"zero_junction");
	WriteDatabase(root/"zero_junction.ntiga", database_ranks);
	WriteZeroDClockGraph(root);
	if (!grouped) return;
	auto config = ReadFile(root/"simulation_config.json");
	const auto position = config.find("  \"domains\":[");
	if (position == std::string::npos) throw std::runtime_error("zero-D domain catalog marker missing");
	config.insert(position,
		"  \"resources\":{\"mode\":\"domain_groups\",\"groups\":["
		"{\"id\":\"zero_d\",\"ranks\":1,\"domains\":[\"source_0d\",\"terminal_0d\"]},"
		"{\"id\":\"junction\",\"ranks\":2,\"domains\":[\"junction\"]}]},\n");
	std::ofstream output(root/"simulation_config.json", std::ios::trunc);
	output << config;
}

void Run(const fs::path& root, const fs::path& output, MPI_Comm communicator)
{
	std::vector<std::string> arguments{"domain_group_test", "--graph-case", root.string(),
		"--output-dir", output.string()};
	std::vector<char*> argv;
	for (auto& argument : arguments) argv.push_back(argument.data());
	argv.push_back(nullptr);
	if (iga::RunMultidomainFlow(static_cast<int>(arguments.size()), argv.data(), communicator))
		throw std::runtime_error("multidomain domain-group run failed");
}

void ValidateResourceManifest(const fs::path& output)
{
	const auto manifest = ReadFile(output/"graph_binding_manifest.json");
	if (manifest.find("\"resource_mode\": \"domain_groups\"") == std::string::npos)
		throw std::runtime_error("completion manifest lacks domain-group resource mode");
	if (manifest.find("\"resource_reconstruction\": \"contiguous_manifest_group_order\"")
		== std::string::npos)
		throw std::runtime_error("completion manifest lacks resource reconstruction rule");
}

double Compare(const fs::path& reference, const fs::path& actual)
{
	double maximum = 0.0;
	for (const auto& table : {"pressure_flow_edges.csv", "pressure_flow_ports.csv",
		"pressure_flow_domain_balances.csv", "zero_d_flow_history.csv"}) {
		const auto first = ReadCsv(reference/table);
		const auto second = ReadCsv(actual/table);
		if (first.size() != second.size()) throw std::runtime_error("grouped row count differs");
		for (std::size_t row = 0; row < first.size(); ++row) {
			if (first[row].size() != second[row].size())
				throw std::runtime_error("grouped column count differs");
			for (const auto& item : first[row]) {
				const auto found = second[row].find(item.first);
				if (found == second[row].end()) throw std::runtime_error("grouped column missing");
				char* end = nullptr;
				const double expected = std::strtod(item.second.c_str(), &end);
				if (!item.second.empty() && end == item.second.c_str()+item.second.size()) {
					const double value = Value(second[row], item.first);
					const double scale = std::max({1.0, std::abs(expected), std::abs(value)});
					maximum = std::max(maximum, std::abs(value-expected)/scale);
				} else if (item.second != found->second)
					throw std::runtime_error("grouped identity differs");
			}
		}
	}
	if (maximum > 1.0e-6) throw std::runtime_error("grouped numerical result differs");
	return maximum;
}

double CompareSpecies(const fs::path& reference, const fs::path& actual)
{
	double maximum = 0.0;
	for (const auto& table : {"species_hydraulic_steps.csv", "species_edge_amounts.csv",
		"species_domain_accounting.csv", "species_global_balance.csv",
		"species_logical_ports.csv"}) {
		const auto first = ReadCsv(reference/table);
		const auto second = ReadCsv(actual/table);
		if (first.size() != second.size()) throw std::runtime_error("grouped species row count differs");
		for (std::size_t row = 0; row < first.size(); ++row) for (const auto& item : first[row]) {
			const auto found = second[row].find(item.first);
			if (found == second[row].end()) throw std::runtime_error("grouped species column missing");
			char* end = nullptr;
			const double expected = std::strtod(item.second.c_str(), &end);
			if (!item.second.empty() && end == item.second.c_str()+item.second.size()) {
				const double value = Value(second[row], item.first);
				maximum = std::max(maximum, std::abs(value-expected)
					/std::max({1.0, std::abs(value), std::abs(expected)}));
			} else if (item.second != found->second)
				throw std::runtime_error("grouped species identity differs");
		}
	}
	if (maximum > 1.0e-6) throw std::runtime_error("grouped species result differs");
	return maximum;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		if (argc != 2 || ranks != 5)
			throw std::runtime_error("usage: mpiexec -np 5 multidomain_domain_groups_test OUTPUT");
		const fs::path root(argv[1]);
		if (rank == 0) {
			fs::create_directories(root/"serial");
			fs::create_directories(root/"grouped");
			Prepare(root/"serial", 1, false);
			Prepare(root/"grouped", 2, true);
			Run(root/"serial", root/"serial-output", PETSC_COMM_SELF);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		Run(root/"grouped", root/"grouped-output", PETSC_COMM_WORLD);
		if (rank == 0) {
			ValidateResourceManifest(root/"grouped-output");
			std::cout << "multidomain domain groups maximum_relative_error="
				<< Compare(root/"serial-output", root/"grouped-output") << " passed\n";
			fs::create_directories(root/"species-serial");
			fs::create_directories(root/"species-grouped");
			PrepareSpecies(root/"species-serial", 1, false);
			PrepareSpecies(root/"species-grouped", 2, true);
			Run(root/"species-serial", root/"species-serial-output", PETSC_COMM_SELF);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		Run(root/"species-grouped", root/"species-grouped-output", PETSC_COMM_WORLD);
		if (rank == 0) std::cout << "multidomain staged domain groups maximum_relative_error="
			<< CompareSpecies(root/"species-serial-output", root/"species-grouped-output")
			<< " passed\n";
		if (rank == 0) {
			fs::create_directories(root/"zero-serial");
			fs::create_directories(root/"zero-grouped");
			PrepareZeroD(root/"zero-serial", 1, false);
			PrepareZeroD(root/"zero-grouped", 2, true);
			Run(root/"zero-serial", root/"zero-serial-output", PETSC_COMM_SELF);
		}
		MPI_Barrier(PETSC_COMM_WORLD);
		Run(root/"zero-grouped", root/"zero-grouped-output", PETSC_COMM_WORLD);
		if (rank == 0) std::cout << "multidomain zero-D owner groups maximum_relative_error="
			<< Compare(root/"zero-serial-output", root/"zero-grouped-output") << " passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
	}
	PetscFinalize();
	return 0;
}

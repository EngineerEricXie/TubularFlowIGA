#include "MultidomainRunner.hpp"

// Reuse the established physical fixture and its conservation gates.
#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

namespace {

void RunGraph(const fs::path& root, const fs::path& output, MPI_Comm communicator)
{
	std::vector<std::string> arguments{"multidomain_subcommunicator_test",
		"--graph-case", root.string(), "--output-dir", output.string()};
	std::vector<char*> argv;
	for (auto& argument : arguments) argv.push_back(argument.data());
	argv.push_back(nullptr);
	if (iga::RunMultidomainFlow(static_cast<int>(arguments.size()), argv.data(), communicator))
		throw std::runtime_error("graph failed on supplied communicator");
	if (!fs::exists(output/"graph_binding_manifest.json"))
		throw std::runtime_error("successful graph did not publish completion manifest");
}

double CompareColumn(const std::vector<CsvRow>& reference,
	const std::vector<CsvRow>& current, const std::string& name)
{
	double norm = 0.0, error = 0.0;
	for (std::size_t row = 0; row < reference.size(); ++row) {
		const double expected = Value(reference[row], name);
		const double actual = Value(current[row], name);
		if (!std::isfinite(expected) || !std::isfinite(actual))
			throw std::runtime_error("nonfinite graph output: "+name);
		norm = std::hypot(norm, expected);
		error = std::hypot(error, actual-expected);
	}
	if (norm == 0.0 ? error > 1e-12 : error/norm > 1e-6)
		throw std::runtime_error("graph differs from serial reference: "+name);
	return norm > 0.0 ? error/norm : error;
}

double CompareGraph(const fs::path& reference, const fs::path& current)
{
	Validate(reference, true);
	Validate(current, true);
	double maximum = 0.0;
	const auto one = ReadCsv(reference/"pressure_flow_edges.csv");
	const auto many = ReadCsv(current/"pressure_flow_edges.csv");
	if (one.size() != many.size()) throw std::runtime_error("graph edge count differs");
	for (std::size_t row = 0; row < one.size(); ++row)
		for (const auto& key : {"step", "time_s", "edge_id"})
			if (one[row].at(key) != many[row].at(key))
				throw std::runtime_error("graph edge identity differs");
	for (const auto& name : {"applied_pressure_pa", "measured_pressure_pa",
		"first_outward_flow_m3_s", "second_outward_flow_m3_s"})
		maximum = std::max(maximum, CompareColumn(one, many, name));
	const auto steps = ReadCsv(reference/"pressure_flow_steps.csv");
	const auto distributed_steps = ReadCsv(current/"pressure_flow_steps.csv");
	if (steps.size() != distributed_steps.size()) throw std::runtime_error("graph step count differs");
	for (std::size_t row = 0; row < steps.size(); ++row)
		if (steps[row].at("iterations") != distributed_steps[row].at("iterations"))
			throw std::runtime_error("graph acceptance iteration differs");
	return maximum;
}

double CompareSpeciesGraph(const fs::path& reference, const fs::path& current)
{
	ValidateSpeciesMultidomain(reference);
	ValidateSpeciesMultidomain(current);
	double maximum = 0.0;
	const std::map<std::string, std::vector<std::string>> tables{
		{"species_edge_amounts.csv", {"first_outward_amount", "second_outward_amount"}},
		{"species_domain_accounting.csv", {"M0_amount", "M1_amount", "source_amount"}},
		{"species_global_balance.csv", {"M0_amount", "M1_amount", "source_amount"}}};
	for (const auto& table : tables) {
		const auto one = ReadCsv(reference/table.first);
		const auto many = ReadCsv(current/table.first);
		if (one.size() != many.size()) throw std::runtime_error("species row count differs");
		std::map<std::string, std::vector<CsvRow>> expected, actual;
		for (std::size_t row = 0; row < one.size(); ++row) {
			for (const auto& key : {"step", "time_s", "edge_id", "domain_id", "species_id",
				"first_domain_id", "first_port_id", "second_domain_id", "second_port_id", "donor"})
				if (one[row].count(key) && one[row].at(key) != many[row].at(key))
					throw std::runtime_error("species identity or routing differs");
			const auto& species = one[row].at("species_id");
			expected[species].push_back(one[row]);
			actual[species].push_back(many[row]);
		}
		// Each norm contains only one species and one physical quantity.
		for (const auto& species : expected) {
			for (const auto& column : table.second)
				maximum = std::max(maximum, CompareColumn(species.second, actual.at(species.first), column));
			// In this source-free, constant-concentration fixture, each net
			// amount is analytically zero. Check BOTH solves against that zero
			// with the declared 1e-12 absolute gate; dividing one cancellation
			// roundoff by another is not a relative field comparison.
			if (table.first != "species_edge_amounts.csv") {
				const std::string net = table.first == "species_domain_accounting.csv"
					? "total_outward_amount" : "outward_amount";
				auto zero = species.second;
				for (auto& row : zero) row[net] = "0";
				CompareColumn(zero, species.second, net);
				CompareColumn(zero, actual.at(species.first), net);
			}
		}
	}
	return maximum;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		if ((argc != 2 && argc != 3) || ranks != 3
			|| (argc == 3 && std::string(argv[2]) != "species"))
			throw std::runtime_error("usage: mpiexec -np 3 multidomain_subcommunicator_test OUTPUT_PARENT [species]");
		const bool species = argc == 3;
		const int color = rank == 0 ? 0 : 1;
		MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
		int local_rank = 0, local_size = 1;
		MPI_Comm_rank(group, &local_rank);
		MPI_Comm_size(group, &local_size);
		// Deliberately distinct solver-option databases across disjoint groups.
		// PREONLY reference solves preserve the existing physical comparison.
		PetscOptionsSetValue(nullptr, "-ksp_rtol", color == 0 ? "1e-9" : "1e-10");
		const auto root = fs::path(argv[1])/("graph-group-"+std::to_string(color));
		// Different coupling schemes execute different collective sequences.
		// A stray world collective cannot be satisfied by the other graph.
		const std::string scheme = species ? "species_explicit" : color == 0 ? "explicit" : "aitken";
		const auto write_graph = [&](const std::string& database) {
			if (species) WriteSpeciesGraph(root, database, color == 1);
			else WriteGraph(root, database, scheme);
		};
		if (local_rank == 0) {
			if (!fs::create_directories(root)) throw std::runtime_error("fixture directory already exists");
			for (const auto& domain : {"three_d", "source", "branch_a", "branch_b",
				"junction", "species_source", "species_leaf_a", "species_leaf_b"})
				fs::create_directory(root/domain);
			if (species) {
				WriteThreeDTransportCase(root/"junction");
				WriteOneDTransportCase(root/"species_source", 1e-3, "source_red", "source_blue");
				WriteOneDTransportCase(root/"species_leaf_a", 6e-4, "leaf_a_red", "leaf_a_blue");
				WriteOneDTransportCase(root/"species_leaf_b", 4e-4, "leaf_b_red", "leaf_b_blue");
			} else {
				WriteThreeDCase(root/"three_d");
				WriteOneDCase(root/"source", 1e-3);
				WriteOneDCase(root/"branch_a", 6e-4);
				WriteOneDCase(root/"branch_b", 4e-4);
			}
			WriteDatabase(root/"serial.ntiga", 1);
			write_graph("serial.ntiga");
			RunGraph(root, root/"serial", PETSC_COMM_SELF);
			WriteDatabase(root/"group.ntiga", static_cast<std::uint32_t>(local_size));
			write_graph("group.ntiga");
		}
		MPI_Barrier(group);
		RunGraph(root, root/"distributed", group);
		if (local_rank == 0) {
			const double error = species ? CompareSpeciesGraph(root/"serial", root/"distributed")
				: CompareGraph(root/"serial", root/"distributed");
			std::cout << "multidomain_subcommunicator group=" << color
				<< " ranks=" << local_size << " scheme=" << scheme
				<< " maximum_relative_l2=" << error << " passed\n";
		}
		// The synchronous runner destroys every runtime before returning.
		MPI_Comm_free(&group);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

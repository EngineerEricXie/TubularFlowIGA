#include "CollectiveFailure.hpp"
#include "CollectivePetscOptions.hpp"
#include "MultidomainRunner.hpp"
#include "MultidomainConfig.hpp"
#include "OneDConfig.hpp"
#include "SimulationConfig.hpp"
#include "TemporalFunction.hpp"
#include "CaseInput.hpp"
#include "BoundarySupport.hpp"
#include "OneDNetwork.hpp"
#include <sys/stat.h>

#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1, status = 2;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		if (argc != 3 || ranks != 3)
			throw std::runtime_error("usage: mpiexec -np 3 multidomain_failure_test OUTPUT_PARENT MODE");
		const fs::path root(argv[1]);
		const std::string mode(argv[2]);
		const std::set<std::string> modes{"arguments", "input", "database", "output", "healthy",
			"stop", "newton", "injection", "manifest", "one-d-config", "three-d-config",
			"replica", "zero-d-config", "zero-d-healthy", "petsc-ksp", "petsc-prefix",
			"petsc-flag", "petsc-used", "petsc-file", "petsc-file-different",
			"asset-database", "asset-network", "asset-mesh", "asset-velocity",
			"asset-waveform-1d", "asset-waveform-3d", "asset-missing", "asset-fifo", "asset-replica"};
		if (!modes.count(mode))
			throw std::runtime_error("unknown test mode");
		const bool zero_d = mode == "zero-d-config" || mode == "zero-d-healthy";
		const bool asset_mode = mode.rfind("asset-", 0) == 0;
		const bool healthy = mode == "healthy" || mode == "replica" || mode == "zero-d-healthy" || mode == "asset-replica"
			|| mode == "petsc-used" || mode == "petsc-file";
		const bool replica = mode == "database" || mode == "manifest" || mode == "one-d-config"
			|| mode == "three-d-config" || mode == "replica" || mode == "zero-d-config" || asset_mode;
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "failure fixture", [&] {
			if (rank != 0) return;
			if (!fs::create_directories(root)) throw std::runtime_error("test directory already exists");
			for (const auto& domain : {"three_d", "source", "branch_a", "branch_b"})
				fs::create_directory(root/domain);
			WriteThreeDCase(root/"three_d");
			WriteOneDCase(root/"source", 1e-3);
			WriteOneDCase(root/"branch_a", 6e-4);
			WriteOneDCase(root/"branch_b", 4e-4);
			WriteDatabase(root/"graph.ntiga", 3);
			WriteGraph(root, "graph.ntiga", "explicit");
			const auto replace = [&](const fs::path& path, const std::string& before,
				const std::string& after) {
				auto text = ReadFile(path);
				const auto position = text.find(before);
				if (position == std::string::npos) throw std::runtime_error("missing mutation target");
				text.replace(position, before.size(), after);
				std::ofstream(path, std::ios::trunc) << text;
			};
			if (asset_mode) {
				replace(root/"source/simulation_config.json",
					"\"kind\":\"constant\",\"units\":\"m3/s\",\"value\":"+Number(1e-3),
					"\"kind\":\"periodic_table\",\"units\":\"m3/s\",\"period\":1,\"file\":\"inlet.csv\",\"interpolation\":\"linear\"");
				replace(root/"three_d/simulation_config.json", "\"fields\":[",
					"\"temporal_functions\":[{\"name\":\"inlet_scale\",\"kind\":\"periodic_table\","
					"\"units\":\"dimensionless\",\"period\":1,\"file\":\"scale.csv\",\"interpolation\":\"linear\"}],\"fields\":[");
				replace(root/"three_d/simulation_config.json", "\"scale\":1}",
					"\"scale\":1,\"waveform\":\"inlet_scale\"}");
				std::ofstream(root/"source/inlet.csv") << "time,value\n0,0.001\n0.5,0.001\n";
				std::ofstream(root/"three_d/scale.csv") << "time,value\n0,1\n0.5,1\n";
				(void)iga::ParseOneDConfiguration(ReadFile(root/"source/simulation_config.json"));
				(void)iga::ReadSimulationConfiguration((root/"three_d/simulation_config.json").string());
			}
			if (zero_d) {
				for (const auto& domain : {"zero_source", "zero_terminal", "zero_junction"})
					fs::create_directory(root/domain);
				WriteZeroDModel(root/"zero_source", true);
				WriteZeroDModel(root/"zero_terminal", false);
				WriteZeroDThreeDCase(root/"zero_junction");
				WriteDatabase(root/"zero_junction.ntiga", 3);
				WriteZeroDClockGraph(root);
			}
			if (replica) {
				const auto bad = root/"bad";
				// List before making a child, so recursive copying cannot include itself.
				std::vector<fs::path> entries;
				for (const auto& entry : fs::directory_iterator(root)) entries.push_back(entry.path());
				fs::create_directory(bad);
				for (const auto& entry : entries)
					fs::copy(entry, bad/entry.filename(), fs::copy_options::recursive);
				if (mode == "database") std::ofstream(bad/"graph.ntiga") << "truncated";
				if (mode == "manifest") WriteGraph(bad, "graph.ntiga", "aitken");
				if (mode == "one-d-config") replace(bad/"source/simulation_config.json",
					"\"cells_per_segment\":1", "\"cells_per_segment\":2");
				if (mode == "three-d-config") replace(bad/"three_d/simulation_config.json",
					"\"dt\":"+Number(kDt), "\"dt\":"+Number(kDt/2));
				if (mode == "zero-d-config") replace(bad/"zero_source/zero_d_model.json",
					"\"initial_pressure_pa\":0.1", "\"initial_pressure_pa\":0.2");
				if (mode == "asset-database") {
					std::fstream database(bad/"graph.ntiga", std::ios::in | std::ios::out | std::ios::binary);
					database.seekp(72); iga::Write(database, 1.125); database.close();
					iga::Database verified((bad/"graph.ntiga").string());
					if (verified.header().nodes != 64 || verified.header().elements != 1
						|| verified.header().ranks != 3 || verified.Load(0).connectivity.size() != 64)
						throw std::runtime_error("mutated database failed standalone validation");
				}
				if (mode == "asset-network") {
					replace(bad/"source/tree.swc", "2 2 1 0 0 ", "2 2 1 0.1 0 ");
					(void)iga::ReadOneDNetwork(bad/"source/tree.swc", 1.0, 1, kViscosity);
				}
				if (mode == "asset-mesh") {
					replace(bad/"three_d/controlmesh.vtk", "POINTS 64 double\n0 0 0\n", "POINTS 64 double\n0.001 0 0\n");
					(void)iga::ReadLabeledHexMesh((bad/"three_d/controlmesh.vtk").string(), 64, 1);
				}
				if (mode == "asset-velocity") {
					replace(bad/"three_d/initial_velocityfield.txt", "0 0 0\n", "1 0 0\n");
					(void)iga::ReadVelocity((bad/"three_d/initial_velocityfield.txt").string(), 64);
				}
				if (mode == "asset-waveform-1d") {
					replace(bad/"source/inlet.csv", "0.5,0.001", "0.5,0.002");
					(void)iga::ReadTemporalCsv((bad/"source/inlet.csv").string(), 1.0);
				}
				if (mode == "asset-waveform-3d") {
					replace(bad/"three_d/scale.csv", "0.5,1", "0.5,2");
					(void)iga::ReadTemporalCsv((bad/"three_d/scale.csv").string(), 1.0);
				}
				if (mode == "asset-missing") fs::remove(bad/"source/inlet.csv");
				if (mode == "asset-fifo") {
					fs::remove(bad/"source/inlet.csv");
					if (::mkfifo((bad/"source/inlet.csv").c_str(), 0600) != 0)
						throw std::runtime_error("cannot create asset FIFO fixture");
				}
				// Both inputs must parse independently; otherwise this is a local
				// validation failure, not a test of differing accepted settings.
				(void)iga::ReadMultidomainConfiguration((bad/"simulation_config.json").string());
				if (mode == "one-d-config")
					(void)iga::ParseOneDConfiguration(ReadFile(bad/"source/simulation_config.json"));
				if (mode == "three-d-config")
					(void)iga::ReadSimulationConfiguration((bad/"three_d/simulation_config.json").string());
			}
			if (mode == "output") std::ofstream(root/"blocked-parent") << "regular file";
			if (mode == "petsc-file" || mode == "petsc-file-different")
				for (int peer = 0; peer < ranks; ++peer) {
					std::ofstream file(root/("rank "+std::to_string(peer)+" options.txt"));
					file << "-ksp_type preonly\n-pc_type lu\n-pc_factor_mat_solver_type mumps\n-ksp_rtol "
						<< (mode == "petsc-file-different" && peer == 1 ? "1e-2" : "1e-10") << '\n';
				}
		});
		const auto graph = rank == 1 && mode == "input" ? root/"missing"
			: rank == 1 && replica ? root/"bad" : root;
		const auto output = mode == "output" ? root/"blocked-parent"/"results" : root/"results";
		std::vector<std::string> arguments{"multidomain_failure_test", "--graph-case",
			graph.string(), "--output-dir", output.string()};
		if (rank == 1 && mode == "arguments") arguments.push_back("--invalid-local-option");
		if (rank == 1 && mode == "stop") arguments.insert(arguments.end(), {"--stop-after-step", "1"});
		if (rank == 1 && mode == "newton") arguments.insert(arguments.end(), {"--three-d-max-newton", "1"});
		if (rank == 1 && mode == "injection")
			setenv("TUBULARFLOWIGA_INJECT_EXPLICIT_COUPLING_FAILURE_STEP", "1", 1);
		if (rank == 1 && mode == "petsc-ksp") PetscOptionsSetValue(nullptr, "-ksp_type", "cg");
		if (rank == 1 && mode == "petsc-prefix") PetscOptionsSetValue(nullptr, "-fieldsplit_0_pc_type", "jacobi");
		if (rank == 1 && mode == "petsc-flag") PetscOptionsSetValue(nullptr, "-ksp_monitor", nullptr);
		if (mode == "petsc-used") {
			PetscOptionsSetValue(nullptr, "-hpc_note", "quoted \"value\" -fake_option\nline");
			if (rank == 1) {
				char text[128]{};
				PetscOptionsGetString(nullptr, nullptr, "-hpc_note", text, sizeof(text), nullptr);
			}
		}
		if (mode == "petsc-file" || mode == "petsc-file-different") {
			const auto path = root/("rank "+std::to_string(rank)+" options.txt");
			// Each file is loaded on SELF, outside a group-local-stage callback.
			iga::petsc_options_detail::Check(PetscOptionsInsertFile(PETSC_COMM_SELF, nullptr,
				path.string().c_str(), PETSC_TRUE), "PetscOptionsInsertFile");
			PetscOptionsSetValue(nullptr, "-options_file", path.string().c_str());
		}
		if (mode == "replica") {
			// Also populate PETSc's copy of the application arguments, as native
			// PetscInitialize does. Its rank-local path must not force rejection.
			PetscOptionsSetValue(nullptr, "--graph-case", graph.string().c_str());
			PetscOptionsSetValue(nullptr, "--output-dir", output.string().c_str());
		}
		std::vector<char*> raw;
		for (auto& argument : arguments) raw.push_back(argument.data());
		raw.push_back(nullptr);
		status = iga::RunMultidomainFlow(static_cast<int>(arguments.size()), raw.data(), PETSC_COMM_WORLD);
		int minimum = 0, maximum = 0;
		MPI_Allreduce(&status, &minimum, 1, MPI_INT, MPI_MIN, PETSC_COMM_WORLD);
		MPI_Allreduce(&status, &maximum, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
		if (minimum != maximum || status != (healthy ? 0 : 1))
			throw std::runtime_error("runner exit status was not the expected common result");
		iga::CollectiveLocalStage(PETSC_COMM_WORLD, "output verification", [&] {
			if (mode == "petsc-used") {
				iga::petsc_options_detail::Buffers buffers(nullptr);
				iga::petsc_options_detail::Check(PetscOptionsLeftGet(nullptr, &buffers.count,
					&buffers.names, &buffers.values), "PetscOptionsLeftGet");
				buffers.have_unused = true;
				bool unused_note = false;
				for (PetscInt i = 0; i < buffers.count; ++i)
					unused_note = unused_note || std::string(buffers.names[i]) == "hpc_note";
				if (unused_note != (rank != 1))
					throw std::runtime_error("graph option agreement changed unused-option diagnostics");
			}
			if (healthy) {
				if (zero_d) ValidateZeroDClockRun(output);
				else Validate(output, true);
			}
			else if (fs::exists(output/"graph_binding_manifest.json"))
				throw std::runtime_error("failed graph published a success marker");
		});
		std::cout << "failure test mode=" << mode << " rank=" << rank << " runner_status="
			<< status << " verified\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 2);
	}
	PetscFinalize();
	// Fault modes deliberately return the runner's nonzero status to mpiexec.
	return status;
}

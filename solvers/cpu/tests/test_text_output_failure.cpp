#include "PhysiologyOutput.hpp"
#include "VelocitySeries.hpp"
#include "BezierVisualization.hpp"
#include "CouplingHistory.hpp"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <csignal>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace {

std::string Read(const fs::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot read test output");
	return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}

void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

} // namespace

int main(int argc, char** argv)
{
#if defined(__unix__) || defined(__APPLE__)
	try {
		Require(argc == 2, "usage: text_output_failure_test NEW_OUTPUT_DIRECTORY");
		const fs::path directory(argv[1]);
		Require(fs::create_directory(directory), "test directory must be new");
		const auto mesh = directory/"mesh.vtk";
		{
			std::ofstream output(mesh);
			output << "# vtk DataFile Version 3.0\ntest\nASCII\n"
				<< "DATASET UNSTRUCTURED_GRID\n"
				<< "POINTS 4 double\n0 0 0 1 0 0 0 1 0 0 0 1\n"
				<< "CELLS 1 5\n4 0 1 2 3\nCELL_TYPES 1\n10\n";
			output.close();
			Require(bool(output), "cannot write test mesh");
		}
		struct Writer {
			std::string name;
			std::string diagnostic;
			std::function<void(const fs::path&, bool)> write;
			bool supports_large = true;
		};
		const std::vector<Writer> writers{
			{"velocity", "cannot write velocity manifest", [&](const fs::path& path, bool large) {
				std::vector<iga::VelocitySnapshot> snapshots;
				for (int i = 0; i < (large ? 2048 : 5); ++i)
					snapshots.push_back({i*0.25, "field.step000001.txt"});
				iga::WriteVelocityManifest(path, snapshots);
			}},
			{"vtu", "cannot write VTU output", [&](const fs::path& path, bool large) {
				iga::WriteVtu(mesh, path, {{large ? std::string(32768, 'v') : "scalar",
					1, {1, 2, 3, 4}}}, 0.25);
			}},
			{"pvd", "cannot write PVD output", [&](const fs::path& path, bool large) {
				std::vector<std::pair<double, fs::path>> snapshots(
					large ? 1024 : 1, {0.25, "field.vtu"});
				iga::WritePvd(path, snapshots);
			}},
			{"physiology", "cannot write physiology manifest", [&](const fs::path& path, bool large) {
				iga::PhysiologyDefinition physiology;
				physiology.enabled = true;
				physiology.derived_fields = {"pO2"};
				iga::WritePhysiologyManifest(path, physiology,
					{large ? std::string(32768, 's') : "oxygen"});
			}},
			{"geometry", "cannot write Bezier geometry report", [&](const fs::path& path, bool) {
				iga::BezierGeometryValidation validation;
				validation.elements = 1;
				validation.minimum_jacobian = 1.0;
				iga::WriteBezierGeometryReport(path, validation);
			}, false},
			{"history", "cannot write coupling manifest", [&](const fs::path& path, bool large) {
				// The public writer fixes its basename. Redirect only this test's
				// owned fixture link to exercise the same output-path fault cases.
				const auto link = directory/"coupling_manifest.json";
				fs::remove(link);
				fs::create_symlink(path.filename(), link);
				iga::CouplingHistoryWriter writer(directory, iga::SimulationScopeMode::VcaClosedLoop);
				writer.Write(large ? std::string(32768, 'b') : "test");
			}}
		};
		int failures = 0;
		int injected = 0;
		int retries = 0;
		for (const auto& writer : writers) {
			const auto healthy = directory/(writer.name+"-healthy");
			writer.write(healthy, false);
			const auto reference = Read(healthy);
			Require(reference.size() > 64 && reference.size() < 4096,
				"small fixture must exceed file limit while fitting in a stream buffer");
			for (const bool large : {false, true}) {
				if (large && !writer.supports_large) continue;
				++injected;
				const auto path = directory/(writer.name+(large ? "-stream-failure" : "-close-failure"));
				std::cout.flush();
				const auto child = fork();
				Require(child >= 0, "fork failed");
				if (child == 0) {
					// A real regular file becomes unwritable after 64 bytes. Ignore
					// SIGXFSZ so the writer must observe the returned I/O failure.
					if (std::signal(SIGXFSZ, SIG_IGN) == SIG_ERR) _exit(3);
					struct rlimit limit{};
					if (getrlimit(RLIMIT_FSIZE, &limit) != 0) _exit(3);
					limit.rlim_cur = 64;
					if (setrlimit(RLIMIT_FSIZE, &limit) != 0) _exit(3);
					try { writer.write(path, large); }
					catch (const std::exception& error) {
						_exit(std::string(error.what()).find(writer.diagnostic) != std::string::npos ? 0 : 2);
					}
					_exit(1); // Silent success despite incomplete output.
				}
				int status = 0;
				Require(waitpid(child, &status, 0) == child, "waitpid failed");
				Require(fs::is_regular_file(path) && fs::file_size(path) == 64,
					"fault did not produce the expected truncated regular file");
				const bool passed = WIFEXITED(status) && WEXITSTATUS(status) == 0;
				if (!passed) ++failures;
				std::cout << writer.name << (large ? " stream" : " close")
					<< " failure: " << (passed ? "PASS" : "FAIL") << '\n';
				// A fresh call after the failed write must produce the full bytes.
				writer.write(path, false);
				Require(Read(path) == reference, "retry changed healthy output bytes");
				++retries;
			}
			const auto invalid = directory/(writer.name+"-directory");
			fs::create_directory(invalid);
			bool rejected = false;
			try { writer.write(invalid, false); }
			catch (const std::exception&) { rejected = true; }
			Require(rejected, "directory output was accepted");
			++injected;
			std::cout << writer.name << " open failure and retries: PASS\n";
		}
		std::cout << "writers=" << writers.size() << " injected_failures=" << injected
			<< " retries=" << retries << " unexpected=" << failures << '\n';
		return failures ? 1 : 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
#else
	(void)argc;
	(void)argv;
	std::cerr << "SKIP: file-limit fault injection requires POSIX fork and RLIMIT_FSIZE\n";
	return 77;
#endif
}

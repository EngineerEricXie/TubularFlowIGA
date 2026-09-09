#include "CudaExecution.hpp"
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace {
int checks = 0;
const std::vector<const char*> names = {"OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE", "PMI_RANK", "PMI_SIZE",
	"SLURM_PROCID", "SLURM_NTASKS", "SLURM_STEP_NUM_TASKS", "SLURM_STEP_ID", "SLURM_STEPID",
	"OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"};

void Require(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
	++checks;
}

void Environment(const std::map<std::string, std::string>& values)
{
	for (const auto name : names) if (unsetenv(name)) throw std::runtime_error("unsetenv failed");
	for (const auto& value : values)
		if (setenv(value.first.c_str(), value.second.c_str(), 1)) throw std::runtime_error("setenv failed");
}

void Reject(const std::function<void()>& work, const char* expected)
{
	try { work(); }
	catch (const std::runtime_error& error) {
		Require(std::string(error.what()).find(expected) != std::string::npos, "wrong rejection diagnostic");
		return;
	}
	throw std::runtime_error("invalid execution accepted");
}
}

int main()
{
	try {
		const auto report = [] { std::ostringstream out; iga::cuda::RequireExecutionEnvironment(out); return out.str(); };
		Environment({});
		Require(report().find("launcher=direct") != std::string::npos, "missing direct launcher");
		for (const auto& pair : {std::pair<const char*, const char*>{"OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"},
			{"PMI_RANK", "PMI_SIZE"}, {"SLURM_PROCID", "SLURM_NTASKS"}}) {
			Environment({{pair.first, "0"}, {pair.second, "1"}});
			Require(report().find("processes=1") != std::string::npos, "single process rejected");
			for (const auto* rank : {"0", "1"}) {
				Environment({{pair.first, rank}, {pair.second, "2"}});
				Reject(report, "requires one task");
			}
			for (const auto* bad : {"", "-1", "+1", "1x", "1,2", "2147483648"}) {
				Environment({{pair.first, "0"}, {pair.second, bad}}); Reject(report, pair.second);
				Environment({{pair.first, bad}, {pair.second, "1"}}); Reject(report, pair.first);
			}
			Environment({{pair.first, "0"}, {pair.second, "0"}}); Reject(report, "positive");
			Environment({{pair.first, "2"}, {pair.second, "2"}}); Reject(report, "outside");
			Environment({{pair.first, "0"}}); Reject(report, pair.second);
		}
		Environment({{"SLURM_NTASKS", "8"}});
		Require(report().find("launcher=direct") != std::string::npos, "allocation mistaken for launch");
		for (const auto* step : {"batch", "extern"}) {
			Environment({{"SLURM_PROCID", "0"}, {"SLURM_NTASKS", "8"}, {"SLURM_STEP_ID", step}});
			Require(report().find("slurm_allocation") != std::string::npos, "batch shell rejected");
		}
		Environment({{"SLURM_PROCID", "0"}, {"SLURM_NTASKS", "8"}, {"SLURM_STEP_NUM_TASKS", "1"}});
		Require(report().find("slurm_step") != std::string::npos, "step count precedence lost");
		Environment({{"OMPI_COMM_WORLD_RANK", "0"}, {"OMPI_COMM_WORLD_SIZE", "1"},
			{"SLURM_PROCID", "2"}, {"SLURM_NTASKS", "8"}});
		Require(report().find("launcher=openmpi") != std::string::npos, "MPI precedence lost");
		Environment({{"OMP_NUM_THREADS", "4,2"}, {"OMP_THREAD_LIMIT", "8"}, {"OPENBLAS_NUM_THREADS", "0"}});
		Require(report().find("omp_requested=4") != std::string::npos, "nested thread request rejected");
		for (const auto* name : {"OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS"}) {
			Environment({{name, "bad"}}); Reject(report, name);
		}
		Environment({});
		std::ostringstream failed; failed.setstate(std::ios::badbit);
		Reject([&] { iga::cuda::RequireExecutionEnvironment(failed); }, "cannot write complete text stream");
		iga::cuda::RequireMeshIndexCapacity(INT_MAX, INT_MAX);
		Reject([] { iga::cuda::RequireMeshIndexCapacity(0, 1); }, "node count");
		Reject([] { iga::cuda::RequireMeshIndexCapacity(1, 0); }, "element count");
		Reject([] { iga::cuda::RequireMeshIndexCapacity(std::uint64_t(INT_MAX)+1, 1); }, "node count");
		Reject([] { iga::cuda::RequireMeshIndexCapacity(1, std::uint64_t(INT_MAX)+1); }, "element count");
		for (unsigned fields = 1; fields <= 8; ++fields) {
			const std::uint64_t maximum = INT_MAX/std::max(fields, 3u);
			iga::cuda::RequireFieldIndexCapacity(maximum, fields);
			Reject([&] { iga::cuda::RequireFieldIndexCapacity(maximum+1, fields); }, "node/field rows");
		}
		Reject([] { iga::cuda::RequireFieldIndexCapacity(1, 0); }, "field count");
		Reject([] { iga::cuda::RequireFieldIndexCapacity(1, 9); }, "field count");
		// Block value addresses use size_t, so a valid int32 block count must
		// not be rejected just because its scalar-value count exceeds INT_MAX.
		if (sizeof(std::size_t) >= 8)
			Require(iga::cuda::CheckedBlockValueCount(INT_MAX, 8) == std::size_t(INT_MAX)*64, "block offset range narrowed");
		Reject([] { iga::cuda::CheckedBlockValueCount(std::size_t(INT_MAX)+1, 4); }, "block matrix");
		Reject([] { iga::cuda::CheckedBlockValueCount(1, 0); }, "block matrix");
		std::cout << "cuda execution environment/index checks=" << checks << " passed\n";
		return 0;
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

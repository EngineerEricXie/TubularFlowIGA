#ifndef IGA_MEMORY_REPORT_HPP
#define IGA_MEMORY_REPORT_HPP

#include <petscsys.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

class DistributedMemoryRecorder {
public:
	DistributedMemoryRecorder(MPI_Comm communicator, const std::filesystem::path& path)
		: communicator_(communicator), path_(path)
	{
		MPI_Comm_rank(communicator_, &rank_);
		MPI_Comm_size(communicator_, &size_);
		if (path_.empty()) return;
		PetscCallThrow(PetscMemorySetGetMaximumUsage(), "PetscMemorySetGetMaximumUsage");
		if (rank_ == 0) {
			output_.open(path_, std::ios::trunc);
			if (!output_) throw std::runtime_error("cannot create memory report: "+path_.string());
			output_ << std::setprecision(17);
		}
	}

	void Record(const std::string& stage, int step = -1)
	{
		if (path_.empty()) return;
		PetscLogDouble petsc_current = 0.0, petsc_peak = 0.0;
		PetscLogDouble process_current = 0.0, process_peak = 0.0;
		PetscCallThrow(PetscMallocGetCurrentUsage(&petsc_current), "PetscMallocGetCurrentUsage");
		PetscCallThrow(PetscMallocGetMaximumUsage(&petsc_peak), "PetscMallocGetMaximumUsage");
		PetscCallThrow(PetscMemoryGetCurrentUsage(&process_current), "PetscMemoryGetCurrentUsage");
		PetscCallThrow(PetscMemoryGetMaximumUsage(&process_peak), "PetscMemoryGetMaximumUsage");
		const double local[6] = {ReadStatusBytes("VmRSS:"), ReadStatusBytes("VmHWM:"),
			static_cast<double>(petsc_current), static_cast<double>(petsc_peak),
			static_cast<double>(process_current), static_cast<double>(process_peak)};
		std::vector<double> gathered;
		if (rank_ == 0) gathered.resize(static_cast<std::size_t>(size_)*6);
		MPI_Gather(local, 6, MPI_DOUBLE, rank_ == 0 ? gathered.data() : nullptr,
			6, MPI_DOUBLE, 0, communicator_);
		if (rank_ != 0) return;

		const auto rss_max = Maximum(gathered, 0);
		const auto rss_sum = Sum(gathered, 0);
		const auto petsc_max = Maximum(gathered, 2);
		const auto petsc_sum = Sum(gathered, 2);
		std::cout << "memory stage=" << stage;
		if (step >= 0) std::cout << " step=" << step;
		std::cout << " rss_max_gib=" << rss_max/kGiB
			<< " rss_sum_gib=" << rss_sum/kGiB
			<< " petsc_alloc_max_gib=" << petsc_max/kGiB
			<< " petsc_alloc_sum_gib=" << petsc_sum/kGiB << '\n';
		output_ << "{\"stage\":\"" << stage << "\",\"step\":" << step
			<< ",\"mpi_ranks\":" << size_;
		WriteMetric("rss_current_bytes", gathered, 0);
		WriteMetric("rss_peak_bytes", gathered, 1);
		WriteMetric("petsc_alloc_current_bytes", gathered, 2);
		WriteMetric("petsc_alloc_peak_bytes", gathered, 3);
		WriteMetric("petsc_process_current_bytes", gathered, 4);
		WriteMetric("petsc_process_peak_bytes", gathered, 5);
		output_ << "}\n";
		output_.flush();
		if (!output_) throw std::runtime_error("cannot write memory report: "+path_.string());
	}

private:
	static constexpr double kGiB = 1024.0*1024.0*1024.0;

	static double ReadStatusBytes(const std::string& key)
	{
		std::ifstream input("/proc/self/status");
		std::string name;
		double kib = 0.0;
		std::string unit;
		while (input >> name) {
			if (name == key) {
				input >> kib >> unit;
				return kib*1024.0;
			}
			std::string rest;
			std::getline(input, rest);
		}
		return 0.0;
	}

	static double Maximum(const std::vector<double>& values, std::size_t metric)
	{
		double result = 0.0;
		for (std::size_t rank = 0; rank < values.size()/6; ++rank)
			result = std::max(result, values[rank*6+metric]);
		return result;
	}

	static double Sum(const std::vector<double>& values, std::size_t metric)
	{
		double result = 0.0;
		for (std::size_t rank = 0; rank < values.size()/6; ++rank)
			result += values[rank*6+metric];
		return result;
	}

	void WriteMetric(const char* name, const std::vector<double>& values, std::size_t metric)
	{
		output_ << ",\"" << name << "_per_rank\":[";
		for (std::size_t rank = 0; rank < values.size()/6; ++rank) {
			if (rank) output_ << ',';
			output_ << values[rank*6+metric];
		}
		output_ << "],\"" << name << "_max\":" << Maximum(values, metric)
			<< ",\"" << name << "_sum\":" << Sum(values, metric);
	}

	static void PetscCallThrow(PetscErrorCode code, const char* operation)
	{
		if (code != 0)
			throw std::runtime_error(std::string(operation)+" failed with PETSc error "+std::to_string(code));
	}

	MPI_Comm communicator_;
	std::filesystem::path path_;
	std::ofstream output_;
	int rank_ = 0;
	int size_ = 1;
};

} // namespace iga

#endif

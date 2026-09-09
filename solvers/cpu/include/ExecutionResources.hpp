#ifndef IGA_EXECUTION_RESOURCES_HPP
#define IGA_EXECUTION_RESOURCES_HPP

#include "CollectiveFailure.hpp"
#include "CheckedText.hpp"
#include <petscksp.h>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cctype>
#include <limits>
#include <ostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace iga {

// Parse documented numeric resource settings without changing environment or
// runtime policy. -1 means unset; a BLAS zero is retained as a library default.
inline int ResourceThreadSetting(const char* name, bool list = false, bool allow_zero = false)
{
	const auto* value = std::getenv(name);
	if (!value) return -1;
	int first = -1;
	for (;;) {
		while (std::isspace(static_cast<unsigned char>(*value))) ++value;
		char* end = nullptr;
		errno = 0;
		const auto number = std::strtol(value, &end, 10);
		if (end == value || errno == ERANGE || number < (allow_zero ? 0 : 1) || number > INT_MAX)
			throw std::runtime_error(std::string(name)+" requires "+(allow_zero ? "nonnegative" : "positive")+" integer thread counts");
		if (first < 0) first = static_cast<int>(number);
		while (std::isspace(static_cast<unsigned char>(*end))) ++end;
		if (!*end) return first;
		if (!list || *end != ',') throw std::runtime_error(std::string(name)+" has an invalid thread count list");
		value = end+1;
	}
}

inline void RequirePetscRealDouble()
{
#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE)
	throw std::runtime_error("this CPU backend requires a real, double-precision PETSc build");
#endif
}

inline void ValidatePackedExecution(std::uint64_t partitions, std::uint64_t nodes,
	std::uint64_t fields, int ranks)
{
	RequirePetscRealDouble();
	if (ranks < 1 || partitions != static_cast<std::uint64_t>(ranks))
		throw std::runtime_error("database partition count "+std::to_string(partitions)
			+" does not match communicator size "+std::to_string(ranks));
	if (!fields || nodes > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max())/fields)
		throw std::runtime_error("database node/field rows exceed this PETSc index capacity");
}

namespace execution_resources_detail {

inline std::string FormatReport(int ranks, const std::array<int, 6>& minimum,
	const std::array<int, 6>& maximum)
{
	std::ostringstream line;
	line.exceptions(std::ios::badbit | std::ios::failbit);
	line << "execution_resources ranks=" << ranks << " petsc_int_bits=" << sizeof(PetscInt)*CHAR_BIT
		<< " petsc_scalar=real" << sizeof(PetscScalar)*CHAR_BIT;
	const char* names[] = {"mpi_thread", "openmp_compiled", "omp_max_threads",
		"openblas_requested", "mkl_requested", "blis_requested"};
	for (std::size_t i = 0; i < minimum.size(); ++i)
		line << ' ' << names[i] << "_min=" << minimum[i] << ' ' << names[i] << "_max=" << maximum[i];
#ifdef PETSC_HAVE_MUMPS
	line << " petsc_mumps=1";
#else
	line << " petsc_mumps=0";
#endif
#ifdef PETSC_HAVE_HYPRE
	line << " petsc_hypre=1";
#else
	line << " petsc_hypre=0";
#endif
	return line.str();
}

} // namespace execution_resources_detail

// Collective, called from the MPI initialization thread before numerical work.
// Unequal but legal rank-local thread counts are allowed and reported as ranges.
// The report distinguishes compiled OpenMP/effective limits from BLAS requests;
// neither proves actual worker counts, core placement, or speedup.
inline void RequireExecutionResources(MPI_Comm communicator, std::ostream* report = nullptr)
{
	std::array<int, 6> local{};
	std::string abi;
	CollectiveLocalStage(communicator, "execution resource preflight", [&] {
		RequirePetscRealDouble();
		if (MPI_Query_thread(&local[0]) != MPI_SUCCESS) throw std::runtime_error("cannot query MPI thread support");
		(void)ResourceThreadSetting("OMP_NUM_THREADS", true);
		(void)ResourceThreadSetting("OMP_THREAD_LIMIT");
#ifdef _OPENMP
		local[1] = 1;
		local[2] = omp_get_max_threads();
		if (local[2] > 1 && local[0] < MPI_THREAD_FUNNELED)
			throw std::runtime_error("OpenMP with multiple threads requires MPI_THREAD_FUNNELED or stronger; initialize MPI accordingly");
#else
		local[1] = 0;
		local[2] = 1;
#endif
		local[3] = ResourceThreadSetting("OPENBLAS_NUM_THREADS", false, true);
		local[4] = ResourceThreadSetting("MKL_NUM_THREADS", false, true);
		local[5] = ResourceThreadSetting("BLIS_NUM_THREADS", false, true);
		abi = std::to_string(sizeof(PetscInt)*CHAR_BIT)+" "+std::to_string(sizeof(PetscScalar)*CHAR_BIT);
	});
	RequireCollectiveSameText(communicator, "execution PETSc ABI agreement", abi);
	auto minimum = local, maximum = local;
	MPI_Allreduce(local.data(), minimum.data(), static_cast<int>(local.size()), MPI_INT, MPI_MIN, communicator);
	MPI_Allreduce(local.data(), maximum.data(), static_cast<int>(local.size()), MPI_INT, MPI_MAX, communicator);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank); MPI_Comm_size(communicator, &ranks);
	CollectiveLocalStage(communicator, "execution resource report", [&] {
		if (rank != 0 || !report) return;
		*report << execution_resources_detail::FormatReport(ranks, minimum, maximum) << '\n';
		if (!*report) throw std::runtime_error("cannot write execution resource report");
		FlushCheckedText(*report);
	});
}

// Query actual matrix/backend compatibility, not just a compile-time package
// flag. This local PETSc query precedes collective factorization/solve.
inline void RequireFactorBackend(Mat matrix, MPI_Comm communicator, const char* backend,
	MatFactorType factor = MAT_FACTOR_LU)
{
	std::string description;
	CollectiveLocalStage(communicator, "factor backend preparation", [&] {
		if (!backend || !*backend) throw std::runtime_error("factor backend name is empty");
		description = std::string(backend)+" "+std::to_string(static_cast<int>(factor));
	});
	RequireCollectiveSameText(communicator, "factor backend agreement", description);
	CollectiveLocalStage(communicator, "factor backend availability", [&] {
		PetscBool available = PETSC_FALSE;
		const auto error = MatGetFactorAvailable(matrix, backend, factor, &available);
		if (error || !available)
			throw std::runtime_error("factor backend "+std::string(backend)+" is unavailable for this matrix/build; configure a compatible PETSc solver backend");
	});
}

inline void RequireKspFactorBackend(KSP solver, Mat matrix, MPI_Comm communicator)
{
	std::string backend, description;
	MatFactorType factor = MAT_FACTOR_LU;
	CollectiveLocalStage(communicator, "KSP backend preparation", [&] {
		PC pc = nullptr;
		PCType type = nullptr;
		if (KSPGetPC(solver, &pc) || PCGetType(pc, &type) || !type)
			throw std::runtime_error("cannot query KSP preconditioner");
		description = type;
		if (description == PCLU || description == PCCHOLESKY) {
			MatSolverType selected = nullptr;
			if (PCFactorGetMatSolverType(pc, &selected)) throw std::runtime_error("cannot query factor backend");
			if (selected) backend = selected;
			if (description == PCCHOLESKY) factor = MAT_FACTOR_CHOLESKY;
		}
		description += " "+backend;
	});
	RequireCollectiveSameText(communicator, "KSP backend agreement", description);
	// A null backend asks PETSc to choose later. Do not replace that policy.
	if (!backend.empty()) RequireFactorBackend(matrix, communicator, backend.c_str(), factor);
}

} // namespace iga

#endif

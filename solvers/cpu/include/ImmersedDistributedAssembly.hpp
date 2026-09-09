#ifndef IGA_IMMERSED_DISTRIBUTED_ASSEMBLY_HPP
#define IGA_IMMERSED_DISTRIBUTED_ASSEMBLY_HPP

#include "CollectiveFailure.hpp"
#include "ExecutionResources.hpp"
#include "PetscReadArray.hpp"
#include "RuntimeCleanup.hpp"
#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <set>
#include <sstream>
#include <vector>

namespace iga {

// Cell volume/wall terms use Dense; cubic ghost penalties couple the same
// field only; a gauge or port stencil couples its last (scalar) row to the
// other rows. Gauge/port stencils are split by cell so no owner needs the
// complete pressure/velocity field just to assemble one scalar constraint.
enum class ImmersedStencilPattern { Dense, SameField, Scalar };
struct ImmersedAssemblyStencil {
	std::uint64_t id = 0;
	int owner = 0;
	ImmersedStencilPattern pattern = ImmersedStencilPattern::Dense;
	std::vector<PetscInt> rows;

	bool Couples(std::size_t row, std::size_t column) const
	{
		if (pattern == ImmersedStencilPattern::Dense) return true;
		if (pattern == ImmersedStencilPattern::SameField) return rows[row]%4 == rows[column]%4;
		return row+1 == rows.size() || column+1 == rows.size();
	}
};

// Symbolic topology is initially replicated, numerical integration is unique
// to the stencil owner, and Mat/Vec storage is distributed by row offsets.
// The borrowed communicator outlives this object. Construction and every
// public operation that performs MPI must be entered by all its members.
class ImmersedDistributedAssembly {
public:
	ImmersedDistributedAssembly(MPI_Comm communicator, const std::vector<PetscInt>& row_offsets,
		const std::vector<ImmersedAssemblyStencil>& stencils)
		: communicator_(communicator)
	{
		MPI_Comm_rank(communicator_, &rank_);
		MPI_Comm_size(communicator_, &size_);
		std::string signature;
		std::vector<std::vector<PetscInt>> adjacency;
		std::vector<PetscInt> diagonal, remote;
		CollectiveLocalStage(communicator_, "immersed distributed topology", [&] {
			RequirePetscRealDouble();
			if (row_offsets.size() != static_cast<std::size_t>(size_)+1 || row_offsets.front() != 0
				|| row_offsets.back() <= 0 || !std::is_sorted(row_offsets.begin(), row_offsets.end()))
				throw std::invalid_argument("invalid immersed owned-row offsets");
			begin_ = row_offsets[rank_]; end_ = row_offsets[rank_+1]; rows_ = row_offsets.back();
			adjacency.resize(static_cast<std::size_t>(end_-begin_));
			std::ostringstream text;
			text.exceptions(std::ios::badbit | std::ios::failbit);
			for (const auto offset : row_offsets) text << offset << ' ';
			text << '|';
			std::set<std::uint64_t> ids;
			for (const auto& stencil : stencils) {
				if (stencil.owner < 0 || stencil.owner >= size_ || stencil.rows.empty()
					|| stencil.rows.size() > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max())
					|| !ids.insert(stencil.id).second)
					throw std::invalid_argument("invalid or duplicate immersed stencil");
				if (stencil.pattern != ImmersedStencilPattern::Dense && stencil.pattern != ImmersedStencilPattern::SameField
					&& stencil.pattern != ImmersedStencilPattern::Scalar)
					throw std::invalid_argument("invalid immersed stencil pattern");
				std::set<PetscInt> unique;
				text << stencil.id << ':' << stencil.owner << ':' << static_cast<int>(stencil.pattern) << ':';
				for (const auto row : stencil.rows) {
					if (row < 0 || row >= rows_ || !unique.insert(row).second)
						throw std::invalid_argument("invalid or duplicate immersed stencil row");
					text << row << ',';
				}
				text << ';';
				if (stencil.owner == rank_) {
					owned_.push_back(stencil);
					required_.insert(required_.end(), stencil.rows.begin(), stencil.rows.end());
				}
				for (std::size_t i = 0; i < stencil.rows.size(); ++i) {
					const auto row = stencil.rows[i];
					if (row < begin_ || row >= end_) continue;
					for (std::size_t j = 0; j < stencil.rows.size(); ++j)
						if (stencil.Couples(i, j)) adjacency[static_cast<std::size_t>(row-begin_)].push_back(stencil.rows[j]);
				}
			}
			signature = text.str();
			std::sort(required_.begin(), required_.end());
			required_.erase(std::unique(required_.begin(), required_.end()), required_.end());
			if (required_.size() > static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()))
				throw std::overflow_error("immersed halo exceeds PetscInt capacity");
			halo_values_.resize(required_.size());
			for (std::size_t i = 0; i < adjacency.size(); ++i) {
				auto& columns = adjacency[i];
				columns.push_back(begin_+static_cast<PetscInt>(i));
				std::sort(columns.begin(), columns.end());
				columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
				const auto local = std::count_if(columns.begin(), columns.end(), [&](PetscInt column) { return column >= begin_ && column < end_; });
				diagonal.push_back(static_cast<PetscInt>(local));
				remote.push_back(static_cast<PetscInt>(columns.size()-static_cast<std::size_t>(local)));
				symbolic_entries_ += columns.size();
			}
		});
		RequireCollectiveSameText(communicator_, "immersed topology agreement", signature);
		try {
			Check("immersed matrix create", MatCreateAIJ(communicator_, end_-begin_, end_-begin_, rows_, rows_,
				0, diagonal.data(), 0, remote.data(), &matrix_));
			Check("immersed vector create", VecCreateMPI(communicator_, end_-begin_, rows_, &state_));
			Check("immersed residual create", VecDuplicate(state_, &residual_));
			Check("immersed initial state", VecSet(state_, 0.0));
			Check("immersed structural zeros", MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE));
			CollectiveLocalStage(communicator_, "immersed symbolic insertion", [&] {
				for (std::size_t i = 0; i < adjacency.size(); ++i) {
					const auto row = begin_+static_cast<PetscInt>(i);
					std::vector<PetscScalar> zeros(adjacency[i].size(), 0.0);
					if (MatSetValues(matrix_, 1, &row, static_cast<PetscInt>(zeros.size()), adjacency[i].data(), zeros.data(), INSERT_VALUES))
						throw std::runtime_error("cannot seed immersed matrix structure");
				}
			});
			Check("immersed symbolic begin", MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY));
			Check("immersed symbolic end", MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY));
			Check("immersed exact allocation", MatSetOption(matrix_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
			Check("immersed fixed structure", MatSetOption(matrix_, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE));
			Check("immersed numerical zeros", MatSetOption(matrix_, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE));
			Check("immersed halo vector", VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(required_.size()), &halo_));
			Check("immersed halo rows", ISCreateGeneral(PETSC_COMM_SELF, static_cast<PetscInt>(required_.size()), required_.data(), PETSC_COPY_VALUES, &source_rows_));
			Check("immersed halo offsets", ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(required_.size()), 0, 1, &target_rows_));
			Check("immersed halo create", VecScatterCreate(state_, source_rows_, halo_, target_rows_, &scatter_));
		} catch (...) { Release(); throw; }
	}

	~ImmersedDistributedAssembly() { Release(); }
	ImmersedDistributedAssembly(const ImmersedDistributedAssembly&) = delete;
	ImmersedDistributedAssembly& operator=(const ImmersedDistributedAssembly&) = delete;
	PetscInt RowBegin() const noexcept { return begin_; }
	PetscInt RowEnd() const noexcept { return end_; }
	PetscInt Rows() const noexcept { return rows_; }
	std::size_t SymbolicEntries() const noexcept { return symbolic_entries_; }
	const std::vector<PetscInt>& RequiredRows() const noexcept { return required_; }
	const std::vector<ImmersedAssemblyStencil>& OwnedStencils() const noexcept { return owned_; }
	Vec State() const { RequireOpen(); return state_; }
	Mat Matrix() const { RequireAssembled(); return matrix_; }
	Vec Residual() const { RequireAssembled(); return residual_; }

	// Valid inside the local assembly callback, after the required-state scatter.
	double StateAt(PetscInt row) const
	{
		const auto found = std::lower_bound(required_.begin(), required_.end(), row);
		if (found == required_.end() || *found != row) throw std::out_of_range("row is outside immersed halo");
		return halo_values_[static_cast<std::size_t>(found-required_.begin())];
	}

	// callback(stencil, matrix, negative_residual) performs only local work.
	// An interrupted callback still drains all Mat/Vec stashes collectively;
	// its partial operator is inaccessible, and the next call starts from zero.
	template <class Function>
	void Assemble(Function&& function)
	{
		CollectiveLocalStage(communicator_, "immersed assembly preparation", [&] { RequireOpen(); });
		assembled_ = false;
		Check("immersed halo begin", VecScatterBegin(scatter_, state_, halo_, INSERT_VALUES, SCATTER_FORWARD));
		Check("immersed halo end", VecScatterEnd(scatter_, state_, halo_, INSERT_VALUES, SCATTER_FORWARD));
		CollectiveLocalStage(communicator_, "immersed halo values", [&] {
			PetscReadArray view;
			view.Acquire(halo_);
			for (std::size_t i = 0; i < required_.size(); ++i) {
				halo_values_[i] = PetscRealPart(view.Data()[i]);
				if (!std::isfinite(halo_values_[i])) throw std::runtime_error("nonfinite immersed halo state");
			}
			view.Restore();
		});
		Check("immersed matrix reset", MatZeroEntries(matrix_));
		Check("immersed residual reset", VecSet(residual_, 0.0));
		std::exception_ptr error;
		try {
			CollectiveLocalStage(communicator_, "immersed owned integration", [&] {
				std::vector<PetscScalar> matrix, residual;
				for (const auto& stencil : owned_) {
					matrix.clear(); residual.clear();
					function(stencil, matrix, residual);
					const auto n = stencil.rows.size();
					if (n > std::numeric_limits<std::size_t>::max()/n || matrix.size() != n*n || residual.size() != n)
						throw std::invalid_argument("immersed local block has an invalid size");
					for (std::size_t i = 0; i < n; ++i) {
						if (!std::isfinite(PetscRealPart(residual[i]))) throw std::runtime_error("nonfinite immersed local residual");
						for (std::size_t j = 0; j < n; ++j)
							if (!std::isfinite(PetscRealPart(matrix[i*n+j])) || (matrix[i*n+j] != 0.0 && !stencil.Couples(i, j)))
								throw std::runtime_error("immersed local matrix violates finite-value or stencil pattern contract");
					}
					if (MatSetValues(matrix_, static_cast<PetscInt>(n), stencil.rows.data(), static_cast<PetscInt>(n), stencil.rows.data(), matrix.data(), ADD_VALUES)
						|| VecSetValues(residual_, static_cast<PetscInt>(n), stencil.rows.data(), residual.data(), ADD_VALUES))
						throw std::runtime_error("cannot insert immersed local block");
				}
			});
		} catch (...) { error = std::current_exception(); }
		Check("immersed matrix assembly begin", MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY));
		Check("immersed matrix assembly end", MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY));
		Check("immersed residual assembly begin", VecAssemblyBegin(residual_));
		Check("immersed residual assembly end", VecAssemblyEnd(residual_));
		if (error) std::rethrow_exception(error);
		assembled_ = true;
	}

	void Close()
	{
		Release();
		cleanup_.Check(communicator_, "immersed assembly close");
	}
private:
	void RequireOpen() const { if (closed_) throw std::logic_error("immersed assembly is closed"); }
	void RequireAssembled() const { RequireOpen(); if (!assembled_) throw std::logic_error("immersed operator is not assembled"); }
	void Check(const char* stage, PetscErrorCode code) const { RequireCollectivePetscSuccess(communicator_, stage, code); }
	void Release() noexcept
	{
		if (closed_) return;
		closed_ = true;
		if (scatter_) { cleanup_.Observe("VecScatterDestroy", VecScatterDestroy(&scatter_)); scatter_ = nullptr; }
		if (source_rows_) { cleanup_.Observe("ISDestroy source", ISDestroy(&source_rows_)); source_rows_ = nullptr; }
		if (target_rows_) { cleanup_.Observe("ISDestroy target", ISDestroy(&target_rows_)); target_rows_ = nullptr; }
		if (halo_) { cleanup_.Observe("VecDestroy halo", VecDestroy(&halo_)); halo_ = nullptr; }
		if (residual_) { cleanup_.Observe("VecDestroy residual", VecDestroy(&residual_)); residual_ = nullptr; }
		if (state_) { cleanup_.Observe("VecDestroy state", VecDestroy(&state_)); state_ = nullptr; }
		if (matrix_) { cleanup_.Observe("MatDestroy", MatDestroy(&matrix_)); matrix_ = nullptr; }
	}
	MPI_Comm communicator_;
	int rank_ = 0, size_ = 1;
	PetscInt begin_ = 0, end_ = 0, rows_ = 0;
	std::size_t symbolic_entries_ = 0;
	std::vector<ImmersedAssemblyStencil> owned_;
	std::vector<PetscInt> required_;
	std::vector<double> halo_values_;
	Mat matrix_ = nullptr;
	Vec state_ = nullptr, residual_ = nullptr, halo_ = nullptr;
	IS source_rows_ = nullptr, target_rows_ = nullptr;
	VecScatter scatter_ = nullptr;
	RuntimeCleanupResult cleanup_;
	bool assembled_ = false, closed_ = false;
};

} // namespace iga
#endif

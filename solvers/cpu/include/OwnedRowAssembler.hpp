#ifndef OWNED_ROW_ASSEMBLER_HPP
#define OWNED_ROW_ASSEMBLER_HPP

#include "IgaDatabase.hpp"

#include <petscmat.h>
#include <petscvec.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace iga {

class OwnedRowAssembler {
public:
	OwnedRowAssembler(Database& database, MPI_Comm communicator, PetscInt fields)
		: database_(database), communicator_(communicator), fields_(fields)
	{
		if (fields_ <= 0) throw std::invalid_argument("fields must be positive");
		MPI_Comm_rank(communicator_, &rank_);
		MPI_Comm_size(communicator_, &size_);
		if (database_.header().ranks != static_cast<std::uint32_t>(size_))
			throw std::runtime_error("database rank count does not match MPI size");
		auto range = database_.NodeRange(rank_);
		node_begin_ = range.first;
		node_end_ = range.second;
		local_elements_ = database_.LoadRequired(rank_);
	}

	PetscInt fields() const { return fields_; }
	std::uint64_t node_begin() const { return node_begin_; }
	std::uint64_t node_end() const { return node_end_; }
	PetscInt local_rows() const { return CheckedPetsc((node_end_ - node_begin_) * fields_); }
	PetscInt global_rows() const { return CheckedPetsc(database_.header().nodes * fields_); }
	const std::vector<Element>& elements() const { return local_elements_; }

	Mat CreateMatrix(bool keep_nonzero_pattern = false) const
	{
		const auto local_nodes = static_cast<std::size_t>(node_end_ - node_begin_);
		std::vector<std::vector<std::int32_t>> adjacency(local_nodes);
		for (const auto& element : local_elements_) {
			for (const auto row_node : element.connectivity) {
				if (!Owns(row_node)) continue;
				auto& neighbors = adjacency[static_cast<std::size_t>(row_node - node_begin_)];
				neighbors.insert(neighbors.end(), element.connectivity.begin(), element.connectivity.end());
			}
		}

		std::vector<PetscInt> diagonal(static_cast<std::size_t>(local_rows()), 0);
		std::vector<PetscInt> off_diagonal(static_cast<std::size_t>(local_rows()), 0);
		for (std::size_t local_node = 0; local_node < adjacency.size(); ++local_node) {
			auto& neighbors = adjacency[local_node];
			std::sort(neighbors.begin(), neighbors.end());
			neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
			PetscInt local_count = 0;
			for (auto node : neighbors) if (Owns(node)) ++local_count;
			const auto remote_count = static_cast<PetscInt>(neighbors.size()) - local_count;
			for (PetscInt field = 0; field < fields_; ++field) {
				const auto row = local_node * static_cast<std::size_t>(fields_) + static_cast<std::size_t>(field);
				diagonal[row] = local_count * fields_;
				off_diagonal[row] = remote_count * fields_;
			}
		}

		Mat matrix = nullptr;
		PetscCallThrow(MatCreateAIJ(communicator_, local_rows(), local_rows(), global_rows(), global_rows(),
			0, diagonal.data(), 0, off_diagonal.data(), &matrix), "MatCreateAIJ");
		PetscCallThrow(MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE), "MatSetOption");
		if (keep_nonzero_pattern)
			PetscCallThrow(MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE), "MatSetOption");
		PetscCallThrow(MatSetOption(matrix, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE), "MatSetOption");
		return matrix;
	}

	Vec CreateVector() const
	{
		Vec vector = nullptr;
		PetscCallThrow(VecCreateMPI(communicator_, local_rows(), global_rows(), &vector), "VecCreateMPI");
		return vector;
	}

	void AddElementMatrix(Mat matrix, const Element& element, const std::vector<PetscScalar>& values) const
	{
		const auto nodes = element.connectivity.size();
		const auto element_rows = nodes * static_cast<std::size_t>(fields_);
		if (values.size() != element_rows * element_rows)
			throw std::invalid_argument("element matrix has the wrong size");
		std::vector<PetscInt> columns(element_rows);
		for (std::size_t node = 0; node < nodes; ++node)
			for (PetscInt field = 0; field < fields_; ++field)
				columns[node * fields_ + field] = GlobalRow(element.connectivity[node], field);

		std::vector<PetscInt> rows;
		std::vector<PetscScalar> owned_values;
		rows.reserve(element_rows);
		owned_values.reserve(element_rows * element_rows);
		for (std::size_t node = 0; node < nodes; ++node) {
			if (!Owns(element.connectivity[node])) continue;
			for (PetscInt field = 0; field < fields_; ++field) {
				const auto local_element_row = node * fields_ + field;
				rows.push_back(columns[local_element_row]);
				const auto begin = values.begin() + static_cast<std::ptrdiff_t>(local_element_row * element_rows);
				owned_values.insert(owned_values.end(), begin, begin + static_cast<std::ptrdiff_t>(element_rows));
			}
		}
		if (!rows.empty())
			PetscCallThrow(MatSetValues(matrix, static_cast<PetscInt>(rows.size()), rows.data(),
				static_cast<PetscInt>(columns.size()), columns.data(), owned_values.data(), ADD_VALUES), "MatSetValues");
	}

	void AddElementVector(Vec vector, const Element& element, const std::vector<PetscScalar>& values) const
	{
		const auto nodes = element.connectivity.size();
		if (values.size() != nodes * static_cast<std::size_t>(fields_))
			throw std::invalid_argument("element vector has the wrong size");
		std::vector<PetscInt> rows;
		std::vector<PetscScalar> owned_values;
		for (std::size_t node = 0; node < nodes; ++node) {
			if (!Owns(element.connectivity[node])) continue;
			for (PetscInt field = 0; field < fields_; ++field) {
				rows.push_back(GlobalRow(element.connectivity[node], field));
				owned_values.push_back(values[node * fields_ + field]);
			}
		}
		if (!rows.empty())
			PetscCallThrow(VecSetValues(vector, static_cast<PetscInt>(rows.size()), rows.data(),
				owned_values.data(), ADD_VALUES), "VecSetValues");
	}

	static void Assemble(Mat matrix)
	{
		PetscCallThrow(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
		PetscCallThrow(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
	}

	static void Assemble(Vec vector)
	{
		PetscCallThrow(VecAssemblyBegin(vector), "VecAssemblyBegin");
		PetscCallThrow(VecAssemblyEnd(vector), "VecAssemblyEnd");
	}

private:
	bool Owns(std::int64_t node) const
	{
		return node >= static_cast<std::int64_t>(node_begin_) && node < static_cast<std::int64_t>(node_end_);
	}

	PetscInt GlobalRow(std::int32_t node, PetscInt field) const
	{
		return CheckedPetsc(static_cast<std::uint64_t>(node) * fields_ + field);
	}

	static PetscInt CheckedPetsc(std::uint64_t value)
	{
		if (value > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max()))
			throw std::overflow_error("problem exceeds PetscInt range; rebuild PETSc with 64-bit indices");
		return static_cast<PetscInt>(value);
	}

	static void PetscCallThrow(PetscErrorCode code, const char* operation)
	{
		if (code != 0) throw std::runtime_error(std::string(operation) + " failed with PETSc error " + std::to_string(code));
	}

	Database& database_;
	MPI_Comm communicator_;
	PetscInt fields_;
	int rank_ = 0;
	int size_ = 1;
	std::uint64_t node_begin_ = 0;
	std::uint64_t node_end_ = 0;
	std::vector<Element> local_elements_;
};

} // namespace iga

#endif

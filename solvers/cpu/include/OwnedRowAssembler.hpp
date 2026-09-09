#ifndef OWNED_ROW_ASSEMBLER_HPP
#define OWNED_ROW_ASSEMBLER_HPP

#include "CollectiveFailure.hpp"
#include "RuntimeConstruction.hpp"
#include "ExecutionResources.hpp"
#include "FieldCouplingPattern.hpp"
#include "IgaDatabase.hpp"
#include "ParallelOwnershipValidation.hpp"

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
	// Borrowed communicator: caller keeps it valid until this assembler and
	// every PETSc object it creates have been destroyed. No duplication/free.
	OwnedRowAssembler(Database& database, MPI_Comm communicator, PetscInt fields)
		: database_(database), communicator_(communicator), fields_(fields)
	{
		MPI_Comm_rank(communicator_, &rank_);
		MPI_Comm_size(communicator_, &size_);
		CollectiveLocalStage(communicator_, "owned-row database input", [&] {
			RequirePetscRealDouble();
			if (fields_ <= 0) throw std::invalid_argument("fields must be positive");
			if (database_.header().ranks != static_cast<std::uint32_t>(size_))
				throw std::runtime_error("database rank count does not match MPI size");
			(void)global_rows();
			auto range = database_.NodeRange(rank_);
			node_begin_ = range.first;
			node_end_ = range.second;
			local_elements_ = database_.LoadRequired(rank_);
		});
		const std::array<std::uint64_t, 3> shape{{database_.header().nodes,
			database_.header().elements, static_cast<std::uint64_t>(fields_)}};
		auto root_shape = shape;
		{
			PhaseScope communication_phase(ProfilePhase::Communication);
			MPI_Bcast(root_shape.data(), static_cast<int>(root_shape.size()), MPI_UINT64_T,
				0, communicator_);
		}
		CollectiveLocalStage(communicator_, "owned-row global shape", [&] {
			if (shape != root_shape)
				throw std::runtime_error("node, element or field count differs from communicator rank 0");
		});
	}

	PetscInt fields() const { return fields_; }
	std::uint64_t node_begin() const { return node_begin_; }
	std::uint64_t node_end() const { return node_end_; }
	PetscInt local_rows() const { return RowCount(node_end_-node_begin_); }
	PetscInt global_rows() const { return RowCount(database_.header().nodes); }
	const std::vector<Element>& elements() const { return local_elements_; }
	OwnedNodeRange NodeOwnership() const { return {node_begin_, node_end_, database_.header().nodes}; }
	PetscInt GlobalRow(std::int32_t node, PetscInt field) const
	{
		return static_cast<PetscInt>(CheckedFieldRow(node, field, database_.header().nodes,
			fields_, static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max())));
	}
	std::vector<PetscInt> RequiredRows(const std::vector<std::int32_t>& nodes) const
	{
		ValidateRequiredNodeIds(nodes, database_.header().nodes);
		std::vector<PetscInt> rows;
		rows.reserve(static_cast<std::size_t>(RowCount(nodes.size())));
		for (auto node : nodes)
			for (PetscInt field = 0; field < fields_; ++field) rows.push_back(GlobalRow(node, field));
		return rows;
	}

	// Optional collective audit for validation runs. Assembly ownership and
	// physical integration ownership are checked separately.
	void ValidateOwnership() const
	{
		ValidateDistributedNodeRanges(NodeOwnership(), communicator_);
		std::vector<std::uint64_t> required, partition_owned, minimum_node_owned;
		std::vector<OwnershipIncidence> authoritative_rows, consumed_rows;
		bool valid = true;
		try {
		for (const auto& element : local_elements_) {
			required.push_back(element.id);
			valid = valid && element.id < database_.header().elements
				&& element.owner >= 0 && element.owner < size_ && !element.connectivity.empty();
			for (auto node : element.connectivity) {
				valid = valid && node >= 0 && static_cast<std::uint64_t>(node) < database_.header().nodes;
				if (Owns(node)) consumed_rows.emplace_back(element.id, node);
			}
			if (!element.connectivity.empty() && OwnsElementByMinimumNode(element))
				minimum_node_owned.push_back(element.id);
		}
		// The packed required-element index guarantees touching elements,
		// not all METIS-owned elements. Read the independent partition catalog
		// for the authoritative physical incidences, one record at a time.
		for (std::uint64_t index = 0; index < database_.header().elements; ++index) {
			if (database_.owners().at(static_cast<std::size_t>(index)) != rank_) continue;
			const auto element = database_.Load(index);
			valid = valid && element.id == index && element.owner == rank_
				&& !element.connectivity.empty();
			partition_owned.push_back(element.id);
			for (auto node : element.connectivity) {
				valid = valid && node >= 0 && static_cast<std::uint64_t>(node) < database_.header().nodes;
				authoritative_rows.emplace_back(element.id, node);
			}
		}
		} catch (const std::exception&) {
			valid = false;
		}
		std::sort(required.begin(), required.end());
		valid = valid && std::adjacent_find(required.begin(), required.end()) == required.end();
		RequireOwnershipCondition(valid, communicator_, "invalid, unreadable or duplicate element catalog");
		ValidateUniqueEntityCoverage(partition_owned, database_.header().elements,
			communicator_, "database-owned integration");
		ValidateUniqueEntityCoverage(minimum_node_owned, database_.header().elements,
			communicator_, "minimum-node integration");
		ValidateOwnershipIncidences(authoritative_rows, consumed_rows, communicator_, "owned-row assembly");
	}
	bool OwnsElementByMinimumNode(const Element& element) const
	{
		if (element.connectivity.empty())
			throw std::invalid_argument("cannot assign ownership to an empty element");
		return Owns(*std::min_element(element.connectivity.begin(), element.connectivity.end()));
	}

	Mat CreateMatrix(bool keep_nonzero_pattern = false) const
	{
		FieldCouplingPattern pattern;
		CollectiveLocalStage(communicator_, "owned-row dense field pattern", [&] {
			pattern = FieldCouplingPattern::Dense(static_cast<std::size_t>(fields_));
		});
		return CreateMatrix(pattern, keep_nonzero_pattern);
	}

	Mat CreateMatrix(const FieldCouplingPattern& pattern,
		bool keep_nonzero_pattern = false) const
	{
		RequireCollectiveSameInt(communicator_, "owned-row matrix pattern policy", keep_nonzero_pattern ? 1 : 0);
		std::vector<PetscInt> diagonal, off_diagonal;
		CollectiveLocalStage(communicator_, "owned-row matrix preallocation", [&] {
			if (pattern.fields() != static_cast<std::size_t>(fields_))
				throw std::invalid_argument("matrix field-coupling pattern has the wrong field count");
			const auto local_nodes = static_cast<std::size_t>(node_end_ - node_begin_);
			std::vector<std::vector<std::int32_t>> adjacency(local_nodes);
			for (const auto& element : local_elements_) {
				for (const auto row_node : element.connectivity) {
					if (!Owns(row_node)) continue;
					auto& neighbors = adjacency[static_cast<std::size_t>(row_node - node_begin_)];
					neighbors.insert(neighbors.end(), element.connectivity.begin(), element.connectivity.end());
				}
			}

			diagonal.assign(static_cast<std::size_t>(local_rows()), 0);
			off_diagonal.assign(static_cast<std::size_t>(local_rows()), 0);
			for (std::size_t local_node = 0; local_node < adjacency.size(); ++local_node) {
				auto& neighbors = adjacency[local_node];
				std::sort(neighbors.begin(), neighbors.end());
				neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
				PetscInt local_count = 0;
				for (auto node : neighbors) if (Owns(node)) ++local_count;
				const auto remote_count = static_cast<PetscInt>(neighbors.size()) - local_count;
				for (PetscInt field = 0; field < fields_; ++field) {
					const auto row = local_node * static_cast<std::size_t>(fields_) + static_cast<std::size_t>(field);
					const auto active_trials = CheckedPetsc(pattern.ActiveTrials(static_cast<std::size_t>(field)));
					diagonal[row] = local_count * active_trials;
					off_diagonal[row] = remote_count * active_trials;
				}
			}
		});

		Mat matrix = nullptr;
		try {
			RequireCollectivePetscSuccess(communicator_, "owned-row matrix creation",
				MatCreateAIJ(communicator_, local_rows(), local_rows(), global_rows(), global_rows(),
					0, diagonal.data(), 0, off_diagonal.data(), &matrix));
			ObserveConstructedObject(communicator_, "owned-row matrix created", reinterpret_cast<PetscObject>(matrix));
			RequireCollectivePetscSuccess(communicator_, "owned-row new nonzero policy",
				MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
			if (keep_nonzero_pattern)
				RequireCollectivePetscSuccess(communicator_, "owned-row nonzero pattern policy",
					MatSetOption(matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));
			RequireCollectivePetscSuccess(communicator_, "owned-row off-process policy",
				MatSetOption(matrix, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_FALSE));
			return matrix;
		} catch (...) {
			MatDestroy(&matrix);
			throw;
		}
	}

	Vec CreateVector() const
	{
		Vec vector = nullptr;
		try {
			RequireCollectivePetscSuccess(communicator_, "owned-row vector creation",
				VecCreateMPI(communicator_, local_rows(), global_rows(), &vector));
			ObserveConstructedObject(communicator_, "owned-row vector created", reinterpret_cast<PetscObject>(vector));
			return vector;
		} catch (...) {
			VecDestroy(&vector);
			throw;
		}
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

	void AddElementMatrix(Mat matrix, const Element& element,
		const FieldBlockElementMatrix& blocks) const
	{
		const auto nodes = element.connectivity.size();
		if (blocks.nodes() != nodes
			|| blocks.pattern().fields() != static_cast<std::size_t>(fields_))
			throw std::invalid_argument("element field blocks have the wrong size");
		sparse_rows_.clear();
		sparse_columns_.resize(nodes);
		sparse_values_.clear();
		sparse_rows_.reserve(nodes);
		sparse_values_.reserve(nodes*nodes);
		for (std::size_t pair = 0; pair < blocks.pattern().pairs(); ++pair) {
			const auto coupling = blocks.pattern().active_pairs()[pair];
			sparse_rows_.clear();
			sparse_values_.clear();
			for (std::size_t node = 0; node < nodes; ++node)
				sparse_columns_[node] = GlobalRow(element.connectivity[node],
					static_cast<PetscInt>(coupling.second));
			const auto* values = blocks.Block(pair);
			for (std::size_t node = 0; node < nodes; ++node) {
				if (!Owns(element.connectivity[node])) continue;
				sparse_rows_.push_back(GlobalRow(element.connectivity[node],
					static_cast<PetscInt>(coupling.first)));
				const auto begin = values+static_cast<std::ptrdiff_t>(node*nodes);
				sparse_values_.insert(sparse_values_.end(), begin,
					begin+static_cast<std::ptrdiff_t>(nodes));
			}
			if (!sparse_rows_.empty())
				PetscCallThrow(MatSetValues(matrix, static_cast<PetscInt>(sparse_rows_.size()),
					sparse_rows_.data(), static_cast<PetscInt>(sparse_columns_.size()),
					sparse_columns_.data(), sparse_values_.data(), ADD_VALUES),
					"MatSetValues");
		}
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

	// Coordinate each returned status before any peer enters the next PETSc
	// operation. The communicator must be the object's borrowed runtime group.
	static void Assemble(Mat matrix, MPI_Comm communicator)
	{
		RequireCollectivePetscSuccess(communicator, "matrix assembly begin", MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY));
		RequireCollectivePetscSuccess(communicator, "matrix assembly end", MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY));
	}

	static void Assemble(Vec vector, MPI_Comm communicator)
	{
		RequireCollectivePetscSuccess(communicator, "vector assembly begin", VecAssemblyBegin(vector));
		RequireCollectivePetscSuccess(communicator, "vector assembly end", VecAssemblyEnd(vector));
	}

private:
	bool Owns(std::int64_t node) const
	{
		return node >= static_cast<std::int64_t>(node_begin_) && node < static_cast<std::int64_t>(node_end_);
	}

	PetscInt RowCount(std::uint64_t nodes) const
	{
		return static_cast<PetscInt>(CheckedFieldRowCount(nodes, fields_,
			static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max())));
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
	mutable std::vector<PetscInt> sparse_rows_;
	mutable std::vector<PetscInt> sparse_columns_;
	mutable std::vector<PetscScalar> sparse_values_;
};

} // namespace iga

#endif

#include "ImmersedDistributedAssembly.hpp"
#include <iostream>
#include <memory>

namespace {
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

double State(PetscInt row) { return 0.03*std::sin(0.31*(row+1)); }
double Coefficient(const iga::ImmersedAssemblyStencil& stencil, std::size_t i, std::size_t j)
{
	return stencil.Couples(i,j) ? (i == j ? 1.0 : 0.125)/(1.0+stencil.id+i+j) : 0.0;
}

template <class Function>
void Reject(MPI_Comm comm, Function&& function, const char* expected)
{
	std::string message;
	try { function(); }
	catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm, "distributed assembly test rejection", [&] {
		Require(message.find(expected) != std::string::npos, "wrong or missing rejection");
	});
	iga::RequireCollectiveSameText(comm, "distributed assembly test failure agreement", message);
}

void Exercise(MPI_Comm comm, bool empty)
{
	int rank = 0, size = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
	constexpr PetscInt nodes = 17, rows = 4*nodes+2;
	std::vector<PetscInt> offsets;
	std::vector<iga::ImmersedAssemblyStencil> stencils;
	iga::CollectiveLocalStage(comm, "distributed assembly test topology", [&] {
		for (int r = 0; r <= size; ++r) offsets.push_back(empty ? (r ? rows : 0) : (r == size ? rows : 4*(nodes*r/size)));
		const auto add = [&](int cell, iga::ImmersedStencilPattern pattern, std::vector<PetscInt> indices) {
			stencils.push_back({stencils.size(), empty ? 0 : cell*size/16, pattern, std::move(indices)});
		};
		for (int cell = 0; cell < 16; ++cell) {
			std::vector<PetscInt> indices;
			for (PetscInt node = cell; node <= cell+1; ++node)
				for (PetscInt field = 0; field < 4; ++field) indices.push_back(4*node+field);
			add(cell, iga::ImmersedStencilPattern::Dense, indices);
			if (cell < 15) {
				for (PetscInt field = 0; field < 4; ++field) indices.push_back(4*(cell+2)+field);
				add(cell, iga::ImmersedStencilPattern::SameField, indices);
			}
			add(cell, iga::ImmersedStencilPattern::Scalar, {4*cell+3, 4*(cell+1)+3, 4*nodes});
			add(cell, iga::ImmersedStencilPattern::Scalar, {4*cell, 4*cell+1, 4*(cell+1), 4*(cell+1)+1, 4*nodes+1});
		}
	});
	iga::ImmersedDistributedAssembly assembly(comm, offsets, stencils);
	iga::CollectiveLocalStage(comm, "distributed assembly test state", [&] {
		for (PetscInt row = assembly.RowBegin(); row < assembly.RowEnd(); ++row) {
			const PetscScalar value = State(row);
			Require(!VecSetValue(assembly.State(), row, value, INSERT_VALUES), "state insertion failed");
		}
	});
	iga::RequireCollectivePetscSuccess(comm, "test state begin", VecAssemblyBegin(assembly.State()));
	iga::RequireCollectivePetscSuccess(comm, "test state end", VecAssemblyEnd(assembly.State()));
	const auto integrate = [&](const iga::ImmersedAssemblyStencil& stencil,
		std::vector<PetscScalar>& matrix, std::vector<PetscScalar>& residual) {
		const auto n = stencil.rows.size(); matrix.resize(n*n); residual.assign(n, 0.0);
		for (std::size_t i = 0; i < n; ++i) {
			Require(assembly.StateAt(stencil.rows[i]) == State(stencil.rows[i]), "halo value/order differs");
			for (std::size_t j = 0; j < n; ++j) {
				matrix[i*n+j] = Coefficient(stencil, i, j);
				residual[i] -= matrix[i*n+j]*assembly.StateAt(stencil.rows[j]);
			}
		}
	};
	std::vector<double> reference(rows, 0.0);
	for (const auto& stencil : stencils)
		for (std::size_t i = 0; i < stencil.rows.size(); ++i)
			for (std::size_t j = 0; j < stencil.rows.size(); ++j)
				reference[static_cast<std::size_t>(stencil.rows[i])] -= Coefficient(stencil, i, j)*State(stencil.rows[j]);
	const auto verify = [&] {
		iga::CollectiveLocalStage(comm, "distributed assembly test residual", [&] {
			iga::PetscReadArray view; view.Acquire(assembly.Residual());
			for (PetscInt row = assembly.RowBegin(); row < assembly.RowEnd(); ++row)
				Require(std::abs(view.Data()[row-assembly.RowBegin()]-reference[static_cast<std::size_t>(row)]) < 1e-13, "owned residual differs from independent global sum");
			view.Restore();
		});
		Vec action = nullptr;
		iga::RequireCollectivePetscSuccess(comm, "test action create", VecDuplicate(assembly.State(), &action));
		iga::RequireCollectivePetscSuccess(comm, "test matrix action", MatMult(assembly.Matrix(), assembly.State(), action));
		iga::RequireCollectivePetscSuccess(comm, "test matrix consistency", VecAXPY(action, 1.0, assembly.Residual()));
		PetscReal error = 0.0;
		iga::RequireCollectivePetscSuccess(comm, "test action norm", VecNorm(action, NORM_2, &error));
		iga::RequireCollectivePetscSuccess(comm, "test action destroy", VecDestroy(&action));
		iga::CollectiveLocalStage(comm, "distributed assembly test action gate", [&] { Require(error < 1e-13, "distributed matrix action differs"); });
	};
	assembly.Assemble(integrate); verify();
	MatInfo info{};
	iga::RequireCollectivePetscSuccess(comm, "test local matrix info", MatGetInfo(assembly.Matrix(), MAT_LOCAL, &info));
	iga::CollectiveLocalStage(comm, "distributed assembly test sparsity", [&] {
		Require(info.nz_used == static_cast<double>(assembly.SymbolicEntries()) && info.nz_allocated == info.nz_used && info.mallocs == 0.0,
			"matrix does not have exact preallocation");
		if (empty && rank) Require(assembly.OwnedStencils().empty() && assembly.RequiredRows().empty() && assembly.RowBegin() == assembly.RowEnd(), "empty rank owns work or state");
	});
	const int failed_rank = empty ? 0 : size-1;
	int visits = 0;
	Reject(comm, [&] { assembly.Assemble([&](const auto& stencil, auto& matrix, auto& residual) {
		if (rank == failed_rank && ++visits == 2) throw std::runtime_error("injected integration failure after insertion");
		integrate(stencil, matrix, residual);
	}); }, "injected integration failure after insertion");
	Reject(comm, [&] { (void)assembly.Matrix(); }, "operator is not assembled");
	assembly.Assemble(integrate); verify();
	Reject(comm, [&] { assembly.Assemble([&](const auto& stencil, auto& matrix, auto& residual) {
		integrate(stencil, matrix, residual);
		if (rank == failed_rank && stencil.pattern == iga::ImmersedStencilPattern::SameField) matrix[1] = 1.0;
	}); }, "stencil pattern contract");
	assembly.Assemble(integrate); verify();
	const auto owned = assembly.OwnedStencils().size(), halo = assembly.RequiredRows().size();
	assembly.Close(); assembly.Close();
	Reject(comm, [&] { (void)assembly.State(); }, "assembly is closed");
	for (int fault = 0; fault < 3; ++fault) {
		auto invalid = stencils;
		if (rank == size-1) {
			if (fault == 0) invalid.back().id = invalid.front().id;
			if (fault == 1) invalid.front().rows[0] = rows;
			if (fault == 2) invalid.front().rows[0] = 32;
		}
		if (fault == 2 && size == 1) continue;
		Reject(comm, [&] { iga::ImmersedDistributedAssembly rejected(comm, offsets, invalid); },
			fault == 0 ? "duplicate immersed stencil" : fault == 1 ? "invalid or duplicate immersed stencil row" : "input differs");
	}
	std::cout << "immersed_distributed_test rank=" << rank << " ranks=" << size << " empty=" << empty
		<< " owned_stencils=" << owned << " halo_rows=" << halo << " global_rows=" << rows
		<< " local_nz=" << info.nz_used << " passed\n";
}
}

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, size = 1, status = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank); MPI_Comm_size(PETSC_COMM_WORLD, &size);
	MPI_Comm comm = MPI_COMM_NULL;
	try {
		const bool split = argc > 1 && std::string(argv[1]) == "split";
		MPI_Comm_split(PETSC_COMM_WORLD, split && rank ? 1 : 0, rank, &comm);
		Exercise(comm, false); Exercise(comm, true);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status, &global, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

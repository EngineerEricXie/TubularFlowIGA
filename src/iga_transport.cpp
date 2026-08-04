#include "CaseInput.hpp"
#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"
#include "TransportElement.hpp"

#include <petscksp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "NeuronTransportIGA v2 transport solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc < 3) throw std::runtime_error("usage: iga_transport DATABASE.ntiga CASE_DIR [STEPS] [OUTPUT] [VELOCITY]");
		iga::Database database(argv[1]);
		const fs::path case_dir(argv[2]);
		auto parameters = iga::ReadTransportParameters((case_dir / "simulation_parameter.txt").string());
		if (argc >= 4) parameters.steps = std::stoi(argv[3]);
		const auto labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
		const auto velocity_path = argc >= 6 ? fs::path(argv[5]) : case_dir / "initial_velocityfield.txt";
		const auto velocity = iga::ReadVelocity(velocity_path.string(), database.header().nodes);
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 2);
		iga::RequireValidGeometry(assembler.elements(), rank, PETSC_COMM_WORLD);
		Mat left = assembler.CreateMatrix();
		Mat previous = assembler.CreateMatrix();
		MatSetOption(previous, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
		const auto assembly_start = std::chrono::steady_clock::now();
		for (const auto& element : assembler.elements()) {
			auto matrices = iga::BuildTransportElement(element, velocity, parameters);
			assembler.AddElementMatrix(left, element, matrices.left);
			assembler.AddElementMatrix(previous, element, matrices.previous);
		}
		iga::OwnedRowAssembler::Assemble(left);
		iga::OwnedRowAssembler::Assemble(previous);

		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
			if (labels[static_cast<std::size_t>(node)] == 1) {
				boundary_rows.push_back(static_cast<PetscInt>(2*node));
				boundary_rows.push_back(static_cast<PetscInt>(2*node+1));
			}
		MatZeroRows(left, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 1.0, nullptr, nullptr);
		MatZeroRows(previous, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 0.0, nullptr, nullptr);

		Vec current = assembler.CreateVector();
		Vec next = assembler.CreateVector();
		Vec rhs = assembler.CreateVector();
		VecSet(current, 0.0);
		VecSet(next, 0.0);
		VecSet(rhs, 0.0);
		PetscScalar* current_array = nullptr;
		VecGetArray(current, &current_array);
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto local = static_cast<std::size_t>(node - assembler.node_begin()) * 2;
			if (labels[static_cast<std::size_t>(node)] == 1) {
				current_array[local] = parameters.n0_bc;
				current_array[local+1] = parameters.nplus_bc;
			}
		}
		VecRestoreArray(current, &current_array);

		KSP solver = nullptr;
		KSPCreate(PETSC_COMM_WORLD, &solver);
		KSPSetOperators(solver, left, left);
		KSPSetType(solver, KSPGMRES);
		KSPGMRESSetRestart(solver, 50);
		KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
		PC preconditioner = nullptr;
		KSPGetPC(solver, &preconditioner);
		PCSetType(preconditioner, PCBJACOBI);
		KSPSetFromOptions(solver);
		KSPSetUp(solver);
		const auto assembly_end = std::chrono::steady_clock::now();

		PetscInt total_iterations = 0;
		const auto solve_start = std::chrono::steady_clock::now();
		for (int step = 0; step < parameters.steps; ++step) {
			MatMult(previous, current, rhs);
			std::vector<PetscScalar> boundary_values(boundary_rows.size());
			for (std::size_t i = 0; i < boundary_values.size(); i += 2) {
				boundary_values[i] = parameters.n0_bc;
				boundary_values[i+1] = parameters.nplus_bc;
			}
			VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(rhs);
			if (step > 0) VecCopy(current, next);
			KSPSetInitialGuessNonzero(solver, step > 0 ? PETSC_TRUE : PETSC_FALSE);
			KSPSolve(solver, rhs, next);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver, &reason);
			if (reason <= 0) throw std::runtime_error("transport KSP did not converge at step " + std::to_string(step));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver, &iterations);
			total_iterations += iterations;
			VecSwap(current, next);
		}
		const auto solve_end = std::chrono::steady_clock::now();
		PetscReal norm = 0.0;
		VecNorm(current, NORM_2, &norm);
		if (argc >= 5) {
			Vec root = nullptr;
			VecScatter scatter = nullptr;
			VecScatterCreateToZero(current, &scatter, &root);
			VecScatterBegin(scatter, current, root, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(scatter, current, root, INSERT_VALUES, SCATTER_FORWARD);
			if (rank == 0) {
				const PetscScalar* values = nullptr;
				VecGetArrayRead(root, &values);
				std::ofstream output(argv[4]);
				if (!output) throw std::runtime_error("cannot create final state output");
				output.precision(17);
				for (std::uint64_t node = 0; node < database.header().nodes; ++node)
					output << node << ' ' << PetscRealPart(values[2*node]) << ' ' << PetscRealPart(values[2*node+1]) << '\n';
				VecRestoreArrayRead(root, &values);
			}
			VecScatterDestroy(&scatter);
			VecDestroy(&root);
		}
		if (rank == 0) {
			const auto assembly_seconds = std::chrono::duration<double>(assembly_end-assembly_start).count();
			const auto solve_seconds = std::chrono::duration<double>(solve_end-solve_start).count();
			std::cout << "transport_v2 nodes=" << database.header().nodes << " elements=" << database.header().elements
				<< " steps=" << parameters.steps << " assembly_s=" << assembly_seconds << " solve_s=" << solve_seconds
				<< " total_iterations=" << total_iterations << " final_l2=" << norm << '\n';
		}
		KSPDestroy(&solver);
		VecDestroy(&rhs); VecDestroy(&next); VecDestroy(&current);
		MatDestroy(&previous); MatDestroy(&left);
	} catch (const std::exception& e) {
		std::cerr << "rank " << rank << ": " << e.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

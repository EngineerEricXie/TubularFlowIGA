#include "CaseInput.hpp"
#include "IgaDatabase.hpp"
#include "NavierStokesElement.hpp"
#include "OwnedRowAssembler.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "NeuronTransportIGA v2 stabilized steady Navier-Stokes solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc < 3) throw std::runtime_error("usage: iga_navier_stokes DATABASE.ntiga CASE_DIR [MAX_NEWTON] [OUTPUT]");
		iga::Database database(argv[1]);
		const fs::path case_dir(argv[2]);
		const int max_newton = argc >= 4 ? std::stoi(argv[3]) : 8;
		const auto parameters = iga::ReadTransportParameters((case_dir / "simulation_parameter.txt").string());
		const auto labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
		const auto boundary_velocity = iga::ReadVelocity((case_dir / "initial_velocityfield.txt").string(), database.header().nodes);
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 4);
		iga::RequireValidGeometry(assembler.elements(), rank, PETSC_COMM_WORLD);
		Mat jacobian = assembler.CreateMatrix(true);
		Vec state = assembler.CreateVector(), update = assembler.CreateVector(), rhs = assembler.CreateVector();
		VecSet(state, 0.0);
		VecSet(update, 0.0);
		VecSet(rhs, 0.0);

		std::vector<std::int32_t> ghost_nodes;
		for (const auto& element : assembler.elements()) ghost_nodes.insert(ghost_nodes.end(), element.connectivity.begin(), element.connectivity.end());
		std::sort(ghost_nodes.begin(), ghost_nodes.end());
		ghost_nodes.erase(std::unique(ghost_nodes.begin(), ghost_nodes.end()), ghost_nodes.end());
		std::unordered_map<std::int32_t, std::size_t> ghost_position;
		std::vector<PetscInt> ghost_rows;
		ghost_rows.reserve(4*ghost_nodes.size());
		for (std::size_t i = 0; i < ghost_nodes.size(); ++i) {
			ghost_position[ghost_nodes[i]] = i;
			for (int field = 0; field < 4; ++field) ghost_rows.push_back(4*ghost_nodes[i]+field);
		}
		IS source_rows = nullptr;
		ISCreateGeneral(PETSC_COMM_WORLD, static_cast<PetscInt>(ghost_rows.size()), ghost_rows.data(), PETSC_COPY_VALUES, &source_rows);
		Vec ghost_state = nullptr;
		VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(ghost_rows.size()), &ghost_state);
		IS destination_rows = nullptr;
		ISCreateStride(PETSC_COMM_SELF, static_cast<PetscInt>(ghost_rows.size()), 0, 1, &destination_rows);
		VecScatter scatter = nullptr;
		VecScatterCreate(state, source_rows, ghost_state, destination_rows, &scatter);

		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto label = labels[static_cast<std::size_t>(node)];
			if (label == 0 || label == 1)
				for (int field = 0; field < 3; ++field) boundary_rows.push_back(static_cast<PetscInt>(4*node+field));
			if (label >= 2) boundary_rows.push_back(static_cast<PetscInt>(4*node+3));
		}

		KSP solver = nullptr;
		KSPCreate(PETSC_COMM_WORLD, &solver);
		KSPSetType(solver, KSPFGMRES);
		KSPSetTolerances(solver, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
		PC pc = nullptr;
		KSPGetPC(solver, &pc);
		PCSetType(pc, PCBJACOBI);
		KSPSetFromOptions(solver);

		const double viscosity = 0.1;
		PetscInt total_linear_iterations = 0;
		PetscReal initial_residual = -1.0;
		bool converged = false;
		const auto start = std::chrono::steady_clock::now();
		for (int nonlinear = 0; nonlinear < max_newton; ++nonlinear) {
			const auto iteration_start = std::chrono::steady_clock::now();
			MatZeroEntries(jacobian);
			VecSet(rhs, 0.0);
			VecScatterBegin(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(scatter, state, ghost_state, INSERT_VALUES, SCATTER_FORWARD);
			const PetscScalar* ghost_values = nullptr;
			VecGetArrayRead(ghost_state, &ghost_values);
			for (const auto& element : assembler.elements()) {
				std::vector<std::array<double,4>> nodal(element.connectivity.size());
				for (std::size_t a = 0; a < element.connectivity.size(); ++a) {
					const auto position = ghost_position.at(element.connectivity[a]);
					for (int field = 0; field < 4; ++field) nodal[a][field] = PetscRealPart(ghost_values[4*position+field]);
				}
				auto local = iga::BuildNavierStokesElement(element, nodal, viscosity);
				assembler.AddElementMatrix(jacobian, element, local.jacobian);
				assembler.AddElementVector(rhs, element, local.negative_residual);
			}
			VecRestoreArrayRead(ghost_state, &ghost_values);
			iga::OwnedRowAssembler::Assemble(jacobian);
			iga::OwnedRowAssembler::Assemble(rhs);

			const PetscScalar* owned_state = nullptr;
			VecGetArrayRead(state, &owned_state);
			std::vector<PetscScalar> boundary_update;
			boundary_update.reserve(boundary_rows.size());
			for (auto row : boundary_rows) {
				const auto node = static_cast<std::uint64_t>(row / 4);
				const auto field = row % 4;
				double target = 0.0;
				if (labels[static_cast<std::size_t>(node)] == 1 && field < 3)
					target = parameters.vplus * boundary_velocity[static_cast<std::size_t>(node)][field];
				const auto local_row = static_cast<std::size_t>(row - 4*assembler.node_begin());
				boundary_update.push_back(target - PetscRealPart(owned_state[local_row]));
			}
			VecRestoreArrayRead(state, &owned_state);
			MatZeroRows(jacobian, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), 1.0, nullptr, nullptr);
			VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_update.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(rhs);
			PetscReal residual_norm = 0.0;
			VecNorm(rhs, NORM_2, &residual_norm);
			if (initial_residual < 0.0) initial_residual = residual_norm;
			const auto nonlinear_tolerance = std::max<PetscReal>(1e-10, 1e-5*initial_residual);
			if (residual_norm <= nonlinear_tolerance) {
				if (rank == 0) std::cout << "converged newton=" << nonlinear << " residual_l2=" << residual_norm
					<< " tolerance=" << nonlinear_tolerance << " assembly_s="
					<< std::chrono::duration<double>(std::chrono::steady_clock::now()-iteration_start).count() << '\n';
				converged = true;
				break;
			}
			const auto linear_start = std::chrono::steady_clock::now();
			KSPSetOperators(solver, jacobian, jacobian);
			KSPSolve(solver, rhs, update);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver, &reason);
			if (reason <= 0) throw std::runtime_error("Navier-Stokes linear solve failed at nonlinear iteration " + std::to_string(nonlinear));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver, &iterations);
			total_linear_iterations += iterations;
			PetscReal linear_residual = 0.0;
			KSPGetResidualNorm(solver, &linear_residual);
			PetscReal update_norm = 0.0;
			VecNorm(update, NORM_2, &update_norm);
			VecAXPY(state, 1.0, update);
			if (rank == 0) std::cout << "newton=" << nonlinear << " residual_l2=" << residual_norm << " update_l2=" << update_norm
				<< " linear_iterations=" << iterations << " linear_residual=" << linear_residual
				<< " assembly_s=" << std::chrono::duration<double>(linear_start-iteration_start).count()
				<< " linear_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-linear_start).count() << '\n';
			if (update_norm < 1e-10) { converged = true; break; }
		}
		if (!converged) throw std::runtime_error("Navier-Stokes nonlinear solve reached MAX_NEWTON without convergence");
		const auto end = std::chrono::steady_clock::now();
		PetscReal state_norm = 0.0;
		VecNorm(state, NORM_2, &state_norm);
		const PetscScalar* owned_values = nullptr;
		VecGetArrayRead(state, &owned_values);
		double local_velocity_squared = 0.0, local_pressure_squared = 0.0;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto local = static_cast<std::size_t>(node - assembler.node_begin()) * 4;
			for (int field = 0; field < 3; ++field) {
				const auto value = PetscRealPart(owned_values[local+field]);
				local_velocity_squared += value*value;
			}
			const auto pressure = PetscRealPart(owned_values[local+3]);
			local_pressure_squared += pressure*pressure;
		}
		VecRestoreArrayRead(state, &owned_values);
		double velocity_squared = 0.0, pressure_squared = 0.0;
		MPI_Allreduce(&local_velocity_squared, &velocity_squared, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		MPI_Allreduce(&local_pressure_squared, &pressure_squared, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		if (rank == 0) std::cout << "navier_stokes_v2 seconds=" << std::chrono::duration<double>(end-start).count()
			<< " total_linear_iterations=" << total_linear_iterations << " state_l2=" << state_norm
			<< " velocity_l2=" << std::sqrt(velocity_squared) << " pressure_l2=" << std::sqrt(pressure_squared) << '\n';

		if (argc >= 5) {
			Vec root = nullptr; VecScatter root_scatter = nullptr;
			VecScatterCreateToZero(state, &root_scatter, &root);
			VecScatterBegin(root_scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(root_scatter, state, root, INSERT_VALUES, SCATTER_FORWARD);
			if (rank == 0) {
				const PetscScalar* values = nullptr; VecGetArrayRead(root, &values);
				std::ofstream output(argv[4]);
				std::ofstream pressure(std::string(argv[4]) + ".pressure");
				if (!output || !pressure) throw std::runtime_error("cannot create Navier-Stokes output");
				output.precision(17); pressure.precision(17);
				for (std::uint64_t node = 0; node < database.header().nodes; ++node) {
					output << PetscRealPart(values[4*node]) << ' ' << PetscRealPart(values[4*node+1]) << ' '
						<< PetscRealPart(values[4*node+2]) << '\n';
					pressure << PetscRealPart(values[4*node+3]) << '\n';
				}
				VecRestoreArrayRead(root, &values);
			}
			VecScatterDestroy(&root_scatter); VecDestroy(&root);
		}
		KSPDestroy(&solver); VecScatterDestroy(&scatter); ISDestroy(&destination_rows); VecDestroy(&ghost_state); ISDestroy(&source_rows);
		VecDestroy(&rhs); VecDestroy(&update); VecDestroy(&state); MatDestroy(&jacobian);
	} catch (const std::exception& e) {
		std::cerr << "rank " << rank << ": " << e.what() << '\n'; status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

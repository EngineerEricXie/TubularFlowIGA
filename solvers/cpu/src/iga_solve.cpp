#include "CaseInput.hpp"
#include "GenericCaseInput.hpp"
#include "GenericTransportElement.hpp"
#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"
#include "SimulationConfig.hpp"
#include "TemporalFunction.hpp"

#include <petscksp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, "TubularFlowIGA configured PDE solver\n");
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc < 3) throw std::runtime_error(
			"usage: iga_solve DATABASE.ntiga CASE_DIR [SYSTEM] [OUTPUT] [VELOCITY]");
		iga::Database database(argv[1]);
		const fs::path case_dir(argv[2]);
		const auto configuration = iga::ReadSimulationConfiguration(
			(case_dir / "simulation_config.json").string());
		const auto system_name = argc >= 4 ? std::string(argv[3])
			: iga::FirstLinearTransportSystem(configuration);
		const auto system = iga::CompileLinearSystem(configuration, system_name);
		const auto labels = iga::ReadPointLabels((case_dir / "controlmesh.vtk").string(), database.header().nodes);
		const auto velocity_path = argc >= 6 ? fs::path(argv[5]) : case_dir / "initial_velocityfield.txt";
		const auto velocity = iga::ReadVelocity(velocity_path.string(), database.header().nodes);
		const auto initial_configuration = iga::MaterializeBoundaryWaveforms(configuration, case_dir, 0.0);
		const auto boundaries = iga::ResolveScalarBoundaries(initial_configuration, system, labels);
		const auto initial = iga::InitialScalarValues(configuration, system);
		const auto fields = system.fields.size();
		if (rank == 0) {
			std::cout << "configuration=simulation_config.json system=" << system.name
				<< " fields=" << fields << " dirichlet_dofs=" << boundaries.constrained_dofs << '\n';
			for (std::size_t i = 0; i < system.fields.size(); ++i)
				std::cout << "field[" << i << "]=" << system.fields[i] << '\n';
		}

		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, fields);
		iga::RequireValidGeometry(assembler.elements(), rank, PETSC_COMM_WORLD);
		std::map<int, long long> surface_faces;
		for (const auto& boundary : configuration.boundaries)
			for (const auto& condition : boundary.conditions)
				if (system.field_index.count(condition.field)
					&& (condition.kind == iga::FieldBoundaryKind::Flux
						|| condition.kind == iga::FieldBoundaryKind::Robin))
					surface_faces.emplace(boundary.label, 0);
		for (const auto& element : assembler.elements()) {
			if (element.owner != rank) continue;
			for (const auto label : element.boundary_labels) {
				auto found = surface_faces.find(label);
				if (found != surface_faces.end()) ++found->second;
			}
		}
		for (auto& entry : surface_faces) {
			long long global_faces = 0;
			MPI_Allreduce(&entry.second, &global_faces, 1, MPI_LONG_LONG, MPI_SUM, PETSC_COMM_WORLD);
			entry.second = global_faces;
			if (global_faces == 0)
				throw std::runtime_error("configured scalar surface boundary label "
					+ std::to_string(entry.first) + " has no boundary faces in the .ntiga database; repack with iga_pack");
		}
		Mat left = assembler.CreateMatrix();
		Mat previous = assembler.CreateMatrix();
		MatSetOption(previous, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
		Vec forcing = assembler.CreateVector();
		VecSet(forcing, 0.0);
		const auto assembly_start = std::chrono::steady_clock::now();
		for (const auto& element : assembler.elements()) {
			const auto matrices = iga::BuildGenericTransportElement(element, velocity, system, configuration);
			assembler.AddElementMatrix(left, element, matrices.left);
			assembler.AddElementMatrix(previous, element, matrices.previous);
			assembler.AddElementVector(forcing, element, matrices.source);
		}
		iga::OwnedRowAssembler::Assemble(left);
		iga::OwnedRowAssembler::Assemble(previous);
		iga::OwnedRowAssembler::Assemble(forcing);

		std::vector<PetscInt> boundary_rows;
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
			for (std::size_t field = 0; field < fields; ++field)
				if (boundaries.constrained[static_cast<std::size_t>(node)*fields+field])
					boundary_rows.push_back(static_cast<PetscInt>(node*fields+field));
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
		for (std::uint64_t node = assembler.node_begin(); node < assembler.node_end(); ++node)
			for (std::size_t field = 0; field < fields; ++field) {
				const auto local = static_cast<std::size_t>(node-assembler.node_begin())*fields+field;
				const auto global = static_cast<std::size_t>(node)*fields+field;
				current_array[local] = boundaries.constrained[global] ? boundaries.value[global] : initial[field];
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
		for (int step = 0; step < system.steps; ++step) {
			const auto step_configuration = iga::MaterializeBoundaryWaveforms(
				configuration, case_dir, (step+1)*system.dt);
			const auto step_boundaries = iga::ResolveScalarBoundaries(step_configuration, system, labels);
			MatMult(previous, current, rhs);
			VecAXPY(rhs, 1.0, forcing);
			std::vector<PetscScalar> boundary_values;
			boundary_values.reserve(boundary_rows.size());
			for (auto row : boundary_rows)
				boundary_values.push_back(step_boundaries.value[static_cast<std::size_t>(row)]);
			VecSetValues(rhs, static_cast<PetscInt>(boundary_rows.size()), boundary_rows.data(), boundary_values.data(), INSERT_VALUES);
			iga::OwnedRowAssembler::Assemble(rhs);
			if (step > 0) VecCopy(current, next);
			KSPSetInitialGuessNonzero(solver, step > 0 ? PETSC_TRUE : PETSC_FALSE);
			KSPSolve(solver, rhs, next);
			KSPConvergedReason reason;
			KSPGetConvergedReason(solver, &reason);
			if (reason <= 0) throw std::runtime_error("KSP did not converge at step " + std::to_string(step));
			PetscInt iterations = 0;
			KSPGetIterationNumber(solver, &iterations);
			total_iterations += iterations;
			VecSwap(current, next);
		}
		const auto solve_end = std::chrono::steady_clock::now();

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
				std::ofstream metadata(std::string(argv[4]) + ".fields");
				if (!output || !metadata) throw std::runtime_error("cannot create configured PDE output");
				output.precision(17);
				for (const auto& name : system.fields) metadata << name << '\n';
				for (std::uint64_t node = 0; node < database.header().nodes; ++node) {
					output << node;
					for (std::size_t field = 0; field < fields; ++field)
						output << ' ' << PetscRealPart(values[node*fields+field]);
					output << '\n';
				}
				VecRestoreArrayRead(root, &values);
			}
			VecScatterDestroy(&scatter);
			VecDestroy(&root);
		}
		PetscReal norm = 0.0;
		VecNorm(current, NORM_2, &norm);
		if (rank == 0) std::cout << "iga_solve system=" << system.name << " steps=" << system.steps
			<< " assembly_s=" << std::chrono::duration<double>(assembly_end-assembly_start).count()
			<< " solve_s=" << std::chrono::duration<double>(solve_end-solve_start).count()
			<< " total_iterations=" << total_iterations << " final_l2=" << norm << '\n';
		KSPDestroy(&solver);
		VecDestroy(&rhs); VecDestroy(&next); VecDestroy(&current); VecDestroy(&forcing);
		MatDestroy(&previous); MatDestroy(&left);
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

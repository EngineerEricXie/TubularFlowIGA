#include "IgaDatabase.hpp"
#include "OwnedRowAssembler.hpp"
#include "TransportBoundaryPreflight.hpp"
#include "TransportElement.hpp"

#include <petscsys.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

void WriteOwnerMismatchDatabase(const fs::path& path)
{
	std::ofstream output(path, std::ios::binary | std::ios::trunc);
	if (!output) throw std::runtime_error("cannot create boundary-preflight database");
	constexpr std::uint64_t header_size = 88;
	const std::uint64_t element_offset = header_size+2*sizeof(std::uint64_t)
		+sizeof(std::int32_t);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, std::uint32_t{2});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{128});
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto rank_index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	for (int axis = 0; axis < 3; ++axis) iga::Write(output, 0.0);
	iga::Write(output, 1.0);
	iga::Write(output, 1.0);
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{1});

	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{1});
	iga::Write(output, std::uint32_t{64});
	const std::array<std::int32_t, 6> boundary_labels{{0, 4, -1, -1, -1, -1}};
	output.write(reinterpret_cast<const char*>(boundary_labels.data()),
		static_cast<std::streamsize>(sizeof(boundary_labels)));
	for (std::int32_t node = 0; node < 64; ++node) iga::Write(output, node);
	for (std::uint8_t row = 0; row < 64; ++row) {
		iga::Write(output, std::uint8_t{1});
		iga::Write(output, row);
		iga::Write(output, 1.0);
	}
	for (int k = 0; k < 4; ++k)
		for (int j = 0; j < 4; ++j)
			for (int i = 0; i < 4; ++i) {
				const std::array<double, 3> point{{i/3.0, j/3.0, k/3.0}};
				output.write(reinterpret_cast<const char*>(point.data()),
					static_cast<std::streamsize>(sizeof(point)));
			}
	const auto element_end = static_cast<std::uint64_t>(output.tellp());
	const auto rank_index_offset = element_end;
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{0});
	output.seekp(header_size+sizeof(std::uint64_t));
	iga::Write(output, element_end);
	output.seekp(rank_index_position);
	iga::Write(output, rank_index_offset);
	if (!output) throw std::runtime_error("failed to finalize boundary-preflight database");
}

iga::SimulationConfiguration MakeConfiguration()
{
	iga::SimulationConfiguration configuration;
	configuration.fields = {{"oxygen", iga::FieldKind::Scalar, 0.0},
		{"drug", iga::FieldKind::Scalar, 0.0}};
	configuration.time = {0.1, 1};
	iga::EquationSystemDefinition system;
	system.name = "surface_transport";
	system.kind = iga::EquationKind::LinearTransport;
	system.unknowns = {"oxygen", "drug"};
	system.terms = {{iga::TermKind::TimeDerivative, "oxygen", "oxygen", 1.0, ""},
		{iga::TermKind::TimeDerivative, "drug", "drug", 1.0, ""}};
	configuration.equation_systems.push_back(system);
	iga::FieldBoundaryCondition flux;
	flux.field = "oxygen";
	flux.kind = iga::FieldBoundaryKind::Flux;
	flux.value = {1.0};
	iga::FieldBoundaryCondition robin;
	robin.field = "drug";
	robin.kind = iga::FieldBoundaryKind::Robin;
	robin.coefficient = 2.0;
	robin.exterior_value = 3.0;
	configuration.boundaries = {{0, "flux", {flux}}, {4, "robin", {robin}}};
	return configuration;
}

std::map<int, long long> CountPackedFaces(iga::Database& database)
{
	std::map<int, long long> result{{0, 0}, {4, 0}};
	for (std::uint64_t index = 0; index < database.header().elements; ++index)
		for (const auto label : database.Load(index).boundary_labels) {
			auto found = result.find(label);
			if (found != result.end()) ++found->second;
		}
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0;
	int size = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	int status = 0;
	const auto path = fs::temp_directory_path()/"tubularflowiga-boundary-preflight.ntiga";
	try {
		Require(size == 2, "boundary-preflight regression requires exactly two MPI ranks");
		if (rank == 0) WriteOwnerMismatchDatabase(path);
		MPI_Barrier(PETSC_COMM_WORLD);
		iga::Database database(path.string());
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 2);
		iga::RequireValidGeometry(assembler.elements(),
			[&assembler](const iga::Element& element) {
				return assembler.OwnsElementByMinimumNode(element);
			}, PETSC_COMM_WORLD);
		Require(database.RequiredElementIndices(0).size() == 1,
			"row owner rank 0 must load the fixture element");
		Require(database.RequiredElementIndices(1).empty(),
			"METIS owner rank 1 must not load the fixture element");
		int local_mismatch = 0;
		for (const auto& element : assembler.elements())
			if (assembler.OwnsElementByMinimumNode(element)) {
				Require(element.owner != rank,
					"fixture did not separate METIS and PETSc row ownership");
				++local_mismatch;
			}
		int global_mismatch = 0;
		MPI_Allreduce(&local_mismatch, &global_mismatch, 1, MPI_INT, MPI_SUM,
			PETSC_COMM_WORLD);
		Require(global_mismatch == 1, "canonical row owner did not count the element once");
		const auto configuration = MakeConfiguration();
		const auto system = iga::CompileLinearSystem(configuration, "surface_transport");
		const auto distributed = iga::CountConfiguredScalarSurfaceFaces(
			configuration, system, assembler, PETSC_COMM_WORLD);
		const auto packed = CountPackedFaces(database);
		Require(distributed == packed,
			"distributed Flux/Robin face counts differ from independent packed count");
		iga::RequireConfiguredScalarSurfaceFaces(
			configuration, system, assembler, PETSC_COMM_WORLD);
		MPI_Barrier(PETSC_COMM_WORLD);
		if (rank == 0) {
			fs::remove(path);
			std::cout << "distributed transport boundary preflight regression passed\n";
		}
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		status = 1;
	}
	int global_status = 0;
	MPI_Allreduce(&status, &global_status, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD);
	PetscFinalize();
	return global_status;
}

#include "ParallelOwnershipValidation.hpp"
#include "OwnedRowAssembler.hpp"
#include "SurfaceOwnershipValidation.hpp"
#include "DynamicWeightedAitkenRelaxation.hpp"
#include "CollectiveFailure.hpp"
#include "OwnedScalarContributions.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unistd.h>
#include <vector>

namespace {

void Require(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

template <class Function>
void Reject(Function function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	Require(rejected, "invalid local ownership input was accepted");
}

template <class Function>
void RejectCollectively(Function function, MPI_Comm communicator)
{
	int rejected = 0, total = 0, ranks = 1;
	try { function(); } catch (const std::runtime_error&) { rejected = 1; }
	MPI_Comm_size(communicator, &ranks);
	MPI_Allreduce(&rejected, &total, 1, MPI_INT, MPI_SUM, communicator);
	Require(total == ranks, "ownership error did not reach every group member");
}

iga::OwnedNodeRange Range(std::uint64_t nodes, int rank, int ranks)
{
	const auto quotient = nodes/static_cast<std::uint64_t>(ranks);
	const auto remainder = nodes%ranks;
	const auto r = static_cast<std::uint64_t>(rank);
	const auto begin = r*quotient+std::min(r, remainder);
	return {begin, begin+quotient+(r < remainder ? 1 : 0), nodes};
}

void CheckLocalMappings()
{
	const iga::OwnedNodeRange range{5, 8, 10};
	range.Validate();
	Require(range.LocalNode(6) == 1 && range.GlobalNode(1) == 6, "global/local node mapping");
	Reject([&] { range.LocalNode(4); });
	Reject([&] { range.GlobalNode(3); });
	Reject([] { iga::OwnedNodeRange{2, 1, 4}.Validate(); });
	Reject([] { iga::OwnedNodeRange{0, 5, 4}.LocalNode(4); });
	const iga::OwnedNodeRange empty{2, 2, 2};
	empty.Validate();
	Reject([&] { empty.GlobalNode(0); });
	constexpr std::uint64_t int32_limit = 2147483647;
	constexpr std::uint64_t int64_limit = 9223372036854775807ULL;
	Require(iga::CheckedFieldRow(536870914, 3, 536870915, 4, int64_limit)
		== 2147483659ULL, "field row overflowed a 32-bit intermediate");
	Reject([&] { iga::CheckedFieldRow(536870914, 3, 536870915, 4, int32_limit); });
	Reject([&] { iga::CheckedFieldRowCount(int64_limit/2+1, 2, int64_limit); });
	Reject([&] { iga::CheckedFieldRowCount(1, 0, int32_limit); });
	Reject([&] { iga::CheckedFieldRow(-1, 0, 10, 4, int32_limit); });
	Reject([&] { iga::CheckedFieldRow(10, 0, 10, 4, int32_limit); });
	Reject([&] { iga::CheckedFieldRow(0, 4, 10, 4, int32_limit); });
	iga::ValidateRequiredNodeIds({0, 3, 7}, 8);
	iga::ValidateRequiredNodeIds({}, 0);
	Reject([] { iga::ValidateRequiredNodeIds({0, 3, 3}, 8); });
	Reject([] { iga::ValidateRequiredNodeIds({3, 0}, 8); });
	Reject([] { iga::ValidateRequiredNodeIds({-1, 0}, 8); });
	Reject([] { iga::ValidateRequiredNodeIds({0, 8}, 8); });
}

void CheckDistributed(MPI_Comm communicator, std::uint64_t nodes)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank);
	MPI_Comm_size(communicator, &ranks);
	const auto range = Range(nodes, rank, ranks);
	iga::ValidateDistributedNodeRanges(range, communicator);
	std::vector<std::uint64_t> owned;
	for (auto node = range.begin; node < range.end; ++node) owned.push_back(node);
	iga::ValidateUniqueEntityCoverage(owned, nodes, communicator, "valid nodes");
	iga::ValidateUniqueEntityCoverage({}, 0, communicator, "empty catalog");
	auto bad = owned;
	if (rank == 0) bad.push_back(nodes);
	RejectCollectively([&] { iga::ValidateUniqueEntityCoverage(bad, nodes,
		communicator, "out-of-range ID"); }, communicator);
	bad = owned;
	if (rank == 0 && !bad.empty()) bad.pop_back();
	if (nodes) RejectCollectively([&] { iga::ValidateUniqueEntityCoverage(bad, nodes,
		communicator, "missing ID"); }, communicator);
	if (ranks > 1) {
		RejectCollectively([&] { iga::ValidateUniqueEntityCoverage(owned,
			nodes+(rank == 0 ? 1 : 0), communicator, "inconsistent global count"); }, communicator);
		auto overlap = range;
		if (rank == 1) overlap.begin = 0;
		if (nodes > 1) RejectCollectively([&] { iga::ValidateDistributedNodeRanges(overlap,
			communicator); }, communicator);
	}
}

void WriteAuditDatabase(const std::filesystem::path& path, bool omit_rank_one)
{
	// Algebraic fixture: ranks 0/1 own the two rows; rank 2 owns the
	// element but no rows and has an empty required-element catalog.
	std::ofstream output(path, std::ios::binary);
	constexpr std::uint64_t header_size = 88;
	constexpr std::uint64_t element_offset = header_size+2*sizeof(std::uint64_t)+sizeof(std::int32_t);
	output.write(iga::kMagic.data(), iga::kMagic.size());
	iga::Write(output, iga::kVersion);
	iga::Write(output, std::uint32_t{3});
	iga::Write(output, std::uint64_t{1});
	iga::Write(output, std::uint64_t{2});
	iga::Write(output, iga::kBezierPointCount);
	iga::Write(output, std::uint32_t{0});
	const auto index_position = output.tellp();
	iga::Write(output, std::uint64_t{0});
	for (int axis = 0; axis < 3; ++axis) iga::Write(output, 0.0);
	iga::Write(output, 1.0); iga::Write(output, 1.0);
	iga::Write(output, element_offset);
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{2});
	iga::Write(output, std::uint64_t{0});
	iga::Write(output, std::int32_t{0});
	iga::Write(output, std::int32_t{2});
	iga::Write(output, std::uint32_t{2});
	for (int face = 0; face < 6; ++face) iga::Write(output, std::int32_t{-1});
	iga::Write(output, std::int32_t{0}); iga::Write(output, std::int32_t{1});
	for (std::uint8_t column : {0, 63}) {
		iga::Write(output, std::uint8_t{1});
		iga::Write(output, column); iga::Write(output, 1.0);
	}
	for (int k = 0; k < 4; ++k) for (int j = 0; j < 4; ++j) for (int i = 0; i < 4; ++i) {
		iga::Write(output, i/3.0); iga::Write(output, j/3.0); iga::Write(output, k/3.0);
	}
	const auto index_offset = static_cast<std::uint64_t>(output.tellp());
	for (const auto offset : {std::uint64_t{0}, std::uint64_t{1},
		std::uint64_t{omit_rank_one ? 1U : 2U}, std::uint64_t{omit_rank_one ? 1U : 2U}})
		iga::Write(output, offset);
	iga::Write(output, std::uint64_t{0});
	if (!omit_rank_one) iga::Write(output, std::uint64_t{0});
	output.seekp(header_size+sizeof(std::uint64_t));
	iga::Write(output, index_offset);
	output.seekp(index_position);
	iga::Write(output, index_offset);
}

void CheckIndependentCatalogs(int rank)
{
	int root_pid = rank == 0 ? static_cast<int>(getpid()) : 0;
	MPI_Bcast(&root_pid, 1, MPI_INT, 0, PETSC_COMM_WORLD);
	const auto directory = std::filesystem::temp_directory_path()
		/("tubularflow-ownership-catalogs-"+std::to_string(root_pid));
	if (rank == 0) {
		Require(std::filesystem::create_directory(directory), "fixture directory exists");
		WriteAuditDatabase(directory/"valid.ntiga", false);
		WriteAuditDatabase(directory/"missing-row.ntiga", true);
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	{
		iga::Database database((directory/"valid.ntiga").string());
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 1);
		assembler.ValidateOwnership();
		if (rank == 2) Require(assembler.local_rows() == 0 && assembler.elements().empty(),
			"fixture did not exercise an element owner with no owned rows");
		Mat matrix = assembler.CreateMatrix();
		for (const auto& element : assembler.elements())
			assembler.AddElementMatrix(matrix, element, std::vector<PetscScalar>{2.0, -1.0, -1.0, 2.0});
		iga::OwnedRowAssembler::Assemble(matrix);
		for (auto node = assembler.node_begin(); node < assembler.node_end(); ++node) {
			const auto row = static_cast<PetscInt>(node);
			const PetscInt columns[2] = {0, 1};
			PetscScalar values[2]{};
			Require(MatGetValues(matrix, 1, &row, 2, columns, values) == 0, "cannot read owned matrix row");
			for (int column = 0; column < 2; ++column)
				Require(PetscRealPart(values[column]) == (column == row ? 2.0 : -1.0),
					"shared element matrix contribution was omitted or duplicated");
		}
		MatDestroy(&matrix);
	}
	{
		iga::Database database((directory/"missing-row.ntiga").string());
		iga::OwnedRowAssembler assembler(database, PETSC_COMM_WORLD, 1);
		RejectCollectively([&] { assembler.ValidateOwnership(); }, PETSC_COMM_WORLD);
	}
	if (rank == 0) std::cout << "independent partition/row catalogs and missing-row rejection passed fixture="
		<< directory << '\n';
}

void CheckScalarContributions(int rank)
{
	const auto id = std::numeric_limits<std::uint64_t>::max()-1;
	const std::vector<std::uint64_t> owned = rank == 1 ? std::vector<std::uint64_t>{id, 0} : std::vector<std::uint64_t>{};
	const std::vector<std::pair<std::uint64_t,double>> values{{id, static_cast<double>(rank+1)}, {id, -.5}};
	const auto expected = rank == 1 ? std::vector<double>{4.5, 0.} : std::vector<double>{};
	Require(iga::SumOwnedScalarContributions(PETSC_COMM_WORLD, owned, values) == expected, "wrong owner scalar sum");
	Require(iga::SumOwnedScalarContributions(PETSC_COMM_WORLD, {}, {}).empty(), "nonempty all-empty scalar sum");
	for(int mode=0;mode<6;++mode) {
		auto ids=owned;auto data=values;iga::PointIdentityLimits limits;
		if(rank==2) {
			if(mode==0)ids.push_back(id);
			if(mode==1)data.push_back({17, 1.});
			if(mode==2)data.push_back({id, std::numeric_limits<double>::infinity()});
			if(mode==3)limits.max_wire_bytes=8;
			if(mode==4)limits.max_local_occurrences=1;
			if(mode==5) {data.push_back({id, std::numeric_limits<double>::max()});data.push_back({id, std::numeric_limits<double>::max()});}
		}
		RejectCollectively([&] { (void)iga::SumOwnedScalarContributions(PETSC_COMM_WORLD, ids, data, limits); }, PETSC_COMM_WORLD);
		Require(iga::SumOwnedScalarContributions(PETSC_COMM_WORLD, owned, values) == expected, "scalar retry failed");
	}
	if(rank==0)std::cout << "owned_scalar_contributions=passed failures=6 retries=6\n";
}

void CheckEmptyAitken(int rank)
{
	const std::string identity(64, static_cast<char>('a'+rank));
	const std::vector<double> weights = rank == 1 ? std::vector<double>{1., 3.} : std::vector<double>{};
	iga::DynamicWeightedAitkenRelaxation distributed(identity, weights, 4.);
	iga::DynamicWeightedAitkenRelaxation serial(std::string(64, 'f'), {1., 3.}, 4.);
	for (const auto& full : {std::vector<double>{2., -1.}, std::vector<double>{1., 2.}}) {
		const auto local = rank == 1 ? full : std::vector<double>{};
		const auto control = distributed.ControlStateIdentitySha256();
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "empty Aitken control", control);
		const auto terms = distributed.LocalContributions(local, identity);
		double contributions[2]{terms.numerator, terms.denominator}, global[2]{};
		MPI_Allreduce(contributions, global, 2, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
		const double scale = distributed.LocalRequiredResidualScale(local, 1., identity);
		double global_scale = 0.;
		MPI_Allreduce(&scale, &global_scale, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
		const auto proposal = distributed.ProposeWithGlobalReduction(std::vector<double>(local.size(), 0.),
			local, 1., global_scale, global[0], global[1], control, identity);
		const auto reference = serial.Propose({0., 0.}, full, 1., std::string(64, 'f'));
		Require(proposal.relaxation_factor == reference.relaxation_factor, "empty Aitken changed global relaxation");
		if (rank != 1) Require(proposal.next.empty(), "empty Aitken published values");
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "empty Aitken pending", distributed.ControlStateIdentitySha256());
		distributed.AcceptApplied(proposal, local, proposal.relaxation_factor, identity);
		serial.AcceptApplied(reference, full, reference.relaxation_factor, std::string(64, 'f'));
		iga::RequireCollectiveSameText(PETSC_COMM_WORLD, "empty Aitken accepted", distributed.ControlStateIdentitySha256());
	}
	if (rank == 0) std::cout << "empty_aitken_reduction=passed iterations=2 empty_ranks=2\n";
}

void CheckSurfaceOwnership(int rank)
{
	iga::SurfaceInterfaceRef reference{"fluid", "flow", "wall"};
	iga::DistributedSurfaceLayout layout;
	layout.reference_mesh_identity_sha256 = std::string(64, 'a');
	layout.global_node_count = 3;
	layout.partition_count = 3;
	layout.partition_rank = static_cast<std::uint64_t>(rank);
	layout.owned_global_node_ids = {10*static_cast<std::uint64_t>(rank+1)};
	layout.reference_positions = {{10, {{0, 0, 0}}}, {20, {{1, 0, 0}}}, {30, {{0, 1, 0}}}};
	layout.reference_triangles = {{{10, 20, 30}}};
	layout.owned_reference_lumped_areas_m2 = {1.0/6.0};
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	iga::ValidateSurfacePublicationOwnership(reference, layout, PETSC_COMM_WORLD);
	// Both the catalog root and another peer may publish no nodes.
	auto sparse = layout;
	sparse.owned_global_node_ids = rank == 1 ? std::vector<std::uint64_t>{10, 20, 30} : std::vector<std::uint64_t>{};
	sparse.owned_reference_lumped_areas_m2.assign(sparse.owned_global_node_ids.size(), 1.0/6.0);
	iga::ValidateSurfacePublicationOwnership(reference, sparse, PETSC_COMM_WORLD);
	auto missing = sparse;
	if (rank == 1) {
		missing.owned_global_node_ids.pop_back();
		missing.owned_reference_lumped_areas_m2.pop_back();
	}
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(reference, missing,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
	iga::ValidateSurfacePublicationOwnership(reference, sparse, PETSC_COMM_WORLD);

	auto duplicate = layout;
	if (rank == 2) duplicate.owned_global_node_ids = {20};
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(reference, duplicate,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
	auto wrong_rank = layout;
	if (rank == 1) wrong_rank.partition_rank = 0;
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(reference, wrong_rank,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
	auto geometry = layout;
	if (rank == 1) {
		geometry.reference_positions[1].position_m[0] = 2.0;
		geometry.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(geometry);
	}
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(reference, geometry,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
	auto other_interface = reference;
	if (rank == 1) other_interface.interface_id = "another_wall";
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(other_interface, layout,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
	auto unknown = layout;
	if (rank == 1) unknown.owned_global_node_ids = {999};
	RejectCollectively([&] { iga::ValidateSurfacePublicationOwnership(reference, unknown,
		PETSC_COMM_WORLD); }, PETSC_COMM_WORLD);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	try {
		Require(ranks == 3, "parallel_ownership_test requires three ranks");
		CheckLocalMappings();
		CheckDistributed(PETSC_COMM_WORLD, 7);
		CheckDistributed(PETSC_COMM_WORLD, 2); // Rank 2 owns no nodes.
		CheckDistributed(PETSC_COMM_WORLD, 0);
		CheckIndependentCatalogs(rank);
		CheckSurfaceOwnership(rank);
		CheckEmptyAitken(rank);
		CheckScalarContributions(rank);
		std::vector<std::uint64_t> balanced;
		for (std::uint64_t id = rank; id < 6; id += ranks) balanced.push_back(id);
		// Replace 1,4 with 2,3: global count AND ID sum remain unchanged.
		// Exact coverage must still detect the paired duplicates and omissions.
		if (rank == 1) balanced = {2, 3};
		RejectCollectively([&] { iga::ValidateUniqueEntityCoverage(balanced, 6,
			PETSC_COMM_WORLD, "duplicate/omission with equal count and sum"); }, PETSC_COMM_WORLD);
		std::vector<iga::OwnershipIncidence> expected, actual;
		if (rank == 0) expected = {{0, 0}, {0, 1}, {1, 1}, {1, 2}};
		if (rank == 0) actual = {{0, 0}};
		if (rank == 1) actual = {{0, 1}, {1, 1}};
		if (rank == 2) actual = {{1, 2}};
		iga::ValidateOwnershipIncidences(expected, actual, PETSC_COMM_WORLD, "shared node in two elements");
		auto missing = actual;
		if (rank == 1) missing.pop_back();
		RejectCollectively([&] { iga::ValidateOwnershipIncidences(expected, missing,
			PETSC_COMM_WORLD, "missing remote row contribution"); }, PETSC_COMM_WORLD);
		auto duplicate = actual;
		if (rank == 2) duplicate.push_back({1, 2});
		RejectCollectively([&] { iga::ValidateOwnershipIncidences(expected, duplicate,
			PETSC_COMM_WORLD, "duplicate row contribution"); }, PETSC_COMM_WORLD);
		auto incorrect = actual;
		if (rank == 1) incorrect[1] = {1, 0};
		RejectCollectively([&] { iga::ValidateOwnershipIncidences(expected, incorrect,
			PETSC_COMM_WORLD, "wrong node with unchanged contribution count"); }, PETSC_COMM_WORLD);
		MPI_Comm group = MPI_COMM_NULL;
		MPI_Comm_split(PETSC_COMM_WORLD, rank == 0 ? 0 : 1, rank, &group);
		CheckDistributed(group, rank == 0 ? 1 : 5);
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "parallel ownership mappings, exact coverage, empty ranks and collective rejection passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

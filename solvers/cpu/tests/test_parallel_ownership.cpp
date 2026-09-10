#include "ParallelOwnershipValidation.hpp"
#include "OwnedRowAssembler.hpp"
#include "SurfaceOwnershipValidation.hpp"
#include "DynamicWeightedAitkenRelaxation.hpp"
#include "CollectiveFailure.hpp"
#include "OwnedScalarContributions.hpp"
#include "DistributedSurfaceAreas.hpp"
#include "OwnedPointValues.hpp"
#include "SurfaceGhostKinematics.hpp"
#include "DistributedSurfaceTriangleKinematics.hpp"

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

void CheckOwnedPointValues(int rank)
{
	const auto id=std::numeric_limits<std::uint64_t>::max();
	const auto owned=rank==1?std::vector<std::uint64_t>{id,0}:std::vector<std::uint64_t>{};
	const auto values=rank==1?std::vector<double>{1.,-0.,3.,4.,5.,6.}:std::vector<double>{};
	const auto queries=rank==2?std::vector<std::uint64_t>{}:std::vector<std::uint64_t>{0,id,0};
	const auto expected=rank==2?std::vector<double>{}:std::vector<double>{4.,5.,6.,1.,-0.,3.,4.,5.,6.};
	const auto check=[&] {
		const auto result=iga::FetchOwnedPointValues(PETSC_COMM_WORLD,owned,values,queries,3);
		Require(result==expected,"wrong fetched owner tuples");
		if(!result.empty())Require(std::signbit(result[4]),"point exchange lost signed zero");
	};
	check();
	Require(iga::FetchOwnedPointValues(PETSC_COMM_WORLD,{},{},{},3).empty(),"all-empty point fetch failed");
	for(int mode=0;mode<6;++mode) {
		auto ids=owned;auto data=values;auto requested=queries;int components=3;iga::PointIdentityLimits limits;
		if(rank==2) {
			if(mode==0){ids.push_back(id);data={7.,8.,9.};}
			if(mode==1)requested.push_back(17);
			if(mode==2){ids.push_back(18);data={1.,std::numeric_limits<double>::infinity(),3.};}
			if(mode==3)components=2;
			if(mode==4){requested.push_back(id);limits.max_wire_bytes=8;}
			if(mode==5){requested={id,0};limits.max_local_occurrences=1;}
		}
		RejectCollectively([&] {(void)iga::FetchOwnedPointValues(PETSC_COMM_WORLD,ids,data,requested,components,limits);},PETSC_COMM_WORLD);
		check();
	}
	if(rank==0)std::cout << "owned_point_values=passed failures=6 retries=6\n";
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

	iga::SurfaceKinematics publication;
	publication.interface=reference;
	publication.stamp.time_s=.25;publication.stamp.step=4;publication.stamp.coupling_iteration=2;
	publication.stamp.reference_mesh_identity_sha256=sparse.reference_mesh_identity_sha256;
	publication.stamp.layout_identity_sha256=sparse.layout_identity_sha256;
	publication.stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(sparse);
	publication.stamp.producer_state_identity_sha256=std::string(64,static_cast<char>('a'+rank));
	for(auto id:sparse.owned_global_node_ids) {
		publication.displacement_m.push_back({{static_cast<double>(id),0.,0.}});
		publication.velocity_m_per_s.push_back({{0.,static_cast<double>(id),0.}});
	}
	const auto expected_stamp=publication.stamp;
	const auto requests=rank==0?std::vector<std::uint64_t>{30,10}:
		(rank==1?std::vector<std::uint64_t>{20}:std::vector<std::uint64_t>{});
	const auto check_ghosts=[&] {
		const auto ghosts=iga::FetchSurfaceGhostKinematics(PETSC_COMM_WORLD,sparse,publication,reference,expected_stamp,requests);
		Require(ghosts.requested_node_ids==requests,"ghost IDs reordered");
		for(std::size_t row=0;row<requests.size();++row) {
			Require(ghosts.displacement_m[row]==std::array<double,3>{{static_cast<double>(requests[row]),0.,0.}},"wrong ghost displacement");
			Require(ghosts.velocity_m_per_s[row]==std::array<double,3>{{0.,static_cast<double>(requests[row]),0.}},"wrong ghost velocity");
		}
	};
	check_ghosts();
	for(int mode=0;mode<7;++mode) {
		auto incoming=publication;auto expected=expected_stamp;auto query=requests;
		if(rank==2) {
			if(mode==0)incoming.stamp.time_s=.5;
			if(mode==1)incoming.stamp.producer_state_identity_sha256=std::string(64,'f');
			if(mode==2)incoming.interface.subsystem_id="wrong";
			if(mode==3)incoming.stamp.partition_identity_sha256=std::string(64,'f');
			if(mode==4){incoming.stamp.coupling_iteration=3;expected.coupling_iteration=3;}
			if(mode==5)query.push_back(999);
		}
		if(mode==6&&rank==1)incoming.displacement_m[0][0]=std::numeric_limits<double>::infinity();
		RejectCollectively([&] {(void)iga::FetchSurfaceGhostKinematics(PETSC_COMM_WORLD,sparse,incoming,reference,expected,query);},PETSC_COMM_WORLD);
		check_ghosts();
	}
	if(rank==0)std::cout << "surface_ghost_kinematics=passed failures=7 retries=7\n";

	// Reference triangle owner need not own any of its material nodes.
	const std::vector<std::uint64_t> triangles = rank == 2 ? std::vector<std::uint64_t>{0} : std::vector<std::uint64_t>{};
	const auto check_triangles=[&] {
		const auto local=iga::BuildDistributedSurfaceTriangleKinematics(PETSC_COMM_WORLD,sparse,publication,reference,expected_stamp,triangles);
		Require(local.size()==triangles.size(),"wrong local triangle count");
		if(rank==2) {
			Require(local[0].reference_triangle_index==0 && local[0].node_ids==std::array<std::uint64_t,3>{{10,20,30}},"triangle identity changed");
			Require(local[0].positions_m==std::array<std::array<double,3>,3>{{{{10,0,0}},{{21,0,0}},{{30,1,0}}}},"wrong deformed triangle positions");
			Require(local[0].velocities_m_per_s==std::array<std::array<double,3>,3>{{{{0,10,0}},{{0,20,0}},{{0,30,0}}}},"wrong triangle velocities");
		}
	};
	check_triangles();
	for(int mode=0;mode<4;++mode) {
		auto selected=triangles;auto incoming=publication;iga::PointIdentityLimits limits;
		if(mode==0 && rank==0)selected.push_back(0);
		if(mode==1)selected.clear();
		if(mode==2 && rank==2)limits.max_local_occurrences=2;
		if(mode==3 && rank==0)incoming.stamp.time_s=.5;
		RejectCollectively([&] {(void)iga::BuildDistributedSurfaceTriangleKinematics(PETSC_COMM_WORLD,sparse,incoming,reference,expected_stamp,selected,limits);},PETSC_COMM_WORLD);
		check_triangles();
	}
	if(rank==0)std::cout << "surface_triangle_kinematics=passed failures=4 retries=4\n";
	for(const auto& partition : {layout, sparse}) {
		const auto before=iga::BuildDistributedSurfacePartitionIdentitySha256(partition);
		const auto areas=iga::ComputeDistributedSurfaceAreas(PETSC_COMM_WORLD, reference, partition, triangles);
		Require(areas.global_area_m2==.5, "wrong reference surface area");
		Require(areas.owned_lumped_areas_m2==std::vector<double>(partition.owned_global_node_ids.size(),1./6.), "wrong nodal areas");
		Require(iga::BuildDistributedSurfacePartitionIdentitySha256(partition)==before, "area calculation mutated partition");
	}
	for(int mode=0;mode<3;++mode) {
		auto selected=triangles;iga::PointIdentityLimits limits;
		if(mode==0&&rank==0)selected.push_back(0);
		if(mode==1)selected.clear();
		if(mode==2&&rank==2)limits.max_local_occurrences=2;
		RejectCollectively([&] { (void)iga::ComputeDistributedSurfaceAreas(PETSC_COMM_WORLD, reference, sparse, selected, limits); },PETSC_COMM_WORLD);
		Require(iga::ComputeDistributedSurfaceAreas(PETSC_COMM_WORLD, reference, sparse, triangles).global_area_m2==.5,"surface area retry failed");
	}
	auto patch=layout;
	patch.global_node_count=4;
	patch.reference_positions.push_back({40,{{2.,1.,0.}}});
	patch.reference_triangles.push_back({{20,40,30}});
	patch.owned_global_node_ids=rank==0?std::vector<std::uint64_t>{10}:
		(rank==1?std::vector<std::uint64_t>{20,30}:std::vector<std::uint64_t>{40});
	patch.owned_reference_lumped_areas_m2.assign(patch.owned_global_node_ids.size(),1.);
	patch.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(patch);
	const auto patch_triangles=rank==0?std::vector<std::uint64_t>{1}:
		(rank==2?std::vector<std::uint64_t>{0}:std::vector<std::uint64_t>{});
	const auto patch_areas=iga::ComputeDistributedSurfaceAreas(PETSC_COMM_WORLD,reference,patch,patch_triangles);
	Require(patch_areas.global_area_m2==1.5,"wrong unequal triangle area total");
	const auto expected_areas=rank==0?std::vector<double>{1./6.}:
		(rank==1?std::vector<double>{.5,.5}:std::vector<double>{1./3.});
	Require(patch_areas.owned_lumped_areas_m2==expected_areas,"shared node triangle areas were not accumulated");
	if(rank==0)std::cout << "distributed_surface_areas=passed failures=3 retries=3\n";
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
		CheckOwnedPointValues(rank);
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

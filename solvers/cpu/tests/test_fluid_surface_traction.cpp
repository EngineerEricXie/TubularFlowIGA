#include "FluidSurfaceTraction.hpp"
#include "OwnedFluidSurfaceTractionPoints.hpp"
#include "DistributedFluidSurfaceTraction.hpp"
#include "PrescribedSurfaceMotion.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kTolerance = 2.0e-9;

bool Near(double left, double right, double tolerance = kTolerance)
{ return std::abs(left-right) <= tolerance*std::max({1.0, std::abs(left), std::abs(right)}); }

void CheckVector(const std::array<double, 3>& actual, const std::array<double, 3>& expected,
	double tolerance = kTolerance)
{ for (int component = 0; component < 3; ++component) assert(Near(actual[component], expected[component], tolerance)); }

template <class Function> void Reject(Function&& function)
{ bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

template <class Function> void RejectContaining(const char* expected, Function&& function)
{
	bool rejected = false;
	try { function(); }
	catch (const std::exception& error) {
		rejected = true;
		assert(std::string(error.what()).find(expected) != std::string::npos);
	}
	assert(rejected);
}

void CheckP1Contribution()
{
	// Degree-two triangle rule integrates P1 traction and consistent mass.
	std::array<std::array<double,3>,3> mass{},force{};
	const std::array<std::array<double,3>,3> nodal_traction{{{{1,2,3}},{{4,5,6}},{{7,8,9}}}};
	for(int point=0;point<3;++point) {
		std::array<double,3> shape{{1./6.,1./6.,1./6.}},traction{};shape[point]=2./3.;
		for(int node=0;node<3;++node)for(int axis=0;axis<3;++axis)traction[axis]+=shape[node]*nodal_traction[node][axis];
		const auto value=iga::BuildSurfaceP1TractionContribution(shape,traction,1./6.);
		for(int left=0;left<3;++left)for(int right=0;right<3;++right) {
			mass[left][right]+=value.consistent_mass_m2[left][right];
			force[left][right]+=value.corner_force_n[left][right];
		}
	}
	for(int node=0;node<3;++node)for(int axis=0;axis<3;++axis) {
		assert(Near(mass[node][axis],node==axis?1./12.:1./24.,1.e-14));
		double exact=0.;for(int other=0;other<3;++other)exact+=(node==other?1./12.:1./24.)*nodal_traction[other][axis];
		assert(Near(force[node][axis],exact,1.e-14));
	}
	const double inf=std::numeric_limits<double>::infinity(),large=std::numeric_limits<double>::max();
	Reject([&]{iga::BuildSurfaceP1TractionContribution({{1,0,0}},{{1,2,3}},inf);});
	Reject([&]{iga::BuildSurfaceP1TractionContribution({{1,inf,0}},{{1,2,3}},1);});
	Reject([&]{iga::BuildSurfaceP1TractionContribution({{1,0,0}},{{1,2,inf}},1);});
	Reject([&]{iga::BuildSurfaceP1TractionContribution({{1,0,0}},{{large,0,0}},2);});
	Reject([&]{iga::BuildSurfaceP1TractionContribution({{large,0,0}},{{0,0,0}},1);});
}

iga::RawSurfaceTriangle Face(int a, int b, int c, std::uint32_t label)
{ iga::RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = label; return result; }

iga::MaterialSurfaceKinematics Material(double time=.5)
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	// A-C-B has outward -z; this is the selected planar material patch.
	soup.triangles = {Face(0,2,1,7), Face(0,1,3,9), Face(0,3,2,9), Face(1,2,3,9)};
	iga::PrescribedSurfaceFrame first{0.0,soup}, second{1.0,soup};
	return iga::PrescribedSurfaceMotion({first,second}).Evaluate(time,0.0,time);
}

iga::DistributedSurfaceInterface StructuralInterface(const std::string& reference_identity)
{
	iga::DistributedSurfaceInterface result;
	result.id = {"structure", "membrane", "patch"}; result.subsystem_id = "membrane";
	result.boundary_labels = {7}; result.reference_mesh_identity_sha256 = reference_identity;
	result.provides = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	result.requires = {iga::SurfaceFieldQuantity::TractionOnStructure};
	return result;
}

iga::DistributedSurfaceLayout Layout(const iga::MaterialSurfaceKinematics& material)
{
	iga::DistributedSurfaceLayout result;
	result.reference_mesh_identity_sha256 = std::string(64, '0'); result.global_node_count = 3;
	result.partition_count = 1; result.partition_rank = 0; result.owned_global_node_ids = {11,20,30};
	// The source patch is A-C-B; layout IDs intentionally do not coincide with source IDs.
	result.reference_positions = {{11,material.ReferenceMaterialVerticesM()[0]}, {20,material.ReferenceMaterialVerticesM()[1]},
		{30,material.ReferenceMaterialVerticesM()[2]}};
	result.reference_triangles = {{{{11,30,20}}}};
	result.owned_reference_lumped_areas_m2 = {1.0/6.0,1.0/6.0,1.0/6.0};
	result.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(result);
	return result;
}

iga::MaterialSurfacePatchMap PatchMap(const iga::MaterialSurfaceKinematics& material)
{
	auto layout = Layout(material);
	const std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> vertices{{11,0},{20,1},{30,2}};
	const std::vector<std::uint32_t> triangles{0}; const std::vector<std::uint64_t> clamps{11,20,30};
	const auto reference = iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(material, 7, vertices,
		layout.reference_triangles, triangles, clamps);
	layout.reference_mesh_identity_sha256 = reference;
	layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	return iga::MaterialSurfacePatchMap::Create(StructuralInterface(reference), std::move(layout), material, 7,
		vertices, triangles, clamps);
}

iga::DistributedSurfaceInterface Interface(const iga::MaterialSurfacePatchMap& patch_map)
{
	iga::DistributedSurfaceInterface result;
	result.id = {"fluid", "immersed", "patch"}; result.subsystem_id = "immersed";
	result.boundary_labels = {7}; result.reference_mesh_identity_sha256 = patch_map.ReferenceIdentitySha256();
	result.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
	result.requires = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::Velocity};
	return result;
}

iga::CartesianDomainClassification Domain(const iga::MaterialSurfaceKinematics& material)
{
	const iga::CubicCartesianGridSpec grid{{{-.1,-.1,-.1}},{{1.1,1.1,1.1}},{{2,2,2}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(grid),
		iga::SurfaceSpatialIndex(material.Surface()));
}

iga::MaterialSurfaceKinematics TranslatedMaterial(double x_offset)
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{x_offset,0,0}},{{x_offset+1,0,0}},{{x_offset,1,0}},{{x_offset,0,1}}}};
	soup.triangles = {Face(0,2,1,7), Face(0,1,3,9), Face(0,3,2,9), Face(1,2,3,9)};
	return iga::PrescribedSurfaceMotion({{0.0,soup},{1.0,soup}}).Evaluate(.5,0.0,.5);
}

iga::CartesianDomainClassification TranslatedDomain(const iga::MaterialSurfaceKinematics& material, double x_offset)
{
	const iga::CubicCartesianGridSpec grid{{{x_offset-.1,-.1,-.1}},{{x_offset+1.1,1.1,1.1}},{{2,2,2}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(grid), iga::SurfaceSpatialIndex(material.Surface()));
}

std::vector<iga::FluidSurfaceElementState> State(const iga::CartesianDomainClassification& domain,
	const iga::ImmersedSurfaceQuadratureCatalog& catalog, double pressure,
	const std::array<std::array<double, 3>, 3>& gradient,
	const std::vector<std::uint32_t>& labels = {7})
{
	std::vector<iga::FluidSurfaceElementState> result;
	for (std::uint64_t id = 0; id < domain.Cells().size(); ++id) {
		if (domain.Cells()[static_cast<std::size_t>(id)].classification != iga::CellClassification::Cut) continue;
		const auto& rule = catalog.UsableRule(domain, id);
		if (!std::any_of(rule.Points().begin(), rule.Points().end(), [&labels](const iga::SurfaceQuadraturePoint& point) {
			return std::find(labels.begin(), labels.end(), point.boundary_id) != labels.end();
		})) continue;
		const auto element = domain.Background().MaterializeElement(id);
		iga::FluidSurfaceElementState item; item.cell_id = id; item.nodal_state.resize(element.connectivity.size());
		const auto basis_x = domain.Background().Spec().cells[0]+3, basis_y = domain.Background().Spec().cells[1]+3;
		for (std::size_t node = 0; node < element.connectivity.size(); ++node) {
			const auto global = static_cast<std::uint32_t>(element.connectivity[node]);
			const std::uint32_t i = global%basis_x, j = (global/basis_x)%basis_y, k = global/(basis_x*basis_y);
			const auto x = domain.Background().Greville(i,j,k);
			for (int velocity = 0; velocity < 3; ++velocity)
				for (int direction = 0; direction < 3; ++direction)
					item.nodal_state[node][velocity] += gradient[velocity][direction]*x[direction];
			item.nodal_state[node][3] = pressure;
		}
		result.push_back(std::move(item));
	}
	return result;
}

iga::SurfaceFieldStamp Stamp(const iga::MaterialSurfaceKinematics& material,
	const iga::CartesianDomainClassification& domain, const iga::ImmersedSurfaceQuadratureCatalog& catalog,
	const iga::DistributedSurfaceLayout& layout, const iga::MaterialSurfacePatchMap& patch_map,
	const std::vector<iga::FluidSurfaceElementState>& state, double viscosity)
{
	iga::SurfaceFieldStamp result; result.time_s = material.EvaluatedTimeS(); result.step = 4; result.coupling_iteration = 2;
	result.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256; result.layout_identity_sha256 = layout.layout_identity_sha256;
	result.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	result.producer_state_identity_sha256 = iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, material, patch_map, viscosity, state);
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		CheckP1Contribution();
		const auto material = Material(); const auto patch_map = PatchMap(material);
		const auto interface = Interface(patch_map); const auto& layout = patch_map.Layout();
		assert(!(interface.id == patch_map.Interface().id));
		const auto domain = Domain(material); const iga::ImmersedSurfaceQuadratureCatalog catalog(domain);
		const std::array<std::array<double,3>,3> zero{}; constexpr double viscosity = 2.0;
		const auto pressure_state = State(domain, catalog, 20.0, zero); const auto pressure_stamp = Stamp(material, domain, catalog, layout, patch_map, pressure_state, viscosity);
		const auto pressure = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, patch_map, pressure_stamp, viscosity, pressure_state);
		// Selected face has n_f=-z, so -sigma n_f=p n_f=(0,0,-20) Pa.
		for (const auto& traction : pressure.traction.traction_on_structure_pa) CheckVector(traction, {{0,0,-20}});
		CheckVector(pressure.diagnostics.quadrature_resultant_n, {{0,0,-10}});
		CheckVector(pressure.diagnostics.nodal_resultant_n, {{0,0,-10}});
		assert(pressure.diagnostics.retained_quadrature_points > 0);

		const std::array<std::array<double,3>,3> affine{{{{0,0,1}},{{0,0,0}},{{2,0,3}}}};
		const auto affine_state = State(domain, catalog, 4.0, affine); const auto affine_stamp = Stamp(material, domain, catalog, layout, patch_map, affine_state, viscosity);
		const auto affine_result = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, patch_map, affine_stamp, viscosity, affine_state);
		// p n - mu(grad u + grad u^T)n = (6,0,8) Pa on this patch.
		for (const auto& traction : affine_result.traction.traction_on_structure_pa) CheckVector(traction, {{6,0,8}});
		CheckVector(affine_result.diagnostics.quadrature_resultant_n, {{3,0,4}});
		CheckVector(affine_result.diagnostics.nodal_resultant_n, {{3,0,4}});
		CheckVector(affine_result.diagnostics.quadrature_moment_n_m, affine_result.diagnostics.nodal_moment_n_m);

		int rank=0,ranks=1;MPI_Comm_rank(PETSC_COMM_WORLD,&rank);MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
		auto distributed_layout=layout;distributed_layout.partition_count=ranks;distributed_layout.partition_rank=rank;
		distributed_layout.owned_global_node_ids.clear();distributed_layout.owned_reference_lumped_areas_m2.clear();
		for(std::size_t i=0;i<layout.owned_global_node_ids.size();++i)if(static_cast<int>(i%ranks)==rank) {
			distributed_layout.owned_global_node_ids.push_back(layout.owned_global_node_ids[i]);
			distributed_layout.owned_reference_lumped_areas_m2.push_back(layout.owned_reference_lumped_areas_m2[i]);
		}
		distributed_layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(distributed_layout);
		std::vector<std::uint64_t> owned_cells;
		for(std::uint64_t id=0;id<domain.Cells().size();++id)if(static_cast<int>(id%ranks)==rank)owned_cells.push_back(id);
		for(const auto* all_state:{&pressure_state,&affine_state}) {
			std::vector<iga::FluidSurfaceElementState> local_state;
			for(const auto& item:*all_state)if(static_cast<int>(item.cell_id%ranks)==rank)local_state.push_back(item);
			std::vector<iga::SurfaceCellTractionPoints> points;
			points=iga::BuildDistributedFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,local_state);
			if(ranks>1) {
				auto conflicting=local_state;
				const auto nx=domain.Background().Spec().cells[0]+3,ny=domain.Background().Spec().cells[1]+3;
				const auto shared=2+2*nx+2*nx*ny;
				int present=0,total_present=0;
				for(auto& item:conflicting) {
					const auto element=domain.Background().MaterializeElement(item.cell_id);
					for(std::size_t node=0;node<element.connectivity.size();++node)if(static_cast<std::uint64_t>(element.connectivity[node])==shared) {
						present=1;if(rank==0)item.nodal_state[node][3]+=1.;
					}
				}
				MPI_Allreduce(&present,&total_present,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(total_present>1);
				int rejected=0,total=0;
				try { (void)iga::BuildDistributedFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,conflicting); }
				catch(const std::runtime_error&) { rejected=1; }
				MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(total==ranks);
				points=iga::BuildDistributedFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,local_state);
				for(int mode=0;mode<2;++mode) {
					const auto alternate=rank==0&&mode==1?Material(.75):material;
					assert(alternate.Surface().CanonicalSha256()==material.Surface().CanonicalSha256());
					const double mu=rank==0&&mode==0?2*viscosity:viscosity;
					int failure=0,failures=0;
					try { (void)iga::BuildDistributedFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,alternate,patch_map,mu,owned_cells,local_state); }
					catch(const std::runtime_error&) { failure=1; }
					MPI_Allreduce(&failure,&failures,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(failures==ranks);
					points=iga::BuildDistributedFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,local_state);
				}
			}
			for(int mode=0;mode<3;++mode) {
				auto bad=local_state;
				if(rank==static_cast<int>(all_state->front().cell_id%ranks)) {
					if(mode==0)bad.erase(bad.begin());
					if(mode==1)bad.front().nodal_state.front()[3]=std::numeric_limits<double>::infinity();
					if(mode==2)bad.front().nodal_state.pop_back();
				}
				int rejected=0,total=0;
				try { iga::CollectiveLocalStage(PETSC_COMM_WORLD,"invalid owned IGA traction",[&] {
					(void)iga::BuildOwnedFluidSurfaceTractionPoints(domain,catalog,material,patch_map,viscosity,owned_cells,bad);
				}); } catch(const std::runtime_error&) { rejected=1; }
				MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(total==ranks);
				iga::CollectiveLocalStage(PETSC_COMM_WORLD,"owned IGA traction retry",[&] {
					points=iga::BuildOwnedFluidSurfaceTractionPoints(domain,catalog,material,patch_map,viscosity,owned_cells,local_state);
				});
			}
			iga::FsiTrialContext expected_trial;expected_trial.step=1;expected_trial.start_time_s=0.;expected_trial.dt_s=.5;expected_trial.coupling_iteration=2;
			const auto expected_material=material.ContentIdentitySha256();
			const auto expected_fluid=iga::BuildFluidSurfaceTractionStateIdentitySha256(domain,catalog,material,patch_map,viscosity,local_state);
			points=iga::BuildTrialFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,local_state,expected_trial,expected_material,expected_fluid);
			for(int mode=0;mode<3;++mode) {
				auto context=expected_trial;auto altered=local_state;
				const auto trial_material=mode==0?Material(.75):material;
				if(mode==1)for(auto& item:altered)for(auto& coefficient:item.nodal_state)coefficient[3]+=1.;
				if(mode==2)context.dt_s=.75;
				int rejected=0,total=0;
				try { (void)iga::BuildTrialFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,trial_material,patch_map,viscosity,owned_cells,altered,context,expected_material,expected_fluid); }
				catch(const std::runtime_error&) { rejected=1; }
				MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(total==ranks);
				points=iga::BuildTrialFluidSurfaceTractionPoints(PETSC_COMM_WORLD,domain,catalog,material,patch_map,viscosity,owned_cells,local_state,expected_trial,expected_material,expected_fluid);
			}
			const auto distributed=iga::AssembleDistributedSurfaceTraction(PETSC_COMM_WORLD,interface.id,distributed_layout,domain.Cells().size(),points,0);
			const auto distributed_map=iga::MaterialSurfacePatchMap::Create(patch_map.Interface(),distributed_layout,patch_map.FullReference(),patch_map.PatchLabel(),
				patch_map.GlobalToSourceVertices(),patch_map.LayoutTriangleToSourceTriangles(),patch_map.ConfiguredClampedGlobalNodeIds());
			const auto producer=iga::BuildFluidSurfaceTractionStateIdentitySha256(domain,catalog,material,distributed_map,viscosity,local_state);
			const auto publish=[&] {
				return iga::BuildDistributedFluidSurfaceTraction(PETSC_COMM_WORLD,domain,catalog,material,distributed_map,interface,viscosity,
					owned_cells,local_state,expected_trial,expected_material,producer,0);
			};
			const auto publication=publish();iga::ValidateSurfaceTraction(publication,distributed_layout);
			assert(publication.stamp.step==expected_trial.step&&publication.stamp.coupling_iteration==expected_trial.coupling_iteration);
			assert(publication.stamp.producer_state_identity_sha256==producer);
			for(std::size_t row=0;row<publication.consistent_nodal_force_n.size();++row) {
				CheckVector(publication.consistent_nodal_force_n[row],distributed.owned_force_n[row],1.e-12);
				CheckVector(publication.traction_on_structure_pa[row],distributed.owned_traction_pa[row],1.e-12);
			}
			for(int mode=0;mode<2;++mode) {
				auto incoming=interface;auto expected_producer=producer;
				if(rank==0&&mode==0)incoming.boundary_labels={9};
				if(mode==1)expected_producer=std::string(64,'0');
				int rejected=0,total=0;
				try { (void)iga::BuildDistributedFluidSurfaceTraction(PETSC_COMM_WORLD,domain,catalog,material,distributed_map,incoming,viscosity,
					owned_cells,local_state,expected_trial,expected_material,expected_producer,0); }
				catch(const std::runtime_error&) { rejected=1; }
				MPI_Allreduce(&rejected,&total,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);assert(total==ranks);
				const auto retry=publish();assert(retry.projection_identity_sha256==publication.projection_identity_sha256);
			}
			const auto& oracle=all_state==&pressure_state?pressure:affine_result;
			for(std::size_t row=0;row<distributed_layout.owned_global_node_ids.size();++row) {
				const auto id=distributed_layout.owned_global_node_ids[row];
				const auto index=std::lower_bound(layout.owned_global_node_ids.begin(),layout.owned_global_node_ids.end(),id)-layout.owned_global_node_ids.begin();
				CheckVector(distributed.owned_force_n[row],oracle.traction.consistent_nodal_force_n[index],1.e-12);
				CheckVector(distributed.owned_traction_pa[row],oracle.traction.traction_on_structure_pa[index],1.e-12);
			}
		}
		if(rank==0)std::cout << "owned_iga_traction=passed pressure_and_viscosity ranks=" << ranks << '\n';

		const auto zero_state = State(domain, catalog, 0.0, zero); const auto zero_stamp = Stamp(material, domain, catalog, layout, patch_map, zero_state, viscosity);
		const auto zero_result = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, patch_map, zero_stamp, viscosity, zero_state);
		CheckVector(zero_result.diagnostics.quadrature_resultant_n, {{0,0,0}});
		for (const auto& traction : zero_result.traction.traction_on_structure_pa) CheckVector(traction, {{0,0,0}});

		// Label 9 has three further faces; the selected resultant is still only label 7.
		assert(pressure.diagnostics.retained_quadrature_points < catalog.Diagnostics().output_points);
		assert(pressure.traction.projection_identity_sha256 == iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, patch_map, pressure_stamp, viscosity, pressure_state).traction.projection_identity_sha256);
		assert(pressure.traction.projection_identity_sha256 != affine_result.traction.projection_identity_sha256);
		auto changed_state = pressure_state; changed_state.front().nodal_state.front()[3] += 1.0;
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, material, patch_map, viscosity, changed_state) != pressure_stamp.producer_state_identity_sha256);
		auto changed_map_layout = layout;
		for (auto& id : changed_map_layout.owned_global_node_ids) id += 100;
		for (auto& position : changed_map_layout.reference_positions) position.global_node_id += 100;
		for (auto& triangle : changed_map_layout.reference_triangles) for (auto& id : triangle) id += 100;
		const std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> changed_map_vertices{{111,0},{120,1},{130,2}};
		const std::vector<std::uint32_t> changed_map_triangles{0}; const std::vector<std::uint64_t> changed_map_clamps{111,120,130};
		const auto changed_map_reference = iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(material, 7,
			changed_map_vertices, changed_map_layout.reference_triangles, changed_map_triangles, changed_map_clamps);
		changed_map_layout.reference_mesh_identity_sha256 = changed_map_reference;
		changed_map_layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(changed_map_layout);
		const auto changed_map = iga::MaterialSurfacePatchMap::Create(StructuralInterface(changed_map_reference),
			std::move(changed_map_layout), material, 7, changed_map_vertices, changed_map_triangles, changed_map_clamps);
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, material, changed_map, viscosity, pressure_state)
			!= pressure_stamp.producer_state_identity_sha256);
		RejectContaining("mapped patch authority", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, changed_map, pressure_stamp, viscosity, pressure_state); });

		const auto changed_material = TranslatedMaterial(.1); const auto changed_patch_map = PatchMap(changed_material);
		const auto changed_domain = Domain(changed_material); const iga::ImmersedSurfaceQuadratureCatalog changed_catalog(changed_domain);
		const auto changed_geometry_state = State(changed_domain, changed_catalog, 20.0, zero);
		const auto changed_stamp = Stamp(changed_material, changed_domain, changed_catalog, changed_patch_map.Layout(),
			changed_patch_map, changed_geometry_state, viscosity);
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(changed_domain, changed_catalog, changed_material,
			changed_patch_map, viscosity, changed_geometry_state) != pressure_stamp.producer_state_identity_sha256);
		const auto changed_result = iga::BuildFluidSurfaceTraction(changed_domain, changed_catalog, changed_material,
			Interface(changed_patch_map), changed_patch_map.Layout(), changed_patch_map, changed_stamp, viscosity, changed_geometry_state);
		assert(changed_result.traction.projection_identity_sha256 != pressure.traction.projection_identity_sha256);

		const iga::CubicCartesianGridSpec alternate_grid{{{-.1,-.1,-.1}},{{1.1,1.1,1.1}},{{3,2,2}}};
		const iga::CartesianDomainClassification alternate_domain(iga::CubicCartesianBackground(alternate_grid), iga::SurfaceSpatialIndex(material.Surface()));
		const iga::ImmersedSurfaceQuadratureCatalog alternate_catalog(alternate_domain);
		const auto alternate_state = State(alternate_domain, alternate_catalog, 20.0, zero);
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(alternate_domain, alternate_catalog, material,
			patch_map, viscosity, alternate_state) != pressure_stamp.producer_state_identity_sha256);

		assert(pressure_state.size() > 1);
		auto missing_state = pressure_state; missing_state.pop_back();
		const auto missing_stamp = Stamp(material, domain, catalog, layout, patch_map, missing_state, viscosity);
		RejectContaining("does not cover exactly", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, missing_stamp, viscosity, missing_state); });
		auto duplicate_state = pressure_state; duplicate_state[1] = duplicate_state[0];
		const auto duplicate_stamp = Stamp(material, domain, catalog, layout, patch_map, duplicate_state, viscosity);
		RejectContaining("cell ownership is ambiguous", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, duplicate_stamp, viscosity, duplicate_state); });
		auto reordered_state = pressure_state; std::swap(reordered_state[0], reordered_state[1]);
		const auto reordered_stamp = Stamp(material, domain, catalog, layout, patch_map, reordered_state, viscosity);
		RejectContaining("cell ownership is ambiguous", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, reordered_stamp, viscosity, reordered_state); });
		auto inconsistent_state = pressure_state; std::map<std::int32_t, std::pair<std::size_t, std::size_t>> owner;
		bool found_shared_node = false;
		for (std::size_t item = 0; item < inconsistent_state.size() && !found_shared_node; ++item) {
			const auto element = domain.Background().MaterializeElement(inconsistent_state[item].cell_id);
			for (std::size_t local = 0; local < element.connectivity.size(); ++local) {
				const auto previous = owner.find(element.connectivity[local]);
				if (previous == owner.end()) owner.emplace(element.connectivity[local], std::make_pair(item, local));
				else { inconsistent_state[item].nodal_state[local][0] += .125; found_shared_node = true; break; }
			}
		}
		assert(found_shared_node);
		const auto inconsistent_stamp = Stamp(material, domain, catalog, layout, patch_map, inconsistent_state, viscosity);
		RejectContaining("shared IGA node coefficients are inconsistent", [&] { (void)iga::BuildFluidSurfaceTraction(domain,
			catalog, material, interface, layout, patch_map, inconsistent_stamp, viscosity, inconsistent_state); });

		auto invalid_options = iga::FluidSurfaceTractionProjectionOptions{};
		invalid_options.conservation_absolute_force_tolerance_n = -1.0;
		RejectContaining("parameters are invalid", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, pressure_stamp, viscosity, pressure_state, invalid_options); });
		invalid_options = {}; invalid_options.conservation_absolute_moment_tolerance_n_m = std::numeric_limits<double>::quiet_NaN();
		RejectContaining("parameters are invalid", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, pressure_stamp, viscosity, pressure_state, invalid_options); });
		invalid_options = {}; invalid_options.conservation_absolute_force_tolerance_n = 0.0;
		invalid_options.conservation_absolute_moment_tolerance_n_m = 0.0; invalid_options.conservation_relative_tolerance = 0.0;
		RejectContaining("parameters are invalid", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, pressure_stamp, viscosity, pressure_state, invalid_options); });

		auto cancellation_state = State(domain, catalog, 0.0, zero);
		for (auto& item : cancellation_state) {
			const auto element = domain.Background().MaterializeElement(item.cell_id);
			for (std::size_t node = 0; node < item.nodal_state.size(); ++node) {
				const auto global = static_cast<std::uint32_t>(element.connectivity[node]);
				const auto basis_x = domain.Background().Spec().cells[0]+3, basis_y = domain.Background().Spec().cells[1]+3;
				const auto x = domain.Background().Greville(global%basis_x, (global/basis_x)%basis_y, global/(basis_x*basis_y));
				item.nodal_state[node][3] = x[0]-1.0/3.0;
			}
		}
		const auto cancellation_stamp = Stamp(material, domain, catalog, layout, patch_map, cancellation_state, viscosity);
		const auto cancellation = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, patch_map,
			cancellation_stamp, viscosity, cancellation_state);
		CheckVector(cancellation.diagnostics.quadrature_resultant_n, {{0,0,0}});
		CheckVector(cancellation.diagnostics.quadrature_resultant_n, cancellation.diagnostics.nodal_resultant_n);

		const auto translated_material = TranslatedMaterial(1.0e9); const auto translated_patch_map = PatchMap(translated_material);
		const auto translated_domain = TranslatedDomain(translated_material, 1.0e9);
		const iga::ImmersedSurfaceQuadratureCatalog translated_catalog(translated_domain);
		const auto translated_state = State(translated_domain, translated_catalog, 20.0, zero);
		const auto translated_stamp = Stamp(translated_material, translated_domain, translated_catalog, translated_patch_map.Layout(),
			translated_patch_map, translated_state, viscosity);
		const auto translated = iga::BuildFluidSurfaceTraction(translated_domain, translated_catalog, translated_material,
			Interface(translated_patch_map), translated_patch_map.Layout(), translated_patch_map, translated_stamp, viscosity, translated_state);
		CheckVector(translated.diagnostics.quadrature_resultant_n, translated.diagnostics.nodal_resultant_n);
		CheckVector(translated.diagnostics.quadrature_moment_n_m, translated.diagnostics.nodal_moment_n_m, 1.0e-8);

		auto no_traction_interface = interface;
		no_traction_interface.provides = {iga::SurfaceFieldQuantity::Velocity};
		no_traction_interface.requires = {iga::SurfaceFieldQuantity::Displacement, iga::SurfaceFieldQuantity::TractionOnStructure};
		RejectContaining("exact traction/displacement/velocity role", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog,
			material, no_traction_interface, layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto extra_role_interface = interface;
		extra_role_interface.provides = {iga::SurfaceFieldQuantity::Velocity, iga::SurfaceFieldQuantity::TractionOnStructure};
		extra_role_interface.requires = {iga::SurfaceFieldQuantity::Displacement};
		RejectContaining("exact traction/displacement/velocity role", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog,
			material, extra_role_interface, layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto bad_role = interface; bad_role.requires = {iga::SurfaceFieldQuantity::Displacement};
		RejectContaining("exact traction/displacement/velocity role", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog,
			material, bad_role, layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto bad_label = interface; bad_label.boundary_labels = {9};
		RejectContaining("mapped patch authority", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			bad_label, layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto bad_layout = layout; bad_layout.partition_count = 2;
		bad_layout.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(bad_layout);
		RejectContaining("complete single-rank partition", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, bad_layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto bad_interface = interface; bad_interface.reference_mesh_identity_sha256 = material.ContentIdentitySha256();
		RejectContaining("mapped patch authority", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			bad_interface, layout, patch_map, pressure_stamp, viscosity, pressure_state); });
		auto bad_stamp = pressure_stamp; bad_stamp.producer_state_identity_sha256 = affine_stamp.producer_state_identity_sha256;
		RejectContaining("producer state identity", [&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material,
			interface, layout, patch_map, bad_stamp, viscosity, pressure_state); });
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; status = 1;
	}
	PetscFinalize(); return status;
}

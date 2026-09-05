#include "FluidSurfaceTraction.hpp"
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

void CheckVector(const std::array<double, 3>& actual, const std::array<double, 3>& expected)
{ for (int component = 0; component < 3; ++component) assert(Near(actual[component], expected[component])); }

template <class Function> void Reject(Function&& function)
{ bool rejected = false; try { function(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

iga::RawSurfaceTriangle Face(int a, int b, int c, std::uint32_t label)
{ iga::RawSurfaceTriangle result; result.indices = {{a,b,c}}; result.boundary_id = label; return result; }

iga::MaterialSurfaceKinematics Material()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}};
	// A-C-B has outward -z; this is the selected planar material patch.
	soup.triangles = {Face(0,2,1,7), Face(0,1,3,9), Face(0,3,2,9), Face(1,2,3,9)};
	iga::PrescribedSurfaceFrame first{0.0,soup}, second{1.0,soup};
	return iga::PrescribedSurfaceMotion({first,second}).Evaluate(0.5,0.0,0.5);
}

iga::DistributedSurfaceInterface Interface(const iga::MaterialSurfaceKinematics& material)
{
	iga::DistributedSurfaceInterface result;
	result.id = {"fluid", "immersed", "membrane"}; result.subsystem_id = "immersed";
	result.boundary_labels = {7}; result.reference_mesh_identity_sha256 = material.MaterialIdentitySha256();
	result.provides = {iga::SurfaceFieldQuantity::TractionOnStructure};
	result.requires = {iga::SurfaceFieldQuantity::Displacement};
	return result;
}

iga::DistributedSurfaceLayout Layout(const iga::MaterialSurfaceKinematics& material)
{
	iga::DistributedSurfaceLayout result;
	result.reference_mesh_identity_sha256 = material.MaterialIdentitySha256(); result.global_node_count = 3;
	result.partition_count = 1; result.partition_rank = 0; result.owned_global_node_ids = {11,20,30};
	// The source patch is A-C-B; layout IDs intentionally do not coincide with source IDs.
	result.reference_positions = {{11,material.ReferenceMaterialVerticesM()[0]}, {20,material.ReferenceMaterialVerticesM()[1]},
		{30,material.ReferenceMaterialVerticesM()[2]}};
	result.reference_triangles = {{{{11,20,30}}}};
	result.owned_reference_lumped_areas_m2 = {1.0/6.0,1.0/6.0,1.0/6.0};
	result.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(result);
	return result;
}

iga::DistributedSurfaceLayout AllSurfaceLayout(const iga::MaterialSurfaceKinematics& material)
{
	iga::DistributedSurfaceLayout result;
	result.reference_mesh_identity_sha256 = material.MaterialIdentitySha256(); result.global_node_count = 4;
	result.partition_count = 1; result.partition_rank = 0; result.owned_global_node_ids = {11,20,30,40};
	for (std::size_t node = 0; node < 4; ++node)
		result.reference_positions.push_back({result.owned_global_node_ids[node], material.ReferenceMaterialVerticesM()[node]});
	result.reference_triangles = {{{{11,20,30}}, {{11,30,40}}, {{11,20,40}}, {{20,30,40}}}};
	result.owned_reference_lumped_areas_m2 = {.5,.5,.5,.5};
	result.layout_identity_sha256 = iga::BuildDistributedSurfaceLayoutIdentitySha256(result);
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
	const iga::DistributedSurfaceLayout& layout, const std::vector<iga::FluidSurfaceElementState>& state, double viscosity)
{
	iga::SurfaceFieldStamp result; result.time_s = material.EvaluatedTimeS(); result.step = 4; result.coupling_iteration = 2;
	result.reference_mesh_identity_sha256 = layout.reference_mesh_identity_sha256; result.layout_identity_sha256 = layout.layout_identity_sha256;
	result.partition_identity_sha256 = iga::BuildDistributedSurfacePartitionIdentitySha256(layout);
	result.producer_state_identity_sha256 = iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, material, viscosity, state);
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		const auto material = Material(); const auto interface = Interface(material); const auto layout = Layout(material);
		const auto domain = Domain(material); const iga::ImmersedSurfaceQuadratureCatalog catalog(domain);
		const std::array<std::array<double,3>,3> zero{}; constexpr double viscosity = 2.0;
		const auto pressure_state = State(domain, catalog, 20.0, zero); const auto pressure_stamp = Stamp(material, domain, catalog, layout, pressure_state, viscosity);
		const auto pressure = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, pressure_state);
		// Selected face has n_f=-z, so -sigma n_f=p n_f=(0,0,-20) Pa.
		for (const auto& traction : pressure.traction.traction_on_structure_pa) CheckVector(traction, {{0,0,-20}});
		CheckVector(pressure.diagnostics.quadrature_resultant_n, {{0,0,-10}});
		CheckVector(pressure.diagnostics.nodal_resultant_n, {{0,0,-10}});
		assert(pressure.diagnostics.retained_quadrature_points > 0);

		const std::array<std::array<double,3>,3> affine{{{{0,0,1}},{{0,0,0}},{{2,0,3}}}};
		const auto affine_state = State(domain, catalog, 4.0, affine); const auto affine_stamp = Stamp(material, domain, catalog, layout, affine_state, viscosity);
		const auto affine_result = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, affine_stamp, viscosity, affine_state);
		// p n - mu(grad u + grad u^T)n = (6,0,8) Pa on this patch.
		for (const auto& traction : affine_result.traction.traction_on_structure_pa) CheckVector(traction, {{6,0,8}});
		CheckVector(affine_result.diagnostics.quadrature_resultant_n, {{3,0,4}});
		CheckVector(affine_result.diagnostics.nodal_resultant_n, {{3,0,4}});
		CheckVector(affine_result.diagnostics.quadrature_moment_n_m, affine_result.diagnostics.nodal_moment_n_m);

		const auto zero_state = State(domain, catalog, 0.0, zero); const auto zero_stamp = Stamp(material, domain, catalog, layout, zero_state, viscosity);
		const auto zero_result = iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, zero_stamp, viscosity, zero_state);
		CheckVector(zero_result.diagnostics.quadrature_resultant_n, {{0,0,0}});
		for (const auto& traction : zero_result.traction.traction_on_structure_pa) CheckVector(traction, {{0,0,0}});

		// Label 9 has three further faces; the selected resultant is still only label 7.
		assert(pressure.diagnostics.retained_quadrature_points < catalog.Diagnostics().output_points);
		assert(pressure.traction.projection_identity_sha256 == iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, pressure_state).traction.projection_identity_sha256);
		assert(pressure.traction.projection_identity_sha256 != affine_result.traction.projection_identity_sha256);
		auto changed_state = pressure_state; changed_state.front().nodal_state.front()[3] += 1.0;
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, material, viscosity, changed_state) != pressure_stamp.producer_state_identity_sha256);
		auto moved_soup = iga::RawSurfaceSoup{}; moved_soup.vertices = {{{{.1,0,0}},{{1.1,0,0}},{{.1,1,0}},{{.1,0,1}}}};
		moved_soup.triangles = {Face(0,2,1,7), Face(0,1,3,9), Face(0,3,2,9), Face(1,2,3,9)};
		const auto changed_material = iga::PrescribedSurfaceMotion({{0.0,moved_soup},{1.0,moved_soup}}).Evaluate(.5,0.0,.5);
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(domain, catalog, changed_material, viscosity, pressure_state) != pressure_stamp.producer_state_identity_sha256);
		const auto changed_domain = Domain(changed_material); const iga::ImmersedSurfaceQuadratureCatalog changed_catalog(changed_domain);
		const auto changed_layout = Layout(changed_material); const auto changed_interface = Interface(changed_material);
		const auto changed_geometry_state = State(changed_domain, changed_catalog, 20.0, zero);
		const auto changed_geometry_stamp = Stamp(changed_material, changed_domain, changed_catalog, changed_layout, changed_geometry_state, viscosity);
		const auto changed_geometry = iga::BuildFluidSurfaceTraction(changed_domain, changed_catalog, changed_material,
			changed_interface, changed_layout, changed_geometry_stamp, viscosity, changed_geometry_state);
		assert(changed_geometry.traction.projection_identity_sha256 != pressure.traction.projection_identity_sha256);
		const iga::CubicCartesianGridSpec alternate_grid{{{-.1,-.1,-.1}},{{1.1,1.1,1.1}},{{3,2,2}}};
		const iga::CartesianDomainClassification alternate_domain(iga::CubicCartesianBackground(alternate_grid), iga::SurfaceSpatialIndex(material.Surface()));
		const iga::ImmersedSurfaceQuadratureCatalog alternate_catalog(alternate_domain);
		const auto alternate_state = State(alternate_domain, alternate_catalog, 20.0, zero);
		assert(iga::BuildFluidSurfaceTractionStateIdentitySha256(alternate_domain, alternate_catalog, material, viscosity, alternate_state)
			!= pressure_stamp.producer_state_identity_sha256);

		assert(pressure_state.size() > 1);
		auto missing_state = pressure_state; missing_state.pop_back();
		const auto missing_stamp = Stamp(material, domain, catalog, layout, missing_state, viscosity);
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, missing_stamp, viscosity, missing_state); });
		auto duplicate_state = pressure_state; duplicate_state[1] = duplicate_state[0];
		const auto duplicate_stamp = Stamp(material, domain, catalog, layout, duplicate_state, viscosity);
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, duplicate_stamp, viscosity, duplicate_state); });
		auto reordered_state = pressure_state; std::swap(reordered_state[0], reordered_state[1]);
		const auto reordered_stamp = Stamp(material, domain, catalog, layout, reordered_state, viscosity);
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, reordered_stamp, viscosity, reordered_state); });
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
		const auto inconsistent_stamp = Stamp(material, domain, catalog, layout, inconsistent_state, viscosity);
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, inconsistent_stamp, viscosity, inconsistent_state); });

		auto no_traction_interface = interface; no_traction_interface.provides = {iga::SurfaceFieldQuantity::Velocity};
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, no_traction_interface, layout, pressure_stamp, viscosity, pressure_state); });
		auto invalid_options = iga::FluidSurfaceTractionProjectionOptions{}; invalid_options.conservation_absolute_force_tolerance_n = -1.0;
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, pressure_state, invalid_options); });
		invalid_options = {}; invalid_options.conservation_absolute_moment_tolerance_n_m = std::numeric_limits<double>::quiet_NaN();
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, pressure_state, invalid_options); });
		invalid_options = {}; invalid_options.conservation_absolute_force_tolerance_n = 0.0;
		invalid_options.conservation_absolute_moment_tolerance_n_m = 0.0; invalid_options.conservation_relative_tolerance = 0.0;
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, pressure_state, invalid_options); });

		auto all_interface = interface; all_interface.boundary_labels = {7,9};
		const auto all_layout = AllSurfaceLayout(material); const auto all_state = State(domain, catalog, 20.0, zero, {7,9});
		const auto all_stamp = Stamp(material, domain, catalog, all_layout, all_state, viscosity);
		const auto cancellation = iga::BuildFluidSurfaceTraction(domain, catalog, material, all_interface, all_layout, all_stamp, viscosity, all_state);
		CheckVector(cancellation.diagnostics.quadrature_resultant_n, {{0,0,0}});
		CheckVector(cancellation.diagnostics.nodal_resultant_n, {{0,0,0}});
		const auto translated_material = TranslatedMaterial(1.0e9);
		const auto translated_domain = TranslatedDomain(translated_material, 1.0e9);
		const iga::ImmersedSurfaceQuadratureCatalog translated_catalog(translated_domain);
		const auto translated_layout = Layout(translated_material); const auto translated_interface = Interface(translated_material);
		const auto translated_state = State(translated_domain, translated_catalog, 20.0, zero);
		const auto translated_stamp = Stamp(translated_material, translated_domain, translated_catalog, translated_layout, translated_state, viscosity);
		const auto translated = iga::BuildFluidSurfaceTraction(translated_domain, translated_catalog, translated_material,
			translated_interface, translated_layout, translated_stamp, viscosity, translated_state);
		CheckVector(translated.diagnostics.quadrature_resultant_n, translated.diagnostics.nodal_resultant_n);
		CheckVector(translated.diagnostics.quadrature_moment_n_m, translated.diagnostics.nodal_moment_n_m);

		auto bad_stamp = pressure_stamp; bad_stamp.producer_state_identity_sha256 = affine_stamp.producer_state_identity_sha256;
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, bad_stamp, viscosity, pressure_state); });
		auto bad_layout = layout; bad_layout.partition_count = 2;
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, bad_layout, pressure_stamp, viscosity, pressure_state); });
		auto bad_interface = interface; bad_interface.reference_mesh_identity_sha256 = material.ContentIdentitySha256();
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, bad_interface, layout, pressure_stamp, viscosity, pressure_state); });
		Reject([&] { (void)iga::BuildFluidSurfaceTraction(domain, catalog, material, interface, layout, pressure_stamp, viscosity, affine_state); });
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n'; status = 1;
	}
	PetscFinalize(); return status;
}

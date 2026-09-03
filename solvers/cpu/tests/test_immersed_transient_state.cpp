#include "ImmersedTransientState.hpp"
#include "NavierStokesElement.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>

namespace {
iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	iga::RawSurfaceTriangle r; r.indices={{a,b,c}}; r.boundary_id=0; return r;
}
iga::RawSurfaceSoup Cube()
{
	iga::RawSurfaceSoup r; r.vertices={{{{.2,.2,.2}},{{.8,.2,.2}},{{.8,.8,.2}},{{.2,.8,.2}},{{.2,.2,.8}},{{.8,.2,.8}},{{.8,.8,.8}},{{.2,.8,.8}}}};
	r.triangles={Face(0,2,1),Face(0,3,2),Face(4,5,6),Face(4,6,7),Face(0,1,5),Face(0,5,4),Face(1,2,6),Face(1,6,5),Face(2,3,7),Face(2,7,6),Face(3,0,4),Face(3,4,7)}; return r;
}
iga::RawSurfaceSoup WideCube()
{
	iga::RawSurfaceSoup r; r.vertices={{{{.1,.2,.2}},{{1.9,.2,.2}},{{1.9,.8,.2}},{{.1,.8,.2}},{{.1,.2,.8}},{{1.9,.2,.8}},{{1.9,.8,.8}},{{.1,.8,.8}}}};
	r.triangles={Face(0,2,1),Face(0,3,2),Face(4,5,6),Face(4,6,7),Face(0,1,5),Face(0,5,4),Face(1,2,6),Face(1,6,5),Face(2,3,7),Face(2,7,6),Face(3,0,4),Face(3,4,7)}; return r;
}
template<class F> void Reject(F&& f) { bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}assert(rejected); }
}
int main()
{
	const iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{1,1,1}}};
	const iga::CartesianDomainClassification domain(iga::CubicCartesianBackground(spec),iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(Cube())));
	const iga::CutCellVolumeQuadratureCatalog volume(domain,{4,300000,300000,3000000});
	const auto layout=iga::ImmersedActiveLayout::Build(domain,volume,"fixed-geometry",{3,9},true);
	const auto repeat=iga::ImmersedActiveLayout::Build(domain,volume,"fixed-geometry",{3,9},true);
	assert(layout.Valid()&&layout.HashSha256()==repeat.HashSha256()&&layout.ControllerRow(3)<layout.ControllerRow(9)&&layout.GaugeRow()+1==layout.Rows());
	std::vector<std::array<double,4>> fields(layout.NodeIds().size());
	for(std::size_t i=0;i<fields.size();++i)fields[i]={{double(i),double(i)+.25,double(i)-.5,1000.0+double(i)}};
	const iga::ImmersedGlobalFlowState state(.1,2,layout,fields,{.3,.9},true,1.2);
	const auto history=iga::BuildIdentityImmersedVelocityHistory(state,layout,.2);
	assert(state.Valid()&&history.Valid()&&history.Provenance().size()==fields.size()); history.ValidateCoverage(layout);
	for(std::size_t i=0;i<fields.size();++i) { assert(history.Provenance()[i]==iga::ImmersedVelocityHistoryProvenance::Committed); for(int c=0;c<3;++c)assert(std::memcmp(&history.Velocities()[i][c],&fields[i][c],sizeof(double))==0); }
	Reject([&]{iga::ImmersedActiveLayout::Build(domain,volume,"fixed-geometry",{9,3});});
	Reject([&]{iga::ImmersedGlobalFlowState(std::numeric_limits<double>::quiet_NaN(),0,layout,fields,{.3,.9},true,0.);});
	Reject([&]{iga::ImmersedVelocityHistory(.1,.2,"a","b",{1,1},{{{0,0,0}},{{0,0,0}}},{iga::ImmersedVelocityHistoryProvenance::Committed,iga::ImmersedVelocityHistoryProvenance::Committed});});
	const auto layout_without_gauge=iga::ImmersedActiveLayout::Build(domain,volume,"fixed-geometry",{3,9},false);
	Reject([&]{iga::ImmersedGlobalFlowState(.1,2,layout_without_gauge,fields,{.3,.9},false,1.0);});
	const iga::ImmersedGlobalFlowState zero_gauge(.1,2,layout_without_gauge,fields,{.3,.9},false,-0.0);
	assert(zero_gauge.GaugeMultiplier() == 0.0 && zero_gauge.HashSha256()
		== iga::ImmersedGlobalFlowState(.1,2,layout_without_gauge,fields,{.3,.9},false,0.0).HashSha256());
	const auto element=domain.Background().MaterializeElement(0);
	const auto transient=iga::BuildTransientNavierStokesElement(element,fields,history,layout,.2,
		{1.0,.1,.1},volume.UsableRule(domain,0));
	assert(!transient.jacobian.empty() && !transient.negative_residual.empty());
	// The named global-ID path retains its body-fitted default for compatibility,
	// while the immersed runtime explicitly selects the conservative pair.  This
	// is a residual-level discriminator, independent of the FD tangent checks.
	const auto conservative=iga::BuildTransientNavierStokesElement(element,fields,history,layout,.2,
		{1.0,.1,.1},volume.UsableRule(domain,0),[](const std::array<double,3>&){return std::array<double,3>{{0,0,0}};},
		iga::NavierStokesResolvedMixedForm::Conservative);
	double mixed_form_difference=0.0;
	for(std::size_t i=0;i<transient.negative_residual.size();++i)
		mixed_form_difference=std::max(mixed_form_difference,std::abs(PetscRealPart(conservative.negative_residual[i]-transient.negative_residual[i])));
	assert(mixed_form_difference>1e-10);
	const std::vector<iga::ImmersedVelocityHistoryProvenance> committed(fields.size(),
		iga::ImmersedVelocityHistoryProvenance::Committed);
	const auto tiny_history=iga::ImmersedVelocityHistory(0.0,1e-16,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	const auto tiny_transient=iga::BuildTransientNavierStokesElement(element,fields,tiny_history,layout,1e-16,
		{1.0,.1,1e-16},volume.UsableRule(domain,0));
	assert(!tiny_transient.jacobian.empty() && !tiny_transient.negative_residual.empty());
	const double epsilon=std::numeric_limits<double>::epsilon();
	const double source_time=1.0;
	const double expected_target=source_time+epsilon;
	assert(expected_target > source_time);
	const auto exact_history=iga::ImmersedVelocityHistory(source_time,expected_target,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	const auto exact_transient=iga::BuildTransientNavierStokesElement(element,fields,exact_history,layout,expected_target,
		{1.0,.1,epsilon},volume.UsableRule(domain,0));
	assert(!exact_transient.jacobian.empty() && !exact_transient.negative_residual.empty());
	const auto two_ulp_late_history=iga::ImmersedVelocityHistory(source_time,source_time+2.0*epsilon,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,two_ulp_late_history,layout,source_time+2.0*epsilon,
		{1.0,.1,epsilon},volume.UsableRule(domain,0));});
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,exact_history,layout,
		std::nextafter(expected_target,std::numeric_limits<double>::infinity()),{1.0,.1,epsilon},volume.UsableRule(domain,0));});
	const auto same_time_history=iga::ImmersedVelocityHistory(source_time,source_time,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,same_time_history,layout,source_time,
		{1.0,.1,epsilon},volume.UsableRule(domain,0));});
	const double large_source=1e20;
	const auto nonadvancing_history=iga::ImmersedVelocityHistory(large_source,large_source,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,nonadvancing_history,layout,large_source,
		{1.0,.1,1e-16},volume.UsableRule(domain,0));});
	Reject([&]{iga::CheckedTransientTargetTime(std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max());});
	const auto wrong_source_geometry=iga::ImmersedVelocityHistory(.1,.2,"other-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),committed);
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,wrong_source_geometry,layout,.2,
		{1.0,.1,.1},volume.UsableRule(domain,0));});
	const auto extended_history=iga::ImmersedVelocityHistory(.1,.2,"fixed-geometry","fixed-geometry",
		layout.NodeIds(),history.Velocities(),std::vector<iga::ImmersedVelocityHistoryProvenance>(fields.size(),
			iga::ImmersedVelocityHistoryProvenance::Extended));
	Reject([&]{iga::BuildTransientNavierStokesElement(element,fields,extended_history,layout,.2,
		{1.0,.1,.1},volume.UsableRule(domain,0));});
	const iga::CubicCartesianGridSpec two_cell_spec{{{0,0,0}},{{2,1,1}},{{2,1,1}}};
	const iga::CartesianDomainClassification two_cell_domain(iga::CubicCartesianBackground(two_cell_spec),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(WideCube())));
	const iga::CutCellVolumeQuadratureCatalog two_cell_volume(two_cell_domain,{4,300000,300000,3000000});
	const auto two_cell_layout=iga::ImmersedActiveLayout::Build(two_cell_domain,two_cell_volume,"two-cell");
	std::vector<std::array<double,3>> global_velocity(two_cell_layout.NodeIds().size());
	for(std::size_t i=0;i<global_velocity.size();++i) global_velocity[i]={{double(two_cell_layout.NodeIds()[i]),double(two_cell_layout.NodeIds()[i])+10.0,double(two_cell_layout.NodeIds()[i])-10.0}};
	const iga::ImmersedVelocityHistory global_history(0.0,.1,"two-cell","two-cell",two_cell_layout.NodeIds(),global_velocity,
		std::vector<iga::ImmersedVelocityHistoryProvenance>(global_velocity.size(),iga::ImmersedVelocityHistoryProvenance::Committed));
	for(std::uint64_t cell_id=0;cell_id<2;++cell_id) {
		const auto local=iga::LocalizeImmersedVelocityHistory(two_cell_domain.Background().MaterializeElement(cell_id),two_cell_layout,global_history,.1);
		const auto local_element=two_cell_domain.Background().MaterializeElement(cell_id);
		for(std::size_t a=0;a<local.size();++a) assert(local[a][0] == double(local_element.connectivity[a]));
	}
	std::vector<std::int32_t> missing_ids=two_cell_layout.NodeIds(); missing_ids.pop_back();
	std::vector<std::array<double,3>> missing_velocity=global_velocity; missing_velocity.pop_back();
	const iga::ImmersedVelocityHistory missing_history(0.0,.1,"two-cell","two-cell",missing_ids,missing_velocity,
		std::vector<iga::ImmersedVelocityHistoryProvenance>(missing_velocity.size(),iga::ImmersedVelocityHistoryProvenance::Committed));
	Reject([&]{iga::LocalizeImmersedVelocityHistory(two_cell_domain.Background().MaterializeElement(0),two_cell_layout,missing_history,.1);});
	Reject([&]{iga::LocalizeImmersedVelocityHistory(two_cell_domain.Background().MaterializeElement(0),two_cell_layout,global_history,.2);});
	auto invalid_element=two_cell_domain.Background().MaterializeElement(0); invalid_element.connectivity[0]=std::numeric_limits<std::int32_t>::max();
	Reject([&]{iga::LocalizeImmersedVelocityHistory(invalid_element,two_cell_layout,global_history,.1);});
	std::cout<<"immersed transient state tests passed\n";
}

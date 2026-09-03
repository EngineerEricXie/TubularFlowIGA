#include "MovingImmersedFlowSnapshot.hpp"

#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace iga;
constexpr double kFieldTolerance = 3.0e-11;
constexpr double kMetricTolerance = 1.0e-12;

void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}
void CheckClose(double actual, double expected, const char* message, double tolerance = kMetricTolerance)
{
	Check(std::isfinite(actual) && std::isfinite(expected) && std::abs(actual-expected) <= tolerance*std::max({1.0,std::abs(actual),std::abs(expected)}), message);
}
void ExpectUnpublished(const std::function<void()>& action, const char* message)
{
	bool published = false;
	try { action(); published = true; }
	catch (const std::exception& error) {
		Check(error.what() == std::string(message), "wrong rejection gate");
		Check(!published, "invalid snapshot was published");
		return;
	}
	throw std::runtime_error("missing rejection gate");
}

RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c)
{
	RawSurfaceTriangle value; value.indices={{a,b,c}}; value.boundary_id=1; return value;
}
RawSurfaceSoup Tetra(double shift=0.0)
{
	RawSurfaceSoup value;
	value.vertices={{{{shift+.18,.25,.25}},{{shift+.78,.25,.25}},{{shift+.18,.85,.25}},{{shift+.18,.25,.85}}}};
	value.triangles={Face(0,2,1),Face(0,1,3),Face(0,3,2),Face(1,2,3)};
	return value;
}
MovingCutGeometryOptions Options(CutCellVolumeQuadratureStorageMode storage)
{
	MovingCutGeometryOptions value; value.volume_storage=storage; value.volume.max_depth=3;
	value.volume.max_nodes=value.volume.max_leaves=value.volume.max_points=value.volume.max_records=value.volume.max_logical_points=1000000;
	value.volume.max_retained_bytes=100000000; return value;
}
struct Fixture { CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}}; std::unique_ptr<MovingCutGeometry> geometry; ImmersedActiveLayout layout; };
Fixture Make(CutCellVolumeQuadratureStorageMode storage=CutCellVolumeQuadratureStorageMode::Expanded)
{
	Fixture value; PrescribedSurfaceMotion motion({{0.,Tetra()},{1.,Tetra()}});
	value.geometry=MovingCutGeometry::Build(value.grid,motion.Evaluate(0,0,1),Options(storage));
	value.layout=ImmersedActiveLayout::Build(value.geometry->Domain(),value.geometry->Volume(),value.geometry->GeometryIdentitySha256(),{7,9}); return value;
}
Fixture MakeTranslated()
{
	Fixture value; PrescribedSurfaceMotion motion({{0.,Tetra()},{1.,Tetra(.05)}});
	value.geometry=MovingCutGeometry::Build(value.grid,motion.Evaluate(1,0,1),Options(CutCellVolumeQuadratureStorageMode::Expanded));
	value.layout=ImmersedActiveLayout::Build(value.geometry->Domain(),value.geometry->Volume(),value.geometry->GeometryIdentitySha256(),{7,9}); return value;
}
double Knot(const CubicCartesianGridSpec& grid, int axis, unsigned index)
{
	if(index<=3) return grid.lower_m[axis];
	if(index>=grid.cells[axis]+3) return grid.upper_m[axis];
	return grid.lower_m[axis]+(grid.upper_m[axis]-grid.lower_m[axis])*(index-3)/grid.cells[axis];
}
std::array<unsigned,3> NodeIndex(std::int32_t id,const CubicCartesianGridSpec& grid)
{
	const unsigned nx=grid.cells[0]+3,ny=grid.cells[1]+3; return {{unsigned(id)%nx,(unsigned(id)/nx)%ny,unsigned(id)/(nx*ny)}};
}
std::array<double,3> GrevillePoint(std::int32_t id,const CubicCartesianGridSpec& grid)
{
	const auto index=NodeIndex(id,grid); std::array<double,3> point{};
	for(int axis=0;axis<3;++axis)
		point[axis]=(Knot(grid,axis,index[axis]+1)+Knot(grid,axis,index[axis]+2)+Knot(grid,axis,index[axis]+3))/3.;
	return point;
}
std::array<double,3> AffineVelocity(const std::array<double,3>& x)
{
	return {{.11+.40*x[0]-.30*x[1]+.20*x[2],-.07+.10*x[0]+.25*x[1]-.40*x[2],.50-.20*x[0]+.30*x[1]+.10*x[2]}};
}
double AffinePressure(const std::array<double,3>& x) { return .80+.90*x[0]-.60*x[1]+.40*x[2]; }
std::vector<std::array<double,4>> Affine(const ImmersedActiveLayout& layout,const CubicCartesianGridSpec& grid)
{
	std::vector<std::array<double,4>> value; value.reserve(layout.NodeIds().size()); for(auto id:layout.NodeIds()){const auto x=GrevillePoint(id,grid);const auto u=AffineVelocity(x);value.push_back({{u[0],u[1],u[2],AffinePressure(x)}});} return value;
}
std::vector<std::array<double,4>> Rotation(const ImmersedActiveLayout& layout,const CubicCartesianGridSpec& grid,double omega)
{
	std::vector<std::array<double,4>> value; value.reserve(layout.NodeIds().size()); for(auto id:layout.NodeIds()){const auto x=GrevillePoint(id,grid);value.push_back({{-omega*x[1],omega*x[0],0.,3.25}});} return value;
}
std::vector<std::array<double,4>> SymmetricStrain(const ImmersedActiveLayout& layout,const CubicCartesianGridSpec& grid,double a)
{
	std::vector<std::array<double,4>> value; value.reserve(layout.NodeIds().size()); for(auto id:layout.NodeIds()){const auto x=GrevillePoint(id,grid);value.push_back({{a*x[0],-a*x[1],0.,0.}});} return value;
}
std::vector<std::array<double,4>> SignThresholdField(const ImmersedActiveLayout& layout,const CubicCartesianGridSpec& grid)
{
	std::vector<std::array<double,4>> value; value.reserve(layout.NodeIds().size()); for(auto id:layout.NodeIds()){const auto x=GrevillePoint(id,grid);value.push_back({{x[1],-3.*(x[0]-.5)*(x[0]-.5),0.,0.}});} return value;
}
std::vector<std::array<double,4>> Constant(const ImmersedActiveLayout& layout,std::array<double,3> u)
{
	std::vector<std::array<double,4>> value(layout.NodeIds().size());for(auto& q:value)q={{u[0],u[1],u[2],-1.}};return value;
}
MovingImmersedFlowSnapshotRequest Request(const Fixture& fixture,double dt=.25)
{
	MovingImmersedFlowSnapshotRequest value;value.time_s=fixture.geometry->Evaluation().EvaluatedTimeS();value.index=value.time_s==0.?4:5;value.transition_available=true;value.dt_s=dt;value.port_flows={{7,-.03},{9,.01}};value.wall_labels={1};return value;
}
double ExpectedWallArea(const MovingCutGeometry& geometry,std::uint32_t label)
{
	long double area=0.;for(std::uint64_t id=0;id<geometry.Surface().Cells().size();++id){const auto& rule=geometry.Surface().UsableRule(geometry.Domain(),id);for(const auto& point:rule.Points())if(point.boundary_id>=0&&static_cast<std::uint64_t>(point.boundary_id)==label)area+=point.weight;}return static_cast<double>(area);
}
std::array<double,3> ExpectedPhysical(const MovingCutGeometry& geometry, std::uint64_t cell_id,
	const std::array<double,3>& parametric)
{
	const auto cell = geometry.Domain().Background().Cell(cell_id);
	std::array<double,3> physical{};
	for (std::size_t d = 0; d < 3; ++d)
		physical[d] = cell.lower_m[d]+(cell.upper_m[d]-cell.lower_m[d])*parametric[d];
	return physical;
}
void CheckEquivalentMetrics(const MovingImmersedFlowSnapshotMetrics& a,const MovingImmersedFlowSnapshotMetrics& b)
{
	for(const auto& pair:std::vector<std::pair<double,double>>{{a.quadrature_volume_m3,b.quadrature_volume_m3},{a.enstrophy_integral_m3_per_s2,b.enstrophy_integral_m3_per_s2},{a.mean_enstrophy_per_s2,b.mean_enstrophy_per_s2},{a.mean_q_criterion_per_s2,b.mean_q_criterion_per_s2}})CheckClose(pair.first,pair.second,"expanded/compact polynomial-moment metric mismatch");
}
void CheckIndicatorEstimates(const MovingImmersedFlowSnapshotMetrics& metrics)
{
	for(const auto& pair:std::vector<std::pair<double,double>>{{metrics.q_positive_volume_m3,metrics.q_positive_volume_fraction},{metrics.stagnant_volume_m3,metrics.stagnant_volume_fraction}}){
		Check(std::isfinite(pair.first)&&std::isfinite(pair.second)&&pair.first>=0.&&pair.first<=metrics.quadrature_volume_m3&&pair.second>=0.&&pair.second<=1.,"indicator estimate is outside its valid range");
	}
}
}

int main()
{
	try {
		auto expanded=Make(); MovingImmersedFlowSnapshotOptions options; options.stagnant_speed_threshold_m_per_s=0.;
		const auto affine_state=ImmersedGlobalFlowState(0.,4,expanded.layout,Affine(expanded.layout,expanded.grid),{0.,0.},false,0.);const auto affine=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,affine_state,Request(expanded),options);Check(!affine.Points().empty(),"affine snapshot is empty");
		for(const auto& point:affine.Points()){const auto physical=ExpectedPhysical(*expanded.geometry,point.background_cell_id,point.parametric);for(std::size_t d=0;d<3;++d)CheckClose(point.physical_m[d],physical[d],"Cartesian physical coordinate mismatch",kFieldTolerance);const auto expected=AffineVelocity(physical);for(std::size_t d=0;d<3;++d)CheckClose(point.velocity_m_per_s[d],expected[d],"affine velocity interpolation mismatch",kFieldTolerance);CheckClose(point.pressure,AffinePressure(physical),"affine pressure interpolation mismatch",kFieldTolerance);}

		constexpr double omega=1.75;const auto rotation_state=ImmersedGlobalFlowState(0.,4,expanded.layout,Rotation(expanded.layout,expanded.grid,omega),{0.,0.},false,0.);const auto request=Request(expanded);const auto snapshot=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,request,options);const auto repeated=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,request,options);const auto& metrics=snapshot.Metrics();Check(snapshot.ContentHashSha256()==repeated.ContentHashSha256()&&snapshot.SnapshotIdentitySha256()==repeated.SnapshotIdentitySha256(),"snapshot hashes are not deterministic");
		for(const auto& point:snapshot.Points()){CheckClose(point.vorticity_per_s[0],0.,"rotation curl-x mismatch",kFieldTolerance);CheckClose(point.vorticity_per_s[1],0.,"rotation curl-y mismatch",kFieldTolerance);CheckClose(point.vorticity_per_s[2],2.*omega,"rotation curl-z mismatch",kFieldTolerance);CheckClose(point.q_criterion_per_s2,omega*omega,"rotation Q mismatch",kFieldTolerance);CheckClose(point.enstrophy_density_per_s2,2.*omega*omega,"rotation enstrophy mismatch",kFieldTolerance);}
		CheckClose(metrics.inlet_flow_m3_s,.03,"inlet port map mismatch");CheckClose(metrics.outlet_flow_m3_s,.01,"outlet port map mismatch");const double expected_rate=.03/metrics.quadrature_volume_m3,expected_replacement=-std::expm1(-.03*request.dt_s/metrics.quadrature_volume_m3);CheckClose(metrics.endpoint_turnover_rate_per_s,expected_rate,"turnover-rate formula changed");Check(metrics.endpoint_turnover_time_s.has_value(),"positive inflow must have a turnover time");CheckClose(*metrics.endpoint_turnover_time_s,metrics.quadrature_volume_m3/.03,"turnover-time formula changed");CheckClose(metrics.well_mixed_replacement_fraction_over_step,expected_replacement,"replacement formula changed");

		constexpr double strain=.83;const auto strain_state=ImmersedGlobalFlowState(0.,4,expanded.layout,SymmetricStrain(expanded.layout,expanded.grid,strain),{0.,0.},false,0.);const auto strain_snapshot=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,strain_state,request,options);for(const auto& point:strain_snapshot.Points()){for(double curl:point.vorticity_per_s)CheckClose(curl,0.,"symmetric-strain curl mismatch",kFieldTolerance);CheckClose(point.q_criterion_per_s2,-strain*strain,"symmetric-strain Q mismatch",kFieldTolerance);}

		auto compact=Make(CutCellVolumeQuadratureStorageMode::Compact);const auto compact_state=ImmersedGlobalFlowState(0.,4,compact.layout,Rotation(compact.layout,compact.grid,omega),{0.,0.},false,0.);const auto compact_snapshot=MovingImmersedFlowSnapshot::Build(*compact.geometry,compact.layout,compact_state,Request(compact),options);CheckEquivalentMetrics(metrics,compact_snapshot.Metrics());Check(snapshot.ContentHashSha256()!=compact_snapshot.ContentHashSha256(),"storage-mode-bound content hashes must differ");
		MovingImmersedFlowSnapshotOptions indicator_options;indicator_options.stagnant_speed_threshold_m_per_s=.55;const auto indicator_state=ImmersedGlobalFlowState(0.,4,expanded.layout,SignThresholdField(expanded.layout,expanded.grid),{0.,0.},false,0.);const auto indicator=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,indicator_state,Request(expanded),indicator_options);const auto compact_indicator_state=ImmersedGlobalFlowState(0.,4,compact.layout,SignThresholdField(compact.layout,compact.grid),{0.,0.},false,0.);const auto compact_indicator=MovingImmersedFlowSnapshot::Build(*compact.geometry,compact.layout,compact_indicator_state,Request(compact),indicator_options);CheckIndicatorEstimates(indicator.Metrics());CheckIndicatorEstimates(compact_indicator.Metrics());Check(indicator.Metrics().q_positive_volume_fraction>0.&&indicator.Metrics().q_positive_volume_fraction<1.,"sign-changing Q indicator is not nontrivial");Check(indicator.Metrics().stagnant_volume_fraction>0.&&indicator.Metrics().stagnant_volume_fraction<1.,"threshold indicator is not nontrivial");

		const auto zero_state=ImmersedGlobalFlowState(0.,4,expanded.layout,Constant(expanded.layout,{{0.,0.,0.}}),{0.,0.},false,0.);auto zero_request=Request(expanded);zero_request.port_flows={{7,0.},{9,0.}};const auto zero=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,options);Check(zero.Metrics().stagnant_volume_fraction==0.&&zero.Metrics().mean_q_criterion_per_s2==0.,"strict zero-field edges failed");Check(zero.Metrics().endpoint_turnover_rate_per_s==0.&&!zero.Metrics().endpoint_turnover_time_s&&zero.Metrics().well_mixed_replacement_fraction_over_step==0.,"zero-inflow transition metrics failed");options.stagnant_speed_threshold_m_per_s=1.;const auto stagnant=MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,options);Check(stagnant.Metrics().stagnant_volume_fraction==1.,"strict stagnant threshold failed");Check(stagnant.Metrics().wall_relative_velocity_rms_m_per_s<=256.*std::numeric_limits<double>::epsilon()&&stagnant.Metrics().wall_relative_velocity_max_m_per_s<=256.*std::numeric_limits<double>::epsilon(),"zero wall-relative velocity failed");

		auto translated=MakeTranslated();auto translated_request=Request(translated);const auto rigid_state=ImmersedGlobalFlowState(1.,5,translated.layout,Constant(translated.layout,{{.05,0.,0.}}),{0.,0.},false,0.);const auto rigid=MovingImmersedFlowSnapshot::Build(*translated.geometry,translated.layout,rigid_state,translated_request,options);Check(rigid.Metrics().wall_relative_velocity_rms_m_per_s<=256.*std::numeric_limits<double>::epsilon()&&rigid.Metrics().wall_relative_velocity_max_m_per_s<=256.*std::numeric_limits<double>::epsilon(),"rigid wall-relative velocity failed");constexpr double delta=.125;const auto mismatch_state=ImmersedGlobalFlowState(1.,5,translated.layout,Constant(translated.layout,{{.05,delta,0.}}),{0.,0.},false,0.);const auto mismatch=MovingImmersedFlowSnapshot::Build(*translated.geometry,translated.layout,mismatch_state,translated_request,options);const double wall_area=ExpectedWallArea(*translated.geometry,1);CheckClose(mismatch.Metrics().wall_area_m2,wall_area,"wall area mismatch");CheckClose(mismatch.Metrics().wall_relative_velocity_squared_area_integral_m4_per_s2,delta*delta*wall_area,"wall squared-area integral mismatch");CheckClose(mismatch.Metrics().wall_relative_velocity_rms_m_per_s,std::abs(delta),"wall RMS mismatch");CheckClose(mismatch.Metrics().wall_relative_velocity_max_m_per_s,std::abs(delta),"wall maximum mismatch");

		ExpectUnpublished([&]{auto bad=options;bad.maximum_points=0;(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,bad);},"moving immersed snapshot options are invalid");ExpectUnpublished([&]{auto bad=options;bad.maximum_output_bytes=1;(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,bad);},"moving immersed snapshot output byte cap exceeded");ExpectUnpublished([&]{auto bad=options;bad.maximum_points=1;(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,bad);},"moving immersed snapshot point cap exceeded");ExpectUnpublished([&]{auto bad=options;bad.stagnant_speed_threshold_m_per_s=-1.;(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,bad);},"moving immersed snapshot options are invalid");ExpectUnpublished([&]{auto bad=options;bad.stagnant_speed_threshold_m_per_s=std::numeric_limits<double>::quiet_NaN();(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,zero_state,zero_request,bad);},"moving immersed snapshot options are invalid");ExpectUnpublished([&]{auto bad=Request(expanded);bad.port_flows={{9,.1},{7,.1}};(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,bad,options);},"moving immersed snapshot port labels are invalid");ExpectUnpublished([&]{auto bad=Request(expanded);bad.port_flows[0].outward_flow_m3_s=std::numeric_limits<double>::quiet_NaN();(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,bad,options);},"moving immersed snapshot port labels are invalid");ExpectUnpublished([&]{auto bad=Request(expanded);bad.dt_s=0.;(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,bad,options);},"moving immersed snapshot transition is invalid");ExpectUnpublished([&]{auto bad=Request(expanded);bad.wall_labels={1,999};(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,rotation_state,bad,options);},"moving immersed snapshot requested wall label has no finite positive surface area");
		auto stale_layout=ImmersedActiveLayout::Build(expanded.geometry->Domain(),expanded.geometry->Volume(),"stale",{7,9});auto stale_state=ImmersedGlobalFlowState(0.,4,stale_layout,Constant(stale_layout,{{0.,0.,0.}}),{0.,0.},false,0.);ExpectUnpublished([&]{(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,stale_layout,stale_state,zero_request,options);},"moving immersed snapshot geometry, layout, and state identities do not match");auto stale_time=ImmersedGlobalFlowState(.125,4,expanded.layout,Constant(expanded.layout,{{0.,0.,0.}}),{0.,0.},false,0.);ExpectUnpublished([&]{(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,stale_time,zero_request,options);},"moving immersed snapshot time or index does not match state");const auto gauge_layout=ImmersedActiveLayout::Build(expanded.geometry->Domain(),expanded.geometry->Volume(),expanded.geometry->GeometryIdentitySha256(),{7,9},true);const auto gauge_state=ImmersedGlobalFlowState(0.,4,gauge_layout,Constant(gauge_layout,{{0.,0.,0.}}),{0.,0.},true,0.);ExpectUnpublished([&]{(void)MovingImmersedFlowSnapshot::Build(*expanded.geometry,expanded.layout,gauge_state,zero_request,options);},"moving immersed snapshot geometry, layout, and state identities do not match");

		std::cout<<std::setprecision(17)<<"moving immersed snapshot evidence: points="<<snapshot.Points().size()<<" volume="<<metrics.quadrature_volume_m3<<" mean_q="<<metrics.mean_q_criterion_per_s2<<" mean_enstrophy="<<metrics.mean_enstrophy_per_s2<<" compact_tolerance="<<kMetricTolerance<<" content_hash="<<snapshot.ContentHashSha256()<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}

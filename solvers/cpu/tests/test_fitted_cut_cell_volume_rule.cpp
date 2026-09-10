// Reuse the analytic nonconvex fixture and its scalar moment oracle.
#define main PolyhedralMomentOracleMain
#include "test_polyhedral_volume_moments.cpp"
#undef main
#include "FittedCutCellVolumeRule.hpp"
#include "CutCellVolumeQuadrature.hpp"
#include <cstring>

int main()
{
	try {
		iga::SurfaceSpatialIndex surface(iga::ClosedTriangulatedSurface::Build(ConcavePrism()));
		const std::array<double,3> lower{{0,0,0}},upper{{2,2,1}};
		std::vector<iga::VolumeQuadraturePoint> seed_points;
		const auto gauss=iga::polyhedral_moment_detail::Gauss(8);
		for(const auto& x:gauss)for(const auto& y:gauss)for(const auto& z:gauss) {
			const std::array<double,3> point{{static_cast<double>(x.first),static_cast<double>(y.first),static_cast<double>(z.first)}};
			if(surface.LocatePoint({{2*point[0],2*point[1],point[2]}})==iga::PointLocation::Inside)
				seed_points.push_back({point,static_cast<double>(x.second*y.second*z.second)});
		}
		const iga::VolumeQuadratureRule seed(seed_points);
		const auto result=iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed);
		Check(!result.rule.Points().empty()&&result.relative_moment_residual<1e-13L,"nonconvex fitted rule failed");
		iga::CompactCutCellVolumeRule compact;compact.fitted_points=result.rule.Points();
		iga::ValidateCompactCutCellVolumeRule(compact);
		Check(iga::CompactCutCellVolumeLogicalPointCount(compact)==result.rule.Points().size(),"fitted compact point count differs");
		std::size_t emitted=0;
		iga::ForEachVolumePoint(compact,[&](const iga::VolumeQuadraturePoint& point){
			const auto& original=result.rule.Points().at(emitted++);
			Check(std::memcmp(point.parametric.data(),original.parametric.data(),3*sizeof(double))==0
				&&std::memcmp(&point.weight,&original.weight,sizeof(double))==0,"fitted compact point bits or ordering changed");
		});
		Check(emitted==result.rule.Points().size(),"fitted compact iterator lost points");
		Reject([&]{auto invalid=compact;invalid.certified_blocks.push_back({{{0,0,0}},{{1,1,1}}});iga::ValidateCompactCutCellVolumeRule(invalid);});
		Reject([&]{auto invalid=compact;invalid.fitted_points[0].weight=0;iga::ForEachVolumePoint(invalid,[](const auto&){});});
		Reject([&]{auto invalid=compact;invalid.fitted_points[0].parametric[0]=std::numeric_limits<double>::quiet_NaN();iga::CompactCutCellVolumeLogicalPointCount(invalid);});
		Reject([&]{auto invalid=compact;invalid.fitted_points[1]=invalid.fitted_points[0];iga::ValidateCompactCutCellVolumeRule(invalid);});
		Reject([&]{auto invalid=compact;std::reverse(invalid.fitted_points.begin(),invalid.fitted_points.end());iga::ValidateCompactCutCellVolumeRule(invalid);});
		for(const auto& point:result.rule.Points()) {
			Check(std::isfinite(point.weight)&&point.weight>0,"nonpositive fitted rule weight");
			Check(surface.LocatePoint({{2*point.parametric[0],2*point.parametric[1],point.parametric[2]}})==iga::PointLocation::Inside,"fitted point is not inside the nonconvex surface");
		}
		long double maximum_error=0;
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			long double actual=0;
			for(const auto& point:result.rule.Points())actual+=point.weight*iga::polyhedral_moment_detail::Power(point.parametric[0],a)
				*iga::polyhedral_moment_detail::Power(point.parametric[1],b)*iga::polyhedral_moment_detail::Power(point.parametric[2],c);
			const long double x=std::ldexp(1.L,a+1),y=std::ldexp(1.L,b+1);
			const long double expected=(x*y-(x-1)*(y-1))/((a+1)*(b+1)*(c+1))/std::ldexp(4.L,a+b);
			maximum_error=std::max(maximum_error,std::abs(actual-expected));
			Check(std::abs(actual-expected)<2e-13L,"fitted rule differs from analytic nonconvex moments");
			Check(std::abs(iga::CompactCutCellVolumeMoment(compact,a,b,c)-expected)<2e-13L,"compact fitted rule differs from analytic nonconvex moments");
		}
		Reject([&]{auto options=iga::FittedCutCellVolumeRuleOptions();options.max_point_queries=1;iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed,options);});
		Reject([&]{auto options=iga::FittedCutCellVolumeRuleOptions();options.max_seed_points=1;iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed,options);});
		Reject([&]{auto options=iga::FittedCutCellVolumeRuleOptions();options.fit.max_columns=1;iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed,options);});
		Reject([&]{auto options=iga::FittedCutCellVolumeRuleOptions();options.candidate_orders={2};iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed,options);});
		Reject([&]{auto options=iga::FittedCutCellVolumeRuleOptions();options.candidate_orders={16,8};iga::BuildFittedCutCellVolumeRule(surface,lower,upper,seed,options);});
		Reject([&]{iga::BuildFittedCutCellVolumeRule(surface,lower,upper,iga::VolumeQuadratureRule());});
		std::cout<<"fitted_cut_cell_volume_rule_test: PASS nonconvex geometry, 343 analytic moments, caps and exhausted candidates; max_error="<<maximum_error<<'\n';
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}

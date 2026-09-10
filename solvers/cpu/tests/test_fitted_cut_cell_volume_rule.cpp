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
		Check(result.candidates<=iga::FittedCutCellVolumeRuleOptions{}.fit.max_columns,"fitted candidate reservoir exceeded its matrix cap");
		// One interior sample has no coordinate spread, but this nonconvex
		// domain still has positive volume and must be fitted, not rejected.
		const iga::VolumeQuadratureRule point_seed({{{{.25,.25,.5}},.75}});
		const auto point_result=iga::BuildFittedCutCellVolumeRule(surface,lower,upper,point_seed);
		Check(point_result.relative_moment_residual<1e-13L,"zero-spread seed moment audit failed");
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			long double measured=0;
			for(const auto& point:point_result.rule.Points()) {
				Check(point.weight>0,"zero-spread fit has nonpositive weight");
				measured+=point.weight*std::pow(static_cast<long double>(point.parametric[0]),a)
					*std::pow(static_cast<long double>(point.parametric[1]),b)*std::pow(static_cast<long double>(point.parametric[2]),c);
			}
			const long double expected=(1-(1-std::pow(.5L,a+1))*(1-std::pow(.5L,b+1)))/((a+1.L)*(b+1.L)*(c+1.L));
			Check(std::abs(measured-expected)<2e-13L,"zero-spread fitted moment differs from analytic prism");
		}
		Reject([&]{auto limited=iga::FittedCutCellVolumeRuleOptions();limited.max_point_queries=1;
			iga::BuildFittedCutCellVolumeRule(surface,lower,upper,point_seed,limited);});

		// Reproduce the moving endpoint's partial cell without changing its
		// clipped volume: the true local box is [0.8,1] x [0.25,1]^2.
		auto shifted=ConcavePrism();for(auto& point:shifted.vertices) {
			point[0]+=.84;point[1]+=.1;point[2]+=.1;
		}
		iga::SurfaceSpatialIndex corner(iga::ClosedTriangulatedSurface::Build(shifted));
		const iga::VolumeQuadratureRule corner_seed({
			{{{.83250236955189294,.26735796105074344,.26735796105074344}},1.},
			{{{.9826420389492565,.9826420389492565,.9826420389492565}},1.}});
		iga::FittedCutCellVolumeRuleOptions corner_options;corner_options.support_expansion=3.;
		const auto corner_rule=iga::BuildFittedCutCellVolumeRule(corner,{{.6,0,0}},{{.9,.4,.4}},corner_seed,corner_options);
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			long double measured=0;
			for(const auto& point:corner_rule.rule.Points()) {
				Check(point.weight>0&&point.parametric[0]>=.8-1e-14&&point.parametric[0]<=1
					&&point.parametric[1]>=.25-1e-14&&point.parametric[2]>=.25-1e-14,"corner rule point or weight is invalid");
				measured+=point.weight*std::pow(static_cast<long double>(point.parametric[0]),a)
					*std::pow(static_cast<long double>(point.parametric[1]),b)*std::pow(static_cast<long double>(point.parametric[2]),c);
			}
			const long double expected=(1-std::pow(.8L,a+1))/(a+1)*(1-std::pow(.25L,b+1))/(b+1)*(1-std::pow(.25L,c+1))/(c+1);
			Check(std::abs(measured-expected)<2e-13L,"clipped corner tensor moment differs from analytic box");
		}
		auto thin_soup=ConcavePrism();for(auto& point:thin_soup.vertices) {
			point[0]+=.285;point[1]+=.1;point[2]+=.1;
		}
		iga::SurfaceSpatialIndex thin_surface(iga::ClosedTriangulatedSurface::Build(thin_soup));
		const iga::VolumeQuadratureRule thin_seed({
			{{{.9826420389492565,.26735796105074344,.26735796105074344}},1.},
			{{{.9826420389492565,.9826420389492565,.9826420389492565}},1.}});
		const auto thin_rule=iga::BuildFittedCutCellVolumeRule(thin_surface,{{0,0,0}},{{.3,.4,.4}},thin_seed,corner_options);
		for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
			long double measured=0;
			for(const auto& point:thin_rule.rule.Points()) {
				Check(point.weight>0&&point.parametric[0]>=.95-1e-14,"thin fitted point or weight is invalid");
				measured+=point.weight*std::pow(static_cast<long double>(point.parametric[0]),a)
					*std::pow(static_cast<long double>(point.parametric[1]),b)*std::pow(static_cast<long double>(point.parametric[2]),c);
			}
			const long double expected=(1-std::pow(.95L,a+1))/(a+1)*(1-std::pow(.25L,b+1))/(b+1)*(1-std::pow(.25L,c+1))/(c+1);
			Check(std::abs(measured-expected)<2e-13L,"zero-spread thin cut moment differs from analytic box");
		}

		iga::CompactCutCellVolumeRule compact;compact.fitted_points=result.rule.Points();
		iga::ValidateCompactCutCellVolumeRule(compact);
		Check(iga::CompactCutCellVolumeRecordCount(compact)==compact.fitted_points.size(),"compact fitted records were not counted");
		Check(iga::CompactCutCellVolumeCapacityBytes(compact)==compact.fitted_points.capacity()*sizeof(iga::VolumeQuadraturePoint),"compact fitted storage was not counted");
		Reject([&]{iga::CompactCutCellVolumeCapacityBytes(std::numeric_limits<std::size_t>::max(),0,0);});
		const auto digest=[](const iga::CompactCutCellVolumeRule& rule){iga::Sha256 hash;iga::AppendCompactCutCellVolumeRuleHash(hash,rule);return hash.Hex();};
		auto changed=compact;changed.fitted_points[0].weight=std::nextafter(changed.fitted_points[0].weight,1.);
		Check(digest(changed)!=digest(compact),"compact fitted weight was not bound by hash");
		changed=compact;changed.fitted_points[0].parametric[0]=std::nextafter(changed.fitted_points[0].parametric[0],1.);
		Check(digest(changed)!=digest(compact),"compact fitted coordinate was not bound by hash");
		iga::CompactCutCellVolumeRule legacy;legacy.max_depth=1;legacy.certified_blocks.push_back({{{0,0,0}},{{2,2,2}}});
		iga::Sha256 legacy_bytes;legacy_bytes.AppendLittleEndian32(1);legacy_bytes.AppendLittleEndian64(1);
		for(unsigned q=0;q<3;++q)legacy_bytes.AppendLittleEndian32(0);
		for(unsigned q=0;q<3;++q)legacy_bytes.AppendLittleEndian32(2);
		legacy_bytes.AppendLittleEndian64(0);
		Check(digest(legacy)==legacy_bytes.Hex(),"legacy compact hash byte stream changed");
		legacy.certified_blocks.clear();legacy.sample_leaves.push_back({{{1,0,1}},1,0x55u});
		iga::Sha256 sample_bytes;sample_bytes.AppendLittleEndian32(1);sample_bytes.AppendLittleEndian64(0);sample_bytes.AppendLittleEndian64(1);
		sample_bytes.AppendLittleEndian32(1);sample_bytes.AppendLittleEndian32(0);sample_bytes.AppendLittleEndian32(1);
		sample_bytes.AppendLittleEndian32(1);sample_bytes.AppendLittleEndian64(0x55u);
		Check(digest(legacy)==sample_bytes.Hex(),"legacy compact sample hash byte stream changed");
		Check(iga::CompactCutCellVolumeRecordCount(legacy)==1,"legacy sample record count changed");
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
		iga::CartesianDomainClassification domain(iga::CubicCartesianBackground({lower,upper,{{1,1,1}}}),
			iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(ConcavePrism())));
		iga::OctreeCutQuadratureOptions quadrature;quadrature.max_depth=2;
		iga::FittedCutCellVolumeRuleOptions fitting;fitting.fit.max_columns=32768;
		const iga::CutCellVolumeQuadratureCatalog seed_catalog(domain,quadrature,iga::CutCellVolumeQuadratureStorageMode::Compact);
		const iga::CutCellVolumeQuadratureCatalog expanded(domain,quadrature,iga::CutCellVolumeQuadratureStorageMode::Expanded,&fitting);
		auto packed_options=quadrature;packed_options.max_points=1;
		const iga::CutCellVolumeQuadratureCatalog packed(domain,packed_options,iga::CutCellVolumeQuadratureStorageMode::Compact,&fitting);
		expanded.ValidateUsableRule(domain,0);packed.ValidateUsableCompactRule(domain,0);
		Check(expanded.Cell(0).moment_fitted&&packed.Cell(0).moment_fitted,"catalog did not publish fitted rule");
		Check(std::abs(expanded.Cell(0).diagnostics.estimated_physical_volume-3)<2e-13,"fitted catalog volume is inaccurate");
		const auto& published=expanded.UsableRule(domain,0).Points();emitted=0;
		iga::ForEachVolumePoint(packed.UsableCompactRule(domain,0),[&](const iga::VolumeQuadraturePoint& point){
			const auto& original=published.at(emitted++);
			Check(std::memcmp(point.parametric.data(),original.parametric.data(),3*sizeof(double))==0
				&&std::memcmp(&point.weight,&original.weight,sizeof(double))==0,"catalog storage modes produced different fitted points");
		});
		Check(emitted==published.size(),"catalog storage point counts differ");
		Check(packed.Cell(0).rule.Points().empty(),"compact fitted catalog retained expanded points");
		Check(packed.Cell(0).diagnostics.observed_retained_bytes==iga::CompactCutCellVolumeCapacityBytes(packed.Cell(0).compact_rule),"fitted catalog capacity was not recorded");
		Reject([&]{auto limited=quadrature;limited.max_points=1;iga::CutCellVolumeQuadratureCatalog failed(domain,limited,iga::CutCellVolumeQuadratureStorageMode::Expanded,&fitting);});
		Reject([&]{auto limited=fitting;limited.max_point_queries=0;iga::CutCellVolumeQuadratureCatalog failed(domain,quadrature,iga::CutCellVolumeQuadratureStorageMode::Compact,&limited);});
		Reject([&]{auto limited=fitting;limited.max_seed_points=1;iga::CutCellVolumeQuadratureCatalog failed(domain,quadrature,iga::CutCellVolumeQuadratureStorageMode::Compact,&limited);});
		// A valid but exhausted query budget must fail during fitting, after
		// seed construction, without publishing the seed as a fallback rule.
		for(const auto mode:{iga::CutCellVolumeQuadratureStorageMode::Expanded,iga::CutCellVolumeQuadratureStorageMode::Compact}) {
			const auto& diagnostics=(mode==iga::CutCellVolumeQuadratureStorageMode::Expanded?expanded:packed).Cell(0).diagnostics;
			const auto& seed_diagnostics=seed_catalog.Cell(0).diagnostics;
			Check(diagnostics.attempted_record_attempts==seed_diagnostics.attempted_record_attempts+published.size(),"fitted publication omitted cumulative records");
			Check(diagnostics.attempted_logical_output_points==seed_diagnostics.attempted_logical_output_points+published.size(),"fitted publication omitted cumulative logical points");
			Check(diagnostics.attempted_output_points==(mode==iga::CutCellVolumeQuadratureStorageMode::Expanded?published.size():0),"fitted publication expanded work count is incorrect");
			Reject([&]{auto limited=quadrature;limited.max_logical_points=seed_diagnostics.attempted_logical_output_points+published.size()-1;iga::CutCellVolumeQuadratureCatalog failed(domain,limited,mode,&fitting);});
			Reject([&]{auto limited=fitting;limited.max_point_queries=1;iga::CutCellVolumeQuadratureCatalog failed(domain,quadrature,mode,&limited);});
			Reject([&]{auto limited=quadrature;limited.max_retained_bytes=2*343*sizeof(iga::VolumeQuadraturePoint)-1;iga::CutCellVolumeQuadratureCatalog failed(domain,limited,mode,&fitting);});
			// This cell touches the prism at its x=2 face but has no interior
			// intersection. Its certified empty seed must remain empty.
			iga::CartesianDomainClassification tangent_domain(iga::CubicCartesianBackground({{{0,0,0}},{{3,3,1}},{{3,3,1}}}),
				iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(ConcavePrism())));
			const iga::CutCellVolumeQuadratureCatalog tangent(tangent_domain,quadrature,mode,&fitting);
			const auto& empty=tangent.Cell(2);
			Check(empty.classification==iga::CellClassification::Cut&&empty.usable&&!empty.moment_fitted,"tangent fitted catalog classification changed");
			Check(empty.rule.Points().empty()&&iga::CompactCutCellVolumeLogicalPointCount(empty.compact_rule)==0,"tangent fitted catalog is not empty");
			Check(empty.diagnostics.estimated_physical_volume==0&&empty.diagnostics.unresolved_physical_volume==0&&empty.fitting_queries==0,"tangent fitted catalog performed fitting or retained unresolved volume");
			if(mode==iga::CutCellVolumeQuadratureStorageMode::Expanded)tangent.ValidateUsableRule(tangent_domain,2);
			else tangent.ValidateUsableCompactRule(tangent_domain,2);
			auto small=ConcavePrism();
			for(auto& point:small.vertices)for(auto& value:point)value=.49+.001*value;
			iga::CartesianDomainClassification small_domain(iga::CubicCartesianBackground({{{0,0,0}},{{1,1,1}},{{1,1,1}}}),
				iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(small)));
			auto shallow=quadrature;shallow.max_depth=0;shallow.empty_rule_rescue_max_depth=0;
			const iga::CutCellVolumeQuadratureCatalog unresolved(small_domain,shallow,mode);
			Check(unresolved.Cell(0).diagnostics.logical_output_points==0&&unresolved.Cell(0).diagnostics.unresolved_reference_volume>0,"small fixture is not provisionally empty");
			Reject([&]{iga::CutCellVolumeQuadratureCatalog failed(small_domain,shallow,mode,&fitting);});
		}
		std::cout<<"fitted_cut_cell_volume_rule_test: PASS nonconvex geometry, 343 analytic moments, caps and exhausted candidates; max_error="<<maximum_error<<'\n';
		std::cout<<"fitted_catalog: PASS expanded/compact point identity, physical volume, diagnostics and caps\n";
	} catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
	return 0;
}

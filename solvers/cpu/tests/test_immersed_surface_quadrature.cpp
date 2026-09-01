#include "ImmersedSurfaceQuadrature.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t label=7)
{ iga::RawSurfaceTriangle t; t.indices={{a,b,c}}; t.boundary_id=label; return t; }
iga::RawSurfaceSoup Cube(double lo, double hi, std::int64_t label=7)
{
	iga::RawSurfaceSoup s; s.vertices={{{{lo,lo,lo}},{{hi,lo,lo}},{{hi,hi,lo}},{{lo,hi,lo}},{{lo,lo,hi}},{{hi,lo,hi}},{{hi,hi,hi}},{{lo,hi,hi}}}};
	s.triangles={Face(0,2,1,label),Face(0,3,2,label),Face(4,5,6,label),Face(4,6,7,label),Face(0,1,5,label),Face(0,5,4,label),Face(1,2,6,label),Face(1,6,5,label),Face(2,3,7,label),Face(2,7,6,label),Face(3,0,4,label),Face(3,4,7,label)}; return s;
}
iga::RawSurfaceSoup Tetra(double scale=1.0, std::array<double,3> offset={{0,0,0}})
{
	iga::RawSurfaceSoup s; s.vertices={{{{offset[0],offset[1],offset[2]}},{{offset[0]+scale,offset[1],offset[2]}},{{offset[0],offset[1]+scale,offset[2]}},{{offset[0],offset[1],offset[2]+scale}}}};
	s.triangles={Face(0,2,1,3),Face(0,1,3,4),Face(0,3,2,5),Face(1,2,3,6)}; return s;
}
iga::RawSurfaceSoup AnisotropicTetra(std::array<double,3> origin, std::array<double,3> length)
{
	iga::RawSurfaceSoup s; s.vertices={{{{origin[0],origin[1],origin[2]}},{{origin[0]+length[0],origin[1],origin[2]}},{{origin[0],origin[1]+length[1],origin[2]}},{{origin[0],origin[1],origin[2]+length[2]}}}};
	s.triangles={Face(0,2,1,3),Face(0,1,3,4),Face(0,3,2,5),Face(1,2,3,6)}; return s;
}
iga::CartesianDomainClassification Domain(iga::RawSurfaceSoup s, iga::CubicCartesianGridSpec spec)
{ return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec),iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(std::move(s)))); }
bool Near(double a,double b,double tol=3e-10) { return std::abs(a-b)<=tol*std::max({1.0,std::abs(a),std::abs(b)}); }
bool AreaNear(double a, double b, double factor=1024.0)
{ return std::abs(a-b) <= factor*std::numeric_limits<double>::epsilon()*(std::abs(a)+std::abs(b)); }
template<class F> void Reject(F&& f) { bool yes=false; try { f(); } catch(const std::exception&) { yes=true; } assert(yes); }
double Moment(const iga::ImmersedSurfaceQuadratureCatalog& q,const iga::CartesianDomainClassification& domain,int x,int y,int z)
{ long double r=0; for(const auto& c:q.Cells()) for(const auto& p:q.UsableRule(domain,c.id).Points()) r+=static_cast<long double>(p.weight)*std::pow(p.physical[0],x)*std::pow(p.physical[1],y)*std::pow(p.physical[2],z); return static_cast<double>(r); }
void AssertAudit(const iga::ImmersedSurfaceQuadratureCatalog& q, double expected)
{
	const auto& d=q.Diagnostics(); assert(q.Usable());
	assert(AreaNear(d.source_total_area_m2,expected) && AreaNear(d.accumulated_total_area_m2,expected));
	assert(AreaNear(d.total_area_m2,d.accumulated_total_area_m2));
	assert(std::abs(d.total_area_residual_m2)<=1024.0*std::numeric_limits<double>::epsilon()*expected);
	assert(std::abs(d.total_area_absolute_residual_m2)<=1024.0*std::numeric_limits<double>::epsilon()*expected);
	assert(d.triangle_area_m2.size()==d.triangle_area_residual_m2.size());
	assert(d.triangle_area_m2.size()==d.triangle_area_absolute_residual_m2.size());
	for(std::size_t i=0;i<d.triangle_area_m2.size();++i) {
		assert(std::isfinite(d.triangle_area_residual_m2[i])); assert(d.triangle_area_absolute_residual_m2[i]>=0.0);
		assert(AreaNear(d.triangle_area_absolute_residual_m2[i],std::abs(d.triangle_area_residual_m2[i])));
	}
	for(const auto& item:d.area_by_boundary_id) {
		const auto label=item.first; assert(d.source_area_by_boundary_id.count(label));
		assert(d.area_residual_by_boundary_id.count(label) && d.area_absolute_residual_by_boundary_id.count(label));
		assert(std::abs(d.area_residual_by_boundary_id.at(label))<=1024.0*std::numeric_limits<double>::epsilon()*d.source_area_by_boundary_id.at(label));
		assert(std::abs(d.area_absolute_residual_by_boundary_id.at(label))<=1024.0*std::numeric_limits<double>::epsilon()*d.source_area_by_boundary_id.at(label));
	}
}
void AssertSameRule(const iga::ImmersedSurfaceQuadratureCatalog& first, const iga::CartesianDomainClassification& first_domain, const iga::ImmersedSurfaceQuadratureCatalog& second, const iga::CartesianDomainClassification& second_domain, double physical_scale=1.0)
{
	assert(first.Cells().size()==second.Cells().size());
	for(std::size_t cell=0;cell<first.Cells().size();++cell) {
		const auto& a=first.Cells()[cell]; const auto& b=second.Cells()[cell]; const auto& ar=first.UsableRule(first_domain,a.id); const auto& br=second.UsableRule(second_domain,b.id); assert(a.id==b.id && a.usable==b.usable && a.point_count==b.point_count);
		for(std::size_t point=0;point<ar.Points().size();++point) {
			const auto& x=ar.Points()[point]; const auto& y=br.Points()[point]; assert(x.boundary_id==y.boundary_id);
			for(int axis=0;axis<3;++axis) { assert(Near(x.parametric[axis],y.parametric[axis])); assert(Near(x.normal[axis],y.normal[axis])); assert(Near(y.physical[axis],physical_scale*x.physical[axis])); }
			assert(Near(y.weight,physical_scale*physical_scale*x.weight));
		}
	}
}
bool SameDoubleBits(double a, double b)
{
	std::uint64_t left=0, right=0; std::memcpy(&left,&a,sizeof(left)); std::memcpy(&right,&b,sizeof(right)); return left==right;
}
void AssertBitIdenticalRule(const iga::ImmersedSurfaceQuadratureCatalog& first, const iga::CartesianDomainClassification& first_domain,
	const iga::ImmersedSurfaceQuadratureCatalog& second, const iga::CartesianDomainClassification& second_domain)
{
	assert(first.Cells().size()==second.Cells().size());
	for(std::size_t cell=0;cell<first.Cells().size();++cell) {
		const auto& a=first.Cells()[cell]; const auto& b=second.Cells()[cell];
		assert(a.id==b.id && a.usable==b.usable && a.ambiguous==b.ambiguous && SameDoubleBits(a.area_m2,b.area_m2) && a.point_count==b.point_count);
		const auto& ar=first.UsableRule(first_domain,a.id); const auto& br=second.UsableRule(second_domain,b.id); assert(ar.Points().size()==br.Points().size());
		for(std::size_t point=0;point<ar.Points().size();++point) {
			const auto& x=ar.Points()[point]; const auto& y=br.Points()[point]; assert(x.boundary_id==y.boundary_id && SameDoubleBits(x.weight,y.weight));
			for(int axis=0;axis<3;++axis) assert(SameDoubleBits(x.parametric[axis],y.parametric[axis]) && SameDoubleBits(x.physical[axis],y.physical[axis]) && SameDoubleBits(x.normal[axis],y.normal[axis]));
		}
	}
}
void AssertTranslationInvariantRule(const iga::ImmersedSurfaceQuadratureCatalog& local, const iga::CartesianDomainClassification& local_domain,
	const iga::ImmersedSurfaceQuadratureCatalog& translated, const iga::CartesianDomainClassification& translated_domain,
	const std::array<double,3>& translated_origin)
{
	assert(local.Cells().size()==translated.Cells().size());
	for(std::size_t cell=0;cell<local.Cells().size();++cell) {
		const auto& a=local.Cells()[cell]; const auto& b=translated.Cells()[cell];
		assert(a.id==b.id && a.usable==b.usable && a.ambiguous==b.ambiguous && SameDoubleBits(a.area_m2,b.area_m2) && a.point_count==b.point_count);
		const auto& ar=local.UsableRule(local_domain,a.id); const auto& br=translated.UsableRule(translated_domain,b.id); assert(ar.Points().size()==br.Points().size());
		for(std::size_t point=0;point<ar.Points().size();++point) {
			const auto& x=ar.Points()[point]; const auto& y=br.Points()[point]; assert(x.boundary_id==y.boundary_id && SameDoubleBits(x.weight,y.weight));
			bool on_face=false;
			for(int axis=0;axis<3;++axis) {
				assert(SameDoubleBits(x.parametric[axis],y.parametric[axis]) && SameDoubleBits(x.normal[axis],y.normal[axis]));
				const double recovered_local=y.physical[axis]-translated_origin[axis];
				assert(std::abs(recovered_local-x.physical[axis])<=4.0*std::max(std::abs(std::nextafter(y.physical[axis],std::numeric_limits<double>::infinity())-y.physical[axis]),std::numeric_limits<double>::denorm_min()));
				on_face=on_face || x.parametric[axis]==0.0 || x.parametric[axis]==1.0;
			}
			if (!on_face) for(int axis=0;axis<3;++axis) assert(x.parametric[axis]>0.0 && x.parametric[axis]<1.0 && y.parametric[axis]>0.0 && y.parametric[axis]<1.0);
		}
	}
}
#ifdef IGA_EXACT_DYADIC_TESTING
struct ResetTestControls {
	~ResetTestControls() { iga::SetSurfaceQuadratureTestRejectPositiveArea(false); iga::SetSurfaceQuadratureTestRejectRuleValidation(false); iga::exact_dyadic::SetTestMagnitudeCap(iga::exact_dyadic::kMaximumMagnitudeLimbs); }
};
#endif
}

int main()
{
#ifdef IGA_EXACT_DYADIC_TESTING
	iga::exact_dyadic::SetTestMagnitudeCap(iga::exact_dyadic::kMaximumMagnitudeLimbs);
#endif
	const iga::CubicCartesianGridSpec root{{{0,0,0}},{{1,1,1}},{{1,1,1}}};
	const auto cube_domain=Domain(Cube(0,1),root); const iga::ImmersedSurfaceQuadratureCatalog cube(cube_domain);
	assert(cube.Cell(0).usable && cube.Cell(0).point_count==144); AssertAudit(cube,6.0); assert(Near(cube.Diagnostics().area_by_boundary_id.at(7),6));
	assert(Near(Moment(cube,cube_domain,1,0,0),3) && Near(Moment(cube,cube_domain,0,1,0),3) && Near(Moment(cube,cube_domain,0,0,1),3));
	assert(Near(Moment(cube,cube_domain,2,0,0),7.0/3.0) && Near(Moment(cube,cube_domain,3,0,0),2.0) && Near(Moment(cube,cube_domain,6,0,0),11.0/7.0));
	assert(Near(Moment(cube,cube_domain,2,1,0),7.0/6.0));
	assert(Near(Moment(cube,cube_domain,2,2,2),1.0/3.0));
	std::array<double,3> normal_sum{}; double divergence=0; for(const auto& point:cube.UsableRule(cube_domain,0).Points()) { for(int a=0;a<3;++a) normal_sum[a]+=point.weight*point.normal[a]; divergence+=point.weight*(point.physical[0]*point.normal[0]+point.physical[1]*point.normal[1]+point.physical[2]*point.normal[2]); }
	assert(Near(normal_sum[0],0) && Near(normal_sum[1],0) && Near(normal_sum[2],0) && Near(divergence,3)); cube.ValidateUsableRule(cube_domain,0);

	const auto tetra_domain=Domain(Tetra(),root); const iga::ImmersedSurfaceQuadratureCatalog tetra_root(tetra_domain); const double tetra_area=1.5+.5*std::sqrt(3.0);
	assert(tetra_root.Cell(0).usable && tetra_root.Cell(0).point_count==48); AssertAudit(tetra_root,tetra_area);
	assert(Near(tetra_root.Diagnostics().area_by_boundary_id.at(3),.5) && Near(tetra_root.Diagnostics().area_by_boundary_id.at(4),.5) && Near(tetra_root.Diagnostics().area_by_boundary_id.at(5),.5) && Near(tetra_root.Diagnostics().area_by_boundary_id.at(6),.5*std::sqrt(3.0)));
	normal_sum={}; divergence=0; for(const auto& point:tetra_root.UsableRule(tetra_domain,0).Points()) { for(int a=0;a<3;++a) normal_sum[a]+=point.weight*point.normal[a]; divergence+=point.weight*(point.physical[0]*point.normal[0]+point.physical[1]*point.normal[1]+point.physical[2]*point.normal[2]); }
	assert(Near(normal_sum[0],0) && Near(normal_sum[1],0) && Near(normal_sum[2],0) && Near(divergence,.5));

	const iga::CubicCartesianGridSpec fine{{{0,0,0}},{{1,1,1}},{{2,2,2}}}; const auto fine_domain=Domain(Cube(0,1),fine); const iga::ImmersedSurfaceQuadratureCatalog refined(fine_domain);
	AssertAudit(refined,6.0); assert(Near(Moment(refined,fine_domain,1,0,0),3)); for (const auto& cell:refined.Cells()) refined.ValidateUsableRule(fine_domain,cell.id);
	std::array<double,3> degree8_error{};
	for(std::uint32_t n=1;n<=4;n*=2) { const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{n,n,n}}}; const auto grid_domain=Domain(Cube(0,1),grid); const iga::ImmersedSurfaceQuadratureCatalog q(grid_domain); assert(q.Usable()); degree8_error[static_cast<std::size_t>(std::log2(n))]=std::abs(Moment(q,grid_domain,8,0,0)-13.0/9.0); }
	assert(degree8_error[1]<degree8_error[0] && degree8_error[2]<degree8_error[1]);

	const iga::CubicCartesianGridSpec oblique_grid{{{0,0,0}},{{1,1,1}},{{4,4,4}}}; const auto oblique_domain=Domain(Tetra(),oblique_grid); const iga::ImmersedSurfaceQuadratureCatalog oblique(oblique_domain);
	AssertAudit(oblique,tetra_area); assert(oblique.Diagnostics().candidate_attempts>tetra_root.Diagnostics().candidate_attempts && oblique.Diagnostics().positive_fragments>tetra_root.Diagnostics().positive_fragments);
	for(const std::uint32_t n:std::array<std::uint32_t,2>{{8,16}}) {
		const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{n,n,n}}}; const auto dense_domain=Domain(Tetra(),grid); const iga::ImmersedSurfaceQuadratureCatalog dense(dense_domain,{1000000,1000000,12000000,128});
		AssertAudit(dense,tetra_area); assert(dense.Diagnostics().peak_exact_limbs<=128); for(const auto& cell:dense.Cells()) dense.ValidateUsableRule(dense_domain,cell.id);
		if(n==16) { iga::RawSurfaceSoup p=Tetra(); std::reverse(p.triangles.begin(),p.triangles.end()); for(auto& triangle:p.triangles) std::rotate(triangle.indices.begin(),triangle.indices.begin()+1,triangle.indices.end()); const auto permuted_dense_domain=Domain(p,grid); const iga::ImmersedSurfaceQuadratureCatalog permuted_dense(permuted_dense_domain,{1000000,1000000,12000000,128}); AssertBitIdenticalRule(dense,dense_domain,permuted_dense,permuted_dense_domain); }
	}

	iga::RawSurfaceSoup permuted=Tetra(); std::reverse(permuted.triangles.begin(),permuted.triangles.end()); for(auto& triangle:permuted.triangles) std::rotate(triangle.indices.begin(),triangle.indices.begin()+1,triangle.indices.end());
	const auto tetra_fine_domain=Domain(Tetra(),fine); const auto permuted_domain=Domain(permuted,fine); const iga::ImmersedSurfaceQuadratureCatalog tetra(tetra_fine_domain); const iga::ImmersedSurfaceQuadratureCatalog canonical_permuted(permuted_domain); AssertSameRule(tetra,tetra_fine_domain,canonical_permuted,permuted_domain);
	const double scale=8.0; const iga::CubicCartesianGridSpec scaled_grid{{{0,0,0}},{{scale,scale,scale}},{{2,2,2}}}; const auto scaled_domain=Domain(Tetra(scale),scaled_grid); const iga::ImmersedSurfaceQuadratureCatalog scaled(scaled_domain); AssertSameRule(tetra,tetra_fine_domain,scaled,scaled_domain,scale); AssertAudit(scaled,tetra_area*scale*scale);
	assert(Near(scaled.Diagnostics().source_total_area_m2,scale*scale*tetra.Diagnostics().source_total_area_m2)); assert(Near(scaled.Diagnostics().accumulated_total_area_m2,scale*scale*tetra.Diagnostics().accumulated_total_area_m2));

	const std::array<double,3> origin{{10,-3,7}}, length{{2,3,4}}, interior{{10.2,-2.7,7.4}}; const iga::CubicCartesianGridSpec translated_grid{{{10,-3,7}},{{12,0,11}},{{2,2,2}}};
	const auto translated_domain=Domain(AnisotropicTetra(origin,length),translated_grid); const iga::ImmersedSurfaceQuadratureCatalog translated(translated_domain);
	for(const auto& cell:translated.Cells()) for(const auto& point:translated.UsableRule(translated_domain,cell.id).Points()) { const auto geometry=iga::EvaluateElementGeometry(translated_domain.Background().MaterializeElement(cell.id),point.parametric); for(int axis=0;axis<3;++axis) assert(Near(geometry.physical[axis],point.physical[axis])); double outward=0; for(int axis=0;axis<3;++axis) outward+=point.normal[axis]*(point.physical[axis]-interior[axis]); assert(outward>0.0); }

	const double large_offset=std::ldexp(1.0,53), large_width=1024.0, large_inset=16.0;
	const iga::CubicCartesianGridSpec local_large_grid{{{0,0,0}},{{large_width,large_width,large_width}},{{1,1,1}}};
	const auto local_large_domain=Domain(Tetra(128.0,{{large_inset,large_inset,large_inset}}),local_large_grid);
	const iga::ImmersedSurfaceQuadratureCatalog local_large(local_large_domain);
	const std::array<double,3> large_origin{{large_offset+large_inset,large_offset+large_inset,large_offset+large_inset}};
	const iga::CubicCartesianGridSpec large_grid{{{large_offset,large_offset,large_offset}},{{large_offset+large_width,large_offset+large_width,large_offset+large_width}},{{1,1,1}}};
	const auto large_domain=Domain(AnisotropicTetra(large_origin,{{128,128,128}}),large_grid); const iga::ImmersedSurfaceQuadratureCatalog large(large_domain);
	AssertTranslationInvariantRule(local_large,local_large_domain,large,large_domain,{{large_offset,large_offset,large_offset}});
	assert(large.Usable() && large.Cell(0).point_count==48); AssertAudit(large,1.5*128.0*128.0+.5*std::sqrt(3.0)*128.0*128.0);
	large.ValidateUsableRule(large_domain,0);
	for(const auto& point:large.UsableRule(large_domain,0).Points()) {
		for(int axis=0;axis<3;++axis) { const double mapped=large_offset+point.parametric[axis]*large_width; assert(std::abs(mapped-point.physical[axis])<=8.0*std::max(std::abs(std::nextafter(point.physical[axis],std::numeric_limits<double>::infinity())-point.physical[axis]),std::abs(std::nextafter(mapped,std::numeric_limits<double>::infinity())-mapped))); }
		if(point.boundary_id==5) { assert(point.physical[0]==large_origin[0]); assert(point.parametric[0]==large_inset/large_width); }
	}

	const iga::CubicCartesianGridSpec plane_grid{{{0,0,0}},{{1,1,1}},{{2,2,2}}}; const auto high_domain=Domain(Cube(.5,1),plane_grid); const iga::ImmersedSurfaceQuadratureCatalog high(high_domain); AssertAudit(high,1.5); for(const auto& cell:high.Cells()) { if(cell.id==7) assert(cell.point_count>0); else assert(cell.point_count==0); }
	const auto low_domain=Domain(Cube(0,.5),plane_grid); const iga::ImmersedSurfaceQuadratureCatalog low(low_domain); AssertAudit(low,1.5); for(const auto& cell:low.Cells()) { if(cell.id==0) assert(cell.point_count>0); else assert(cell.point_count==0); }

	Reject([&]{ cube.UsableRule(Domain(Cube(0,1,8),root),0); }); Reject([&]{ cube.UsableRule(tetra_domain,0); });
	const iga::CubicCartesianGridSpec wrong_grid{{{0,0,0}},{{2,1,1}},{{1,1,1}}}; const auto wrong_domain=Domain(Cube(0,1),wrong_grid); Reject([&]{ cube.UsableRule(wrong_domain,0); });
	Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog bad_label(Domain(Cube(0,1,static_cast<std::int64_t>(std::numeric_limits<int>::max())+1),root)); });
	Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog bad(cube_domain,{1,100,100,512}); }); Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog bad(cube_domain,{100,1,100,512}); }); Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog bad(cube_domain,{100,100,1,512}); });
	const double micro_scale=std::ldexp(1.0,-200); const iga::CubicCartesianGridSpec micro_grid{{{0,0,0}},{{micro_scale,micro_scale,micro_scale}},{{1,1,1}}}; const iga::ImmersedSurfaceQuadratureCatalog micro(Domain(Cube(0,micro_scale),micro_grid)); AssertAudit(micro,6.0*micro_scale*micro_scale);

#ifdef IGA_EXACT_DYADIC_TESTING
	{ ResetTestControls reset; const iga::ImmersedSurfaceQuadratureCatalog exact_capped(fine_domain,{1000,1000,12000,1}); assert(!exact_capped.Usable() && !exact_capped.Diagnostics().catalog_unusable_reason.empty()); assert(exact_capped.Diagnostics().peak_exact_limbs>0 && exact_capped.Diagnostics().peak_exact_limbs<=1); for(const auto& cell:exact_capped.Cells()) Reject([&]{ exact_capped.UsableRule(fine_domain,cell.id); }); }
	{ ResetTestControls reset; iga::SetSurfaceQuadratureTestRejectPositiveArea(true); const iga::ImmersedSurfaceQuadratureCatalog failed(fine_domain); assert(!failed.Usable() && failed.Diagnostics().candidate_attempts>failed.Diagnostics().candidate_pairs); assert(failed.Diagnostics().positive_fragments==0 && failed.Diagnostics().output_points==0); for(const auto& cell:failed.Cells()) { assert(!cell.usable && cell.point_count==0); Reject([&]{ failed.UsableRule(fine_domain,cell.id); }); } }
	{ ResetTestControls reset; iga::SetSurfaceQuadratureTestRejectPositiveArea(true); Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog candidate_capped(fine_domain,{1,1000,12000,512}); }); Reject([&]{ iga::ImmersedSurfaceQuadratureCatalog fragment_capped(fine_domain,{1000,1,12000,512}); }); }
	{ ResetTestControls reset; iga::SetSurfaceQuadratureTestRejectRuleValidationAfter(1); const iga::ImmersedSurfaceQuadratureCatalog failed_late(fine_domain,{1000000,1000000,12000000,512}); assert(!failed_late.Usable() && !failed_late.Diagnostics().catalog_unusable_reason.empty()); assert(failed_late.Cells().size()==fine_domain.Cells().size() && failed_late.Diagnostics().candidate_attempts>0 && failed_late.Diagnostics().positive_fragments>0 && failed_late.Diagnostics().output_points>0); bool pre_failure_completed=false, rolled_back_failure=false, post_failure_completed=false; for(const auto& cell:failed_late.Cells()) { if (!cell.usable) rolled_back_failure=rolled_back_failure || (cell.point_count==0 && cell.area_m2==0.0); else if (cell.point_count>0) { if (rolled_back_failure) post_failure_completed=true; else pre_failure_completed=true; } Reject([&]{ failed_late.UsableRule(fine_domain,cell.id); }); } assert(pre_failure_completed && rolled_back_failure && post_failure_completed); assert(iga::SurfaceQuadratureTestControls().rule_validation_attempts==2); }
	{ ResetTestControls reset; iga::exact_dyadic::SetTestMagnitudeCap(1); const iga::ImmersedSurfaceQuadratureCatalog owned(cube_domain,{1000000,1000000,12000000,128}); assert(owned.Usable() && owned.Diagnostics().peak_exact_limbs<=128); }
#endif
	std::cout << "immersed surface quadrature tests passed\n";
}

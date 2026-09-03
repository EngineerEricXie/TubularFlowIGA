#ifndef IGA_CUT_CELL_GHOST_PENALTY_HPP
#define IGA_CUT_CELL_GHOST_PENALTY_HPP

// A bounded per-face cubic ghost penalty.  The catalog borrows immutable
// caller-owned domain/volume catalogs; callers must keep both alive. Binding
// validation compares stored object identities only (O(1), no stale dereference).
#include "CutCellVolumeQuadrature.hpp"
#include "NavierStokesElement.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
namespace iga {
struct CutCellGhostPenaltyOptions { double gamma_u=.01, gamma_p=.01; std::size_t max_faces=1000000, max_quadrature_points=16000000, max_trace_entries=128; };
struct CutCellGhostPenaltyFace { std::uint64_t minus_cell=0, plus_cell=0; std::uint8_t axis=0; double h_normal_m=0., area_m2=0.; };
struct CutCellGhostPenaltyDiagnostics {
	std::size_t active_cells=0, full_cells=0, cut_cells=0, active_nodes=0, candidate_faces=0, selected_faces=0;
	std::array<std::size_t,3> selected_by_axis{{0,0,0}}; std::size_t cut_inside_faces=0, cut_cut_faces=0, uncovered_cells=0, uncovered_components=0, quadrature_points=0;
	double minimum_h_normal_m=0., maximum_h_normal_m=0., minimum_face_area_m2=0., maximum_face_area_m2=0.;
	// Geometry-normalized, mu=1 values; never a placeholder zero for selected faces.
	double maximum_abs_trace_jump=0., maximum_abs_coefficient=0.;
};
struct CutCellGhostPenaltyTrace { std::vector<std::int32_t> connectivity; std::vector<double> coefficients; };
struct CutCellGhostPenaltyAssembly { std::vector<std::int32_t> connectivity; std::vector<PetscScalar> jacobian, negative_residual; double maximum_abs_contribution=0.; };

class CutCellGhostPenaltyCatalog {
public:
	static constexpr std::uint32_t kCatalogVersion=2, kSplineDegree=3, kContinuity=2, kDerivativeOrder=3, kTangentQuadratureOrder=4;
	CutCellGhostPenaltyCatalog(const CutCellGhostPenaltyCatalog&)=delete;
	CutCellGhostPenaltyCatalog& operator=(const CutCellGhostPenaltyCatalog&)=delete;
	CutCellGhostPenaltyCatalog(const CartesianDomainClassification& domain, const CutCellVolumeQuadratureCatalog& volume, CutCellGhostPenaltyOptions options={})
		: bound_domain_(&domain), bound_volume_(&volume), grid_spec_(domain.Background().Spec()), surface_hash_(domain.SurfaceCanonicalHash()), cell_count_(domain.Cells().size()), options_(ValidateOptions(options))
	{
		if (!SameGrid(grid_spec_,volume.GridSpec()) || surface_hash_!=volume.SurfaceCanonicalHash()) throw std::invalid_argument("ghost penalty volume catalog does not match domain binding");
		const auto& cells=domain.Cells(); const auto& vcells=volume.Cells();
		if(cells.size()!=vcells.size() || cells.size()!=domain.Background().ElementCount()) throw std::invalid_argument("ghost penalty requires complete ordered domain and volume catalogs");
		std::vector<bool> active(cells.size()), positive_cut(cells.size()), sensitive(cells.size()), full(cells.size());
		for(std::size_t i=0;i<cells.size();++i) {
			if(cells[i].id!=i || vcells[i].id!=i || cells[i].classification!=vcells[i].classification) throw std::invalid_argument("ghost penalty catalog cell ordering or classification differs");
			const double f=vcells[i].diagnostics.estimated_reference_volume;
			if(!std::isfinite(f)||f<0.||f>1.) throw std::invalid_argument("ghost penalty cut fraction is not finite in [0,1]");
			if(cells[i].classification==CellClassification::Inside) {active[i]=full[i]=true;++diagnostics_.active_cells;++diagnostics_.full_cells;}
			else if(cells[i].classification==CellClassification::Cut && vcells[i].usable && f>0.) {
				active[i]=positive_cut[i]=true;++diagnostics_.active_cells;++diagnostics_.cut_cells;
				const double tol=2e-12*std::max(1.,std::abs(f)); full[i]=vcells[i].diagnostics.lower_reference_volume>=1.-tol && vcells[i].diagnostics.upper_reference_volume<=1.+tol;
				if(full[i]) ++diagnostics_.full_cells; else sensitive[i]=true;
			}
		}
		for(std::uint64_t id=0;id<cells.size();++id) if(active[id]) {
			const auto e=domain.Background().MaterializeElement(id); active_nodes_.insert(active_nodes_.end(),e.connectivity.begin(),e.connectivity.end()); const auto cell=domain.Background().Cell(id);
			for(std::uint8_t axis=0;axis<3;++axis) {
				const auto slot=axis==0?CubicCartesianBackground::XPlus:axis==1?CubicCartesianBackground::YPlus:CubicCartesianBackground::ZPlus; const auto plus=cell.neighbor[slot];
				if(plus==CubicCartesianBackground::kNoNeighbor || !active[plus]) continue;
				++diagnostics_.candidate_faces;
				// A certified-full Cut remains a Cut face selector; it only becomes an anchor.
				if(!SelectFace(positive_cut[id], positive_cut[plus])) continue;
				if(faces_.size()>=options_.max_faces) throw std::runtime_error("ghost penalty face cap exceeded");
				const auto nq=CheckedAdd(diagnostics_.quadrature_points,16,"ghost penalty quadrature count overflows"); if(nq>options_.max_quadrature_points) throw std::runtime_error("ghost penalty quadrature cap exceeded");
				CutCellGhostPenaltyFace face{id,plus,axis,cell.upper_m[axis]-cell.lower_m[axis],1.}; for(int t=0;t<3;++t) if(t!=axis) face.area_m2*=cell.upper_m[t]-cell.lower_m[t];
				if(!std::isfinite(face.h_normal_m)||!(face.h_normal_m>0.)||!std::isfinite(face.area_m2)||!(face.area_m2>0.)) throw std::runtime_error("ghost penalty face has invalid Cartesian metric");
				const auto nodes=FaceConnectivity(domain,face); if(nodes.size()>options_.max_trace_entries) throw std::runtime_error("ghost penalty trace-entry cap exceeded");
				faces_.push_back(face); diagnostics_.quadrature_points=nq; ++diagnostics_.selected_by_axis[axis]; if(positive_cut[id]&&positive_cut[plus]) ++diagnostics_.cut_cut_faces; else ++diagnostics_.cut_inside_faces;
				Range(diagnostics_.minimum_h_normal_m,diagnostics_.maximum_h_normal_m,face.h_normal_m,faces_.size()==1); Range(diagnostics_.minimum_face_area_m2,diagnostics_.maximum_face_area_m2,face.area_m2,faces_.size()==1); Precompute(domain,faces_.size()-1);
			}
		}
		diagnostics_.selected_faces=faces_.size();
		if (!faces_.empty() && (!std::isfinite(diagnostics_.maximum_abs_coefficient) || !(diagnostics_.maximum_abs_coefficient > 0.)))
			throw std::overflow_error("ghost penalty selected-face coefficient diagnostic is not finite and positive");
		std::sort(active_nodes_.begin(),active_nodes_.end()); active_nodes_.erase(std::unique(active_nodes_.begin(),active_nodes_.end()),active_nodes_.end()); diagnostics_.active_nodes=active_nodes_.size(); AuditCoverage(active,sensitive,full); usable_=true;
	}
	const std::vector<CutCellGhostPenaltyFace>& Faces() const noexcept{return faces_;} const CutCellGhostPenaltyDiagnostics& Diagnostics() const noexcept{return diagnostics_;} const CutCellGhostPenaltyOptions& Options() const noexcept{return options_;} const std::vector<std::int32_t>& ActiveNodes() const noexcept{return active_nodes_;} const std::vector<std::uint64_t>& UncoveredCells() const noexcept{return uncovered_cells_;} bool Usable() const noexcept{return usable_;} std::uint32_t Version() const noexcept{return kCatalogVersion;}
	// Kept separately testable because a classified Cut can be certified full:
	// full status establishes coverage, never suppresses its selected face.
	static bool SelectFace(bool positive_cut_minus, bool positive_cut_plus) noexcept
	{
		return positive_cut_minus || positive_cut_plus;
	}
	bool Covered(std::uint64_t id) const {if(id>=covered_.size())throw std::out_of_range("ghost penalty coverage cell id is out of range");return covered_[id];}
	bool Matches(const CartesianDomainClassification& d,const CutCellVolumeQuadratureCatalog& v) const noexcept{return usable_&&&d==bound_domain_&&&v==bound_volume_;}
	void ValidateBinding(const CartesianDomainClassification& d,const CutCellVolumeQuadratureCatalog& v) const {if(!Matches(d,v))throw std::invalid_argument("ghost penalty catalog binding does not match immutable caller-owned domain/volume catalogs");}
	CutCellGhostPenaltyTrace EvaluateFaceJump(std::size_t fi,std::uint32_t order,double t0,double t1,const CartesianDomainClassification& d,const CutCellVolumeQuadratureCatalog& v) const {
		ValidateBinding(d,v); if(fi>=faces_.size())throw std::out_of_range("ghost penalty face index is out of range"); if(order>3||!std::isfinite(t0)||!std::isfinite(t1)||t0<0.||t0>1.||t1<0.||t1>1.)throw std::invalid_argument("ghost penalty trace evaluation arguments are invalid");
		CutCellGhostPenaltyTrace r; r.connectivity=FaceConnectivity(d,faces_[fi]); if(r.connectivity.size()>options_.max_trace_entries)throw std::runtime_error("ghost penalty trace-entry cap exceeded"); const auto& f=faces_[fi]; r.coefficients=TraceJump(d.Background().MaterializeElement(f.minus_cell),d.Background().MaterializeElement(f.plus_cell),f.axis,order,t0,t1,f.h_normal_m,r.connectivity); return r;
	}
	template <class StateAt> CutCellGhostPenaltyAssembly AssembleFaceLocal(std::size_t fi,const CartesianDomainClassification& d,const CutCellVolumeQuadratureCatalog& v,StateAt&& state_at,double mu) const {
		ValidateBinding(d,v); if(fi>=faces_.size())throw std::out_of_range("ghost penalty face index is out of range"); if(!std::isfinite(mu)||!(mu>0.))throw std::invalid_argument("ghost penalty viscosity must be finite and positive");
		CutCellGhostPenaltyAssembly r; r.connectivity=FaceConnectivity(d,faces_[fi]); const auto n=r.connectivity.size(); if(n>80||n>options_.max_trace_entries)throw std::runtime_error("ghost penalty face trace union exceeds configured bound"); const auto ndof=Mul(4,n,"ghost penalty face dof count overflows"), entries=Mul(ndof,ndof,"ghost penalty face dense allocation overflows"); r.jacobian.assign(entries,PetscScalar(0));r.negative_residual.assign(ndof,PetscScalar(0)); const auto& f=faces_[fi]; const double ku=options_.gamma_u*mu*std::pow(f.h_normal_m,5),kp=options_.gamma_p/mu*std::pow(f.h_normal_m,7);if(!std::isfinite(ku)||!(ku>0.)||!std::isfinite(kp)||!(kp>0.))throw std::overflow_error("ghost penalty coefficient is not finite and positive");
		static constexpr std::array<double,4> x{{-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526}},w{{.3478548451374538,.6521451548625461,.6521451548625461,.3478548451374538}};
		for(std::size_t q0=0;q0<4;++q0)for(std::size_t q1=0;q1<4;++q1){const double wt=w[q0]*w[q1]*f.area_m2/4.;if(!std::isfinite(wt)||!(wt>0.))throw std::overflow_error("ghost penalty face weight is invalid");const auto j=EvaluateFaceJump(fi,3,(x[q0]+1.)/2.,(x[q1]+1.)/2.,d,v).coefficients;for(std::size_t a=0;a<n;++a)for(std::size_t b=0;b<n;++b){const double uv=ku*wt*j[a]*j[b],pp=kp*wt*j[a]*j[b];if(!std::isfinite(uv)||!std::isfinite(pp))throw std::overflow_error("ghost penalty contribution is not finite");r.maximum_abs_contribution=std::max(r.maximum_abs_contribution,std::max(std::abs(uv),std::abs(pp)));for(int c=0;c<3;++c)Add(r.jacobian[(4*a+c)*ndof+4*b+c],uv);Add(r.jacobian[(4*a+3)*ndof+4*b+3],pp);}}
		for(std::size_t a=0;a<ndof;++a)
			for(std::size_t b=0;b<ndof;++b)
				Add(r.negative_residual[a],-PetscRealPart(r.jacobian[a*ndof+b])*state_at(r.connectivity[b/4],static_cast<int>(b%4)));
		return r;
	}
	CutCellGhostPenaltyAssembly AssembleFace(std::size_t fi,const CartesianDomainClassification& d,const CutCellVolumeQuadratureCatalog& v,const std::vector<std::array<double,4>>& state,double mu) const {
		if(state.size()!=d.Background().NodeCount()) throw std::invalid_argument("ghost penalty state must contain four fields for every background node");
		for(const auto& q:state) for(double x:q) if(!std::isfinite(x)) throw std::invalid_argument("ghost penalty state is not finite");
		return AssembleFaceLocal(fi,d,v,[&state](std::int32_t node,int field) { return state[static_cast<std::size_t>(node)][field]; },mu);
	}
private:
	static bool SameGrid(const CubicCartesianGridSpec&a,const CubicCartesianGridSpec&b)noexcept{return a.lower_m==b.lower_m&&a.upper_m==b.upper_m&&a.cells==b.cells;}
	static std::size_t CheckedAdd(std::size_t a,std::size_t b,const char*msg){if(a>std::numeric_limits<std::size_t>::max()-b)throw std::overflow_error(msg);return a+b;} static std::size_t Mul(std::size_t a,std::size_t b,const char*msg){if(a&&b>std::numeric_limits<std::size_t>::max()/a)throw std::overflow_error(msg);return a*b;}
	static CutCellGhostPenaltyOptions ValidateOptions(CutCellGhostPenaltyOptions o){if(!std::isfinite(o.gamma_u)||!(o.gamma_u>0.)||!std::isfinite(o.gamma_p)||!(o.gamma_p>0.)||!o.max_faces||!o.max_quadrature_points||!o.max_trace_entries)throw std::invalid_argument("ghost penalty options must be finite positive");return o;}
	static void Add(PetscScalar& a,double b){if(!std::isfinite(b)||!std::isfinite(PetscRealPart(a)))throw std::overflow_error("ghost penalty accumulation is not finite");const auto n=a+b;if(!std::isfinite(PetscRealPart(n)))throw std::overflow_error("ghost penalty accumulation overflows");a=n;} static void Range(double&lo,double&hi,double x,bool first){if(first)lo=hi=x;else{lo=std::min(lo,x);hi=std::max(hi,x);}}
	static std::vector<std::int32_t> FaceConnectivity(const CartesianDomainClassification&d,const CutCellGhostPenaltyFace&f){auto r=d.Background().MaterializeElement(f.minus_cell).connectivity;const auto p=d.Background().MaterializeElement(f.plus_cell);r.insert(r.end(),p.connectivity.begin(),p.connectivity.end());std::sort(r.begin(),r.end());r.erase(std::unique(r.begin(),r.end()),r.end());if(r.size()>80)throw std::runtime_error("ghost penalty adjacent cubic face union exceeds 80 nodes");return r;}
	static std::array<double,4> B(std::uint32_t o,double t){const double u=1.-t;if(!o)return{{u*u*u,3*u*u*t,3*u*t*t,t*t*t}};if(o==1)return{{-3*u*u,3*u*u-6*u*t,6*u*t-3*t*t,3*t*t}};if(o==2)return{{6*u,-12+18*t,6-18*t,6*t}};return{{-6,18,-18,6}};}
	static std::vector<double> TraceJump(const Element&minus,const Element&plus,std::uint8_t axis,std::uint32_t order,double t0,double t1,double h,const std::vector<std::int32_t>&nodes){const auto trace=[&](const Element&e,bool lower){std::vector<double>r(nodes.size());const auto n=B(order,lower?1.:0.),b0=B(0,t0),b1=B(0,t1);const double scale=std::pow(h,-static_cast<int>(order));for(std::size_t row=0;row<e.connectivity.size();++row){double z=0.;for(int c=0;c<4;++c)for(int b=0;b<4;++b)for(int a=0;a<4;++a){const int idx=a+4*(b+4*c),ni=axis==0?a:axis==1?b:c,ti0=axis==0?b:a,ti1=axis==2?b:c;z+=e.extraction[row][idx]*n[ni]*b0[ti0]*b1[ti1];}const auto pos=static_cast<std::size_t>(std::lower_bound(nodes.begin(),nodes.end(),e.connectivity[row])-nodes.begin());r[pos]+=(lower?1.:-1.)*scale*z;}return r;};auto r=trace(minus,true);const auto p=trace(plus,false);for(std::size_t i=0;i<r.size();++i)r[i]+=p[i];return r;}
	static double NormalizedCoefficient(double gamma,double h,int exponent){const double h_power=std::pow(h,exponent);if(!std::isfinite(h_power)||!(h_power>0.))throw std::overflow_error("ghost penalty metric power is not finite and positive");const double coefficient=gamma*h_power;if(!std::isfinite(coefficient)||!(coefficient>0.))throw std::overflow_error("ghost penalty normalized coefficient is not finite and positive");return coefficient;}
	void Precompute(const CartesianDomainClassification&d,std::size_t fi){const auto&f=faces_[fi];const double ku=NormalizedCoefficient(options_.gamma_u,f.h_normal_m,5),kp=NormalizedCoefficient(options_.gamma_p,f.h_normal_m,7);const auto nodes=FaceConnectivity(d,f);static constexpr std::array<double,4>x{{-.8611363115940526,-.3399810435848563,.3399810435848563,.8611363115940526}};for(double a:x)for(double b:x)for(double z:TraceJump(d.Background().MaterializeElement(f.minus_cell),d.Background().MaterializeElement(f.plus_cell),f.axis,3,(a+1.)/2.,(b+1.)/2.,f.h_normal_m,nodes)){if(!std::isfinite(z))throw std::overflow_error("ghost penalty trace diagnostic is not finite");diagnostics_.maximum_abs_trace_jump=std::max(diagnostics_.maximum_abs_trace_jump,std::abs(z));}diagnostics_.maximum_abs_coefficient=std::max(diagnostics_.maximum_abs_coefficient,std::max(ku,kp));}
	void AuditCoverage(const std::vector<bool>&active,const std::vector<bool>&sensitive,const std::vector<bool>&full){covered_.assign(active.size(),false);std::vector<std::vector<std::uint64_t>>a(active.size());for(const auto&f:faces_){a[f.minus_cell].push_back(f.plus_cell);a[f.plus_cell].push_back(f.minus_cell);}std::queue<std::uint64_t>q;for(std::uint64_t i=0;i<active.size();++i)if(active[i]&&full[i]){covered_[i]=true;q.push(i);}while(!q.empty()){auto i=q.front();q.pop();for(auto j:a[i])if(!covered_[j]){covered_[j]=true;q.push(j);}}std::vector<bool>mark(active.size());for(std::uint64_t i=0;i<active.size();++i)if(sensitive[i]&&!covered_[i]){uncovered_cells_.push_back(i);if(mark[i])continue;++diagnostics_.uncovered_components;mark[i]=true;q.push(i);while(!q.empty()){auto k=q.front();q.pop();for(auto j:a[k])if(sensitive[j]&&!covered_[j]&&!mark[j]){mark[j]=true;q.push(j);}}}diagnostics_.uncovered_cells=uncovered_cells_.size();}
	const CartesianDomainClassification* bound_domain_=nullptr;const CutCellVolumeQuadratureCatalog* bound_volume_=nullptr;CubicCartesianGridSpec grid_spec_{};std::string surface_hash_;std::size_t cell_count_=0;CutCellGhostPenaltyOptions options_{};CutCellGhostPenaltyDiagnostics diagnostics_{};std::vector<CutCellGhostPenaltyFace>faces_;std::vector<std::int32_t>active_nodes_;std::vector<std::uint64_t>uncovered_cells_;std::vector<bool>covered_;bool usable_=false;
};
} // namespace iga
#endif

#ifndef IGA_IMMERSED_SURFACE_QUADRATURE_HPP
#define IGA_IMMERSED_SURFACE_QUADRATURE_HPP

// Surface rules for Cartesian cut cells.  Clipping uses homogeneous exact
// dyadics reconstructed from source and box constraints; doubles are only
// used after all topological choices (including duplicate removal).
#include "CartesianDomainClassification.hpp"
#include "Quadrature.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

struct ImmersedSurfaceQuadratureOptions {
	std::size_t max_candidates = 1000000;
	std::size_t max_fragments = 1000000;
	std::size_t max_points = 12000000;
	std::size_t max_exact_limbs = 512;
};

struct ImmersedSurfaceQuadratureDiagnostics {
	// Attempt counters are deliberately monotone.  They include work later
	// discarded with an unusable cell, so a retry cannot evade a resource cap.
	std::size_t candidate_attempts = 0;
	std::size_t fragment_attempts = 0;
	std::size_t point_attempts = 0;
	std::size_t candidate_pairs = 0;
	std::size_t clipped_polygons = 0;
	std::size_t positive_fragments = 0;
	std::size_t zero_fragments = 0;
	std::size_t precision_limited_fragments = 0;
	std::size_t output_points = 0;
	std::size_t coplanar_owner_skips = 0;
	std::size_t predicate_ambiguities = 0;
	std::size_t peak_exact_limbs = 0;
	std::size_t unusable_cells = 0;
	double total_area_m2 = 0.0;
	double source_total_area_m2 = 0.0;
	double accumulated_total_area_m2 = 0.0;
	double total_area_residual_m2 = 0.0;
	double total_area_absolute_residual_m2 = 0.0;
	std::map<std::uint32_t, double> area_by_boundary_id;
	std::map<std::uint32_t, double> source_area_by_boundary_id;
	std::map<std::uint32_t, double> area_residual_by_boundary_id;
	std::map<std::uint32_t, double> area_absolute_residual_by_boundary_id;
	std::vector<double> triangle_area_m2;
	std::vector<double> triangle_area_residual_m2;
	std::vector<double> triangle_area_absolute_residual_m2;
	double max_mapping_residual_m = 0.0;
	double max_surface_residual_m = 0.0;
	bool catalog_usable = true;
	std::string catalog_unusable_reason;
};

#ifdef IGA_EXACT_DYADIC_TESTING
// Test-only fault injection is intentionally local to this catalog.  It
// exercises the post-exact, pre-publication paths without weakening the
// production exact predicate capacity.
struct ImmersedSurfaceQuadratureTestControls {
	bool reject_positive_area = false;
	bool reject_rule_validation = false;
	std::size_t reject_rule_validation_after = 0;
	std::size_t rule_validation_attempts = 0;
};
inline ImmersedSurfaceQuadratureTestControls& SurfaceQuadratureTestControls()
{ static ImmersedSurfaceQuadratureTestControls controls; return controls; }
inline void SetSurfaceQuadratureTestRejectPositiveArea(bool value)
{ SurfaceQuadratureTestControls().reject_positive_area = value; }
inline void SetSurfaceQuadratureTestRejectRuleValidation(bool value)
{
	SurfaceQuadratureTestControls().reject_rule_validation = value;
	SurfaceQuadratureTestControls().reject_rule_validation_after = 0;
	SurfaceQuadratureTestControls().rule_validation_attempts = 0;
}
inline void SetSurfaceQuadratureTestRejectRuleValidationAfter(std::size_t successful_attempts)
{
	SurfaceQuadratureTestControls().reject_rule_validation = true;
	SurfaceQuadratureTestControls().reject_rule_validation_after = successful_attempts;
	SurfaceQuadratureTestControls().rule_validation_attempts = 0;
}
#endif

struct ImmersedSurfaceQuadratureCell {
	std::uint64_t id = 0;
	bool usable = true;
	bool ambiguous = false;
	double area_m2 = 0.0;
	std::size_t point_count = 0;
};

// Kept parallel to SurfaceQuadratureRule so body-fitted users retain the
// generic point type and its aggregate construction.  The triangle index is
// the canonical index in the bound ClosedTriangulatedSurface.
struct ImmersedSurfaceQuadraturePointProvenance {
	std::uint32_t canonical_triangle = 0;
	std::array<double, 3> canonical_barycentric{{0.0, 0.0, 0.0}};
};

class ImmersedSurfaceQuadratureCatalog {
public:
	ImmersedSurfaceQuadratureCatalog(const ImmersedSurfaceQuadratureCatalog&) = delete;
	ImmersedSurfaceQuadratureCatalog& operator=(const ImmersedSurfaceQuadratureCatalog&) = delete;

	explicit ImmersedSurfaceQuadratureCatalog(const CartesianDomainClassification& domain,
		ImmersedSurfaceQuadratureOptions options = {})
		: grid_spec_(domain.Background().Spec()), surface_canonical_hash_(domain.SurfaceCanonicalHash()),
		options_(ValidateOptions(options))
	{
		const auto& source_cells = domain.Cells();
		const auto& surface = domain.SurfaceIndex().Surface();
		if (source_cells.size() != domain.Background().ElementCount())
			throw std::runtime_error("Cartesian domain catalog is incomplete");
		ValidateBoundaryIds(surface);
		diagnostics_.triangle_area_m2.assign(surface.Triangles().size(), 0.0);
		diagnostics_.triangle_area_residual_m2.assign(surface.Triangles().size(), 0.0);
		diagnostics_.triangle_area_absolute_residual_m2.assign(surface.Triangles().size(), 0.0);
		triangle_area_accumulated_.assign(surface.Triangles().size(), {});
		triangle_fragment_counts_.assign(surface.Triangles().size(), 0);
		provenance_.resize(source_cells.size());
		cells_.reserve(source_cells.size());
		for (const auto& source : source_cells) {
			if (source.id != cells_.size()) throw std::runtime_error("Cartesian domain cell ids are not x-fast");
			ImmersedSurfaceQuadratureCell result; SurfaceQuadratureRule rule; result.id = source.id;
			BuildCell(domain, source, result, rule);
			if (!result.usable) { ++diagnostics_.unusable_cells; MarkCatalogUnusable("a cell could not certify its surface measure"); }
			cells_.push_back(std::move(result));
			rules_.push_back(std::move(rule));
		}
		FinalizeAudit(surface);
	}

	ImmersedSurfaceQuadratureCatalog(ImmersedSurfaceQuadratureCatalog&&) noexcept = default;
	ImmersedSurfaceQuadratureCatalog& operator=(ImmersedSurfaceQuadratureCatalog&&) noexcept = default;
	const CubicCartesianGridSpec& GridSpec() const noexcept { return grid_spec_; }
	const std::string& SurfaceCanonicalHash() const noexcept { return surface_canonical_hash_; }
	const ImmersedSurfaceQuadratureOptions& Options() const noexcept { return options_; }
	const ImmersedSurfaceQuadratureDiagnostics& Diagnostics() const noexcept { return diagnostics_; }
	bool Usable() const noexcept { return diagnostics_.catalog_usable; }
	const std::vector<ImmersedSurfaceQuadratureCell>& Cells() const noexcept { return cells_; }
	const ImmersedSurfaceQuadratureCell& Cell(std::uint64_t id) const
	{
		if (id >= cells_.size()) throw std::out_of_range("immersed surface quadrature id is out of range");
		return cells_[static_cast<std::size_t>(id)];
	}
	const SurfaceQuadratureRule& UsableRule(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		ValidateBinding(domain, id); const auto& cell = Cell(id);
		if (!diagnostics_.catalog_usable || !cell.usable) throw std::runtime_error("immersed surface quadrature catalog is unusable");
		return rules_[static_cast<std::size_t>(id)];
	}
	const std::vector<ImmersedSurfaceQuadraturePointProvenance>& UsableProvenance(
		const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		const auto& rule = UsableRule(domain, id);
		const auto& provenance = provenance_[static_cast<std::size_t>(id)];
		if (provenance.size() != rule.Points().size())
			throw std::runtime_error("immersed surface quadrature provenance point count does not match rule");
		return provenance;
	}
	void ValidateUsableRule(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		const auto& rule = UsableRule(domain, id);
		ValidateSurfaceQuadratureRule(domain.Background().MaterializeElement(id), rule);
	}

private:
	using Number = exact_dyadic::Number;
	struct CompensatedArea {
		long double sum = 0.0L;
		long double correction = 0.0L;
		void Add(long double value)
		{
			const long double next = sum+value;
			if (std::abs(sum) >= std::abs(value)) correction += (sum-next)+value;
			else correction += (value-next)+sum;
			sum = next;
		}
		long double Value() const { return sum+correction; }
	};
	struct HomogeneousPoint {
		std::array<Number, 4> value{};
		std::array<Number, 4> material{};
		std::uint16_t active = 0;
	};
	struct BoxPlane { std::size_t axis = 0; double value = 0.0; };

	static ImmersedSurfaceQuadratureOptions ValidateOptions(ImmersedSurfaceQuadratureOptions options)
	{
		if (!options.max_candidates || !options.max_fragments || !options.max_points
			|| !options.max_exact_limbs || options.max_exact_limbs > 512)
			throw std::invalid_argument("immersed surface quadrature caps must be positive");
		return options;
	}
	Number Add(const Number& a, const Number& b) { return exact_dyadic::Add(a,b,options_.max_exact_limbs); }
	Number Subtract(const Number& a, const Number& b) { return exact_dyadic::Subtract(a,b,options_.max_exact_limbs); }
	Number Multiply(const Number& a, const Number& b) { return exact_dyadic::Multiply(a,b,options_.max_exact_limbs); }
	void Observe(const Number& value)
	{ diagnostics_.peak_exact_limbs = std::max(diagnostics_.peak_exact_limbs, value.magnitude.limbs.size()); }
	Number A(const Number& a, const Number& b) { auto value=Add(a,b); Observe(value); return value; }
	Number S(const Number& a, const Number& b) { auto value=Subtract(a,b); Observe(value); return value; }
	Number M(const Number& a, const Number& b) { auto value=Multiply(a,b); Observe(value); return value; }
	Number E(double value) { auto exact=D(value); Observe(exact); return exact; }
	static Number D(double value) { return exact_dyadic::FromDouble(value); }
	static int Sign(const Number& value) { return exact_dyadic::Sign(value); }
	static Number Negate(Number value) { return exact_dyadic::Negate(std::move(value)); }
	static std::size_t EdgeStart(std::size_t edge) { return edge; }
	static std::size_t EdgeEnd(std::size_t edge) { return (edge+1)%3; }
	Number PlaneValue(const HomogeneousPoint& point, std::size_t axis, double plane)
	{ return S(point.value[axis], M(D(plane),point.value[3])); }
	bool Equal(const HomogeneousPoint& a, const HomogeneousPoint& b)
	{
		for (std::size_t axis=0; axis<3; ++axis) if (Sign(S(M(a.value[axis],b.value[3]),M(b.value[axis],a.value[3])))) return false;
		return true;
	}
	static double ToDouble(const Number& numerator, const Number& denominator)
	{
		if (!numerator.sign) return 0.0;
		// Scale each magnitude by its leading limb before division.  The full
		// numerator and denominator can exceed long-double range even when their
		// ratio is an ordinary physical coordinate.
		auto normalized = [](const Number& number, long double& mantissa, long long& exponent) {
			if (!number.sign || number.magnitude.limbs.empty()) throw std::logic_error("invalid rational dyadic");
			const std::uint64_t leading = number.magnitude.limbs.back(); std::size_t bits = 0;
			for (std::uint64_t value = leading; value; value >>= 1) ++bits;
			mantissa = static_cast<long double>(leading);
			for (std::size_t i = number.magnitude.limbs.size()-1; i != 0; --i)
				mantissa += std::ldexp(static_cast<long double>(number.magnitude.limbs[i-1]),-static_cast<int>(64*(number.magnitude.limbs.size()-i)));
			mantissa = std::ldexp(mantissa, -static_cast<int>(bits));
			exponent = static_cast<long long>(number.exponent)+static_cast<long long>((number.magnitude.limbs.size()-1)*64+bits);
		};
		long double numerator_mantissa = 0.0L, denominator_mantissa = 0.0L; long long numerator_exponent = 0, denominator_exponent = 0;
		normalized(numerator,numerator_mantissa,numerator_exponent); normalized(denominator,denominator_mantissa,denominator_exponent);
		const long long exponent = numerator_exponent-denominator_exponent;
		if (exponent > std::numeric_limits<int>::max() || exponent < std::numeric_limits<int>::min())
			throw std::runtime_error("rational clipped coordinate is not representable as double");
		const long double result = std::ldexp(numerator_mantissa/denominator_mantissa,static_cast<int>(exponent))
			*(numerator.sign == denominator.sign ? 1.0L : -1.0L);
		if (!std::isfinite(result) || !std::isfinite(static_cast<double>(result)))
			throw std::runtime_error("rational clipped coordinate is not representable as double");
		return static_cast<double>(result);
	}
	HomogeneousPoint Point(const std::array<double,3>& point, const std::array<std::array<double,3>,3>& original,
		const std::array<BoxPlane,6>& planes, std::size_t corner)
	{
		HomogeneousPoint result; for(std::size_t a=0;a<3;++a) { result.value[a]=D(point[a]); result.material[a]=D(a==corner ? 1.0 : 0.0); } result.value[3]=D(1.0); result.material[3]=D(1.0); RecomputeActive(result,original,planes); return result;
	}
	std::array<Number,3> TriangleNormal(const std::array<std::array<double,3>,3>& original)
	{
		std::array<Number,3> u{},v{},n{}; for(std::size_t a=0;a<3;++a) { u[a]=S(D(original[1][a]),D(original[0][a])); v[a]=S(D(original[2][a]),D(original[0][a])); }
		n[0]=S(M(u[1],v[2]),M(u[2],v[1])); n[1]=S(M(u[2],v[0]),M(u[0],v[2])); n[2]=S(M(u[0],v[1]),M(u[1],v[0])); return n;
	}
	void RecomputeActive(HomogeneousPoint& point, const std::array<std::array<double,3>,3>& original,
		const std::array<BoxPlane,6>& planes)
	{
		if(Sign(point.value[3])<=0) throw std::runtime_error("surface clipping homogeneous weight is not positive");
		point.active=0;
		for(std::size_t edge=0;edge<3;++edge) {
			const auto begin=EdgeStart(edge), end=EdgeEnd(edge); std::array<Number,3> direction{}, offset{}, cross{};
			for(std::size_t a=0;a<3;++a) { direction[a]=S(D(original[end][a]),D(original[begin][a])); offset[a]=S(point.value[a],M(D(original[begin][a]),point.value[3])); }
			cross[0]=S(M(direction[1],offset[2]),M(direction[2],offset[1])); cross[1]=S(M(direction[2],offset[0]),M(direction[0],offset[2])); cross[2]=S(M(direction[0],offset[1]),M(direction[1],offset[0]));
			if(!Sign(cross[0])&&!Sign(cross[1])&&!Sign(cross[2])) point.active|=static_cast<std::uint16_t>(1u<<edge);
		}
		for(std::size_t plane=0;plane<planes.size();++plane) if(!Sign(PlaneValue(point,planes[plane].axis,planes[plane].value))) point.active|=static_cast<std::uint16_t>(1u<<(3+plane));
	}
	HomogeneousPoint IntersectPlane(const HomogeneousPoint& left, const HomogeneousPoint& right, std::size_t new_plane,
		const std::array<std::array<double,3>,3>& original, const std::array<BoxPlane,6>& planes)
	{
		const auto axis=planes[new_plane].axis; const Number p=D(planes[new_plane].value); const auto common=static_cast<std::uint16_t>(left.active&right.active);
		const Number left_plane=PlaneValue(left,axis,planes[new_plane].value), right_plane=PlaneValue(right,axis,planes[new_plane].value);
		auto material_intersection = [&]() { const Number left_scale=M(left_plane,right.value[3]), right_scale=M(right_plane,left.value[3]); std::array<Number,4> material{}; for (std::size_t i=0;i<3;++i) material[i]=S(M(M(left_scale,right.material[i]),left.material[3]),M(M(right_scale,left.material[i]),right.material[3])); material[3]=S(M(M(left_scale,right.material[3]),left.material[3]),M(M(right_scale,left.material[3]),right.material[3])); if (Sign(material[3]) < 0) for (auto& value : material) value=Negate(std::move(value)); if (Sign(material[3]) <= 0) throw std::runtime_error("surface clipping material intersection has invalid weight"); return material; };
		for(std::size_t edge=0;edge<3;++edge) if(common&(1u<<edge)) {
			const auto begin=EdgeStart(edge),end=EdgeEnd(edge); const Number d=S(D(original[end][axis]),D(original[begin][axis])); const int sign=Sign(d); if(!sign) continue;
			HomogeneousPoint result; result.value[3]=sign>0?d:Negate(d); result.value[axis]=M(p,result.value[3]);
			for(std::size_t coordinate=0;coordinate<3;++coordinate) if(coordinate!=axis) { const Number term=A(M(D(original[begin][coordinate]),d),M(S(p,D(original[begin][axis])),S(D(original[end][coordinate]),D(original[begin][coordinate])))); result.value[coordinate]=sign>0?term:Negate(term); }
			result.material=material_intersection(); RecomputeActive(result,original,planes); return result;
		}
		const auto normal=TriangleNormal(original); Number k{}; for(std::size_t a=0;a<3;++a) k=A(k,M(normal[a],D(original[0][a])));
		for(std::size_t old=0;old<6;++old) if((common&(1u<<(3+old)))&&old!=new_plane) {
			const auto old_axis=planes[old].axis; if(old_axis==axis) continue; std::size_t remaining=0; while(remaining==axis||remaining==old_axis) ++remaining; const int sign=Sign(normal[remaining]); if(!sign) continue;
			HomogeneousPoint result; result.value[3]=sign>0?normal[remaining]:Negate(normal[remaining]); for(std::size_t coordinate=0;coordinate<3;++coordinate) result.value[coordinate]=D(0.0); result.value[axis]=M(p,result.value[3]); result.value[old_axis]=M(D(planes[old].value),result.value[3]); Number remainder=k; remainder=S(remainder,M(normal[axis],p)); remainder=S(remainder,M(normal[old_axis],D(planes[old].value))); result.value[remaining]=sign>0?remainder:Negate(remainder); result.material=material_intersection(); RecomputeActive(result,original,planes); return result;
		}
		throw std::runtime_error("surface clipping intersection lacks unique constraint provenance");
	}
	std::vector<HomogeneousPoint> ClipPlane(const std::vector<HomogeneousPoint>& input, std::size_t plane, bool keep_greater,
		const std::array<std::array<double,3>,3>& original, const std::array<BoxPlane,6>& planes)
	{
		std::vector<HomogeneousPoint> output; if(input.empty()) return output; for(std::size_t i=0;i<input.size();++i) { const auto& a=input[i]; const auto& b=input[(i+1)%input.size()]; const int sa=Sign(PlaneValue(a,planes[plane].axis,planes[plane].value)),sb=Sign(PlaneValue(b,planes[plane].axis,planes[plane].value)); const bool ina=keep_greater?sa>=0:sa<=0,inb=keep_greater?sb>=0:sb<=0; if(ina) output.push_back(a); if(ina!=inb&&sa!=0&&sb!=0) output.push_back(IntersectPlane(a,b,plane,original,planes)); }
		std::vector<HomogeneousPoint> deduped; for(const auto& point:output) if(deduped.empty()||!Equal(deduped.back(),point)) deduped.push_back(point); if(deduped.size()>1&&Equal(deduped.front(),deduped.back())) deduped.pop_back(); return deduped;
	}
	Number ProjectedDet(const HomogeneousPoint& a, const HomogeneousPoint& b, const HomogeneousPoint& c, std::size_t drop)
	{
		std::array<std::size_t,2> keep{}; for(std::size_t axis=0,n=0;axis<3;++axis) if(axis!=drop) keep[n++]=axis; return A(S(M(a.value[keep[0]],S(M(b.value[keep[1]],c.value[3]),M(b.value[3],c.value[keep[1]]))),M(a.value[keep[1]],S(M(b.value[keep[0]],c.value[3]),M(b.value[3],c.value[keep[0]])))),M(a.value[3],S(M(b.value[keep[0]],c.value[keep[1]]),M(b.value[keep[1]],c.value[keep[0]]))));
	}
	void RemoveCollinear(std::vector<HomogeneousPoint>& polygon, std::size_t projection)
	{ bool changed=true; while(changed&&polygon.size()>=3) { changed=false; for(std::size_t i=0;i<polygon.size();++i) if(!Sign(ProjectedDet(polygon[(i+polygon.size()-1)%polygon.size()],polygon[i],polygon[(i+1)%polygon.size()],projection))) { polygon.erase(polygon.begin()+static_cast<std::ptrdiff_t>(i)); changed=true; break; } } }
	enum class CoplanarDecision { Include, Skip, Invalid };
	CoplanarDecision CoplanarOwnership(const std::array<std::array<double,3>,3>& vertices,
		const CubicCartesianBackground& background, const CartesianDomainCell& cell,
		const ClosedSurfaceTriangle&)
	{
		for (std::size_t axis = 0; axis < 3; ++axis) for (int side = 0; side < 2; ++side) {
			const std::uint32_t plane_index = cell.index[axis]+static_cast<std::uint32_t>(side);
			const double plane = background.Plane(axis,plane_index);
			if (Sign(S(E(vertices[0][axis]),E(plane))) == 0
				&& Sign(S(E(vertices[1][axis]),E(plane))) == 0
				&& Sign(S(E(vertices[2][axis]),E(plane))) == 0) {
				std::array<Number,3> u{}, v{}, exact_cross{};
				for (std::size_t coordinate=0; coordinate<3; ++coordinate) {
					u[coordinate]=S(E(vertices[1][coordinate]),E(vertices[0][coordinate]));
					v[coordinate]=S(E(vertices[2][coordinate]),E(vertices[0][coordinate]));
				}
				exact_cross[0]=S(M(u[1],v[2]),M(u[2],v[1]));
				exact_cross[1]=S(M(u[2],v[0]),M(u[0],v[2]));
				exact_cross[2]=S(M(u[0],v[1]),M(u[1],v[0]));
				const int normal_sign = Sign(exact_cross[axis]);
				if (!normal_sign) return CoplanarDecision::Invalid;
				if (plane_index == 0) return normal_sign < 0 ? CoplanarDecision::Include : CoplanarDecision::Invalid;
				if (plane_index == background.Spec().cells[axis]) return normal_sign > 0 ? CoplanarDecision::Include : CoplanarDecision::Invalid;
				// Positive outward component belongs to the lower-index cell;
				// its contact is its high side. Negative is the converse.
				return normal_sign > 0 ? (side == 1 ? CoplanarDecision::Include : CoplanarDecision::Skip)
					: (side == 0 ? CoplanarDecision::Include : CoplanarDecision::Skip);
			}
		}
		return CoplanarDecision::Include;
	}
	static double TriangleArea(const std::array<double,3>& a, const std::array<double,3>& b,
		const std::array<double,3>& c)
	{
		const long double ux=static_cast<long double>(b[0])-a[0], uy=static_cast<long double>(b[1])-a[1], uz=static_cast<long double>(b[2])-a[2], vx=static_cast<long double>(c[0])-a[0], vy=static_cast<long double>(c[1])-a[1], vz=static_cast<long double>(c[2])-a[2];
		const long double nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
		const long double area=.5L*std::sqrt(static_cast<long double>(nx)*nx+static_cast<long double>(ny)*ny+static_cast<long double>(nz)*nz);
		if (!std::isfinite(area)) throw std::runtime_error("clipped triangle area is nonfinite");
		return area;
	}
	static double Ulp(double value)
	{
		if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();
		return std::abs(std::nextafter(value,std::numeric_limits<double>::infinity())-value);
	}
	static constexpr std::array<std::array<double,4>,12> kDunavant{{
		{{.501426509658179,.249286745170910,.249286745170910,.116786275726379}},
		{{.249286745170910,.501426509658179,.249286745170910,.116786275726379}},
		{{.249286745170910,.249286745170910,.501426509658179,.116786275726379}},
		{{.873821971016996,.063089014491502,.063089014491502,.050844906370207}},
		{{.063089014491502,.873821971016996,.063089014491502,.050844906370207}},
		{{.063089014491502,.063089014491502,.873821971016996,.050844906370207}},
		{{.053145049844817,.310352451033784,.636502499121399,.082851075618374}},
		{{.053145049844817,.636502499121399,.310352451033784,.082851075618374}},
		{{.310352451033784,.053145049844817,.636502499121399,.082851075618374}},
		{{.310352451033784,.636502499121399,.053145049844817,.082851075618374}},
		{{.636502499121399,.053145049844817,.310352451033784,.082851075618374}},
		{{.636502499121399,.310352451033784,.053145049844817,.082851075618374}} }};
	static bool SameGridSpec(const CubicCartesianGridSpec& a, const CubicCartesianGridSpec& b)
	{ return a.lower_m == b.lower_m && a.upper_m == b.upper_m && a.cells == b.cells; }
	static void ValidateBoundaryIds(const ClosedTriangulatedSurface& surface)
	{
		for (const auto& triangle : surface.Triangles()) if (triangle.boundary_id > static_cast<std::uint32_t>(std::numeric_limits<int>::max()))
			throw std::invalid_argument("surface boundary id exceeds int range");
	}
	void MarkCatalogUnusable(const std::string& reason)
	{
		diagnostics_.catalog_usable = false;
		if (diagnostics_.catalog_unusable_reason.empty()) diagnostics_.catalog_unusable_reason = reason;
	}
	static void CountAttempt(std::size_t& attempts, std::size_t cap, const char* message)
	{
		if (attempts == std::numeric_limits<std::size_t>::max())
			throw std::overflow_error("immersed surface quadrature attempt count overflows");
		++attempts;
		if (attempts > cap) throw std::runtime_error(message);
	}
	void RollbackCell(const ImmersedSurfaceQuadratureDiagnostics& diagnostics_before,
		const std::vector<CompensatedArea>& triangle_area_before,
		const std::map<std::uint32_t, CompensatedArea>& label_area_before,
		const std::vector<std::size_t>& triangle_counts_before,
		const std::map<std::uint32_t, std::size_t>& label_counts_before,
		const CompensatedArea& total_area_before, ImmersedSurfaceQuadratureCell& result,
		bool predicate_failure, const std::string& failure_reason = {})
	{
		const std::size_t candidate_attempts = diagnostics_.candidate_attempts;
		const std::size_t fragment_attempts = diagnostics_.fragment_attempts;
		const std::size_t point_attempts = diagnostics_.point_attempts;
		// Exact work performed by a failed cell is still part of this catalog's
		// resource history.  In particular, a retry must not make a high-water
		// mark disappear.
		const std::size_t peak_exact_limbs = diagnostics_.peak_exact_limbs;
		diagnostics_ = diagnostics_before;
		diagnostics_.candidate_attempts = candidate_attempts;
		diagnostics_.fragment_attempts = fragment_attempts;
		diagnostics_.point_attempts = point_attempts;
		diagnostics_.peak_exact_limbs = std::max(diagnostics_.peak_exact_limbs,peak_exact_limbs);
		triangle_area_accumulated_ = triangle_area_before; label_area_accumulated_ = label_area_before;
		triangle_fragment_counts_ = triangle_counts_before; label_fragment_counts_ = label_counts_before;
		total_area_accumulated_ = total_area_before;
		result.usable = false; result.ambiguous = true; result.area_m2 = 0.0; result.point_count = 0;
		if (predicate_failure) ++diagnostics_.predicate_ambiguities;
		else ++diagnostics_.precision_limited_fragments;
		MarkCatalogUnusable(predicate_failure ? "exact surface clipping predicate failed"
			: (failure_reason.empty() ? "a cell failed surface-rule validation" : failure_reason));
	}
	void AddArea(std::size_t triangle_id, std::uint32_t boundary_id, long double area, ImmersedSurfaceQuadratureCell& cell)
	{
		triangle_area_accumulated_[triangle_id].Add(area); total_area_accumulated_.Add(area);
		label_area_accumulated_[boundary_id].Add(area); ++triangle_fragment_counts_[triangle_id]; ++label_fragment_counts_[boundary_id]; cell.area_m2 = static_cast<double>(static_cast<long double>(cell.area_m2)+area);
	}
	void ValidateBinding(const CartesianDomainClassification& domain, std::uint64_t id) const
	{
		if (!SameGridSpec(grid_spec_,domain.Background().Spec()) || cells_.size() != domain.Cells().size()
			|| surface_canonical_hash_ != domain.SurfaceCanonicalHash()) throw std::runtime_error("immersed surface quadrature domain binding does not match catalog");
		if (id >= cells_.size() || cells_[static_cast<std::size_t>(id)].id != id || domain.Cells()[static_cast<std::size_t>(id)].id != id)
			throw std::runtime_error("immersed surface quadrature id is not x-fast in bound domain");
	}
	void BuildCell(const CartesianDomainClassification& domain, const CartesianDomainCell& source,
		ImmersedSurfaceQuadratureCell& result, SurfaceQuadratureRule& rule)
	{
		if (source.classification != CellClassification::Cut) return;
		if (source.ambiguous) { result.usable = false; result.ambiguous = true; ++diagnostics_.predicate_ambiguities; return; }
		const auto& background = domain.Background(); const auto cell = background.Cell(source.id);
		const auto& surface = domain.SurfaceIndex().Surface(); std::vector<SurfaceQuadraturePoint> points;
		std::vector<ImmersedSurfaceQuadraturePointProvenance> provenance;
		const ImmersedSurfaceQuadratureDiagnostics diagnostics_before = diagnostics_;
		const auto triangle_area_before = triangle_area_accumulated_; const auto label_area_before = label_area_accumulated_; const auto triangle_counts_before = triangle_fragment_counts_; const auto label_counts_before = label_fragment_counts_;
		const CompensatedArea total_area_before = total_area_accumulated_;
		try {
			for (const std::size_t triangle_id : source.triangle_ids) {
				if (triangle_id >= surface.Triangles().size()) throw std::runtime_error("domain surface triangle id is out of range");
				CountAttempt(diagnostics_.candidate_attempts, options_.max_candidates, "immersed surface quadrature candidate cap reached");
				++diagnostics_.candidate_pairs; const auto& triangle = surface.Triangles()[triangle_id];
				std::array<std::array<double,3>,3> original{{surface.Vertices()[triangle.indices[0]],surface.Vertices()[triangle.indices[1]],surface.Vertices()[triangle.indices[2]]}};
					const CoplanarDecision ownership = CoplanarOwnership(original,background,source,triangle);
					if (ownership == CoplanarDecision::Invalid) throw std::runtime_error("coplanar surface triangle has no outward boundary owner");
					if (ownership == CoplanarDecision::Skip) { ++diagnostics_.coplanar_owner_skips; continue; }
					const std::array<BoxPlane,6> planes{{{0,cell.lower_m[0]},{0,cell.upper_m[0]},{1,cell.lower_m[1]},{1,cell.upper_m[1]},{2,cell.lower_m[2]},{2,cell.upper_m[2]}}};
					std::vector<HomogeneousPoint> polygon{{Point(original[0],original,planes,0),Point(original[1],original,planes,1),Point(original[2],original,planes,2)}};
					for (std::size_t plane=0; plane<planes.size(); ++plane) polygon=ClipPlane(polygon,plane,(plane%2)==0,original,planes);
					if (polygon.empty()) continue;
					++diagnostics_.clipped_polygons;
					const auto normal=TriangleNormal(original); std::size_t projection=0; while(projection<3&&!Sign(normal[projection])) ++projection; if(projection==3) throw std::runtime_error("surface triangle has zero exact normal");
					RemoveCollinear(polygon,projection);
				if (polygon.size() < 3) { ++diagnostics_.zero_fragments; continue; }
					for (std::size_t fan=1; fan+1<polygon.size(); ++fan) {
						if (!Sign(ProjectedDet(polygon[0],polygon[fan],polygon[fan+1],projection))) { ++diagnostics_.zero_fragments; continue; }
						CountAttempt(diagnostics_.fragment_attempts, options_.max_fragments, "immersed surface quadrature fragment cap reached");
						// Convert exact homogeneous vertices directly to cell-local
						// coordinates.  Constructing an absolute double here loses
						// sub-ULP local position near a large grid origin, and recovering
						// it later would incorrectly change the reference coordinates.
						std::array<std::array<double,3>,3> local{};
						std::array<const HomogeneousPoint*,3> vertices{{&polygon[0],&polygon[fan],&polygon[fan+1]}};
						for (std::size_t n=0;n<3;++n) for(std::size_t a=0;a<3;++a) {
							const std::size_t lower_plane=2*a, upper_plane=lower_plane+1;
							const double h=cell.upper_m[a]-cell.lower_m[a];
							if (vertices[n]->active&(1u<<(3+lower_plane))) local[n][a]=0.0;
							else if (vertices[n]->active&(1u<<(3+upper_plane))) local[n][a]=h;
							else local[n][a]=ToDouble(S(vertices[n]->value[a],M(D(cell.lower_m[a]),vertices[n]->value[3])),vertices[n]->value[3]);
							if(local[n][a]==0.0) local[n][a]=0.0;
						}
						std::array<bool,3> fragment_on_lower{}, fragment_on_upper{};
						for (std::size_t a=0; a<3; ++a) {
							fragment_on_lower[a] = (vertices[0]->active&(1u<<(3+2*a))) && (vertices[1]->active&(1u<<(3+2*a))) && (vertices[2]->active&(1u<<(3+2*a)));
							fragment_on_upper[a] = (vertices[0]->active&(1u<<(3+2*a+1))) && (vertices[1]->active&(1u<<(3+2*a+1))) && (vertices[2]->active&(1u<<(3+2*a+1)));
						}
						double area=TriangleArea(local[0],local[1],local[2]);
#ifdef IGA_EXACT_DYADIC_TESTING
					if (SurfaceQuadratureTestControls().reject_positive_area) area = 0.0;
#endif
					if (!(area>0.0)) throw std::runtime_error("positive clipped surface fragment is not representable");
						++diagnostics_.positive_fragments; AddArea(triangle_id,triangle.boundary_id,area,result);
					for (const auto& q : kDunavant) {
						CountAttempt(diagnostics_.point_attempts, options_.max_points, "immersed surface quadrature point cap reached");
						std::array<double,3> local_point{}, parametric{}, physical{};
						for(std::size_t a=0;a<3;++a)
							local_point[a]=local[0][a]+q[1]*(local[1][a]-local[0][a])+q[2]*(local[2][a]-local[0][a]);
						std::array<std::array<double,3>,3> original_local{};
						for (std::size_t n=0;n<3;++n) for (std::size_t a=0;a<3;++a)
							original_local[n][a]=ToDouble(S(D(original[n][a]),M(D(cell.lower_m[a]),D(1.0))),D(1.0));
						double surface_residual = 0.0, surface_scale = 0.0;
						for (std::size_t a=0; a<3; ++a) {
							surface_residual += (local_point[a]-original_local[0][a])*triangle.outward_unit_normal[a];
							surface_scale = std::max(surface_scale,std::max(std::abs(original_local[1][a]-original_local[0][a]),std::abs(original_local[2][a]-original_local[0][a])));
						}
						surface_residual=std::abs(surface_residual);
						diagnostics_.max_surface_residual_m=std::max(diagnostics_.max_surface_residual_m,surface_residual);
						double surface_ulp = 0.0;
						for (std::size_t a=0; a<3; ++a)
							surface_ulp = std::max(surface_ulp,std::max(Ulp(local_point[a]),Ulp(original_local[0][a])));
						if (surface_residual > 128*std::numeric_limits<double>::epsilon()*surface_scale+8.0*surface_ulp)
							throw std::runtime_error("immersed surface quadrature point leaves source triangle plane");
						for(std::size_t a=0;a<3;++a) {
							const double h=cell.upper_m[a]-cell.lower_m[a];
							parametric[a]=local_point[a]/h;
							if (fragment_on_lower[a]) parametric[a]=0.0;
							else if (fragment_on_upper[a]) parametric[a]=1.0;
							const double mapping_residual=std::abs(local_point[a]-parametric[a]*h);
							const double mapping_scale=std::max(std::abs(h),std::numeric_limits<double>::denorm_min());
							if (mapping_residual>128.0*std::numeric_limits<double>::epsilon()*mapping_scale+8.0*Ulp(local_point[a]))
								throw std::runtime_error("immersed surface local point does not map through Cartesian cell");
							physical[a]=cell.lower_m[a]+local_point[a];
						}
						const double weight=area*q[3]; if (!std::isfinite(weight) || !(weight>0.0)) throw std::runtime_error("immersed surface quadrature weight is not representable");
						std::array<double,3> barycentric{{0.0,0.0,0.0}};
						for (std::size_t corner=0; corner<3; ++corner)
							barycentric[corner]=q[0]*ToDouble(vertices[0]->material[corner],vertices[0]->material[3])
								+q[1]*ToDouble(vertices[1]->material[corner],vertices[1]->material[3])
								+q[2]*ToDouble(vertices[2]->material[corner],vertices[2]->material[3]);
						ValidateProvenancePoint(surface, triangle_id, triangle, cell.lower_m, local_point, barycentric);
						points.push_back({parametric,physical,triangle.outward_unit_normal,weight,static_cast<int>(triangle.boundary_id)});
						provenance.push_back({static_cast<std::uint32_t>(triangle_id), barycentric}); ++diagnostics_.output_points;
					}
				}
			}
			const Element element=background.MaterializeElement(source.id);
				rule=SurfaceQuadratureRule(std::move(points)); result.point_count=rule.Points().size();
				if (provenance.size() != result.point_count) throw std::runtime_error("immersed surface provenance is not one-to-one");
				provenance_[static_cast<std::size_t>(source.id)] = std::move(provenance);
				ValidateCellRule(element,result,rule);
		} catch (const std::overflow_error&) { provenance_[static_cast<std::size_t>(source.id)].clear(); rule=SurfaceQuadratureRule{}; RollbackCell(diagnostics_before,triangle_area_before,label_area_before,triangle_counts_before,label_counts_before,total_area_before,result,true); return; }
		catch (const std::exception& error) {
			const std::string message(error.what());
			if (message.find("cap reached") != std::string::npos) throw;
			provenance_[static_cast<std::size_t>(source.id)].clear(); rule=SurfaceQuadratureRule{}; RollbackCell(diagnostics_before,triangle_area_before,label_area_before,triangle_counts_before,label_counts_before,total_area_before,result,false,message); return;
		}
	}
	static void ValidateProvenancePoint(const ClosedTriangulatedSurface& surface, std::size_t triangle_id,
		const ClosedSurfaceTriangle& triangle, const std::array<double,3>& cell_lower,
		const std::array<double,3>& local_point, const std::array<double, 3>& barycentric)
	{
		if (triangle_id >= surface.Triangles().size() || triangle_id > std::numeric_limits<std::uint32_t>::max()) throw std::out_of_range("immersed surface provenance triangle is out of range");
		double sum = 0.0, scale = 0.0, residual = 0.0;
		std::array<double, 3> reconstruction{{0.0, 0.0, 0.0}};
		for (std::size_t corner = 0; corner < 3; ++corner) {
			if (!std::isfinite(barycentric[corner]) || barycentric[corner] < -1.0e-10 || barycentric[corner] > 1.0+1.0e-10) throw std::runtime_error("immersed surface provenance barycentric coordinate is outside simplex");
			sum += barycentric[corner];
			const auto& vertex = surface.Vertices()[triangle.indices[corner]];
			for (std::size_t axis = 0; axis < 3; ++axis) { const double local_vertex=vertex[axis]-cell_lower[axis]; reconstruction[axis] += barycentric[corner]*local_vertex; scale = std::max(scale, std::abs(local_vertex)); }
		}
		for (std::size_t axis = 0; axis < 3; ++axis) residual = std::max(residual, std::abs(reconstruction[axis]-local_point[axis]));
		if (!std::isfinite(sum) || std::abs(sum-1.0) > 2.0e-10 || residual > 2.0e-10*std::max(std::numeric_limits<double>::min(), scale)
			|| triangle.boundary_id > static_cast<std::uint32_t>(std::numeric_limits<int>::max()))
			throw std::runtime_error("immersed surface provenance reconstruction or label is inconsistent");
	}
	void ValidateCellRule(const Element& element, const ImmersedSurfaceQuadratureCell& cell, const SurfaceQuadratureRule& rule)
	{
		#ifdef IGA_EXACT_DYADIC_TESTING
		auto& controls=SurfaceQuadratureTestControls();
		if (controls.reject_rule_validation
			&& controls.rule_validation_attempts++ == controls.reject_rule_validation_after) {
			controls.reject_rule_validation=false; // exactly one late failure
			throw std::runtime_error("immersed surface quadrature test validation rejection");
		}
		#endif
			long double sum=0.0L; for (const auto& point:rule.Points()) {
			const auto geometry=EvaluateElementGeometry(element,point.parametric); if (!std::isfinite(geometry.raw_determinant) || !(geometry.raw_determinant > 0.0)) throw std::runtime_error("immersed surface point has a non-positive raw element Jacobian"); double residual=0.0, scale=0.0, ulp=0.0;
			for(std::size_t a=0;a<3;++a) {
				const double lower=element.bezier_points[0][a], upper=element.bezier_points[63][a], expected=lower+point.parametric[a]*(upper-lower);
				residual=std::max(residual,std::abs(expected-point.physical[a]));
				scale=std::max(scale,std::abs(upper-lower));
				ulp=std::max(ulp,std::max(Ulp(expected),Ulp(point.physical[a])));
			}
			diagnostics_.max_mapping_residual_m=std::max(diagnostics_.max_mapping_residual_m,residual); if (residual>128*std::numeric_limits<double>::epsilon()*scale+8.0*ulp) throw std::runtime_error("immersed surface point does not map through Cartesian element");
			sum+=point.weight;
		}
			const long double scale=std::abs(static_cast<long double>(cell.area_m2))+std::abs(sum);
			if (std::abs(sum-static_cast<long double>(cell.area_m2))>64.0L*std::numeric_limits<double>::epsilon()*std::max(scale,static_cast<long double>(std::numeric_limits<double>::denorm_min()))) throw std::runtime_error("immersed surface quadrature weights do not sum to clipped area");
			ValidateSurfaceQuadratureRule(element,rule);
	}
	void FinalizeAudit(const ClosedTriangulatedSurface& surface)
	{
		std::map<std::uint32_t,CompensatedArea> source_by_label;
		CompensatedArea source_total;
		for(std::size_t i=0;i<surface.Triangles().size();++i) {
			const long double source=surface.Triangles()[i].area_m2, accumulated=triangle_area_accumulated_[i].Value(), residual=accumulated-source;
			diagnostics_.triangle_area_m2[i]=static_cast<double>(accumulated);
			diagnostics_.triangle_area_residual_m2[i]=static_cast<double>(residual);
			diagnostics_.triangle_area_absolute_residual_m2[i]=static_cast<double>(std::abs(residual));
			source_by_label[surface.Triangles()[i].boundary_id].Add(source); source_total.Add(source);
			const long double scale=std::abs(source)+std::abs(accumulated);
			const long double operations = 64.0L+16.0L*triangle_fragment_counts_[i];
			if (std::abs(residual)>operations*std::numeric_limits<double>::epsilon()*std::max(scale,static_cast<long double>(std::numeric_limits<double>::denorm_min()))) MarkCatalogUnusable("per-triangle surface-area audit failed");
		}
		const long double accumulated_total=total_area_accumulated_.Value(), source_total_value=source_total.Value(), total_residual=accumulated_total-source_total_value;
		diagnostics_.source_total_area_m2=static_cast<double>(source_total_value);
		diagnostics_.accumulated_total_area_m2=static_cast<double>(accumulated_total);
		diagnostics_.total_area_m2=diagnostics_.accumulated_total_area_m2;
		diagnostics_.total_area_residual_m2=static_cast<double>(total_residual);
		diagnostics_.total_area_absolute_residual_m2=static_cast<double>(std::abs(total_residual));
		const long double total_scale=std::abs(source_total_value)+std::abs(accumulated_total);
		const long double total_operations=64.0L+16.0L*std::accumulate(triangle_fragment_counts_.begin(),triangle_fragment_counts_.end(),std::size_t{0});
		if (std::abs(total_residual)>total_operations*std::numeric_limits<double>::epsilon()*std::max(total_scale,static_cast<long double>(std::numeric_limits<double>::denorm_min()))) MarkCatalogUnusable("total surface-area audit failed");
		for (const auto& item:source_by_label) {
			const long double source=item.second.Value(), actual=label_area_accumulated_[item.first].Value(), residual=actual-source, scale=std::abs(actual)+std::abs(source);
			diagnostics_.area_by_boundary_id[item.first]=static_cast<double>(actual);
			diagnostics_.source_area_by_boundary_id[item.first]=static_cast<double>(source);
			diagnostics_.area_residual_by_boundary_id[item.first]=static_cast<double>(residual);
			diagnostics_.area_absolute_residual_by_boundary_id[item.first]=static_cast<double>(std::abs(residual));
			const long double operations = 64.0L+16.0L*label_fragment_counts_[item.first];
			if (std::abs(residual)>operations*std::numeric_limits<double>::epsilon()*std::max(scale,static_cast<long double>(std::numeric_limits<double>::denorm_min()))) MarkCatalogUnusable("per-label surface-area audit failed");
		}
	}

	CubicCartesianGridSpec grid_spec_{}; std::string surface_canonical_hash_; ImmersedSurfaceQuadratureOptions options_{};
	std::vector<ImmersedSurfaceQuadratureCell> cells_; std::vector<SurfaceQuadratureRule> rules_; std::vector<std::vector<ImmersedSurfaceQuadraturePointProvenance>> provenance_; ImmersedSurfaceQuadratureDiagnostics diagnostics_;
	std::vector<CompensatedArea> triangle_area_accumulated_; std::map<std::uint32_t,CompensatedArea> label_area_accumulated_; std::vector<std::size_t> triangle_fragment_counts_; std::map<std::uint32_t,std::size_t> label_fragment_counts_; CompensatedArea total_area_accumulated_;
};

} // namespace iga

#endif

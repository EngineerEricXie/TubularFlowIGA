#ifndef IGA_FITTED_CUT_CELL_VOLUME_RULE_HPP
#define IGA_FITTED_CUT_CELL_VOLUME_RULE_HPP

#include "NonnegativeLeastSquares.hpp"
#include "PolyhedralBoxVolumeMoments.hpp"
#include "Quadrature.hpp"
#include "SurfaceSpatialIndex.hpp"
#include <cstdint>
#include <iomanip>
#include <sstream>

namespace iga {

struct FittedCutCellVolumeRuleOptions {
	std::vector<unsigned> candidate_orders{16,24,32};
	double support_expansion=1.25;
	std::size_t max_point_queries=65536;
	std::size_t max_seed_points=1000000;
	NonnegativeLeastSquaresOptions fit;
	PolyhedralBoxMomentOptions moments;
};

struct FittedCutCellVolumeRuleResult {
	VolumeQuadratureRule rule;
	std::size_t point_queries=0,candidates=0,attempts=0;
	long double reference_volume=0,relative_moment_residual=0;
	NonnegativeLeastSquaresResult fit;
};

inline void ValidateFittedCutCellVolumeRuleOptions(const FittedCutCellVolumeRuleOptions& options)
{
	if(options.candidate_orders.empty()||!options.max_point_queries||!options.max_seed_points
		||!std::isfinite(options.support_expansion)||options.support_expansion<1
		||options.fit.max_rows<343||!options.fit.max_columns||!options.fit.max_workspace_bytes||!options.fit.max_iterations
		||!std::isfinite(options.fit.relative_tolerance)||!(options.fit.relative_tolerance>0)
		||!std::isfinite(options.fit.absolute_tolerance)||options.fit.absolute_tolerance<0
		||!options.moments.max_triangles||!options.moments.max_quadrature_points)
		throw std::invalid_argument("invalid fitted cut-cell options");
	unsigned previous=0;
	for(auto order:options.candidate_orders) {
		if(order<2||order>64||order<=previous)throw std::invalid_argument("fitted cut-cell candidate orders must increase within [2,64]");
		previous=order;
	}
}

// Build a positive tensor-degree-six rule on a cut cell. The seed supplies
// only the support frame, never an acceptance fallback. Geometry predicates
// select every new node. Catalog publication and its diagnostics/hash updates
// remain the caller's responsibility; no catalog is mutated here.
inline FittedCutCellVolumeRuleResult BuildFittedCutCellVolumeRule(
	const SurfaceSpatialIndex& surface,const std::array<double,3>& lower,
	const std::array<double,3>& upper,const VolumeQuadratureRule& seed,
	FittedCutCellVolumeRuleOptions options={})
{
	ValidateFittedCutCellVolumeRuleOptions(options);
	if(seed.Points().empty())throw std::invalid_argument("fitted cut-cell seed is empty");
	if(seed.Points().size()>options.max_seed_points)throw std::runtime_error("fitted cut-cell seed cap reached");
	std::array<double,3> minimum{{1,1,1}},maximum{{0,0,0}};
	std::array<long double,3> extent;
	long double determinant=1;
	for(unsigned q=0;q<3;++q) {
		if(!std::isfinite(lower[q])||!std::isfinite(upper[q])||!(upper[q]>lower[q]))throw std::invalid_argument("invalid fitted cut-cell bounds");
		extent[q]=static_cast<long double>(upper[q])-lower[q];determinant*=extent[q];
	}
	if(!std::isfinite(determinant)||!(determinant>0))throw std::overflow_error("invalid fitted cut-cell determinant");
	for(const auto& point:seed.Points()) {
		if(!std::isfinite(point.weight)||!(point.weight>0))throw std::invalid_argument("invalid fitted cut-cell seed weight");
		for(unsigned q=0;q<3;++q) {
			if(!std::isfinite(point.parametric[q])||point.parametric[q]<0||point.parametric[q]>1)throw std::invalid_argument("invalid fitted cut-cell seed point");
			minimum[q]=std::min(minimum[q],point.parametric[q]);maximum[q]=std::max(maximum[q],point.parametric[q]);
		}
	}
	FittedCutCellVolumeRuleResult result;
	if(minimum[0]==maximum[0]||minimum[1]==maximum[1]||minimum[2]==maximum[2]) {
		// A collapsed sample cloud is not a zero-volume certificate. Recover
		// support from clipped boundary facets and retained box corners. Every
		// extremum of a polyhedron/box intersection is represented by these.
		std::array<long double,3> clipped_min{{1,1,1}},clipped_max{{0,0,0}};
		const auto include=[&](const std::array<long double,3>& point) {
			for(unsigned q=0;q<3;++q) {
				clipped_min[q]=std::min(clipped_min[q],point[q]);
				clipped_max[q]=std::max(clipped_max[q],point[q]);
			}
		};
		if(surface.Surface().Triangles().size()>options.moments.max_triangles)
			throw std::runtime_error("fitted support triangle cap reached");
		for(const auto& triangle:surface.Surface().Triangles()) {
			polyhedral_box_moment_detail::Polygon polygon;
			for(auto vertex:triangle.indices) {
				std::array<long double,3> point;
				for(unsigned q=0;q<3;++q) {
					point[q]=(static_cast<long double>(surface.Surface().Vertices()[vertex][q])-lower[q])/extent[q];
					if(!std::isfinite(point[q]))throw std::overflow_error("fitted support coordinate is nonfinite");
				}
				polygon.push_back(point);
			}
			for(unsigned q=0;q<3;++q) {
				polygon=polyhedral_box_moment_detail::Clip(polygon,q,0,true);
				polygon=polyhedral_box_moment_detail::Clip(polygon,q,1,false);
			}
			for(const auto& point:polygon)include(point);
		}
		for(unsigned corner=0;corner<8;++corner) {
			std::array<long double,3> point;std::array<double,3> physical;
			for(unsigned q=0;q<3;++q) { point[q]=(corner>>q)&1;physical[q]=point[q]?upper[q]:lower[q]; }
			if(result.point_queries==options.max_point_queries)throw std::runtime_error("fitted cut-cell point-query cap reached");
			++result.point_queries;
			const auto location=surface.LocatePoint(physical);
			if(location==PointLocation::Ambiguous)throw std::runtime_error("fitted support corner classification is ambiguous");
			if(location!=PointLocation::Outside)include(point);
		}
		for(unsigned q=0;q<3;++q)if(minimum[q]==maximum[q]) {
			minimum[q]=static_cast<double>(clipped_min[q]);maximum[q]=static_cast<double>(clipped_max[q]);
		}
	}
	for(unsigned q=0;q<3;++q) {
		if(!(maximum[q]>minimum[q]))throw std::runtime_error("fitted cut-cell clipped support has zero width");
		options.moments.coordinate_origin[q]=(maximum[q]+minimum[q])/2;
		options.moments.coordinate_scale[q]=(maximum[q]-minimum[q])/2;
	}

	std::vector<std::array<unsigned,3>> exponents;
	std::vector<double> target;
	std::vector<long double> precise_target;
	for(unsigned a=0;a<=6;++a)for(unsigned b=0;b<=6;++b)for(unsigned c=0;c<=6;++c) {
		exponents.push_back({{a,b,c}});
		const auto moment=PolyhedralBoxVolumeMoment(surface.Surface(),lower,upper,{{a,b,c}},options.moments)/determinant;
		if(!std::isfinite(moment)||std::abs(moment)>std::numeric_limits<double>::max())throw std::overflow_error("fitted cut-cell moment is not representable");
		precise_target.push_back(moment);target.push_back(static_cast<double>(moment));
	}
	if(precise_target[0]>1+options.fit.absolute_tolerance+options.fit.relative_tolerance)throw std::runtime_error("fitted cut-cell target volume exceeds its box");
	if(!(precise_target[0]>0))throw std::runtime_error("fitted cut-cell target volume is not positive");
	// Clip the candidate interval BEFORE placing Gauss nodes. Dropping nodes
	// from an overhanging tensor grid can leave too few distinct coordinates
	// near a cell face to reproduce degree-six moments, even with enough points.
	std::array<long double,3> candidate_lower,candidate_upper;
	for(unsigned q=0;q<3;++q) {
		const long double radius=static_cast<long double>(options.support_expansion)*options.moments.coordinate_scale[q];
		candidate_lower[q]=std::max(0.L,options.moments.coordinate_origin[q]-radius);
		candidate_upper[q]=std::min(1.L,options.moments.coordinate_origin[q]+radius);
	}
	for(auto order:options.candidate_orders) {
		++result.attempts;
		const auto gauss=polyhedral_moment_detail::Gauss(order);
		std::vector<VolumeQuadraturePoint> candidates;
		std::vector<double> matrix;
		std::uint64_t inside_candidates=0;
		for(const auto& gx:gauss)for(const auto& gy:gauss)for(const auto& gz:gauss) {
			const std::array<long double,3> g{{gx.first,gy.first,gz.first}};
			VolumeQuadraturePoint point;std::array<double,3> physical;
			bool in_box=true;
			for(unsigned q=0;q<3;++q) {
				point.parametric[q]=static_cast<double>(candidate_lower[q]+g[q]*(candidate_upper[q]-candidate_lower[q]));
				if(point.parametric[q]<0||point.parametric[q]>1)in_box=false;
				physical[q]=static_cast<double>(lower[q]+extent[q]*point.parametric[q]);
			}
			if(!in_box)continue;
			if(result.point_queries==options.max_point_queries)throw std::runtime_error("fitted cut-cell point-query cap reached");
			++result.point_queries;
			const auto location=surface.LocatePoint(physical);
			if(location==PointLocation::Ambiguous)throw std::runtime_error("fitted cut-cell point classification is ambiguous");
			if(location!=PointLocation::Inside)continue;
			++inside_candidates;
			if(candidates.size()<options.fit.max_columns)candidates.push_back(point);
			else {
				// A deterministic reservoir keeps the matrix within max_columns
				// without truncating one end of the tensor grid. Query limits still
				// count every classified proposal, and the full moment audit is
				// mandatory for the retained subset.
				std::uint64_t sample=inside_candidates+0x9e3779b97f4a7c15ULL;
				sample=(sample^(sample>>30))*0xbf58476d1ce4e5b9ULL;
				sample=(sample^(sample>>27))*0x94d049bb133111ebULL;
				sample^=sample>>31;
				const auto slot=sample%inside_candidates;
				if(slot<candidates.size())candidates[static_cast<std::size_t>(slot)]=point;
			}
		}
		if(candidates.empty())continue;
		// Reservoir replacements disturb tensor traversal order. Restore the
		// canonical point order required by compact quadrature publication.
		std::sort(candidates.begin(),candidates.end(),[](const VolumeQuadraturePoint& a,const VolumeQuadraturePoint& b) {
			return a.parametric<b.parametric;
		});
		// Allocate the matrix once, so vector capacity growth cannot exceed
		// the checked caller-owned matrix bound.
		const auto entries=nonnegative_least_squares_detail::Product(candidates.size(),exponents.size());
		if(nonnegative_least_squares_detail::Product(entries,sizeof(double))>options.fit.max_workspace_bytes)
			throw std::runtime_error("fitted cut-cell matrix cap reached");
		matrix.resize(entries);
		for(std::size_t j=0;j<candidates.size();++j)for(std::size_t i=0;i<exponents.size();++i) {
			long double value=1;
			for(unsigned q=0;q<3;++q)value*=polyhedral_moment_detail::Power(
				(static_cast<long double>(candidates[j].parametric[q])-options.moments.coordinate_origin[q])/options.moments.coordinate_scale[q],exponents[i][q]);
			matrix[j*exponents.size()+i]=static_cast<double>(value);
		}
		result.fit=FitNonnegativeLeastSquares(matrix,target,options.fit);
		result.candidates=candidates.size();
		if(!result.fit.within_tolerance)continue;
		std::vector<VolumeQuadraturePoint> retained;
		retained.reserve(exponents.size());
		for(std::size_t j=0;j<candidates.size();++j)if(result.fit.coefficients[j]>0) {
			candidates[j].weight=result.fit.coefficients[j];retained.push_back(candidates[j]);
		}
		std::vector<long double> residual=precise_target;
		for(const auto& point:retained)for(std::size_t i=0;i<exponents.size();++i) {
			long double value=point.weight;
			for(unsigned q=0;q<3;++q)value*=polyhedral_moment_detail::Power(
				(static_cast<long double>(point.parametric[q])-options.moments.coordinate_origin[q])/options.moments.coordinate_scale[q],exponents[i][q]);
			residual[i]-=value;
		}
		const auto norm=nonnegative_least_squares_detail::Norm(residual),target_norm=nonnegative_least_squares_detail::Norm(precise_target);
		result.relative_moment_residual=norm/target_norm;
		if(!std::isfinite(result.relative_moment_residual))throw std::runtime_error("fitted cut-cell moment audit is not finite");
		if(norm>options.fit.absolute_tolerance+options.fit.relative_tolerance*target_norm)continue;
		result.reference_volume=precise_target[0];result.rule=VolumeQuadratureRule(std::move(retained));
		return result;
	}
	std::ostringstream message;message<<std::setprecision(std::numeric_limits<long double>::max_digits10)
		<<"fitted cut-cell candidate orders exhausted without a positive accurate rule; bounds=";
	for(unsigned q=0;q<3;++q)message<<'['<<lower[q]<<','<<upper[q]<<']';
	message<<" support=";for(unsigned q=0;q<3;++q)message<<'['<<minimum[q]<<','<<maximum[q]<<']';
	message<<" candidates="<<result.candidates<<" iterations="<<result.fit.iterations
		<<" positive="<<result.fit.positive_columns<<" rank_rejections="<<result.fit.rank_rejections
		<<" nnls_relative_residual="<<result.fit.relative_residual
		<<" moment_relative_residual="<<result.relative_moment_residual;
	throw std::runtime_error(message.str());
}

}

#endif

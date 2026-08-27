#include "BSpline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace tubular {
namespace {

using Value = std::array<double,4>;

Value Add(Value a, const Value& b)
{
	for (int i=0; i<4; ++i) a[i] += b[i];
	return a;
}

Value Scale(Value a, double s)
{
	for (double& value : a) value *= s;
	return a;
}

Value Evaluate(const std::vector<double>& knots, const std::vector<Value>& control, int degree, double u)
{
	const int n = static_cast<int>(control.size())-1;
	if (n < degree || static_cast<int>(knots.size()) != n+degree+2)
		throw std::runtime_error("invalid B-spline representation");
	int span = degree;
	if (u >= knots[n+1]) span = n;
	else {
		for (int i=degree; i<=n; ++i)
			if (u >= knots[i] && u < knots[i+1]) { span = i; break; }
	}
	std::vector<Value> d(degree+1);
	for (int j=0; j<=degree; ++j) d[j] = control[span-degree+j];
	for (int r=1; r<=degree; ++r) {
		for (int j=degree; j>=r; --j) {
			const int i = span-degree+j;
			const double denominator = knots[i+degree-r+1]-knots[i];
			const double alpha = denominator > 0.0 ? (u-knots[i])/denominator : 0.0;
			d[j] = Add(Scale(d[j-1], 1.0-alpha), Scale(d[j], alpha));
		}
	}
	return d[degree];
}

struct SplineRepresentation
{
	int degree = 0;
	std::vector<double> knots;
	std::vector<Value> control;
	std::vector<double> derivative_knots;
	std::vector<Value> derivative_control;
};

SplineRepresentation BuildSpline(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	double polyline_length)
{
	const int count = static_cast<int>(points.size());
	SplineRepresentation result;
	result.degree = count > 3 ? 3 : count-1;
	std::vector<double> chord(count, 0.0);
	for (int i=1; i<count; ++i)
		chord[i] = chord[i-1]+Norm(points[i]-points[i-1])/polyline_length;
	chord.back() = 1.0;
	result.knots.assign(result.degree+1, 0.0);
	if (count > 3)
		for (int i=2; i<=count-3; ++i) result.knots.push_back(chord[i]);
	result.knots.insert(result.knots.end(), result.degree+1, 1.0);
	result.control.resize(count);
	for (int i=0; i<count; ++i)
		result.control[i] = {points[i].x, points[i].y, points[i].z, std::log(diameters[i])};
	result.derivative_control.resize(count-1);
	for (int i=0; i<count-1; ++i) {
		const double denominator = result.knots[i+result.degree+1]-result.knots[i+1];
		if (denominator <= 0.0)
			throw std::runtime_error("degenerate B-spline derivative knot span");
		for (int component=0; component<4; ++component)
			result.derivative_control[i][component]
				= result.degree*(result.control[i+1][component]-result.control[i][component])/denominator;
	}
	result.derivative_knots.assign(result.knots.begin()+1, result.knots.end()-1);
	return result;
}

CurveSample EvaluateCurve(const SplineRepresentation& spline, double u)
{
	u = std::max(0.0, std::min(1.0, u));
	const Value value = Evaluate(spline.knots, spline.control, spline.degree, u);
	const Value tangent = Evaluate(
		spline.derivative_knots, spline.derivative_control, spline.degree-1, u);
	const double diameter = std::exp(value[3]);
	return {{value[0],value[1],value[2]},diameter,{tangent[0],tangent[1],tangent[2]}};
}

} // namespace

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	const BranchSamplingOptions& options,
	int mode,
	const std::string& context)
{
	if (points.size() != diameters.size() || points.size() < 2)
		throw std::runtime_error("SampleBranch requires matching point and diameter arrays");
	double polyline_length = 0.0;
	for (std::size_t i=1; i<points.size(); ++i)
		polyline_length += Norm(points[i]-points[i-1]);
	if (!std::isfinite(polyline_length) || polyline_length <= 0.0
		|| !std::isfinite(options.target_spacing) || options.target_spacing <= 0.0
		|| !std::isfinite(options.max_spacing_over_diameter)
		|| options.max_spacing_over_diameter <= 0.0
		|| !std::isfinite(options.max_turn_degrees) || options.max_turn_degrees <= 0.0
		|| !std::isfinite(options.max_diameter_change_fraction)
		|| options.max_diameter_change_fraction <= 0.0)
		throw std::runtime_error("SampleBranch requires positive curve and segment lengths");
	for (double value : diameters)
		if (!std::isfinite(value) || value <= 0.0) throw std::runtime_error("SampleBranch received an invalid diameter");
	if (mode < 1 || mode > 4)
		throw std::runtime_error("unknown branch sampling mode");

	const auto spline = BuildSpline(points,diameters,polyline_length);
	const int dense_count = std::min(200000, std::max({256,
		static_cast<int>(points.size()*64),
		static_cast<int>(std::ceil(polyline_length/options.target_spacing))*16}));
	std::vector<double> dense_arc(dense_count+1,0.0);
	std::vector<CurveSample> dense(dense_count+1);
	for (int i=0;i<=dense_count;++i) {
		dense[i]=EvaluateCurve(spline,static_cast<double>(i)/dense_count);
		if (!IsFinite(dense[i].point) || !std::isfinite(dense[i].diameter)
			|| dense[i].diameter<=0.0 || Norm(dense[i].tangent)<=1.0e-14)
			throw std::runtime_error(context+": B-spline evaluation produced invalid geometry");
		if(i>0) dense_arc[i]=dense_arc[i-1]+Norm(dense[i].point-dense[i-1].point);
	}
	const double length=dense_arc.back();
	if(!std::isfinite(length)||length<=0.0)
		throw std::runtime_error(context+": smoothed branch has invalid arc length");

	auto at_arc=[&](double s) {
		s=std::max(0.0,std::min(length,s));
		auto upper=std::lower_bound(dense_arc.begin(),dense_arc.end(),s);
		if(upper==dense_arc.begin()) return dense.front();
		if(upper==dense_arc.end()) return dense.back();
		const int hi=static_cast<int>(upper-dense_arc.begin());
		const int lo=hi-1;
		const double span=dense_arc[hi]-dense_arc[lo];
		const double alpha=span>0.0?(s-dense_arc[lo])/span:0.0;
		const double u=(lo+alpha)/dense_count;
		return EvaluateCurve(spline,u);
	};

	const double start_clearance=(mode==1||mode==2)
		? options.downstream_clearance_over_diameter*diameters.front():0.0;
	const double end_clearance=(mode==1||mode==3)
		? options.upstream_clearance_over_diameter*diameters.back():0.0;
	const double usable_begin=start_clearance;
	const double usable_end=length-end_clearance;
	const double tolerance=1.0e-10*std::max(1.0,length);
	if(usable_begin>=usable_end-tolerance) {
		std::ostringstream message;
		message<<context<<": insufficient bifurcation clearance; arc_length="<<length
			<<" required_start="<<start_clearance<<" required_end="<<end_clearance;
		throw std::runtime_error(message.str());
	}

	auto local_spacing=[&](double s) {
		const auto center=at_arc(s);
		double spacing=std::min(options.target_spacing,
			options.max_spacing_over_diameter*center.diameter);
		const double delta=std::max(length/dense_count*2.0,length*1.0e-6);
		const double left=std::max(0.0,s-delta),right=std::min(length,s+delta);
		if(right-left>tolerance) {
			const auto a=at_arc(left),b=at_arc(right);
			const double angle=std::acos(ClampUnit(Dot(
				Normalized(a.tangent,"adaptive tangent"),Normalized(b.tangent,"adaptive tangent"))));
			const double curvature=angle/(right-left);
			if(curvature>1.0e-14)
				spacing=std::min(spacing,
					(options.max_turn_degrees*std::acos(-1.0)/180.0)/curvature);
			const double log_rate=std::abs(std::log(b.diameter/a.diameter))/(right-left);
			if(log_rate>1.0e-14)
				spacing=std::min(spacing,
					std::log1p(options.max_diameter_change_fraction)/log_rate);
		}
		if(!std::isfinite(spacing)||spacing<=tolerance)
			throw std::runtime_error(context+": adaptive sampling requested a degenerate spacing");
		return spacing;
	};

	std::vector<double> arc_samples{0.0};
	if(usable_begin>tolerance) arc_samples.push_back(usable_begin);
	double current=usable_begin;
	while(current<usable_end-tolerance) {
		const double spacing=local_spacing(current);
		double next=std::min(usable_end,current+spacing);
		if(usable_end-next<0.25*spacing) next=usable_end;
		if(next<=current+tolerance)
			throw std::runtime_error(context+": adaptive sampling did not advance");
		arc_samples.push_back(next);
		current=next;
		if(arc_samples.size()>100000000)
			throw std::runtime_error(context+": adaptive sampling generated too many points");
	}
	if(length-arc_samples.back()>tolerance) arc_samples.push_back(length);
	const double maximum_turn=options.max_turn_degrees*std::acos(-1.0)/180.0;
	auto interval_is_valid=[&](double begin,double end) {
		if(begin<usable_begin-tolerance||end>usable_end+tolerance) return true;
		const auto first=at_arc(begin),middle=at_arc((begin+end)/2.0),last=at_arc(end);
		const double minimum_diameter=std::min({first.diameter,middle.diameter,last.diameter});
		const double maximum_diameter=std::max({first.diameter,middle.diameter,last.diameter});
		if(end-begin>options.target_spacing*(1.0+1.0e-8)
			||end-begin>options.max_spacing_over_diameter*minimum_diameter*(1.0+1.0e-8)
			||(maximum_diameter-minimum_diameter)/minimum_diameter
				>options.max_diameter_change_fraction*(1.0+1.0e-8)) return false;
		auto tangent_angle=[](const CurveSample& a,const CurveSample& b) {
			return std::acos(ClampUnit(Dot(Normalized(a.tangent,"adaptive interval tangent"),
				Normalized(b.tangent,"adaptive interval tangent"))));
		};
		return tangent_angle(first,middle)<=maximum_turn*(1.0+1.0e-8)
			&&tangent_angle(middle,last)<=maximum_turn*(1.0+1.0e-8)
			&&tangent_angle(first,last)<=maximum_turn*(1.0+1.0e-8);
	};
	for(std::size_t interval=0;interval+1<arc_samples.size();) {
		if(interval_is_valid(arc_samples[interval],arc_samples[interval+1])) {
			++interval;
			continue;
		}
		const double middle=(arc_samples[interval]+arc_samples[interval+1])/2.0;
		if(middle<=arc_samples[interval]+tolerance||middle>=arc_samples[interval+1]-tolerance)
			throw std::runtime_error(context+": adaptive interval refinement reached numerical resolution");
		arc_samples.insert(arc_samples.begin()+interval+1,middle);
		if(arc_samples.size()>100000000)
			throw std::runtime_error(context+": adaptive sampling generated too many points");
	}
	while(arc_samples.size()<4) {
		std::size_t best=0;
		double best_length=-1.0;
		for(std::size_t i=0;i+1<arc_samples.size();++i) {
			const double a=std::max(arc_samples[i],usable_begin);
			const double b=std::min(arc_samples[i+1],usable_end);
			if(b-a>best_length) {best_length=b-a;best=i;}
		}
		if(best_length<=tolerance)
			throw std::runtime_error(context+": branch cannot provide four clearance-respecting samples");
		const double middle=(std::max(arc_samples[best],usable_begin)
			+std::min(arc_samples[best+1],usable_end))/2.0;
		arc_samples.insert(arc_samples.begin()+best+1,middle);
	}

	std::vector<CurveSample> result;
	result.reserve(arc_samples.size());
	for(double s:arc_samples) result.push_back(at_arc(s));
	return result;
}

std::vector<CurveSample> SampleBranch(
	const std::vector<Vec3>& points,
	const std::vector<double>& diameters,
	double segment_length,
	int mode)
{
	BranchSamplingOptions options;
	options.target_spacing=segment_length;
	return SampleBranch(points,diameters,options,mode,"branch");
}

} // namespace tubular

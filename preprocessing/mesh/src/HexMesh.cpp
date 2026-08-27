#include "HexMesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace tubular {
namespace {

QualityResult EvaluateOne(const std::vector<Vec3>& points, const std::vector<Hex>& elements, const std::vector<int>& indices)
{
	static const std::array<double,6> samples{
		-1.0, -0.861136311594053, -0.339981043584856,
		0.339981043584856, 0.861136311594053, 1.0};
	static const std::array<std::array<int,3>,8> signs{{
		{{-1,-1,-1}}, {{1,-1,-1}}, {{1,1,-1}}, {{-1,1,-1}},
		{{-1,-1,1}}, {{1,-1,1}}, {{1,1,1}}, {{-1,1,1}}
	}};
	QualityResult result;
	for (int element_index : indices) {
		if (element_index < 0 || static_cast<std::size_t>(element_index) >= elements.size())
			throw std::runtime_error("hex quality subset contains an invalid element index");
		const auto& element = elements[element_index];
		for (int node : element)
			if (node < 0 || static_cast<std::size_t>(node) >= points.size())
				throw std::runtime_error("hex element contains an invalid point index");
		bool bad = false;
		for (double r : samples) for (double s : samples) for (double t : samples) {
			std::array<Vec3,3> jacobian{};
			for (int a=0; a<8; ++a) {
				const double dr = 0.125*signs[a][0]*(1.0+signs[a][1]*s)*(1.0+signs[a][2]*t);
				const double ds = 0.125*signs[a][1]*(1.0+signs[a][0]*r)*(1.0+signs[a][2]*t);
				const double dt = 0.125*signs[a][2]*(1.0+signs[a][0]*r)*(1.0+signs[a][1]*s);
				const Vec3& p = points[element[a]];
				jacobian[0] += p*dr;
				jacobian[1] += p*ds;
				jacobian[2] += p*dt;
			}
			const double determinant = Determinant(jacobian);
			const double scale = Norm(jacobian[0])*Norm(jacobian[1])*Norm(jacobian[2]);
			const double scaled = std::isfinite(scale) && scale > 0.0
				? determinant/scale : -std::numeric_limits<double>::infinity();
			result.minimum_determinant = std::min(result.minimum_determinant, determinant);
			result.minimum_scaled_jacobian = std::min(result.minimum_scaled_jacobian, scaled);
			if (!std::isfinite(determinant) || determinant <= 0.0) bad = true;
		}
		if (bad) {
			++result.bad_elements;
			if (result.first_bad_element < 0) result.first_bad_element = element_index;
		}
	}
	return result;
}

using Quad = std::array<int,4>;
using Triangle = std::array<int,3>;

struct BoundaryFace
{
	Quad nodes;
	int element = -1;
};

struct SurfaceTriangle
{
	Triangle nodes;
	int face = -1;
	Vec3 lower;
	Vec3 upper;
};

std::vector<BoundaryFace> BoundaryFaces(const std::vector<Hex>& elements)
{
	static const std::array<std::array<int,4>,6> local{{
		{{0,1,2,3}},{{4,5,6,7}},{{0,1,5,4}},
		{{1,2,6,5}},{{2,3,7,6}},{{3,0,4,7}}
	}};
	struct Entry {int count=0;BoundaryFace face;};
	std::map<Quad,Entry> entries;
	for(std::size_t element=0;element<elements.size();++element) for(const auto& indices:local) {
		Quad face{{elements[element][indices[0]],elements[element][indices[1]],
			elements[element][indices[2]],elements[element][indices[3]]}};
		Quad key=face;std::sort(key.begin(),key.end());
		auto& entry=entries[key];
		++entry.count;entry.face={face,static_cast<int>(element)};
	}
	std::vector<BoundaryFace> result;
	for(const auto& item:entries) {
		if(item.second.count==1) result.push_back(item.second.face);
		else if(item.second.count!=2)
			throw std::runtime_error("non-manifold control mesh face detected");
	}
	return result;
}

bool ShareNode(const Triangle& first,const Triangle& second)
{
	for(int a:first) for(int b:second) if(a==b) return true;
	return false;
}

bool SegmentTriangleIntersection(
	const Vec3& p,const Vec3& q,const Vec3& a,const Vec3& b,const Vec3& c,double epsilon)
{
	const Vec3 direction=q-p,edge1=b-a,edge2=c-a;
	const Vec3 h=Cross(direction,edge2);
	const double determinant=Dot(edge1,h);
	if(std::abs(determinant)<=epsilon) return false;
	const double inverse=1.0/determinant;
	const Vec3 s=p-a;
	const double u=inverse*Dot(s,h);
	if(u<-epsilon||u>1.0+epsilon) return false;
	const double v=inverse*Dot(direction,Cross(s,edge1));
	if(v<-epsilon||u+v>1.0+epsilon) return false;
	const double t=inverse*Dot(edge2,Cross(s,edge1));
	return t>=-epsilon&&t<=1.0+epsilon;
}

struct Vec2 {double x=0.0,y=0.0;};

Vec2 Project2(const Vec3& point,int drop)
{
	if(drop==0) return {point.y,point.z};
	if(drop==1) return {point.x,point.z};
	return {point.x,point.y};
}

double Orient2(const Vec2& a,const Vec2& b,const Vec2& c)
{
	return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
}

bool OnSegment2(const Vec2& a,const Vec2& b,const Vec2& p,double epsilon)
{
	return std::abs(Orient2(a,b,p))<=epsilon
		&&p.x>=std::min(a.x,b.x)-epsilon&&p.x<=std::max(a.x,b.x)+epsilon
		&&p.y>=std::min(a.y,b.y)-epsilon&&p.y<=std::max(a.y,b.y)+epsilon;
}

bool SegmentIntersection2(const Vec2& a,const Vec2& b,const Vec2& c,const Vec2& d,double epsilon)
{
	const double o1=Orient2(a,b,c),o2=Orient2(a,b,d),o3=Orient2(c,d,a),o4=Orient2(c,d,b);
	if(((o1>epsilon&&o2<-epsilon)||(o1<-epsilon&&o2>epsilon))
		&&((o3>epsilon&&o4<-epsilon)||(o3<-epsilon&&o4>epsilon))) return true;
	return OnSegment2(a,b,c,epsilon)||OnSegment2(a,b,d,epsilon)
		||OnSegment2(c,d,a,epsilon)||OnSegment2(c,d,b,epsilon);
}

bool PointInTriangle2(const Vec2& p,const std::array<Vec2,3>& triangle,double epsilon)
{
	const double a=Orient2(triangle[0],triangle[1],p);
	const double b=Orient2(triangle[1],triangle[2],p);
	const double c=Orient2(triangle[2],triangle[0],p);
	return (a>=-epsilon&&b>=-epsilon&&c>=-epsilon)
		||(a<=epsilon&&b<=epsilon&&c<=epsilon);
}

bool CoplanarTriangleIntersection(
	const std::array<Vec3,3>& first,const std::array<Vec3,3>& second,const Vec3& normal,double epsilon)
{
	const std::array<double,3> magnitude{{std::abs(normal.x),std::abs(normal.y),std::abs(normal.z)}};
	const int drop=static_cast<int>(std::max_element(magnitude.begin(),magnitude.end())-magnitude.begin());
	std::array<Vec2,3> a,b;
	for(int i=0;i<3;++i) {a[i]=Project2(first[i],drop);b[i]=Project2(second[i],drop);}
	for(int i=0;i<3;++i) for(int j=0;j<3;++j)
		if(SegmentIntersection2(a[i],a[(i+1)%3],b[j],b[(j+1)%3],epsilon)) return true;
	return PointInTriangle2(a[0],b,epsilon)||PointInTriangle2(b[0],a,epsilon);
}

bool TriangleIntersection(
	const std::vector<Vec3>& points,const Triangle& first_nodes,const Triangle& second_nodes)
{
	std::array<Vec3,3> first{{points[first_nodes[0]],points[first_nodes[1]],points[first_nodes[2]]}};
	std::array<Vec3,3> second{{points[second_nodes[0]],points[second_nodes[1]],points[second_nodes[2]]}};
	const Vec3 origin=first[0];
	double coordinate_scale=0.0;
	for(const auto& point:first) coordinate_scale=std::max(coordinate_scale,Norm(point-origin));
	for(const auto& point:second) coordinate_scale=std::max(coordinate_scale,Norm(point-origin));
	if(!(coordinate_scale>0.0)||!std::isfinite(coordinate_scale)) return true;
	for(auto& point:first) point=(point-origin)/coordinate_scale;
	for(auto& point:second) point=(point-origin)/coordinate_scale;
	const Vec3 first_normal=Cross(first[1]-first[0],first[2]-first[0]);
	const Vec3 second_normal=Cross(second[1]-second[0],second[2]-second[0]);
	const double epsilon=1.0e-11;
	if(Norm(first_normal)<=epsilon||Norm(second_normal)<=epsilon) return true;
	const bool parallel=Norm(Cross(first_normal,second_normal))
		<=epsilon*Norm(first_normal)*Norm(second_normal);
	if(parallel) {
		const double plane_distance=std::abs(Dot(second[0]-first[0],Normalized(first_normal,"surface normal")));
		if(plane_distance>epsilon) return false;
		return CoplanarTriangleIntersection(first,second,first_normal,epsilon);
	}
	for(int i=0;i<3;++i)
		if(SegmentTriangleIntersection(first[i],first[(i+1)%3],second[0],second[1],second[2],epsilon)) return true;
	for(int i=0;i<3;++i)
		if(SegmentTriangleIntersection(second[i],second[(i+1)%3],first[0],first[1],first[2],epsilon)) return true;
	return false;
}

} // namespace

QualityResult EvaluateHexQuality(
	const std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<int>* subset)
{
	std::vector<int> all;
	if (!subset) {
		all.resize(elements.size());
		for (std::size_t i=0; i<elements.size(); ++i) all[i] = static_cast<int>(i);
		subset = &all;
	}
	return EvaluateOne(points, elements, *subset);
}

SurfaceIntersectionResult EvaluateBoundarySelfIntersections(
	const std::vector<Vec3>& points,
	const std::vector<Hex>& elements)
{
	const auto faces=BoundaryFaces(elements);
	std::vector<SurfaceTriangle> triangles;
	triangles.reserve(2*faces.size());
	for(std::size_t face=0;face<faces.size();++face) {
		for(const Triangle nodes:{Triangle{{faces[face].nodes[0],faces[face].nodes[1],faces[face].nodes[2]}},
			Triangle{{faces[face].nodes[0],faces[face].nodes[2],faces[face].nodes[3]}}}) {
			SurfaceTriangle triangle;triangle.nodes=nodes;triangle.face=static_cast<int>(face);
			triangle.lower=triangle.upper=points[nodes[0]];
			for(int node:nodes) {
				const auto& p=points[node];
				triangle.lower.x=std::min(triangle.lower.x,p.x);triangle.upper.x=std::max(triangle.upper.x,p.x);
				triangle.lower.y=std::min(triangle.lower.y,p.y);triangle.upper.y=std::max(triangle.upper.y,p.y);
				triangle.lower.z=std::min(triangle.lower.z,p.z);triangle.upper.z=std::max(triangle.upper.z,p.z);
			}
			triangles.push_back(triangle);
		}
	}
	std::vector<int> order(triangles.size());std::iota(order.begin(),order.end(),0);
	std::sort(order.begin(),order.end(),[&](int a,int b){return triangles[a].lower.x<triangles[b].lower.x;});
	for(std::size_t oi=0;oi<order.size();++oi) {
		const auto& first=triangles[order[oi]];
		for(std::size_t oj=oi+1;oj<order.size();++oj) {
			const auto& second=triangles[order[oj]];
			if(second.lower.x>first.upper.x) break;
			if(first.face==second.face||ShareNode(first.nodes,second.nodes)
				||second.lower.y>first.upper.y||second.upper.y<first.lower.y
				||second.lower.z>first.upper.z||second.upper.z<first.lower.z) continue;
			if(TriangleIntersection(points,first.nodes,second.nodes))
				return {1,first.face,second.face};
		}
	}
	return {};
}

std::vector<std::vector<int>> BuildIncidentElements(
	std::size_t point_count,
	const std::vector<Hex>& elements)
{
	std::vector<std::vector<int>> incident(point_count);
	for (std::size_t e=0; e<elements.size(); ++e)
		for (int point : elements[e]) {
			if (point < 0 || static_cast<std::size_t>(point) >= point_count)
				throw std::runtime_error("cannot build incidence for invalid connectivity");
			incident[point].push_back(static_cast<int>(e));
		}
	for (auto& list : incident) {
		std::sort(list.begin(), list.end());
		list.erase(std::unique(list.begin(), list.end()), list.end());
	}
	return incident;
}

double MovePointSafely(
	std::vector<Vec3>& points,
	const std::vector<Hex>& elements,
	const std::vector<std::vector<int>>& incident_elements,
	int point_index,
	const Vec3& candidate,
	double required_scaled_jacobian)
{
	if (point_index < 0 || static_cast<std::size_t>(point_index) >= points.size())
		throw std::runtime_error("point move index is out of range");
	if (!IsFinite(candidate)) throw std::runtime_error("point move candidate is not finite");
	const auto& incident = incident_elements.at(point_index);
	if (incident.empty()) throw std::runtime_error("point move has no incident volume element");
	const Vec3 original = points[point_index];
	static const std::array<double,8> alphas{{1.0,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0}};
	for (double alpha : alphas) {
		points[point_index] = original+(candidate-original)*alpha;
		const auto quality = EvaluateHexQuality(points, elements, &incident);
		if (quality.bad_elements == 0 && quality.minimum_determinant > 0.0
			&& quality.minimum_scaled_jacobian >= required_scaled_jacobian)
			return alpha;
	}
	points[point_index] = original;
	throw std::runtime_error("no valid position found during point backtracking");
}

} // namespace tubular

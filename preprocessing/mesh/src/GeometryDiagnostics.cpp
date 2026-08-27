#include "GeometryDiagnostics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace tubular {
namespace {

struct CenterlineSegment
{
	int parent = -1;
	int child = -1;
	Vec3 a;
	Vec3 b;
	double ra = 0.0;
	double rb = 0.0;
	double length = 0.0;
	Vec3 lower;
	Vec3 upper;
};

double Degrees(double radians)
{
	return radians*180.0/std::acos(-1.0);
}

double VectorAngle(const Vec3& a, const Vec3& b)
{
	return std::acos(ClampUnit(Dot(
		Normalized(a,"diagnostic angle a"),Normalized(b,"diagnostic angle b"))));
}

double NodeDistance(const SwcGraph& graph, int first, int second)
{
	double distance = 0.0;
	auto depth = [&](int node) {
		int result = 0;
		while(graph.nodes[node].parent>=0) {++result;node=graph.nodes[node].parent;}
		return result;
	};
	int first_depth=depth(first),second_depth=depth(second);
	while(first_depth>second_depth) {
		const int parent=graph.nodes[first].parent;
		distance+=Norm(graph.nodes[first].position-graph.nodes[parent].position);
		first=parent;--first_depth;
	}
	while(second_depth>first_depth) {
		const int parent=graph.nodes[second].parent;
		distance+=Norm(graph.nodes[second].position-graph.nodes[parent].position);
		second=parent;--second_depth;
	}
	while(first!=second) {
		const int first_parent=graph.nodes[first].parent;
		const int second_parent=graph.nodes[second].parent;
		distance+=Norm(graph.nodes[first].position-graph.nodes[first_parent].position);
		distance+=Norm(graph.nodes[second].position-graph.nodes[second_parent].position);
		first=first_parent;second=second_parent;
	}
	return distance;
}

double SegmentDistance(
	const Vec3& p1,const Vec3& q1,const Vec3& p2,const Vec3& q2,
	double& first_fraction,double& second_fraction)
{
	const Vec3 d1=q1-p1,d2=q2-p2,r=p1-p2;
	const double a=NormSquared(d1),e=NormSquared(d2),f=Dot(d2,r);
	const double epsilon=1.0e-14;
	if(a<=epsilon&&e<=epsilon) {
		first_fraction=second_fraction=0.0;
		return Norm(p1-p2);
	}
	if(a<=epsilon) {
		first_fraction=0.0;
		second_fraction=std::max(0.0,std::min(1.0,f/e));
	} else {
		const double c=Dot(d1,r);
		if(e<=epsilon) {
			second_fraction=0.0;
			first_fraction=std::max(0.0,std::min(1.0,-c/a));
		} else {
			const double b=Dot(d1,d2),denominator=a*e-b*b;
			first_fraction=denominator!=0.0
				? std::max(0.0,std::min(1.0,(b*f-c*e)/denominator)):0.0;
			second_fraction=(b*first_fraction+f)/e;
			if(second_fraction<0.0) {
				second_fraction=0.0;
				first_fraction=std::max(0.0,std::min(1.0,-c/a));
			} else if(second_fraction>1.0) {
				second_fraction=1.0;
				first_fraction=std::max(0.0,std::min(1.0,(b-c)/a));
			}
		}
	}
	return Norm((p1+d1*first_fraction)-(p2+d2*second_fraction));
}

double PathDistanceBetweenSegmentPoints(
	const SwcGraph& graph,
	const CenterlineSegment& first,double first_fraction,
	const CenterlineSegment& second,double second_fraction)
{
	const std::array<int,2> first_nodes{{first.parent,first.child}};
	const std::array<int,2> second_nodes{{second.parent,second.child}};
	const std::array<double,2> first_offsets{{first_fraction*first.length,
		(1.0-first_fraction)*first.length}};
	const std::array<double,2> second_offsets{{second_fraction*second.length,
		(1.0-second_fraction)*second.length}};
	double result=std::numeric_limits<double>::infinity();
	for(int i=0;i<2;++i) for(int j=0;j<2;++j)
		result=std::min(result,first_offsets[i]+NodeDistance(graph,first_nodes[i],second_nodes[j])
			+second_offsets[j]);
	return result;
}

std::string NodeName(const SwcGraph& graph,int node)
{
	return std::to_string(graph.nodes[node].id>0?graph.nodes[node].id:node+1);
}

void Mark(std::vector<int>& risk,int node,int value)
{
	if(node>=0) risk.at(static_cast<std::size_t>(node))
		=std::max(risk.at(static_cast<std::size_t>(node)),value);
}

std::string JsonEscape(const std::string& value)
{
	std::string result;
	for(char c:value) {
		if(c=='\\'||c=='"') {result.push_back('\\');result.push_back(c);}
		else if(c=='\n') result+="\\n";
		else result.push_back(c);
	}
	return result;
}

} // namespace

GeometryDiagnostics AnalyzeSkeletonGeometry(
	const SwcGraph& graph,
	const MeshParameters& parameters)
{
	graph.Validate();
	parameters.Validate();
	GeometryDiagnostics result;
	const std::size_t count=graph.nodes.size();
	result.node_min_length_over_diameter.assign(count,std::numeric_limits<double>::infinity());
	result.node_curvature_radius_product.assign(count,0.0);
	result.node_junction_min_angle_degrees.assign(count,0.0);
	result.node_risk.assign(count,0);
	std::vector<CenterlineSegment> centerline;
	centerline.reserve(count-1);

	for(std::size_t child=0;child<count;++child) {
		const int parent=graph.nodes[child].parent;
		if(parent<0) continue;
		const double length=Norm(graph.nodes[child].position-graph.nodes[parent].position);
		const double minimum_diameter=std::min(graph.nodes[child].diameter,graph.nodes[parent].diameter);
		const double maximum_diameter=std::max(graph.nodes[child].diameter,graph.nodes[parent].diameter);
		const double ratio=length/maximum_diameter;
		const double change=std::abs(graph.nodes[child].diameter-graph.nodes[parent].diameter)
			/minimum_diameter;
		result.segments.push_back({parent,static_cast<int>(child),length,ratio,change});
		result.node_min_length_over_diameter[parent]
			=std::min(result.node_min_length_over_diameter[parent],ratio);
		result.node_min_length_over_diameter[child]
			=std::min(result.node_min_length_over_diameter[child],ratio);
		if(change>parameters.max_diameter_change_fraction*(1.0+1.0e-4)
			&& !graph.is_branch(parent) && !graph.is_branch(static_cast<int>(child))) {
			std::ostringstream message;
			message<<"segment "<<NodeName(graph,parent)<<"->"<<NodeName(graph,static_cast<int>(child))
				<<" changes diameter by "<<change<<", above adaptive limit "
				<<parameters.max_diameter_change_fraction;
			result.errors.push_back(message.str());
			Mark(result.node_risk,parent,2);Mark(result.node_risk,static_cast<int>(child),2);
		}
		CenterlineSegment segment;
		segment.parent=parent;segment.child=static_cast<int>(child);
		segment.a=graph.nodes[parent].position;segment.b=graph.nodes[child].position;
		segment.ra=graph.nodes[parent].diameter/2.0;segment.rb=graph.nodes[child].diameter/2.0;
		segment.length=length;
		const double expansion=parameters.collision_safety_factor*std::max(segment.ra,segment.rb);
		segment.lower={std::min(segment.a.x,segment.b.x)-expansion,
			std::min(segment.a.y,segment.b.y)-expansion,
			std::min(segment.a.z,segment.b.z)-expansion};
		segment.upper={std::max(segment.a.x,segment.b.x)+expansion,
			std::max(segment.a.y,segment.b.y)+expansion,
			std::max(segment.a.z,segment.b.z)+expansion};
		centerline.push_back(segment);
	}

	for(std::size_t node=0;node<count;++node) {
		const int parent=graph.nodes[node].parent;
		if(parent>=0&&graph.nodes[node].children.size()==1) {
			const int child=graph.nodes[node].children.front();
			const Vec3 incoming=graph.nodes[node].position-graph.nodes[parent].position;
			const Vec3 outgoing=graph.nodes[child].position-graph.nodes[node].position;
			const double turn=VectorAngle(incoming,outgoing);
			const double scale=std::min(Norm(incoming),Norm(outgoing));
			const double product=2.0*std::sin(turn/2.0)/scale*graph.nodes[node].diameter/2.0;
			result.node_curvature_radius_product[node]=product;
			if(product>=parameters.maximum_curvature_radius_product) {
				std::ostringstream message;
				message<<"node "<<NodeName(graph,static_cast<int>(node))
					<<" has curvature*radius="<<product<<", limit="
					<<parameters.maximum_curvature_radius_product;
				result.errors.push_back(message.str());Mark(result.node_risk,static_cast<int>(node),2);
			}
		}
		if(!graph.is_branch(static_cast<int>(node))) continue;
		const int branch=static_cast<int>(node);
		const int branch_parent=graph.nodes[node].parent;
		if(branch_parent<0) continue;
		const int child0=graph.nodes[node].children[0],child1=graph.nodes[node].children[1];
		const std::array<Vec3,3> directions{{
			graph.nodes[branch_parent].position-graph.nodes[node].position,
			graph.nodes[child0].position-graph.nodes[node].position,
			graph.nodes[child1].position-graph.nodes[node].position}};
		const double minimum_angle=Degrees(std::min({VectorAngle(directions[0],directions[1]),
			VectorAngle(directions[0],directions[2]),VectorAngle(directions[1],directions[2])}));
		const std::array<double,3> radii{{graph.nodes[branch_parent].diameter/2.0,
			graph.nodes[child0].diameter/2.0,graph.nodes[child1].diameter/2.0}};
		const double radius_ratio=*std::max_element(radii.begin(),radii.end())
			/ *std::min_element(radii.begin(),radii.end());
		const double upstream=Norm(directions[0])
			/std::max(graph.nodes[node].diameter,graph.nodes[branch_parent].diameter);
		const double downstream0=Norm(directions[1])
			/std::max(graph.nodes[node].diameter,graph.nodes[child0].diameter);
		const double downstream1=Norm(directions[2])
			/std::max(graph.nodes[node].diameter,graph.nodes[child1].diameter);
		const double downstream=std::min(downstream0,downstream1);
		result.junctions.push_back({branch,minimum_angle,radius_ratio,upstream,downstream});
		result.node_junction_min_angle_degrees[node]=minimum_angle;
		auto junction_error=[&](bool invalid,const std::string& detail) {
			if(!invalid) return;
			result.errors.push_back("bifurcation node "+NodeName(graph,branch)+" "+detail);
			Mark(result.node_risk,branch,2);Mark(result.node_risk,branch_parent,2);
			Mark(result.node_risk,child0,2);Mark(result.node_risk,child1,2);
		};
		junction_error(minimum_angle+1.0e-8<parameters.minimum_bifurcation_angle_degrees,
			"minimum angle "+std::to_string(minimum_angle)+" degrees is below limit "
			+std::to_string(parameters.minimum_bifurcation_angle_degrees));
		junction_error(radius_ratio>parameters.maximum_junction_radius_ratio*(1.0+1.0e-8),
			"radius ratio "+std::to_string(radius_ratio)+" exceeds limit "
			+std::to_string(parameters.maximum_junction_radius_ratio));
		junction_error(upstream+1.0e-6<0.98*parameters.upstream_clearance_over_diameter,
			"upstream clearance/diameter "+std::to_string(upstream)+" is below required "
			+std::to_string(parameters.upstream_clearance_over_diameter));
		junction_error(downstream+1.0e-6<0.98*parameters.downstream_clearance_over_diameter,
			"downstream clearance/diameter "+std::to_string(downstream)+" is below required "
			+std::to_string(parameters.downstream_clearance_over_diameter));
	}

	if(parameters.check_self_intersection) {
		std::vector<int> order(centerline.size());
		std::iota(order.begin(),order.end(),0);
		std::sort(order.begin(),order.end(),[&](int a,int b) {
			return centerline[a].lower.x<centerline[b].lower.x;
		});
		for(std::size_t oi=0;oi<order.size();++oi) {
			const auto& first=centerline[order[oi]];
			for(std::size_t oj=oi+1;oj<order.size();++oj) {
				const auto& second=centerline[order[oj]];
				if(second.lower.x>first.upper.x) break;
				if(second.lower.y>first.upper.y||second.upper.y<first.lower.y
					||second.lower.z>first.upper.z||second.upper.z<first.lower.z) continue;
				if(first.parent==second.parent||first.parent==second.child
					||first.child==second.parent||first.child==second.child) continue;
				double sf=0.0,tf=0.0;
				const double distance=SegmentDistance(first.a,first.b,second.a,second.b,sf,tf);
				const double first_radius=first.ra+(first.rb-first.ra)*sf;
				const double second_radius=second.ra+(second.rb-second.ra)*tf;
				const double required=parameters.collision_safety_factor*(first_radius+second_radius);
				if(distance>=required) continue;
				const double path_distance=PathDistanceBetweenSegmentPoints(graph,first,sf,second,tf);
				if(path_distance<=1.05*required) continue;
				result.collisions.push_back({first.parent,first.child,second.parent,second.child,
					distance,required});
				std::ostringstream message;
				message<<"centerline tubes "<<NodeName(graph,first.parent)<<"->"<<NodeName(graph,first.child)
					<<" and "<<NodeName(graph,second.parent)<<"->"<<NodeName(graph,second.child)
					<<" overlap: distance="<<distance<<" required="<<required;
				result.warnings.push_back(message.str());
				for(int node:{first.parent,first.child,second.parent,second.child}) Mark(result.node_risk,node,1);
			}
		}
	}

	for(double& value:result.node_min_length_over_diameter)
		if(!std::isfinite(value)) value=0.0;
	return result;
}

void WriteGeometryDiagnosticsJson(
	const GeometryDiagnostics& diagnostics,
	const SwcGraph& graph,
	const MeshParameters& parameters,
	const std::filesystem::path& path)
{
	std::ofstream output(path);
	if(!output) throw std::runtime_error("cannot write geometry diagnostics: "+path.string());
	output<<std::setprecision(17)<<"{\n  \"schema_version\": 1,\n  \"valid\": "
		<<(diagnostics.valid()?"true":"false")<<",\n  \"nodes\": "<<graph.nodes.size()
		<<",\n  \"effective_limits\": {"
		<<"\"max_spacing_over_diameter\": "<<parameters.max_spacing_over_diameter
		<<", \"max_turn_degrees\": "<<parameters.max_turn_degrees
		<<", \"max_diameter_change_fraction\": "<<parameters.max_diameter_change_fraction
		<<", \"maximum_curvature_radius_product\": "<<parameters.maximum_curvature_radius_product
		<<", \"minimum_bifurcation_angle_degrees\": "<<parameters.minimum_bifurcation_angle_degrees
		<<", \"maximum_junction_radius_ratio\": "<<parameters.maximum_junction_radius_ratio
		<<", \"minimum_scaled_jacobian\": "<<parameters.minimum_scaled_jacobian
		<<", \"collision_safety_factor\": "<<parameters.collision_safety_factor
		<<"},\n  \"segments\": [\n";
	for(std::size_t i=0;i<diagnostics.segments.size();++i) {
		const auto& row=diagnostics.segments[i];
		output<<"    {\"parent_id\": "<<graph.nodes[row.parent].id
			<<", \"child_id\": "<<graph.nodes[row.child].id<<", \"length\": "<<row.length
			<<", \"length_over_diameter\": "<<row.length_over_diameter
			<<", \"diameter_change_fraction\": "<<row.diameter_change_fraction<<"}"
			<<(i+1==diagnostics.segments.size()?"\n":",\n");
	}
	output<<"  ],\n  \"junctions\": [\n";
	for(std::size_t i=0;i<diagnostics.junctions.size();++i) {
		const auto& row=diagnostics.junctions[i];
		output<<"    {\"node_id\": "<<graph.nodes[row.node].id
			<<", \"minimum_angle_degrees\": "<<row.minimum_angle_degrees
			<<", \"radius_ratio\": "<<row.radius_ratio
			<<", \"upstream_clearance_over_diameter\": "<<row.upstream_clearance_over_diameter
			<<", \"minimum_downstream_clearance_over_diameter\": "
			<<row.minimum_downstream_clearance_over_diameter<<"}"
			<<(i+1==diagnostics.junctions.size()?"\n":",\n");
	}
	output<<"  ],\n  \"collisions\": [\n";
	for(std::size_t i=0;i<diagnostics.collisions.size();++i) {
		const auto& row=diagnostics.collisions[i];
		output<<"    {\"first\": ["<<graph.nodes[row.first_parent].id<<","<<graph.nodes[row.first_child].id
			<<"], \"second\": ["<<graph.nodes[row.second_parent].id<<","<<graph.nodes[row.second_child].id
			<<"], \"distance\": "<<row.distance<<", \"required_distance\": "<<row.required_distance<<"}"
			<<(i+1==diagnostics.collisions.size()?"\n":",\n");
	}
	output<<"  ],\n";
	auto write_messages=[&](const char* name,const std::vector<std::string>& messages) {
		output<<"  \""<<name<<"\": [";
		for(std::size_t i=0;i<messages.size();++i)
			output<<(i?", ":"")<<"\""<<JsonEscape(messages[i])<<"\"";
		output<<"]";
	};
	write_messages("warnings",diagnostics.warnings);
	output<<",\n";
	write_messages("errors",diagnostics.errors);
	output<<"\n}\n";
}

void WriteGeometryDiagnosticsVtp(
	const GeometryDiagnostics& diagnostics,
	const SwcGraph& graph,
	const std::filesystem::path& path)
{
	std::ofstream output(path);
	if(!output) throw std::runtime_error("cannot write diagnostic VTP: "+path.string());
	output<<std::setprecision(17)
		<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		<<"<PolyData><Piece NumberOfPoints=\""<<graph.nodes.size()<<"\" NumberOfLines=\""
		<<diagnostics.segments.size()<<"\">\n<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for(const auto& node:graph.nodes) output<<node.position.x<<' '<<node.position.y<<' '<<node.position.z<<'\n';
	output<<"</DataArray></Points>\n<Lines><DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for(const auto& row:diagnostics.segments) output<<row.parent<<' '<<row.child<<'\n';
	output<<"</DataArray><DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for(std::size_t i=0;i<diagnostics.segments.size();++i) output<<2*(i+1)<<'\n';
	output<<"</DataArray></Lines>\n<PointData>\n";
	auto point_array=[&](const char* name,const auto& values,const char* type) {
		output<<"<DataArray type=\""<<type<<"\" Name=\""<<name<<"\" format=\"ascii\">\n";
		for(const auto& value:values) output<<value<<'\n';
		output<<"</DataArray>\n";
	};
	point_array("min_length_over_diameter",diagnostics.node_min_length_over_diameter,"Float64");
	point_array("curvature_radius_product",diagnostics.node_curvature_radius_product,"Float64");
	point_array("junction_min_angle_degrees",diagnostics.node_junction_min_angle_degrees,"Float64");
	point_array("risk",diagnostics.node_risk,"Int32");
	std::vector<int> ids;ids.reserve(graph.nodes.size());for(const auto& node:graph.nodes) ids.push_back(node.id);
	point_array("node_id",ids,"Int32");
	output<<"</PointData>\n<CellData>\n";
	output<<"<DataArray type=\"Float64\" Name=\"length_over_diameter\" format=\"ascii\">\n";
	for(const auto& row:diagnostics.segments) output<<row.length_over_diameter<<'\n';
	output<<"</DataArray>\n<DataArray type=\"Float64\" Name=\"diameter_change_fraction\" format=\"ascii\">\n";
	for(const auto& row:diagnostics.segments) output<<row.diameter_change_fraction<<'\n';
	output<<"</DataArray>\n</CellData>\n</Piece></PolyData></VTKFile>\n";
}

void RequireValidGeometry(const GeometryDiagnostics& diagnostics)
{
	if(diagnostics.valid()) return;
	throw std::runtime_error("geometry preflight failed with "+std::to_string(diagnostics.errors.size())
		+" error(s): "+diagnostics.errors.front());
}

} // namespace tubular

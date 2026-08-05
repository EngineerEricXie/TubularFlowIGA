#include "MeshGenerator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace tubular {
namespace {

using Face = std::array<int,4>;
using FaceList = std::vector<Face>;

struct Templates
{
	std::vector<Vec3> circle;
	std::vector<Vec3> merge;
	FaceList circle_faces;
	FaceList merge_faces;
	FaceList branch_bottom;
	FaceList branch_left;
	FaceList branch_right;
	std::vector<Vec3> branch_bottom_points;
	std::vector<Vec3> branch_left_points;
	std::vector<Vec3> branch_right_points;
	std::vector<int> boundary_circle;
	std::vector<int> boundary_merge;
	std::vector<int> bottom_boundary;
	std::vector<int> left_boundary;
	std::vector<int> right_boundary;
};

struct Layer
{
	std::vector<Vec3> points;
	FaceList primary;
	FaceList alternate;
	FaceList right;
	Vec3 reference;
	bool has_reference = false;
};

struct Bifurcation
{
	int node = -1;
	int parent = -1;
	int child1 = -1;
	int child2 = -1;
	std::array<int,4> offsets{{-1,-1,-1,-1}};
};

std::vector<Vec3> ReadPoints(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open template points: "+path.string());
	std::vector<Vec3> result;
	Vec3 point;
	while (input >> point.x >> point.y >> point.z) result.push_back(point);
	if (result.empty()) throw std::runtime_error("template point file is empty: "+path.string());
	return result;
}

FaceList ReadFaces(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open template elements: "+path.string());
	FaceList result;
	int count;
	Face face;
	while (input >> count >> face[0] >> face[1] >> face[2] >> face[3]) {
		if (count != 4) throw std::runtime_error("template surface element is not a quad");
		result.push_back(face);
	}
	if (result.empty()) throw std::runtime_error("template element file is empty: "+path.string());
	return result;
}

Vec3 RotateRowX(const Vec3& p, double degrees)
{
	const double angle = degrees*std::acos(-1.0)/180.0;
	const double c = std::cos(angle);
	const double s = std::sin(angle);
	return {p.x, p.y*c+p.z*s, -p.y*s+p.z*c};
}

Templates ReadTemplates(const std::filesystem::path& directory)
{
	Templates t;
	t.circle = ReadPoints(directory/"template_circle_points90.txt");
	t.merge = ReadPoints(directory/"template_merge120_points90.txt");
	t.circle_faces = ReadFaces(directory/"template_circle_elements.txt");
	t.merge_faces = ReadFaces(directory/"template_merge_elements.txt");
	if (t.circle.size() != 201 || t.circle_faces.size() != 180
		|| t.merge.size() != 294 || t.merge_faces.size() != 270)
		throw std::runtime_error("unexpected mesh template dimensions");

	for (std::size_t i=0; i<t.circle.size(); ++i)
		if (Norm(t.circle[i]) > 0.95) t.boundary_circle.push_back(static_cast<int>(i));
	for (std::size_t i=0; i<t.merge.size(); ++i) {
		if (Norm(t.merge[i]) > 0.95) {
			t.boundary_merge.push_back(static_cast<int>(i));
			if (t.merge[i].y > 0.0) t.left_boundary.push_back(static_cast<int>(i));
			else if (t.merge[i].y < 0.0) t.right_boundary.push_back(static_cast<int>(i));
			if (t.merge[i].z > 0.0) t.bottom_boundary.push_back(static_cast<int>(i));
		}
	}
	if (t.bottom_boundary.size() != t.left_boundary.size()
		|| t.left_boundary.size() != t.right_boundary.size())
		throw std::runtime_error("bifurcation boundary template lists have inconsistent sizes");

	FaceList reordered = t.merge_faces;
	for (std::size_t i=2*reordered.size()/3; i<reordered.size(); ++i)
		reordered[i] = {reordered[i][0], reordered[i][3], reordered[i][2], reordered[i][1]};
	std::vector<int> a{90,91,92,93,95,96,97,100,101,105};
	for (int value=125; value<=134; ++value) a.push_back(value);
	std::vector<int> b{114,109,104,99,113,108,103,112,107,111,124,122,120,118,116,123,121,119,117,115};
	if (a.size() != b.size()) throw std::runtime_error("internal bifurcation reorder table mismatch");
	for (std::size_t i=0; i<a.size(); ++i) {
		const int ai = a[i]+1+static_cast<int>(t.circle_faces.size()/2)-1;
		const int bi = b[i]+1+static_cast<int>(t.circle_faces.size()/2)-1;
		std::swap(reordered.at(ai), reordered.at(bi));
		std::swap(reordered.at(ai+static_cast<int>(t.circle_faces.size()/4)),
			reordered.at(bi+static_cast<int>(t.circle_faces.size()/4)));
	}
	for (int i=1; i<=static_cast<int>(t.circle_faces.size()/4); ++i) {
		const int first = i+static_cast<int>(t.circle_faces.size())-1;
		const int second = i+static_cast<int>(5*t.circle_faces.size()/4)-1;
		std::swap(reordered.at(first), reordered.at(second));
	}

	t.branch_bottom = t.circle_faces;
	t.branch_left.insert(t.branch_left.end(), t.circle_faces.begin(), t.circle_faces.begin()+t.circle_faces.size()/2);
	t.branch_left.insert(t.branch_left.end(), t.merge_faces.begin()+t.circle_faces.size(), t.merge_faces.end());
	t.branch_right.insert(t.branch_right.end(), reordered.begin()+t.circle_faces.size(), reordered.end());
	t.branch_right.insert(t.branch_right.end(), t.circle_faces.begin()+t.circle_faces.size()/2, t.circle_faces.end());
	if (t.branch_left.size() != t.circle_faces.size() || t.branch_right.size() != t.circle_faces.size())
		throw std::runtime_error("branch template face count mismatch");

	t.branch_bottom_points.assign(t.merge.begin(), t.merge.begin()+t.circle.size());
	t.branch_left_points = t.branch_bottom_points;
	t.branch_right_points = t.branch_bottom_points;
	for (std::size_t i=0; i<t.circle.size(); ++i) {
		if (t.merge[i].y < 0.0) t.branch_left_points[i] = RotateRowX(t.merge[i], 120.0);
		if (t.merge[i].y > 0.0) t.branch_right_points[i] = RotateRowX(t.merge[i], -120.0);
	}
	return t;
}

std::vector<Vec3> Scale(const std::vector<Vec3>& input, double scale)
{
	std::vector<Vec3> output = input;
	for (auto& point : output) point *= scale;
	return output;
}

void Translate(std::vector<Vec3>& points, const Vec3& shift)
{
	for (auto& point : points) point += shift;
}

FaceList OffsetFaces(const FaceList& faces, int offset)
{
	FaceList result = faces;
	for (auto& face : result) for (int& point : face) point += offset;
	return result;
}

int AddPoints(
	ControlMesh& mesh,
	const std::vector<Vec3>& points,
	int interior_label,
	const std::vector<int>& boundary,
	const Vec3& velocity)
{
	const int offset = static_cast<int>(mesh.points.size());
	mesh.points.insert(mesh.points.end(), points.begin(), points.end());
	mesh.labels.insert(mesh.labels.end(), points.size(), interior_label);
	mesh.velocity.insert(mesh.velocity.end(), points.size(), velocity);
	for (int index : boundary) mesh.labels.at(offset+index) = 0;
	return offset;
}

Vec3 VelocityAt(const Vec3& direction, const Vec3& template_point)
{
	return Normalized(direction, "velocity direction")*std::abs(1.0-NormSquared(template_point));
}

std::vector<Vec3> MakeCircle(
	const Templates& t,
	double radius,
	const Vec3& first_axis,
	const Vec3& second_axis,
	const Vec3& center)
{
	std::vector<Vec3> result(t.circle.size());
	for (std::size_t i=0; i<t.circle.size(); ++i)
		result[i] = center+(first_axis*t.circle[i].x+second_axis*t.circle[i].y)*radius;
	return result;
}

void ConnectFaces(ControlMesh& mesh, const FaceList& first, const FaceList& second)
{
	if (first.size() != second.size()) throw std::runtime_error("cannot connect layers with different face counts");
	for (std::size_t i=0; i<first.size(); ++i)
		mesh.elements.push_back({first[i][0],first[i][1],first[i][2],first[i][3],
			second[i][0],second[i][1],second[i][2],second[i][3]});
}

double AngleBetween(const Vec3& a, const Vec3& b)
{
	return std::acos(ClampUnit(Dot(Normalized(a, "angle vector a"), Normalized(b, "angle vector b"))));
}

Vec3 RotateIfNeeded(const Vec3& point, const Vec3& axis, double angle)
{
	if (std::abs(angle) <= 1.0e-14) return point;
	if (Norm(axis) <= 1.0e-14) throw std::runtime_error("nonzero rotation has a degenerate axis");
	return RotateAroundAxis(point, axis, angle);
}

Face ReorderFace(const Face& f)
{
	return {f[2],f[3],f[0],f[1]};
}

void AdjustEndFaces(FaceList& faces, double& angle, double sign)
{
	const double pi = std::acos(-1.0);
	const std::size_t quarter = faces.size()/4;
	FaceList old = faces;
	if (angle > pi/4.0 && angle <= 3.0*pi/4.0) {
		angle -= pi/2.0;
		for (std::size_t q=0; q<quarter; ++q) {
			if (sign > 0.0) {
				faces[q] = ReorderFace(old[q+3*quarter]);
				faces[q+quarter] = old[q];
				faces[q+2*quarter] = ReorderFace(old[q+quarter]);
				faces[q+3*quarter] = old[q+2*quarter];
			} else if (sign < 0.0) {
				faces[q] = old[q+quarter];
				faces[q+quarter] = ReorderFace(old[q+2*quarter]);
				faces[q+2*quarter] = old[q+3*quarter];
				faces[q+3*quarter] = ReorderFace(old[q]);
			}
		}
	} else if (angle > 3.0*pi/4.0 && angle <= pi+1.0e-12) {
		angle -= pi;
		for (std::size_t q=0; q<quarter; ++q) {
			faces[q] = ReorderFace(old[q+2*quarter]);
			faces[q+quarter] = ReorderFace(old[q+3*quarter]);
			faces[q+2*quarter] = ReorderFace(old[q]);
			faces[q+3*quarter] = ReorderFace(old[q+quarter]);
		}
	}
}

std::vector<int> TrimSection(const SwcGraph& graph, const std::vector<int>& section)
{
	std::size_t begin = graph.is_branch(section.front()) ? 1 : 0;
	std::size_t end = section.size()-(graph.is_branch(section.back()) ? 1 : 0);
	if (begin >= end) throw std::runtime_error("section has no non-branch endpoint layer");
	return std::vector<int>(section.begin()+begin, section.begin()+end);
}

} // namespace

ControlMesh GenerateControlMesh(
	const SwcGraph& skeleton,
	const MeshParameters& parameters,
	const std::filesystem::path& template_directory,
	double minimum_scaled_jacobian)
{
	skeleton.Validate();
	if (!std::isfinite(minimum_scaled_jacobian) || minimum_scaled_jacobian <= 0.0)
		throw std::runtime_error("minimum scaled Jacobian must be positive");
	const Templates t = ReadTemplates(template_directory);
	const int root = skeleton.root();
	const auto sections = skeleton.Sections();
	const std::size_t node_count = skeleton.nodes.size();

	std::vector<Vec3> segment(node_count);
	std::vector<int> layer_count(node_count, 0);
	std::vector<int> child_arm(node_count, 0);
	std::vector<int> branch_nodes;
	for (std::size_t i=0; i<node_count; ++i) {
		if (skeleton.is_branch(static_cast<int>(i))) branch_nodes.push_back(static_cast<int>(i));
		if (skeleton.nodes[i].parent >= 0) {
			const int parent = skeleton.nodes[i].parent;
			segment[i] = skeleton.nodes[i].position-skeleton.nodes[parent].position;
			const double length = Norm(segment[i]);
			const double diameter = std::min(skeleton.nodes[i].diameter, skeleton.nodes[parent].diameter);
			if (!std::isfinite(length) || length <= 0.0 || !std::isfinite(diameter) || diameter <= 0.0)
				throw std::runtime_error("invalid skeleton segment at node "+std::to_string(i+1));
			layer_count[i] = std::max(1, static_cast<int>(std::ceil(length/(1.1*diameter))));
			if (skeleton.is_branch(static_cast<int>(i)) || skeleton.is_branch(parent)) {
				if (parameters.bifurcation_refinement <= 0.0)
					throw std::runtime_error("bifurcation refinement must be positive for branched skeletons");
				const int branch = skeleton.is_branch(static_cast<int>(i)) ? static_cast<int>(i) : parent;
				layer_count[i] = std::max(layer_count[i], static_cast<int>(std::ceil(
					length/(skeleton.nodes[branch].diameter*parameters.bifurcation_refinement))));
			}
		}
	}

	ControlMesh mesh;
	std::vector<Layer> layers(node_count);
	std::vector<std::array<std::vector<Vec3>,3>> branch_points(node_count);
	std::vector<Bifurcation> bifurcations;
	int tip_label = 1;
	const Vec3 root_position = skeleton.nodes[root].position;

	for (int branch : branch_nodes) {
		const int parent = skeleton.nodes[branch].parent;
		if (parent < 0 || skeleton.nodes[parent].parent < 0)
			throw std::runtime_error("a bifurcation must have at least two upstream segments");
		const int child1 = skeleton.nodes[branch].children[0];
		const int child2 = skeleton.nodes[branch].children[1];
		child_arm[child1] = 1;
		child_arm[child2] = 2;
		Bifurcation info;
		info.node=branch; info.parent=parent; info.child1=child1; info.child2=child2;

		const Vec3 sv_parent_parent = segment[parent];
		const Vec3 sv_parent = segment[branch];
		const Vec3 sv_child1 = segment[child1];
		const Vec3 sv_child2 = segment[child2];
		const double ri = skeleton.nodes[parent].diameter/2.0;
		const double rj = skeleton.nodes[child1].diameter/2.0;
		const double rk = skeleton.nodes[child2].diameter/2.0;
		const double average_radius = (ri+rj+rk)/3.0;

		const Vec3 vi = -Normalized(sv_parent, "bifurcation parent direction");
		const Vec3 vj = Normalized(sv_child1, "bifurcation child1 direction");
		const Vec3 vk = Normalized(sv_child2, "bifurcation child2 direction");
		const double aij=AngleBetween(vi,vj), aik=AngleBetween(vi,vk), akj=AngleBetween(vk,vj);
		Vec3 kij=Normalized(ri*vj+rj*vi, "bifurcation ij separator");
		Vec3 kik=Normalized(ri*vk+rk*vi, "bifurcation ik separator");
		Vec3 kkj=Normalized(rk*vj+rj*vk, "bifurcation kj separator");
		if (Dot(Cross(vi,kij),Cross(vi,kik)) > 0.0) {
			if (aij>aik) kij=-kij; else kik=-kik;
		}
		const double pi=std::acos(-1.0);
		Vec3 spij = aij<=pi/2.0 ? kij*ri/std::sin(std::atan(ri/rj)) : kij*(ri+rj)/2.0;
		Vec3 spik = aik<=pi/2.0 ? kik*ri/std::sin(std::atan(ri/rk)) : kik*(ri+rk)/2.0;
		Vec3 spkj = akj<=pi/2.0 ? kkj*rk/std::sin(std::atan(rk/rj)) : kkj*(rk+rj)/2.0;
		const Vec3 cpn=Normalized(Cross(spik-spkj,spij-spkj), "bifurcation plane");
		std::cout << "bifurcation_node=" << branch+1 << " plane_normal=[" << cpn.x << "," << cpn.y << "," << cpn.z << "]\n";
		const Vec3 cp1=cpn*average_radius;
		const Vec3 w_parent=Normalized(Cross(sv_parent,cp1), "bifurcation parent frame");
		const Vec3 w_child1=Normalized(Cross(sv_child1,cp1), "bifurcation child1 frame");
		const Vec3 w_child2=Normalized(Cross(sv_child2,cp1), "bifurcation child2 frame");

		auto parent_terminal=Scale(t.circle,ri);
		auto child1_terminal=Scale(t.circle,rj);
		auto child2_terminal=Scale(t.circle,rk);
		for (std::size_t k=0;k<t.circle.size();++k) {
			parent_terminal[k]=Normalized(cp1,"bifurcation cp1")*parent_terminal[k].x+w_parent*parent_terminal[k].y;
			child1_terminal[k]=Normalized(cp1,"bifurcation cp1")*child1_terminal[k].x+w_child1*child1_terminal[k].y;
			child2_terminal[k]=Normalized(cp1,"bifurcation cp1")*child2_terminal[k].x+w_child2*child2_terminal[k].y;
		}

		Vec3 ciji=w_parent*ri-sv_parent, ciki=-w_parent*ri-sv_parent;
		Vec3 cijj=w_child1*rj+sv_child1, cjkj=-w_child1*rj+sv_child1;
		Vec3 cjkk=w_child2*rk+sv_child2, cikk=-w_child2*rk+sv_child2;
		const int ni=layer_count[branch], nj=layer_count[child1], nk=layer_count[child2];
		ciji=(spij*(ni-1)+ciji)/ni; ciki=(spik*(ni-1)+ciki)/ni;
		cijj=(spij*(nj-1)+cijj)/nj; cjkj=(spkj*(nj-1)+cjkj)/nj;
		cikk=(spik*(nk-1)+cikk)/nk; cjkk=(spkj*(nk-1)+cjkk)/nk;
		auto balance=[](const Vec3& split,const Vec3& a,const Vec3& b) {
			const double da=Norm(split-a), db=Norm(split-b);
			if (da+db<=1.0e-14) throw std::runtime_error("degenerate bifurcation separator balance");
			return (a*db+b*da)/(da+db);
		};
		spij=balance(spij,ciji,cijj); spik=balance(spik,ciki,cikk); spkj=balance(spkj,cjkj,cjkk);

		auto bottom=Scale(t.branch_bottom_points,average_radius);
		auto left=Scale(t.branch_left_points,average_radius);
		auto right=Scale(t.branch_right_points,average_radius);
		for (std::size_t k=0;k<t.circle.size();++k) {
			const Vec3& pb=t.branch_bottom_points[k];
			bottom[k]=cp1*pb.x+(pb.y<0.0?spik*(std::hypot(pb.y,pb.z)):pb.y>0.0?spij*(std::hypot(pb.y,pb.z)):Vec3{});
			const Vec3& pl=t.branch_left_points[k];
			left[k]=cp1*pl.x+(pl.z<0.0?spij*std::hypot(pl.y,pl.z):pl.z>0.0?spkj*pl.z:Vec3{});
			const Vec3& pr=t.branch_right_points[k];
			right[k]=cp1*pr.x+(pr.z<0.0?spik*std::hypot(pr.y,pr.z):pr.z>0.0?spkj*pr.z:Vec3{});
		}
		auto merge=Scale(t.merge,average_radius);
		for (std::size_t k=0;k<t.merge.size();++k) {
			const Vec3& p=t.merge[k];
			if (p.z>0.0) merge[k]=cp1*p.x+spkj*p.z;
			else if (p.y>0.0) merge[k]=cp1*p.x+spij*std::hypot(p.y,p.z);
			else if (p.y<0.0) merge[k]=cp1*p.x+spik*std::hypot(p.y,p.z);
			else merge[k]=cp1*p.x;
		}
		const Vec3 center=skeleton.nodes[branch].position-root_position;
		Translate(merge,center); Translate(bottom,center); Translate(left,center); Translate(right,center);
		Translate(parent_terminal,center-sv_parent);
		Translate(child1_terminal,center+sv_child1);
		Translate(child2_terminal,center+sv_child2);

		const int central_offset=AddPoints(mesh,merge,-1,t.boundary_merge,{0,0,0});
		info.offsets[0]=central_offset;
		layers[branch].points=merge;
		layers[branch].primary=OffsetFaces(t.branch_bottom,central_offset);
		layers[branch].alternate=OffsetFaces(t.branch_left,central_offset);
		layers[branch].right=OffsetFaces(t.branch_right,central_offset);
		branch_points[branch]={bottom,left,right};

		auto add_terminal=[&](int node,std::vector<Vec3> points,const Vec3& direction,int arm) {
			int label=-1;
			if ((node==root || skeleton.is_terminal(node))) label=tip_label++;
			std::vector<Vec3> velocities(points.size());
			for(std::size_t q=0;q<points.size();++q) velocities[q]=VelocityAt(direction,t.circle[q]);
			const int offset=static_cast<int>(mesh.points.size());
			mesh.points.insert(mesh.points.end(),points.begin(),points.end());
			mesh.labels.insert(mesh.labels.end(),points.size(),label);
			mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
			for(int q:t.boundary_circle) mesh.labels[offset+q]=0;
			layers[node].points=std::move(points);
			layers[node].primary=OffsetFaces(t.circle_faces,offset);
			layers[node].alternate=layers[node].primary;
			layers[node].reference=cp1;
			layers[node].has_reference=true;
			info.offsets[arm]=offset;
		};
		add_terminal(parent,std::move(parent_terminal),sv_parent_parent,1);
		add_terminal(child1,std::move(child1_terminal),sv_child1,2);
		add_terminal(child2,std::move(child2_terminal),sv_child2,3);
		bifurcations.push_back(info);
	}

	if (branch_nodes.empty()) {
		const auto& section=sections.at(0);
		Vec3 w,reference;
		for(std::size_t q=0;q<section.size();++q) {
			int node=section[q];
			Vec3 direction;
			if(q==0) direction=segment[section[1]]; else direction=segment[node];
			if(q==0) {
				reference=Normalized(RotateSurface({1,0,0},{0,0,1},direction),"pipe initial reference");
				w=Normalized(Cross(direction,reference),"pipe initial frame");
			} else {
				reference=Normalized(Cross(w,direction),"pipe propagated reference");
				w=Normalized(Cross(direction,reference),"pipe propagated frame");
			}
			std::vector<Vec3> points;
			if(q==0) {
				points=Scale(t.circle,skeleton.nodes[node].diameter/2.0);
				for(auto& p:points) p=RotateSurface(p,{0,0,1},direction);
				Translate(points,skeleton.nodes[node].position-root_position);
			} else points=MakeCircle(t,skeleton.nodes[node].diameter/2.0,reference,w,skeleton.nodes[node].position-root_position);
			std::vector<Vec3> velocities(points.size());
			for(std::size_t k=0;k<points.size();++k) velocities[k]=VelocityAt(direction,t.circle[k]);
			const int offset=static_cast<int>(mesh.points.size());
			mesh.points.insert(mesh.points.end(),points.begin(),points.end());
			mesh.labels.insert(mesh.labels.end(),points.size(),(q==0||q+1==section.size())?tip_label++:-1);
			mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
			for(int k:t.boundary_circle) mesh.labels[offset+k]=0;
			layers[node].points=std::move(points);
			layers[node].primary=OffsetFaces(t.circle_faces,offset);
			layers[node].alternate=layers[node].primary;
		}
	} else {
		for(const auto& section:sections) {
			const auto trim=TrimSection(skeleton,section);
			const int start=trim.front(), end=trim.back();
			if(start==root) {
				Vec3 direction_end=segment[end];
				Vec3 reference=layers[end].reference;
				Vec3 w=Normalized(Cross(direction_end,reference),"root-bifurcation frame");
				reference=Normalized(reference,"root-bifurcation reference");
				for(std::size_t rev=trim.size()-1;rev>0;--rev) {
					const int child=trim[rev], node=trim[rev-1];
					const Vec3 direction=segment[child];
					reference=Normalized(Cross(w,direction),"root-bifurcation propagated reference");
					auto points=MakeCircle(t,skeleton.nodes[node].diameter/2.0,reference,w,skeleton.nodes[node].position-root_position);
					w=Normalized(Cross(direction,reference),"root-bifurcation propagated frame");
					std::vector<Vec3> velocities(points.size());
					for(std::size_t k=0;k<points.size();++k) velocities[k]=VelocityAt(direction,t.circle[k]);
					const int offset=static_cast<int>(mesh.points.size());
					mesh.points.insert(mesh.points.end(),points.begin(),points.end());
					mesh.labels.insert(mesh.labels.end(),points.size(),rev==1?tip_label++:-1);
					mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
					for(int k:t.boundary_circle) mesh.labels[offset+k]=0;
					layers[node].points=std::move(points);
					layers[node].primary=OffsetFaces(t.circle_faces,offset);
					layers[node].alternate=layers[node].primary;
				}
			} else if(skeleton.is_terminal(end)) {
				Vec3 reference=layers[start].reference;
				Vec3 w=Normalized(Cross(segment[start],reference),"bifurcation-terminal frame");
				reference=Normalized(reference,"bifurcation-terminal reference");
				for(std::size_t q=1;q<trim.size();++q) {
					const int node=trim[q];
					const Vec3 direction=segment[node];
					reference=Normalized(Cross(w,direction),"bifurcation-terminal propagated reference");
					w=Normalized(Cross(direction,reference),"bifurcation-terminal propagated frame");
					auto points=MakeCircle(t,skeleton.nodes[node].diameter/2.0,reference,w,skeleton.nodes[node].position-root_position);
					std::vector<Vec3> velocities(points.size());
					for(std::size_t k=0;k<points.size();++k) velocities[k]=VelocityAt(direction,t.circle[k]);
					const int offset=static_cast<int>(mesh.points.size());
					mesh.points.insert(mesh.points.end(),points.begin(),points.end());
					mesh.labels.insert(mesh.labels.end(),points.size(),q+1==trim.size()?tip_label++:-1);
					mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
					for(int k:t.boundary_circle) mesh.labels[offset+k]=0;
					layers[node].points=std::move(points);
					layers[node].primary=OffsetFaces(t.circle_faces,offset);
					layers[node].alternate=layers[node].primary;
				}
			} else {
				const Vec3 sec_start=segment[start];
				const int end_branch=skeleton.nodes[end].children.front();
				const Vec3 sec_end=segment[end_branch];
				const Vec3 ref_start=layers[start].reference, ref_end=layers[end].reference;
				const Vec3 n_start=Normalized(Cross(ref_start,Cross(sec_start,ref_start)),"bif-bif start normal");
				const Vec3 n_end=Normalized(Cross(ref_end,Cross(sec_end,ref_end)),"bif-bif end normal");
				const Vec3 mapped=RotateSurface(ref_start,n_start,n_end);
				Vec3 rotation_axis=Cross(mapped,ref_end);
				double angle=AngleBetween(mapped,ref_end);
				const double rotation_sign=Norm(rotation_axis)>1.0e-14?Dot(Normalized(rotation_axis,"bif-bif rotation axis"),n_end):0.0;
				std::cout << "section_start=" << start+1 << " section_end=" << end+1 << " twist_degrees=" << angle*180.0/std::acos(-1.0) << " sign=" << rotation_sign << "\n";
				AdjustEndFaces(layers[end].alternate,angle,rotation_sign);
				const int inserted=static_cast<int>(trim.size())-2;
				const double angle_per=angle/(inserted+1);
				const Vec3 reference_template=RotateSurface({1,0,0},{0,0,1},n_start);
				const Vec3 start_axis=Cross(reference_template,ref_start);
				const double start_angle=AngleBetween(reference_template,ref_start);
				for(std::size_t q=1;q+1<trim.size();++q) {
					const int node=trim[q];
					const Vec3 direction=segment[node];
					auto points=Scale(t.circle,skeleton.nodes[node].diameter/2.0);
					for(auto& p:points) {
						p=RotateSurface(p,{0,0,1},n_start);
						p=RotateIfNeeded(p,start_axis,start_angle);
						if(rotation_sign>0.0) p=RotateIfNeeded(p,n_start,angle_per*q);
						else if(rotation_sign<0.0) p=RotateIfNeeded(p,-n_start,angle_per*q);
						p=RotateSurface(p,n_start,direction);
						p+=skeleton.nodes[node].position-root_position;
					}
					std::vector<Vec3> velocities(points.size());
					for(std::size_t k=0;k<points.size();++k) velocities[k]=VelocityAt(direction,t.circle[k]);
					const int offset=static_cast<int>(mesh.points.size());
					mesh.points.insert(mesh.points.end(),points.begin(),points.end());
					mesh.labels.insert(mesh.labels.end(),points.size(),-1);
					mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
					for(int k:t.boundary_circle) mesh.labels[offset+k]=0;
					layers[node].points=std::move(points);
					layers[node].primary=OffsetFaces(t.circle_faces,offset);
					layers[node].alternate=layers[node].primary;
				}
			}
		}
	}

	for(const auto& section:sections) for(std::size_t q=0;q+1<section.size();++q) {
		const int parent=section[q], child=section[q+1];
		const int count=layer_count[child];
		const std::vector<Vec3>* start_points;
		FaceList start_faces;
		if(skeleton.is_branch(parent)) {
			if(child_arm[child]==1) {start_points=&branch_points[parent][1];start_faces=layers[parent].alternate;}
			else if(child_arm[child]==2) {start_points=&branch_points[parent][2];start_faces=layers[parent].right;}
			else throw std::runtime_error("missing bifurcation child arm label");
		} else {start_points=&layers[parent].points;start_faces=layers[parent].primary;}
		const std::vector<Vec3>* end_points;
		FaceList end_faces;
		if(skeleton.is_branch(child)) {end_points=&branch_points[child][0];end_faces=layers[child].primary;}
		else {end_points=&layers[child].points;end_faces=layers[child].alternate;}
		if(start_points->size()!=t.circle.size()||end_points->size()!=t.circle.size())
			throw std::runtime_error("missing cross-section layer while connecting skeleton");
		FaceList previous=start_faces;
		for(int k=1;k<count;++k) {
			std::vector<Vec3> points(t.circle.size());
			for(std::size_t i=0;i<points.size();++i)
				points[i]=((*end_points)[i]*k+(*start_points)[i]*(count-k))/count;
			std::vector<Vec3> velocities(points.size());
			for(std::size_t i=0;i<points.size();++i) velocities[i]=VelocityAt(segment[child],t.circle[i]);
			const int offset=static_cast<int>(mesh.points.size());
			mesh.points.insert(mesh.points.end(),points.begin(),points.end());
			mesh.labels.insert(mesh.labels.end(),points.size(),-1);
			mesh.velocity.insert(mesh.velocity.end(),velocities.begin(),velocities.end());
			for(int i:t.boundary_circle) mesh.labels[offset+i]=0;
			FaceList current=OffsetFaces(t.circle_faces,offset);
			ConnectFaces(mesh,previous,current);
			previous=std::move(current);
			if(skeleton.is_branch(child)&&k==count-1)
				for(auto& b:bifurcations)if(b.node==child)b.offsets[1]=offset;
			if(skeleton.is_branch(parent)&&k==1)
				for(auto& b:bifurcations)if(b.node==parent)b.offsets[child_arm[child]+1]=offset;
		}
		ConnectFaces(mesh,previous,end_faces);
	}

	const auto incident=BuildIncidentElements(mesh.points.size(),mesh.elements);
	const std::array<std::array<int,4>,2> extra{{{{1,24,123,216}},{{3,31,124,217}}}};
	for(std::size_t bi=0;bi<bifurcations.size();++bi) {
		const auto& b=bifurcations[bi];
		for(std::size_t q=0;q<t.left_boundary.size();++q) {
			std::array<std::pair<int,Vec3>,3> moves{{
				{b.offsets[0]+t.bottom_boundary[q],ProjectToLine(mesh.points[b.offsets[2]+t.right_boundary[q]],mesh.points[b.offsets[3]+t.left_boundary[q]],mesh.points[b.offsets[0]+t.bottom_boundary[q]])},
				{b.offsets[0]+t.left_boundary[q],ProjectToLine(mesh.points[b.offsets[1]+t.left_boundary[q]],mesh.points[b.offsets[2]+t.left_boundary[q]],mesh.points[b.offsets[0]+t.left_boundary[q]])},
				{b.offsets[0]+t.right_boundary[q],ProjectToLine(mesh.points[b.offsets[1]+t.right_boundary[q]],mesh.points[b.offsets[3]+t.right_boundary[q]],mesh.points[b.offsets[0]+t.right_boundary[q]])}
			}};
			for(const auto& move:moves) {
				const double alpha=MovePointSafely(mesh.points,mesh.elements,incident,move.first,move.second,minimum_scaled_jacobian);
				if(alpha<1.0) std::cout<<"limited_bifurcation="<<bi<<" point="<<move.first<<" alpha="<<alpha<<'\n';
			}
		}
		for(const auto& row:extra) {
			const Vec3 a=mesh.points[b.offsets[0]+row[1]], c=mesh.points[b.offsets[0]+row[2]], d=mesh.points[b.offsets[0]+row[3]];
			for(int offset:b.offsets) {
				const int point=offset+row[0];
				const Vec3 candidate=ProjectToPlane(a,c,d,mesh.points[point]);
				const double alpha=MovePointSafely(mesh.points,mesh.elements,incident,point,candidate,minimum_scaled_jacobian);
				if(alpha<1.0) std::cout<<"limited_bifurcation="<<bi<<" point="<<point<<" alpha="<<alpha<<'\n';
			}
		}
	}

	if(mesh.points.size()!=mesh.labels.size()||mesh.points.size()!=mesh.velocity.size())
		throw std::runtime_error("control mesh point data arrays are inconsistent");
	mesh.quality=EvaluateHexQuality(mesh.points,mesh.elements);
	if(mesh.quality.bad_elements>0)
		throw std::runtime_error("generated mesh has invalid elements; first="+std::to_string(mesh.quality.first_bad_element));
	if(mesh.quality.minimum_scaled_jacobian<minimum_scaled_jacobian)
		throw std::runtime_error("generated mesh is below the scaled-Jacobian quality floor");
	return mesh;
}

void WriteControlMeshVtk(const ControlMesh& mesh,const std::filesystem::path& path)
{
	std::ofstream out(path);
	if(!out)throw std::runtime_error("cannot write control mesh: "+path.string());
	out<<"# vtk DataFile Version 3.1\nTubularFlowIGA control mesh\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	out<<"POINTS "<<mesh.points.size()<<" FLOAT\n"<<std::fixed<<std::setprecision(6);
	for(const auto&p:mesh.points)out<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
	out<<"CELLS "<<mesh.elements.size()<<' '<<9*mesh.elements.size()<<'\n';
	for(const auto&e:mesh.elements)out<<"8 "<<e[0]<<' '<<e[1]<<' '<<e[2]<<' '<<e[3]<<' '<<e[4]<<' '<<e[5]<<' '<<e[6]<<' '<<e[7]<<'\n';
	out<<"CELL_TYPES "<<mesh.elements.size()<<'\n';
	for(std::size_t i=0;i<mesh.elements.size();++i)out<<"12\n";
	out<<"POINT_DATA "<<mesh.points.size()<<"\nSCALARS label float 1\nLOOKUP_TABLE default\n";
	for(int label:mesh.labels)out<<label<<'\n';
}

void WriteVelocity(const ControlMesh& mesh,const std::filesystem::path& path)
{
	std::ofstream out(path);
	if(!out)throw std::runtime_error("cannot write initial velocity: "+path.string());
	out<<std::fixed<<std::setprecision(6);
	for(const auto&v:mesh.velocity)out<<v.x<<' '<<v.y<<' '<<v.z<<'\n';
}

} // namespace tubular

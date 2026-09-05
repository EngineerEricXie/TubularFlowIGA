#include "MaterialSurfacePatchKinematics.hpp"
#include "PrescribedSurfaceMotion.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {

template <class F> void Reject(F&& f) { bool rejected = false; try { f(); } catch (const std::exception&) { rejected = true; } assert(rejected); }
iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::int64_t label) { iga::RawSurfaceTriangle t; t.indices = {{a,b,c}}; t.boundary_id = label; return t; }
std::vector<iga::RawSurfaceTriangle> Topology()
{ return {Face(0,2,4,7), Face(2,1,4,7), Face(1,0,4,7), Face(0,1,3,3), Face(0,3,2,3), Face(1,2,3,3)}; }
std::vector<std::array<double,3>> Reference()
{ return {{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}},{{1.0/3.0,1.0/3.0,0}}}}; }
iga::MaterialSurfaceKinematics Full(double evaluated = 0.0, double begin = -.5, double end = 0.0)
{ auto x = Reference(); return iga::MaterialSurfaceKinematics::CreateFromSourceTopology(Reference(), x, std::vector<std::array<double,3>>(x.size(), {{0,0,0}}), Topology(), evaluated, begin, end); }
iga::MaterialSurfaceKinematics FullWithChangedNonpatchLabel()
{ auto topology=Topology(); for (std::size_t i=3;i<topology.size();++i) topology[i].boundary_id=9; auto x=Reference(); return iga::MaterialSurfaceKinematics::CreateFromSourceTopology(Reference(),x,std::vector<std::array<double,3>>(x.size(),{{0,0,0}}),topology,0.,-.5,0.); }
iga::MaterialSurfaceKinematics SmallFull()
{ auto reference=Reference(); for (auto& vertex:reference) for (double& x:vertex) x*=1e-12; auto x=reference; return iga::MaterialSurfaceKinematics::CreateFromSourceTopology(reference,x,std::vector<std::array<double,3>>(x.size(),{{0,0,0}}),Topology(),0.,-.5,0.); }
using Mapping = iga::MaterialSurfacePatchMap::GlobalToSourceVertex;
std::vector<Mapping> MappingData() { return {{10,0},{11,1},{12,2},{14,4}}; }
std::vector<std::uint32_t> TriangleData() { return {0,1,2}; }
std::vector<std::uint64_t> Clamps() { return {10,11,12}; }
std::vector<std::array<std::uint64_t,3>> PatchTriangles() { return {{{10,12,14}},{{12,11,14}},{{11,10,14}}}; }
std::string ReferenceDigest(const iga::MaterialSurfaceKinematics& full, const std::vector<Mapping>& mapping = MappingData(), const std::vector<std::uint32_t>& triangles = TriangleData(), const std::vector<std::uint64_t>& clamps = Clamps(), const std::vector<std::array<std::uint64_t,3>>& layout_triangles = PatchTriangles())
{ return iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(full, 7, mapping, layout_triangles, triangles, clamps); }
iga::DistributedSurfaceLayout Layout(const iga::MaterialSurfaceKinematics& full, const std::vector<Mapping>& mapping = MappingData(), const std::vector<std::uint32_t>& triangles = TriangleData(), const std::vector<std::uint64_t>& clamps = Clamps())
{ iga::DistributedSurfaceLayout l; l.global_node_count=4; l.partition_count=1; l.partition_rank=0; l.owned_global_node_ids={10,11,12,14}; for (const auto& m:mapping) l.reference_positions.push_back({m.first,full.ReferenceMaterialVerticesM()[m.second]}); l.reference_triangles=PatchTriangles(); l.reference_mesh_identity_sha256 = ReferenceDigest(full,mapping,triangles,clamps,l.reference_triangles); l.owned_reference_lumped_areas_m2={.1,.1,.1,.1}; l.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(l); return l; }
iga::DistributedSurfaceInterface Interface(const std::string& digest)
{ iga::DistributedSurfaceInterface i; i.id={"structure","membrane","patch"}; i.subsystem_id="membrane"; i.boundary_labels={7}; i.reference_mesh_identity_sha256=digest; i.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity}; i.requires={iga::SurfaceFieldQuantity::TractionOnStructure}; return i; }
iga::MaterialSurfacePatchMap Map(const iga::MaterialSurfaceKinematics& full)
{ const auto layout=Layout(full); return iga::MaterialSurfacePatchMap::Create(Interface(layout.reference_mesh_identity_sha256),layout,full,7,MappingData(),TriangleData(),Clamps()); }
iga::SurfaceKinematics Trial(const iga::MaterialSurfacePatchMap& map, const iga::FsiTrialContext& c)
{ iga::SurfaceKinematics p; p.interface=map.Interface().id; p.stamp.time_s=c.EndTime(); p.stamp.step=c.step; p.stamp.coupling_iteration=c.coupling_iteration; p.stamp.reference_mesh_identity_sha256=map.Layout().reference_mesh_identity_sha256; p.stamp.layout_identity_sha256=map.Layout().layout_identity_sha256; p.stamp.partition_identity_sha256=iga::BuildDistributedSurfacePartitionIdentitySha256(map.Layout()); iga::Sha256 h; h.Append("patch",5); p.stamp.producer_state_identity_sha256=h.Hex(); p.displacement_m.assign(4,{{0,0,0}}); p.velocity_m_per_s.assign(4,{{0,0,0}}); p.displacement_m[3][2]=.1; p.velocity_m_per_s[3][2]=.2; return p; }

void FactoryRoundTrip()
{ iga::RawSurfaceSoup soup; soup.vertices={{{{0,0,0}},{{1,0,0}},{{0,1,0}},{{0,0,1}}}}; soup.triangles={Face(0,2,1,7),Face(0,1,3,3),Face(0,3,2,3),Face(1,2,3,3)}; iga::PrescribedSurfaceMotion motion({{0.,soup},{1.,soup}}); const auto p=motion.Evaluate(.5,0.,1.); std::vector<iga::RawSurfaceTriangle> t; for(const auto& s:p.SourceTriangles()){iga::RawSurfaceTriangle r; for(std::size_t i=0;i<3;++i)r.indices[i]=s.source_vertex_indices[i]; r.boundary_id=s.boundary_id;t.push_back(r);} const auto n=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(p.ReferenceMaterialVerticesM(),p.SourceVerticesM(),p.SourceVertexVelocitiesMPerS(),t,p.EvaluatedTimeS(),p.StepStartS(),p.StepEndS()); assert(n.MaterialIdentitySha256()==p.MaterialIdentitySha256()&&n.TopologyIdentitySha256()==p.TopologyIdentitySha256()&&n.ContentIdentitySha256()==p.ContentIdentitySha256()); assert(p.IdentitySha256()!=n.IdentitySha256()); }
}

int main()
{
	FactoryRoundTrip(); const auto full=Full(); const auto map=Map(full); assert(map.SourceVertexForGlobalNode(14)==4&&map.CanonicalTriangleForLayoutTriangle(1)<full.Surface().Triangles().size());
	// Public digest construction is a validated API: malformed mapping inputs
	// must throw before any source vertex or triangle is dereferenced.
	Reject([&]{auto m=MappingData();m[3].second=99;ReferenceDigest(full,m);});
	Reject([&]{auto t=TriangleData();t[1]=99;ReferenceDigest(full,MappingData(),t);});
	Reject([&]{auto t=TriangleData();t[1]=0;ReferenceDigest(full,MappingData(),t);});
	Reject([&]{ReferenceDigest(full,MappingData(),TriangleData(),{10,999});});
	// Explicit mapping mutations are rejected without any coordinate search fallback.
	Reject([&]{auto m=MappingData();m[1].first=10;auto l=Layout(full,m);iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,m,TriangleData(),Clamps());});
	Reject([&]{auto l=Layout(full);l.reference_positions[3].position_m[0]+=.01;l.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(l);iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,MappingData(),TriangleData(),Clamps());});
	Reject([&]{auto l=Layout(full);l.reference_triangles[0]={{10,14,12}};l.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(l);iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,MappingData(),TriangleData(),Clamps());});
	{ auto l=Layout(full);l.reference_triangles[0]={{12,14,10}};l.reference_mesh_identity_sha256=ReferenceDigest(full,MappingData(),TriangleData(),Clamps(),l.reference_triangles);l.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(l);const auto cyclic=iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,MappingData(),TriangleData(),Clamps());assert(cyclic.ReferenceIdentitySha256()!=map.ReferenceIdentitySha256()); }
	Reject([&]{auto l=Layout(full);auto i=Interface(l.reference_mesh_identity_sha256);i.boundary_labels={7,8};iga::MaterialSurfacePatchMap::Create(i,l,full,7,MappingData(),TriangleData(),Clamps());});
	Reject([&]{auto l=Layout(full);auto i=Interface(l.reference_mesh_identity_sha256);i.provides={iga::SurfaceFieldQuantity::Displacement};i.requires={iga::SurfaceFieldQuantity::Velocity,iga::SurfaceFieldQuantity::TractionOnStructure};iga::MaterialSurfacePatchMap::Create(i,l,full,7,MappingData(),TriangleData(),Clamps());});
	Reject([&]{auto clamps=std::vector<std::uint64_t>{10,11};auto l=Layout(full,MappingData(),TriangleData(),clamps);iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,MappingData(),TriangleData(),clamps);});
	Reject([&]{auto l=Layout(full);l.partition_count=2;l.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(l);iga::MaterialSurfacePatchMap::Create(Interface(l.reference_mesh_identity_sha256),l,full,7,MappingData(),TriangleData(),Clamps());});
	// Patch reference authority is local to the selected mesh; whole-material
	// and interface/mapping bindings remain visible in the map identity.
	const auto changed_nonpatch=FullWithChangedNonpatchLabel(); const auto changed_map=Map(changed_nonpatch);
	assert(changed_map.ReferenceIdentitySha256()==map.ReferenceIdentitySha256());
	assert(changed_map.IdentitySha256()!=map.IdentitySha256());
	{ auto i=Interface(map.ReferenceIdentitySha256());i.id.interface_id="other-patch";const auto other=iga::MaterialSurfacePatchMap::Create(i,Layout(full),full,7,MappingData(),TriangleData(),Clamps());assert(other.IdentitySha256()!=map.IdentitySha256()); }
	{ auto m=MappingData();std::swap(m[0].second,m[1].second);const std::vector<std::array<std::uint64_t,3>> remapped={{{11,12,14}},{{12,10,14}},{{10,11,14}}};assert(ReferenceDigest(full,m,TriangleData(),Clamps(),remapped)!=ReferenceDigest(full)); }

	const iga::FsiTrialContext context{9,0.,.5,2}; const auto patch=Trial(map,context); const auto result=iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,patch,context);
	assert(result.target.SourceVerticesM()[4][2]==.1 && result.target.SourceVertexVelocitiesMPerS()[4][2]==.2);
	const std::array<double,3> zero{{0,0,0}};
	assert(result.target.SourceVerticesM()[3]==Reference()[3] && result.target.SourceVertexVelocitiesMPerS()[3]==zero);
	assert(result.target.MaterialIdentitySha256()==full.MaterialIdentitySha256()&&result.target.TopologyIdentitySha256()==full.TopologyIdentitySha256());
	assert(iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,patch,context).composition_identity_sha256==result.composition_identity_sha256);
	{ auto changed=patch;iga::Sha256 h;h.Append("patch-mutation",14);changed.stamp.producer_state_identity_sha256=h.Hex();assert(iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,changed,context).composition_identity_sha256!=result.composition_identity_sha256); }
	{ const iga::FsiTrialContext next{10,.5,.5,2};auto next_patch=Trial(map,next);next_patch.displacement_m[3][2]=.2;const auto next_result=iga::MaterialSurfacePatchKinematics::ComposeTarget(map,result.target,next_patch,next);assert(next_result.target.SourceVerticesM()[4][2]==.2&&next_result.composition_identity_sha256!=result.composition_identity_sha256); }
	{ const iga::FsiTrialContext other_iteration{9,0.,.5,3};const auto retried=iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,Trial(map,other_iteration),other_iteration);assert(retried.composition_identity_sha256!=result.composition_identity_sha256); }
	Reject([&]{auto bad=patch;bad.velocity_m_per_s[3][2]=.1;iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,bad,context);});
	Reject([&]{auto bad=patch;bad.displacement_m[0][0]=1e-3;iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,bad,context);});
	Reject([&]{auto bad=patch;bad.stamp.coupling_iteration=3;iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,bad,context);});
	Reject([&]{auto bad=context;bad.dt_s=.25;iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,patch,bad);});
	Reject([&]{auto bad=patch;bad.displacement_m[3]={{-1.0/3.0,-1.0/3.0,0}};bad.velocity_m_per_s[3]={{-2.0/3.0,-2.0/3.0,0}};iga::MaterialSurfacePatchKinematics::ComposeTarget(map,full,bad,context);});
	{ const auto small=SmallFull();const auto small_map=Map(small);auto small_patch=Trial(small_map,context);small_patch.displacement_m.assign(4,{{0,0,0}});small_patch.velocity_m_per_s.assign(4,{{0,0,0}});small_patch.displacement_m[3][0]=2e-16;Reject([&]{iga::MaterialSurfacePatchKinematics::ComposeTarget(small_map,small,small_patch,context);}); }
	return 0;
}

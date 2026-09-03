#include "MovingImmersedFlowSnapshotPublisher.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Check(bool value,const char* message){if(!value)throw std::runtime_error(message);}
std::string Hash(char c){return std::string(64,c);}
iga::MovingImmersedFlowSnapshotPublicationIdentity Identity(){return {Hash('c'),Hash('d'),Hash('e'),Hash('f'),Hash('0')};}
std::string Read(const std::filesystem::path& path){std::ifstream input(path);return {std::istreambuf_iterator<char>(input),{}};}
iga::RawSurfaceTriangle Face(std::int64_t a,std::int64_t b,std::int64_t c,int label){iga::RawSurfaceTriangle value;value.indices={{a,b,c}};value.boundary_id=label;return value;}
iga::RawSurfaceSoup Tetra(){iga::RawSurfaceSoup value;value.vertices={{{{.18,.25,.25}},{{.78,.25,.25}},{{.18,.85,.25}},{{.18,.25,.85}}}};value.triangles={Face(0,2,1,1),Face(0,1,3,7),Face(0,3,2,9),Face(1,2,3,1)};return value;}
struct Fixture {iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};std::unique_ptr<iga::MovingCutGeometry> geometry;iga::ImmersedActiveLayout layout;};
Fixture MakeFixture(){Fixture value;iga::MovingCutGeometryOptions options;options.volume.max_depth=3;options.volume.max_nodes=options.volume.max_leaves=options.volume.max_points=1000000;options.volume.max_records=options.volume.max_logical_points=1000000;options.volume.max_retained_bytes=100000000;iga::PrescribedSurfaceMotion motion({{0.,Tetra()},{1.,Tetra()}});value.geometry=iga::MovingCutGeometry::Build(value.grid,motion.Evaluate(0.,0.,1.),options);value.layout=iga::ImmersedActiveLayout::Build(value.geometry->Domain(),value.geometry->Volume(),value.geometry->GeometryIdentitySha256(),{9});return value;}
iga::MovingImmersedFlowSnapshot Snapshot(const Fixture& fixture,std::uint64_t index,bool flows_available,double flow){std::vector<std::array<double,4>> coefficients(fixture.layout.NodeIds().size());iga::ImmersedGlobalFlowState state(0.,index,fixture.layout,coefficients,{0.},false,0.);iga::MovingImmersedFlowSnapshotRequest request;request.time_s=0.;request.index=index;request.transition_available=true;request.dt_s=.25;request.wall_labels={1};request.port_labels={7,9};request.port_flows_available=flows_available;if(flows_available)request.port_flows={{7,flow},{9,0.}};return iga::MovingImmersedFlowSnapshot::Build(*fixture.geometry,fixture.layout,state,request,{});}
void ReplaceOnce(const std::filesystem::path& path,const std::string& from,const std::string& to){std::string text=Read(path);const auto position=text.find(from);Check(position!=std::string::npos,"test manifest field was absent");text.replace(position,from.size(),to);std::ofstream output(path,std::ios::trunc);output<<text;}
}

int main()
{
	try {
		using iga::MovingImmersedFlowSnapshotPublisher;
		const auto stem=std::filesystem::path("/tmp/moving snapshot & <test>/field");
		Check(MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,0).filename()=="field.epoch00000000000000000000","zero epoch naming changed");
		Check(MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem,UINT64_MAX).filename()=="field.epoch18446744073709551615","uint64 epoch naming changed");
		Check(MovingImmersedFlowSnapshotPublisher::CollectionPath(stem).filename()=="field.pvd","PVD naming changed");
		bool rejected=false;try{(void)MovingImmersedFlowSnapshotPublisher::EpochDirectory({},0);}catch(const std::invalid_argument&){rejected=true;}Check(rejected,"empty output stem was accepted");
		const auto root=std::filesystem::temp_directory_path()/"moving_snapshot_publisher_test";std::filesystem::remove_all(root);const auto output=root/"quoted & <name>";
		const auto fixture=MakeFixture();const auto unavailable=Snapshot(fixture,4,false,0.);const auto known_zero=Snapshot(fixture,5,true,0.);const auto positive=Snapshot(fixture,3,true,-.25);
		Check(unavailable.ContentHashSha256()!=known_zero.ContentHashSha256()&&unavailable.SnapshotIdentitySha256()!=known_zero.SnapshotIdentitySha256(),"unavailable and known-zero snapshots share an identity");
		const auto first=MovingImmersedFlowSnapshotPublisher::Publish(unavailable,output,Identity());const auto zero_publication=MovingImmersedFlowSnapshotPublisher::Publish(known_zero,output,Identity());const auto positive_publication=MovingImmersedFlowSnapshotPublisher::Publish(positive,output,Identity());
		Check(std::filesystem::is_regular_file(first.vtu_path)&&std::filesystem::is_regular_file(first.metrics_path),"epoch files missing");
		const auto vtu=Read(first.vtu_path),unavailable_json=Read(first.metrics_path),zero_json=Read(zero_publication.metrics_path),positive_json=Read(positive_publication.metrics_path),pvd=Read(first.pvd_path);
		Check(vtu.find("moving_snapshot_schema_version\" format=\"ascii\">\n          3")!=std::string::npos,"VTU schema v3 missing");
		Check(unavailable_json.find("\"schema_version\": 3")!=std::string::npos&&unavailable_json.find("\"transition_available\": true")!=std::string::npos&&unavailable_json.find("\"dt_s\": 0.25")!=std::string::npos&&unavailable_json.find("\"wall_labels\": [1]")!=std::string::npos&&unavailable_json.find("\"port_labels\": [7, 9]")!=std::string::npos&&unavailable_json.find("\"port_flows_available\": false")!=std::string::npos&&unavailable_json.find("\"endpoint_turnover_available\": false")!=std::string::npos&&unavailable_json.find("\"port_label_outward_flow_m3_s\": []")!=std::string::npos,"unavailable v3 request semantics missing");
		Check(zero_json.find("\"port_flows_available\": true")!=std::string::npos&&zero_json.find("\"endpoint_turnover_available\": true")!=std::string::npos&&zero_json.find("{\"label\": 7, \"outward_flow_m3_s\": 0}")!=std::string::npos,"known-zero v3 request semantics missing");
		Check(positive_json.find("{\"label\": 7, \"outward_flow_m3_s\": -0.25}")!=std::string::npos&&first.metrics_sha256!=zero_publication.metrics_sha256,"available port-flow serialization or identity did not differ");
		Check(pvd.find("quoted &amp; &lt;name&gt;.epoch00000000000000000004/fields.vtu")!=std::string::npos,"PVD did not escape name");
		const auto retry=MovingImmersedFlowSnapshotPublisher::Publish(unavailable,output,Identity());Check(first.vtu_sha256==retry.vtu_sha256&&first.metrics_sha256==retry.metrics_sha256,"idempotent retry changed epoch");
		auto conflict=Snapshot(fixture,4,true,0.);rejected=false;try{(void)MovingImmersedFlowSnapshotPublisher::Publish(conflict,output,Identity());}catch(const std::runtime_error&){rejected=true;}Check(rejected,"index with different semantic identity was accepted");
		ReplaceOnce(positive_publication.metrics_path,"\"endpoint_turnover_available\": true","\"endpoint_turnover_available\": false");rejected=false;try{MovingImmersedFlowSnapshotPublisher::RebuildCollection(output);}catch(const std::runtime_error&){rejected=true;}Check(rejected,"recovery accepted mismatched endpoint availability");
		std::filesystem::remove_all(root);
	} catch(const std::exception&) {return 1;}
	return 0;
}

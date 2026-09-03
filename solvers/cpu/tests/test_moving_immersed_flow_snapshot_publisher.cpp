#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#define private public
#include "MovingImmersedFlowSnapshotPublisher.hpp"
#undef private

namespace {
void Check(bool value, const char* message)
{
	if (!value) throw std::runtime_error(message);
}

std::string Hash(char c) { return std::string(64, c); }
iga::MovingImmersedFlowSnapshot Snapshot(std::uint64_t index, double time)
{
	iga::MovingImmersedFlowSnapshot value;
	value.options_.maximum_points=8; value.options_.maximum_output_bytes=1024*1024;
	value.request_.index=index; value.request_.time_s=time;
	value.request_.port_flows={{7,-.25},{9,.125}};
	iga::MovingImmersedFieldPoint point;
	point.physical_m={{1.,2.,3.}}; point.parametric={{.1,.2,.3}};
	point.velocity_m_per_s={{.4,.5,.6}}; point.vorticity_per_s={{.7,.8,.9}};
	point.pressure=2.; point.speed_m_per_s=1.; point.q_criterion_per_s2=-.1;
	point.enstrophy_density_per_s2=.2; point.physical_integration_weight_m3=.3;
	point.background_cell_id=12; point.cell_quadrature_ordinal=4;
	point.cell_kind=iga::MovingImmersedFieldCellKind::Cut; point.stagnant=true;
	value.points_={point,point}; value.points_[1].background_cell_id=13;
	value.metrics_.quadrature_volume_m3=.6; value.metrics_.audited_volume_m3=.6;
	value.metrics_.inlet_flow_m3_s=.25; value.metrics_.outlet_flow_m3_s=.125;
	value.metrics_.endpoint_turnover_rate_per_s=.25/.6;
	value.content_hash_sha256_=Hash('a'); value.snapshot_identity_sha256_=Hash('b'); return value;
}
iga::MovingImmersedFlowSnapshotPublicationIdentity Identity() { return {Hash('c'),Hash('d'),Hash('e'),Hash('f'),Hash('0')}; }
std::string Read(const std::filesystem::path& path) { std::ifstream input(path); return {std::istreambuf_iterator<char>(input), {}}; }
}

int main()
{
	try {
		using iga::MovingImmersedFlowSnapshotPublisher;
		const auto stem = std::filesystem::path("/tmp/moving snapshot & <test>/field");
		Check(MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem, 0).filename()
			== "field.epoch00000000000000000000", "zero epoch naming changed");
		Check(MovingImmersedFlowSnapshotPublisher::EpochDirectory(stem, UINT64_MAX).filename()
			== "field.epoch18446744073709551615", "uint64 epoch naming changed");
		Check(MovingImmersedFlowSnapshotPublisher::CollectionPath(stem).filename()=="field.pvd", "PVD naming changed");
		bool rejected=false; try { (void)MovingImmersedFlowSnapshotPublisher::EpochDirectory({}, 0); } catch (const std::invalid_argument&) { rejected=true; }
		Check(rejected, "empty output stem was accepted");
		const auto root=std::filesystem::temp_directory_path()/"moving_snapshot_publisher_test";
		std::filesystem::remove_all(root); const auto output=root/"quoted & <name>";
		const auto snapshot=Snapshot(4,.5); const auto before=snapshot.ContentHashSha256();
		const auto first=MovingImmersedFlowSnapshotPublisher::Publish(snapshot,output,Identity());
		Check(std::filesystem::is_regular_file(first.vtu_path)&&std::filesystem::is_regular_file(first.metrics_path),"epoch files missing");
		const auto vtu=Read(first.vtu_path), json=Read(first.metrics_path), pvd=Read(first.pvd_path);
		Check(vtu.find("NumberOfPoints=\"2\" NumberOfCells=\"2\"")!=std::string::npos&&vtu.find("Name=\"velocity\" NumberOfComponents=\"3\"")!=std::string::npos,"VTU schema missing");
		Check(vtu.find("Name=\"types\"")!=std::string::npos&&vtu.find("1 1")!=std::string::npos,"VTK_VERTEX cells missing");
		Check(json.find("\"endpoint_turnover_time_s\": null")!=std::string::npos&&json.find("\"vtu_sha256\": ")!=std::string::npos,"metrics schema missing");
		Check(pvd.find("quoted &amp; &lt;name&gt;.epoch00000000000000000004/fields.vtu")!=std::string::npos,"PVD did not escape name");
		const auto retry=MovingImmersedFlowSnapshotPublisher::Publish(snapshot,output,Identity());
		Check(first.vtu_sha256==retry.vtu_sha256&&first.metrics_sha256==retry.metrics_sha256,"idempotent retry changed epoch");
		auto conflict=Snapshot(4,.5); conflict.snapshot_identity_sha256_=Hash('9');
		bool failed=false;try{(void)MovingImmersedFlowSnapshotPublisher::Publish(conflict,output,Identity());}catch(const std::runtime_error&){failed=true;}
		Check(failed,"duplicate epoch index with a different identity was accepted");
		(void)MovingImmersedFlowSnapshotPublisher::Publish(Snapshot(3,.4),output,Identity());
		const auto sorted=Read(MovingImmersedFlowSnapshotPublisher::CollectionPath(output));
		Check(sorted.find(".epoch00000000000000000003/fields.vtu") < sorted.find(".epoch00000000000000000004/fields.vtu"),"collection was not sorted by index");
		MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::BeforeDirectoryRename);
		failed=false;try{(void)MovingImmersedFlowSnapshotPublisher::Publish(Snapshot(6,.7),output,Identity());}catch(const std::runtime_error&){failed=true;}
		Check(failed&&!std::filesystem::exists(MovingImmersedFlowSnapshotPublisher::EpochDirectory(output,6)),"pre-rename failure published a partial epoch");
		MovingImmersedFlowSnapshotPublisher::SetTestFault(iga::MovingImmersedFlowSnapshotPublicationFault::AfterDirectoryRename);
		failed=false;try{(void)MovingImmersedFlowSnapshotPublisher::Publish(Snapshot(7,.8),output,Identity());}catch(const std::runtime_error&){failed=true;}
		Check(failed&&std::filesystem::is_directory(MovingImmersedFlowSnapshotPublisher::EpochDirectory(output,7)),"post-rename failure did not leave recoverable epoch");
		MovingImmersedFlowSnapshotPublisher::RebuildCollection(output);
		Check(Read(MovingImmersedFlowSnapshotPublisher::CollectionPath(output)).find(".epoch00000000000000000007/fields.vtu")!=std::string::npos,"orphan recovery failed");
		std::ofstream(first.vtu_path,std::ios::trunc) << "truncated";
		failed=false;try{MovingImmersedFlowSnapshotPublisher::RebuildCollection(output);}catch(const std::runtime_error&){failed=true;}
		Check(failed,"corrupt complete epoch was silently included");
		Check(snapshot.ContentHashSha256()==before,"publisher mutated source snapshot");
		std::filesystem::remove_all(root);
	} catch (const std::exception&) { return 1; }
	return 0;
}

#ifndef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_HPP
#define IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_HPP

// Transactional publication for an already-built MovingImmersedFlowSnapshot.
// This is deliberately independent of the moving runtime: it only reads its
// immutable input and makes each completed epoch recoverable from disk.
#include "MovingImmersedFlowSnapshot.hpp"
#include "Sha256.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace iga {

struct MovingImmersedFlowSnapshotPublicationIdentity {
	std::string geometry_identity_sha256;
	std::string publication_identity_sha256;
	std::string layout_identity_sha256;
	std::string state_identity_sha256;
	std::string runtime_identity_sha256;
};

// A caller that has a committed transition conservation audit may attach it
// without making the publisher inspect mutable runtime state.
struct MovingImmersedFlowSnapshotTransitionConservation {
	std::string identity_sha256;
	double reynolds_residual_m3_s = 0.0;
	double moving_residual_m3_s = 0.0;
	double volume_change_rate_m3_s = 0.0;
};

struct MovingImmersedFlowSnapshotPublication {
	std::filesystem::path epoch_directory;
	std::filesystem::path vtu_path;
	std::filesystem::path metrics_path;
	std::filesystem::path pvd_path;
	std::string vtu_sha256;
	std::string metrics_sha256;
};

#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
enum class MovingImmersedFlowSnapshotPublicationFault {
	None, AfterTemporaryVtuWrite, BeforeDirectoryRename, AfterDirectoryRename
};
#endif

class MovingImmersedFlowSnapshotPublisher {
public:
	static constexpr std::uint32_t SchemaVersion() noexcept { return 1; }

	static std::filesystem::path EpochDirectory(const std::filesystem::path& stem, std::uint64_t index)
	{
		ValidateStem(stem);
		std::ostringstream name;
		name << stem.filename().string() << ".epoch" << std::setw(20) << std::setfill('0') << index;
		return stem.parent_path()/name.str();
	}
	static std::filesystem::path CollectionPath(const std::filesystem::path& stem)
	{
		ValidateStem(stem); return stem.parent_path()/(stem.filename().string()+".pvd");
	}

	static MovingImmersedFlowSnapshotPublication Publish(const MovingImmersedFlowSnapshot& snapshot,
		const std::filesystem::path& stem, const MovingImmersedFlowSnapshotPublicationIdentity& identity,
		const std::optional<MovingImmersedFlowSnapshotTransitionConservation>& conservation = std::nullopt)
	{
		ValidateSnapshot(snapshot, identity, conservation); ValidateStem(stem);
		const auto final_directory = EpochDirectory(stem, snapshot.Request().index);
		const auto parent = final_directory.parent_path();
		std::filesystem::create_directories(parent);
		if (!std::filesystem::is_directory(parent)) throw std::runtime_error("moving snapshot output parent is not a directory");
		if (std::filesystem::exists(final_directory)) {
			const auto published = ValidateEpoch(final_directory, stem, snapshot.Request().index);
			if (published.snapshot_identity != snapshot.SnapshotIdentitySha256()
				|| published.content_hash != snapshot.ContentHashSha256())
				throw std::runtime_error("moving snapshot epoch index conflicts with a different identity");
			const auto manifest = ReadAll(final_directory/"metrics.json");
			if (JsonString(manifest,"geometry_identity_sha256") != identity.geometry_identity_sha256
				|| JsonString(manifest,"publication_identity_sha256") != identity.publication_identity_sha256
				|| JsonString(manifest,"layout_identity_sha256") != identity.layout_identity_sha256
				|| JsonString(manifest,"state_identity_sha256") != identity.state_identity_sha256
				|| JsonString(manifest,"runtime_identity_sha256") != identity.runtime_identity_sha256)
				throw std::runtime_error("moving snapshot epoch index conflicts with a different publication identity");
			RebuildCollection(stem);
			return Result(final_directory, CollectionPath(stem), published.vtu_sha256, published.metrics_sha256);
		}
		const auto temporary = TemporaryDirectory(final_directory);
		bool renamed = false;
		try {
			std::filesystem::create_directory(temporary);
			const auto vtu = temporary/"fields.vtu";
			WriteVtu(vtu, snapshot, identity);
			const auto vtu_hash = FileHash(vtu);
#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
			ThrowFault(MovingImmersedFlowSnapshotPublicationFault::AfterTemporaryVtuWrite);
#endif
			const auto metrics = temporary/"metrics.json";
			WriteMetrics(metrics, snapshot, identity, conservation, vtu_hash);
			const auto metrics_hash = FileHash(metrics);
			// Validate the exact on-disk manifest before publishing the directory.
			const auto checked = ValidateEpoch(temporary, stem, snapshot.Request().index);
			if (checked.vtu_sha256 != vtu_hash || checked.metrics_sha256 != metrics_hash)
				throw std::runtime_error("moving snapshot temporary epoch hash mismatch");
#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
			ThrowFault(MovingImmersedFlowSnapshotPublicationFault::BeforeDirectoryRename);
#endif
			std::filesystem::rename(temporary, final_directory);
			renamed = true;
#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
			ThrowFault(MovingImmersedFlowSnapshotPublicationFault::AfterDirectoryRename);
#endif
			RebuildCollection(stem);
			return Result(final_directory, CollectionPath(stem), vtu_hash, metrics_hash);
		} catch (...) {
			if (!renamed) { std::error_code ignored; std::filesystem::remove_all(temporary, ignored); }
			throw;
		}
	}

	// This is intentionally public so a process can recover a complete orphan
	// after a crash between the directory and PVD renames.
	static void RebuildCollection(const std::filesystem::path& stem)
	{
		ValidateStem(stem); const auto parent = stem.parent_path();
		if (!std::filesystem::exists(parent)) return;
		std::vector<Epoch> epochs;
		const std::string prefix = stem.filename().string()+".epoch";
		for (const auto& entry : std::filesystem::directory_iterator(parent)) {
			if (!entry.is_directory()) continue;
			const auto name = entry.path().filename().string();
			if (name.compare(0, prefix.size(), prefix) != 0) continue;
			const auto suffix = name.substr(prefix.size());
			if (suffix.size() != 20 || !std::all_of(suffix.begin(), suffix.end(), [](unsigned char c){ return std::isdigit(c) != 0; }))
				throw std::runtime_error("moving snapshot epoch directory name is invalid: "+name);
			std::uint64_t index = 0;
			for (const char c : suffix) {
				if (index > (std::numeric_limits<std::uint64_t>::max()-static_cast<unsigned>(c-'0'))/10)
					throw std::runtime_error("moving snapshot epoch index overflows");
				index = 10*index+static_cast<unsigned>(c-'0');
			}
			epochs.push_back(ValidateEpoch(entry.path(), stem, index));
		}
		std::sort(epochs.begin(), epochs.end(), [](const Epoch& a, const Epoch& b) { return a.index < b.index; });
		for (std::size_t i=1; i<epochs.size(); ++i) if (epochs[i-1].index == epochs[i].index)
			throw std::runtime_error("moving snapshot collection has duplicate epoch index");
		const auto output = CollectionPath(stem);
		const std::filesystem::path temporary = output.string()+".tmp."+Nonce();
		try {
			std::ofstream file(temporary, std::ios::binary|std::ios::trunc);
			if (!file) throw std::runtime_error("cannot create moving snapshot PVD");
			file << std::setprecision(17) << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <Collection>\n";
			for (const auto& epoch : epochs) file << "    <DataSet timestep=\"" << epoch.time_s << "\" group=\"\" part=\"0\" file=\""
				<< EscapeXml(epoch.directory.filename().string()+"/fields.vtu") << "\"/>\n";
			file << "  </Collection>\n</VTKFile>\n"; file.close();
			if (!file) throw std::runtime_error("cannot write moving snapshot PVD");
			std::filesystem::rename(temporary, output);
		} catch (...) { std::error_code ignored; std::filesystem::remove(temporary, ignored); throw; }
	}

#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
	static void SetTestFault(MovingImmersedFlowSnapshotPublicationFault fault) { TestFault() = fault; }
#endif

private:
	struct Epoch { std::filesystem::path directory; std::uint64_t index = 0; double time_s = 0.; std::string snapshot_identity, content_hash, vtu_sha256, metrics_sha256; };
	static bool ValidHash(const std::string& value)
	{
		return value.size() == 64 && std::all_of(value.begin(), value.end(), [](unsigned char c) { return std::isxdigit(c) != 0; });
	}
	static void RequireHash(const std::string& value, const char* what) { if (!ValidHash(value)) throw std::invalid_argument(std::string("moving snapshot ")+what+" must be a SHA-256 hex digest"); }
	static void ValidateStem(const std::filesystem::path& stem)
	{
		if (stem.empty() || stem.filename().empty() || stem.filename() == "." || stem.filename() == ".." || stem.filename().string().find('/') != std::string::npos || stem.filename().string().find('\\') != std::string::npos)
			throw std::invalid_argument("moving snapshot output stem is invalid");
	}
	static void ValidateSnapshot(const MovingImmersedFlowSnapshot& snapshot, const MovingImmersedFlowSnapshotPublicationIdentity& id,
		const std::optional<MovingImmersedFlowSnapshotTransitionConservation>& conservation)
	{
		RequireHash(snapshot.ContentHashSha256(), "content identity"); RequireHash(snapshot.SnapshotIdentitySha256(), "identity");
		for (const auto* item : {&id.geometry_identity_sha256, &id.publication_identity_sha256, &id.layout_identity_sha256, &id.state_identity_sha256, &id.runtime_identity_sha256}) RequireHash(*item, "publication identity");
		if (!std::isfinite(snapshot.Request().time_s) || snapshot.Points().empty() || snapshot.Points().size() > snapshot.Options().maximum_points)
			throw std::invalid_argument("moving snapshot publication input is invalid");
		for (const auto& p : snapshot.Points()) for (double x : {p.physical_m[0],p.physical_m[1],p.physical_m[2],p.parametric[0],p.parametric[1],p.parametric[2],p.velocity_m_per_s[0],p.velocity_m_per_s[1],p.velocity_m_per_s[2],p.vorticity_per_s[0],p.vorticity_per_s[1],p.vorticity_per_s[2],p.pressure,p.speed_m_per_s,p.q_criterion_per_s2,p.enstrophy_density_per_s2,p.physical_integration_weight_m3}) if (!std::isfinite(x)) throw std::invalid_argument("moving snapshot publication has nonfinite point data");
		if (conservation) { RequireHash(conservation->identity_sha256, "conservation identity"); for (double x : {conservation->reynolds_residual_m3_s, conservation->moving_residual_m3_s, conservation->volume_change_rate_m3_s}) if (!std::isfinite(x)) throw std::invalid_argument("moving snapshot conservation value is nonfinite"); }
	}
	static std::string EscapeXml(const std::string& input) { std::string out; for (unsigned char c : input) { if(c=='&')out+="&amp;"; else if(c=='<')out+="&lt;"; else if(c=='>')out+="&gt;"; else if(c=='\"')out+="&quot;"; else if(c=='\'')out+="&apos;"; else out.push_back(static_cast<char>(c)); } return out; }
	static std::string EscapeJson(const std::string& input) { std::ostringstream out; for(unsigned char c:input) { switch(c){case '\"':out<<"\\\\\"";break;case '\\':out<<"\\\\\\\\";break;case '\b':out<<"\\\\b";break;case '\f':out<<"\\\\f";break;case '\n':out<<"\\\\n";break;case '\r':out<<"\\\\r";break;case '\t':out<<"\\\\t";break;default:if(c<0x20)out<<"\\\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c)<<std::dec<<std::setfill(' ');else out<<static_cast<char>(c);}} return out.str(); }
	static std::string FileHash(const std::filesystem::path& path) { std::ifstream file(path, std::ios::binary); if(!file)throw std::runtime_error("cannot hash moving snapshot file"); Sha256 hash; std::array<char,8192> buffer{}; while(file.read(buffer.data(),buffer.size()) || file.gcount()) hash.Append(buffer.data(),static_cast<std::size_t>(file.gcount())); if(!file.eof())throw std::runtime_error("cannot read moving snapshot file"); return hash.Hex(); }
	static void Data(std::ostream& out, const char* type, const char* name, unsigned components, const std::string& values) { out << "        <DataArray type=\""<<type<<"\" Name=\""<<name<<"\""; if(components!=1)out<<" NumberOfComponents=\""<<components<<"\""; out<<" format=\"ascii\">\n          "<<values<<"\n        </DataArray>\n"; }
	static void WriteVtu(const std::filesystem::path& path, const MovingImmersedFlowSnapshot& s, const MovingImmersedFlowSnapshotPublicationIdentity& id)
	{
		std::ofstream out(path, std::ios::binary|std::ios::trunc); if(!out)throw std::runtime_error("cannot create moving snapshot VTU"); out<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <UnstructuredGrid>\n    <FieldData>\n";
		Data(out,"String","moving_snapshot_schema",1,"MovingImmersedFlowSnapshotPublication"); Data(out,"UInt32","moving_snapshot_schema_version",1,std::to_string(SchemaVersion())); Data(out,"Float64","time_s",1,Number(s.Request().time_s)); Data(out,"UInt64","index",1,std::to_string(s.Request().index));
		Data(out,"String","geometry_identity_sha256",1,id.geometry_identity_sha256); Data(out,"String","publication_identity_sha256",1,id.publication_identity_sha256); Data(out,"String","layout_identity_sha256",1,id.layout_identity_sha256); Data(out,"String","state_identity_sha256",1,id.state_identity_sha256); Data(out,"String","runtime_identity_sha256",1,id.runtime_identity_sha256); Data(out,"String","content_hash_sha256",1,s.ContentHashSha256()); Data(out,"String","snapshot_identity_sha256",1,s.SnapshotIdentitySha256());
		out<<"    </FieldData>\n    <Piece NumberOfPoints=\""<<s.Points().size()<<"\" NumberOfCells=\""<<s.Points().size()<<"\">\n      <PointData>\n";
		auto fields=[&](const char* type,const char* name,unsigned n,auto getter){std::ostringstream values;values<<std::setprecision(17);for(const auto& p:s.Points())getter(values,p);Data(out,type,name,n,values.str());};
		fields("Float64","velocity",3,[](auto&v,const auto&p){for(double x:p.velocity_m_per_s)v<<x<<' ';}); fields("Float64","pressure",1,[](auto&v,const auto&p){v<<p.pressure;}); fields("Float64","vorticity",3,[](auto&v,const auto&p){for(double x:p.vorticity_per_s)v<<x<<' ';}); fields("Float64","q_criterion",1,[](auto&v,const auto&p){v<<p.q_criterion_per_s2;}); fields("Float64","enstrophy_density",1,[](auto&v,const auto&p){v<<p.enstrophy_density_per_s2;}); fields("Float64","speed",1,[](auto&v,const auto&p){v<<p.speed_m_per_s;}); fields("Float64","integration_weight",1,[](auto&v,const auto&p){v<<p.physical_integration_weight_m3;}); fields("Float64","reference_parametric",3,[](auto&v,const auto&p){for(double x:p.parametric)v<<x<<' ';}); fields("UInt64","background_cell_id",1,[](auto&v,const auto&p){v<<p.background_cell_id;}); fields("UInt64","quadrature_ordinal",1,[](auto&v,const auto&p){v<<p.cell_quadrature_ordinal;}); fields("UInt8","cell_kind",1,[](auto&v,const auto&p){v<<static_cast<unsigned>(p.cell_kind);}); fields("UInt8","stagnant",1,[](auto&v,const auto&p){v<<(p.stagnant?1:0);});
		out<<"      </PointData>\n      <Points>\n"; {std::ostringstream values;values<<std::setprecision(17);for(const auto&p:s.Points())for(double x:p.physical_m)values<<x<<' ';Data(out,"Float64","points",3,values.str());} out<<"      </Points>\n      <Cells>\n";
		{std::ostringstream v;for(std::size_t i=0;i<s.Points().size();++i)v<<i<<' ';Data(out,"Int64","connectivity",1,v.str());} {std::ostringstream v;for(std::size_t i=0;i<s.Points().size();++i)v<<(i+1)<<' ';Data(out,"Int64","offsets",1,v.str());} {std::ostringstream v;for(std::size_t i=0;i<s.Points().size();++i)v<<1<<' ';Data(out,"UInt8","types",1,v.str());}
		out<<"      </Cells>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n"; out.close();if(!out)throw std::runtime_error("cannot write moving snapshot VTU");
	}
	static std::string Number(double value) { if(!std::isfinite(value))throw std::invalid_argument("moving snapshot JSON value is nonfinite");std::ostringstream out;out<<std::setprecision(17)<<value;return out.str(); }
	static void WriteMetrics(const std::filesystem::path& path, const MovingImmersedFlowSnapshot& s, const MovingImmersedFlowSnapshotPublicationIdentity& id, const std::optional<MovingImmersedFlowSnapshotTransitionConservation>& c, const std::string& vtu_hash)
	{
		std::ofstream out(path,std::ios::binary|std::ios::trunc);if(!out)throw std::runtime_error("cannot create moving snapshot metrics"); const auto& m=s.Metrics(); const auto q=[&](const std::string&x){return std::string("\"")+EscapeJson(x)+"\"";}; const auto n=[](double x){return Number(x);};
		out<<"{\n  \"schema\": \"MovingImmersedFlowSnapshotPublication\",\n  \"schema_version\": "<<SchemaVersion()<<",\n  \"time_s\": "<<n(s.Request().time_s)<<",\n  \"index\": "<<s.Request().index<<",\n  \"geometry_identity_sha256\": "<<q(id.geometry_identity_sha256)<<",\n  \"publication_identity_sha256\": "<<q(id.publication_identity_sha256)<<",\n  \"layout_identity_sha256\": "<<q(id.layout_identity_sha256)<<",\n  \"state_identity_sha256\": "<<q(id.state_identity_sha256)<<",\n  \"runtime_identity_sha256\": "<<q(id.runtime_identity_sha256)<<",\n  \"content_hash_sha256\": "<<q(s.ContentHashSha256())<<",\n  \"snapshot_identity_sha256\": "<<q(s.SnapshotIdentitySha256())<<",\n  \"vtu_sha256\": "<<q(vtu_hash)<<",\n  \"stagnant_speed_threshold_m_per_s\": "<<n(s.Options().stagnant_speed_threshold_m_per_s)<<",\n  \"units\": {\"length\": \"m\", \"time\": \"s\", \"volume\": \"m3\", \"velocity\": \"m/s\", \"pressure\": \"Pa\"},\n  \"metrics\": {\n";
		auto metric=[&](const char*name,double value,bool last=false){out<<"    \""<<name<<"\": "<<n(value)<<(last?"\n":" ,\n");}; metric("quadrature_volume_m3",m.quadrature_volume_m3);metric("audited_volume_m3",m.audited_volume_m3);metric("enstrophy_integral_m3_per_s2",m.enstrophy_integral_m3_per_s2);metric("mean_enstrophy_per_s2",m.mean_enstrophy_per_s2);metric("mean_q_criterion_per_s2",m.mean_q_criterion_per_s2);metric("q_positive_volume_m3",m.q_positive_volume_m3);metric("q_positive_volume_fraction",m.q_positive_volume_fraction);metric("stagnant_volume_m3",m.stagnant_volume_m3);metric("stagnant_volume_fraction",m.stagnant_volume_fraction);metric("inlet_flow_m3_s",m.inlet_flow_m3_s);metric("outlet_flow_m3_s",m.outlet_flow_m3_s);metric("endpoint_turnover_rate_per_s",m.endpoint_turnover_rate_per_s);out<<"    \"endpoint_turnover_time_s\": ";if(m.endpoint_turnover_time_s)out<<n(*m.endpoint_turnover_time_s);else out<<"null";out<<",\n";metric("well_mixed_replacement_fraction_over_step",m.well_mixed_replacement_fraction_over_step);metric("wall_relative_velocity_squared_area_integral_m4_per_s2",m.wall_relative_velocity_squared_area_integral_m4_per_s2);metric("wall_relative_velocity_rms_m_per_s",m.wall_relative_velocity_rms_m_per_s);metric("wall_relative_velocity_max_m_per_s",m.wall_relative_velocity_max_m_per_s);metric("wall_area_m2",m.wall_area_m2,true);out<<"  },\n  \"port_label_outward_flow_m3_s\": [";for(std::size_t i=0;i<s.Request().port_flows.size();++i){if(i)out<<", ";out<<"{\"label\": "<<s.Request().port_flows[i].label<<", \"outward_flow_m3_s\": "<<n(s.Request().port_flows[i].outward_flow_m3_s)<<"}";}out<<"]";if(c)out<<",\n  \"transition_conservation\": {\"identity_sha256\": "<<q(c->identity_sha256)<<", \"reynolds_residual_m3_s\": "<<n(c->reynolds_residual_m3_s)<<", \"moving_residual_m3_s\": "<<n(c->moving_residual_m3_s)<<", \"volume_change_rate_m3_s\": "<<n(c->volume_change_rate_m3_s)<<"}";out<<"\n}\n";out.close();if(!out)throw std::runtime_error("cannot write moving snapshot metrics");
	}
	static std::string ReadAll(const std::filesystem::path&p){std::ifstream in(p,std::ios::binary);if(!in)throw std::runtime_error("moving snapshot epoch file is missing");std::ostringstream out;out<<in.rdbuf();if(!in)throw std::runtime_error("cannot read moving snapshot epoch file");return out.str();}
	static std::string JsonString(const std::string& text,const std::string& key){const auto needle="\""+key+"\": \"";const auto pos=text.find(needle);if(pos==std::string::npos)throw std::runtime_error("moving snapshot metrics field is missing: "+key);const auto begin=pos+needle.size(),end=text.find('"',begin);if(end==std::string::npos)throw std::runtime_error("moving snapshot metrics string is truncated");return text.substr(begin,end-begin);}
	static std::string XmlField(const std::string& text, const std::string& name)
	{
		const auto marker = "Name=\""+name+"\""; const auto field = text.find(marker);
		if (field == std::string::npos) throw std::runtime_error("moving snapshot VTU field is missing: "+name);
		const auto open = text.find('>', field); const auto close = text.find('<', open == std::string::npos ? field : open+1);
		if (open == std::string::npos || close == std::string::npos) throw std::runtime_error("moving snapshot VTU field is truncated");
		std::string value = text.substr(open+1, close-open-1); const auto first=value.find_first_not_of(" \t\r\n"), last=value.find_last_not_of(" \t\r\n");
		if (first == std::string::npos)
			throw std::runtime_error("moving snapshot VTU field is empty");
		return value.substr(first,last-first+1);
	}
	static std::uint64_t JsonUnsigned(const std::string& text, const std::string& key){const auto needle="\""+key+"\": ";const auto pos=text.find(needle);if(pos==std::string::npos)throw std::runtime_error("moving snapshot metrics integer is missing: "+key);std::size_t used=0;try{const auto v=std::stoull(text.substr(pos+needle.size()),&used);if(!used)return 0;return v;}catch(...){throw std::runtime_error("moving snapshot metrics integer is invalid: "+key);}}
	static std::uint64_t JsonIndex(const std::string& text){return JsonUnsigned(text,"index");}
	static double JsonTime(const std::string& text){const auto needle=std::string("\"time_s\": ");const auto pos=text.find(needle);if(pos==std::string::npos)throw std::runtime_error("moving snapshot metrics time is missing");std::size_t used=0;try{const auto v=std::stod(text.substr(pos+needle.size()),&used);if(!used||!std::isfinite(v))throw std::runtime_error("x");return v;}catch(...){throw std::runtime_error("moving snapshot metrics time is invalid");}}
	static Epoch ValidateEpoch(const std::filesystem::path& directory,const std::filesystem::path& stem,std::uint64_t expected_index)
	{
		if (!std::filesystem::is_directory(directory))
			throw std::runtime_error("moving snapshot epoch directory is invalid");
		const auto vtu=directory/"fields.vtu", metrics=directory/"metrics.json";
		if (!std::filesystem::is_regular_file(vtu) || !std::filesystem::is_regular_file(metrics))
			throw std::runtime_error("moving snapshot epoch is incomplete");
		const auto json=ReadAll(metrics), xml=ReadAll(vtu);
		if (JsonString(json,"schema")!="MovingImmersedFlowSnapshotPublication"
			|| JsonUnsigned(json,"schema_version")!=SchemaVersion()
			|| XmlField(xml,"moving_snapshot_schema")!="MovingImmersedFlowSnapshotPublication"
			|| XmlField(xml,"moving_snapshot_schema_version")!=std::to_string(SchemaVersion()))
			throw std::runtime_error("moving snapshot epoch schema is invalid");
		const auto index=JsonIndex(json); std::uint64_t xml_index=0;
		try { xml_index=std::stoull(XmlField(xml,"index")); }
		catch (...) { throw std::runtime_error("moving snapshot VTU index is invalid"); }
		if (index!=expected_index || xml_index!=expected_index)
			throw std::runtime_error("moving snapshot epoch index is inconsistent");
		const auto time=JsonTime(json);
		if (!std::isfinite(time) || JsonString(json,"snapshot_identity_sha256")!=XmlField(xml,"snapshot_identity_sha256")
			|| JsonString(json,"content_hash_sha256")!=XmlField(xml,"content_hash_sha256")
			|| JsonString(json,"geometry_identity_sha256")!=XmlField(xml,"geometry_identity_sha256")
			|| JsonString(json,"publication_identity_sha256")!=XmlField(xml,"publication_identity_sha256")
			|| JsonString(json,"layout_identity_sha256")!=XmlField(xml,"layout_identity_sha256")
			|| JsonString(json,"state_identity_sha256")!=XmlField(xml,"state_identity_sha256")
			|| JsonString(json,"runtime_identity_sha256")!=XmlField(xml,"runtime_identity_sha256"))
			throw std::runtime_error("moving snapshot epoch identity is inconsistent");
		if (XmlField(xml,"time_s") != Number(time))
			throw std::runtime_error("moving snapshot epoch time is inconsistent");
		for (const auto& hash : {JsonString(json,"snapshot_identity_sha256"), JsonString(json,"content_hash_sha256"),
			JsonString(json,"geometry_identity_sha256"), JsonString(json,"publication_identity_sha256"),
			JsonString(json,"layout_identity_sha256"), JsonString(json,"state_identity_sha256"), JsonString(json,"runtime_identity_sha256")})
			if (!ValidHash(hash)) throw std::runtime_error("moving snapshot epoch hash is invalid");
		const auto vtu_hash=FileHash(vtu);
		if (JsonString(json,"vtu_sha256")!=vtu_hash)
			throw std::runtime_error("moving snapshot VTU hash is invalid");
		if (xml.find("NumberOfPoints=\"")==std::string::npos || xml.find("Name=\"types\"")==std::string::npos)
			throw std::runtime_error("moving snapshot VTU is incomplete");
		(void)stem;
		return {directory,index,time,JsonString(json,"snapshot_identity_sha256"),JsonString(json,"content_hash_sha256"),vtu_hash,FileHash(metrics)};
	}
	static MovingImmersedFlowSnapshotPublication Result(const std::filesystem::path& dir,const std::filesystem::path&pvd,const std::string&vtu,const std::string&metrics){return {dir,dir/"fields.vtu",dir/"metrics.json",pvd,vtu,metrics};}
	static std::string Nonce(){std::ostringstream out;out<<std::chrono::steady_clock::now().time_since_epoch().count()<<'.'<<std::random_device{}();return out.str();}
	static std::filesystem::path TemporaryDirectory(const std::filesystem::path& final){return final.string()+".tmp."+Nonce();}
#ifdef IGA_MOVING_IMMERSED_FLOW_SNAPSHOT_PUBLISHER_TESTING
	static MovingImmersedFlowSnapshotPublicationFault& TestFault(){static MovingImmersedFlowSnapshotPublicationFault value=MovingImmersedFlowSnapshotPublicationFault::None;return value;}
	static void ThrowFault(MovingImmersedFlowSnapshotPublicationFault point){if(TestFault()==point){TestFault()=MovingImmersedFlowSnapshotPublicationFault::None;throw std::runtime_error("injected moving snapshot publication failure");}}
#endif
};

} // namespace iga

#endif

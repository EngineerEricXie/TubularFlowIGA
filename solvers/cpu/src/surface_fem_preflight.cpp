#include "SurfaceReaders.hpp"
#include "Sha256.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Arguments
{
	std::string input;
	std::string output;
	std::string manifest;
	std::string boundary_array = "boundary_id";
	iga::SurfaceValidationOptions validation;
};

double ParseDouble(const std::string& text,const std::string& name)
{
	std::size_t consumed=0;
	const double value=std::stod(text,&consumed);
	if(consumed!=text.size()) throw std::invalid_argument(name+" is not a number");
	return value;
}

std::int64_t ParseInteger(const std::string& text,const std::string& name)
{
	std::size_t consumed=0;
	const auto value=std::stoll(text,&consumed);
	if(consumed!=text.size()) throw std::invalid_argument(name+" is not an integer");
	return value;
}

Arguments ParseArguments(int argc,char** argv)
{
	if(argc<3) throw std::invalid_argument(
		"usage: surface_fem_preflight INPUT.{vtp,stl} OUTPUT.msh [--manifest FILE.json] "
		"[--boundary-array NAME] [--length-scale-to-m VALUE] [--weld-tolerance-m VALUE] "
		"[--default-boundary-id ID] [--max-triangles COUNT]");
	Arguments result;
	result.input=argv[1];result.output=argv[2];
	for(int i=3;i<argc;++i) {
		const std::string option=argv[i];
		if(i+1>=argc) throw std::invalid_argument("missing value for "+option);
		const std::string value=argv[++i];
		if(option=="--manifest") result.manifest=value;
		else if(option=="--boundary-array") result.boundary_array=value;
		else if(option=="--length-scale-to-m") result.validation.length_scale_to_m=ParseDouble(value,option);
		else if(option=="--weld-tolerance-m") result.validation.weld_tolerance_m=ParseDouble(value,option);
		else if(option=="--default-boundary-id") result.validation.default_boundary_id=ParseInteger(value,option);
		else if(option=="--max-triangles") result.validation.max_triangles=ParseInteger(value,option);
		else throw std::invalid_argument("unknown option "+option);
	}
	if(result.manifest.empty()) result.manifest=result.output+".json";
	return result;
}

std::string LowerExtension(const std::string& path)
{
	std::string result=std::filesystem::path(path).extension().string();
	std::transform(result.begin(),result.end(),result.begin(),[](unsigned char value){return static_cast<char>(std::tolower(value));});
	return result;
}

std::string FileSha256(const std::string& path)
{
	std::ifstream input(path,std::ios::binary);
	if(!input) throw std::runtime_error("cannot open input surface: "+path);
	iga::Sha256 hash;
	std::vector<char> buffer(64*1024);
	while(input) {
		input.read(buffer.data(),static_cast<std::streamsize>(buffer.size()));
		const auto count=input.gcount();
		if(count>0) hash.Append(buffer.data(),static_cast<std::size_t>(count));
	}
	if(!input.eof()) throw std::runtime_error("cannot read input surface: "+path);
	return hash.Hex();
}

std::string JsonEscape(const std::string& value)
{
	std::ostringstream result;
	for(const unsigned char character:value) {
		if(character=='"'||character=='\\') result<<'\\'<<character;
		else if(character=='\n') result<<"\\n";
		else if(character<0x20) result<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')
			<<static_cast<unsigned>(character)<<std::dec;
		else result<<character;
	}
	return result.str();
}

std::map<std::uint32_t,std::size_t> LabelCounts(const iga::ClosedTriangulatedSurface& surface)
{
	std::map<std::uint32_t,std::size_t> result;
	for(const auto& triangle:surface.Triangles()) ++result[triangle.boundary_id];
	return result;
}

std::map<std::uint32_t,double> LabelAreas(const iga::ClosedTriangulatedSurface& surface)
{
	std::map<std::uint32_t,double> result;
	for(const auto& triangle:surface.Triangles()) result[triangle.boundary_id]+=triangle.area_m2;
	return result;
}

struct ShapeMetrics
{
	double bounds_diagonal_m=0.0;
	double volume_equivalent_radius_m=0.0;
	double area_equivalent_radius_m=0.0;
	double integrated_absolute_mean_curvature_m=0.0;
	double area_average_absolute_mean_curvature_per_m=0.0;
	double maximum_edge_dihedral_rad=0.0;
};

ShapeMetrics ComputeShapeMetrics(const iga::ClosedTriangulatedSurface& surface)
{
	ShapeMetrics result;
	const auto& diagnostics=surface.Diagnostics();
	double diagonal_squared=0.0;
	for(std::size_t axis=0;axis<3;++axis) {
		const double extent=diagnostics.bounds.maximum[axis]-diagnostics.bounds.minimum[axis];
		diagonal_squared+=extent*extent;
	}
	result.bounds_diagonal_m=std::sqrt(diagonal_squared);
	const double pi=std::acos(-1.0);
	result.volume_equivalent_radius_m=std::cbrt(3.0*diagnostics.volume_m3/(4.0*pi));
	result.area_equivalent_radius_m=std::sqrt(diagnostics.area_m2/(4.0*pi));
	using Edge=std::array<std::uint32_t,2>;
	std::map<Edge,std::vector<std::array<double,3>>> edge_normals;
	for(const auto& triangle:surface.Triangles())
		for(const auto edge:std::array<std::array<int,2>,3>{{{{0,1}},{{1,2}},{{2,0}}}}) {
			Edge key{{triangle.indices[edge[0]],triangle.indices[edge[1]]}};
			if(key[1]<key[0]) std::swap(key[0],key[1]);
			edge_normals[key].push_back(triangle.outward_unit_normal);
		}
	for(const auto& entry:edge_normals) {
		if(entry.second.size()!=2) throw std::runtime_error("curvature edge is not manifold");
		double cosine=0.0;
		for(std::size_t axis=0;axis<3;++axis) cosine+=entry.second[0][axis]*entry.second[1][axis];
		cosine=std::max(-1.0,std::min(1.0,cosine));
		const double angle=std::acos(cosine);
		const auto& first=surface.Vertices()[entry.first[0]];
		const auto& second=surface.Vertices()[entry.first[1]];
		double length_squared=0.0;
		for(std::size_t axis=0;axis<3;++axis) {
			const double delta=second[axis]-first[axis];
			length_squared+=delta*delta;
		}
		result.integrated_absolute_mean_curvature_m+=0.5*std::sqrt(length_squared)*angle;
		result.maximum_edge_dihedral_rad=std::max(result.maximum_edge_dihedral_rad,angle);
	}
	result.area_average_absolute_mean_curvature_per_m=
		result.integrated_absolute_mean_curvature_m/diagnostics.area_m2;
	return result;
}

void WriteMsh(const iga::ClosedTriangulatedSurface& surface,const std::string& path)
{
	const auto counts=LabelCounts(surface);
	std::map<std::uint32_t,int> tags;
	int next=1;for(const auto& entry:counts) tags.emplace(entry.first,next++);
	std::ofstream output(path);
	if(!output) throw std::runtime_error("cannot write canonical Gmsh surface: "+path);
	output<<std::setprecision(17)<<"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n"
		<<tags.size()<<'\n';
	for(const auto& entry:tags) output<<"2 "<<entry.second<<" \"boundary_label_"<<entry.first<<"\"\n";
	output<<"$EndPhysicalNames\n$Nodes\n"<<surface.Vertices().size()<<'\n';
	for(std::size_t i=0;i<surface.Vertices().size();++i) {
		const auto& point=surface.Vertices()[i];
		output<<i+1<<' '<<point[0]<<' '<<point[1]<<' '<<point[2]<<'\n';
	}
	output<<"$EndNodes\n$Elements\n"<<surface.Triangles().size()<<'\n';
	for(std::size_t i=0;i<surface.Triangles().size();++i) {
		const auto& triangle=surface.Triangles()[i];
		const int tag=tags.at(triangle.boundary_id);
		output<<i+1<<" 2 2 "<<tag<<' '<<tag;
		for(const auto node:triangle.indices) output<<' '<<node+1;
		output<<'\n';
	}
	output<<"$EndElements\n";
	output.close();
	if(!output) throw std::runtime_error("cannot finalize canonical Gmsh surface: "+path);
}

void WriteManifest(const Arguments& arguments,const iga::ClosedTriangulatedSurface& surface,
	const std::string& source_sha256,const std::string& output_sha256)
{
	const auto counts=LabelCounts(surface);
	const auto areas=LabelAreas(surface);
	const auto& diagnostics=surface.Diagnostics();
	const auto shape=ComputeShapeMetrics(surface);
	std::ofstream output(arguments.manifest);
	if(!output) throw std::runtime_error("cannot write surface preflight manifest: "+arguments.manifest);
	output<<std::setprecision(17)<<"{\n  \"schema_version\": 1,\n  \"kind\": \"canonical_closed_surface\",\n"
		<<"  \"source\": {\"path\": \""<<JsonEscape(arguments.input)<<"\", \"sha256\": \""
		<<source_sha256<<"\"},\n  \"canonical_sha256\": \""<<surface.CanonicalSha256()<<"\",\n"
		<<"  \"canonical_msh_sha256\": \""<<output_sha256<<"\",\n"
		<<"  \"boundary_array\": \""<<JsonEscape(arguments.boundary_array)<<"\",\n"
		<<"  \"default_boundary_id\": "<<arguments.validation.default_boundary_id<<",\n"
		<<"  \"max_triangles\": "<<arguments.validation.max_triangles<<",\n"
		<<"  \"applied_length_scale_to_m\": "<<diagnostics.applied_length_scale_to_m<<",\n"
		<<"  \"applied_weld_tolerance_m\": "<<diagnostics.applied_weld_tolerance_m<<",\n"
		<<"  \"input_vertices\": "<<diagnostics.input_vertex_count<<",\n"
		<<"  \"canonical_vertices\": "<<diagnostics.vertex_count<<",\n"
		<<"  \"welded_vertices\": "<<diagnostics.welded_vertex_count<<",\n"
		<<"  \"triangles\": "<<diagnostics.triangle_count<<",\n"
		<<"  \"area_m2\": "<<diagnostics.area_m2<<",\n  \"volume_m3\": "<<diagnostics.volume_m3<<",\n"
		<<"  \"bounds_m\": {\"minimum\": ["<<diagnostics.bounds.minimum[0]<<", "
		<<diagnostics.bounds.minimum[1]<<", "<<diagnostics.bounds.minimum[2]
		<<"], \"maximum\": ["<<diagnostics.bounds.maximum[0]<<", "
		<<diagnostics.bounds.maximum[1]<<", "<<diagnostics.bounds.maximum[2]
		<<"], \"diagonal_m\": "<<shape.bounds_diagonal_m<<"},\n"
		<<"  \"equivalent_radii_m\": {\"enclosed_volume_sphere\": "
		<<shape.volume_equivalent_radius_m<<", \"surface_area_sphere\": "
		<<shape.area_equivalent_radius_m<<"},\n"
		<<"  \"discrete_curvature\": {\"definition\": \"one_half_sum_edge_length_times_absolute_dihedral\", "
		<<"\"integrated_absolute_mean_curvature_m\": "
		<<shape.integrated_absolute_mean_curvature_m
		<<", \"area_average_absolute_mean_curvature_per_m\": "
		<<shape.area_average_absolute_mean_curvature_per_m
		<<", \"maximum_edge_dihedral_rad\": "<<shape.maximum_edge_dihedral_rad<<"},\n"
		<<"  \"minimum_edge_length_m\": "<<diagnostics.minimum_edge_length_m<<",\n"
		<<"  \"minimum_triangle_area_m2\": "<<diagnostics.minimum_triangle_area_m2<<",\n"
		<<"  \"flipped_inward_shell\": "<<(diagnostics.flipped_inward_shell?"true":"false")<<",\n"
		<<"  \"boundary_labels\": {\n";
	std::size_t index=0;
	for(const auto& entry:counts) {
		output<<"    \""<<entry.first<<"\": {\"triangles\": "<<entry.second
			<<", \"area_m2\": "<<areas.at(entry.first)<<"}"<<(++index==counts.size()?"\n":",\n");
	}
	output<<"  },\n  \"geometry_change_policy\": {\"implicit_repair\": false, "
		<<"\"implicit_remesh\": false, \"welding_is_explicitly_parameterized\": true},\n"
		<<"  \"checks\": {\"closed\": true, \"oriented\": true, \"manifold\": true, "
		<<"\"connected\": true, \"self_intersection_free\": true, \"positive_volume\": true}\n}\n";
	output.close();
	if(!output) throw std::runtime_error("cannot finalize surface preflight manifest: "+arguments.manifest);
}

} // namespace

int main(int argc,char** argv)
{
	try {
		const auto arguments=ParseArguments(argc,argv);
		const auto extension=LowerExtension(arguments.input);
		iga::ClosedTriangulatedSurface surface=extension==".vtp"
			? iga::SurfaceReaders::ReadVtpPath(arguments.input,arguments.validation,arguments.boundary_array)
			: extension==".stl"?iga::SurfaceReaders::ReadStlPath(arguments.input,arguments.validation)
			: throw std::invalid_argument("surface input extension must be .vtp or .stl");
		const auto source_sha256=FileSha256(arguments.input);
		WriteMsh(surface,arguments.output);
		WriteManifest(arguments,surface,source_sha256,FileSha256(arguments.output));
		std::cout<<"surface_fem_preflight: PASS vertices="<<surface.Vertices().size()
			<<" triangles="<<surface.Triangles().size()<<" canonical_sha256="
			<<surface.CanonicalSha256()<<'\n';
		return 0;
	} catch(const std::exception& error) {
		std::cerr<<"surface_fem_preflight: ERROR: "<<error.what()<<'\n';
		return 2;
	}
}

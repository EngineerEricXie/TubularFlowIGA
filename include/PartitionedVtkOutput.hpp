#ifndef IGA_PARTITIONED_VTK_OUTPUT_HPP
#define IGA_PARTITIONED_VTK_OUTPUT_HPP

#include "VtkOutput.hpp"
#include <cmath>
#include <set>
#include <limits>

namespace iga {

struct VtkArraySchema {
	std::string name;
	int components=1;
};

// A piece owns its cells. Shared points may appear in several pieces with the
// same global id; ids must be assigned by the caller, never by coordinates.
struct VtkPartition {
	LegacyUnstructuredGrid grid;
	std::vector<std::int64_t> point_ids,cell_ids;
	std::vector<VtkPointArray> point_arrays,cell_arrays;
};

namespace partitioned_vtk_detail {
inline void ValidateName(const std::string& name)
{
	if(name.empty()||name=="GlobalPointIds"||name=="GlobalCellIds")
		throw std::invalid_argument("invalid or reserved partitioned VTU array name");
	for(unsigned char c:name)if(c<32||c==127)
		throw std::invalid_argument("control character in partitioned VTU array name");
}
inline void ValidateSchema(const std::vector<VtkArraySchema>& schema)
{
	std::set<std::string> names;
	for(const auto& array:schema) {
		ValidateName(array.name);
		if(array.components<1||!names.insert(array.name).second)
			throw std::invalid_argument("invalid partitioned VTU array schema");
	}
}
inline std::vector<VtkArraySchema> Schema(const std::vector<VtkPointArray>& arrays,std::size_t count)
{
	std::vector<VtkArraySchema> schema;
	for(const auto& array:arrays) {
		if(array.components<1||count>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(array.components)
			||array.values.size()!=count*static_cast<std::size_t>(array.components))
			throw std::invalid_argument("partitioned VTU tuple count mismatch");
		for(double value:array.values)if(!std::isfinite(value))throw std::invalid_argument("nonfinite partitioned VTU field");
		schema.push_back({array.name,array.components});
	}
	ValidateSchema(schema);return schema;
}
template<class T>
inline void Array(std::ostream& output,const char* type,const std::string& name,int components,const std::vector<T>& values)
{
	output<<"<DataArray type=\""<<type<<"\" Name=\""<<EscapeVtkXml(name)<<"\" NumberOfComponents=\""<<components<<"\" format=\"ascii\">\n";
	for(const auto& value:values)output<<value<<' ';
	output<<"\n</DataArray>\n";
}
inline void Fields(std::ostream& output,const char* tag,const char* id_name,const std::vector<std::int64_t>& ids,const std::vector<VtkPointArray>& arrays)
{
	output<<'<'<<tag<<" GlobalIds=\""<<id_name<<"\">\n";
	Array(output,"Int64",id_name,1,ids);
	for(const auto& array:arrays)Array(output,"Float64",array.name,array.components,array.values);
	output<<"</"<<tag<<">\n";
}
inline void ParallelFields(std::ostream& output,const char* tag,const char* id_name,const std::vector<VtkArraySchema>& schema)
{
	output<<'<'<<tag<<" GlobalIds=\""<<id_name<<"\">\n<PDataArray type=\"Int64\" Name=\""<<id_name<<"\" NumberOfComponents=\"1\"/>\n";
	for(const auto& array:schema)output<<"<PDataArray type=\"Float64\" Name=\""<<EscapeVtkXml(array.name)<<"\" NumberOfComponents=\""<<array.components<<"\"/>\n";
	output<<"</"<<tag<<">\n";
}
}

inline void ValidateVtkPartition(const VtkPartition& piece,double physical_time)
{
	const auto& grid=piece.grid;const auto points=grid.points.size()/3,cells=grid.offsets.size();
	if(!std::isfinite(physical_time)||grid.points.size()%3||grid.types.size()!=cells
		||piece.point_ids.size()!=points||piece.cell_ids.size()!=cells)
		throw std::invalid_argument("invalid partitioned VTU geometry dimensions or time");
	for(double value:grid.points)if(!std::isfinite(value))throw std::invalid_argument("nonfinite partitioned VTU coordinate");
	std::int64_t previous=0;
	for(std::size_t cell=0;cell<cells;++cell) {
		const auto end=grid.offsets[cell];
		if(end<=previous||static_cast<std::uint64_t>(end)>grid.connectivity.size()||grid.types[cell]==0||grid.types[cell]>255)
			throw std::invalid_argument("invalid partitioned VTU cell");
		previous=end;
	}
	if(static_cast<std::uint64_t>(previous)!=grid.connectivity.size())throw std::invalid_argument("partitioned VTU trailing connectivity");
	for(auto node:grid.connectivity)if(node<0||static_cast<std::uint64_t>(node)>=points)throw std::invalid_argument("partitioned VTU connectivity out of range");
	for(const auto* ids:{&piece.point_ids,&piece.cell_ids}) {
		std::set<std::int64_t> unique;
		for(auto id:*ids)if(id<0||!unique.insert(id).second)throw std::invalid_argument("invalid or duplicate partitioned VTU global id");
	}
	partitioned_vtk_detail::Schema(piece.point_arrays,points);
	partitioned_vtk_detail::Schema(piece.cell_arrays,cells);
}

inline void WriteVtuPartition(const std::filesystem::path& path,const VtkPartition& piece,double physical_time)
{
	ValidateVtkPartition(piece,physical_time);
	if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
	std::ofstream output(path);if(!output)throw std::runtime_error("cannot create VTU partition: "+path.string());
	output<<std::setprecision(17)<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid>\n"
		<<"<FieldData><DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"<<physical_time<<"</DataArray></FieldData>\n"
		<<"<Piece NumberOfPoints=\""<<piece.point_ids.size()<<"\" NumberOfCells=\""<<piece.cell_ids.size()<<"\">\n";
	partitioned_vtk_detail::Fields(output,"PointData","GlobalPointIds",piece.point_ids,piece.point_arrays);
	partitioned_vtk_detail::Fields(output,"CellData","GlobalCellIds",piece.cell_ids,piece.cell_arrays);
	output<<"<Points>\n";partitioned_vtk_detail::Array(output,"Float64","Points",3,piece.grid.points);output<<"</Points><Cells>\n";
	partitioned_vtk_detail::Array(output,"Int64","connectivity",1,piece.grid.connectivity);
	partitioned_vtk_detail::Array(output,"Int64","offsets",1,piece.grid.offsets);
	partitioned_vtk_detail::Array(output,"UInt8","types",1,piece.grid.types);
	output<<"</Cells></Piece></UnstructuredGrid></VTKFile>\n";
	output.close();if(!output)throw std::runtime_error("cannot write VTU partition: "+path.string());
}

// The collective caller publishes this index only after every piece succeeds.
// This format primitive does not coordinate MPI or validate files on other ranks.
inline void WritePvtu(const std::filesystem::path& path,const std::vector<std::filesystem::path>& pieces,
	const std::vector<VtkArraySchema>& point_schema,const std::vector<VtkArraySchema>& cell_schema={})
{
	partitioned_vtk_detail::ValidateSchema(point_schema);partitioned_vtk_detail::ValidateSchema(cell_schema);
	if(pieces.empty())throw std::invalid_argument("PVTU requires at least one piece");
	std::set<std::filesystem::path> unique;
	for(const auto& piece:pieces) {
		if(piece.empty()||piece.is_absolute()||piece!=piece.filename()||!unique.insert(piece).second)
			throw std::invalid_argument("PVTU piece names must be distinct local filenames");
		for(unsigned char c:piece.string())if(c<32||c==127)throw std::invalid_argument("control character in PVTU filename");
	}
	if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
	std::ofstream output(path);if(!output)throw std::runtime_error("cannot create PVTU index: "+path.string());
	output<<"<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><PUnstructuredGrid GhostLevel=\"0\">\n";
	partitioned_vtk_detail::ParallelFields(output,"PPointData","GlobalPointIds",point_schema);
	partitioned_vtk_detail::ParallelFields(output,"PCellData","GlobalCellIds",cell_schema);
	output<<"<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n";
	for(const auto& piece:pieces)output<<"<Piece Source=\""<<EscapeVtkXml(piece.string())<<"\"/>\n";
	output<<"</PUnstructuredGrid></VTKFile>\n";
	output.close();if(!output)throw std::runtime_error("cannot write PVTU index: "+path.string());
}

} // namespace iga
#endif

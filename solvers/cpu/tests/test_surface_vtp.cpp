#include "SurfaceReaders.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

template <class Function>
void RequireRejected(Function&& function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	assert(rejected);
}

using Point = std::array<double, 3>;
using Face = std::array<int, 3>;

const std::array<Point, 4> kTetraPoints{{{{0,0,0}}, {{0,1,0}}, {{1,0,0}}, {{0,0,1}}}};
const std::array<Face, 4> kTetraFaces{{{{0,1,2}}, {{0,2,3}}, {{0,3,1}}, {{2,1,3}}}};
const std::array<Point, 8> kCubePoints{{{{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}}, {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}}};
const std::array<Face, 12> kCubeFaces{{{{0,2,1}}, {{0,3,2}}, {{4,5,6}}, {{4,6,7}}, {{0,1,5}}, {{0,5,4}}, {{1,2,6}}, {{1,6,5}}, {{2,3,7}}, {{2,7,6}}, {{3,0,4}}, {{3,4,7}}}};

template <std::size_t Points, std::size_t Faces>
std::string Vtp(const std::array<Point, Points>& points, const std::array<Face, Faces>& faces,
	const char* point_type = "Float64", const char* index_type = "Int32", const char* label_type = "UInt32",
	const std::vector<unsigned>& labels = {}, const std::string& label_name = "boundary_id")
{
	std::ostringstream result;
	result << "<?xml version=\"1.0\"?>\n<!-- strict ASCII fixture -->\n"
		<< "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
		<< "<PolyData><Piece NumberOfPoints=\"" << Points << "\" NumberOfPolys=\"" << Faces << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\">\n"
		<< "<Points><DataArray type=\"" << point_type << "\" NumberOfComponents=\"3\" format=\"ascii\">";
	for (const auto& point : points) result << point[0] << ' ' << point[1] << ' ' << point[2] << ' ';
	result << "</DataArray></Points><Polys><DataArray type=\"" << index_type << "\" Name=\"connectivity\" format=\"ascii\">";
	for (const auto& face : faces) result << face[0] << ' ' << face[1] << ' ' << face[2] << ' ';
	result << "</DataArray><DataArray type=\"" << index_type << "\" Name=\"offsets\" format=\"ascii\">";
	for (std::size_t face = 1; face <= Faces; ++face) result << 3*face << ' ';
	result << "</DataArray></Polys>";
	if (!labels.empty()) {
		assert(labels.size() == Faces);
		result << "<CellData Scalars=\"" << label_name << "\"><DataArray type=\"" << label_type << "\" Name=\"" << label_name << "\" format=\"ascii\">";
		for (const auto label : labels) result << label << ' ';
		result << "</DataArray></CellData>";
	}
	return result.str()+"</Piece></PolyData></VTKFile>\n";
}

std::string Replace(std::string text, const std::string& before, const std::string& after)
{
	const auto at = text.find(before); assert(at != std::string::npos); text.replace(at, before.size(), after); return text;
}

} // namespace

int main()
{
	const auto tetra64 = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces), {}, "tetra64");
	const auto tetra32 = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces, "Float32", "Int64"), {}, "tetra32");
	assert(tetra64.CanonicalSha256() == tetra32.CanonicalSha256());
	assert(std::fabs(tetra64.Diagnostics().volume_m3-1.0/6.0) < 1.0e-12);
	auto decimal_points = kTetraPoints; decimal_points[2][0] = 0.1;
	const auto decimal32 = iga::SurfaceReaders::ReadVtpText(Vtp(decimal_points, kTetraFaces, "Float32"));
	const auto decimal64 = iga::SurfaceReaders::ReadVtpText(Vtp(decimal_points, kTetraFaces, "Float64"));
	const double parsed_float32 = static_cast<double>(static_cast<float>(0.1)); bool found_float32 = false;
	for (const auto& point : decimal32.Vertices()) if (point[0] == parsed_float32) found_float32 = true;
	assert(found_float32 && decimal32.CanonicalSha256() != decimal64.CanonicalSha256());
	iga::SurfaceValidationOptions default_label; default_label.default_boundary_id = 5;
	assert(iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces), default_label).Diagnostics().boundary_id_histogram.at(5) == 4U);
	const auto cube = iga::SurfaceReaders::ReadVtpText(Vtp(kCubePoints, kCubeFaces, "Float64", "Int64", "UInt64", std::vector<unsigned>(12, 7)), {}, "cube");
	assert(std::fabs(cube.Diagnostics().area_m2-6.0) < 1.0e-12);
	assert(std::fabs(cube.Diagnostics().volume_m3-1.0) < 1.0e-12);
	const auto labelled_a = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", {1, 1, 1, 1}));
	const auto labelled_b = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "Int64", {2, 1, 1, 1}));
	assert(labelled_a.CanonicalSha256() != labelled_b.CanonicalSha256());
	assert(labelled_b.Diagnostics().boundary_id_histogram.at(2) == 1U);
	const auto named = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "Int32", {9, 9, 9, 9}, "wall_id"), {}, "named", "wall_id");
	assert(named.Diagnostics().boundary_id_histogram.at(9) == 4U);
	iga::SurfaceValidationOptions millimetres; millimetres.length_scale_to_m = 0.001;
	auto millimetre_points = kTetraPoints; for (auto& point : millimetre_points) for (auto& coordinate : point) coordinate *= 1000.0;
	assert(iga::SurfaceReaders::ReadVtpText(Vtp(millimetre_points, kTetraFaces), millimetres).CanonicalSha256() == tetra64.CanonicalSha256());

	const auto valid = Vtp(kTetraPoints, kTetraFaces);
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "type=\"PolyData\"", "type=\"UnstructuredGrid\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "byte_order=\"LittleEndian\"", "byte_order=\"BigEndian\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "header_type=\"UInt64\"", "compressor=\"vtkZLibDataCompressor\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"binary\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"appended\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"base64\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "0 1 0 1 0 0", "0 1 0 nan 0 0")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "0 1 0 1 0 0", "0 1 0 1x 0 0")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "NumberOfPolys=\"4\"", "NumberOfPolys=\"3\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "NumberOfPolys=\"4\"", "NumberOfPolys=\"3074457345618258603\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "0 1 2 0 2 3", "0 1 9 0 2 3")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "3 6 9 12", "3 7 9 12")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "NumberOfVerts=\"0\"", "NumberOfVerts=\"1\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "<Polys>", "<Polys><DataArray type=\"Int32\" Name=\"unknown\" format=\"ascii\">0</DataArray>")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "</PolyData>", "<Piece NumberOfPoints=\"0\" NumberOfPolys=\"0\"/></PolyData>")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText("<!DOCTYPE VTKFile>"+valid); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "<PolyData>", "<?bad instruction?><PolyData>")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "type=\"PolyData\"", "type=\"PolyData\" type=\"PolyData\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "type=\"PolyData\"", "type=\"Poly&amp;Data\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "</VTKFile>", "&amp;</VTKFile>")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "</Polys>", "</Points>")); });
	std::string deeply_nested = "<VTKFile type=\"PolyData\">";
	for (int level = 0; level < 40; ++level) deeply_nested += "<node>";
	for (int level = 0; level < 40; ++level) deeply_nested += "</node>";
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(deeply_nested+"</VTKFile>"); });
	const auto labels = Vtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", {1, 1, 1, 1});
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(labels, "1 1 1 1", "1 1 1 -1")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(labels, "UInt32", "Float64")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(labels, "boundary_id", "other")); });
	iga::SurfaceValidationOptions maximum; maximum.max_triangles = 3; RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(valid, maximum); });
	iga::SurfaceValidationOptions negative_maximum; negative_maximum.max_triangles = -1; RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(valid, negative_maximum); });
	iga::SurfaceValidationOptions empty_cell_default; empty_cell_default.default_boundary_id = 8;
	assert(iga::SurfaceReaders::ReadVtpText(Replace(valid, "</Polys>", "</Polys><CellData/>"), empty_cell_default).Diagnostics().boundary_id_histogram.at(8) == 4U);
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "</Polys>", "</Polys><CellData><DataArray type=\"Int32\" Name=\"unknown\" format=\"ascii\">1 1 1 1</DataArray></CellData>")); });

	const std::string path = "/tmp/tubularflowiga_surface_vtp_test.vtp";
	{ std::ofstream output(path, std::ios::binary); output << valid; }
	assert(iga::SurfaceReaders::ReadVtpPath(path).CanonicalSha256() == tetra64.CanonicalSha256());
	assert(std::remove(path.c_str()) == 0);
	RequireRejected([] { iga::SurfaceReaders::ReadVtpPath("/tmp/tubularflowiga_missing_surface.vtp"); });
	std::cout << "surface VTP tests passed\n";
}

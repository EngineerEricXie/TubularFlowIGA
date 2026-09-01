#include "SurfaceReaders.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
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

void AppendLe(std::vector<std::uint8_t>& bytes, std::uint64_t value, std::size_t size)
{
	for (std::size_t shift = 0; shift < size; ++shift) bytes.push_back(static_cast<std::uint8_t>(value>>(8U*shift)));
}

void AppendFloating(std::vector<std::uint8_t>& bytes, double value, const std::string& type)
{
	if (type == "Float32") { const float converted = static_cast<float>(value); std::uint32_t bits = 0; std::memcpy(&bits, &converted, sizeof(bits)); AppendLe(bytes, bits, 4); }
	else { std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits)); AppendLe(bytes, bits, 8); }
}

void AppendInteger(std::vector<std::uint8_t>& bytes, std::int64_t value, const std::string& type)
{
	if (type == "Int32" || type == "UInt32") AppendLe(bytes, static_cast<std::uint32_t>(value), 4);
	else AppendLe(bytes, static_cast<std::uint64_t>(value), 8);
}

std::vector<std::uint8_t> VtkBlock(const std::vector<std::uint8_t>& payload, bool header64)
{
	std::vector<std::uint8_t> result; AppendLe(result, payload.size(), header64 ? 8 : 4); result.insert(result.end(), payload.begin(), payload.end()); return result;
}

std::string Base64(const std::vector<std::uint8_t>& bytes)
{
	static const char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	std::string result;
	for (std::size_t at = 0; at < bytes.size(); at += 3) {
		const auto first = bytes[at]; const auto second = at+1 < bytes.size() ? bytes[at+1] : 0U; const auto third = at+2 < bytes.size() ? bytes[at+2] : 0U;
		result.push_back(alphabet[first>>2]); result.push_back(alphabet[((first&3U)<<4)|(second>>4)]);
		result.push_back(at+1 < bytes.size() ? alphabet[((second&15U)<<2)|(third>>6)] : '='); result.push_back(at+2 < bytes.size() ? alphabet[third&63U] : '=');
	}
	return result;
}

template <std::size_t Points, std::size_t Faces>
std::string EncodedVtp(const std::array<Point, Points>& points, const std::array<Face, Faces>& faces,
	const std::string& point_type, const std::string& index_type, const std::string& label_type,
	const std::vector<unsigned>& labels, bool appended, bool header64, bool leading_appended_gap = false)
{
	std::vector<std::uint8_t> point_bytes, connectivity_bytes, offsets_bytes, label_bytes;
	for (const auto& point : points) for (const auto coordinate : point) AppendFloating(point_bytes, coordinate, point_type);
	for (const auto& face : faces) for (const auto index : face) AppendInteger(connectivity_bytes, index, index_type);
	for (std::size_t face = 1; face <= Faces; ++face) AppendInteger(offsets_bytes, static_cast<std::int64_t>(3*face), index_type);
	for (const auto label : labels) AppendInteger(label_bytes, label, label_type);
	const auto point_block = VtkBlock(point_bytes, header64), connectivity_block = VtkBlock(connectivity_bytes, header64), offsets_block = VtkBlock(offsets_bytes, header64), label_block = VtkBlock(label_bytes, header64);
	std::string stream = leading_appended_gap ? "AAAA" : ""; const std::array<std::vector<std::uint8_t>, 4> blocks{{point_block, connectivity_block, offsets_block, label_block}}; std::array<std::size_t, 4> offsets{};
	for (std::size_t block = 0; block < blocks.size(); ++block) { offsets[block] = stream.size(); stream += Base64(blocks[block]); }
	auto data = [&](std::size_t ordinal, const std::string& type, const std::string& name, const std::string& components = "") {
		std::ostringstream result; result << "<DataArray type=\"" << type << "\""; if (!name.empty()) result << " Name=\"" << name << "\""; if (!components.empty()) result << " NumberOfComponents=\"" << components << "\"";
		if (appended) result << " format=\"appended\" offset=\"" << offsets[ordinal] << "\"/>";
		else result << " format=\"binary\">" << Base64(blocks[ordinal]) << "</DataArray>";
		return result.str();
	};
	std::ostringstream result;
	result << "<VTKFile type=\"PolyData\" byte_order=\"LittleEndian\" header_type=\"" << (header64 ? "UInt64" : "UInt32") << "\"><PolyData><Piece NumberOfPoints=\"" << Points << "\" NumberOfPolys=\"" << Faces << "\">"
		<< "<Points>" << data(0, point_type, "", "3") << "</Points><Polys>" << data(1, index_type, "connectivity") << data(2, index_type, "offsets") << "</Polys>";
	if (!labels.empty()) result << "<CellData Scalars=\"boundary_id\">" << data(3, label_type, "boundary_id") << "</CellData>";
	result << "</Piece></PolyData>";
	if (appended) result << "<AppendedData encoding=\"base64\">_" << stream << "</AppendedData>";
	return result.str()+"</VTKFile>";
}

std::string WithOffset(std::string text, std::size_t ordinal, const std::string& value)
{
	std::size_t at = 0; for (std::size_t found = 0; found <= ordinal; ++found) { at = text.find("offset=\"", at); assert(at != std::string::npos); if (found == ordinal) break; at += 8; }
	const auto begin = at+8; const auto end = text.find('\"', begin); assert(end != std::string::npos); text.replace(begin, end-begin, value); return text;
}

std::string ReplaceFirstBinaryText(std::string text, const std::string& value)
{
	const auto begin = text.find("format=\"binary\">"); assert(begin != std::string::npos); const auto content = begin+16; const auto end = text.find("</DataArray>", content); assert(end != std::string::npos); text.replace(content, end-content, value); return text;
}

std::string ReplaceBinaryText(std::string text, std::size_t ordinal, const std::string& value)
{
	std::size_t at = 0; for (std::size_t found = 0; found <= ordinal; ++found) { at = text.find("format=\"binary\">", at); assert(at != std::string::npos); if (found == ordinal) break; at += 16; }
	const auto content = at+16; const auto end = text.find("</DataArray>", content); assert(end != std::string::npos); text.replace(content, end-content, value); return text;
}

std::string WithoutAppendedData(std::string text)
{
	const auto begin = text.find("<AppendedData"); assert(begin != std::string::npos); const auto end = text.find("</AppendedData>", begin); assert(end != std::string::npos); text.erase(begin, end+15-begin); return text;
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
	const std::vector<unsigned> tetra_labels{3, 3, 3, 3};
	const auto tetra_ascii_labelled = iga::SurfaceReaders::ReadVtpText(Vtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", tetra_labels));
	const auto tetra_inline_text = EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", tetra_labels, false, false);
	const auto tetra_inline = iga::SurfaceReaders::ReadVtpText(tetra_inline_text);
	const auto tetra_appended = iga::SurfaceReaders::ReadVtpText(EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int64", "UInt64", tetra_labels, true, true));
	assert(tetra_inline.CanonicalSha256() == tetra_ascii_labelled.CanonicalSha256());
	assert(tetra_appended.CanonicalSha256() == tetra_ascii_labelled.CanonicalSha256());
	assert(iga::SurfaceReaders::ReadVtpText(Replace(tetra_inline_text, " header_type=\"UInt32\"", "")).CanonicalSha256() == tetra_ascii_labelled.CanonicalSha256());
	const std::vector<unsigned> cube_labels(12, 7);
	const auto cube_ascii_binary_equivalent = iga::SurfaceReaders::ReadVtpText(Vtp(kCubePoints, kCubeFaces, "Float32", "Int64", "UInt32", cube_labels));
	const auto cube_inline = iga::SurfaceReaders::ReadVtpText(EncodedVtp(kCubePoints, kCubeFaces, "Float32", "Int64", "UInt32", cube_labels, false, true));
	const auto cube_appended = iga::SurfaceReaders::ReadVtpText(EncodedVtp(kCubePoints, kCubeFaces, "Float32", "Int32", "UInt64", cube_labels, true, false));
	assert(cube_inline.CanonicalSha256() == cube_ascii_binary_equivalent.CanonicalSha256());
	assert(cube_appended.CanonicalSha256() == cube_ascii_binary_equivalent.CanonicalSha256());

	const auto valid = Vtp(kTetraPoints, kTetraFaces);
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "type=\"PolyData\"", "type=\"UnstructuredGrid\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "byte_order=\"LittleEndian\"", "byte_order=\"BigEndian\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "header_type=\"UInt64\"", "compressor=\"vtkZLibDataCompressor\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"binary\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"appended\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "format=\"ascii\"", "format=\"base64\"")); });
	const auto binary_fixture = EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", tetra_labels, false, false);
	const auto appended_fixture = EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", tetra_labels, true, false);
	assert(appended_fixture.find('=') != appended_fixture.rfind('=')); // independently padded appended blocks
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, "!!!!")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, "A===")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, "AB==")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, Base64(std::vector<std::uint8_t>{0, 0, 0, 0}))); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, Base64(std::vector<std::uint8_t>{0, 0, 0}))); });
	std::vector<std::uint8_t> point_payload; for (const auto& point : kTetraPoints) for (const auto coordinate : point) AppendFloating(point_payload, coordinate, "Float64"); auto trailing_point_block = VtkBlock(point_payload, false); trailing_point_block.push_back(0);
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, Base64(trailing_point_block))); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(WithoutAppendedData(appended_fixture)); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(appended_fixture, "encoding=\"base64\"", "encoding=\"raw\"")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(WithOffset(appended_fixture, 0, "999999999")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(WithOffset(appended_fixture, 1, "1")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(WithOffset(appended_fixture, 1, "0")); });
	auto decoded_byte_offset_fixture = WithOffset(appended_fixture, 1, "100"); decoded_byte_offset_fixture = WithOffset(decoded_byte_offset_fixture, 2, "152"); decoded_byte_offset_fixture = WithOffset(decoded_byte_offset_fixture, 3, "172");
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(decoded_byte_offset_fixture); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(appended_fixture, ">_", ">_A ")); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "UInt32", tetra_labels, true, false, true)); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(valid, "</VTKFile>", "<AppendedData encoding=\"base64\">_AAAA</AppendedData></VTKFile>")); });
	std::vector<std::uint8_t> nan_payload = point_payload; nan_payload[0] = 0; nan_payload[1] = 0; nan_payload[2] = 0; nan_payload[3] = 0; nan_payload[4] = 0; nan_payload[5] = 0; nan_payload[6] = 0xf8; nan_payload[7] = 0x7f;
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceFirstBinaryText(binary_fixture, Base64(VtkBlock(nan_payload, false)))); });
	std::vector<std::uint8_t> bad_connectivity; for (const auto& face : kTetraFaces) for (std::size_t corner = 0; corner < 3; ++corner) AppendInteger(bad_connectivity, face == kTetraFaces[0] && corner == 0 ? 9 : face[corner], "Int32");
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceBinaryText(binary_fixture, 1, Base64(VtkBlock(bad_connectivity, false)))); });
	std::vector<std::uint8_t> bad_offsets; for (const auto offset : {3, 7, 9, 12}) AppendInteger(bad_offsets, offset, "Int32");
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceBinaryText(binary_fixture, 2, Base64(VtkBlock(bad_offsets, false)))); });
	const auto signed_label_fixture = EncodedVtp(kTetraPoints, kTetraFaces, "Float64", "Int32", "Int32", tetra_labels, false, false);
	std::vector<std::uint8_t> bad_labels; for (const auto label : {3, 3, 3, -1}) AppendInteger(bad_labels, label, "Int32");
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(ReplaceBinaryText(signed_label_fixture, 3, Base64(VtkBlock(bad_labels, false)))); });
	RequireRejected([&] { iga::SurfaceReaders::ReadVtpText(Replace(binary_fixture, "Float64", "Float16")); });
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

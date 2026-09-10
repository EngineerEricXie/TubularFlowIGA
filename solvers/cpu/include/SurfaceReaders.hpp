#ifndef IGA_SURFACE_READERS_HPP
#define IGA_SURFACE_READERS_HPP

#include "SurfaceGeometry.hpp"

#include <cerrno>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

class SurfaceReaders {
public:
	static ClosedTriangulatedSurface ReadStlBuffer(const std::vector<std::uint8_t>& bytes,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "buffer")
	{
		try {
			return ClosedTriangulatedSurface::Build(IsExactBinaryStl(bytes)
				? ParseBinary(bytes, options) : ParseAscii(bytes, options), options);
		} catch (const std::exception& error) {
			throw std::invalid_argument("STL " + context + ": " + error.what());
		}
	}

	static ClosedTriangulatedSurface ReadStlText(const std::string& text,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "text")
	{
		return ReadStlBuffer(std::vector<std::uint8_t>(text.begin(), text.end()), options, context);
	}

	static ClosedTriangulatedSurface ReadStlPath(const std::string& path,
		const SurfaceValidationOptions& options = SurfaceValidationOptions())
	{
		std::ifstream input(path, std::ios::binary);
		if (!input) throw std::invalid_argument("STL " + path + ": cannot open file");
		input.seekg(0, std::ios::end);
		const auto size = input.tellg();
		if (size < 0) throw std::invalid_argument("STL " + path + ": cannot determine file size");
		if (static_cast<std::uintmax_t>(size) > std::numeric_limits<std::size_t>::max()
			|| static_cast<std::uintmax_t>(size) > std::numeric_limits<std::streamsize>::max())
			throw std::invalid_argument("STL " + path + ": file is too large");
		input.seekg(0, std::ios::beg);
		std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
		if (!bytes.empty()) input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
		if (!input && !bytes.empty()) throw std::invalid_argument("STL " + path + ": cannot read file");
		return ReadStlBuffer(bytes, options, path);
	}

	// This supports the bounded, uncompressed VTK XML PolyData subset: ASCII,
	// inline base64 blocks, and base64-appended blocks.
	static ClosedTriangulatedSurface ReadVtpBuffer(const std::vector<std::uint8_t>& bytes,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "buffer",
		const std::string& boundary_array_name = "boundary_id")
	{
		try {
			if (boundary_array_name.empty()) throw std::invalid_argument("boundary array name is empty");
			if (options.max_triangles < 0) throw std::invalid_argument("configured maximum triangle count is negative");
			return ClosedTriangulatedSurface::Build(ParseVtp(bytes, options, boundary_array_name), options);
		} catch (const std::exception& error) {
			throw std::invalid_argument("VTP " + context + ": " + error.what());
		}
	}

	static ClosedTriangulatedSurface ReadVtpText(const std::string& text,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "text",
		const std::string& boundary_array_name = "boundary_id")
	{
		return ReadVtpBuffer(std::vector<std::uint8_t>(text.begin(), text.end()), options, context, boundary_array_name);
	}

	static ClosedTriangulatedSurface ReadVtpPath(const std::string& path,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& boundary_array_name = "boundary_id")
	{
		return ReadVtpBuffer(ReadVtpBytes(path), options, path, boundary_array_name);
	}

	// Material frames retain source vertex/facet order and directed connectivity.
	// Validate the closed surface without publishing its canonical reordering.
	// Scaling and welding would change the material coordinates or identities.
	static RawSurfaceSoup ReadMaterialVtpBuffer(const std::vector<std::uint8_t>& bytes,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "buffer",
		const std::string& boundary_array_name = "boundary_id")
	{
		try {
			if (boundary_array_name.empty()) throw std::invalid_argument("boundary array name is empty");
			if (options.max_triangles < 0) throw std::invalid_argument("configured maximum triangle count is negative");
			if (options.length_scale_to_m != 1.0 || options.weld_tolerance_m != 0.0)
				throw std::invalid_argument("material VTP requires metre coordinates without welding");
			auto source = ParseVtp(bytes, options, boundary_array_name);
			// Canonical construction normally resolves absent labels. Material
			// consumers need those effective labels on the retained source too.
			for (auto& triangle : source.triangles)
				if (triangle.boundary_id < 0) triangle.boundary_id = options.default_boundary_id;
			(void)ClosedTriangulatedSurface::Build(source, options);
			return source;
		} catch (const std::exception& error) {
			throw std::invalid_argument("VTP " + context + ": " + error.what());
		}
	}

	static RawSurfaceSoup ReadMaterialVtpText(const std::string& text,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& context = "text",
		const std::string& boundary_array_name = "boundary_id")
	{
		return ReadMaterialVtpBuffer(std::vector<std::uint8_t>(text.begin(), text.end()), options, context, boundary_array_name);
	}

	static RawSurfaceSoup ReadMaterialVtpPath(const std::string& path,
		const SurfaceValidationOptions& options = SurfaceValidationOptions(), const std::string& boundary_array_name = "boundary_id")
	{
		return ReadMaterialVtpBuffer(ReadVtpBytes(path), options, path, boundary_array_name);
	}

private:
	static std::vector<std::uint8_t> ReadVtpBytes(const std::string& path)
	{
		std::ifstream input(path, std::ios::binary);
		if (!input) throw std::invalid_argument("VTP " + path + ": cannot open file");
		input.seekg(0, std::ios::end);
		const auto size = input.tellg();
		if (size < 0) throw std::invalid_argument("VTP " + path + ": cannot determine file size");
		if (static_cast<std::uintmax_t>(size) > kVtpMaximumBytes
			|| static_cast<std::uintmax_t>(size) > std::numeric_limits<std::size_t>::max()
			|| static_cast<std::uintmax_t>(size) > std::numeric_limits<std::streamsize>::max())
			throw std::invalid_argument("VTP " + path + ": file is too large");
		input.seekg(0, std::ios::beg);
		std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
		if (!bytes.empty()) input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
		if (!input && !bytes.empty()) throw std::invalid_argument("VTP " + path + ": cannot read file");
		return bytes;
	}

	static constexpr std::size_t kVtpMaximumBytes = 64U*1024U*1024U;
	static constexpr std::size_t kVtpMaximumTokenBytes = 1024U*1024U;
	static constexpr std::size_t kVtpMaximumDepth = 32U;
	static constexpr std::size_t kVtpMaximumNodes = 128U;

	struct XmlElement {
		std::string name;
		std::map<std::string, std::string> attributes;
		std::string text;
		std::vector<std::unique_ptr<XmlElement>> children;
	};

	static bool IsXmlWhitespace(char value)
	{
		return value == ' ' || value == '\t' || value == '\r' || value == '\n';
	}
	static bool IsXmlNameStart(char value)
	{
		return std::isalpha(static_cast<unsigned char>(value)) || value == '_';
	}
	static bool IsXmlNameChar(char value)
	{
		return IsXmlNameStart(value) || std::isdigit(static_cast<unsigned char>(value)) || value == '-' || value == '.' || value == ':';
	}
	static void SkipXmlWhitespace(const std::string& text, std::size_t& at)
	{
		while (at < text.size() && IsXmlWhitespace(text[at])) ++at;
	}
	static std::string ReadXmlName(const std::string& text, std::size_t& at)
	{
		if (at == text.size() || !IsXmlNameStart(text[at])) throw std::invalid_argument("XML name is malformed");
		const auto begin = at++;
		while (at < text.size() && IsXmlNameChar(text[at])) ++at;
		if (at-begin > 128U) throw std::invalid_argument("XML name exceeds parser limit");
		return text.substr(begin, at-begin);
	}
	static void AppendXmlText(XmlElement& element, const std::string& source, std::size_t begin, std::size_t end)
	{
		const auto limit = element.name == "DataArray" || element.name == "AppendedData" ? kVtpMaximumBytes : kVtpMaximumTokenBytes;
		if (end < begin || end-begin > limit || element.text.size() > limit-(end-begin))
			throw std::invalid_argument("XML text token exceeds parser limit");
		for (std::size_t i = begin; i < end; ++i)
			if (source[i] == '&' || (static_cast<unsigned char>(source[i]) < 0x20U && !IsXmlWhitespace(source[i])))
				throw std::invalid_argument("XML entities or control characters are unsupported");
		element.text.append(source, begin, end-begin);
	}
	static std::unique_ptr<XmlElement> ParseXml(const std::vector<std::uint8_t>& bytes)
	{
		if (bytes.empty()) throw std::invalid_argument("XML document is empty");
		if (bytes.size() > kVtpMaximumBytes) throw std::invalid_argument("VTP exceeds parser limit");
		std::string text; text.reserve(bytes.size());
		for (const auto byte : bytes) {
			if (byte > 0x7fU) throw std::invalid_argument("XML is not ASCII");
			text.push_back(static_cast<char>(byte));
		}
		std::unique_ptr<XmlElement> root;
		std::vector<XmlElement*> stack;
		std::size_t nodes = 0, at = 0; bool declaration_seen = false, root_seen = false;
		while (at < text.size()) {
			if (text[at] != '<') {
				const auto begin = at; while (at < text.size() && text[at] != '<') ++at;
				if (stack.empty()) { for (std::size_t i = begin; i < at; ++i) if (!IsXmlWhitespace(text[i])) throw std::invalid_argument("XML has text outside its root"); }
				else AppendXmlText(*stack.back(), text, begin, at);
				continue;
			}
			if (text.compare(at, 4, "<!--") == 0) {
				const auto end = text.find("-->", at+4); if (end == std::string::npos || end-(at+4) > kVtpMaximumTokenBytes) throw std::invalid_argument("XML comment is malformed or too large");
				if (text.find("--", at+4) < end) throw std::invalid_argument("XML comment is malformed");
				at = end+3;
				continue;
			}
			if (text.compare(at, 5, "<?xml") == 0) {
				if (declaration_seen || root_seen || at != 0 || (at+5 < text.size() && !IsXmlWhitespace(text[at+5]))) throw std::invalid_argument("XML declaration is not leading");
				const auto end = text.find("?>", at+5); if (end == std::string::npos || end-(at+5) > kVtpMaximumTokenBytes || text.find('&', at+5) < end) throw std::invalid_argument("XML declaration is malformed");
				declaration_seen = true; at = end+2; continue;
			}
			if (text.compare(at, 2, "<?") == 0 || text.compare(at, 2, "<!") == 0)
				throw std::invalid_argument("XML processing instructions, DTDs, and entities are unsupported");
			if (text.compare(at, 2, "</") == 0) {
				at += 2; const auto name = ReadXmlName(text, at); SkipXmlWhitespace(text, at);
				if (at == text.size() || text[at++] != '>' || stack.empty() || stack.back()->name != name) throw std::invalid_argument("XML nesting is malformed");
				stack.pop_back(); continue;
			}
			++at; auto element = std::make_unique<XmlElement>(); element->name = ReadXmlName(text, at);
			for (;;) {
				SkipXmlWhitespace(text, at); if (at == text.size()) throw std::invalid_argument("XML start tag is truncated");
				if (text[at] == '>') { ++at; break; }
				if (text[at] == '/' && at+1 < text.size() && text[at+1] == '>') { at += 2; break; }
				const auto key = ReadXmlName(text, at); SkipXmlWhitespace(text, at);
				if (at == text.size() || text[at++] != '=') throw std::invalid_argument("XML attribute is malformed");
				SkipXmlWhitespace(text, at); if (at == text.size() || (text[at] != '\'' && text[at] != '\"')) throw std::invalid_argument("XML attribute quote is malformed");
				const char quote = text[at++]; const auto begin = at; while (at < text.size() && text[at] != quote) { if (text[at] == '<' || text[at] == '&') throw std::invalid_argument("XML attribute escape is unsupported"); ++at; }
				if (at == text.size() || at-begin > 4096U) throw std::invalid_argument("XML attribute is malformed or too large");
				if (!element->attributes.emplace(key, text.substr(begin, at-begin)).second) throw std::invalid_argument("XML has duplicate attribute");
				++at;
			}
			if (++nodes > kVtpMaximumNodes) throw std::invalid_argument("XML has too many elements");
			XmlElement* current = element.get();
			if (stack.empty()) { if (root) throw std::invalid_argument("XML has multiple roots"); root = std::move(element); root_seen = true; }
			else { stack.back()->children.push_back(std::move(element)); }
			// A self-closing tag is complete only when it was not placed on the stack.
			if (at >= 2 && text[at-2] == '/' && text[at-1] == '>') continue;
			if (stack.size() >= kVtpMaximumDepth) throw std::invalid_argument("XML nesting exceeds parser limit");
			stack.push_back(current);
		}
		if (!root || !stack.empty()) throw std::invalid_argument("XML document is truncated");
		return root;
	}
	static const std::string* Attribute(const XmlElement& element, const char* key)
	{
		const auto found = element.attributes.find(key); return found == element.attributes.end() ? nullptr : &found->second;
	}
	static const std::string& RequireAttribute(const XmlElement& element, const char* key)
	{
		const auto value = Attribute(element, key); if (!value) throw std::invalid_argument("XML element '"+element.name+"' lacks attribute '"+key+"'"); return *value;
	}
	static void CheckAttributes(const XmlElement& element, std::initializer_list<const char*> allowed)
	{
		for (const auto& attribute : element.attributes) {
			bool known = false; for (const auto* name : allowed) if (attribute.first == name) { known = true; break; }
			if (!known) throw std::invalid_argument("XML element '"+element.name+"' has unsupported attribute '"+attribute.first+"'");
		}
	}
	static void RequireWhitespace(const XmlElement& element)
	{
		for (const auto value : element.text) if (!IsXmlWhitespace(value)) throw std::invalid_argument("XML element '"+element.name+"' has unexpected text");
	}
	static std::uint64_t ParseUnsigned(const std::string& token, const char* description)
	{
		if (token.empty()) throw std::invalid_argument(std::string(description)+" is empty");
		std::uint64_t value = 0;
		for (const auto character : token) {
			if (character < '0' || character > '9') throw std::invalid_argument(std::string(description)+" is invalid");
			const auto digit = static_cast<unsigned>(character-'0'); if (value > (std::numeric_limits<std::uint64_t>::max()-digit)/10U) throw std::invalid_argument(std::string(description)+" overflows"); value = 10U*value+digit;
		}
		return value;
	}
	static std::size_t ParseCount(const std::string& token, const char* description)
	{
		const auto value = ParseUnsigned(token, description); if (value > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument(std::string(description)+" exceeds platform range"); return static_cast<std::size_t>(value);
	}
	static std::vector<std::string> SplitAsciiValues(const std::string& text, std::size_t expected, const char* description)
	{
		if (expected > kVtpMaximumTokenBytes || expected > text.size()) throw std::invalid_argument(std::string(description)+" has an invalid number of values");
		std::vector<std::string> values; values.reserve(expected); std::size_t at = 0;
		while (at < text.size()) {
			while (at < text.size() && IsXmlWhitespace(text[at])) ++at;
			if (at == text.size()) break;
			const auto begin = at; while (at < text.size() && !IsXmlWhitespace(text[at])) ++at;
			if (at-begin > 256U || values.size() == expected) throw std::invalid_argument(std::string(description)+" has an invalid number of values");
			values.emplace_back(text, begin, at-begin);
		}
		if (values.size() != expected) throw std::invalid_argument(std::string(description)+" has an invalid number of values");
		return values;
	}
	static double ParseVtpFloat(const std::string& token, bool float32)
	{
		errno = 0; char* end = nullptr;
		if (float32) {
			const float value = std::strtof(token.c_str(), &end);
			if (end != token.c_str()+token.size() || errno == ERANGE || !std::isfinite(value)) throw std::invalid_argument("VTP floating-point value is invalid");
			return static_cast<double>(value);
		}
		const double value = std::strtod(token.c_str(), &end);
		if (end != token.c_str()+token.size() || errno == ERANGE || !std::isfinite(value)) throw std::invalid_argument("VTP floating-point value is invalid");
		return value;
	}
	static std::int64_t ParseVtpInteger(const std::string& token, const std::string& type)
	{
		const bool unsigned_type = type == "UInt32" || type == "UInt64";
		if (type != "Int32" && type != "Int64" && !unsigned_type) throw std::invalid_argument("VTP integral type is unsupported");
		if (token.empty()) throw std::invalid_argument("VTP integer value is invalid");
		if (unsigned_type) {
			const auto value = ParseUnsigned(token, "VTP integer value");
			const auto maximum = type == "UInt32" ? static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) : std::numeric_limits<std::uint64_t>::max();
			if (value > maximum || value > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) throw std::invalid_argument("VTP integer value is out of range");
			return static_cast<std::int64_t>(value);
		}
		std::size_t at = token[0] == '-' ? 1U : 0U; if (at == token.size()) throw std::invalid_argument("VTP integer value is invalid");
		for (std::size_t i = at; i < token.size(); ++i) if (token[i] < '0' || token[i] > '9') throw std::invalid_argument("VTP integer value is invalid");
		errno = 0; char* end = nullptr; const auto value = std::strtoll(token.c_str(), &end, 10);
		if (end != token.c_str()+token.size() || errno == ERANGE || (type == "Int32" && (value < std::numeric_limits<std::int32_t>::min() || value > std::numeric_limits<std::int32_t>::max()))) throw std::invalid_argument("VTP integer value is out of range");
		return static_cast<std::int64_t>(value);
	}
	struct VtpBlockRange { std::size_t begin = 0; std::size_t end = 0; };
	static std::size_t VtpTypeSize(const std::string& type)
	{
		if (type == "Float32" || type == "Int32" || type == "UInt32") return 4U;
		if (type == "Float64" || type == "Int64" || type == "UInt64") return 8U;
		throw std::invalid_argument("VTP DataArray type is unsupported");
	}
	static std::size_t CheckedVtpBytes(std::size_t values, std::size_t value_size)
	{
		if (values > kVtpMaximumBytes/value_size) throw std::invalid_argument("VTP binary payload exceeds parser limit");
		return values*value_size;
	}
	static int Base64Value(char value)
	{
		if (value >= 'A' && value <= 'Z') return value-'A';
		if (value >= 'a' && value <= 'z') return value-'a'+26;
		if (value >= '0' && value <= '9') return value-'0'+52;
		if (value == '+') return 62;
		if (value == '/') return 63;
		return -1;
	}
	static std::vector<std::uint8_t> DecodeBase64(const std::string& text, const char* description)
	{
		std::string compact; compact.reserve(text.size());
		for (const auto value : text) {
			if (IsXmlWhitespace(value)) continue;
			if (value == '=' || Base64Value(value) >= 0) compact.push_back(value);
			else throw std::invalid_argument(std::string(description)+" has invalid base64 data");
		}
		if (compact.empty() || compact.size()%4U != 0U) throw std::invalid_argument(std::string(description)+" has malformed base64 padding");
		const std::size_t padding = compact[compact.size()-1] == '=' ? (compact[compact.size()-2] == '=' ? 2U : 1U) : 0U;
		if (padding > 2U || compact.size()/4U > (kVtpMaximumBytes+2U)/3U) throw std::invalid_argument(std::string(description)+" exceeds decoded-size limit");
		for (std::size_t at = 0; at < compact.size(); at += 4U) {
			const bool last = at+4U == compact.size();
			if (compact[at] == '=' || compact[at+1] == '=' || (!last && (compact[at+2] == '=' || compact[at+3] == '=')))
				throw std::invalid_argument(std::string(description)+" has malformed base64 padding");
			if ((!last && (compact[at+2] == '=' || compact[at+3] == '=')) || (last && (compact[at+2] == '=' && compact[at+3] != '=')))
				throw std::invalid_argument(std::string(description)+" has malformed base64 padding");
		}
		std::vector<std::uint8_t> result; result.reserve(3U*(compact.size()/4U)-padding);
		for (std::size_t at = 0; at < compact.size(); at += 4U) {
			const int first = Base64Value(compact[at]), second = Base64Value(compact[at+1]);
			const int third = compact[at+2] == '=' ? 0 : Base64Value(compact[at+2]);
			const int fourth = compact[at+3] == '=' ? 0 : Base64Value(compact[at+3]);
			if (first < 0 || second < 0 || third < 0 || fourth < 0) throw std::invalid_argument(std::string(description)+" has invalid base64 data");
			if (at+4U == compact.size() && ((compact[at+2] == '=' && (second&15) != 0) || (compact[at+3] == '=' && compact[at+2] != '=' && (third&3) != 0)))
				throw std::invalid_argument(std::string(description)+" has noncanonical base64 padding");
			result.push_back(static_cast<std::uint8_t>((first<<2)|(second>>4)));
			if (compact[at+2] != '=') result.push_back(static_cast<std::uint8_t>((second<<4)|(third>>2)));
			if (compact[at+3] != '=') result.push_back(static_cast<std::uint8_t>((third<<6)|fourth));
		}
		return result;
	}
	static std::uint64_t ReadLe64(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		std::uint64_t result = 0;
		for (std::size_t shift = 0; shift < 8U; ++shift) result |= static_cast<std::uint64_t>(bytes[offset+shift]) << (8U*shift);
		return result;
	}
	static std::size_t ReadVtpBlock(const std::vector<std::uint8_t>& stream, std::size_t offset, std::size_t header_size,
		std::size_t expected_payload, std::vector<std::uint8_t>& payload)
	{
		if (offset > stream.size() || header_size > stream.size()-offset) throw std::invalid_argument("VTP binary block header is truncated");
		const auto length = header_size == 4U ? static_cast<std::uint64_t>(ReadLe32(stream, offset)) : ReadLe64(stream, offset);
		if (length != expected_payload || length > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("VTP binary block length is invalid");
		const auto payload_size = static_cast<std::size_t>(length);
		if (payload_size > stream.size()-offset-header_size) throw std::invalid_argument("VTP binary block payload is truncated");
		const auto end = offset+header_size+payload_size;
		payload.assign(stream.begin()+static_cast<std::ptrdiff_t>(offset+header_size), stream.begin()+static_cast<std::ptrdiff_t>(end));
		return end;
	}
	static const std::string& DataArrayFormat(const XmlElement& array)
	{
		const auto& format = RequireAttribute(array, "format");
		if (format != "ascii" && format != "binary" && format != "appended") throw std::invalid_argument("VTP DataArray format is unsupported");
		return format;
	}
	static void ValidateDataArray(const XmlElement& array, const std::string& expected_type, const char* expected_name, bool components_three)
	{
		CheckAttributes(array, {"type", "Name", "NumberOfComponents", "format", "offset"});
		if (array.name != "DataArray" || RequireAttribute(array, "type") != expected_type) throw std::invalid_argument("VTP DataArray type is unsupported");
		if (expected_name && RequireAttribute(array, "Name") != expected_name) throw std::invalid_argument("VTP DataArray name is unexpected");
		const auto components = Attribute(array, "NumberOfComponents");
		if (components_three ? (!components || *components != "3") : (components && *components != "1")) throw std::invalid_argument("VTP DataArray component count is unsupported");
		if (!array.children.empty()) throw std::invalid_argument("VTP DataArray cannot contain XML elements");
		const auto& format = DataArrayFormat(array);
		if (format == "appended") {
			(void)ParseUnsigned(RequireAttribute(array, "offset"), "VTP appended offset");
			RequireWhitespace(array);
		} else if (Attribute(array, "offset")) throw std::invalid_argument("VTP inline DataArray has an offset");
	}
	static std::vector<std::uint8_t> BinaryDataArrayPayload(const XmlElement& array, std::size_t values,
		std::size_t header_size, const std::string* appended, std::vector<VtpBlockRange>& appended_ranges)
	{
		const auto expected_payload = CheckedVtpBytes(values, VtpTypeSize(RequireAttribute(array, "type")));
		const auto& format = DataArrayFormat(array); std::vector<std::uint8_t> payload;
		if (format == "binary") {
			const auto encoded = DecodeBase64(array.text, "VTP inline DataArray");
			const auto end = ReadVtpBlock(encoded, 0, header_size, expected_payload, payload);
			if (end != encoded.size()) throw std::invalid_argument("VTP inline DataArray has trailing decoded bytes");
			return payload;
		}
		if (format != "appended" || !appended) throw std::invalid_argument("VTP DataArray has no binary payload source");
		const auto offset_value = ParseUnsigned(RequireAttribute(array, "offset"), "VTP appended offset");
		if (offset_value > std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("VTP appended offset exceeds platform range");
		const auto offset = static_cast<std::size_t>(offset_value);
		if (expected_payload > std::numeric_limits<std::size_t>::max()-header_size) throw std::invalid_argument("VTP binary block size overflows");
		const auto raw_size = header_size+expected_payload;
		if (raw_size > (std::numeric_limits<std::size_t>::max()-2U)/4U) throw std::invalid_argument("VTP appended encoded span overflows");
		const auto encoded_size = 4U*((raw_size+2U)/3U);
		if (offset > appended->size() || encoded_size > appended->size()-offset) throw std::invalid_argument("VTP appended offset is out of range");
		const auto decoded = DecodeBase64(appended->substr(offset, encoded_size), "VTP appended DataArray");
		if (decoded.size() != raw_size) throw std::invalid_argument("VTP appended DataArray has an invalid encoded span");
		const auto end = ReadVtpBlock(decoded, 0, header_size, expected_payload, payload);
		if (end != decoded.size()) throw std::invalid_argument("VTP appended DataArray has trailing decoded bytes");
		appended_ranges.push_back({offset, offset+encoded_size});
		return payload;
	}
	static double ReadBinaryVtpFloat(const std::vector<std::uint8_t>& bytes, std::size_t offset, const std::string& type)
	{
		if (type == "Float32") {
			static_assert(sizeof(float) == 4U, "VTP Float32 requires 32-bit float");
			const auto bits = ReadLe32(bytes, offset); float value = 0.0F; static_assert(std::numeric_limits<float>::is_iec559, "VTP Float32 requires IEEE-754"); std::memcpy(&value, &bits, sizeof(value));
			if (!std::isfinite(value)) throw std::invalid_argument("VTP floating-point value is invalid");
			return static_cast<double>(value);
		}
		if (type == "Float64") {
			static_assert(sizeof(double) == 8U, "VTP Float64 requires 64-bit double");
			const auto bits = ReadLe64(bytes, offset); double value = 0.0; static_assert(std::numeric_limits<double>::is_iec559, "VTP Float64 requires IEEE-754"); std::memcpy(&value, &bits, sizeof(value));
			if (!std::isfinite(value)) throw std::invalid_argument("VTP floating-point value is invalid");
			return value;
		}
		throw std::invalid_argument("VTP point type is unsupported");
	}
	static std::int64_t ReadBinaryVtpInteger(const std::vector<std::uint8_t>& bytes, std::size_t offset, const std::string& type)
	{
		if (type == "Int32") { const auto value = ReadLe32(bytes, offset); return value <= static_cast<std::uint32_t>(std::numeric_limits<std::int32_t>::max()) ? static_cast<std::int64_t>(value) : static_cast<std::int64_t>(value)-4294967296LL; }
		if (type == "UInt32") return static_cast<std::int64_t>(ReadLe32(bytes, offset));
		const auto value = ReadLe64(bytes, offset);
		if (type == "UInt64") { if (value > static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) throw std::invalid_argument("VTP integer value is out of range"); return static_cast<std::int64_t>(value); }
		if (type == "Int64") return value <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) ? static_cast<std::int64_t>(value) : std::numeric_limits<std::int64_t>::min()+static_cast<std::int64_t>(value & static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()));
		throw std::invalid_argument("VTP integral type is unsupported");
	}
	static void ValidateAppendedRanges(std::vector<VtpBlockRange>& ranges, std::size_t stream_size)
	{
		if (ranges.empty()) throw std::invalid_argument("VTP AppendedData is unused");
		std::sort(ranges.begin(), ranges.end(), [](const VtpBlockRange& left, const VtpBlockRange& right) { return left.begin < right.begin; });
		std::size_t cursor = 0;
		for (const auto& range : ranges) { if (range.begin != cursor || range.end < range.begin) throw std::invalid_argument("VTP appended blocks overlap or leave a gap"); cursor = range.end; }
		if (cursor != stream_size) throw std::invalid_argument("VTP appended data has trailing bytes");
	}
	static std::string DecodeAppendedData(const XmlElement& appended)
	{
		CheckAttributes(appended, {"encoding"});
		if (appended.name != "AppendedData" || RequireAttribute(appended, "encoding") != "base64" || !appended.children.empty()) throw std::invalid_argument("VTP AppendedData encoding is unsupported");
		std::size_t at = 0; while (at < appended.text.size() && IsXmlWhitespace(appended.text[at])) ++at;
		if (at == appended.text.size() || appended.text[at++] != '_') throw std::invalid_argument("VTP AppendedData lacks its base64 sentinel");
		std::string encoded; encoded.reserve(appended.text.size()-at); bool trailing_whitespace = false;
		for (; at < appended.text.size(); ++at) {
			const auto value = appended.text[at];
			if (IsXmlWhitespace(value)) { trailing_whitespace = true; continue; }
			if (trailing_whitespace || (value != '=' && Base64Value(value) < 0)) throw std::invalid_argument("VTP AppendedData has invalid encoded text");
			encoded.push_back(value);
		}
		if (encoded.empty()) throw std::invalid_argument("VTP AppendedData is empty");
		return encoded;
	}
	static RawSurfaceSoup ParseVtp(const std::vector<std::uint8_t>& bytes, const SurfaceValidationOptions& options, const std::string& boundary_array_name)
	{
		auto root = ParseXml(bytes); CheckAttributes(*root, {"type", "version", "byte_order", "header_type", "compressor"});
		if (root->name != "VTKFile" || RequireAttribute(*root, "type") != "PolyData") throw std::invalid_argument("VTP root is not VTK PolyData");
		if (const auto value = Attribute(*root, "byte_order"); value && *value != "LittleEndian") throw std::invalid_argument("VTP byte order is unsupported");
		if (const auto value = Attribute(*root, "header_type"); value && *value != "UInt32" && *value != "UInt64") throw std::invalid_argument("VTP header type is unsupported");
		if (Attribute(*root, "compressor")) throw std::invalid_argument("VTP compression is unsupported");
		RequireWhitespace(*root);
		const XmlElement *polydata_element = nullptr, *appended_element = nullptr;
		for (const auto& child : root->children) {
			if (child->name == "PolyData" && !polydata_element) polydata_element = child.get();
			else if (child->name == "AppendedData" && !appended_element) appended_element = child.get();
			else throw std::invalid_argument("VTP root has duplicate or unsupported content");
		}
		if (!polydata_element) throw std::invalid_argument("VTP requires exactly one PolyData element");
		const auto& polydata = *polydata_element; CheckAttributes(polydata, {}); RequireWhitespace(polydata);
		if (polydata.children.size() != 1U || polydata.children[0]->name != "Piece") throw std::invalid_argument("VTP requires exactly one Piece element");
		const auto& piece = *polydata.children[0]; CheckAttributes(piece, {"NumberOfPoints", "NumberOfPolys", "NumberOfVerts", "NumberOfLines", "NumberOfStrips"}); RequireWhitespace(piece);
		const auto points_count = ParseCount(RequireAttribute(piece, "NumberOfPoints"), "NumberOfPoints");
		const auto polygons_count = ParseCount(RequireAttribute(piece, "NumberOfPolys"), "NumberOfPolys");
		if (points_count > std::numeric_limits<std::uint32_t>::max() || polygons_count > static_cast<std::uintmax_t>(options.max_triangles)
			|| polygons_count > static_cast<std::uintmax_t>(std::numeric_limits<std::int64_t>::max())/3U) throw std::invalid_argument("VTP count exceeds configured maximum");
		for (const auto* name : {"NumberOfVerts", "NumberOfLines", "NumberOfStrips"}) if (const auto value = Attribute(piece, name); value && ParseUnsigned(*value, name) != 0U) throw std::invalid_argument("VTP supports no verts, lines, or strips");
		const XmlElement *points = nullptr, *polys = nullptr, *cell_data = nullptr;
		for (const auto& child : piece.children) {
			if (child->name == "Points" && !points) points = child.get(); else if (child->name == "Polys" && !polys) polys = child.get(); else if (child->name == "CellData" && !cell_data) cell_data = child.get();
			else if (child->name == "Verts" || child->name == "Lines" || child->name == "Strips") { CheckAttributes(*child, {}); RequireWhitespace(*child); if (!child->children.empty()) throw std::invalid_argument("VTP supports no verts, lines, or strips"); }
			else throw std::invalid_argument("VTP Piece has duplicate or unsupported container");
		}
		if (!points || !polys) throw std::invalid_argument("VTP Piece lacks Points or Polys");
		std::string appended_stream;
		if (appended_element) appended_stream = DecodeAppendedData(*appended_element);
		const std::string* appended = appended_element ? &appended_stream : nullptr;
		std::vector<VtpBlockRange> appended_ranges; bool uses_appended = false;
		CheckAttributes(*points, {}); RequireWhitespace(*points); if (points->children.size() != 1U) throw std::invalid_argument("VTP Points requires exactly one DataArray");
		const auto& point_array = *points->children[0]; const auto& point_type = RequireAttribute(point_array, "type");
		if (point_type != "Float32" && point_type != "Float64") throw std::invalid_argument("VTP point type is unsupported");
		ValidateDataArray(point_array, point_type, nullptr, true);
		if (points_count > std::numeric_limits<std::size_t>::max()/3U) throw std::invalid_argument("VTP point count overflows");
		const auto point_value_count = 3U*points_count; RawSurfaceSoup soup; soup.vertices.reserve(points_count);
		if (DataArrayFormat(point_array) == "ascii") {
			const auto point_values = SplitAsciiValues(point_array.text, point_value_count, "VTP point DataArray");
			for (std::size_t point = 0; point < points_count; ++point) soup.vertices.push_back({{ParseVtpFloat(point_values[3*point], point_type == "Float32"), ParseVtpFloat(point_values[3*point+1], point_type == "Float32"), ParseVtpFloat(point_values[3*point+2], point_type == "Float32")}});
		} else {
			uses_appended = uses_appended || DataArrayFormat(point_array) == "appended";
			const auto point_values = BinaryDataArrayPayload(point_array, point_value_count, Attribute(*root, "header_type") && *Attribute(*root, "header_type") == "UInt64" ? 8U : 4U, appended, appended_ranges);
			const auto value_size = VtpTypeSize(point_type);
			for (std::size_t point = 0; point < points_count; ++point) soup.vertices.push_back({{ReadBinaryVtpFloat(point_values, value_size*(3*point), point_type), ReadBinaryVtpFloat(point_values, value_size*(3*point+1), point_type), ReadBinaryVtpFloat(point_values, value_size*(3*point+2), point_type)}});
		}
		CheckAttributes(*polys, {}); RequireWhitespace(*polys); if (polys->children.size() != 2U) throw std::invalid_argument("VTP Polys requires connectivity and offsets");
		const XmlElement *connectivity = nullptr, *offsets = nullptr;
		for (const auto& child : polys->children) {
			if (child->name != "DataArray") throw std::invalid_argument("VTP Polys has unsupported content");
			const auto& name = RequireAttribute(*child, "Name");
			if (name == "connectivity" && !connectivity) connectivity = child.get(); else if (name == "offsets" && !offsets) offsets = child.get(); else throw std::invalid_argument("VTP Polys has duplicate or unknown DataArray");
		}
		if (!connectivity || !offsets) throw std::invalid_argument("VTP Polys lacks connectivity or offsets");
		const auto& connectivity_type = RequireAttribute(*connectivity, "type"); const auto& offsets_type = RequireAttribute(*offsets, "type");
		if ((connectivity_type != "Int32" && connectivity_type != "Int64") || (offsets_type != "Int32" && offsets_type != "Int64")) throw std::invalid_argument("VTP polygon type is unsupported");
		ValidateDataArray(*connectivity, connectivity_type, "connectivity", false); ValidateDataArray(*offsets, offsets_type, "offsets", false);
		if (polygons_count > std::numeric_limits<std::size_t>::max()/3U) throw std::invalid_argument("VTP polygon count overflows");
		const auto connectivity_count = 3U*polygons_count;
		const auto header_size = Attribute(*root, "header_type") && *Attribute(*root, "header_type") == "UInt64" ? 8U : 4U;
		const auto connectivity_format = DataArrayFormat(*connectivity), offsets_format = DataArrayFormat(*offsets);
		if (connectivity_format == "appended" || offsets_format == "appended") uses_appended = true;
		const auto connectivity_values = connectivity_format == "ascii" ? std::vector<std::string>(SplitAsciiValues(connectivity->text, connectivity_count, "VTP connectivity DataArray")) : std::vector<std::string>();
		const auto offset_values = offsets_format == "ascii" ? std::vector<std::string>(SplitAsciiValues(offsets->text, polygons_count, "VTP offsets DataArray")) : std::vector<std::string>();
		const auto connectivity_binary = connectivity_format == "ascii" ? std::vector<std::uint8_t>() : BinaryDataArrayPayload(*connectivity, connectivity_count, header_size, appended, appended_ranges);
		const auto offsets_binary = offsets_format == "ascii" ? std::vector<std::uint8_t>() : BinaryDataArrayPayload(*offsets, polygons_count, header_size, appended, appended_ranges);
		soup.triangles.reserve(polygons_count); std::int64_t previous = 0;
		for (std::size_t polygon = 0; polygon < polygons_count; ++polygon) {
			const auto offset = offsets_format == "ascii" ? ParseVtpInteger(offset_values[polygon], offsets_type) : ReadBinaryVtpInteger(offsets_binary, VtpTypeSize(offsets_type)*polygon, offsets_type); if (offset != previous+3 || offset < 0) throw std::invalid_argument("VTP offsets must describe triangles"); previous = offset;
			RawSurfaceTriangle triangle; for (std::size_t corner = 0; corner < 3; ++corner) { const auto index = connectivity_format == "ascii" ? ParseVtpInteger(connectivity_values[3*polygon+corner], connectivity_type) : ReadBinaryVtpInteger(connectivity_binary, VtpTypeSize(connectivity_type)*(3*polygon+corner), connectivity_type); if (index < 0 || static_cast<std::uintmax_t>(index) >= points_count) throw std::invalid_argument("VTP connectivity index is out of range"); triangle.indices[corner] = index; } soup.triangles.push_back(triangle);
		}
		if (previous != static_cast<std::int64_t>(3U*polygons_count)) throw std::invalid_argument("VTP final polygon offset is invalid");
		if (cell_data) {
			CheckAttributes(*cell_data, {"Scalars"}); RequireWhitespace(*cell_data);
			if (!cell_data->children.empty()) {
				if (const auto scalar_name = Attribute(*cell_data, "Scalars"); scalar_name && *scalar_name != boundary_array_name) throw std::invalid_argument("VTP CellData scalar is unexpected");
				if (cell_data->children.size() != 1U || cell_data->children[0]->name != "DataArray") throw std::invalid_argument("VTP CellData requires exactly one configured DataArray");
				const auto& labels = *cell_data->children[0]; const auto& type = RequireAttribute(labels, "type"); if (type != "Int32" && type != "UInt32" && type != "Int64" && type != "UInt64") throw std::invalid_argument("VTP boundary id type is unsupported"); ValidateDataArray(labels, type, boundary_array_name.c_str(), false);
				const auto& format = DataArrayFormat(labels); if (format == "appended") uses_appended = true;
				const auto values = format == "ascii" ? std::vector<std::string>(SplitAsciiValues(labels.text, polygons_count, "VTP boundary DataArray")) : std::vector<std::string>();
				const auto binary = format == "ascii" ? std::vector<std::uint8_t>() : BinaryDataArrayPayload(labels, polygons_count, header_size, appended, appended_ranges);
				for (std::size_t polygon = 0; polygon < polygons_count; ++polygon) { const auto label = format == "ascii" ? ParseVtpInteger(values[polygon], type) : ReadBinaryVtpInteger(binary, VtpTypeSize(type)*polygon, type); if (label < 0 || static_cast<std::uintmax_t>(label) > std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("VTP boundary id is out of range"); soup.triangles[polygon].boundary_id = label; }
			}
		}
		if (uses_appended) {
			if (!appended_element) throw std::invalid_argument("VTP appended DataArray lacks AppendedData");
			ValidateAppendedRanges(appended_ranges, appended_stream.size());
		} else if (appended_element) throw std::invalid_argument("VTP AppendedData is unused");
		return soup;
	}
	static std::uint32_t ReadLe32(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		return static_cast<std::uint32_t>(bytes[offset]) | (static_cast<std::uint32_t>(bytes[offset+1])<<8)
			| (static_cast<std::uint32_t>(bytes[offset+2])<<16) | (static_cast<std::uint32_t>(bytes[offset+3])<<24);
	}
	static std::uint16_t ReadLe16(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		return static_cast<std::uint16_t>(bytes[offset]) | static_cast<std::uint16_t>(bytes[offset+1]<<8);
	}
	static float ReadLeFloat(const std::vector<std::uint8_t>& bytes, std::size_t offset)
	{
		const auto bits = ReadLe32(bytes, offset); float value = 0.0F;
		static_assert(sizeof(bits) == sizeof(value), "unsupported float representation");
		static_assert(std::numeric_limits<float>::is_iec559, "binary STL requires IEEE-754 float32");
		std::memcpy(&value, &bits, sizeof(value)); return value;
	}
	static bool IsExactBinaryStl(const std::vector<std::uint8_t>& bytes)
	{
		constexpr std::size_t prefix = 84U, record = 50U;
		if (bytes.size() < prefix) return false;
		const auto count = static_cast<std::size_t>(ReadLe32(bytes, 80));
		return count <= (std::numeric_limits<std::size_t>::max()-prefix)/record
			&& bytes.size() == prefix+record*count;
	}
	static void ValidateCount(std::size_t count, const SurfaceValidationOptions& options)
	{
		if (options.max_triangles < 0) throw std::invalid_argument("configured maximum triangle count is negative");
		if (static_cast<std::uintmax_t>(count) > static_cast<std::uintmax_t>(options.max_triangles))
			throw std::invalid_argument("triangle count exceeds configured maximum");
	}
	static RawSurfaceSoup ParseBinary(const std::vector<std::uint8_t>& bytes, const SurfaceValidationOptions& options)
	{
		const auto count = static_cast<std::size_t>(ReadLe32(bytes, 80)); ValidateCount(count, options);
		if (count > std::numeric_limits<std::size_t>::max()/3U) throw std::overflow_error("vertex count overflows");
		RawSurfaceSoup soup; soup.vertices.reserve(3*count); soup.triangles.reserve(count);
		for (std::size_t face = 0, offset = 84; face < count; ++face, offset += 50) {
			for (std::size_t component = 0; component < 3; ++component)
				if (!std::isfinite(ReadLeFloat(bytes, offset+4*component))) throw std::invalid_argument("binary facet normal is nonfinite");
			RawSurfaceTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) {
				std::array<double, 3> vertex{};
				for (std::size_t component = 0; component < 3; ++component) {
					const float value = ReadLeFloat(bytes, offset+12+12*corner+4*component);
					if (!std::isfinite(value)) throw std::invalid_argument("binary vertex is nonfinite");
					vertex[component] = value;
				}
				triangle.indices[corner] = static_cast<std::int64_t>(soup.vertices.size()); soup.vertices.push_back(vertex);
			}
			if (ReadLe16(bytes, offset+48) != 0) throw std::invalid_argument("binary attribute word must be zero");
			soup.triangles.push_back(triangle);
		}
		return soup;
	}
	static double ParseFinite(const std::string& token)
	{
		errno = 0; char* end = nullptr; const double value = std::strtod(token.c_str(), &end);
		if (end != token.c_str()+token.size() || errno == ERANGE || !std::isfinite(value))
			throw std::invalid_argument("ASCII numeric token is invalid");
		return value;
	}
	static RawSurfaceSoup ParseAscii(const std::vector<std::uint8_t>& bytes, const SurfaceValidationOptions& options)
	{
		if (bytes.size() > 64U*1024U*1024U) throw std::invalid_argument("ASCII STL exceeds parser limit");
		const std::string text(bytes.begin(), bytes.end()); std::istringstream input(text); std::vector<std::string> tokens;
		for (std::string token; input >> token;) { if (tokens.size() == 10000000U) throw std::invalid_argument("ASCII STL has too many tokens"); tokens.push_back(token); }
		std::size_t at = 0;
		auto need = [&]() -> const std::string& { if (at == tokens.size()) throw std::invalid_argument("ASCII STL is truncated"); return tokens[at++]; };
		auto expect = [&](const char* word) { if (need() != word) throw std::invalid_argument(std::string("expected '")+word+"'"); };
		expect("solid");
		if (at < tokens.size() && tokens[at] != "facet" && tokens[at] != "endsolid") ++at;
		RawSurfaceSoup soup;
		while (at < tokens.size() && tokens[at] == "facet") {
			expect("facet"); expect("normal"); for (int i = 0; i < 3; ++i) (void)ParseFinite(need());
			expect("outer"); expect("loop"); RawSurfaceTriangle triangle;
			for (std::size_t corner = 0; corner < 3; ++corner) {
				expect("vertex"); std::array<double, 3> vertex{{ParseFinite(need()), ParseFinite(need()), ParseFinite(need())}};
				triangle.indices[corner] = static_cast<std::int64_t>(soup.vertices.size()); soup.vertices.push_back(vertex);
			}
			expect("endloop"); expect("endfacet"); soup.triangles.push_back(triangle);
			ValidateCount(soup.triangles.size(), options);
		}
		expect("endsolid");
		if (at < tokens.size()) ++at;
		if (at != tokens.size()) throw std::invalid_argument("ASCII STL has trailing tokens");
		return soup;
	}
};

} // namespace iga

#endif

#ifndef IGA_CASE_CONFIG_HPP
#define IGA_CASE_CONFIG_HPP

#include <array>
#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iga {

enum class BoundaryType {
	Interior,
	Wall,
	Inlet,
	Outlet
};

struct BoundaryRule {
	int label = -1;
	std::string name;
	BoundaryType type = BoundaryType::Interior;
	bool has_velocity = false;
	std::array<double, 3> velocity{{0.0, 0.0, 0.0}};
	bool has_velocity_scale = false;
	double velocity_scale = 1.0;
	bool has_pressure = false;
	double pressure = 0.0;
	bool has_transport = false;
	double n0 = 0.0;
	double nplus = 0.0;
};

struct CaseConfiguration {
	bool present = false;
	bool inherit_legacy = true;
	std::vector<BoundaryRule> boundaries;
};

namespace config_detail {

struct JsonValue {
	enum class Type { Null, Boolean, Number, String, Array, Object };

	Type type = Type::Null;
	bool boolean = false;
	double number = 0.0;
	std::string string;
	std::vector<JsonValue> array;
	std::map<std::string, JsonValue> object;
};

class JsonParser {
public:
	explicit JsonParser(std::string input) : input_(std::move(input)) {}

	JsonValue Parse()
	{
		SkipWhitespace();
		auto value = ParseValue();
		SkipWhitespace();
		if (position_ != input_.size()) Fail("unexpected trailing JSON content");
		return value;
	}

private:
	[[noreturn]] void Fail(const std::string& message) const
	{
		throw std::runtime_error("case_config.json: " + message + " at byte "
			+ std::to_string(position_));
	}

	void SkipWhitespace()
	{
		while (position_ < input_.size()) {
			const char c = input_[position_];
			if (c != ' ' && c != '\t' && c != '\n' && c != '\r') break;
			++position_;
		}
	}

	bool Consume(char expected)
	{
		SkipWhitespace();
		if (position_ >= input_.size() || input_[position_] != expected) return false;
		++position_;
		return true;
	}

	JsonValue ParseValue()
	{
		SkipWhitespace();
		if (position_ >= input_.size()) Fail("expected a JSON value");
		const char c = input_[position_];
		if (c == '{') return ParseObject();
		if (c == '[') return ParseArray();
		if (c == '"') {
			JsonValue value;
			value.type = JsonValue::Type::String;
			value.string = ParseString();
			return value;
		}
		if (c == '-' || (c >= '0' && c <= '9')) return ParseNumber();
		if (MatchLiteral("true")) {
			JsonValue value; value.type = JsonValue::Type::Boolean; value.boolean = true; return value;
		}
		if (MatchLiteral("false")) {
			JsonValue value; value.type = JsonValue::Type::Boolean; value.boolean = false; return value;
		}
		if (MatchLiteral("null")) return {};
		Fail("invalid JSON value");
	}

	JsonValue ParseObject()
	{
		if (!Consume('{')) Fail("expected '{'");
		JsonValue result;
		result.type = JsonValue::Type::Object;
		if (Consume('}')) return result;
		while (true) {
			SkipWhitespace();
			if (position_ >= input_.size() || input_[position_] != '"')
				Fail("expected an object key");
			const auto key = ParseString();
			if (!Consume(':')) Fail("expected ':' after object key");
			auto inserted = result.object.emplace(key, ParseValue());
			if (!inserted.second) Fail("duplicate object key '" + key + "'");
			if (Consume('}')) break;
			if (!Consume(',')) Fail("expected ',' or '}' in object");
		}
		return result;
	}

	JsonValue ParseArray()
	{
		if (!Consume('[')) Fail("expected '['");
		JsonValue result;
		result.type = JsonValue::Type::Array;
		if (Consume(']')) return result;
		while (true) {
			result.array.push_back(ParseValue());
			if (Consume(']')) break;
			if (!Consume(',')) Fail("expected ',' or ']' in array");
		}
		return result;
	}

	std::string ParseString()
	{
		if (position_ >= input_.size() || input_[position_] != '"') Fail("expected string");
		++position_;
		std::string result;
		while (position_ < input_.size()) {
			const char c = input_[position_++];
			if (c == '"') return result;
			if (static_cast<unsigned char>(c) < 0x20) Fail("control character in string");
			if (c != '\\') {
				result.push_back(c);
				continue;
			}
			if (position_ >= input_.size()) Fail("unterminated escape sequence");
			const char escaped = input_[position_++];
			switch (escaped) {
			case '"': result.push_back('"'); break;
			case '\\': result.push_back('\\'); break;
			case '/': result.push_back('/'); break;
			case 'b': result.push_back('\b'); break;
			case 'f': result.push_back('\f'); break;
			case 'n': result.push_back('\n'); break;
			case 'r': result.push_back('\r'); break;
			case 't': result.push_back('\t'); break;
			default: Fail("unsupported string escape");
			}
		}
		Fail("unterminated string");
	}

	JsonValue ParseNumber()
	{
		const char* first = input_.c_str() + position_;
		char* last = nullptr;
		errno = 0;
		const double number = std::strtod(first, &last);
		if (last == first || errno == ERANGE || !std::isfinite(number)) Fail("invalid number");
		position_ += static_cast<std::size_t>(last - first);
		JsonValue result;
		result.type = JsonValue::Type::Number;
		result.number = number;
		return result;
	}

	bool MatchLiteral(const std::string& literal)
	{
		if (input_.compare(position_, literal.size(), literal) != 0) return false;
		position_ += literal.size();
		return true;
	}

	std::string input_;
	std::size_t position_ = 0;
};

inline const std::map<std::string, JsonValue>& RequireObject(
	const JsonValue& value, const std::string& context)
{
	if (value.type != JsonValue::Type::Object)
		throw std::runtime_error("case_config.json: " + context + " must be an object");
	return value.object;
}

inline const std::vector<JsonValue>& RequireArray(
	const JsonValue& value, const std::string& context)
{
	if (value.type != JsonValue::Type::Array)
		throw std::runtime_error("case_config.json: " + context + " must be an array");
	return value.array;
}

inline double RequireNumber(const JsonValue& value, const std::string& context)
{
	if (value.type != JsonValue::Type::Number || !std::isfinite(value.number))
		throw std::runtime_error("case_config.json: " + context + " must be a finite number");
	return value.number;
}

inline std::string RequireString(const JsonValue& value, const std::string& context)
{
	if (value.type != JsonValue::Type::String)
		throw std::runtime_error("case_config.json: " + context + " must be a string");
	return value.string;
}

inline bool RequireBoolean(const JsonValue& value, const std::string& context)
{
	if (value.type != JsonValue::Type::Boolean)
		throw std::runtime_error("case_config.json: " + context + " must be a boolean");
	return value.boolean;
}

inline void RequireKnownKeys(const std::map<std::string, JsonValue>& object,
	const std::set<std::string>& allowed, const std::string& context)
{
	for (const auto& item : object)
		if (!allowed.count(item.first))
			throw std::runtime_error("case_config.json: unknown " + context + " key '" + item.first + "'");
}

inline const JsonValue* Find(const std::map<std::string, JsonValue>& object, const std::string& key)
{
	const auto it = object.find(key);
	return it == object.end() ? nullptr : &it->second;
}

inline int RequireInteger(const JsonValue& value, const std::string& context)
{
	const double number = RequireNumber(value, context);
	if (std::floor(number) != number || number < 0.0 || number > 2147483647.0)
		throw std::runtime_error("case_config.json: " + context + " must be a non-negative integer");
	return static_cast<int>(number);
}

inline BoundaryType ParseBoundaryType(const std::string& type)
{
	if (type == "interior") return BoundaryType::Interior;
	if (type == "wall") return BoundaryType::Wall;
	if (type == "inlet") return BoundaryType::Inlet;
	if (type == "outlet") return BoundaryType::Outlet;
	throw std::runtime_error("case_config.json: unsupported boundary type '" + type + "'");
}

inline BoundaryRule ParseBoundaryRule(const JsonValue& value, std::size_t index)
{
	const std::string context = "boundaries.conditions[" + std::to_string(index) + "]";
	const auto& object = RequireObject(value, context);
	RequireKnownKeys(object, {"label", "name", "type", "velocity", "velocity_scale", "pressure", "transport"}, context);
	const auto* label = Find(object, "label");
	const auto* type = Find(object, "type");
	if (!label || !type) throw std::runtime_error("case_config.json: " + context + " requires label and type");

	BoundaryRule rule;
	rule.label = RequireInteger(*label, context + ".label");
	rule.type = ParseBoundaryType(RequireString(*type, context + ".type"));
	if (const auto* name = Find(object, "name")) rule.name = RequireString(*name, context + ".name");
	if (const auto* scale = Find(object, "velocity_scale")) {
		rule.has_velocity_scale = true;
		rule.velocity_scale = RequireNumber(*scale, context + ".velocity_scale");
	}
	if (const auto* velocity = Find(object, "velocity")) {
		const auto& values = RequireArray(*velocity, context + ".velocity");
		if (values.size() != 3) throw std::runtime_error("case_config.json: " + context + ".velocity must have three entries");
		rule.has_velocity = true;
		for (std::size_t i = 0; i < 3; ++i)
			rule.velocity[i] = RequireNumber(values[i], context + ".velocity[" + std::to_string(i) + "]");
	}
	if (const auto* pressure = Find(object, "pressure")) {
		rule.has_pressure = true;
		rule.pressure = RequireNumber(*pressure, context + ".pressure");
	}
	if (const auto* transport = Find(object, "transport")) {
		const auto& values = RequireObject(*transport, context + ".transport");
		RequireKnownKeys(values, {"N0", "Nplus"}, context + ".transport");
		const auto* n0 = Find(values, "N0");
		const auto* nplus = Find(values, "Nplus");
		if (!n0 || !nplus)
			throw std::runtime_error("case_config.json: " + context + ".transport requires N0 and Nplus");
		rule.has_transport = true;
		rule.n0 = RequireNumber(*n0, context + ".transport.N0");
		rule.nplus = RequireNumber(*nplus, context + ".transport.Nplus");
	}

	if (rule.has_velocity && rule.has_velocity_scale)
		throw std::runtime_error("case_config.json: " + context + " cannot set both velocity and velocity_scale");
	if (rule.type == BoundaryType::Interior) {
		if (rule.has_velocity || rule.has_velocity_scale || rule.has_pressure || rule.has_transport)
			throw std::runtime_error("case_config.json: interior boundary rules cannot set values");
	} else if (rule.type == BoundaryType::Wall) {
		if (rule.has_velocity_scale || rule.has_pressure || rule.has_transport)
			throw std::runtime_error("case_config.json: wall rules only support an optional velocity");
	} else if (rule.type == BoundaryType::Inlet) {
		if (rule.has_pressure)
			throw std::runtime_error("case_config.json: inlet rules cannot set pressure");
	} else if (rule.type == BoundaryType::Outlet) {
		if (rule.has_velocity || rule.has_velocity_scale || rule.has_transport)
			throw std::runtime_error("case_config.json: outlet rules only support an optional pressure");
	}
	return rule;
}

} // namespace config_detail

inline CaseConfiguration ParseCaseConfiguration(const std::string& text)
{
	using namespace config_detail;
	const auto root_value = JsonParser(text).Parse();
	const auto& root = RequireObject(root_value, "root");
	RequireKnownKeys(root, {"schema_version", "boundaries"}, "root");
	if (const auto* version = Find(root, "schema_version")) {
		if (RequireInteger(*version, "schema_version") != 1)
			throw std::runtime_error("case_config.json: only schema_version 1 is supported");
	}
	const auto* boundaries = Find(root, "boundaries");
	if (!boundaries) throw std::runtime_error("case_config.json: missing boundaries object");
	const auto& boundary_object = RequireObject(*boundaries, "boundaries");
	RequireKnownKeys(boundary_object, {"inherit_legacy", "conditions"}, "boundaries");

	CaseConfiguration configuration;
	configuration.present = true;
	if (const auto* inherit = Find(boundary_object, "inherit_legacy"))
		configuration.inherit_legacy = RequireBoolean(*inherit, "boundaries.inherit_legacy");
	const auto* conditions = Find(boundary_object, "conditions");
	if (!conditions) throw std::runtime_error("case_config.json: boundaries.conditions is required");
	const auto& rules = RequireArray(*conditions, "boundaries.conditions");
	std::set<int> labels;
	for (std::size_t i = 0; i < rules.size(); ++i) {
		auto rule = ParseBoundaryRule(rules[i], i);
		if (!labels.insert(rule.label).second)
			throw std::runtime_error("case_config.json: duplicate boundary label " + std::to_string(rule.label));
		configuration.boundaries.push_back(std::move(rule));
	}
	return configuration;
}

inline CaseConfiguration ReadCaseConfiguration(const std::string& path)
{
	std::ifstream input(path);
	if (!input) {
		std::error_code error;
		const bool exists = std::filesystem::exists(path, error);
		if (error) throw std::runtime_error("cannot inspect case configuration: " + path
			+ ": " + error.message());
		if (!exists) return {};
		throw std::runtime_error("cannot open case configuration: " + path);
	}
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error("cannot read case configuration: " + path);
	return ParseCaseConfiguration(contents.str());
}

inline const char* BoundaryTypeName(BoundaryType type)
{
	switch (type) {
	case BoundaryType::Interior: return "interior";
	case BoundaryType::Wall: return "wall";
	case BoundaryType::Inlet: return "inlet";
	case BoundaryType::Outlet: return "outlet";
	}
	return "unknown";
}

} // namespace iga

#endif

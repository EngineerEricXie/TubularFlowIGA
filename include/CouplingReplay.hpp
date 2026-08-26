#ifndef IGA_COUPLING_REPLAY_HPP
#define IGA_COUPLING_REPLAY_HPP

#include "CaseConfig.hpp"
#include "VascularCoupling.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

namespace coupling_replay_detail {

inline std::string TrimReplayText(std::string value)
{
	while (!value.empty() && (value.front() == ' ' || value.front() == '\t'
		|| value.front() == '\r')) value.erase(value.begin());
	while (!value.empty() && (value.back() == ' ' || value.back() == '\t'
		|| value.back() == '\r')) value.pop_back();
	return value;
}

inline std::vector<std::string> SplitReplayCsv(const std::string& line)
{
	std::vector<std::string> result;
	std::string current;
	bool quoted = false;
	for (std::size_t i = 0; i < line.size(); ++i) {
		const char character = line[i];
		if (character == '"') {
			if (quoted && i+1 < line.size() && line[i+1] == '"') {
				current.push_back('"');
				++i;
			} else quoted = !quoted;
		} else if (character == ',' && !quoted) {
			result.push_back(TrimReplayText(current));
			current.clear();
		} else current.push_back(character);
	}
	if (quoted) throw std::runtime_error("unterminated quote in coupling replay CSV");
	result.push_back(TrimReplayText(current));
	return result;
}

inline double ReplayNumber(const std::string& text, const std::string& context)
{
	std::size_t used = 0;
	double result = 0.0;
	try { result = std::stod(text, &used); }
	catch (const std::exception&) {
		throw std::runtime_error("invalid coupling replay number for "+context);
	}
	if (used != text.size() || !std::isfinite(result))
		throw std::runtime_error("invalid coupling replay number for "+context);
	return result;
}

inline VascularInletState ParseReplayJsonState(
	const config_detail::JsonValue& value, std::size_t index)
{
	using namespace config_detail;
	const std::string context = "coupling replay states["+std::to_string(index)+"]";
	const auto& object = RequireObject(value, context);
	RequireKnownKeys(object, {"time_s", "flow_m3_s", "pressure_pa", "pressure_Pa",
		"species", "temperature_c", "temperature_C", "hematocrit_percent"}, context);
	VascularInletState result;
	const auto* time = Find(object, "time_s");
	if (!time) throw std::runtime_error(context+" requires time_s");
	result.time_s = RequireNumber(*time, context+".time_s");
	if (const auto* flow = Find(object, "flow_m3_s")) {
		result.has_flow = true;
		result.flow_m3_s = RequireNumber(*flow, context+".flow_m3_s");
	}
	const auto* pressure = Find(object, "pressure_pa");
	if (!pressure) pressure = Find(object, "pressure_Pa");
	if (pressure) {
		result.has_pressure = true;
		result.pressure_pa = RequireNumber(*pressure, context+".pressure_pa");
	}
	if (const auto* species = Find(object, "species")) {
		const auto& items = RequireObject(*species, context+".species");
		for (const auto& item : items)
			result.species.emplace(item.first,
				RequireNumber(item.second, context+".species."+item.first));
	}
	const auto* temperature = Find(object, "temperature_c");
	if (!temperature) temperature = Find(object, "temperature_C");
	if (temperature) {
		result.has_temperature = true;
		result.temperature_c = RequireNumber(*temperature, context+".temperature_c");
	}
	if (const auto* hematocrit = Find(object, "hematocrit_percent")) {
		result.has_hematocrit = true;
		result.hematocrit_percent = RequireNumber(*hematocrit,
			context+".hematocrit_percent");
	}
	ValidateVascularInletState(result);
	return result;
}

inline std::vector<VascularInletState> ReadReplayJson(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open coupling replay file: "+path.string());
	std::ostringstream text;
	text << input.rdbuf();
	const auto root = config_detail::JsonParser(text.str()).Parse();
	const std::vector<config_detail::JsonValue>* rows = nullptr;
	if (root.type == config_detail::JsonValue::Type::Array) rows = &root.array;
	else {
		const auto& object = config_detail::RequireObject(root, "coupling replay root");
		config_detail::RequireKnownKeys(object, {"states"}, "coupling replay root");
		const auto* states = config_detail::Find(object, "states");
		if (!states) throw std::runtime_error("coupling replay JSON requires states");
		rows = &config_detail::RequireArray(*states, "coupling replay states");
	}
	std::vector<VascularInletState> result;
	for (std::size_t i = 0; i < rows->size(); ++i)
		result.push_back(ParseReplayJsonState((*rows)[i], i));
	return result;
}

inline std::vector<VascularInletState> ReadReplayCsv(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open coupling replay file: "+path.string());
	std::string line;
	if (!std::getline(input, line)) throw std::runtime_error("coupling replay CSV is empty");
	const auto header = SplitReplayCsv(line);
	std::map<std::string, std::size_t> columns;
	for (std::size_t i = 0; i < header.size(); ++i)
		if (!columns.emplace(header[i], i).second)
			throw std::runtime_error("coupling replay CSV has duplicate columns");
	if (!columns.count("time_s"))
		throw std::runtime_error("coupling replay CSV requires time_s");
	std::vector<VascularInletState> result;
	std::size_t row = 1;
	while (std::getline(input, line)) {
		++row;
		if (TrimReplayText(line).empty()) continue;
		const auto values = SplitReplayCsv(line);
		if (values.size() != header.size())
			throw std::runtime_error("coupling replay CSV row has the wrong column count");
		auto optional = [&](const std::string& name, double& target) {
			const auto found = columns.find(name);
			if (found == columns.end() || values[found->second].empty()) return false;
			target = ReplayNumber(values[found->second], name+" row "+std::to_string(row));
			return true;
		};
		VascularInletState state;
		state.time_s = ReplayNumber(values[columns.at("time_s")], "time_s");
		state.has_flow = optional("flow_m3_s", state.flow_m3_s);
		state.has_pressure = optional("pressure_pa", state.pressure_pa)
			|| optional("pressure_Pa", state.pressure_pa);
		state.has_temperature = optional("temperature_c", state.temperature_c)
			|| optional("temperature_C", state.temperature_c);
		state.has_hematocrit = optional("hematocrit_percent", state.hematocrit_percent);
		for (const auto& column : columns) {
			const std::string prefix = "species.";
			if (column.first.compare(0, prefix.size(), prefix) != 0) continue;
			if (!values[column.second].empty())
				state.species.emplace(column.first.substr(prefix.size()),
					ReplayNumber(values[column.second], column.first));
		}
		ValidateVascularInletState(state);
		result.push_back(std::move(state));
	}
	return result;
}

} // namespace coupling_replay_detail

class ReplayInletProvider {
public:
	explicit ReplayInletProvider(std::vector<VascularInletState> states)
		: states_(std::move(states))
	{
		if (states_.empty()) throw std::runtime_error("coupling replay requires a state");
		std::sort(states_.begin(), states_.end(), [](const auto& first, const auto& second) {
			return first.time_s < second.time_s;
		});
		for (std::size_t i = 1; i < states_.size(); ++i)
			if (states_[i-1].time_s == states_[i].time_s)
				throw std::runtime_error("coupling replay times must be unique");
	}

	static ReplayInletProvider Read(const std::filesystem::path& path)
	{
		const auto extension = path.extension().string();
		if (extension == ".json")
			return ReplayInletProvider(coupling_replay_detail::ReadReplayJson(path));
		if (extension == ".csv")
			return ReplayInletProvider(coupling_replay_detail::ReadReplayCsv(path));
		throw std::runtime_error("coupling replay file must use .json or .csv");
	}

	VascularInletState StateAt(double time_s) const
	{
		const VascularInletState* selected = &states_.front();
		for (const auto& state : states_) {
			if (state.time_s > time_s) break;
			selected = &state;
		}
		auto result = *selected;
		result.time_s = time_s;
		return result;
	}

private:
	std::vector<VascularInletState> states_;
};

} // namespace iga

#endif

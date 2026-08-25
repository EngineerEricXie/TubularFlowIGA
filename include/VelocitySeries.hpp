#ifndef IGA_VELOCITY_SERIES_HPP
#define IGA_VELOCITY_SERIES_HPP

#include "SimulationConfig.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct VelocitySnapshot {
	double time = 0.0;
	std::filesystem::path file;
};

struct VelocityInterpolation {
	std::size_t lower = 0;
	std::size_t upper = 0;
	double upper_weight = 0.0;
};

inline std::string TrimVelocityToken(const std::string& text)
{
	const auto first = text.find_first_not_of(" \t\r");
	if (first == std::string::npos) return {};
	const auto last = text.find_last_not_of(" \t\r");
	return text.substr(first, last-first+1);
}

inline std::vector<VelocitySnapshot> ParseVelocityManifest(
	const std::string& text, const std::string& context = "velocity manifest")
{
	std::vector<VelocitySnapshot> snapshots;
	std::istringstream input(text);
	std::string line;
	std::size_t line_number = 0;
	while (std::getline(input, line)) {
		++line_number;
		line = TrimVelocityToken(line);
		if (line.empty() || line.front() == '#') continue;
		const auto comma = line.find(',');
		if (comma == std::string::npos || line.find(',', comma+1) != std::string::npos)
			throw std::runtime_error(context + " line " + std::to_string(line_number)
				+ " must contain exactly time,file");
		const auto time_token = TrimVelocityToken(line.substr(0, comma));
		const auto file_token = TrimVelocityToken(line.substr(comma+1));
		if (snapshots.empty() && time_token == "time" && file_token == "file") continue;
		std::size_t used = 0;
		double time = 0.0;
		try {
			time = std::stod(time_token, &used);
		} catch (const std::exception&) {
			throw std::runtime_error(context + " line " + std::to_string(line_number)
				+ " has an invalid time");
		}
		if (used != time_token.size() || !std::isfinite(time) || time < 0.0)
			throw std::runtime_error(context + " times must be finite and nonnegative");
		if (!snapshots.empty() && !(time > snapshots.back().time))
			throw std::runtime_error(context + " times must be strictly increasing");
		const std::filesystem::path file(file_token);
		if (file.empty() || file.is_absolute())
			throw std::runtime_error(context + " files must be non-empty paths relative to the case directory");
		snapshots.push_back({time, file});
	}
	if (snapshots.empty()) throw std::runtime_error(context + " contains no snapshots");
	return snapshots;
}

inline std::vector<VelocitySnapshot> ReadVelocityManifest(const std::filesystem::path& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open velocity manifest: " + path.string());
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof())
		throw std::runtime_error("cannot read velocity manifest: " + path.string());
	return ParseVelocityManifest(contents.str(), path.string());
}

inline VelocityInterpolation ResolveVelocityInterpolation(
	const std::vector<VelocitySnapshot>& snapshots, double time, const std::string& out_of_range)
{
	if (snapshots.empty() || !std::isfinite(time))
		throw std::runtime_error("cannot interpolate an empty velocity series or non-finite time");
	if (time < snapshots.front().time) {
		if (out_of_range == "hold") return {0, 0, 0.0};
		throw std::runtime_error("transport time precedes the first velocity snapshot");
	}
	if (time > snapshots.back().time) {
		if (out_of_range == "hold") return {snapshots.size()-1, snapshots.size()-1, 0.0};
		throw std::runtime_error("transport time exceeds the last velocity snapshot");
	}
	const auto upper = std::lower_bound(snapshots.begin(), snapshots.end(), time,
		[](const VelocitySnapshot& snapshot, double value) { return snapshot.time < value; });
	if (upper == snapshots.end()) return {snapshots.size()-1, snapshots.size()-1, 0.0};
	const auto upper_index = static_cast<std::size_t>(upper-snapshots.begin());
	if (upper->time == time || upper_index == 0) return {upper_index, upper_index, 0.0};
	const auto lower_index = upper_index-1;
	const auto weight = (time-snapshots[lower_index].time)
		/(snapshots[upper_index].time-snapshots[lower_index].time);
	return {lower_index, upper_index, weight};
}

inline const VelocitySourceDefinition& FindVelocitySource(
	const SimulationConfiguration& configuration, const std::string& name)
{
	const auto found = std::find_if(configuration.velocity_sources.begin(),
		configuration.velocity_sources.end(),
		[&](const VelocitySourceDefinition& source) { return source.name == name; });
	if (found == configuration.velocity_sources.end())
		throw std::runtime_error("unknown velocity source '" + name + "'");
	return *found;
}

} // namespace iga

#endif

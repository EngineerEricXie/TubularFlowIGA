#ifndef IGA_TEMPORAL_FUNCTION_HPP
#define IGA_TEMPORAL_FUNCTION_HPP

#include "SimulationConfig.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct TemporalSample {
	double time = 0.0;
	double value = 0.0;
};

inline std::string TrimTemporalToken(const std::string& text)
{
	const auto first = text.find_first_not_of(" \t\r");
	if (first == std::string::npos) return {};
	const auto last = text.find_last_not_of(" \t\r");
	return text.substr(first, last-first+1);
}

inline double ParseTemporalNumber(const std::string& text, const std::string& context)
{
	const auto token = TrimTemporalToken(text);
	std::size_t used = 0;
	double value = 0.0;
	try {
		value = std::stod(token, &used);
	} catch (const std::exception&) {
		throw std::runtime_error(context + " contains a non-numeric value");
	}
	if (used != token.size() || !std::isfinite(value))
		throw std::runtime_error(context + " contains a non-finite or malformed value");
	return value;
}

inline std::vector<TemporalSample> ParseTemporalCsv(
	const std::string& text, double period, const std::string& context = "temporal CSV")
{
	if (!(period > 0.0)) throw std::runtime_error(context + " requires a positive period");
	std::vector<TemporalSample> samples;
	std::istringstream input(text);
	std::string line;
	std::size_t line_number = 0;
	while (std::getline(input, line)) {
		++line_number;
		line = TrimTemporalToken(line);
		if (line.empty() || line.front() == '#') continue;
		const auto comma = line.find(',');
		if (comma == std::string::npos || line.find(',', comma+1) != std::string::npos) {
			throw std::runtime_error(context + " line " + std::to_string(line_number)
				+ " must contain exactly time,value");
		}
		const auto first = TrimTemporalToken(line.substr(0, comma));
		const auto second = TrimTemporalToken(line.substr(comma+1));
		if (samples.empty() && first == "time" && second == "value") continue;
		const auto time = ParseTemporalNumber(first,
			context + " line " + std::to_string(line_number));
		const auto value = ParseTemporalNumber(second,
			context + " line " + std::to_string(line_number));
		if (time < 0.0 || time >= period)
			throw std::runtime_error(context + " sample times must be in [0, period)");
		if (!samples.empty() && !(time > samples.back().time))
			throw std::runtime_error(context + " sample times must be strictly increasing");
		samples.push_back({time, value});
	}
	if (samples.size() < 2 || samples.front().time != 0.0)
		throw std::runtime_error(context + " requires at least two samples beginning at time zero");
	return samples;
}

inline std::vector<TemporalSample> ReadTemporalCsv(const std::string& path, double period)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open temporal CSV: " + path);
	std::ostringstream contents;
	contents << input.rdbuf();
	if (!input.good() && !input.eof()) throw std::runtime_error("cannot read temporal CSV: " + path);
	return ParseTemporalCsv(contents.str(), period, path);
}

inline double PeriodicTime(double time, double period)
{
	if (!(period > 0.0)) throw std::runtime_error("temporal function period must be positive");
	double wrapped = std::fmod(time, period);
	if (wrapped < 0.0) wrapped += period;
	return wrapped;
}

inline double EvaluateTemporalFunction(const TemporalFunctionDefinition& function,
	double time, const std::vector<TemporalSample>* samples = nullptr)
{
	constexpr double pi = 3.141592653589793238462643383279502884;
	if (function.kind == TemporalFunctionKind::Constant) return function.value;
	const double wrapped = PeriodicTime(time, function.period);
	const double theta = 2.0*pi*wrapped/function.period + function.phase;
	if (function.kind == TemporalFunctionKind::Sinusoid)
		return function.mean + function.amplitude*std::sin(theta);
	if (function.kind == TemporalFunctionKind::Fourier) {
		double value = function.mean;
		for (std::size_t harmonic = 0; harmonic < function.cosine.size(); ++harmonic) {
			const double angle = static_cast<double>(harmonic+1)*theta;
			value += function.cosine[harmonic]*std::cos(angle)
				+ function.sine[harmonic]*std::sin(angle);
		}
		return value;
	}
	if (!samples) throw std::runtime_error("periodic_table evaluation requires loaded samples");
	if (samples->size() < 2 || samples->front().time != 0.0)
		throw std::runtime_error("periodic_table samples are invalid");
	const auto upper = std::upper_bound(samples->begin(), samples->end(), wrapped,
		[](double value, const TemporalSample& sample) { return value < sample.time; });
	const TemporalSample* lower_sample = nullptr;
	TemporalSample upper_sample;
	if (upper == samples->end()) {
		lower_sample = &samples->back();
		upper_sample = {function.period, samples->front().value};
	} else {
		if (upper == samples->begin()) throw std::runtime_error("periodic_table samples must begin at zero");
		lower_sample = &*(upper-1);
		upper_sample = *upper;
	}
	const double fraction = (wrapped-lower_sample->time)/(upper_sample.time-lower_sample->time);
	return lower_sample->value + fraction*(upper_sample.value-lower_sample->value);
}

inline const TemporalFunctionDefinition& FindTemporalFunction(
	const SimulationConfiguration& configuration, const std::string& name)
{
	const auto found = std::find_if(configuration.temporal_functions.begin(),
		configuration.temporal_functions.end(),
		[&](const TemporalFunctionDefinition& function) { return function.name == name; });
	if (found == configuration.temporal_functions.end())
		throw std::runtime_error("unknown temporal function '" + name + "'");
	return *found;
}

inline SimulationConfiguration MaterializeBoundaryWaveforms(
	const SimulationConfiguration& configuration, const std::filesystem::path& case_directory,
	double physical_time)
{
	SimulationConfiguration result = configuration;
	for (auto& boundary : result.boundaries)
		for (auto& condition : boundary.conditions) {
			if (condition.waveform.empty()) continue;
			const auto& function = FindTemporalFunction(configuration, condition.waveform);
			std::vector<TemporalSample> samples;
			const std::vector<TemporalSample>* sample_pointer = nullptr;
			if (function.kind == TemporalFunctionKind::PeriodicTable) {
				const std::filesystem::path relative(function.file);
				if (relative.is_absolute())
					throw std::runtime_error("temporal CSV path must be relative to the case directory");
				samples = ReadTemporalCsv((case_directory/relative).string(), function.period);
				sample_pointer = &samples;
			}
			const double factor = EvaluateTemporalFunction(function, physical_time, sample_pointer);
			if (condition.profile.empty())
				for (auto& value : condition.value) value *= factor;
			else
				condition.scale *= factor;
			condition.waveform.clear();
		}
	return result;
}

} // namespace iga

#endif

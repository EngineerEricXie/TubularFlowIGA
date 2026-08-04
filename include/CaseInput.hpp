#ifndef CASE_INPUT_HPP
#define CASE_INPUT_HPP

#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct TransportParameters {
	double diffusion = 0.0;
	double vplus = 0.0;
	double vminus = 0.0;
	double kplus = 0.0;
	double kminus = 0.0;
	double detach_plus = 0.0;
	double detach_minus = 0.0;
	double dt = 0.0;
	int steps = 0;
	double n0_bc = 0.0;
	double nplus_bc = 0.0;
	double nminus_bc = 0.0;
	double artificial_diffusion = 0.0;
};

inline std::vector<int> ReadPointLabels(const std::string& path, std::uint64_t expected_nodes)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open mesh: " + path);
	std::string token;
	std::uint64_t point_data = 0;
	while (in >> token) {
		if (token == "POINT_DATA") {
			if (!(in >> point_data)) throw std::runtime_error("invalid POINT_DATA record");
			break;
		}
	}
	if (point_data != expected_nodes) throw std::runtime_error("VTK POINT_DATA count does not match database nodes");
	while (in >> token) if (token == "LOOKUP_TABLE") { in >> token; break; }
	if (!in) throw std::runtime_error("VTK point labels were not found");
	std::vector<int> labels(static_cast<std::size_t>(expected_nodes));
	for (auto& label : labels) {
		double value = 0.0;
		if (!(in >> value)) throw std::runtime_error("VTK point-label array is truncated");
		label = static_cast<int>(value);
	}
	return labels;
}

inline std::vector<std::array<double, 3>> ReadVelocity(const std::string& path, std::uint64_t expected_nodes)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open velocity field: " + path);
	std::vector<std::array<double, 3>> velocity(static_cast<std::size_t>(expected_nodes));
	for (auto& value : velocity)
		if (!(in >> value[0] >> value[1] >> value[2]))
			throw std::runtime_error("velocity field is truncated");
	std::string extra;
	if (in >> extra) throw std::runtime_error("velocity field contains extra records");
	return velocity;
}

inline TransportParameters ReadTransportParameters(const std::string& path)
{
	std::ifstream in(path);
	if (!in) throw std::runtime_error("cannot open simulation parameters: " + path);
	std::map<std::string, double> values;
	std::string key;
	double value = 0.0;
	while (in >> key >> value) values[key] = value;
	auto required = [&](const std::string& name) {
		auto it = values.find(name);
		if (it == values.end()) throw std::runtime_error("missing transport parameter: " + name);
		return it->second;
	};
	TransportParameters p;
	p.diffusion = required("D");
	p.vplus = required("vplus");
	p.vminus = required("vminus");
	p.kplus = required("kplus");
	p.kminus = required("kminus");
	p.detach_plus = required("k'plus");
	p.detach_minus = required("k'minus");
	p.dt = required("dt");
	p.steps = static_cast<int>(required("nstep"));
	p.n0_bc = required("N0bc");
	p.nplus_bc = required("Nplusbc");
	p.nminus_bc = required("Nminusbc");
	auto artificial = values.find("artificial_diffusion");
	if (artificial != values.end()) p.artificial_diffusion = artificial->second;
	else {
		artificial = values.find("artificial_diffusion_weight");
		if (artificial != values.end()) p.artificial_diffusion = artificial->second;
	}
	if (p.dt <= 0.0 || p.steps < 0) throw std::runtime_error("dt and nstep must be positive");
	return p;
}

} // namespace iga

#endif

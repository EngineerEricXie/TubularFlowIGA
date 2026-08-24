#include "CaseInput.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void Require(bool condition, const std::string& message)
{
	if (!condition) throw std::runtime_error(message);
}

template <class Function>
void RequireFailure(Function function, const std::string& expected)
{
	try {
		function();
	} catch (const std::exception& error) {
		Require(std::string(error.what()).find(expected) != std::string::npos,
			"unexpected error: " + std::string(error.what()));
		return;
	}
	throw std::runtime_error("expected failure containing: " + expected);
}

}

int main()
{
	try {
		iga::TransportParameters parameters;
		parameters.vplus = 2.0;
		parameters.n0_bc = 3.0;
		parameters.nplus_bc = 4.0;
		const std::vector<int> labels{0, 1, 2, 3, -1};
		const std::vector<std::array<double, 3>> velocity{
			{{9.0, 9.0, 9.0}}, {{1.0, 2.0, 3.0}}, {{0.0, 0.0, 0.0}},
			{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}};

		const auto legacy = iga::ResolveBoundaryConditions({}, labels, velocity, parameters);
		Require(legacy.velocity_nodes == 2, "legacy velocity boundary count");
		Require(legacy.pressure_nodes == 2, "legacy pressure boundary count");
		Require(legacy.transport_nodes == 1, "legacy transport boundary count");
		Require(legacy.velocity[1][0] == 2.0 && legacy.velocity[1][2] == 6.0,
			"legacy inlet velocity scale");

		const auto configured = iga::ParseCaseConfiguration(R"json({
			"schema_version": 1,
			"boundaries": {
				"inherit_legacy": true,
				"conditions": [
					{"label": 0, "type": "wall", "name": "moving", "velocity": [0.1, 0.0, 0.0]},
					{"label": 1, "type": "inlet", "velocity_scale": 0.5,
					 "transport": {"N0": 7.0, "Nplus": 8.0}},
					{"label": 3, "type": "outlet", "pressure": 12.0}
				]
			}
		})json");
		const auto resolved = iga::ResolveBoundaryConditions(configured, labels, velocity, parameters);
		Require(std::abs(resolved.velocity[0][0]-0.1) < 1e-14, "moving wall velocity");
		Require(std::abs(resolved.velocity[1][1]-1.0) < 1e-14, "configured inlet scale");
		Require(resolved.n0[1] == 7.0 && resolved.nplus[1] == 8.0, "configured transport inlet");
		Require(resolved.pressure[2] == 0.0 && resolved.pressure[3] == 12.0,
			"per-outlet pressure");

		RequireFailure([] {
			iga::ParseCaseConfiguration(R"json({"boundaries":{"conditions":[]},"unknown":1})json");
		}, "unknown root key");
		RequireFailure([&] {
			const auto missing = iga::ParseCaseConfiguration(
				R"json({"boundaries":{"inherit_legacy":false,"conditions":[{"label":0,"type":"wall"}]}})json");
			iga::ResolveBoundaryConditions(missing, labels, velocity, parameters);
		}, "has no case_config.json rule");
		RequireFailure([&] {
			const auto absent = iga::ParseCaseConfiguration(
				R"json({"boundaries":{"conditions":[{"label":99,"type":"outlet"}]}})json");
			iga::ResolveBoundaryConditions(absent, labels, velocity, parameters);
		}, "is not present");

		std::cout << "case_config_tests=passed\n";
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "case_config_tests: " << error.what() << '\n';
		return 1;
	}
}

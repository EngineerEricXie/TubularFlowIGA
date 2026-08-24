#include "SimulationConfig.hpp"
#include "GenericCaseInput.hpp"
#include "TemporalFunction.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main()
{
	const std::string json = R"json({
		"schema_version": 2,
		"fields": [
			{"name": "velocity", "kind": "vector3"},
			{"name": "pressure", "kind": "pressure"},
			{"name": "oxygen", "kind": "scalar", "initial_value": 0.2},
			{"name": "drug", "kind": "scalar"}
		],
		"time": {"dt": 0.01, "steps": 5},
		"temporal_functions": [
			{"name": "fixed", "kind": "constant", "units": "1", "value": 0.75},
			{"name": "pulse", "kind": "sinusoid", "units": "1",
			 "mean": 2.0, "amplitude": 3.0, "period": 1.0, "phase": 0.0},
			{"name": "measured", "kind": "periodic_table", "units": "1",
			 "file": "pulse.csv", "period": 1.0, "interpolation": "linear"},
			{"name": "harmonics", "kind": "fourier", "units": "1",
			 "mean": 1.0, "period": 2.0, "cosine": [2.0], "sine": [3.0]}
		],
		"equation_systems": [
			{"name": "flow", "kind": "navier_stokes", "unknowns": ["velocity", "pressure"], "viscosity": 0.0035},
			{"name": "species", "kind": "linear_transport", "unknowns": ["oxygen", "drug"],
			 "terms": [
				{"operator": "time_derivative", "equation": "oxygen"},
				{"operator": "diffusion", "equation": "oxygen", "coefficient": 0.12},
				{"operator": "linear_coupling", "equation": "oxygen", "trial": "drug", "coefficient": -2.0},
				{"operator": "time_derivative", "equation": "drug"},
				{"operator": "advection", "equation": "drug", "velocity": "prescribed"}
			 ],
			 "stabilization": [{"equation": "drug", "method": "supg", "velocity": "prescribed"}]}
		],
		"boundaries": [
			{"label": 0, "name": "vessel_wall", "conditions": [
				{"field": "velocity", "type": "dirichlet", "value": [0, 0, 0], "waveform": "fixed"},
				{"field": "oxygen", "type": "no_flux"},
				{"field": "drug", "type": "robin", "coefficient": 0.5, "exterior_value": 0.1}
			]},
			{"label": 1, "name": "inlet", "conditions": [
				{"field": "oxygen", "type": "dirichlet", "value": 1.0}
			]}
		]
	})json";
	const auto configuration = iga::ParseSimulationConfiguration(json);
	assert(configuration.fields.size() == 4);
	assert(configuration.equation_systems.size() == 2);
	assert(configuration.boundaries[0].name == "vessel_wall");
	assert(configuration.temporal_functions.size() == 4);
	const auto& fixed = iga::FindTemporalFunction(configuration, "fixed");
	const auto& pulse = iga::FindTemporalFunction(configuration, "pulse");
	const auto& measured = iga::FindTemporalFunction(configuration, "measured");
	const auto& harmonics = iga::FindTemporalFunction(configuration, "harmonics");
	assert(std::abs(iga::EvaluateTemporalFunction(fixed, 123.0)-0.75) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(pulse, 0.25)-5.0) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(pulse, 1.25)-5.0) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(harmonics, 0.0)-3.0) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(harmonics, 0.5)-4.0) < 1e-14);
	const auto samples = iga::ParseTemporalCsv(
		" time , value \r\n 0, 1 \r\n0.25,3\r\n0.75,-1\r\n", measured.period);
	assert(std::abs(iga::EvaluateTemporalFunction(measured, 0.125, &samples)-2.0) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(measured, 0.875, &samples)) < 1e-14);
	assert(std::abs(iga::EvaluateTemporalFunction(measured, -0.125, &samples)) < 1e-14);
	auto materialized = configuration;
	materialized.boundaries[0].conditions[0].value = {1.0, 2.0, 3.0};
	materialized = iga::MaterializeBoundaryWaveforms(materialized, ".", 7.0);
	assert(materialized.boundaries[0].conditions[0].waveform.empty());
	assert(std::abs(materialized.boundaries[0].conditions[0].value[0]-0.75) < 1e-14);
	assert(std::abs(materialized.boundaries[0].conditions[0].value[1]-1.5) < 1e-14);
	assert(std::abs(materialized.boundaries[0].conditions[0].value[2]-2.25) < 1e-14);
	const auto compiled = iga::CompileLinearSystem(configuration, "species");
	assert(compiled.fields.size() == 2);
	assert(compiled.field_index.at("oxygen") == 0);
	assert(compiled.field_index.at("drug") == 1);
	assert(compiled.terms.size() == 5);
	assert(compiled.dt == 0.01 && compiled.steps == 5);
	auto waveform_scalar = configuration;
	waveform_scalar.boundaries[1].conditions[0].waveform = "fixed";
	bool waveform_rejected = false;
	try {
		iga::ResolveScalarBoundaries(waveform_scalar, compiled, {0, 1});
	} catch (const std::runtime_error& error) {
		waveform_rejected = std::string(error.what()).find("time-dependent") != std::string::npos;
	}
	assert(waveform_rejected);
	waveform_rejected = false;
	try {
		iga::ResolveFlowBoundaries(configuration,
			iga::FirstNavierStokesSystem(configuration), {0, 1},
			{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
	} catch (const std::runtime_error& error) {
		waveform_rejected = std::string(error.what()).find("time-dependent") != std::string::npos;
	}
	assert(waveform_rejected);
	auto unsupported_flow = configuration;
	unsupported_flow.boundaries[0].conditions[0].waveform.clear();
	unsupported_flow.boundaries[0].conditions.push_back(
		{"pressure", iga::FieldBoundaryKind::Flux, {0.0}, 0.0, 0.0, "", 1.0, ""});
	bool flow_boundary_rejected = false;
	try {
		iga::ResolveFlowBoundaries(unsupported_flow,
			iga::FirstNavierStokesSystem(unsupported_flow), {0, 1},
			{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}});
	} catch (const std::runtime_error&) {
		flow_boundary_rejected = true;
	}
	assert(flow_boundary_rejected);

	bool rejected = false;
	try {
		iga::ParseSimulationConfiguration(R"json({"schema_version":2,"fields":[{"name":"c","kind":"scalar"}],"time":{"dt":1,"steps":1},"equation_systems":[{"name":"s","kind":"linear_transport","unknowns":["c"],"terms":[{"operator":"invented","equation":"c"}]}],"boundaries":[]})json");
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	rejected = false;
	try {
		iga::ParseTemporalCsv("0,1\n0.5,2\n0.4,3\n", 1.0);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "simulation configuration tests passed\n";
}

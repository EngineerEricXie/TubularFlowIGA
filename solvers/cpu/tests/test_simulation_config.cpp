#include "SimulationConfig.hpp"

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

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
				{"field": "velocity", "type": "dirichlet", "value": [0, 0, 0]},
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
	const auto compiled = iga::CompileLinearSystem(configuration, "species");
	assert(compiled.fields.size() == 2);
	assert(compiled.field_index.at("oxygen") == 0);
	assert(compiled.field_index.at("drug") == 1);
	assert(compiled.terms.size() == 5);
	assert(compiled.dt == 0.01 && compiled.steps == 5);

	bool rejected = false;
	try {
		iga::ParseSimulationConfiguration(R"json({"schema_version":2,"fields":[{"name":"c","kind":"scalar"}],"time":{"dt":1,"steps":1},"equation_systems":[{"name":"s","kind":"linear_transport","unknowns":["c"],"terms":[{"operator":"invented","equation":"c"}]}],"boundaries":[]})json");
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "simulation configuration tests passed\n";
}

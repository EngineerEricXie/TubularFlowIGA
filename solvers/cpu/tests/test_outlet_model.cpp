#include "OutletModel.hpp"
#include "OutletCheckpoint.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
	const auto configuration = iga::ParseSimulationConfiguration(R"json({
		"schema_version":2,
		"fields":[{"name":"velocity","kind":"vector3"},{"name":"pressure","kind":"pressure"}],
		"time":{"dt":0.1,"steps":2},
		"equation_systems":[{"name":"flow","kind":"navier_stokes",
			"unknowns":["velocity","pressure"],"viscosity":0.1,"density":1,
			"time_integration":"backward_euler"}],
		"boundaries":[
			{"label":0,"name":"wall","conditions":[{"field":"velocity","type":"dirichlet","value":[0,0,0]}]},
			{"label":2,"name":"outlet","conditions":[{"field":"pressure","type":"windkessel_rcr",
				"proximal_resistance":2,"distal_resistance":4,"capacitance":2,
				"reference_pressure":1,"initial_pressure":5}]}
		]
	})json");
	const auto configured_models = iga::InitializeOutletModels(
		configuration, configuration.equation_systems[0]);
	assert(configured_models.size() == 1);
	assert(configured_models[0].label == 2);
	assert(configured_models[0].capacitor_pressure == 5.0);
	const auto materialized = iga::MaterializeOutletPressures(configuration, configured_models);
	assert(materialized.boundaries[1].conditions[0].kind
		== iga::FieldBoundaryKind::PressureTraction);
	assert(materialized.boundaries[1].conditions[0].value[0] == 5.0);

	iga::OutletModelState resistance;
	resistance.kind = iga::FieldBoundaryKind::Resistance;
	resistance.resistance = 3.0;
	resistance.reference_pressure = 4.0;
	const auto resistance_value = iga::EvaluateOutletModel(resistance, 2.0, 0.0);
	assert(std::abs(resistance_value.pressure-10.0) < 1e-14);

	iga::OutletModelState rc;
	rc.kind = iga::FieldBoundaryKind::WindkesselRC;
	rc.resistance = 4.0;
	rc.capacitance = 2.0;
	rc.reference_pressure = 1.0;
	rc.capacitor_pressure = 5.0;
	const auto rc_value = iga::EvaluateOutletModel(rc, 3.0, 0.5);
	assert(std::abs(rc_value.capacitor_pressure-23.25/4.25) < 1e-14);
	assert(std::abs(rc_value.pressure-rc_value.capacitor_pressure) < 1e-14);

	iga::OutletModelState rcr = rc;
	rcr.kind = iga::FieldBoundaryKind::WindkesselRCR;
	rcr.proximal_resistance = 2.0;
	rcr.distal_resistance = 4.0;
	const auto rcr_value = iga::EvaluateOutletModel(rcr, 3.0, 0.5);
	assert(std::abs(rcr_value.capacitor_pressure-rc_value.capacitor_pressure) < 1e-14);
	assert(std::abs(rcr_value.pressure-(rc_value.capacitor_pressure+6.0)) < 1e-14);

	auto checkpoint_model = configured_models[0];
	checkpoint_model.flow = 2.5;
	checkpoint_model.pressure = 8.0;
	checkpoint_model.capacitor_pressure = 6.0;
	iga::FlowCheckpointMetadata metadata;
	iga::AppendOutletCheckpoint({checkpoint_model}, metadata);
	assert(metadata.outlets.size() == 1);
	auto restored_models = configured_models;
	iga::RestoreOutletCheckpoint(metadata, restored_models);
	assert(restored_models[0].flow == 2.5);
	assert(restored_models[0].pressure == 8.0);
	assert(restored_models[0].capacitor_pressure == 6.0);
	const auto coupled = iga::EvaluateOutletCoupling(
		restored_models, {5.0}, {3.0}, 0.5);
	assert(coupled.flow[0] == 3.0);
	assert(std::abs(coupled.capacitor_pressure[0]-rc_value.capacitor_pressure) < 1e-14);
	assert(std::abs(coupled.pressure[0]-rcr_value.pressure) < 1e-14);
	const auto original_pressure = restored_models[0].pressure;
	iga::RelaxOutletCoupling(restored_models, coupled);
	assert(std::abs(restored_models[0].pressure
		-0.5*(original_pressure+coupled.pressure[0])) < 1e-14);
	iga::CommitOutletCoupling(restored_models, coupled);
	assert(restored_models[0].flow == 3.0);
	assert(restored_models[0].pressure == coupled.pressure[0]);
	assert(restored_models[0].capacitor_pressure == coupled.capacitor_pressure[0]);

	bool rejected = false;
	try {
		iga::EvaluateOutletModel(rc, 3.0, 0.0);
	} catch (const std::runtime_error&) {
		rejected = true;
	}
	assert(rejected);
	std::cout << "outlet model tests passed\n";
}

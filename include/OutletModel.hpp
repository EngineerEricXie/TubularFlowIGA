#ifndef IGA_OUTLET_MODEL_HPP
#define IGA_OUTLET_MODEL_HPP

#include "SimulationConfig.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OutletModelState {
	int label = -1;
	FieldBoundaryKind kind = FieldBoundaryKind::Resistance;
	double resistance = 0.0;
	double proximal_resistance = 0.0;
	double distal_resistance = 0.0;
	double capacitance = 0.0;
	double reference_pressure = 0.0;
	double capacitor_pressure = 0.0;
	double flow = 0.0;
	double pressure = 0.0;
};

struct OutletModelEvaluation {
	double pressure = 0.0;
	double capacitor_pressure = 0.0;
};

struct OutletCouplingEvaluation {
	std::vector<double> flow;
	std::vector<double> pressure;
	std::vector<double> capacitor_pressure;
	double maximum_pressure_change = 0.0;
	double pressure_scale = 1.0;
};

inline std::vector<OutletModelState> InitializeOutletModels(
	const SimulationConfiguration& configuration, const EquationSystemDefinition& system)
{
	if (system.kind != EquationKind::NavierStokes || system.unknowns.size() != 2)
		throw std::runtime_error("outlet models require a Navier-Stokes velocity-pressure system");
	const auto& pressure_name = system.unknowns[1];
	std::vector<OutletModelState> models;
	for (const auto& boundary : configuration.boundaries)
		for (const auto& condition : boundary.conditions) {
			if (condition.field != pressure_name
				|| (condition.kind != FieldBoundaryKind::Resistance
					&& condition.kind != FieldBoundaryKind::WindkesselRC
					&& condition.kind != FieldBoundaryKind::WindkesselRCR)) continue;
			OutletModelState model;
			model.label = boundary.label;
			model.kind = condition.kind;
			model.resistance = condition.resistance;
			model.proximal_resistance = condition.proximal_resistance;
			model.distal_resistance = condition.distal_resistance;
			model.capacitance = condition.capacitance;
			model.reference_pressure = condition.reference_pressure;
			model.capacitor_pressure = condition.initial_pressure;
			model.pressure = condition.kind == FieldBoundaryKind::Resistance
				? condition.reference_pressure : condition.initial_pressure;
			models.push_back(model);
		}
	return models;
}

inline OutletModelEvaluation EvaluateOutletModel(
	const OutletModelState& model, double flow, double dt)
{
	if (!std::isfinite(flow)) throw std::runtime_error("outlet flow must be finite");
	if (model.kind == FieldBoundaryKind::Resistance)
		return {model.reference_pressure+model.resistance*flow,
			model.capacitor_pressure};
	if (!(dt > 0.0)) throw std::runtime_error("RC/RCR outlets require positive transient dt");
	double resistance = 0.0;
	if (model.kind == FieldBoundaryKind::WindkesselRC) resistance = model.resistance;
	else if (model.kind == FieldBoundaryKind::WindkesselRCR) resistance = model.distal_resistance;
	else throw std::runtime_error("unsupported outlet model kind");
	const double capacitor_pressure =
		(model.capacitance*model.capacitor_pressure/dt + flow
			+ model.reference_pressure/resistance)
		/(model.capacitance/dt + 1.0/resistance);
	const double pressure = capacitor_pressure
		+ (model.kind == FieldBoundaryKind::WindkesselRCR
			? model.proximal_resistance*flow : 0.0);
	return {pressure, capacitor_pressure};
}

inline OutletCouplingEvaluation EvaluateOutletCoupling(
	const std::vector<OutletModelState>& models,
	const std::vector<double>& previous_capacitor_pressure,
	const std::vector<double>& flow, double dt)
{
	if (models.size() != previous_capacitor_pressure.size()
		|| models.size() != flow.size())
		throw std::runtime_error("outlet coupling vectors have inconsistent sizes");
	OutletCouplingEvaluation result;
	result.flow = flow;
	result.pressure.resize(models.size());
	result.capacitor_pressure.resize(models.size());
	for (std::size_t i = 0; i < models.size(); ++i) {
		auto previous_model = models[i];
		previous_model.capacitor_pressure = previous_capacitor_pressure[i];
		const auto evaluated = EvaluateOutletModel(previous_model, flow[i], dt);
		result.pressure[i] = evaluated.pressure;
		result.capacitor_pressure[i] = evaluated.capacitor_pressure;
		result.maximum_pressure_change = std::max(result.maximum_pressure_change,
			std::abs(evaluated.pressure-models[i].pressure));
		result.pressure_scale = std::max(result.pressure_scale,
			std::abs(evaluated.pressure));
	}
	return result;
}

inline double OutletCouplingTolerance(const OutletCouplingEvaluation& evaluation)
{
	return 1e-8+1e-6*evaluation.pressure_scale;
}

inline void RelaxOutletCoupling(std::vector<OutletModelState>& models,
	const OutletCouplingEvaluation& evaluation, double relaxation = 0.5)
{
	if (models.size() != evaluation.flow.size()
		|| models.size() != evaluation.pressure.size()
		|| !(relaxation > 0.0 && relaxation <= 1.0))
		throw std::runtime_error("invalid outlet coupling relaxation");
	for (std::size_t i = 0; i < models.size(); ++i) {
		models[i].flow = evaluation.flow[i];
		models[i].pressure = (1.0-relaxation)*models[i].pressure
			+relaxation*evaluation.pressure[i];
	}
}

inline void CommitOutletCoupling(std::vector<OutletModelState>& models,
	const OutletCouplingEvaluation& evaluation)
{
	if (models.size() != evaluation.flow.size()
		|| models.size() != evaluation.pressure.size()
		|| models.size() != evaluation.capacitor_pressure.size())
		throw std::runtime_error("invalid outlet coupling commit");
	for (std::size_t i = 0; i < models.size(); ++i) {
		models[i].flow = evaluation.flow[i];
		models[i].pressure = evaluation.pressure[i];
		models[i].capacitor_pressure = evaluation.capacitor_pressure[i];
	}
}

inline SimulationConfiguration MaterializeOutletPressures(
	const SimulationConfiguration& configuration,
	const std::vector<OutletModelState>& models)
{
	SimulationConfiguration result = configuration;
	std::map<int, double> pressures;
	for (const auto& model : models) pressures.emplace(model.label, model.pressure);
	for (auto& boundary : result.boundaries) {
		const auto pressure = pressures.find(boundary.label);
		if (pressure == pressures.end()) continue;
		for (auto& condition : boundary.conditions)
			if (condition.kind == FieldBoundaryKind::Resistance
				|| condition.kind == FieldBoundaryKind::WindkesselRC
				|| condition.kind == FieldBoundaryKind::WindkesselRCR) {
				condition.kind = FieldBoundaryKind::PressureTraction;
				condition.value = {pressure->second};
			}
	}
	return result;
}

} // namespace iga

#endif

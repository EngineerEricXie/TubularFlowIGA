#ifndef IGA_ONE_D_TRANSPORT_HPP
#define IGA_ONE_D_TRANSPORT_HPP

#include "OneDFlow.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace iga {

struct OneDSpeciesState {
	OneDSpeciesDefinition definition;
	std::vector<double> concentration;
	double inlet_value = 0.0;
	std::string inlet_waveform;
	OneDWallBoundaryKind wall_kind = OneDWallBoundaryKind::NoFlux;
	double wall_value = 0.0;
	double wall_coefficient = 0.0;
	double exterior_value = 0.0;
};

struct OneDTransportState {
	std::string name;
	std::vector<OneDSpeciesState> species;
};

inline double OneDSpeciesSourceIntegral(const OneDConfiguration& configuration,
	const OneDNetwork& network, const OneDFlowState& flow,
	const OneDSpeciesState& species)
{
	double result = 0.0;
	const auto metabolism = configuration.physiology.metabolism_rates.find(
		species.definition.field);
	for (const auto& segment : network.segments) {
		const double dx = segment.length/segment.cells;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
			const double concentration = species.concentration[index];
			double volume_source = species.definition.volume_source
				-species.definition.reaction_rate*concentration;
			if (configuration.physiology.enabled
				&& metabolism != configuration.physiology.metabolism_rates.end())
				volume_source += metabolism->second;
			const double perimeter = 2.0*OneDPi*std::sqrt(flow.area[index]/OneDPi);
			double wall_flux = 0.0;
			if (species.wall_kind == OneDWallBoundaryKind::ConstantFlux)
				wall_flux = species.wall_value;
			else if (species.wall_kind == OneDWallBoundaryKind::Robin)
				wall_flux = species.wall_coefficient*(concentration-species.exterior_value);
			result += dx*(flow.area[index]*volume_source-perimeter*wall_flux);
		}
	}
	return result;
}

inline OneDTransportState InitializeOneDTransport(const OneDConfiguration& configuration,
	const OneDTransportSystemDefinition& transport, const OneDNetwork& network)
{
	OneDTransportState state;
	state.name = transport.name;
	for (const auto& definition : transport.species) {
		OneDSpeciesState species;
		species.definition = definition;
		species.concentration.assign(static_cast<std::size_t>(network.cells),
			FindOneDField(configuration, definition.field).initial_value);
		for (const auto& boundary : configuration.boundaries) {
			for (const auto& condition : boundary.conditions) {
				if (condition.field != definition.field) continue;
				if (boundary.role == "inlet" && condition.type == "dirichlet") {
					species.inlet_value = condition.value;
					species.inlet_waveform = condition.waveform;
				} else if (boundary.role == "wall") {
					if (condition.type == "no_flux") species.wall_kind = OneDWallBoundaryKind::NoFlux;
					else if (condition.type == "constant_flux") {
						species.wall_kind = OneDWallBoundaryKind::ConstantFlux;
						species.wall_value = condition.value;
					} else if (condition.type == "robin") {
						species.wall_kind = OneDWallBoundaryKind::Robin;
						species.wall_coefficient = condition.coefficient;
						species.exterior_value = condition.exterior_value;
					}
				}
			}
		}
		state.species.push_back(std::move(species));
	}
	return state;
}

inline double EvaluateOneDSpeciesInlet(const OneDConfiguration& configuration,
	const OneDSpeciesState& species, const std::filesystem::path& case_directory,
	double time)
{
	if (species.inlet_waveform.empty()) return species.inlet_value;
	const auto& function = FindOneDTemporalFunction(configuration, species.inlet_waveform);
	std::vector<TemporalSample> samples;
	const std::vector<TemporalSample>* pointer = nullptr;
	if (function.kind == TemporalFunctionKind::PeriodicTable) {
		samples = ReadTemporalCsv((case_directory/function.file).string(), function.period);
		pointer = &samples;
	}
	return EvaluateTemporalFunction(function, time, pointer);
}

inline double OneDTransportStableDt(const OneDNetwork& network,
	const OneDFlowState& flow, const OneDSpeciesState& species)
{
	double dt = std::numeric_limits<double>::infinity();
	for (const auto& segment : network.segments) {
		const double dx = segment.length/segment.cells;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
			const double velocity = std::abs(flow.flow[index]/flow.area[index]);
			if (velocity > 0.0) dt = std::min(dt, 0.75*dx/velocity);
		}
		if (species.definition.diffusivity > 0.0)
			dt = std::min(dt, 0.45*dx*dx/species.definition.diffusivity);
	}
	return dt;
}

inline void AdvanceOneDSpecies(const OneDConfiguration& configuration,
	const OneDNetwork& network, const OneDFlowState& flow, OneDSpeciesState& species,
	const std::filesystem::path& case_directory, double start_time, double requested_dt)
{
	double remaining = requested_dt;
	std::vector<double> scalar(species.concentration.size());
	std::vector<double> next(species.concentration.size());
	for (std::size_t i = 0; i < scalar.size(); ++i) scalar[i] = flow.area[i]*species.concentration[i];
	while (remaining > 0.0) {
		const double stable = OneDTransportStableDt(network, flow, species);
		const double dt = std::min(remaining,
			std::isfinite(stable) && stable > 0.0 ? stable : remaining);
		const double elapsed = requested_dt-remaining;
		const double inlet_concentration = EvaluateOneDSpeciesInlet(configuration,
			species, case_directory, start_time+elapsed+dt);
		for (const auto& segment : network.segments) {
			const double dx = segment.length/segment.cells;
			std::vector<double> concentration(static_cast<std::size_t>(segment.cells+2));
			std::vector<double> area(static_cast<std::size_t>(segment.cells+2));
			std::vector<double> q(static_cast<std::size_t>(segment.cells+1));
			for (int cell = 0; cell < segment.cells; ++cell) {
				const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
				concentration[static_cast<std::size_t>(cell+1)] = scalar[index]/flow.area[index];
				area[static_cast<std::size_t>(cell+1)] = flow.area[index];
			}
			if (segment.parent == network.root) concentration[0] = inlet_concentration;
			else {
				const int incoming = OneDSegmentIntoNode(network, segment.parent);
				const auto& parent = network.segments[static_cast<std::size_t>(incoming)];
				concentration[0] = scalar[static_cast<std::size_t>(parent.cell_offset+parent.cells-1)]
					/flow.area[static_cast<std::size_t>(parent.cell_offset+parent.cells-1)];
			}
			concentration[static_cast<std::size_t>(segment.cells+1)] = concentration[static_cast<std::size_t>(segment.cells)];
			area[0] = area[1];
			area[static_cast<std::size_t>(segment.cells+1)] = area[static_cast<std::size_t>(segment.cells)];
			for (int face = 0; face <= segment.cells; ++face) {
				const double face_flow = face == 0
					? flow.flow[static_cast<std::size_t>(segment.cell_offset)]
					: flow.flow[static_cast<std::size_t>(segment.cell_offset+std::min(face, segment.cells)-1)];
				const double upwind = face_flow >= 0.0
					? concentration[static_cast<std::size_t>(face)]
					: concentration[static_cast<std::size_t>(face+1)];
				const double face_area = 0.5*(area[static_cast<std::size_t>(face)]+area[static_cast<std::size_t>(face+1)]);
				q[static_cast<std::size_t>(face)] = face_flow*upwind
					-face_area*species.definition.diffusivity
					*(concentration[static_cast<std::size_t>(face+1)]-concentration[static_cast<std::size_t>(face)])/dx;
			}
			for (int cell = 0; cell < segment.cells; ++cell) {
				const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
				const double c = scalar[index]/flow.area[index];
				double source = species.definition.volume_source-species.definition.reaction_rate*c;
				const auto metabolism = configuration.physiology.metabolism_rates.find(species.definition.field);
				if (configuration.physiology.enabled && metabolism != configuration.physiology.metabolism_rates.end())
					source += metabolism->second;
				const double perimeter = 2.0*OneDPi*std::sqrt(flow.area[index]/OneDPi);
				double wall_flux = 0.0;
				if (species.wall_kind == OneDWallBoundaryKind::ConstantFlux) wall_flux = species.wall_value;
				else if (species.wall_kind == OneDWallBoundaryKind::Robin)
					wall_flux = species.wall_coefficient*(c-species.exterior_value);
				next[index] = scalar[index]-dt/dx*(q[static_cast<std::size_t>(cell+1)]-q[static_cast<std::size_t>(cell)])
					+dt*flow.area[index]*source-dt*perimeter*wall_flux;
				if (!std::isfinite(next[index])) throw std::runtime_error("1d species update produced a non-finite state");
			}
		}
		scalar.swap(next);
		remaining -= dt;
	}
	for (std::size_t i = 0; i < scalar.size(); ++i)
		species.concentration[i] = scalar[i]/flow.area[i];
}

inline void AdvanceOneDTransport(const OneDConfiguration& configuration,
	const OneDNetwork& network, const OneDFlowState& flow, OneDTransportState& transport,
	const std::filesystem::path& case_directory, double start_time, double dt)
{
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) if(transport.species.size() >= 4)
	#endif
	for (long long i = 0; i < static_cast<long long>(transport.species.size()); ++i)
		AdvanceOneDSpecies(configuration, network, flow,
			transport.species[static_cast<std::size_t>(i)], case_directory, start_time, dt);
}

inline const OneDSpeciesState* FindOneDSpecies(const std::vector<OneDTransportState>& transports,
	const std::string& name)
{
	for (const auto& transport : transports)
		for (const auto& species : transport.species)
			if (species.definition.field == name) return &species;
	return nullptr;
}

inline OneDSpeciesState* FindOneDSpecies(std::vector<OneDTransportState>& transports,
	const std::string& name)
{
	for (auto& transport : transports)
		for (auto& species : transport.species)
			if (species.definition.field == name) return &species;
	return nullptr;
}

inline std::map<std::string, std::vector<double>> ComputeOneDDerivedFields(
	const OneDConfiguration& configuration, const std::vector<OneDTransportState>& transports)
{
	std::map<std::string, std::vector<double>> result;
	if (!configuration.physiology.enabled) return result;
	const std::string oxygen_source = configuration.physiology.oxygen_state
		== OxygenTransportState::Total ? "total_oxygen" : "oxygen";
	const auto* oxygen = FindOneDSpecies(transports, oxygen_source);
	const auto* carbon_dioxide = FindOneDSpecies(transports, "carbon_dioxide");
	const auto* bicarbonate = FindOneDSpecies(transports, "bicarbonate");
	OxygenCapacityParameters oxygen_parameters;
	oxygen_parameters.hematocrit_percent = configuration.physiology.hematocrit_percent;
	oxygen_parameters.oxygen_solubility_ml_dl_mmhg
		= configuration.physiology.oxygen_solubility;
	oxygen_parameters.hemoglobin_g_dl = configuration.physiology.hemoglobin_g_dl;
	oxygen_parameters.p50_mmhg = configuration.physiology.p50_mmhg;
	oxygen_parameters.hill_exponent = configuration.physiology.hill_exponent;
	oxygen_parameters.gas_molar_volume_ml_mmol
		= configuration.physiology.gas_molar_volume_ml_mmol;
	std::size_t cells = 0;
	for (const auto& transport : transports)
		if (!transport.species.empty()) { cells = transport.species.front().concentration.size(); break; }
	for (const auto& name : configuration.physiology.derived_fields) {
		std::vector<double> values(cells, 0.0);
		if (name == "pO2" && oxygen) {
			for (std::size_t i = 0; i < cells; ++i)
				values[i] = OxygenFromTransported(oxygen->concentration[i],
					configuration.physiology.oxygen_state, oxygen_parameters).po2_mmhg;
		} else if (name == "pCO2" && carbon_dioxide) {
			for (std::size_t i = 0; i < cells; ++i)
				values[i] = carbon_dioxide->concentration[i]/0.0301;
		} else if (name == "pH" && carbon_dioxide && bicarbonate) {
			for (std::size_t i = 0; i < cells; ++i) {
				const double pco2 = std::max(carbon_dioxide->concentration[i]/0.0301, 1.0e-30);
				values[i] = 6.1+std::log10(std::max(bicarbonate->concentration[i], 1.0e-30)/(0.0301*pco2));
			}
		} else if ((name == "SaO2" || name == "SvO2") && oxygen) {
			for (std::size_t i = 0; i < cells; ++i)
				values[i] = OxygenFromTransported(oxygen->concentration[i],
					configuration.physiology.oxygen_state, oxygen_parameters).saturation;
		} else if (name == "dissolved_oxygen" && oxygen) {
			for (std::size_t i = 0; i < cells; ++i)
				values[i] = OxygenFromTransported(oxygen->concentration[i],
					configuration.physiology.oxygen_state,
					oxygen_parameters).dissolved_oxygen_mol_m3;
		} else if ((name == "bound_oxygen" || name == "total_oxygen") && oxygen) {
			for (std::size_t i = 0; i < cells; ++i) {
				const auto equilibrium = OxygenFromTransported(oxygen->concentration[i],
					configuration.physiology.oxygen_state, oxygen_parameters);
				values[i] = name == "bound_oxygen"
					? equilibrium.bound_oxygen_mol_m3
					: equilibrium.total_oxygen_mol_m3;
			}
		} else if (name == "hematocrit") {
			std::fill(values.begin(), values.end(), configuration.physiology.hematocrit_percent);
		} else continue;
		result.emplace(name, std::move(values));
	}
	return result;
}

inline void ApplyOneDVasodilation(const OneDConfiguration& configuration,
	OneDNetwork& network, const std::vector<OneDTransportState>& transports,
	double dt, double dynamic_viscosity)
{
	if (!configuration.physiology.enabled || !configuration.physiology.vasodilation) return;
	const auto* signal = FindOneDSpecies(transports, configuration.physiology.vasodilator_field);
	if (!signal) throw std::runtime_error("vasodilation requires its configured transported field");
	const double relaxation = 1.0-std::exp(-dt/configuration.physiology.relaxation_tau);
	for (auto& segment : network.segments) {
		double mean = 0.0;
		for (int cell = 0; cell < segment.cells; ++cell)
			mean += signal->concentration[static_cast<std::size_t>(segment.cell_offset+cell)];
		mean /= segment.cells;
		const double response = mean/(configuration.physiology.ec50+std::max(mean, 0.0));
		const double target = segment.baseline_radius0
			*(1.0+configuration.physiology.emax_radius_fraction*response);
		segment.radius0 += relaxation*(target-segment.radius0);
		segment.area0 = OneDPi*segment.radius0*segment.radius0;
		segment.resistance = 8.0*dynamic_viscosity*segment.length
			/(OneDPi*std::pow(segment.radius0, 4.0));
	}
}

} // namespace iga

#endif

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
	double root_native_flux = 0.0;
	std::map<int, double> outlet_native_flux;
	bool boundary_flux_valid = false;
	double step_initial_mass = 0.0;
	double step_root_native_amount = 0.0;
	std::map<int, double> step_outlet_native_amount;
	double step_source_amount = 0.0;
	bool step_accounting_valid = false;
};

struct OneDTransportState {
	std::string name;
	std::vector<OneDSpeciesState> species;
};

struct OneDSpeciesStepAccounting {
	double initial_mass = 0.0;
	double final_mass = 0.0;
	double root_outward_amount = 0.0;
	std::map<int, double> outlet_outward_amount;
	double source_amount = 0.0;
	double balance_residual = 0.0;
};

inline double OneDSpeciesFaceFlux(double flow_m3_s, double left_concentration,
	double right_concentration, double face_area_m2, double diffusivity_m2_s,
	double spacing_m)
{
	if (!std::isfinite(flow_m3_s) || !std::isfinite(left_concentration)
		|| !std::isfinite(right_concentration) || !(face_area_m2 > 0.0)
		|| !std::isfinite(face_area_m2) || diffusivity_m2_s < 0.0
		|| !std::isfinite(diffusivity_m2_s) || !(spacing_m > 0.0)
		|| !std::isfinite(spacing_m))
		throw std::runtime_error("1d species face flux requires finite physical inputs");
	const double upwind = flow_m3_s >= 0.0
		? left_concentration : right_concentration;
	return flow_m3_s*upwind-face_area_m2*diffusivity_m2_s
		*(right_concentration-left_concentration)/spacing_m;
}

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

inline double OneDSpeciesMass(const OneDNetwork& network,
	const OneDFlowState& flow, const OneDSpeciesState& species)
{
	double mass = 0.0;
	for (const auto& segment : network.segments) {
		const double dx = segment.length/segment.cells;
		for (int cell = 0; cell < segment.cells; ++cell) {
			const auto index = static_cast<std::size_t>(segment.cell_offset+cell);
			mass += flow.area.at(index)*species.concentration.at(index)*dx;
		}
	}
	return mass;
}

inline void ResetOneDSpeciesStepAccounting(const OneDNetwork& network,
	const OneDFlowState& flow, OneDSpeciesState& species)
{
	species.root_native_flux = 0.0;
	species.outlet_native_flux.clear();
	species.boundary_flux_valid = false;
	species.step_initial_mass = OneDSpeciesMass(network, flow, species);
	species.step_root_native_amount = 0.0;
	species.step_outlet_native_amount.clear();
	species.step_source_amount = 0.0;
	species.step_accounting_valid = false;
}

inline OneDSpeciesStepAccounting GetOneDSpeciesStepAccounting(
	const OneDNetwork& network, const OneDFlowState& flow,
	const OneDSpeciesState& species)
{
	if (!species.step_accounting_valid)
		throw std::runtime_error("1d species step accounting requires a successful trial");
	OneDSpeciesStepAccounting result;
	result.initial_mass = species.step_initial_mass;
	result.final_mass = OneDSpeciesMass(network, flow, species);
	result.root_outward_amount = -species.step_root_native_amount;
	result.outlet_outward_amount = species.step_outlet_native_amount;
	result.source_amount = species.step_source_amount;
	double outward_amount = result.root_outward_amount;
	for (const auto& outlet : result.outlet_outward_amount)
		outward_amount += outlet.second;
	result.balance_residual = result.final_mass-result.initial_mass
		+outward_amount-result.source_amount;
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
	const std::filesystem::path& case_directory, double start_time, double requested_dt,
	const std::vector<double>* initial_area = nullptr,
	const std::map<std::string, double>* root_concentrations = nullptr,
	const std::map<int, std::map<std::string, double>>* outlet_concentrations = nullptr)
{
	double remaining = requested_dt;
	std::vector<double> scalar(species.concentration.size());
	std::vector<double> next(species.concentration.size());
	if (initial_area && initial_area->size() != scalar.size())
		throw std::runtime_error("1d transport initial area size mismatch");
	for (std::size_t i = 0; i < scalar.size(); ++i) {
		const double storage_area = initial_area ? initial_area->at(i) : flow.area.at(i);
		if (!(storage_area > 0.0) || !std::isfinite(storage_area))
			throw std::runtime_error("1d transport initial area must be finite and positive");
		scalar[i] = storage_area*species.concentration[i];
	}
	while (remaining > 0.0) {
		const double stable = OneDTransportStableDt(network, flow, species);
		const double dt = std::min(remaining,
			std::isfinite(stable) && stable > 0.0 ? stable : remaining);
		const double elapsed = requested_dt-remaining;
		double inlet_concentration = EvaluateOneDSpeciesInlet(configuration,
			species, case_directory, start_time+elapsed+dt);
		bool root_concentration_supplied = false;
		if (root_concentrations) {
			const auto supplied = root_concentrations->find(species.definition.field);
			if (supplied != root_concentrations->end()) {
				inlet_concentration = supplied->second;
				root_concentration_supplied = true;
			}
		}
		double root_native_flux = 0.0;
		std::map<int, double> outlet_native_flux;
		double source_amount = 0.0;
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
			if (segment.parent == network.root) {
				const double root_flow = flow.flow.at(
					static_cast<std::size_t>(segment.cell_offset));
				concentration[0] = root_concentrations && !root_concentration_supplied
					&& root_flow < -configuration.coupling.flow_epsilon_m3_s
					? concentration[1] : inlet_concentration;
			}
			else {
				const int incoming = OneDSegmentIntoNode(network, segment.parent);
				const auto& parent = network.segments[static_cast<std::size_t>(incoming)];
				concentration[0] = scalar[static_cast<std::size_t>(parent.cell_offset+parent.cells-1)]
					/flow.area[static_cast<std::size_t>(parent.cell_offset+parent.cells-1)];
			}
			const bool outlet = std::find(network.outlet_nodes.begin(),
				network.outlet_nodes.end(), segment.child) != network.outlet_nodes.end();
			bool outlet_concentration_supplied = false;
			if (outlet && outlet_concentrations) {
				const auto node = outlet_concentrations->find(segment.child);
				if (node != outlet_concentrations->end()) {
					const auto supplied = node->second.find(species.definition.field);
					if (supplied != node->second.end()) {
						concentration[static_cast<std::size_t>(segment.cells+1)]
							= supplied->second;
						outlet_concentration_supplied = true;
					}
				}
			}
			if (!outlet_concentration_supplied)
				concentration[static_cast<std::size_t>(segment.cells+1)]
					= concentration[static_cast<std::size_t>(segment.cells)];
			area[0] = area[1];
			area[static_cast<std::size_t>(segment.cells+1)] = area[static_cast<std::size_t>(segment.cells)];
			for (int face = 0; face <= segment.cells; ++face) {
				const double face_flow = face == 0
					? flow.flow[static_cast<std::size_t>(segment.cell_offset)]
					: flow.flow[static_cast<std::size_t>(segment.cell_offset+std::min(face, segment.cells)-1)];
				const double face_area = 0.5*(area[static_cast<std::size_t>(face)]+area[static_cast<std::size_t>(face+1)]);
				q[static_cast<std::size_t>(face)] = OneDSpeciesFaceFlux(face_flow,
					concentration[static_cast<std::size_t>(face)],
					concentration[static_cast<std::size_t>(face+1)], face_area,
					species.definition.diffusivity, dx);
			}
			if (segment.parent == network.root)
				root_native_flux += q.front();
			if (outlet)
				outlet_native_flux[segment.child] = q.back();
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
				const double source_rate = flow.area[index]*source-perimeter*wall_flux;
				next[index] = scalar[index]-dt/dx*(q[static_cast<std::size_t>(cell+1)]-q[static_cast<std::size_t>(cell)])
					+dt*source_rate;
				if (!std::isfinite(next[index])) throw std::runtime_error("1d species update produced a non-finite state");
				source_amount += dt*dx*source_rate;
			}
		}
		species.root_native_flux = root_native_flux;
		species.outlet_native_flux = std::move(outlet_native_flux);
		species.boundary_flux_valid = true;
		species.step_root_native_amount += dt*root_native_flux;
		for (const auto& outlet : species.outlet_native_flux)
			species.step_outlet_native_amount[outlet.first] += dt*outlet.second;
		species.step_source_amount += source_amount;
		scalar.swap(next);
		remaining -= dt;
	}
	for (std::size_t i = 0; i < scalar.size(); ++i)
		species.concentration[i] = scalar[i]/flow.area[i];
}

inline void AdvanceOneDTransport(const OneDConfiguration& configuration,
	const OneDNetwork& network, const OneDFlowState& flow, OneDTransportState& transport,
	const std::filesystem::path& case_directory, double start_time, double dt,
	const std::vector<double>* initial_area = nullptr,
	const std::map<std::string, double>* root_concentrations = nullptr,
	const std::map<int, std::map<std::string, double>>* outlet_concentrations = nullptr)
{
	if (outlet_concentrations)
		for (const auto& segment : network.segments) {
			if (std::find(network.outlet_nodes.begin(), network.outlet_nodes.end(),
				segment.child) == network.outlet_nodes.end()
				|| flow.flow.at(static_cast<std::size_t>(
					segment.cell_offset+segment.cells-1))
					>= -configuration.coupling.flow_epsilon_m3_s)
				continue;
			const auto node = outlet_concentrations->find(segment.child);
			if (node == outlet_concentrations->end()) continue;
			for (const auto& species : transport.species)
				if (node->second.find(species.definition.field) == node->second.end())
					throw std::runtime_error("reversed 1d outlet flow requires a supplied concentration for species '"
						+species.definition.field+"'");
		}
	#ifdef _OPENMP
	#pragma omp parallel for schedule(static) if(transport.species.size() >= 4)
	#endif
	for (long long i = 0; i < static_cast<long long>(transport.species.size()); ++i)
		AdvanceOneDSpecies(configuration, network, flow,
			transport.species[static_cast<std::size_t>(i)], case_directory, start_time, dt,
			initial_area, root_concentrations, outlet_concentrations);
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

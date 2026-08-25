#include "CaseInput.hpp"
#include "GenericCaseInput.hpp"
#include "IgaDatabase.hpp"
#include "OutletModel.hpp"
#include "TemporalFunction.hpp"

#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	try {
		if (argc != 3) throw std::runtime_error("usage: iga_case_check DATABASE.ntiga CASE_DIR");
		iga::Database database(argv[1]);
		const fs::path case_dir(argv[2]);
		const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
		const auto velocity = iga::ReadVelocity((case_dir/"initial_velocityfield.txt").string(), database.header().nodes);
		std::map<int, std::size_t> label_counts;
		for (int label : labels) ++label_counts[label];
		if (fs::exists(case_dir/"simulation_config.json")) {
			const auto configuration = iga::ReadSimulationConfiguration(
				(case_dir/"simulation_config.json").string());
			std::cout << "case_config=simulation_config.json schema_version="
				<< configuration.schema_version << " nodes=" << labels.size() << '\n';
			for (const auto& item : label_counts)
				std::cout << "label=" << item.first << " nodes=" << item.second << '\n';
			for (const auto& system : configuration.equation_systems) {
				if (system.kind == iga::EquationKind::LinearTransport) {
					const auto compiled = iga::CompileLinearSystem(configuration, system.name);
					const auto at_first_step = iga::MaterializeBoundaryWaveforms(
						configuration, case_dir.string(), configuration.time.dt);
					const auto boundaries = iga::ResolveScalarBoundaries(
						at_first_step, compiled, labels);
					std::cout << "system=" << system.name << " kind=linear_transport fields="
						<< compiled.fields.size() << " constrained_dofs="
						<< boundaries.constrained_dofs << '\n';
				} else {
					auto flow_configuration = configuration;
					if (system.time_integration == "backward_euler")
						flow_configuration = iga::MaterializeBoundaryWaveforms(
							configuration, case_dir.string(), 0.0);
					const auto models = iga::InitializeOutletModels(configuration, system);
					flow_configuration = iga::MaterializeOutletPressures(
						flow_configuration, models);
					const auto& flow = iga::FindEquationSystem(
						flow_configuration, system.name);
					const auto boundaries = iga::ResolveFlowBoundaries(
						flow_configuration, flow, labels, velocity);
					std::cout << "system=" << system.name << " kind=navier_stokes velocity_nodes="
						<< boundaries.velocity_nodes << " pressure_nodes="
						<< boundaries.pressure_nodes << " outlet_models=" << models.size() << '\n';
				}
			}
			for (const auto& boundary : configuration.boundaries)
				std::cout << "boundary label=" << boundary.label << " name=" << boundary.name
					<< " conditions=" << boundary.conditions.size() << '\n';
			return 0;
		}

		const auto parameters = iga::ReadTransportParameters((case_dir/"simulation_parameter.txt").string());
		const auto configuration = iga::ReadCaseConfiguration((case_dir/"case_config.json").string());
		const auto boundaries = iga::ResolveBoundaryConditions(configuration, labels, velocity, parameters);
		std::cout << "case_config=" << (configuration.present ? "case_config.json" : "legacy-defaults")
			<< " inherit_legacy=" << (configuration.inherit_legacy ? "true" : "false")
			<< " nodes=" << labels.size()
			<< " velocity_nodes=" << boundaries.velocity_nodes
			<< " pressure_nodes=" << boundaries.pressure_nodes
			<< " transport_nodes=" << boundaries.transport_nodes << '\n';
		for (const auto& item : label_counts)
			std::cout << "label=" << item.first << " nodes=" << item.second << '\n';
		for (const auto& rule : configuration.boundaries)
			std::cout << "rule label=" << rule.label << " type=" << iga::BoundaryTypeName(rule.type)
				<< (rule.name.empty() ? "" : " name="+rule.name) << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "iga_case_check: " << error.what() << '\n';
		return 1;
	}
}

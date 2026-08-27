#include "SimulationConfig.hpp"
#include "OneDConfig.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

int main(int argc, char** argv)
{
	try {
		if (argc != 2) throw std::runtime_error("usage: iga_config_check SIMULATION_CONFIG.json");
		std::ifstream input(argv[1]);
		if (!input) throw std::runtime_error("cannot open simulation configuration");
		std::ostringstream contents;
		contents << input.rdbuf();
		const auto text = contents.str();
		const auto root = iga::config_detail::RequireObject(
			iga::config_detail::JsonParser(text).Parse(), "root");
		const auto* version_value = iga::config_detail::Find(root, "schema_version");
		if (!version_value) throw std::runtime_error("simulation_config.json requires schema_version");
		const int version = iga::config_detail::RequireInteger(*version_value, "schema_version");
		if (version == 3) {
			const auto* dimension_value = iga::config_detail::Find(root, "dimension");
			if (!dimension_value) throw std::runtime_error("schema_version 3 requires dimension");
			const auto dimension = iga::config_detail::RequireString(*dimension_value, "dimension");
			if (dimension == "1d") {
				const auto configuration = iga::ParseOneDConfiguration(text);
				std::cout << "schema_version=3 dimension=" << configuration.dimension
					<< " fields=" << configuration.fields.size()
					<< " flow_systems=" << configuration.flow_systems.size()
					<< " transport_systems=" << configuration.transport_systems.size()
					<< " boundaries=" << configuration.boundaries.size() << '\n';
				for (const auto& system : configuration.flow_systems)
					std::cout << "system=" << system.name << " kind=network_flow_1d unknowns="
						<< system.unknowns.size() << '\n';
				for (const auto& system : configuration.transport_systems)
					std::cout << "system=" << system.name << " kind=network_transport_1d unknowns="
						<< system.unknowns.size() << " flow_system=" << system.flow_system << '\n';
				return 0;
			}
			if (dimension != "3d") throw std::runtime_error("dimension must be '1d' or '3d'");
		}
		const auto configuration = iga::ParseSimulationConfiguration(text);
		std::cout << "schema_version=" << configuration.schema_version
			<< " dimension=" << configuration.dimension
			<< " fields=" << configuration.fields.size()
			<< " systems=" << configuration.equation_systems.size()
			<< " boundaries=" << configuration.boundaries.size() << '\n';
		if (configuration.has_mesh)
			std::cout << "geometry=" << configuration.geometry.kind
				<< " file=" << configuration.geometry.file
				<< " target_spacing=" << configuration.mesh.centerline.target_spacing
				<< " min_scaled_J=" << configuration.mesh.quality.minimum_scaled_jacobian
				<< " self_intersection="
				<< (configuration.mesh.quality.check_self_intersection ? "on" : "off") << '\n';
		for (const auto& system : configuration.equation_systems) {
			std::cout << "system=" << system.name << " unknowns=" << system.unknowns.size();
			if (system.kind == iga::EquationKind::LinearTransport)
				std::cout << " kind=linear_transport terms=" << system.terms.size();
			else std::cout << " kind=navier_stokes viscosity=" << system.viscosity
				<< " density=" << system.density
				<< " time_integration=" << system.time_integration;
			std::cout << '\n';
		}
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

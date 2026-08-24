#include "SimulationConfig.hpp"

#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
	try {
		if (argc != 2) throw std::runtime_error("usage: iga_config_check SIMULATION_CONFIG.json");
		const auto configuration = iga::ReadSimulationConfiguration(argv[1]);
		std::cout << "schema_version=" << configuration.schema_version
			<< " fields=" << configuration.fields.size()
			<< " systems=" << configuration.equation_systems.size()
			<< " boundaries=" << configuration.boundaries.size() << '\n';
		for (const auto& system : configuration.equation_systems) {
			std::cout << "system=" << system.name << " unknowns=" << system.unknowns.size();
			if (system.kind == iga::EquationKind::LinearTransport)
				std::cout << " kind=linear_transport terms=" << system.terms.size();
			else std::cout << " kind=navier_stokes viscosity=" << system.viscosity;
			std::cout << '\n';
		}
		return 0;
	} catch (const std::exception& error) {
		std::cerr << error.what() << '\n';
		return 1;
	}
}

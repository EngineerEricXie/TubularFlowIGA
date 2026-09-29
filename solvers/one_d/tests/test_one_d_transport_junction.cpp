#include "OneDTransport.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
	for (const bool reversed : {false, true})
		for (const double diffusivity : {0.0, 0.1}) {
			iga::OneDNetwork network;
			network.root = 0;
			network.nodes.resize(4);
			network.outlet_nodes = {2, 3};
			network.cells = 3;
			network.segments = {
				{0, 0, 1, 1.0, 1.0, 1.0, 1.0, 0.0, {}, 0, 1},
				{1, 1, 2, 2.0, 1.0, 1.0, 2.0, 0.0, {}, 1, 1},
				{2, 1, 3, 4.0, 1.0, 1.0, 3.0, 0.0, {}, 2, 1}};
			iga::OneDFlowState flow;
			flow.area = {1.0, 2.0, 3.0};
			flow.flow = reversed ? std::vector<double>{2.0, -1.0, 3.0}
				: std::vector<double>{0.0, 0.0, 0.0};
			iga::OneDSpeciesState species;
			species.definition.field = "oxygen";
			species.definition.diffusivity = diffusivity;
			species.concentration = {1.0, 4.0, 2.0};
			species.inlet_value = 1.0;
			iga::OneDConfiguration configuration;
			iga::ResetOneDSpeciesStepAccounting(network, flow, species);
			iga::AdvanceOneDSpecies(configuration, network, flow, species, ".", 0.0, 0.01);
			species.step_accounting_valid = true;
			const auto accounting = iga::GetOneDSpeciesStepAccounting(network, flow, species);
			assert(std::abs(accounting.balance_residual) < 1.0e-13);
			for (const auto value : species.concentration)
				assert(std::isfinite(value) && value >= 1.0 && value <= 4.0);
		}
	std::cout << "junction mixing and diffusion conservation passed\n";
}

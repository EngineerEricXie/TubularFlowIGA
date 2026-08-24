#include "CaseInput.hpp"
#include "IgaDatabase.hpp"

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
		const auto parameters = iga::ReadTransportParameters((case_dir/"simulation_parameter.txt").string());
		const auto labels = iga::ReadPointLabels((case_dir/"controlmesh.vtk").string(), database.header().nodes);
		const auto velocity = iga::ReadVelocity((case_dir/"initial_velocityfield.txt").string(), database.header().nodes);
		const auto configuration = iga::ReadCaseConfiguration((case_dir/"case_config.json").string());
		const auto boundaries = iga::ResolveBoundaryConditions(configuration, labels, velocity, parameters);
		std::map<int, std::size_t> label_counts;
		for (int label : labels) ++label_counts[label];
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

#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

iga::NativeTetMesh Read(const char* path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error(std::string("cannot open mesh: ")+path);
	return iga::ReadNativeTetMeshGmsh41(input);
}

int Label(const char* value)
{
	std::size_t used = 0;
	const int label = std::stoi(value, &used);
	if (used != std::string(value).size() || label < 0)
		throw std::invalid_argument("invalid boundary label");
	return label;
}

double Mobility(const char* value)
{
	std::size_t used = 0;
	const double mobility = std::stod(value, &used);
	if (used != std::string(value).size() || !(mobility > 0.)
			|| !std::isfinite(mobility))
		throw std::invalid_argument("invalid Darcy mobility in m2/(Pa s)");
	return mobility;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc != 7) throw std::runtime_error(
			"usage: native_tet_matching_darcy_smoke vessel.msh tissue.msh interface_label exit_label_1 exit_label_2 mobility_m2_pa_s");
		const auto vessel = Read(argv[1]);
		const auto tissue = Read(argv[2]);
		const int interface = Label(argv[3]);
		const int exit1 = Label(argv[4]);
		const int exit2 = Label(argv[5]);
		if (interface == exit1 || interface == exit2 || exit1 == exit2)
			throw std::invalid_argument("interface and artificial exit labels must differ");
		const double mobility = Mobility(argv[6]);
		const auto topology = iga::BuildNativeTaylorHoodTopology(vessel);
		const std::size_t velocity_nodes = vessel.points.size()+topology.edges.size();
		std::vector<double> state(3*velocity_nodes+vessel.points.size(), 0.0);
		double minimum_x = vessel.points.front()[0];
		for (const auto& point : vessel.points) minimum_x = std::min(minimum_x, point[0]);
		for (std::size_t node = 0; node < vessel.points.size(); ++node)
			state[3*node] = (vessel.points[node][0]-minimum_x)/0.01;
		for (std::size_t edge = 0; edge < topology.edges.size(); ++edge) {
			const auto& ends = topology.edges[edge];
			const double midpoint_x = 0.5*(vessel.points[ends[0]][0]
				+vessel.points[ends[1]][0]);
			state[3*(vessel.points.size()+edge)] = (midpoint_x-minimum_x)/0.01;
		}
		const iga::NativeTetMatchingInterfacePort port{"manufactured_geometry_probe",
			interface, interface};
		const auto mapped = iga::MapNativeTetMatchingInterfaceToTissueSource(
			vessel, topology, state, tissue, {port});
		const std::vector<double> mobility_cells(tissue.cells.size(), mobility);
		const auto darcy = iga::SolveNativeTetDarcyPetsc(tissue, mobility_cells,
			mapped.source.tissue_source_s_inv, {{exit1, 0.0}, {exit2, 0.0}});
		const double source = mapped.source.vessel_outward_flow_m3_s.at(port.name);
		const double outward = darcy.conservative_outward_boundary_flow_m3_s.at(exit1)
			+darcy.conservative_outward_boundary_flow_m3_s.at(exit2);
		const double tolerance = 1e-8*std::max(1e-12, std::abs(source));
		if (darcy.converged_reason <= 0 || darcy.maximum_cell_balance_defect_m3_s > tolerance
				|| std::abs(source-darcy.volume_source_m3_s) > tolerance
				|| std::abs(source-outward) > tolerance)
			throw std::runtime_error("manufactured interface-Darcy conservation or convergence gate failed");
		if (rank == 0) std::cout << std::setprecision(17)
			<< "native_tet_matching_darcy_smoke: PASS facets="
			<< mapped.matched_facets.at(port.name)
			<< " source_m3_s=" << source
			<< " darcy_outward_m3_s=" << outward
			<< " max_cell_defect_m3_s=" << darcy.maximum_cell_balance_defect_m3_s
			<< " converged_reason=" << darcy.converged_reason << '\n';
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << "native_tet_matching_darcy_smoke: ERROR: "
			<< error.what() << '\n';
		status = 2;
	}
	PetscFinalize();
	return status;
}

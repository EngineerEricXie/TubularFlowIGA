#include "NativeTetMatchingInterfaceSource.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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

} // namespace

int main(int argc, char** argv)
{
	try {
		if (argc != 4)
			throw std::runtime_error("usage: native_tet_matching_mesh_inspect vessel.msh tissue.msh interface_label");
		const auto vessel = Read(argv[1]);
		const auto tissue = Read(argv[2]);
		const int label = std::stoi(argv[3]);
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
		const iga::NativeTetMatchingInterfacePort port{"geometry_probe", label, label};
		const auto result = iga::MapNativeTetMatchingInterfaceToTissueSource(
			vessel, topology, state, tissue, {port});
		const double flow = result.source.vessel_outward_flow_m3_s.at(port.name);
		double deposited = 0.0;
		for (std::size_t cell = 0; cell < tissue.cells.size(); ++cell)
			deposited += result.source.tissue_source_s_inv[cell]
				*iga::EvaluateNativeTetGeometry(tissue, tissue.cells[cell]).determinant/6.0;
		std::cout << std::setprecision(17);
		if (!std::isfinite(flow) || !std::isfinite(deposited)
				|| std::abs(flow-deposited) > 1e-12*std::max(1e-12, std::abs(flow))) {
			std::ostringstream message;
			message << std::setprecision(17)
				<< "manufactured matching source is not conservative: flow=" << flow
				<< " deposited=" << deposited;
			throw std::runtime_error(message.str());
		}
		std::cout << "native_tet_matching_mesh_inspect: PASS facets="
			<< result.matched_facets.at(port.name)
			<< " manufactured_velocity_x_m_s=(x-min_x)/0.01_m"
			<< " vessel_outward_flow_m3_s=" << flow
			<< " tissue_source_flow_m3_s=" << deposited << '\n';
		return 0;
	} catch (const std::exception& error) {
		std::cerr << "native_tet_matching_mesh_inspect: ERROR: " << error.what() << '\n';
		return 2;
	}
}

#include "NativeTetAlePetscRuntime.hpp"
#include "NativeTetBoundaryFlow.hpp"
#include "NativeTetDarcyPetsc.hpp"
#include "NativeTetDarcyVisualization.hpp"
#include "NativeTetHydraulicVisualization.hpp"
#include "NativeTetMatchingInterfaceSource.hpp"
#include "NativeTetMeshComponents.hpp"
#include "ParallelVtkOutput.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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

int Label(const char* token)
{
	std::size_t used = 0;
	const int value = std::stoi(token, &used);
	if (used != std::string(token).size() || value < 0)
		throw std::invalid_argument("boundary label must be a nonnegative integer");
	return value;
}

double Real(const char* token, const char* name)
{
	std::size_t used = 0;
	const double value = std::stod(token, &used);
	if (used != std::string(token).size() || !std::isfinite(value))
		throw std::invalid_argument(std::string(name)+" must be a finite number");
	return value;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	int status = 0;
	try {
		if (argc != 15 && argc != 16) throw std::invalid_argument(
			"usage: native_tet_matching_solved_roi_smoke vessel.msh tissue.msh interface_label inlet_label tissue_exit_label_1 tissue_exit_label_2 inlet_vx_m_s inlet_vy_m_s inlet_vz_m_s density_kg_m3 viscosity_pa_s mobility_m2_pa_s interface_pressure_pa tissue_exit_pressure_pa [field_output_dir]");
		const auto vessel = Read(argv[1]);
		const auto tissue = Read(argv[2]);
		iga::RequireNativeTetSingleFaceComponent(vessel, "fluid");
		iga::RequireNativeTetSingleFaceComponent(tissue, "tissue");
		const int interface = Label(argv[3]);
		const int inlet = Label(argv[4]);
		const int exit1 = Label(argv[5]);
		const int exit2 = Label(argv[6]);
		if (interface == inlet || interface == exit1 || interface == exit2
				|| exit1 == exit2 || inlet == exit1 || inlet == exit2)
			throw std::invalid_argument("inlet, interface, and tissue exits must have distinct labels");
		const std::array<double, 3> inlet_velocity{{
			Real(argv[7], "inlet_vx_m_s"), Real(argv[8], "inlet_vy_m_s"),
			Real(argv[9], "inlet_vz_m_s")}};
		const double density = Real(argv[10], "density_kg_m3");
		const double viscosity = Real(argv[11], "viscosity_pa_s");
		const double mobility = Real(argv[12], "mobility_m2_pa_s");
		const double interface_pressure = Real(argv[13], "interface_pressure_pa");
		const double tissue_pressure = Real(argv[14], "tissue_exit_pressure_pa");
		if (!(density > 0. && viscosity > 0. && mobility > 0.)
				|| std::all_of(inlet_velocity.begin(), inlet_velocity.end(),
					[](double value){return value == 0.;}))
			throw std::invalid_argument("functional fluid/material inputs must be positive and inlet velocity nonzero");
		const auto topology = iga::BuildNativeTaylorHoodTopology(vessel);
		if (!topology.boundary_velocity_nodes.count(inlet)
				|| !topology.boundary_velocity_nodes.count(interface))
			throw std::invalid_argument("fluid inlet or interface label is absent");
		std::map<std::uint32_t, std::array<double, 3>> velocity_bc;
		for (const auto node : topology.boundary_velocity_nodes.at(inlet))
			velocity_bc.emplace(node, inlet_velocity);
		const std::size_t velocity_nodes = vessel.points.size()+topology.edges.size();
		std::vector<double> initial(3*velocity_nodes+vessel.points.size(), 0.);
		std::vector<std::array<double, 3>> grid_velocity(vessel.points.size(), {{0.,0.,0.}});
		iga::NativeTetAleBoundaryConditions fluid_bc;
		fluid_bc.prescribed_pressure_pa.emplace(interface, interface_pressure);
		const auto fluid = iga::SolveNativeTetAlePetscSteady(vessel, grid_velocity,
			initial, velocity_bc, std::numeric_limits<std::uint32_t>::max(),
			{density, viscosity}, 1e-10, 12, fluid_bc);
		if (fluid.final_linear_converged_reason <= 0 || fluid.final_residual_l2 > 1e-10)
			throw std::runtime_error("native fluid convergence gate failed");
		const auto flows = iga::EvaluateNativeTetBoundaryFlows(vessel, topology,
			fluid.replicated_state);
		const double inlet_flow = flows.at(inlet).outward_flow_m3_s;
		const double interface_flow = flows.at(interface).outward_flow_m3_s;
		const double fluid_scale = std::max({std::abs(inlet_flow),
			std::abs(interface_flow), 1e-18});
		if (!(inlet_flow < 0. && interface_flow > 0.)
				|| std::abs(inlet_flow+interface_flow) > 1e-6*fluid_scale)
			throw std::runtime_error("native fluid boundary flow or global mass gate failed");
		const iga::NativeTetMatchingInterfacePort port{"solved_roi_functional",
			interface, interface};
		const auto mapped = iga::MapNativeTetMatchingInterfaceToTissueSource(
			vessel, topology, fluid.replicated_state, tissue, {port});
		const double source = mapped.source.vessel_outward_flow_m3_s.at(port.name);
		double facet_source = 0.;
		for (const auto& facet : mapped.facet_transfers)
			facet_source += facet.vessel_outward_flow_m3_s;
		if (std::abs(source-interface_flow) > 1e-8*fluid_scale)
			throw std::runtime_error("native matching source differs from fluid boundary quadrature");
		if (mapped.facet_transfers.size() != mapped.matched_facets.at(port.name)
				|| std::abs(facet_source-source) > 1e-12*fluid_scale)
			throw std::runtime_error("native matching facet ledger differs from source");
		const std::vector<double> mobility_cells(tissue.cells.size(), mobility);
		const auto darcy = iga::SolveNativeTetDarcyPetsc(tissue, mobility_cells,
			mapped.source.tissue_source_s_inv,
			{{exit1, tissue_pressure}, {exit2, tissue_pressure}});
		const double outward = darcy.conservative_outward_boundary_flow_m3_s.at(exit1)
			+darcy.conservative_outward_boundary_flow_m3_s.at(exit2);
		const double exchange_scale = std::max(std::abs(source), 1e-18);
		if (darcy.converged_reason <= 0
				|| darcy.maximum_cell_balance_defect_m3_s > 1e-8*exchange_scale
				|| std::abs(source-darcy.volume_source_m3_s) > 1e-8*exchange_scale
				|| std::abs(source-outward) > 1e-8*exchange_scale)
			throw std::runtime_error("native Darcy convergence or source/outflow gate failed");
		if (argc == 16) {
			const std::filesystem::path output(argv[15]);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"ROI field output directory",[&] {
				if (output.empty()) throw std::invalid_argument("ROI field output directory is empty");
				if (rank != 0) return;
				if (!output.parent_path().empty())
					std::filesystem::create_directories(output.parent_path());
				if (!std::filesystem::create_directory(output))
					throw std::runtime_error("ROI field output directory already exists");
			});
			int ranks = 0;
			MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
			iga::VtkPartition fluid_piece, tissue_piece;
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"ROI field partition build",[&] {
				fluid_piece = iga::BuildNativeTetHydraulicVtkPartition(vessel,vessel,
					fluid.replicated_state,viscosity,rank,ranks);
				tissue_piece = iga::BuildNativeTetDarcyVtkPartition(tissue,darcy,
					mobility_cells,mapped.source.tissue_source_s_inv,rank,ranks);
			});
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"fluid",
				fluid_piece,0.);
			iga::WriteParallelVtkSnapshot(PETSC_COMM_WORLD,output/"tissue",
				tissue_piece,0.);
			iga::CollectiveLocalStage(PETSC_COMM_WORLD,"ROI interface ledger",[&] {
				if (rank != 0) return;
				std::ofstream ledger(output/"interface_facets.json");
				if (!ledger) throw std::runtime_error("cannot create ROI interface ledger");
				ledger << std::setprecision(17)
					<< "{\n  \"schema_version\": 1,\n"
					<< "  \"kind\": \"native_matching_facet_flow_functional_only\",\n"
					<< "  \"flow_sign\": \"vessel_outward_positive_tissue_source_positive\",\n"
					<< "  \"facets\": [\n";
				for (std::size_t index = 0; index < mapped.facet_transfers.size(); ++index) {
					const auto& facet = mapped.facet_transfers[index];
					ledger << "    {\"port_name\": \"solved_roi_functional\", "
						<< "\"vessel_cell_id\": " << facet.vessel_cell_id
						<< ", \"tissue_cell_id\": " << facet.tissue_cell_id
						<< ", \"vessel_outward_m3_s\": "
						<< facet.vessel_outward_flow_m3_s << ", \"triangle_m\": [";
					for (std::size_t node = 0; node < 3; ++node) {
						if (node) ledger << ", ";
						ledger << '[' << facet.triangle_m[node][0] << ", "
							<< facet.triangle_m[node][1] << ", "
							<< facet.triangle_m[node][2] << ']';
					}
					ledger << "]}" << (index+1 == mapped.facet_transfers.size() ? "\n" : ",\n");
				}
				ledger << "  ]\n}\n";
				ledger.close();
				if (!ledger) throw std::runtime_error("cannot write ROI interface ledger");
			});
		}
		if (rank == 0) std::cout << std::setprecision(17)
			<< "native_tet_matching_solved_roi_smoke: PASS facets="
			<< mapped.matched_facets.at(port.name)
			<< " fluid_inlet_m3_s=" << inlet_flow
			<< " interface_source_m3_s=" << source
			<< " darcy_outward_m3_s=" << outward
			<< " fluid_newton_iterations=" << fluid.newton_iterations
			<< " fluid_linear_reason=" << fluid.final_linear_converged_reason
			<< " darcy_linear_reason=" << darcy.converged_reason
			<< " darcy_max_cell_defect_m3_s=" << darcy.maximum_cell_balance_defect_m3_s
			<< '\n';
	} catch (const std::exception& error) {
		if (rank == 0) std::cerr << "native_tet_matching_solved_roi_smoke: ERROR: "
			<< error.what() << '\n';
		status = 2;
	}
	PetscFinalize();
	return status;
}

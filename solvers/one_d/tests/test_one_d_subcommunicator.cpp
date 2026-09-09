#include "OneDImplicit.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

double Compare(const std::vector<double>& reference, const std::vector<double>& current)
{
	if (reference.size() != current.size()) throw std::runtime_error("field sizes differ");
	double norm = 0.0, error = 0.0;
	for (std::size_t i = 0; i < reference.size(); ++i) {
		if (!std::isfinite(reference[i]) || !std::isfinite(current[i]))
			throw std::runtime_error("nonfinite field");
		norm = std::hypot(norm, reference[i]);
		error = std::hypot(error, current[i]-reference[i]);
	}
	if (norm == 0.0 ? error > 1e-12 : error/norm > 1e-6)
		throw std::runtime_error("subcommunicator field differs from independent serial reference");
	return norm > 0.0 ? error/norm : error;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		if (argc != 2 || ranks < 2)
			throw std::runtime_error("usage: mpiexec -np N>=2 one_d_subcommunicator_test INPUT.swc");
		const int color = rank == 0 ? 0 : 1;
		MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
		int local_rank = 0, local_size = 1;
		MPI_Comm_rank(group, &local_rank); MPI_Comm_size(group, &local_size);
		iga::OneDFlowSystemDefinition flow;
		flow.name = "subcommunicator-flow";
		flow.model = iga::OneDFlowModel::Compliant;
		flow.scheme = iga::OneDFlowScheme::ImplicitPetsc;
		flow.dynamic_viscosity = 0.004*(color+1);
		flow.density = 1060.0;
		flow.discretization.cells_per_segment = 3;
		const auto network = iga::ReadOneDNetwork(argv[1], 1.0, 3, flow.dynamic_viscosity);
		const double inlet = (color+1)*1e-9, dt = 1e-3;
		for (auto formulation : {iga::OneDImplicitFormulation::PressureNetwork,
			iga::OneDImplicitFormulation::LinearizedAQ,
			iga::OneDImplicitFormulation::NonlinearAQ,
			iga::OneDImplicitFormulation::ImplicitPde}) {
			flow.formulation = formulation;
			iga::OneDFlowState initial;
			for (auto node : network.outlet_nodes) {
				iga::OneDOutletState outlet;
				outlet.node = node; outlet.kind = iga::OneDOutletKind::Pressure;
				initial.outlets.push_back(outlet);
			}
			iga::InitializeCompliantOneDFromRigid(network, flow, initial, inlet, dt);
			auto reference = initial, current = initial;
			// Each rank computes its own COMM_SELF oracle. Disjoint groups use
			// different data AND different solve counts, exposing world collectives.
			for (int step = 0; step < color+1; ++step) {
				iga::AdvanceImplicitOneD(network, flow, reference, inlet, dt, PETSC_COMM_SELF);
				iga::AdvanceImplicitOneD(network, flow, current, inlet, dt, group);
			}
			const double area = Compare(reference.area, current.area);
			const double pressure = Compare(reference.pressure, current.pressure);
			const double node_pressure = Compare(reference.node_pressure, current.node_pressure);
			const double flux = Compare(reference.flow, current.flow);
			const double segment_flux = Compare(reference.segment_flow, current.segment_flow);
			if (local_rank == 0)
				std::cout << "group=" << color << " ranks=" << local_size
					<< " formulation=" << static_cast<int>(formulation)
					<< " area_error=" << area << " pressure_error=" << pressure
					<< " node_pressure_error=" << node_pressure << " flow_error=" << flux
					<< " segment_flow_error=" << segment_flux << '\n';
		}
		MPI_Comm_free(&group);
		if (rank == 0) std::cout << "one_d_subcommunicator_test passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		// This is a bounded test harness, not a runtime recovery protocol.
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

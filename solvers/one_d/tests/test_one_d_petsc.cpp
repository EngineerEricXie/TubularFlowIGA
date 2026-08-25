#include "OneDCheckpoint.hpp"

#include <petscsys.h>

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int rank = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	const fs::path swc = fs::temp_directory_path()/"tubularflowiga-one-d-petsc-test.swc";
	if (rank == 0) {
		std::ofstream output(swc);
		output << "1 2 0 0 0 0.001 -1\n2 2 0.01 0 0 0.001 1\n3 2 0.02 0 0 0.001 2\n";
	}
	MPI_Barrier(PETSC_COMM_WORLD);
	iga::OneDFlowSystemDefinition flow;
	flow.name = "flow";
	flow.model = iga::OneDFlowModel::Compliant;
	flow.scheme = iga::OneDFlowScheme::ImplicitPetsc;
	flow.dynamic_viscosity = 0.004;
	flow.density = 1060.0;
	flow.discretization.cells_per_segment = 3;
	auto network = iga::ReadOneDNetwork(swc, 1.0, 3, flow.dynamic_viscosity);
	iga::OneDFlowState state;
	for (const int node : network.outlet_nodes) {
		iga::OneDOutletState outlet;
		outlet.node = node;
		outlet.kind = iga::OneDOutletKind::Pressure;
		state.outlets.push_back(outlet);
	}
	iga::InitializeCompliantOneDFromRigid(network, flow, state, 1.0e-9, 1.0e-3);
	for (const auto formulation : {iga::OneDImplicitFormulation::PressureNetwork,
		iga::OneDImplicitFormulation::LinearizedAQ,
		iga::OneDImplicitFormulation::NonlinearAQ,
		iga::OneDImplicitFormulation::ImplicitPde}) {
		flow.formulation = formulation;
		iga::AdvanceImplicitOneD(network, flow, state, 1.0e-9, 1.0e-3);
		assert(state.area.size() == static_cast<std::size_t>(network.cells));
		for (std::size_t i = 0; i < state.area.size(); ++i) {
			assert(state.area[i] > 0.0);
			assert(std::isfinite(state.flow[i]));
			assert(std::isfinite(state.pressure[i]));
		}
	}
	iga::OneDCheckpointMetadata metadata;
	metadata.completed_step = 1;
	metadata.internal_substeps = 3;
	metadata.physical_time = 0.001;
	metadata.dt = 0.001;
	metadata.inlet_flow = 1.0e-9;
	metadata.cells = network.cells;
	metadata.nodes = static_cast<int>(network.nodes.size());
	metadata.segments = static_cast<int>(network.segments.size());
	metadata.outlets = static_cast<int>(state.outlets.size());
	metadata.config_fingerprint = 18446744073709551557ull;
	metadata.network_fingerprint = iga::OneDNetworkFingerprint(network);
	metadata.state_file = "state.bin";
	const auto packed = iga::PackOneDCheckpointState(state, {}, network);
	auto restored_network = network;
	restored_network.segments.front().radius0 *= 0.5;
	iga::OneDFlowState restored;
	restored.outlets = state.outlets;
	std::vector<iga::OneDTransportState> no_transports;
	iga::UnpackOneDCheckpointState(packed, restored, no_transports, metadata,
		restored_network, flow.dynamic_viscosity);
	assert(restored.area == state.area);
	assert(restored.flow == state.flow);
	assert(restored_network.segments.front().radius0 == network.segments.front().radius0);
	bool rejected = false;
	auto truncated = packed;
	truncated.pop_back();
	try { iga::UnpackOneDCheckpointState(truncated, restored, no_transports, metadata,
		restored_network, flow.dynamic_viscosity); }
	catch (const std::runtime_error&) { rejected = true; }
	assert(rejected);
	const auto parsed = iga::ParseOneDCheckpointMetadata(
		iga::SerializeOneDCheckpointMetadata(metadata));
	assert(parsed.config_fingerprint == metadata.config_fingerprint);
	assert(parsed.network_fingerprint == metadata.network_fingerprint);
	assert(parsed.internal_substeps == metadata.internal_substeps);
	assert(parsed.inlet_flow == metadata.inlet_flow);
	MPI_Barrier(PETSC_COMM_WORLD);
	if (rank == 0) fs::remove(swc);
	if (rank == 0) std::cout << "one-dimensional PETSc tests passed\n";
	PetscFinalize();
}

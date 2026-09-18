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
	iga::OneDFlowState rigid_reference;
	rigid_reference.outlets = state.outlets;
	iga::SolveRigidOneD(network, flow, rigid_reference, 1.0e-9, 0.0);
	iga::OneDFlowSystemDefinition steady = flow;
	steady.model = iga::OneDFlowModel::Lumped;
	steady.scheme = iga::OneDFlowScheme::Petsc;
	steady.formulation = iga::OneDImplicitFormulation::SteadyR;
	iga::OneDFlowState steady_state;
	steady_state.outlets = state.outlets;
	iga::InitializeLumpedOneD(network, steady, steady_state, 1.0e-9);
	iga::AdvanceImplicitOneD(network, steady, steady_state, 1.0e-9, 1.0e-3);
	for (std::size_t i = 0; i < steady_state.node_pressure.size(); ++i)
		assert(std::abs(steady_state.node_pressure[i]-rigid_reference.node_pressure[i]) < 1.0e-10);
	for (std::size_t i = 0; i < steady_state.segment_flow.size(); ++i)
		assert(std::abs(steady_state.segment_flow[i]-rigid_reference.segment_flow[i]) < 1.0e-18);
	for (const auto formulation : {iga::OneDImplicitFormulation::TransientRc,
		iga::OneDImplicitFormulation::TransientRlc,
		iga::OneDImplicitFormulation::NonlinearRlc,
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
	iga::OneDFlowSystemDefinition lumped = flow;
	lumped.model = iga::OneDFlowModel::Lumped;
	lumped.scheme = iga::OneDFlowScheme::Petsc;
	lumped.formulation = iga::OneDImplicitFormulation::TransientRc;
	lumped.lumped.compliance_scale = 0.0;
	lumped.lumped.segment_resistance[2] = 1.0e8;
	lumped.lumped.segment_resistance[3] = 2.0e8;
	iga::OneDFlowState rc_state;
	iga::OneDOutletState rc_outlet;
	rc_outlet.node = network.outlet_nodes.front();
	rc_outlet.kind = iga::OneDOutletKind::WindkesselRc;
	rc_outlet.distal_resistance = 1.0e9;
	rc_outlet.capacitance = 1.0e-10;
	rc_state.outlets.push_back(rc_outlet);
	iga::InitializeLumpedOneD(network, lumped, rc_state, 1.0e-9);
	iga::AdvanceImplicitOneD(network, lumped, rc_state, 1.0e-9, 1.0e-2);
	const double expected_capacitor_pressure = 0.1/1.1;
	assert(std::abs(rc_state.outlets.front().capacitor_pressure
		-expected_capacitor_pressure) < 1.0e-11);
	assert(std::abs(rc_state.node_pressure.back()-expected_capacitor_pressure) < 1.0e-11);
	assert(std::abs(rc_state.node_pressure.front()
		-(expected_capacitor_pressure+0.3)) < 1.0e-11);

	iga::OneDFlowState rcr_state;
	iga::OneDOutletState rcr_outlet = rc_outlet;
	rcr_outlet.kind = iga::OneDOutletKind::WindkesselRcr;
	rcr_outlet.proximal_resistance = 2.0e8;
	rcr_state.outlets.push_back(rcr_outlet);
	iga::InitializeLumpedOneD(network, lumped, rcr_state, 1.0e-9);
	iga::AdvanceImplicitOneD(network, lumped, rcr_state, 1.0e-9, 1.0e-2);
	assert(std::abs(rcr_state.outlets.front().capacitor_pressure
		-expected_capacitor_pressure) < 1.0e-11);
	assert(std::abs(rcr_state.node_pressure.back()
		-(expected_capacitor_pressure+0.2)) < 1.0e-11);
	assert(std::abs(rcr_state.node_pressure.front()
		-(expected_capacitor_pressure+0.5)) < 1.0e-11);
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
	// Independent groups perform different numbers of checkpoint operations.
	// A world collective inside either API would hang this test.
	MPI_Comm group = MPI_COMM_NULL;
	const int color = rank == 0 ? 0 : 1;
	MPI_Comm_split(PETSC_COMM_WORLD, color, rank, &group);
	int local_rank = 0, local_ranks = 1;
	MPI_Comm_rank(group, &local_rank);
	MPI_Comm_size(group, &local_ranks);
	const auto prefix = fs::temp_directory_path()
		/("tubularflowiga-one-d-checkpoint-group-"+std::to_string(color));
	auto group_metadata = metadata;
	group_metadata.state_file = iga::OneDCheckpointStatePath(prefix).filename().string();
	for (int round = 0; round <= color; ++round) {
		auto group_state = state;
		for (auto& pressure : group_state.pressure) pressure += color+round;
		iga::WriteOneDCheckpoint(prefix, group_metadata, group_state, no_transports,
			network, local_rank, group);
		auto loaded = state;
		auto loaded_network = network;
		const auto loaded_metadata = iga::ReadOneDCheckpoint(prefix, loaded, no_transports,
			loaded_network, flow.dynamic_viscosity, group);
		assert(loaded_metadata.completed_step == group_metadata.completed_step);
		assert(iga::PackOneDCheckpointState(loaded, no_transports, loaded_network)
			== iga::PackOneDCheckpointState(group_state, no_transports, network));
	}
	if (local_ranks > 1) {
		int caught = 0;
		try {
			auto loaded = state;
			auto loaded_network = network;
			// The deliberately absent file lives below the fixture SWC file,
			// ensuring that no previous test can accidentally create it.
			iga::ReadOneDCheckpoint(local_rank == 1 ? swc/"absent" : prefix,
				loaded, no_transports, loaded_network, flow.dynamic_viscosity, group);
		} catch (const std::runtime_error& error) {
			caught = std::string(error.what()).find("1d checkpoint metadata read: rank 1:") != std::string::npos;
		}
		int all_caught = 0;
		MPI_Allreduce(&caught, &all_caught, 1, MPI_INT, MPI_MIN, group);
		assert(all_caught == 1);
	}
	MPI_Barrier(group);
	if (local_rank == 0) {
		fs::remove(iga::OneDCheckpointStatePath(prefix));
		fs::remove(iga::OneDCheckpointMetadataPath(prefix));
	}
	MPI_Comm_free(&group);
	MPI_Barrier(PETSC_COMM_WORLD);
	if (rank == 0) fs::remove(swc);
	if (rank == 0) std::cout << "one-dimensional PETSc tests passed\n";
	PetscFinalize();
}

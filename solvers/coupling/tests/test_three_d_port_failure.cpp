#include "BoundarySupport.hpp"
#include "ThreeDBodyFittedFlowDomainAdapter.hpp"
#include "../../cpu/tests/StringStreamFailure.hpp"

// Reuse the established 64-control-point bifurcation geometry. The renamed
// smoke entry is never called; this test runs runtimes directly on its group.
#define main BifurcationFixtureMain
#include "test_bifurcation_coupling_smoke.cpp"
#undef main

#include <limits>

namespace {

using Measurements = std::map<std::string, iga::PortState>;

std::vector<iga::CouplingPort> Ports()
{
	std::vector<iga::CouplingPort> ports;
	for (int label = 1; label <= 3; ++label) {
		iga::CouplingPort port;
		port.id = "port"+std::to_string(label);
		port.subsystem_id = "flow";
		port.locator_kind = "boundary_label";
		port.locator = std::to_string(label);
		port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
			iga::PortQuantity::MeanPressure};
		ports.push_back(port);
	}
	return ports;
}

void CheckValue(double reference, double actual)
{
	if (!std::isfinite(reference) || !std::isfinite(actual)
		|| (reference == 0.0 ? std::abs(actual) > 1e-12
			: std::abs((actual-reference)/reference) > 1e-6))
		throw std::runtime_error("port quantity differs from reference");
}

void Compare(const Measurements& reference, const Measurements& actual)
{
	if (reference.size() != actual.size()) throw std::runtime_error("port count differs");
	for (const auto& entry : reference) {
		const auto& a = entry.second;
		const auto& b = actual.at(entry.first);
		CheckValue(a.time_s, b.time_s);
		CheckValue(a.area_m2.value(), b.area_m2.value());
		CheckValue(a.outward_flow_m3_s.value(), b.outward_flow_m3_s.value());
		CheckValue(a.mean_pressure_pa.value(), b.mean_pressure_pa.value());
		for (const auto& pair : {std::make_pair(&a.concentration, &b.concentration),
			std::make_pair(&a.outward_species_flux, &b.outward_species_flux)}) {
			if (pair.first->size() != pair.second->size())
				throw std::runtime_error("species count differs");
			for (const auto& field : *pair.first) CheckValue(field.second, pair.second->at(field.first));
		}
	}
}

template <class Work>
void ExpectFailure(MPI_Comm communicator, const std::string& stage, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(communicator, "test expected failure", [&] {
		if (message.find(stage) == std::string::npos)
			throw std::runtime_error("expected "+stage+", got: "+message);
	});
	iga::RequireCollectiveSameText(communicator, "test common diagnostic", message);
}

Measurements Run(const fs::path& root, const fs::path& database_path,
	MPI_Comm communicator, bool faults)
{
	int rank = 0, ranks = 1;
	MPI_Comm_rank(communicator, &rank); MPI_Comm_size(communicator, &ranks);
	iga::Database database(database_path.string());
	const auto mesh = iga::ReadLabeledHexMesh((root/"controlmesh.vtk").string(), 64, 1);
	auto velocity = iga::ReadVelocity((root/"initial_velocityfield.txt").string(), 64);
	const auto configuration = iga::ReadSimulationConfiguration((root/"simulation_config.json").string());
	const auto& definition = iga::FirstNavierStokesSystem(configuration);
	const auto boundaries = iga::ResolveFlowBoundaries(configuration, definition, mesh.labels, velocity);
	iga::TransientFlowRuntime flow(database, communicator, true, true,
		{definition.density, definition.viscosity, kDt}, boundaries, mesh.labels,
		std::move(velocity), iga::WallTraceBasis(database, mesh, 0), {});
	// Runtime construction copies its inputs. Mutate the owned test storage,
	// not the caller's now-independent vector.
	auto* reference_velocity = flow.ReferenceVelocityForTesting().data();
	const auto ports = Ports();
	const std::vector<std::string> fields{"red", "blue"};
	std::vector<double> species(2*flow.RequiredNodes().size());
	for (std::size_t i = 0; i < species.size(); ++i) species[i] = i%2 == 0 ? 2.0 : 3.0;
	// Manufactured affine pressure and constant velocity are represented
	// exactly by this cubic element, independent of the partition.
	PetscInt begin = 0, end = 0;
	VecGetOwnershipRange(flow.State(), &begin, &end);
	PetscScalar* values = nullptr;
	VecGetArray(flow.State(), &values);
	for (PetscInt row = begin; row < end; ++row) {
		const auto node = row/4;
		values[row-begin] = row%4 == 3
			? 5.0+(node%4)/3.0+2.0*((node/4)%4)/3.0+3.0*(node/16)/3.0
			: 1.0+row%4;
	}
	VecRestoreArray(flow.State(), &values);
	const auto check_summary = [&] {
		const auto summary = flow.Summary();
		// 64 coefficients: velocity=(1,2,3), pressure=5+x+2y+3z,
		// tensor grid coordinates 0,1/3,2/3,1 have variance 5/36.
		const double velocity_squared = 64*14.0;
		const double pressure_squared = 64*(64.0+14.0*5.0/36.0);
		CheckValue(std::sqrt(velocity_squared), summary.velocity_l2);
		CheckValue(std::sqrt(pressure_squared), summary.pressure_l2);
		CheckValue(std::sqrt(velocity_squared+pressure_squared), summary.state_l2);
		if (summary.linear_iterations != 0) throw std::runtime_error("query changed iteration count");
	};
	check_summary();
	const auto reference_inlet = flow.ReferenceBoundaryFlow(1);
	if (!(reference_inlet < -1e-4)) throw std::runtime_error("reference inlet must be nonzero");
	CheckValue(0.0, flow.ReferenceBoundaryFlow(99));
	const auto measure = [&] { return flow.MeasurePorts(ports, kDt, fields, species); };
	const auto baseline = measure();
	for (int label = 1; label <= 3; ++label) {
		const auto& state = baseline.at("port"+std::to_string(label));
		CheckValue(label == 1 ? -1.0 : label == 2 ? 1.0 : 2.0, state.outward_flow_m3_s.value());
		CheckValue(label == 1 ? 7.5 : label == 2 ? 8.5 : 9.0, state.mean_pressure_pa.value());
	}
	for (const auto& entry : baseline) {
		CheckValue(1.0, entry.second.area_m2.value());
		CheckValue(2.0, entry.second.concentration.at("red"));
		CheckValue(3.0, entry.second.concentration.at("blue"));
		CheckValue(2.0*entry.second.outward_flow_m3_s.value(), entry.second.outward_species_flux.at("red"));
		CheckValue(3.0*entry.second.outward_flow_m3_s.value(), entry.second.outward_species_flux.at("blue"));
	}
	int tested = 0;
	if (faults) {
		iga_test::ExpectStringStreamFailure(communicator, "flow port measurement preparation", measure);
		Compare(baseline, measure());
		tested += 3;
		if (ranks > 1) {
			ExpectFailure(communicator, "flow reference label agreement", [&] {
				flow.ReferenceBoundaryFlow(rank == ranks-1 ? 2 : 1);
			});
			++tested;
		}
		std::int32_t saved_node = 0;
		if (rank == 0) {
			auto& element = const_cast<iga::Element&>(flow.OwnedElements().at(0));
			saved_node = element.connectivity[0]; element.connectivity[0] = -1;
		}
		ExpectFailure(communicator, "flow reference integration", [&] { flow.ReferenceBoundaryFlow(1); });
		if (rank == 0) const_cast<iga::Element&>(flow.OwnedElements().at(0)).connectivity[0] = saved_node;
		CheckValue(reference_inlet, flow.ReferenceBoundaryFlow(1)); ++tested;
		const auto saved_velocity = reference_velocity[0];
		if (rank == 0) reference_velocity[0][0] = std::numeric_limits<double>::quiet_NaN();
		ExpectFailure(communicator, "flow reference result", [&] { flow.ReferenceBoundaryFlow(1); });
		if (rank == 0) reference_velocity[0] = saved_velocity;
		CheckValue(reference_inlet, flow.ReferenceBoundaryFlow(1)); ++tested;
		for (int mode = 0; mode < 2; ++mode) {
			PetscScalar saved = 0.0;
			// Writable array access/restore is logically collective, even though
			// only one rank changes a value. Keep Vec norm cache states aligned.
			VecGetArray(flow.State(), &values);
			if (rank == ranks-1) {
				saved = values[0];
				values[0] = mode == 0 ? std::numeric_limits<double>::quiet_NaN()
					: std::numeric_limits<double>::max();
			}
			VecRestoreArray(flow.State(), &values);
			ExpectFailure(communicator, mode == 0 ? "flow summary local state" : "flow summary result",
				[&] { flow.Summary(); });
			VecGetArray(flow.State(), &values);
			if (rank == ranks-1) values[0] = saved;
			VecRestoreArray(flow.State(), &values);
			check_summary(); ++tested;
		}
		for (int mode = 0; mode < 17; ++mode) {
			if (ranks == 1 && mode >= 5 && mode <= 14) continue;
			auto query_ports = ports;
			auto query_fields = fields;
			auto query_species = species;
			double time = kDt;
			iga::CompiledLinearSystem system;
			system.fields = fields;
			system.terms = {{iga::TermKind::Advection, 0, 0, 1.0, ""},
				{iga::TermKind::Advection, 1, 1, 1.0, ""}};
			if (rank == ranks-1) {
				switch (mode) {
				case 0: time = std::numeric_limits<double>::quiet_NaN(); break;
				case 1: query_ports.clear(); break;
				case 2: query_ports[0].locator = "bad"; break;
				case 3: query_species.push_back(0.0); break;
				case 4: system.fields[0] = "wrong"; break;
				case 5: query_ports.pop_back(); break;
				case 6: std::swap(query_ports[0], query_ports[1]); break;
				case 7: query_ports[0].id = "different"; break;
				case 8: query_ports[0].orientation.native_to_outward_sign = -1; break;
				case 9: time *= 2; break;
				case 10: query_fields[0] = "green"; system.fields = query_fields; break;
				case 11: system.terms[0].coefficient = 2.0; break;
				case 12: query_ports[0].locator = "4"; break;
				case 13:
					query_fields.pop_back(); system.fields = query_fields;
					query_species.resize(flow.RequiredNodes().size()); system.terms.resize(1);
					break;
				case 14: system.terms.push_back({iga::TermKind::Diffusion, 0, 1, 0.1, ""}); break;
				default: break;
				}
			}
			// Force a genuine local integration exception while a PETSc read
			// view is borrowed. Mutate and restore only this test-owned runtime.
			std::int32_t saved = 0;
			if (mode == 15 && rank == 0) {
				auto& element = const_cast<iga::Element&>(flow.OwnedElements().at(0));
				saved = element.connectivity[0]; element.connectivity[0] = -1;
			}
			if (mode == 16) query_ports[0].locator = "99";
			const std::string stage = mode < 5 ? "flow port measurement preparation"
				: mode <= 14 ? "flow port measurement agreement"
				: mode == 15 ? "flow port measurement integration" : "flow port measurement result";
			ExpectFailure(communicator, stage, [&] {
				flow.MeasurePorts(query_ports, time, query_fields, query_species, &system);
			});
			if (mode == 15 && rank == 0)
				const_cast<iga::Element&>(flow.OwnedElements().at(0)).connectivity[0] = saved;
			Compare(baseline, measure());
			++tested;
		}
		// VCA conversion traverses the same collective measurement boundary.
		iga::ThreeDVascularPortDefinition vca;
		vca.inlet_label = 1; vca.outlet_labels = {2, 3};
		const auto measured = flow.MeasurePorts(vca, fields, species);
		for (int label = 1; label <= 3; ++label) {
			const auto& reference = baseline.at("port"+std::to_string(label));
			CheckValue(reference.outward_flow_m3_s.value(), measured.flows.at(label));
			CheckValue(reference.mean_pressure_pa.value(), measured.pressures.at(label));
		}
	}
	// Exercise adapter validation, actual solve, failed lookup, and rollback.
	VecSet(flow.State(), 0.0);
	flow.InitializeState(configuration);
	auto adapter_ports = ports;
	adapter_ports[1].requires = {iga::PortQuantity::MeanPressure};
	iga::ThreeDBodyFittedFlowDomainAdapter adapter("flow", flow, adapter_ports,
		configuration, root, {}, {12, 1e-8, 1e-12, 1e-6});
	adapter.BeginStep({0, 0.0, kDt});
	iga::PortBoundaryData input;
	input.time_s = kDt; input.mean_pressure_pa = 0.0;
	if (faults) {
		auto bad = input;
		if (rank == ranks-1) bad.time_s = 2*kDt;
		ExpectFailure(communicator, "3d flow adapter input", [&] { adapter.SetPortInput("port2", bad); });
		// Valid ranks must not publish the rejected input: no rank now has it.
		ExpectFailure(communicator, "3d flow adapter solve preparation: rank 0", [&] { adapter.SolveTrial(); });
		tested += 2;
	}
	adapter.SetPortInput("port2", input);
	adapter.SolveTrial();
	if (faults) {
		ExpectFailure(communicator, "3d flow adapter port lookup", [&] {
			adapter.GetPortState(rank == ranks-1 ? "unknown" : "port1");
		});
		++tested;
	}
	Measurements solved;
	for (const auto& port : ports) solved.emplace(port.id, adapter.GetPortState(port.id));
	if (std::abs(solved.at("port1").outward_flow_m3_s.value()) < 1e-4)
		throw std::runtime_error("adapter fixture must have nonzero inlet flow");
	adapter.RollbackTrial();
	adapter.SolveTrial();
	for (const auto& port : ports) Compare({{port.id, solved.at(port.id)}}, {{port.id, adapter.GetPortState(port.id)}});
	adapter.PrepareCommitStep(); adapter.FinalizeCommitStep();
	if (rank == 0 && faults)
		std::cout << "three_d_port_failure ranks=" << ranks << " modes=" << tested << " retry=passed\n";
	return solved;
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int world_rank = 0, world_size = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank); MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		if (argc != 2 || world_size != 3)
			throw std::runtime_error("usage: mpiexec -np 3 three_d_port_failure_test FRESH_OUTPUT");
		for (int round = 0; round < 2; ++round) {
			const int color = round == 0 || world_rank == 0 ? 0 : 1;
			MPI_Comm_split(PETSC_COMM_WORLD, color, world_rank, &group);
			int rank = 0, ranks = 0;
			MPI_Comm_rank(group, &rank); MPI_Comm_size(group, &ranks);
			const auto root = fs::path(argv[1])/(std::to_string(round)+"-"+std::to_string(color));
			iga::CollectiveLocalStage(group, "test fixture", [&] {
				if (rank == 0) {
					if (!fs::create_directories(root)) throw std::runtime_error("output already exists");
					WriteThreeDCase(root); WriteDatabase(root/"group.ntiga", ranks); WriteDatabase(root/"serial.ntiga", 1);
				}
			});
			Measurements reference;
			iga::CollectiveLocalStage(group, "test serial reference", [&] {
				// COMM_SELF runtime collectives are local to this process.
				if (rank == 0) reference = Run(root, root/"serial.ntiga", PETSC_COMM_SELF, false);
			});
			const auto actual = Run(root, root/"group.ntiga", group, true);
			iga::CollectiveLocalStage(group, "test reference comparison", [&] {
				if (rank == 0) Compare(reference, actual);
			});
			MPI_Comm_free(&group);
		}
		if (world_rank == 0) std::cout << "three_d_port_failure serial/group comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << world_rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD, 1);
		return 1;
	}
	PetscFinalize();
	return 0;
}

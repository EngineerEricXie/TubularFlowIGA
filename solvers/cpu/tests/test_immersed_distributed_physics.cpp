#include "ImmersedStaticDistributedOperator.hpp"
#include "ImmersedStaticFlowRuntime.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
iga::CartesianDomainClassification Domain(bool padded)
{
	iga::RawSurfaceSoup soup;
	const double a = padded ? 0.3 : 0.1, b = padded ? 0.7 : 0.9;
	soup.vertices = {{{{a,a,a}},{{b,a,a}},{{b,b,a}},{{a,b,a}},{{a,a,b}},{{b,a,b}},{{b,b,b}},{{a,b,b}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	const unsigned cells = padded ? 5 : 3;
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{cells,cells,cells}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(grid),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	Fixture(bool padded,bool expanded)
		: domain(Domain(padded)),volume(domain,{2,500000,500000,3000000},expanded ? iga::CutCellVolumeQuadratureStorageMode::Expanded : iga::CutCellVolumeQuadratureStorageMode::Compact) {}
	iga::CartesianDomainClassification domain;
	iga::CutCellVolumeQuadratureCatalog volume;
	iga::ImmersedSurfaceQuadratureCatalog surface{domain};
	iga::CutCellGhostPenaltyCatalog ghost{domain, volume};
};

void Run(MPI_Comm comm, const std::string& mode,iga::ImmersedWorkPartition partition,int repetitions)
{
	int rank = 0, size = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
	if (mode != "flow" && mode != "pressure" && mode != "traction" && mode != "closed" && mode != "wall-only" && mode != "faults" && mode != "padded" && mode != "expanded") throw std::invalid_argument("unknown physics mode");
	bool inject = false;
	std::size_t force_calls = 0;
	std::unique_ptr<Fixture> fixture;
	std::unique_ptr<iga::ImmersedStaticFlowRuntime> serial;
	iga::ImmersedStaticFlowOptions options;
	std::vector<PetscScalar> state, expected_residual, expected_action;
	// Each validation process builds its own COMM_SELF reference. This is
	// reference-only; the production operator owns numerical integration and
	// never calls the reference runtime to obtain a block.
	std::exception_ptr setup_error;
	try {
		fixture = std::make_unique<Fixture>(mode == "padded",mode == "expanded");
		options.wall_labels = {7};
		options.ports = {{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-0.001},
			{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,0.001}};
		if (mode == "pressure") options.ports[1] = {"outlet",9,iga::ImmersedFlowPortControlMode::Pressure,0.003};
		if (mode == "traction") {
			options.ports[0] = {"inlet",8,iga::ImmersedFlowPortControlMode::MeanNormalTraction,-0.004};
			options.ports[1] = {"outlet",9,iga::ImmersedFlowPortControlMode::Pressure,0.003};
		}
		if (mode == "closed") { options.ports.clear(); options.wall_labels = {7,8,9}; }
		if (mode == "wall-only") options.assemble_volume = options.assemble_ghost = options.assemble_gauge = false;
		options.body_force = [&](const std::array<double,3>& x) {
			if (inject && rank == size-1 && ++force_calls > 3000) throw std::runtime_error("injected owned volume evaluator failure");
			return std::array<double,3>{{0.1+x[0],-0.2*x[1],0.3*x[2]}}; };
		serial = std::make_unique<iga::ImmersedStaticFlowRuntime>(fixture->domain, fixture->volume, fixture->surface, fixture->ghost, options);
		state.resize(serial->Diagnostics().total_dofs);
		for (std::size_t i = 0; i < state.size(); ++i) state[i] = 0.01*std::sin(0.13*(i+1));
		serial->SetCommittedState(state); serial->Assemble();
		expected_residual = serial->AssembledNegativeResidual(); expected_action = serial->AssembledJacobianAction(state);
	} catch (...) { setup_error = std::current_exception(); }
	iga::CollectiveLocalStage(comm, "immersed physics reference setup", [&] { if (setup_error) std::rethrow_exception(setup_error); });
	auto& f = *fixture;
	if (mode == "faults") {
		for (int scenario = 0; scenario < 4; ++scenario) {
			if (size == 1 && scenario%2) continue;
			auto invalid = options;
			auto invalid_partition = partition;
			if (rank == size-1) {
				if (scenario < 2) invalid.wall_gamma0 = scenario == 0 ? -1.0 : 3.0;
				else invalid_partition = scenario == 2 ? static_cast<iga::ImmersedWorkPartition>(-1)
					: partition == iga::ImmersedWorkPartition::CellCount ? iga::ImmersedWorkPartition::WeightedContiguous : iga::ImmersedWorkPartition::CellCount;
			}
			bool rejected = false;
			try { iga::ImmersedStaticDistributedOperator bad(comm,f.domain,f.volume,f.surface,f.ghost,invalid,invalid_partition); }
			catch (const std::exception&) { rejected = true; }
			iga::CollectiveLocalStage(comm,"physics invalid collective setup",[&] { if (!rejected) throw std::runtime_error("invalid or inconsistent controls accepted"); });
		}
	}
	iga::ImmersedStaticDistributedOperator op(comm,f.domain,f.volume,f.surface,f.ghost,options,partition);
	auto& distributed = op.Assembly();
	iga::CollectiveLocalStage(comm, "immersed physics state insertion", [&] {
		for (PetscInt row = distributed.RowBegin(); row < distributed.RowEnd(); ++row)
			if (VecSetValue(distributed.State(),row,state[static_cast<std::size_t>(row)],INSERT_VALUES)) throw std::runtime_error("cannot insert physical state");
	});
	iga::RequireCollectivePetscSuccess(comm,"physics state begin",VecAssemblyBegin(distributed.State()));
	iga::RequireCollectivePetscSuccess(comm,"physics state end",VecAssemblyEnd(distributed.State()));
	if (mode == "faults") {
		inject = true;
		bool rejected = false;
		try { op.Assemble(); } catch (const std::exception&) { rejected = true; }
		inject = false;
		iga::CollectiveLocalStage(comm,"physics evaluator failure agreement",[&] {
			if (!rejected || op.Diagnostics().aggregate_assembly_seconds != 0.0) throw std::runtime_error("failed operator published diagnostics");
			iga::PetscReadArray values; values.Acquire(distributed.State());
			for (PetscInt row = distributed.RowBegin(); row < distributed.RowEnd(); ++row)
				if (values.Data()[row-distributed.RowBegin()] != state[row]) throw std::runtime_error("failed operator changed state");
			values.Restore();
		});
	}
	double assembly_seconds = 0.0;
	for (int repeat = 0; repeat < repetitions; ++repeat) {
		MPI_Barrier(comm);
		const auto start = std::chrono::steady_clock::now();
		op.Assemble();
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
		assembly_seconds += seconds;
		const auto& timing = distributed.LastAssemblyTiming();
		std::cout << std::setprecision(17) << "immersed_work_sample rank=" << rank << " repeat=" << repeat
			<< " assembly_s=" << seconds << " integration_insert_s=" << timing.local_integration_insert_seconds
			<< " halo_s=" << timing.halo_seconds << " stash_s=" << timing.stash_exchange_seconds << '\n';
	}
	assembly_seconds /= repetitions;
	Vec action = nullptr;
	iga::RequireCollectivePetscSuccess(comm,"physics action create",VecDuplicate(distributed.State(),&action));
	iga::RequireCollectivePetscSuccess(comm,"physics matrix action",MatMult(distributed.Matrix(),distributed.State(),action));
	double local[4]{}, global[4]{};
	iga::CollectiveLocalStage(comm,"immersed physics compare",[&] {
		iga::PetscReadArray rhs, product; rhs.Acquire(distributed.Residual()); product.Acquire(action);
		for (PetscInt row = distributed.RowBegin(); row < distributed.RowEnd(); ++row) {
			const auto offset = row-distributed.RowBegin();
			const double r = PetscRealPart(rhs.Data()[offset])-PetscRealPart(expected_residual[row]);
			const double a = PetscRealPart(product.Data()[offset])-PetscRealPart(expected_action[row]);
			local[0] += r*r; local[1] += PetscRealPart(expected_residual[row])*PetscRealPart(expected_residual[row]);
			local[2] += a*a; local[3] += PetscRealPart(expected_action[row])*PetscRealPart(expected_action[row]);
		}
		rhs.Restore(); product.Restore();
	});
	MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_SUM,comm);
	const double residual_error = std::sqrt(global[0]/global[1]), action_error = std::sqrt(global[2]/global[3]);
	iga::CollectiveLocalStage(comm,"immersed physics numerical gates",[&] {
		if (!(residual_error <= 1e-6) || !(action_error <= 1e-6)) throw std::runtime_error("distributed physical operator differs from serial runtime");
	});
	iga::CollectiveLocalStage(comm,"immersed physics diagnostic parity",[&] {
		const auto& actual = op.Diagnostics(); const auto& expected = serial->Diagnostics();
		const auto close = [](double a,double b) {
			if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b) > 1e-11*std::max({1.0,std::abs(a),std::abs(b)}))
				throw std::runtime_error("distributed physical diagnostic differs from serial");
		};
		if (actual.active_nodes != expected.active_nodes || actual.physical_dofs != expected.physical_dofs || actual.total_dofs != expected.total_dofs
			|| actual.volume_cells != expected.volume_cells || actual.surface_cells != expected.surface_cells || actual.ghost_faces != expected.ghost_faces
			|| actual.gauge_present != expected.gauge_present || actual.gauge_row != expected.gauge_row || !actual.scalar_diagonal_structure_verified
			|| actual.wall_selected_points != expected.wall_selected_points || actual.wall_selected_points_by_label != expected.wall_selected_points_by_label)
			throw std::runtime_error("distributed physical topology/count diagnostic differs from serial");
		close(actual.pressure_measure,expected.pressure_measure); close(actual.constant_pressure_defect,expected.constant_pressure_defect);
		for (std::size_t p = 0; p < actual.ports.size(); ++p) {
			const auto& a = actual.ports[p]; const auto& b = expected.ports[p];
			if (a.id != b.id || a.multiplier_row != b.multiplier_row || a.assembled_surface_points != b.assembled_surface_points)
				throw std::runtime_error("distributed port topology differs from serial");
			close(a.measurement.area_m2,b.measurement.area_m2); close(a.measurement.outward_flow_m3_s,b.measurement.outward_flow_m3_s);
			close(a.measurement.mean_pressure_pa,b.measurement.mean_pressure_pa); close(a.measurement.mean_normal_traction_pa,b.measurement.mean_normal_traction_pa);
			close(a.measurement.mean_velocity_squared_m2_s2,b.measurement.mean_velocity_squared_m2_s2);
			close(a.multiplier,b.multiplier); close(a.constraint_residual,b.constraint_residual); close(a.normalized_flow_residual,b.normalized_flow_residual);
			close(a.flow_tolerance_m3_s,b.flow_tolerance_m3_s);
		}
	});
	MatInfo info{};
	iga::RequireCollectivePetscSuccess(comm,"physics matrix info",MatGetInfo(distributed.Matrix(),MAT_LOCAL,&info));
	iga::CollectiveLocalStage(comm,"immersed physics exact pattern",[&] {
		if (info.nz_used != static_cast<double>(distributed.SymbolicEntries()) || info.nz_allocated != info.nz_used || info.mallocs != 0.0)
			throw std::runtime_error("physical matrix allocated outside exact stencil graph");
	});
	iga::RequireCollectivePetscSuccess(comm,"physics action destroy",VecDestroy(&action));
	std::cout << std::setprecision(17) << "immersed_physics rank=" << rank << " ranks=" << size
		<< " cells=" << f.domain.Cells().size() << " ghost_faces=" << f.ghost.Faces().size()
		<< " owned_stencils=" << distributed.OwnedStencils().size() << " owned_rows=" << distributed.RowEnd()-distributed.RowBegin()
		<< " halo_rows=" << distributed.RequiredRows().size()
		<< " remote_halo_rows=" << std::count_if(distributed.RequiredRows().begin(),distributed.RequiredRows().end(),[&](PetscInt row) { return row < distributed.RowBegin() || row >= distributed.RowEnd(); })
		<< " global_rows=" << distributed.Rows()
		<< " estimated_work=" << op.EstimatedOwnedWork() << " weighted=" << (partition == iga::ImmersedWorkPartition::WeightedContiguous)
		<< " local_nz=" << info.nz_used << " assembly_s=" << assembly_seconds
		<< " residual_relative_l2=" << residual_error << " action_relative_l2=" << action_error << " passed\n";
	op.Close(); serial.reset();
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int status = 0;
	try {
		const std::string partition = argc > 2 ? argv[2] : "cell-count";
		if (partition != "cell-count" && partition != "weighted") throw std::invalid_argument("unknown partition");
		const int repetitions = argc > 3 ? std::stoi(argv[3]) : 1;
		if (repetitions < 1 || repetitions > 10) throw std::invalid_argument("invalid assembly repetition count");
		Run(PETSC_COMM_WORLD,argc > 1 ? argv[1] : "flow",partition == "weighted" ? iga::ImmersedWorkPartition::WeightedContiguous : iga::ImmersedWorkPartition::CellCount,repetitions);
	}
	catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	PetscFinalize(); return status;
}

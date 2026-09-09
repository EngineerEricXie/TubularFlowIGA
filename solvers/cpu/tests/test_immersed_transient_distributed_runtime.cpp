#include "ImmersedTransientDistributedRuntime.hpp"
#include "PrescribedSurfaceMotion.hpp"
#include <iostream>
#include <memory>

namespace {
iga::RawSurfaceSoup Cube()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	return soup;
}
void Require(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}
template<class Function> void Reject(MPI_Comm comm, Function&& function)
{
	std::string message;
	try { function(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm,"history test rejection",[&] { Require(!message.empty(),"missing history rejection"); });
	iga::RequireCollectiveSameText(comm,"history rejection agreement",message);
}
std::unique_ptr<iga::MovingCutGeometry> Geometry(double time,bool empty)
{
	const auto soup = Cube(); iga::PrescribedSurfaceMotion motion({{0,soup},{1,soup}});
	iga::MovingCutGeometryOptions options; options.volume.max_depth = 2;
	options.volume_storage = iga::CutCellVolumeQuadratureStorageMode::Compact;
	return iga::MovingCutGeometry::Build({{{0,0,0}},{{1,1,1}},{{empty ? 2u : 4u,1,1}}},motion.Evaluate(time,0,1),options);
}
std::vector<PetscScalar> Owned(Vec vector)
{
	PetscInt size = 0; if (VecGetLocalSize(vector,&size)) throw std::runtime_error("cannot inspect owned vector");
	iga::PetscReadArray values; values.Acquire(vector);
	std::vector<PetscScalar> copy;
	if (size) copy.assign(values.Data(),values.Data()+size);
	values.Restore(); return copy;
}
void Run(MPI_Comm comm,const std::string& mode,bool weighted)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	const bool closed = mode == "closed" || mode == "empty-work",empty = mode == "empty-work";
	std::unique_ptr<iga::MovingCutGeometry> geometry;
	iga::ImmersedTransientFlowOptions options;
	const auto force = [closed](const std::array<double,3>&) { return closed ? std::array<double,3>{{.3,-.2,.1}} : std::array<double,3>{}; };
	bool inject = false,forbid = false; std::size_t calls = 0;
	iga::CollectiveLocalStage(comm,"transient runtime fixture",[&] {
		if (mode != "flow" && mode != "pressure" && mode != "closed" && mode != "inertial" && !empty) throw std::invalid_argument("unknown transient runtime mode");
		geometry = Geometry(0,empty); options.parameters = {1,1,0}; options.wall_labels = closed ? std::vector<int>{7,8,9} : std::vector<int>{7};
		options.nonlinear_absolute_tolerance = 1e-14; options.nonlinear_relative_tolerance = 1e-11;
		options.ksp_relative_tolerance = 1e-11; options.lu_pivot_shift = 1e-12;
		if (!closed) options.ports = {{"outlet",9,mode == "pressure" ? iga::ImmersedFlowPortControlMode::Pressure : iga::ImmersedFlowPortControlMode::FlowRate,mode == "pressure" ? .001 : .0001},
			{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-.0001}};
		if (mode == "inertial") options.wall_inertial_gamma0 = .6;
		options.body_force = [&](const std::array<double,3>& x) {
			++calls;
			if (forbid || (inject && rank == 0 && calls == 2)) throw std::runtime_error("injected transient runtime force failure");
			return force(x);
		};
	});
	auto& g = *geometry;
	iga::ImmersedTransientDistributedRuntime runtime(comm,g.Domain(),g.Volume(),g.Surface(),g.Ghost(),g.GeometryIdentitySha256(),options,
		weighted ? iga::ImmersedWorkPartition::WeightedContiguous : iga::ImmersedWorkPartition::CellCount);
	std::vector<PetscScalar> previous(runtime.Layout().Rows(),0),local;
	iga::CollectiveLocalStage(comm,"transient runtime initial data",[&] {
		if (closed) for (std::size_t row = 0; row < runtime.Layout().NodeFieldRows(); ++row)
			if (row%4 < 3) previous[row] = .0001*std::sin(.37*(row+1));
		local.assign(previous.begin()+runtime.RowBegin(),previous.begin()+runtime.RowEnd());
	});
	runtime.SetCommittedOwnedState(local,0,0);
	Reject(comm,[&] { runtime.BeginTrial(.125,rank == size-1 ? 2 : 1,.125); });
	if (size > 1) Reject(comm,[&] { runtime.BeginTrial(rank == size-1 ? .25 : .125,1,rank == size-1 ? .25 : .125); });
	inject = true; Reject(comm,[&] { runtime.BeginTrial(.125,1,.125); }); inject = false;
	iga::CollectiveLocalStage(comm,"transient failed begin publication",[&] {
		Require(!runtime.Clock().trial_active && !runtime.History().Active() && runtime.Clock().index == 0,"failed begin published transient clock");
		Require(Owned(runtime.CommittedState()) == local,"failed begin changed committed field");
	});
	double maximum_field_error = 0,physical_fraction = 0,history_difference = 0;
	const auto serial_solve = [&](double time,std::uint64_t index,const std::vector<PetscScalar>& old,bool zero_velocity) {
		auto target_geometry = Geometry(time,empty);
		auto settings = options; settings.body_force = force;
		if (!closed) {
			for (auto& p : settings.ports) if (p.control_mode == iga::ImmersedFlowPortControlMode::FlowRate)
				p.value *= index == 1 ? 1 : .5;
		}
		auto reference = std::make_unique<iga::ImmersedTransientFlowRuntime>(*target_geometry,settings);
		const auto& layout = reference->Layout(); std::vector<std::array<double,4>> fields(layout.NodeIds().size());
		for (std::size_t i = 0; i < fields.size(); ++i) for (int c = 0; c < 4; ++c) fields[i][c] = zero_velocity && c < 3 ? 0 : old[4*i+c];
		std::vector<double> multipliers(layout.PortIds().size());
		for (std::size_t i = 0; i < multipliers.size(); ++i) multipliers[i] = old[layout.NodeFieldRows()+i];
		reference->SetCommittedGlobalState(iga::ImmersedGlobalFlowState(time-.125,index-1,layout,fields,multipliers,layout.HasGaugeRow(),layout.HasGaugeRow() ? old[layout.GaugeRow()] : 0));
		reference->BeginTrial(time,index,.125); Require(reference->SolveTrial(),"serial transient solve failed"); reference->Commit();
		// Keep the immutable catalogs alive through reference destruction.
		return std::make_pair(std::move(target_geometry),std::move(reference));
	};
	for (std::uint64_t step = 1; step <= 2; ++step) {
		const double time = .125*step;
		if (step == 2 && !closed) for (const auto& p : options.ports)
			if (p.control_mode == iga::ImmersedFlowPortControlMode::FlowRate) runtime.SetPortControlValue(p.id,.5*p.value);
		std::pair<std::unique_ptr<iga::MovingCutGeometry>,std::unique_ptr<iga::ImmersedTransientFlowRuntime>> serial;
		iga::CollectiveLocalStage(comm,"serial transient step",[&] { serial = serial_solve(time,step,previous,false); });
		calls = 0; forbid = false; runtime.BeginTrial(time,step,.125); forbid = true;
		const auto frozen_calls = calls;
		Reject(comm,[&] { runtime.SetCommittedOwnedState(local,0,0); });
		if (!closed) Reject(comm,[&] { runtime.SetPortControlValue("inlet",0); });
		if (step == 1) {
			if (rank == size-1) runtime.FailNextCandidateForTesting();
			Reject(comm,[&] { runtime.SolveTrial(); });
			iga::CollectiveLocalStage(comm,"transient Newton failure rollback",[&] {
				Require(runtime.Clock().trial_active && runtime.History().Active() && runtime.Clock().index == 0,"Newton failure lost frozen trial or advanced clock");
				Require(Owned(runtime.State()) == local && Owned(runtime.CommittedState()) == local,"Newton failure changed committed field");
				Require(runtime.Diagnostics().rollback_count > 0 && runtime.Diagnostics().ksp_iterations > 0,"candidate failure preceded linear update");
			});
		}
		Require(runtime.SolveTrial(),"distributed transient solve failed");
		if (step == 1) {
			if (rank == size-1) runtime.FailNextPrepareForTesting();
			Reject(comm,[&] { runtime.PrepareCommit(); }); runtime.FinalizeCommit();
			iga::CollectiveLocalStage(comm,"transient failed prepare clock",[&] { Require(runtime.Clock().index == 0 && !runtime.Diagnostics().prepared,"failed prepare published"); });
			runtime.PrepareCommit(); runtime.AbortPrepared(); runtime.FinalizeCommit();
			iga::CollectiveLocalStage(comm,"transient cancelled prepare clock",[&] { Require(runtime.Clock().index == 0,"cancelled preparation advanced clock"); });
			runtime.PrepareCommit(); runtime.AbortTrial(); runtime.FinalizeCommit();
			iga::CollectiveLocalStage(comm,"transient aborted trial",[&] {
				Require(!runtime.Clock().trial_active && !runtime.History().Active() && runtime.Clock().index == 0,"abort published or retained history");
				Require(Owned(runtime.State()) == local,"abort did not restore field");
			});
			forbid = false; runtime.BeginTrial(time,step,.125); forbid = true; Require(runtime.SolveTrial(),"retry after abort failed");
		} else {
			iga::CollectiveLocalStage(comm,"transient frozen callbacks",[&] { Require(calls == frozen_calls,"Newton reevaluated frozen force"); });
		}
		const auto expected = serial.second->CommittedState();
		std::array<double,6> sums{},global{};
		iga::CollectiveLocalStage(comm,"transient field parity",[&] {
			const auto actual = Owned(runtime.State());
			for (PetscInt row = runtime.RowBegin(); row < runtime.RowEnd(); ++row) {
				const int block = static_cast<std::size_t>(row) >= runtime.Layout().NodeFieldRows() ? 2 : row%4 == 3 ? 1 : 0;
				const double difference = actual[row-runtime.RowBegin()]-expected[row];
				sums[2*block] += difference*difference; sums[2*block+1] += expected[row]*expected[row];
			}
			Require(runtime.Diagnostics().converged && runtime.Diagnostics().ksp_reason > 0,"distributed transient convergence missing");
			for (const auto& iteration : runtime.Diagnostics().newton_steps)
				Require(std::isfinite(iteration.linear_relative_residual) && iteration.linear_relative_residual < 1e-6,"true linear residual failed");
		});
		MPI_Allreduce(sums.data(),global.data(),6,MPI_DOUBLE,MPI_SUM,comm);
		iga::CollectiveLocalStage(comm,"transient field gate",[&] {
			for (int block = 0; block < 3; ++block) {
				const double error = std::sqrt(global[2*block]),norm = std::sqrt(global[2*block+1]);
				if (norm < 1e-10) Require(error < 1e-10,"zero transient field differs");
				else { const double relative = error/norm; Require(relative < 1e-6,"transient field relative L2 failed"); maximum_field_error = std::max(maximum_field_error,relative); }
			}
		});
		const auto actual_conservation = runtime.ConservationDiagnostics();
		iga::CollectiveLocalStage(comm,"transient physical parity",[&] {
			const auto check = [&](double actual,double expected_value) {
				const double fraction = std::abs(actual-expected_value)/(1e-11+1e-6*std::abs(expected_value));
				Require(std::isfinite(fraction) && fraction <= 1,"transient physical diagnostic differs"); physical_fraction = std::max(physical_fraction,fraction);
			};
			const auto& ports = serial.second->Diagnostics().ports;
			for (std::size_t i = 0; i < ports.size(); ++i) {
				const auto& a = runtime.Diagnostics().ports[i].measurement; const auto& e = ports[i].measurement;
				check(a.outward_flow_m3_s,e.outward_flow_m3_s); check(a.mean_pressure_pa,e.mean_pressure_pa); check(a.mean_normal_traction_pa,e.mean_normal_traction_pa);
			}
			const auto e = serial.second->ConservationDiagnostics();
			check(actual_conservation.volume_divergence_integral_m3_s,e.endpoint_volume_divergence_m3_s);
			check(actual_conservation.total_surface_outward_flow_m3_s,e.total_surface_outward_flow_m3_s);
			check(actual_conservation.wall_outward_flow_m3_s,e.wall_outward_flow_m3_s);
		});
		runtime.PrepareCommit();
		iga::CollectiveLocalStage(comm,"transient preparation is unpublished",[&] { Require(runtime.Clock().index == step-1 && runtime.Clock().time_s == time-.125,"prepare advanced clock"); });
		runtime.FinalizeCommit(); runtime.FinalizeCommit();
		iga::CollectiveLocalStage(comm,"transient accepted clock",[&] {
			Require(runtime.Clock().index == step && runtime.Clock().time_s == time && !runtime.Clock().trial_active && !runtime.History().Active()
				&& runtime.Diagnostics().commit_count == step,"accepted clock/state publication differs");
		});
		if (step == 2) iga::CollectiveLocalStage(comm,"transient history discriminator",[&] {
			auto cold = serial_solve(time,step,previous,true); const auto reset = cold.second->CommittedState();
			for (std::size_t row = 0; row < runtime.Layout().NodeFieldRows(); ++row) if (row%4 < 3)
				history_difference = std::max(history_difference,std::abs(expected[row]-reset[row]));
			Require(history_difference > 1e-9,"second step does not discriminate committed velocity history");
		});
		previous = expected;
	}
	forbid = false; runtime.BeginTrial(.375,3,.125); forbid = true; Require(runtime.SolveTrial(),"third trial failed"); runtime.PrepareCommit();
	const auto diagnostics = runtime.Diagnostics(); const auto owned = runtime.History().OwnedNodes();
	runtime.Close(); runtime.FinalizeCommit(); runtime.Close();
	iga::CollectiveLocalStage(comm,"transient closed finalize",[&] { Require(runtime.Clock().index == 2 && runtime.Clock().time_s == .25,"close finalized uncommitted clock"); });
	Reject(comm,[&] { runtime.BeginTrial(.375,3,.125); });
	std::cout << "immersed_transient_runtime_mpi rank=" << rank << " ranks=" << size << " steps=2 owned_rows=" << runtime.RowEnd()-runtime.RowBegin()
		<< " owned_history_nodes=" << owned << " field_relative_l2=" << maximum_field_error << " physical_fraction=" << physical_fraction
		<< " history_difference=" << history_difference << " assembly_s=" << diagnostics.aggregate_assembly_seconds
		<< " solve_s=" << diagnostics.aggregate_linear_solve_seconds << " passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0,status = 0; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm comm = MPI_COMM_NULL;
	try {
		bool split = false,weighted = false;
		for (int i = 2; i < argc; ++i) { if (std::string(argv[i]) == "split") split = true; else if (std::string(argv[i]) == "weighted") weighted = true; else throw std::invalid_argument("unknown runtime test flag"); }
		MPI_Comm_split(PETSC_COMM_WORLD,split && rank ? 1 : 0,rank,&comm);
		Run(comm,argc > 1 ? argv[1] : "flow",weighted);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

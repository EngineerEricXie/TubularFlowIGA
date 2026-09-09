#include "ImmersedTransientDistributedOperator.hpp"
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
double Value(std::int32_t id, int field, double scale)
{
	return scale*(0.02+0.003*std::sin(0.31*(id+1)*(field+1)));
}
void Run(MPI_Comm comm,const std::string& mode)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	std::unique_ptr<iga::MovingCutGeometry> geometry;
	std::unique_ptr<iga::ImmersedTransientFlowRuntime> reference;
	iga::ImmersedTransientFlowOptions options;
	std::vector<PetscScalar> committed,trial,expected_residual,expected_action;
	iga::ImmersedTransientFlowConservationDiagnostics conservation;
	iga::CollectiveLocalStage(comm,"transient operator reference",[&] {
		if (mode != "flow" && mode != "pressure" && mode != "traction" && mode != "closed" && mode != "inertial") throw std::invalid_argument("unknown transient operator mode");
		const auto soup = Cube(); iga::PrescribedSurfaceMotion motion({{0,soup},{.125,soup}});
		iga::MovingCutGeometryOptions go; go.volume.max_depth = 2; go.volume_storage = iga::CutCellVolumeQuadratureStorageMode::Compact;
		geometry = iga::MovingCutGeometry::Build({{{0,0,0}},{{1,1,1}},{{4,1,1}}},motion.Evaluate(.125,0,.125),go);
		options.parameters = {1,.1,0}; options.wall_labels = {7};
		// Deliberately reverse declaration order: the existing transient API
		// assigns controllers by numeric label, not this vector order.
		options.ports = {{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,.001},
			{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-.001}};
		if (mode == "pressure" || mode == "traction") options.ports[0].control_mode = mode == "pressure"
			? iga::ImmersedFlowPortControlMode::Pressure : iga::ImmersedFlowPortControlMode::MeanNormalTraction;
		if (mode == "closed") { options.ports.clear(); options.wall_labels = {7,8,9}; }
		if (mode == "inertial") options.wall_inertial_gamma0 = .6;
		options.body_force = [](const std::array<double,3>& x) { return std::array<double,3>{{.1+x[0],-.2*x[1],.3*x[2]}}; };
		reference = std::make_unique<iga::ImmersedTransientFlowRuntime>(*geometry,options);
		const auto& layout = reference->Layout();
		std::vector<std::array<double,4>> fields(layout.NodeIds().size());
		for (std::size_t i = 0; i < fields.size(); ++i)
			for (int c = 0; c < 4; ++c) fields[i][c] = Value(layout.NodeIds()[i],c,1);
		reference->SetCommittedGlobalState(iga::ImmersedGlobalFlowState(0,0,layout,fields,std::vector<double>(layout.PortIds().size()),layout.HasGaugeRow(),0));
		committed = reference->CommittedState(); trial.resize(layout.Rows());
		for (std::size_t i = 0; i < trial.size(); ++i) trial[i] = .01*std::sin(.13*(i+1));
		reference->BeginTrial(.125,1,.125); reference->SetTrialState(trial); reference->Assemble();
		expected_residual = reference->AssembledNegativeResidual(); expected_action = reference->AssembledJacobianAction(trial);
		conservation = reference->ConservationDiagnostics();
		if (mode == "inertial") Require(reference->Diagnostics().wall_penalty.maximum_eta_t > 0,"reference inertial wall term is absent");
	});
	auto& g = *geometry;
	auto invalid = options;
	if (rank == size-1) invalid.wall_inertial_gamma0 = -1;
	Reject(comm,[&] { iga::ImmersedTransientDistributedOperator rejected(comm,g.Domain(),g.Volume(),g.Surface(),g.Ghost(),g.GeometryIdentitySha256(),invalid); });
	if (size > 1) {
		auto different = options;
		if (rank == size-1) different.wall_inertial_gamma0 += .1;
		Reject(comm,[&] { iga::ImmersedTransientDistributedOperator rejected(comm,g.Domain(),g.Volume(),g.Surface(),g.Ghost(),g.GeometryIdentitySha256(),different); });
	}
	iga::ImmersedTransientDistributedOperator op(comm,g.Domain(),g.Volume(),g.Surface(),g.Ghost(),g.GeometryIdentitySha256(),options);
	Vec old = nullptr;
	iga::RequireCollectivePetscSuccess(comm,"transient test history vector",VecDuplicate(op.Assembly().State(),&old));
	try {
	const auto fill = [&](Vec vector,const std::vector<PetscScalar>& values) {
		iga::CollectiveLocalStage(comm,"transient test state insert",[&] {
			for (PetscInt row = op.Assembly().RowBegin(); row < op.Assembly().RowEnd(); ++row)
				Require(!VecSetValue(vector,row,values[row],INSERT_VALUES),"transient vector insertion failed");
		});
		iga::RequireCollectivePetscSuccess(comm,"transient state begin",VecAssemblyBegin(vector));
		iga::RequireCollectivePetscSuccess(comm,"transient state end",VecAssemblyEnd(vector));
	};
	fill(old,committed); fill(op.Assembly().State(),trial);
	Reject(comm,[&] { op.Assemble(); });
	op.Freeze(old,reference->Layout(),0,0,.125,1,.125);
	if (!options.ports.empty()) Reject(comm,[&] { op.SetPortControlValue(options.ports[0].id,0); });
	std::array<double,4> errors{}; double physical_error = 0;
	const auto verify = [&] {
		op.Assemble(); Vec action = nullptr;
		iga::RequireCollectivePetscSuccess(comm,"transient action create",VecDuplicate(op.Assembly().State(),&action));
		try {
		iga::RequireCollectivePetscSuccess(comm,"transient action",MatMult(op.Assembly().Matrix(),op.Assembly().State(),action));
		iga::CollectiveLocalStage(comm,"transient numerical comparison",[&] {
			iga::PetscReadArray residual,product; residual.Acquire(op.Assembly().Residual()); product.Acquire(action);
			for (PetscInt row = op.Assembly().RowBegin(); row < op.Assembly().RowEnd(); ++row) {
				const double r = residual.Data()[row-op.Assembly().RowBegin()],a = product.Data()[row-op.Assembly().RowBegin()];
				Require(std::isfinite(r) && std::isfinite(a),"nonfinite transient operator result");
				errors[0] += std::pow(r-expected_residual[row],2); errors[1] += std::pow(expected_residual[row],2);
				errors[2] += std::pow(a-expected_action[row],2); errors[3] += std::pow(expected_action[row],2);
			}
			residual.Restore(); product.Restore();
			const auto& d = op.Diagnostics(); const auto& serial = reference->Diagnostics();
			Require(d.volume_cells == serial.volume_cells && d.surface_cells == serial.surface_cells && d.ghost_faces == serial.ghost_faces,"transient work counts differ");
			Require(d.scalar_diagonal_structure_verified,"transient scalar structure missing");
			physical_error = std::max(physical_error,std::abs(d.pressure_measure-serial.pressure_measure));
			for (const auto& port : d.ports) {
				std::array<double,5> sums{};
				for (std::uint64_t c = 0; c < g.Domain().Cells().size(); ++c) {
					const auto element = g.Domain().Background().MaterializeElement(c);
					std::vector<std::array<double,4>> nodal(element.connectivity.size());
					for (std::size_t i = 0; i < nodal.size(); ++i)
						for (int f = 0; f < 4; ++f) nodal[i][f] = trial[4*op.Layout().LocalNode(element.connectivity[i])+f];
					const auto m = iga::MeasureImmersedFlowPortElement(element,g.Surface().UsableRule(g.Domain(),c),port.boundary_label,nodal,.1);
					sums[0] += m.area_m2; sums[1] += m.outward_flow_m3_s; sums[2] += m.area_m2*m.mean_pressure_pa;
					sums[3] += m.area_m2*m.mean_normal_traction_pa; sums[4] += m.area_m2*m.mean_velocity_squared_m2_s2;
				}
				const auto& m = port.measurement;
				for (double difference : {m.area_m2-sums[0],m.outward_flow_m3_s-sums[1],m.mean_pressure_pa-sums[2]/sums[0],
					m.mean_normal_traction_pa-sums[3]/sums[0],m.mean_velocity_squared_m2_s2-sums[4]/sums[0]}) physical_error = std::max(physical_error,std::abs(difference));
				if (port.multiplier_row >= 0) Require(port.multiplier_row == static_cast<PetscInt>(op.Layout().ControllerRow(port.boundary_label)),"transient controller order differs");
			}
		});
		} catch (...) { VecDestroy(&action); throw; }
		iga::RequireCollectivePetscSuccess(comm,"transient action destroy",VecDestroy(&action));
		const auto actual = op.ConservationDiagnostics();
		iga::CollectiveLocalStage(comm,"transient conservation comparison",[&] {
			for (double difference : {actual.volume_divergence_integral_m3_s-conservation.endpoint_volume_divergence_m3_s,
				actual.total_surface_outward_flow_m3_s-conservation.total_surface_outward_flow_m3_s,
				actual.open_port_outward_flow_m3_s-conservation.open_port_outward_flow_m3_s,
				actual.wall_outward_flow_m3_s-conservation.wall_outward_flow_m3_s}) physical_error = std::max(physical_error,std::abs(difference));
			for (const auto& item : conservation.surface_flow_by_boundary_label_m3_s)
				physical_error = std::max(physical_error,std::abs(actual.surface_flow_by_boundary_label_m3_s.at(item.first)-item.second));
		});
	};
	verify();
	if (mode == "flow" || mode == "inertial") {
		op.ReleaseTrial(); op.SetPortControlValue("inlet",-.0005);
		Reject(comm,[&] { op.Freeze(old,reference->Layout(),0,0,.125,1,.125); });
		op.SetPortControlValue("outlet",.0005); op.Freeze(old,reference->Layout(),0,0,.125,1,.125);
		expected_residual[op.Layout().ControllerRow(8)] += .0005; expected_residual[op.Layout().ControllerRow(9)] -= .0005;
		verify();
	}
	std::array<double,4> global{}; double global_physical = 0;
	MPI_Allreduce(errors.data(),global.data(),4,MPI_DOUBLE,MPI_SUM,comm);
	MPI_Allreduce(&physical_error,&global_physical,1,MPI_DOUBLE,MPI_MAX,comm);
	const double residual_error = std::sqrt(global[0]/global[1]),action_error = std::sqrt(global[2]/global[3]);
	iga::CollectiveLocalStage(comm,"transient operator gate",[&] {
		Require(std::isfinite(residual_error) && std::isfinite(action_error) && residual_error < 1e-9 && action_error < 1e-9
			&& std::isfinite(global_physical) && global_physical < 1e-12,"transient operator numerical gate failed");
	});
	const auto rows = op.Assembly().RowEnd()-op.Assembly().RowBegin();
	const auto halo = op.Assembly().RequiredRows().size(),work = op.Inputs().OwnedCells().size();
	op.Close(); op.Close(); Reject(comm,[&] { op.Assemble(); });
	Reject(comm,[&] { (void)op.ConservationDiagnostics(); });
	iga::CollectiveLocalStage(comm,"transient operator close state",[&] {
		Require(!op.Inputs().History().Active() && op.Options().parameters.dt == 0,"closed transient input remains active");
	});
	std::cout << "immersed_transient_operator_mpi rank=" << rank << " ranks=" << size << " owned_rows=" << rows
		<< " halo_rows=" << halo << " owned_cells=" << work << " global_rows=" << op.Layout().Rows()
		<< " residual_relative_l2=" << residual_error << " action_relative_l2=" << action_error
		<< " physical_error=" << global_physical << " passed\n";
	} catch (...) { VecDestroy(&old); throw; }
	iga::RequireCollectivePetscSuccess(comm,"transient history destroy",VecDestroy(&old));
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0,status = 0; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm comm = MPI_COMM_NULL;
	try {
		const bool split = argc > 2 && std::string(argv[2]) == "split";
		MPI_Comm_split(PETSC_COMM_WORLD,split && rank ? 1 : 0,rank,&comm);
		Run(comm,argc > 1 ? argv[1] : "flow");
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

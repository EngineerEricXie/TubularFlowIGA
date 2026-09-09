#include "ThreeDImmersedDistributedFlowDomain.hpp"
#include <iomanip>
#include <iostream>

namespace {
iga::CartesianDomainClassification Domain()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{0,0,0}},{{1,0,0}},{{1,1,0}},{{0,1,0}},{{0,0,1}},{{1,0,1}},{{1,1,1}},{{0,1,1}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground({{{0,0,0}},{{1,1,1}},{{4,1,1}}}),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	iga::CartesianDomainClassification domain{Domain()};
	iga::CutCellVolumeQuadratureCatalog volume{domain,{2,500000,500000,3000000},iga::CutCellVolumeQuadratureStorageMode::Compact};
	iga::ImmersedSurfaceQuadratureCatalog surface{domain};
	iga::CutCellGhostPenaltyCatalog ghost{domain,volume};
};
template<class Function> void Reject(MPI_Comm comm,Function&& function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	iga::CollectiveLocalStage(comm,"distributed immersed domain expected rejection",[&] { if (!rejected) throw std::runtime_error("expected collective rejection"); });
}
iga::PortBoundaryData Input(double time,double flow)
{
	iga::PortBoundaryData value; value.time_s = time; value.outward_flow_m3_s = flow; return value;
}
void Run(MPI_Comm comm)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	std::unique_ptr<Fixture> fixture;
	std::unique_ptr<iga::ImmersedStaticFlowRuntime> reference;
	std::unique_ptr<iga::ThreeDImmersedFlowDomain> serial;
	std::vector<iga::CouplingPort> ports;
	iga::ImmersedStaticFlowOptions options;
	bool inject = false; std::size_t calls = 0;
	iga::CollectiveLocalStage(comm,"distributed domain reference setup",[&] {
		fixture = std::make_unique<Fixture>();
		options.wall_labels = {7}; options.ports = {{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,0},
			{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,0}};
		options.nonlinear_absolute_tolerance = 1e-14; options.nonlinear_relative_tolerance = 1e-11; options.lu_pivot_shift = 1e-12;
		options.body_force = [&](const std::array<double,3>&) {
			if (inject && rank == size-1 && ++calls > static_cast<std::size_t>(256/size)) throw std::runtime_error("injected graph candidate failure");
			return std::array<double,3>{};
		};
		for (const auto& definition : options.ports) {
			iga::CouplingPort port; port.id = definition.id; port.subsystem_id = "immersed";
			port.locator_kind = "boundary_label"; port.locator = std::to_string(definition.boundary_label);
			port.requires = {iga::PortQuantity::FlowRate};
			port.provides = {iga::PortQuantity::Area,iga::PortQuantity::FlowRate,iga::PortQuantity::MeanPressure,iga::PortQuantity::MeanNormalTraction};
			ports.push_back(port);
		}
		auto& f = *fixture;
		reference = std::make_unique<iga::ImmersedStaticFlowRuntime>(f.domain,f.volume,f.surface,f.ghost,options);
		serial = std::make_unique<iga::ThreeDImmersedFlowDomain>("immersed",*reference,ports);
	});
	auto& f = *fixture;
	iga::ImmersedStaticDistributedRuntime runtime(comm,f.domain,f.volume,f.surface,f.ghost,options);
	for (int scenario = 0; scenario < (size > 1 ? 2 : 1); ++scenario) {
		auto invalid = ports;
		if (rank == size-1) {
			if (scenario == 0) invalid[0].orientation.native_to_outward_sign = -1;
			else invalid[0].provides.erase(iga::PortQuantity::MeanPressure);
		}
		Reject(comm,[&] { iga::ThreeDImmersedDistributedFlowDomain bad("immersed",runtime,invalid); });
	}
	iga::ThreeDImmersedDistributedFlowDomain domain("immersed",runtime,ports);
	iga::DomainStepContext first; first.start_time_s = 0; first.dt_s = 0.01; first.step_index = 0;
	for (int scenario = 0; scenario < (size > 1 ? 2 : 1); ++scenario) {
		auto bad = first;
		if (rank == size-1) { if (scenario == 0) bad.step_index = 1; else bad.dt_s *= 2; }
		Reject(comm,[&] { domain.BeginStep(bad); });
	}
	domain.BeginStep(first);
	for (int scenario = 0; scenario < (size > 1 ? 2 : 1); ++scenario) {
		auto bad = Input(first.EndTime(),-1e-6);
		if (rank == size-1) bad.outward_flow_m3_s = scenario == 0 ? std::numeric_limits<double>::quiet_NaN() : -2e-6;
		Reject(comm,[&] { domain.SetPortInput("inlet",bad); });
	}
	domain.SetPortInput("inlet",Input(first.EndTime(),-1e-6));
	Reject(comm,[&] { domain.SetPortInput("inlet",Input(first.EndTime(),-1e-6)); });
	Reject(comm,[&] { domain.SolveTrial(); });
	domain.SetPortInput("outlet",Input(first.EndTime(),1e-6));
	inject = true; calls = 0; Reject(comm,[&] { domain.SolveTrial(); }); inject = false;
	domain.SolveTrial();
	if (rank == size-1) runtime.FailNextPrepareForTesting();
	Reject(comm,[&] { domain.PrepareCommitStep(); });
	domain.RollbackTrial();
	const auto unpublished = [&] {
		iga::CollectiveLocalStage(comm,"distributed graph unpublished state",[&] {
			if (runtime.Diagnostics().commit_count || domain.Diagnostics().committed_steps || !domain.CommittedPortStates().empty()) throw std::runtime_error("graph failure published a commit");
			for (auto value : domain.CommittedOwnedBackendState()) if (value != 0.0) throw std::runtime_error("graph rollback changed committed field");
			for (const auto& control : runtime.PortDefinitions()) if (control.value != 0.0) throw std::runtime_error("graph rollback did not restore controls");
		});
	};
	unpublished();
	domain.SetPortInput("inlet",Input(first.EndTime(),-2e-6)); domain.SetPortInput("outlet",Input(first.EndTime(),2e-6));
	domain.SolveTrial(); domain.PrepareCommitStep(); domain.AbortStep(); domain.FinalizeCommitStep(); unpublished();
	double worst = 0.0;
	for (int step = 0; step < 2; ++step) {
		auto context = first; context.start_time_s = 0.01*step; context.step_index = step;
		const double flow = step ? 0.5e-6 : 1e-6;
		domain.BeginStep(context); domain.SetPortInput("inlet",Input(context.EndTime(),-flow)); domain.SetPortInput("outlet",Input(context.EndTime(),flow));
		domain.SolveTrial(); domain.PrepareCommitStep(); domain.FinalizeCommitStep();
		std::vector<PetscScalar> expected;
		iga::CollectiveLocalStage(comm,"distributed graph accepted serial reference",[&] {
			serial->BeginStep(context); serial->SetPortInput("inlet",Input(context.EndTime(),-flow)); serial->SetPortInput("outlet",Input(context.EndTime(),flow));
			serial->SolveTrial(); serial->PrepareCommitStep(); serial->FinalizeCommitStep(); expected = serial->CommittedBackendState();
		});
		double local[6]{},global[6]{};
		iga::CollectiveLocalStage(comm,"distributed graph field comparison",[&] {
			const auto actual = domain.CommittedOwnedBackendState();
			for (PetscInt row = domain.RowBegin(); row < domain.RowEnd(); ++row) {
				const int field = static_cast<std::size_t>(row) >= runtime.Diagnostics().physical_dofs ? 2 : row%4 == 3 ? 1 : 0;
				const double value = actual[row-domain.RowBegin()],ref = expected[row];
				local[2*field] += (value-ref)*(value-ref); local[2*field+1] += ref*ref;
			}
			const auto close = [](double a,double b) { if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b) > 1e-12+1e-6*std::abs(b)) throw std::runtime_error("distributed graph port differs from serial"); };
			for (const auto& entry : domain.CommittedPortStates()) {
				const auto& a = entry.second; const auto& b = serial->CommittedPortStates().at(entry.first);
				close(a.time_s,b.time_s); close(*a.area_m2,*b.area_m2); close(*a.outward_flow_m3_s,*b.outward_flow_m3_s);
				close(*a.mean_pressure_pa,*b.mean_pressure_pa); close(*a.mean_normal_traction_pa,*b.mean_normal_traction_pa);
			}
			if (domain.Diagnostics().committed_time_s != context.EndTime() || domain.Diagnostics().committed_step_index != step
				|| runtime.Diagnostics().commit_count != static_cast<std::size_t>(step+1)) throw std::runtime_error("distributed graph commit clock differs");
		});
		MPI_Allreduce(local,global,6,MPI_DOUBLE,MPI_SUM,comm);
		iga::CollectiveLocalStage(comm,"distributed graph per-field gate",[&] {
			for (int field = 0; field < 3; ++field) {
				const double relative = std::sqrt(global[2*field]/global[2*field+1]); worst = std::max(worst,relative);
				if (!(relative <= 1e-6)) throw std::runtime_error("distributed graph field exceeds relative L2 gate");
			}
		});
	}
	const auto conservation = runtime.ConservationDiagnostics();
	iga::CollectiveLocalStage(comm,"distributed graph conservation",[&] {
		const auto expected = reference->ConservationDiagnostics();
		if (std::abs(conservation.total_surface_outward_flow_m3_s-expected.total_surface_outward_flow_m3_s) > 1e-12
			|| std::abs(conservation.volume_divergence_integral_m3_s-expected.volume_divergence_integral_m3_s) > 1e-12)
			throw std::runtime_error("distributed graph conservation differs");
	});
	Reject(comm,[&] { domain.AbortStep(); });
	auto closing_step = first; closing_step.start_time_s = 0.02; closing_step.step_index = 2;
	domain.BeginStep(closing_step); domain.SetPortInput("inlet",Input(closing_step.EndTime(),-0.5e-6)); domain.SetPortInput("outlet",Input(closing_step.EndTime(),0.5e-6));
	domain.SolveTrial(); domain.PrepareCommitStep(); runtime.Close(); domain.FinalizeCommitStep();
	iga::CollectiveLocalStage(comm,"closed backend does not publish graph metadata",[&] {
		if (domain.Diagnostics().committed_steps != 2 || domain.Diagnostics().committed_time_s != 0.02
			|| runtime.Diagnostics().commit_count != 2) throw std::runtime_error("closed graph backend published a prepared step");
	});
	serial.reset(); reference.reset();
	std::cout << std::setprecision(17) << "immersed_domain_mpi rank=" << rank << " ranks=" << size
		<< " cells=" << f.domain.Cells().size() << " owned_rows=" << domain.RowEnd()-domain.RowBegin()
		<< " steps=" << domain.Diagnostics().committed_steps << " field_relative_l2=" << worst << " passed\n";
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int status = 0; MPI_Comm group = MPI_COMM_NULL;
	try {
		if (argc > 1 && std::string(argv[1]) == "split") {
			int rank = 0,size = 1; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&size);
			if (size != 3) throw std::invalid_argument("split domain test requires three ranks");
			MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group); Run(group);
		} else Run(PETSC_COMM_WORLD);
	} catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	if (group != MPI_COMM_NULL) MPI_Comm_free(&group);
	PetscFinalize(); return status;
}

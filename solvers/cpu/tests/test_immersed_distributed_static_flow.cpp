#include "ImmersedStaticDistributedRuntime.hpp"
#include "ImmersedStaticFlowRuntime.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
iga::CartesianDomainClassification Domain(bool padded,bool empty_work)
{
	iga::RawSurfaceSoup soup;
	const double a = empty_work ? 0.0 : padded ? 0.3 : 0.1, b = empty_work ? 1.0 : padded ? 0.7 : 0.9;
	soup.vertices = {{{{a,a,a}},{{b,a,a}},{{b,b,a}},{{a,b,a}},{{a,a,b}},{{b,a,b}},{{b,b,b}},{{a,b,b}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	const unsigned cells = padded ? 5 : 3;
	const auto shape = empty_work ? std::array<std::uint32_t,3>{{2,1,1}} : std::array<std::uint32_t,3>{{cells,cells,cells}};
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},shape};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(grid),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	Fixture(bool padded,bool expanded,bool empty_work)
		: domain(Domain(padded,empty_work)),volume(domain,{2,500000,500000,3000000},expanded ? iga::CutCellVolumeQuadratureStorageMode::Expanded : iga::CutCellVolumeQuadratureStorageMode::Compact) {}
	iga::CartesianDomainClassification domain;
	iga::CutCellVolumeQuadratureCatalog volume;
	iga::ImmersedSurfaceQuadratureCatalog surface{domain};
	iga::CutCellGhostPenaltyCatalog ghost{domain, volume};
};

template<class Function>
void Reject(MPI_Comm comm,Function&& function)
{
	bool rejected = false;
	try { function(); } catch (const std::exception&) { rejected = true; }
	iga::CollectiveLocalStage(comm,"static MPI expected rejection",[&] { if (!rejected) throw std::runtime_error("expected collective rejection"); });
}
void Run(MPI_Comm comm,const std::string& mode)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	if (mode != "closed" && mode != "flow" && mode != "pressure" && mode != "empty-work") throw std::invalid_argument("unknown static MPI mode");
	const bool closed = mode == "closed" || mode == "empty-work";
	int failure_rank = size-1;
	std::unique_ptr<Fixture> fixture;
	std::unique_ptr<iga::ImmersedStaticFlowRuntime> serial;
	iga::ImmersedStaticFlowOptions options;
	std::size_t calls = 0,throw_on = 0;
	std::vector<PetscScalar> reference;
	iga::ImmersedStaticFlowConservationDiagnostics reference_conservation;
	iga::CollectiveLocalStage(comm,"static MPI serial reference",[&] {
		fixture = std::make_unique<Fixture>(false,false,mode == "empty-work");
		options.wall_labels = {7,8,9}; options.lu_pivot_shift = 1e-12;
		// Both reference and distributed paths use the same tighter tolerance:
		// tiny flow fields must converge beyond the first Newton correction.
		options.nonlinear_absolute_tolerance = 1e-14;
		options.nonlinear_relative_tolerance = 1e-11;
		if (!closed) {
			options.wall_labels = {7};
			options.ports = {{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-1e-6},
				{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,1e-6}};
			if (mode == "pressure") options.ports[1] = {"outlet",9,iga::ImmersedFlowPortControlMode::Pressure,1e-3};
		}
		options.body_force = [&](const std::array<double,3>&) {
			++calls;
			if (rank == failure_rank && throw_on && calls >= throw_on) throw std::runtime_error("injected static candidate failure");
			return closed ? std::array<double,3>{{.3,-.2,.1}} : std::array<double,3>{};
		};
		serial = std::make_unique<iga::ImmersedStaticFlowRuntime>(fixture->domain,fixture->volume,fixture->surface,fixture->ghost,options);
		if (!serial->SolveTrial()) throw std::runtime_error("serial static reference did not converge");
		serial->Commit(); reference = serial->CommittedState();
		reference_conservation = serial->ConservationDiagnostics();
	});
	auto& f = *fixture;
	iga::ImmersedStaticDistributedRuntime runtime(comm,f.domain,f.volume,f.surface,f.ghost,options);
	std::vector<PetscScalar> zero(static_cast<std::size_t>(runtime.RowEnd()-runtime.RowBegin()),0.0);
	runtime.SetCommittedOwnedState(zero);
	calls = 0; runtime.Assemble(); const auto calls_per_assembly = calls;
	const int candidate_rank = calls_per_assembly ? rank : -1;
	MPI_Allreduce(&candidate_rank,&failure_rank,1,MPI_INT,MPI_MAX,comm);
	iga::CollectiveLocalStage(comm,"static MPI integration coverage",[&] {
		if (failure_rank < 0 || (mode != "empty-work" && !calls_per_assembly)) throw std::runtime_error("test rank has no volume work");
		if (mode == "empty-work" && size == 4 && (rank == 1 || rank == 3)
			&& (calls_per_assembly || runtime.OwnedStencilCount() || runtime.RequiredStateRows())) throw std::runtime_error("empty-work rank integrated or retained a halo");
	});
	calls = 0; throw_on = calls_per_assembly+1;
	Reject(comm,[&] { runtime.SolveTrial(); }); throw_on = 0;
	iga::CollectiveLocalStage(comm,"static MPI rollback state",[&] {
		const auto& d = runtime.Diagnostics();
		if (d.trial_active || d.prepared || d.converged || !d.committed || d.rollback_count != 1) throw std::runtime_error("static rollback flags differ");
		iga::PetscReadArray state; state.Acquire(runtime.State());
		for (PetscInt row = runtime.RowBegin(); row < runtime.RowEnd(); ++row)
			if (state.Data()[row-runtime.RowBegin()] != 0.0) throw std::runtime_error("static candidate rollback changed committed values");
		state.Restore();
	});
	if (!runtime.SolveTrial()) throw std::runtime_error("distributed static trial did not converge");
	const auto accepted = runtime.Diagnostics();
	if (rank == size-1) runtime.FailNextPrepareForTesting();
	Reject(comm,[&] { runtime.PrepareCommit(); });
	iga::CollectiveLocalStage(comm,"static MPI unpublished preparation",[&] {
		if (runtime.Diagnostics().prepared || !runtime.Diagnostics().trial_active) throw std::runtime_error("failed prepare published flags");
		iga::PetscReadArray committed; committed.Acquire(runtime.CommittedState());
		for (PetscInt row = runtime.RowBegin(); row < runtime.RowEnd(); ++row)
			if (committed.Data()[row-runtime.RowBegin()] != 0.0) throw std::runtime_error("failed prepare changed committed values");
		committed.Restore();
	});
	runtime.PrepareCommit(); runtime.FinalizeCommit();
	double local[6]{},global[6]{};
	iga::CollectiveLocalStage(comm,"static MPI solution parity",[&] {
		const auto& d = runtime.Diagnostics();
		if (!d.committed || d.trial_active || d.prepared || !d.converged || d.commit_count != 1 || d.ksp_reason <= 0 || d.newton_steps.empty())
			throw std::runtime_error("static committed solve diagnostics invalid");
		for (const auto& port : d.ports) if (port.control_mode == iga::ImmersedFlowPortControlMode::FlowRate)
			if (!(port.normalized_flow_residual <= 1.0)) throw std::runtime_error("static MPI flow controller did not converge");
		for (const auto& step : d.newton_steps) {
			std::cout << std::setprecision(17) << "immersed_static_linear rank=" << rank << " iteration=" << step.iteration
				<< " rhs=" << step.residual_norm << " absolute=" << step.linear_residual_norm << " relative=" << step.linear_relative_residual << '\n';
			if (!std::isfinite(step.linear_relative_residual) || step.linear_relative_residual > 1e-6)
				throw std::runtime_error("static true linear residual exceeds gate");
		}
		iga::PetscReadArray committed; committed.Acquire(runtime.CommittedState());
		for (PetscInt row = runtime.RowBegin(); row < runtime.RowEnd(); ++row) {
			const double value = PetscRealPart(committed.Data()[row-runtime.RowBegin()]),expected = PetscRealPart(reference[row]);
			const std::size_t field = static_cast<std::size_t>(row) >= d.physical_dofs ? 2 : row%4 == 3 ? 1 : 0;
			local[2*field] += (value-expected)*(value-expected); local[2*field+1] += expected*expected;
		}
		committed.Restore();
	});
	MPI_Allreduce(local,global,6,MPI_DOUBLE,MPI_SUM,comm);
	const double relative = std::sqrt((global[0]+global[2]+global[4])/(global[1]+global[3]+global[5]));
	const double velocity_error = std::sqrt(global[0]), pressure_relative = std::sqrt(global[2]/global[3]);
	const double velocity_relative = global[1] > 0.0 ? std::sqrt(global[0]/global[1]) : 0.0;
	const double scalar_error = std::sqrt(global[4]), scalar_relative = global[5] > 0.0 ? std::sqrt(global[4]/global[5]) : 0.0;
	std::cout << std::setprecision(17) << "immersed_static_fields rank=" << rank
		<< " velocity_absolute=" << velocity_error << " velocity_relative=" << velocity_relative
		<< " pressure_relative=" << pressure_relative << " scalar_absolute=" << scalar_error << " scalar_relative=" << scalar_relative << '\n';
	iga::CollectiveLocalStage(comm,"static MPI per-field gates",[&] {
		if (!(pressure_relative <= 1e-6) || !(closed ? velocity_error <= 1e-10 && scalar_error <= 1e-10
			: velocity_relative <= 1e-6 && scalar_relative <= 1e-6)) throw std::runtime_error("static MPI per-field error exceeds gate");
	});
	iga::CollectiveLocalStage(comm,"static MPI field gate",[&] { if (!(relative <= 1e-6)) throw std::runtime_error("static MPI field differs from serial"); });
	const auto conservation = runtime.ConservationDiagnostics();
	iga::CollectiveLocalStage(comm,"static MPI conservation parity",[&] {
		const auto close = [](double a,double b) {
			if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a-b) > 1e-12+1e-6*std::abs(b))
				throw std::runtime_error("distributed static conservation differs from serial");
		};
		close(conservation.volume_divergence_integral_m3_s,reference_conservation.volume_divergence_integral_m3_s);
		close(conservation.total_surface_outward_flow_m3_s,reference_conservation.total_surface_outward_flow_m3_s);
		close(conservation.wall_outward_flow_m3_s,reference_conservation.wall_outward_flow_m3_s);
		close(conservation.open_port_outward_flow_m3_s,reference_conservation.open_port_outward_flow_m3_s);
		if (conservation.surface_flow_by_boundary_label_m3_s.size() != reference_conservation.surface_flow_by_boundary_label_m3_s.size())
			throw std::runtime_error("distributed conservation label set differs");
		for (const auto& entry : reference_conservation.surface_flow_by_boundary_label_m3_s)
			close(conservation.surface_flow_by_boundary_label_m3_s.at(entry.first),entry.second);
	});
	std::vector<PetscScalar> committed_before;
	iga::CollectiveLocalStage(comm,"static MPI committed snapshot",[&] {
		iga::PetscReadArray view; view.Acquire(runtime.CommittedState());
		committed_before.assign(view.Data(),view.Data()+zero.size()); view.Restore();
	});
	for (int invalid = 0; invalid < 2; ++invalid) {
		auto bad = zero;
		if (rank == size-1) {
			if (invalid == 0) bad.pop_back();
			else bad[0] = std::numeric_limits<double>::quiet_NaN();
		}
		Reject(comm,[&] { runtime.SetCommittedOwnedState(bad); });
		iga::CollectiveLocalStage(comm,"static MPI rejected input preserves state",[&] {
			iga::PetscReadArray committed,trial; committed.Acquire(runtime.CommittedState()); trial.Acquire(runtime.State());
			for (std::size_t i = 0; i < committed_before.size(); ++i)
				if (committed.Data()[i] != committed_before[i] || trial.Data()[i] != committed_before[i]) throw std::runtime_error("invalid owned input changed state");
			committed.Restore(); trial.Restore();
		});
	}
	runtime.SetCommittedOwnedState(zero);
	const auto cleared_conservation = runtime.ConservationDiagnostics();
	iga::CollectiveLocalStage(comm,"static MPI conservation refreshes halo",[&] {
		if (cleared_conservation.volume_divergence_integral_m3_s != 0.0 || cleared_conservation.total_surface_outward_flow_m3_s != 0.0)
			throw std::runtime_error("conservation reused a stale state halo");
		for (const auto& entry : cleared_conservation.surface_flow_by_boundary_label_m3_s)
			if (entry.second != 0.0) throw std::runtime_error("boundary conservation reused a stale state halo");
	});
	runtime.SetCommittedOwnedState(committed_before);
	Reject(comm,[&] { runtime.Commit(); });
	if (!runtime.SolveTrial()) throw std::runtime_error("committed static root no longer converges");
	runtime.PrepareCommit(); runtime.Rollback(); runtime.FinalizeCommit();
	iga::CollectiveLocalStage(comm,"static MPI prepared rollback",[&] {
		if (runtime.Diagnostics().commit_count != 1 || runtime.Diagnostics().prepared || runtime.Diagnostics().trial_active || runtime.Diagnostics().rollback_count != 2)
			throw std::runtime_error("prepared static rollback changed publication");
	});
	if (!runtime.SolveTrial()) throw std::runtime_error("static root failed before close");
	runtime.PrepareCommit();
	const auto commits_before_close = runtime.Diagnostics().commit_count;
	runtime.Close(); runtime.FinalizeCommit(); runtime.Close();
	iga::CollectiveLocalStage(comm,"static MPI close discards preparation",[&] {
		if (runtime.Diagnostics().prepared || runtime.Diagnostics().trial_active || runtime.Diagnostics().commit_count != commits_before_close)
			throw std::runtime_error("closed static runtime published a prepared state");
	});
	Reject(comm,[&] { runtime.Rollback(); });
	serial.reset();
	std::cout << std::setprecision(17) << "immersed_static_mpi rank=" << rank << " ranks=" << size
		<< " owned_rows=" << runtime.RowEnd()-runtime.RowBegin() << " owned_stencils=" << runtime.OwnedStencilCount() << " halo_rows=" << runtime.RequiredStateRows() << " field_relative_l2=" << relative
		<< " velocity_absolute_l2=" << velocity_error << " velocity_relative_l2=" << velocity_relative
		<< " volume_divergence=" << conservation.volume_divergence_integral_m3_s << " surface_flux=" << conservation.total_surface_outward_flow_m3_s
		<< " pressure_relative_l2=" << pressure_relative << " scalar_absolute_l2=" << scalar_error << " scalar_relative_l2=" << scalar_relative
		<< " nonlinear_iterations=" << accepted.nonlinear_iterations << " ksp_iterations=" << accepted.ksp_iterations
		<< " residual=" << accepted.residual_norm << " assembly_s=" << accepted.aggregate_assembly_seconds
		<< " linear_s=" << accepted.aggregate_linear_solve_seconds << " passed\n";

}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int status = 0;
	MPI_Comm group = MPI_COMM_NULL;
	try {
		const std::string mode = argc > 1 ? argv[1] : "closed";
		if (argc > 2 && std::string(argv[2]) == "split") {
			int rank = 0,size = 1; MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&size);
			if (size != 3) throw std::invalid_argument("static split test requires three ranks");
			MPI_Comm_split(PETSC_COMM_WORLD,rank == 0 ? 0 : 1,rank,&group);
			Run(group,mode);
		} else Run(PETSC_COMM_WORLD,mode);
	}
	catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	if (group != MPI_COMM_NULL) MPI_Comm_free(&group);
	PetscFinalize(); return status;
}

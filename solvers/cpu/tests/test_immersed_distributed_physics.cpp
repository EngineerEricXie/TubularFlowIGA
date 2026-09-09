#include "ImmersedDistributedAssembly.hpp"
#include "ImmersedStaticFlowRuntime.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>

namespace {
iga::CartesianDomainClassification Domain()
{
	iga::RawSurfaceSoup soup;
	const double a = 0.1, b = 0.9;
	soup.vertices = {{{{a,a,a}},{{b,a,a}},{{b,b,a}},{{a,b,a}},{{a,a,b}},{{b,a,b}},{{b,b,b}},{{a,b,b}}}};
	const std::array<std::array<std::int64_t,3>,12> triangles{{{{0,2,1}},{{0,3,2}},{{4,5,6}},{{4,6,7}},
		{{0,1,5}},{{0,5,4}},{{1,2,6}},{{1,6,5}},{{2,3,7}},{{2,7,6}},{{3,0,4}},{{3,4,7}}}};
	for (std::size_t i = 0; i < triangles.size(); ++i) {
		iga::RawSurfaceTriangle triangle; triangle.indices = triangles[i]; triangle.boundary_id = i < 2 ? 8 : i < 4 ? 9 : 7;
		soup.triangles.push_back(triangle);
	}
	const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(grid),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	iga::CartesianDomainClassification domain{Domain()};
	iga::CutCellVolumeQuadratureCatalog volume{domain, {2,500000,500000,3000000}, iga::CutCellVolumeQuadratureStorageMode::Compact};
	iga::ImmersedSurfaceQuadratureCatalog surface{domain};
	iga::CutCellGhostPenaltyCatalog ghost{domain, volume};
};
struct Work { int kind; std::uint64_t cell; std::size_t port; };

void Run(MPI_Comm comm)
{
	int rank = 0, size = 1;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &size);
	std::unique_ptr<Fixture> fixture;
	std::unique_ptr<iga::ImmersedStaticFlowRuntime> serial;
	iga::ImmersedStaticFlowOptions options;
	std::vector<PetscScalar> state, expected_residual, expected_action;
	std::vector<PetscInt> offsets;
	std::vector<iga::ImmersedAssemblyStencil> stencils;
	std::vector<Work> work;
	// Each validation process builds its own COMM_SELF reference. This is
	// reference-only; the distributed numerical callback below owns each cell
	// and face once and never calls the reference runtime to obtain a block.
	std::exception_ptr setup_error;
	try {
		fixture = std::make_unique<Fixture>();
		options.wall_labels = {7};
		options.ports = {{"inlet",8,iga::ImmersedFlowPortControlMode::FlowRate,-0.001},
			{"outlet",9,iga::ImmersedFlowPortControlMode::FlowRate,0.001}};
		options.body_force = [](const std::array<double,3>& x) { return std::array<double,3>{{0.1+x[0],-0.2*x[1],0.3*x[2]}}; };
		serial = std::make_unique<iga::ImmersedStaticFlowRuntime>(fixture->domain, fixture->volume, fixture->surface, fixture->ghost, options);
		state.resize(serial->Diagnostics().total_dofs);
		for (std::size_t i = 0; i < state.size(); ++i) state[i] = 0.01*std::sin(0.13*(i+1));
		serial->SetCommittedState(state); serial->Assemble();
		expected_residual = serial->AssembledNegativeResidual(); expected_action = serial->AssembledJacobianAction(state);
	} catch (...) { setup_error = std::current_exception(); }
	iga::CollectiveLocalStage(comm, "immersed physics reference setup", [&] { if (setup_error) std::rethrow_exception(setup_error); });
	auto& f = *fixture;
	iga::CollectiveLocalStage(comm, "immersed physics distributed topology", [&] {
		const auto physical = static_cast<PetscInt>(serial->Diagnostics().physical_dofs);
		for (int r = 0; r <= size; ++r) offsets.push_back(r == size ? static_cast<PetscInt>(state.size()) : 4*((physical/4)*r/size));
		const auto add = [&](int owner, iga::ImmersedStencilPattern pattern, std::vector<PetscInt> rows, Work task) {
			stencils.push_back({stencils.size(),owner,pattern,std::move(rows)}); work.push_back(task);
		};
		for (std::uint64_t cell = 0; cell < f.domain.Cells().size(); ++cell) {
			if (!f.volume.Cell(cell).usable || !iga::CompactCutCellVolumeLogicalPointCount(f.volume.Cell(cell).compact_rule)) continue;
			const int owner = static_cast<int>(cell*size/f.domain.Cells().size());
			const auto element = f.domain.Background().MaterializeElement(cell);
			std::vector<PetscInt> rows, pressure, velocity;
			for (auto node : element.connectivity) {
				for (int field = 0; field < 4; ++field) rows.push_back(serial->Dof(node,field));
				pressure.push_back(serial->Dof(node,3));
				for (int field = 0; field < 3; ++field) velocity.push_back(serial->Dof(node,field));
			}
			add(owner,iga::ImmersedStencilPattern::Dense,rows,{0,cell,0});
			pressure.push_back(serial->GaugeDof()); add(owner,iga::ImmersedStencilPattern::Scalar,pressure,{2,cell,0});
			const auto& rule = f.surface.UsableRule(f.domain,cell);
			for (std::size_t port = 0; port < options.ports.size(); ++port) {
				if (!std::any_of(rule.Points().begin(),rule.Points().end(),[&](const auto& point) { return point.boundary_id == options.ports[port].boundary_label; })) continue;
				auto port_rows = velocity; port_rows.push_back(serial->PortMultiplierDof(options.ports[port].id));
				add(owner,iga::ImmersedStencilPattern::Scalar,std::move(port_rows),{3,cell,port});
			}
		}
		for (std::size_t face = 0; face < f.ghost.Faces().size(); ++face) {
			const auto& pair = f.ghost.Faces()[face];
			const auto nodes = iga::CubicCartesianSplineFaceConnectivity(f.domain,pair.minus_cell,pair.plus_cell);
			std::vector<PetscInt> rows;
			for (auto node : nodes) for (int field = 0; field < 4; ++field) rows.push_back(serial->Dof(node,field));
			add(static_cast<int>(pair.minus_cell*size/f.domain.Cells().size()),iga::ImmersedStencilPattern::SameField,std::move(rows),{1,face,0});
		}
		for (std::size_t port = 0; port < options.ports.size(); ++port)
			add(size-1,iga::ImmersedStencilPattern::Scalar,{serial->PortMultiplierDof(options.ports[port].id)},{4,0,port});
	});
	iga::ImmersedDistributedAssembly distributed(comm,offsets,stencils);
	iga::CollectiveLocalStage(comm, "immersed physics state insertion", [&] {
		for (PetscInt row = distributed.RowBegin(); row < distributed.RowEnd(); ++row)
			if (VecSetValue(distributed.State(),row,state[static_cast<std::size_t>(row)],INSERT_VALUES)) throw std::runtime_error("cannot insert physical state");
	});
	iga::RequireCollectivePetscSuccess(comm,"physics state begin",VecAssemblyBegin(distributed.State()));
	iga::RequireCollectivePetscSuccess(comm,"physics state end",VecAssemblyEnd(distributed.State()));
	std::vector<int> seen(stencils.size(),0);
	const auto integrate = [&](const iga::ImmersedAssemblyStencil& stencil, auto& matrix, auto& residual) {
		++seen[stencil.id]; const auto& task = work[stencil.id]; const auto n = stencil.rows.size();
		if (task.kind == 1) {
			auto block = f.ghost.AssembleFaceLocal(task.cell,f.domain,f.volume,[&](std::int32_t node,int field) {
				return distributed.StateAt(serial->Dof(node,field));
			},options.parameters.dynamic_viscosity);
			matrix = std::move(block.jacobian); residual = std::move(block.negative_residual); return;
		}
		matrix.assign(n*n,0.0); residual.assign(n,0.0);
		if (task.kind == 4) { residual[0] = options.ports[task.port].value; return; }
		const auto element = f.domain.Background().MaterializeElement(task.cell);
		const auto& rule = f.surface.UsableRule(f.domain,task.cell);
		std::vector<std::array<double,4>> nodal(element.connectivity.size());
		for (std::size_t a = 0; a < nodal.size(); ++a)
			for (int field = 0; field < 4; ++field) nodal[a][field] = distributed.StateAt(serial->Dof(element.connectivity[a],field));
		const auto points = [&](const auto& consume) { iga::ForEachVolumePoint(f.volume.UsableCompactRule(f.domain,task.cell),consume); };
		if (task.kind == 0) {
			auto volume = iga::BuildNavierStokesElementFromPoints(element,nodal,{},options.parameters,points,options.body_force,iga::NavierStokesResolvedMixedForm::Conservative);
			matrix = volume.jacobian; residual = volume.negative_residual;
			if (f.domain.Cells()[task.cell].classification == iga::CellClassification::Cut) {
				const auto trace = iga::BuildImmersedConservativeMixedTraceElement(element,rule,nodal);
				for (std::size_t i = 0; i < matrix.size(); ++i) matrix[i] += trace.jacobian[i];
				for (std::size_t i = 0; i < residual.size(); ++i) residual[i] += trace.negative_residual[i];
				if (std::any_of(rule.Points().begin(),rule.Points().end(),[](const auto& point) { return point.boundary_id == 7; })) {
					const auto wall = iga::BuildImmersedNitscheWallElementFromVolumeSystem(f.domain,f.volume,f.surface,task.cell,nodal,{},options.parameters,
						options.wall_labels,volume,f.ghost,options.wall_gamma0,options.wall_velocity);
					for (std::size_t i = 0; i < matrix.size(); ++i) matrix[i] += wall.system.jacobian[i]-volume.jacobian[i];
					for (std::size_t i = 0; i < residual.size(); ++i) residual[i] += wall.system.negative_residual[i]-volume.negative_residual[i];
				}
			}
		} else {
			std::vector<double> coefficients(n-1,0.0);
			if (task.kind == 2) {
				points([&](const iga::VolumeQuadraturePoint& point) {
					const auto basis = iga::EvaluateBasis(element,point.parametric[0],point.parametric[1],point.parametric[2],false);
					for (std::size_t a = 0; a < coefficients.size(); ++a) coefficients[a] += basis.value[a]*point.weight*basis.raw_determinant;
				});
			} else {
				const auto& port = options.ports[task.port];
				const auto block = iga::BuildImmersedFlowPortElement(element,rule,port.boundary_label,port.control_mode,port.value,nodal);
				for (std::size_t a = 0; a < element.connectivity.size(); ++a)
					for (int field = 0; field < 3; ++field) coefficients[3*a+field] = block.flow_coefficient[4*a+field];
			}
			const double lambda = distributed.StateAt(stencil.rows.back());
			for (std::size_t a = 0; a < coefficients.size(); ++a) {
				matrix[a*n+n-1] = matrix[(n-1)*n+a] = coefficients[a];
				residual[a] = -coefficients[a]*lambda;
				residual.back() -= coefficients[a]*distributed.StateAt(stencil.rows[a]);
			}
		}
	};
	const auto start = std::chrono::steady_clock::now();
	distributed.Assemble(integrate);
	const double assembly_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
	std::vector<int> global_seen(seen.size());
	MPI_Allreduce(seen.data(),global_seen.data(),static_cast<int>(seen.size()),MPI_INT,MPI_SUM,comm);
	iga::CollectiveLocalStage(comm,"immersed physics unique work",[&] {
		if (!std::all_of(global_seen.begin(),global_seen.end(),[](int count) { return count == 1; })) throw std::runtime_error("physical stencil omitted or repeated");
	});
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
		<< " halo_rows=" << distributed.RequiredRows().size() << " global_rows=" << distributed.Rows()
		<< " local_nz=" << info.nz_used << " assembly_s=" << assembly_seconds
		<< " residual_relative_l2=" << residual_error << " action_relative_l2=" << action_error << " passed\n";
	distributed.Close(); serial.reset();
}
}
int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int status = 0;
	try { Run(PETSC_COMM_WORLD); }
	catch (const std::exception& error) { std::cerr << error.what() << '\n'; status = 1; }
	PetscFinalize(); return status;
}

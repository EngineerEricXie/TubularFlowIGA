#include "ImmersedDistributedTransientVolume.hpp"
#include "ImmersedDistributedAssembly.hpp"
#include <iostream>
#include <memory>

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
	for (auto& vertex : soup.vertices) for (auto& x : vertex) x = .125+.75*x;
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground({{{0,0,0}},{{1,1,1}},{{4,1,1}}}),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	iga::CartesianDomainClassification domain{Domain()};
	iga::CutCellVolumeQuadratureCatalog expanded{domain,{2,500000,500000,3000000}};
	iga::CutCellVolumeQuadratureCatalog compact{domain,{2,500000,500000,3000000},iga::CutCellVolumeQuadratureStorageMode::Compact};
};
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
void Run(MPI_Comm comm, int mode, bool compact)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	std::unique_ptr<Fixture> fixture;
	iga::ImmersedActiveLayout layout;
	std::vector<std::uint64_t> cells;
	std::vector<PetscInt> offsets;
	std::vector<iga::ImmersedAssemblyStencil> stencils;
	iga::CollectiveLocalStage(comm,"transient volume fixture",[&] {
		fixture = std::make_unique<Fixture>();
		const auto& volume = compact ? fixture->compact : fixture->expanded;
		layout = iga::ImmersedActiveLayout::Build(fixture->domain,volume,"fixed",{8,9},true);
		const auto nodes = static_cast<PetscInt>(layout.NodeIds().size()),rows = static_cast<PetscInt>(layout.Rows());
		for (int r = 0; r <= size; ++r)
			offsets.push_back(mode == 1 ? (r ? rows : 0) : mode == 2 ? (r == size ? rows : 0) : r == size ? rows : 4*(nodes*r/size));
		for (std::uint64_t c = 0; c < 4; ++c) {
			iga::ImmersedAssemblyStencil stencil; stencil.id = c; stencil.owner = mode ? 0 : static_cast<int>(c)*size/4;
			for (auto id : fixture->domain.Background().MaterializeElement(c).connectivity)
				for (int f = 0; f < 4; ++f) stencil.rows.push_back(static_cast<PetscInt>(4*layout.LocalNode(id)+f));
			stencils.push_back(stencil);
			if (stencil.owner == rank) cells.push_back(c);
		}
	});
	auto& f = *fixture; const auto& volume = compact ? f.compact : f.expanded;
	iga::ImmersedDistributedAssembly assembly(comm,offsets,stencils);
	const auto fill = [&](double scale) {
		iga::CollectiveLocalStage(comm,"transient volume state",[&] {
			for (PetscInt row = assembly.RowBegin(); row < assembly.RowEnd(); ++row) {
				const double value = static_cast<std::size_t>(row) < layout.NodeFieldRows() ? Value(layout.NodeIds()[row/4],row%4,scale) : 100+row;
				Require(!VecSetValue(assembly.State(),row,value,INSERT_VALUES),"volume state insertion failed");
			}
		});
		iga::RequireCollectivePetscSuccess(comm,"volume state begin",VecAssemblyBegin(assembly.State()));
		iga::RequireCollectivePetscSuccess(comm,"volume state end",VecAssemblyEnd(assembly.State()));
	};
	fill(1);
	// Missing ownership and duplicate ownership must fail before halo creation.
	auto missing = cells;
	missing.erase(std::remove(missing.begin(),missing.end(),0),missing.end());
	Reject(comm,[&] { iga::ImmersedDistributedTransientVolume invalid(comm,f.domain,volume,layout,assembly.State(),missing); });
	if (size > 1) {
		auto duplicate = cells; if (rank == size-1) duplicate.push_back(0);
		Reject(comm,[&] { iga::ImmersedDistributedTransientVolume invalid(comm,f.domain,volume,layout,assembly.State(),duplicate); });
		Reject(comm,[&] { iga::ImmersedDistributedTransientVolume invalid(comm,f.domain,rank == size-1 ? f.expanded : f.compact,layout,assembly.State(),cells); });
	}
	iga::ImmersedDistributedTransientVolume inputs(comm,f.domain,volume,layout,assembly.State(),cells);
	std::size_t calls = 0; bool inject = false,forbid = false; double force_scale = 1;
	const auto field = [](const std::array<double,3>& x) { return std::array<double,3>{{.1+x[0],-.2*x[1],.3*x[2]}}; };
	const iga::NavierStokesBodyForceEvaluator force = [&](const std::array<double,3>& x) {
		++calls;
		if (forbid || (inject && rank == 0 && calls == 2)) throw std::runtime_error("injected owned force callback failure");
		auto value = field(x); for (double& v : value) v *= force_scale; return value;
	};
	const iga::NavierStokesParameters parameters{1,.1,.125};
	inject = true;
	Reject(comm,[&] { inputs.Freeze(assembly.State(),layout,0,0,.125,1,parameters,force); });
	inject = false; calls = 0;
	iga::CollectiveLocalStage(comm,"volume force failure state",[&] {
		Require(!inputs.History().Active() && inputs.FrozenPointCount() == 0,"failed volume inputs published");
	});
	const iga::NavierStokesBodyForceEvaluator nonfinite = [&](const std::array<double,3>& x) {
		auto value = field(x); if (rank == 0) value[0] = std::numeric_limits<double>::quiet_NaN(); return value;
	};
	Reject(comm,[&] { inputs.Freeze(assembly.State(),layout,0,0,.125,1,parameters,nonfinite); });
	if (size > 1) {
		auto different = parameters; if (rank == size-1) different.dynamic_viscosity *= 2;
		Reject(comm,[&] { inputs.Freeze(assembly.State(),layout,0,0,.125,1,different,force); });
	}
	inputs.Freeze(assembly.State(),layout,0,0,.125,1,parameters,force);
	iga::CollectiveLocalStage(comm,"volume force single evaluation",[&] { Require(calls == inputs.FrozenPointCount(),"force was not called once per owned point"); });
	const auto frozen_points = inputs.FrozenPointCount();
	forbid = true; force_scale = 99; fill(.75);
	const auto integrate = [&](const auto& stencil,auto& matrix,auto& residual) {
		const auto element = f.domain.Background().MaterializeElement(stencil.id);
		std::vector<std::array<double,4>> current(element.connectivity.size());
		for (std::size_t a = 0; a < current.size(); ++a)
			for (int c = 0; c < 4; ++c) current[a][c] = assembly.StateAt(stencil.rows[4*a+c]);
		auto block = inputs.BuildVolume(stencil.id,current);
		matrix = std::move(block.jacobian); residual = std::move(block.negative_residual);
	};
	double worst_error = 0,history_effect = 0;
	// The failed-assembly retry has identical inputs. Cache its independent
	// serial baseline; still compare every owned residual/action on every retry.
	std::map<std::pair<double,double>,std::pair<std::vector<double>,std::vector<double>>> references;
	const auto verify = [&](double old_scale,double expected_force_scale) {
		assembly.Assemble(integrate);
		Vec action = nullptr;
		iga::RequireCollectivePetscSuccess(comm,"volume action create",VecDuplicate(assembly.State(),&action));
		try {
		iga::RequireCollectivePetscSuccess(comm,"volume Jacobian action",MatMult(assembly.Matrix(),assembly.State(),action));
		iga::CollectiveLocalStage(comm,"volume serial reference",[&] {
			auto insertion = references.try_emplace(std::make_pair(old_scale,expected_force_scale));
			auto& reference_residual = insertion.first->second.first;
			auto& reference_action = insertion.first->second.second;
			if (insertion.second) {
				reference_residual.assign(layout.Rows(),0); reference_action.assign(layout.Rows(),0);
				std::vector<std::array<double,4>> coefficients(layout.NodeIds().size());
				for (std::size_t i = 0; i < coefficients.size(); ++i)
					for (int c = 0; c < 4; ++c) coefficients[i][c] = Value(layout.NodeIds()[i],c,old_scale);
				const auto source_time = inputs.History().SourceTimeS(),target_time = inputs.History().TargetTimeS();
				const iga::ImmersedGlobalFlowState global(source_time,inputs.History().SourceIndex(),layout,coefficients,{0,0},true,0);
				const auto history = iga::BuildIdentityImmersedVelocityHistory(global,layout,target_time);
				const auto reference_force = [&](const std::array<double,3>& x) { auto value = field(x); for (auto& v : value) v *= expected_force_scale; return value; };
				for (const auto& stencil : stencils) {
					const auto element = f.domain.Background().MaterializeElement(stencil.id);
					std::vector<std::array<double,4>> current(element.connectivity.size());
					for (std::size_t a = 0; a < current.size(); ++a)
						for (int c = 0; c < 4; ++c) current[a][c] = Value(element.connectivity[a],c,.75);
					const auto quadrature = f.expanded.UsableRule(f.domain,stencil.id);
					const auto reference = iga::BuildTransientNavierStokesElement(element,current,history,layout,target_time,parameters,quadrature,reference_force,iga::NavierStokesResolvedMixedForm::Conservative);
					const auto zero = iga::BuildNavierStokesElementFromLocalVelocityHistory(element,current,iga::NavierStokesVelocityHistory(current.size()),parameters,quadrature,reference_force,iga::NavierStokesResolvedMixedForm::Conservative);
					for (std::size_t i = 0; i < stencil.rows.size(); ++i) {
						reference_residual[stencil.rows[i]] += reference.negative_residual[i];
						history_effect = std::max(history_effect,std::abs(reference.negative_residual[i]-zero.negative_residual[i]));
						for (std::size_t j = 0; j < stencil.rows.size(); ++j)
							reference_action[stencil.rows[i]] += reference.jacobian[i*stencil.rows.size()+j]*current[j/4][j%4];
					}
				}
			}
			iga::PetscReadArray residual,product; residual.Acquire(assembly.Residual()); product.Acquire(action);
			for (PetscInt row = assembly.RowBegin(); row < assembly.RowEnd(); ++row) {
				worst_error = std::max(worst_error,std::abs(residual.Data()[row-assembly.RowBegin()]-reference_residual[row]));
				worst_error = std::max(worst_error,std::abs(product.Data()[row-assembly.RowBegin()]-reference_action[row]));
			}
			residual.Restore(); product.Restore();
			Require(calls == frozen_points,"assembly called mutable force callback");
		});
		} catch (...) { VecDestroy(&action); throw; }
		iga::RequireCollectivePetscSuccess(comm,"volume action destroy",VecDestroy(&action));
	};
	verify(1,1);
	// Partial insertion failure is drained by the distributed backend. Retry
	// consumes the same immutable force/history candidate.
	std::size_t visits = 0;
	Reject(comm,[&] { assembly.Assemble([&](const auto& stencil,auto& matrix,auto& residual) {
		integrate(stencil,matrix,residual);
		if (rank == 0 && ++visits == (mode ? 2u : 1u)) throw std::runtime_error("injected transient volume assembly failure");
	}); });
	verify(1,1);
	inputs.ReleaseTrial(); fill(2); calls = 0; forbid = false; force_scale = 2;
	inputs.Freeze(assembly.State(),layout,.125,1,.25,2,parameters,force);
	forbid = true; fill(.75); verify(2,2);
	double global_error = 0,global_effect = 0;
	MPI_Allreduce(&worst_error,&global_error,1,MPI_DOUBLE,MPI_MAX,comm);
	MPI_Allreduce(&history_effect,&global_effect,1,MPI_DOUBLE,MPI_MAX,comm);
	iga::CollectiveLocalStage(comm,"transient volume numerical gate",[&] {
		Require(std::isfinite(global_error) && global_error < 1e-12 && global_effect > 1e-6,"transient volume numerical mismatch");
		Require(inputs.OwnedCells().size() == cells.size(),"transient cell ownership changed");
		if (mode && rank) Require(inputs.FrozenPointCount() == 0,"empty worker stores force data");
	});
	const auto timing = assembly.LastAssemblyTiming();
	inputs.Close(); inputs.Close();
	Reject(comm,[&] { inputs.Freeze(assembly.State(),layout,.125,1,.25,2,parameters,force); });
	assembly.Close();
	std::cout << "immersed_transient_volume_mpi rank=" << rank << " ranks=" << size << " mode=" << mode << " compact=" << compact
		<< " owned_cells=" << cells.size() << " frozen_points=" << frozen_points << " global_rows=" << layout.Rows()
		<< " maximum_error=" << global_error << " inertia_effect=" << global_effect
		<< " integration_s=" << timing.local_integration_insert_seconds << " halo_s=" << timing.halo_seconds
		<< " stash_s=" << timing.stash_exchange_seconds << " passed\n";
}
}
int main(int argc, char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0,status = 0; MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm comm = MPI_COMM_NULL;
	try {
		const bool split = argc > 1 && std::string(argv[1]) == "split";
		MPI_Comm_split(PETSC_COMM_WORLD,split && rank ? 1 : 0,rank,&comm);
		for (int mode = 0; mode < 3; ++mode) for (bool compact : {false,true}) Run(comm,mode,compact);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

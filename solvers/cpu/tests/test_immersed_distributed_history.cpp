#include "ImmersedDistributedVelocityHistory.hpp"
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
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground({{{0,0,0}},{{1,1,1}},{{4,1,1}}}),
		iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}
struct Fixture {
	iga::CartesianDomainClassification domain{Domain()};
	iga::CutCellVolumeQuadratureCatalog volume{domain,{2,500000,500000,3000000}};
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
void Run(MPI_Comm comm, int mode)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	std::unique_ptr<Fixture> fixture;
	iga::ImmersedActiveLayout layout,other;
	std::vector<std::int32_t> required;
	std::vector<std::uint64_t> cells;
	iga::CollectiveLocalStage(comm,"history fixture",[&] {
		fixture = std::make_unique<Fixture>();
		layout = iga::ImmersedActiveLayout::Build(fixture->domain,fixture->volume,"fixed",{8,9},true);
		other = iga::ImmersedActiveLayout::Build(fixture->domain,fixture->volume,"other",{8,9},true);
		for (std::uint64_t cell = 0; cell < 4; ++cell) {
			if ((mode ? 0 : static_cast<int>(cell)*size/4) != rank) continue;
			cells.push_back(cell);
			const auto element = fixture->domain.Background().MaterializeElement(cell);
			required.insert(required.end(),element.connectivity.begin(),element.connectivity.end());
		}
	});
	const auto nodes = static_cast<PetscInt>(layout.NodeIds().size()),rows = static_cast<PetscInt>(layout.Rows());
	const auto offset = [&](int r) -> PetscInt {
		if (mode == 1) return r ? rows : 0;
		if (mode == 2) return r == size ? rows : 0;
		return r == size ? rows : 4*(nodes*r/size);
	};
	const auto begin = offset(rank),end = offset(rank+1);
	Vec state = nullptr;
	iga::RequireCollectivePetscSuccess(comm,"history test vector",VecCreateMPI(comm,end-begin,rows,&state));
	try {
	const auto fill = [&](double scale, bool poison_pressure) {
		iga::CollectiveLocalStage(comm,"history test values",[&] {
			for (PetscInt row = begin; row < end; ++row) {
				const bool velocity = row < 4*nodes && row%4 < 3;
				const double value = velocity ? Value(layout.NodeIds()[row/4],row%4,scale)
					: poison_pressure ? std::numeric_limits<double>::quiet_NaN() : 1000+row;
				Require(!VecSetValue(state,row,value,INSERT_VALUES),"history test insert failed");
			}
		});
		iga::RequireCollectivePetscSuccess(comm,"history test vector begin",VecAssemblyBegin(state));
		iga::RequireCollectivePetscSuccess(comm,"history test vector end",VecAssemblyEnd(state));
	};
	fill(1,true);
	auto invalid_required = required;
	if (rank == size-1) invalid_required.push_back(-1);
	Reject(comm,[&] { iga::ImmersedDistributedVelocityHistory invalid(comm,layout,state,invalid_required); });
	iga::ImmersedDistributedVelocityHistory history(comm,layout,state,required);
	Reject(comm,[&] { history.Freeze(state,rank == size-1 ? other : layout,0,0,.125,1,.125); });
	Reject(comm,[&] { history.Freeze(state,layout,0,0,rank == size-1 ? .25 : .125,1,.125); });
	Reject(comm,[&] { history.Freeze(state,layout,0,rank == size-1 ? UINT64_MAX : 0,.125,1,.125); });
	if (size > 1) Reject(comm,[&] { history.Freeze(state,layout,0,0,rank == size-1 ? .25 : .125,1,rank == size-1 ? .25 : .125); });
	history.Freeze(state,layout,0,0,.125,1,.125);
	Reject(comm,[&] { history.Freeze(state,layout,0,0,.125,1,.125); });
	iga::Element absent; absent.connectivity = {-1};
	Reject(comm,[&] { (void)history.Localize(absent,.125); });
	// Poisoning only pressure/scalars is accepted; mutating every source value
	// after Freeze must not alter even a remotely fetched history node.
	iga::RequireCollectivePetscSuccess(comm,"history mutate source",VecSet(state,-42));
	double error = 0,effect = 0;
	const auto verify = [&](double scale,double source_time,double target_time,std::uint64_t index) {
		iga::CollectiveLocalStage(comm,"history numerical verification",[&] {
			Require(history.Active() && history.SourceTimeS() == source_time && history.TargetTimeS() == target_time
				&& history.SourceIndex() == index && history.TargetIndex() == index+1,"history clock differs");
			std::vector<std::array<double,4>> coefficients(nodes);
			for (PetscInt i = 0; i < nodes; ++i)
				for (int c = 0; c < 4; ++c) coefficients[i][c] = Value(layout.NodeIds()[i],c,scale);
			const iga::ImmersedGlobalFlowState global(source_time,index,layout,coefficients,{7,8},true,9);
			const auto serial = iga::BuildIdentityImmersedVelocityHistory(global,layout,target_time);
			for (const auto cell : cells) {
				const auto element = fixture->domain.Background().MaterializeElement(cell);
				const auto actual = history.Localize(element,target_time);
				const auto expected = iga::LocalizeImmersedVelocityHistory(element,layout,serial,target_time);
				Require(actual == expected,"history differs from global-ID serial reference");
				std::vector<std::array<double,4>> current(element.connectivity.size());
				for (std::size_t i = 0; i < current.size(); ++i)
					for (int c = 0; c < 4; ++c) current[i][c] = Value(element.connectivity[i],c,.75);
				const iga::NavierStokesParameters parameters{1,.1,.125};
				const auto force = [](const std::array<double,3>&) { return std::array<double,3>{{.1,.2,-.1}}; };
				const auto quadrature = fixture->volume.UsableRule(fixture->domain,cell);
				const auto distributed = iga::BuildTransientNavierStokesElement(element,current,history,target_time,parameters,quadrature,force);
				const auto reference = iga::BuildTransientNavierStokesElement(element,current,serial,layout,target_time,parameters,quadrature,force,iga::NavierStokesResolvedMixedForm::Conservative);
				const auto zero = iga::BuildNavierStokesElementFromLocalVelocityHistory(element,current,iga::NavierStokesVelocityHistory(current.size()),parameters,quadrature,force,iga::NavierStokesResolvedMixedForm::Conservative);
				for (std::size_t i = 0; i < distributed.negative_residual.size(); ++i) {
					error = std::max(error,std::abs(distributed.negative_residual[i]-reference.negative_residual[i]));
					effect = std::max(effect,std::abs(distributed.negative_residual[i]-zero.negative_residual[i]));
				}
				for (std::size_t i = 0; i < distributed.jacobian.size(); ++i)
					error = std::max(error,std::abs(distributed.jacobian[i]-reference.jacobian[i]));
				bool rejected = false;
				try { (void)history.Localize(element,std::nextafter(target_time,1.)); } catch (const std::exception&) { rejected = true; }
				Require(rejected,"history accepted wrong assembly time");
			}
		});
	};
	verify(1,0,.125,0);
	history.ReleaseTrial();
	Reject(comm,[&] { (void)history.SourceTimeS(); });
	// A bad velocity on any owner aborts collectively, even with no work there.
	fill(2,false);
	iga::CollectiveLocalStage(comm,"history inject one owner",[&] {
		if (begin == 0 && end > 0)
			Require(!VecSetValue(state,0,std::numeric_limits<double>::quiet_NaN(),INSERT_VALUES),"cannot inject history fault");
	});
	iga::RequireCollectivePetscSuccess(comm,"history fault begin",VecAssemblyBegin(state));
	iga::RequireCollectivePetscSuccess(comm,"history fault end",VecAssemblyEnd(state));
	Reject(comm,[&] { history.Freeze(state,layout,.125,1,.25,2,.125); });
	iga::CollectiveLocalStage(comm,"history failed freeze state",[&] { Require(!history.Active(),"failed history freeze published"); });
	fill(2,false); history.Freeze(state,layout,.125,1,.25,2,.125); verify(2,.125,.25,1);
	double global_error = 0,global_effect = 0;
	MPI_Allreduce(&error,&global_error,1,MPI_DOUBLE,MPI_MAX,comm);
	MPI_Allreduce(&effect,&global_effect,1,MPI_DOUBLE,MPI_MAX,comm);
	iga::CollectiveLocalStage(comm,"history numerical gates",[&] {
		Require(global_error == 0 && global_effect > 1e-6,"history element parity or inertia discriminator failed");
		if (mode == 2 && size > 1 && rank == 0) Require(history.OwnedNodes() == 0 && !required.empty(),"history worker must have only remote nodes");
		if (mode && rank) Require(history.RequiredNodeIds().empty(),"empty worker unexpectedly holds halo");
	});
	const auto owned = history.OwnedNodes(),halo = static_cast<PetscInt>(history.RequiredNodeIds().size());
	history.Close(); history.Close(); history.ReleaseTrial();
	Reject(comm,[&] { history.Freeze(state,layout,.125,1,.25,2,.125); });
	std::cout << "immersed_history_mpi rank=" << rank << " ranks=" << size << " mode=" << mode
		<< " owned_nodes=" << owned << " halo_nodes=" << halo << " global_nodes=" << nodes
		<< " element_error=" << global_error << " inertia_effect=" << global_effect << " passed\n";
	} catch (...) { VecDestroy(&state); throw; }
	iga::RequireCollectivePetscSuccess(comm,"history test destroy",VecDestroy(&state));
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
		for (int mode = 0; mode < 3; ++mode) Run(comm,mode);
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

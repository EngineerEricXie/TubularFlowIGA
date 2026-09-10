#include "ImmersedMovingTransientDistributedOperator.hpp"
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
iga::MaterialSurfacePatchMap BottomPatch(const iga::MaterialSurfaceKinematics& input)
{
	std::vector<iga::RawSurfaceTriangle> topology;
	for (const auto& source : input.SourceTriangles()) {
		iga::RawSurfaceTriangle triangle;triangle.boundary_id=source.boundary_id;
		for (int corner=0;corner<3;++corner)triangle.indices[corner]=source.source_vertex_indices[corner];
		topology.push_back(triangle);
	}
	const auto& reference=input.ReferenceMaterialVerticesM();
	const auto material=iga::MaterialSurfaceKinematics::CreateFromSourceTopology(reference,reference,
		std::vector<std::array<double,3>>(reference.size(),{{0,0,0}}),topology,0.,-.125,0.);
	iga::DistributedSurfaceLayout layout;
	layout.global_node_count=4;layout.partition_count=1;layout.partition_rank=0;
	layout.owned_global_node_ids={10,11,12,13};
	std::vector<iga::MaterialSurfacePatchMap::GlobalToSourceVertex> vertices;
	for (std::uint32_t node=0;node<4;++node) {
		vertices.emplace_back(10+node,node);
		layout.reference_positions.push_back({10+node,material.ReferenceMaterialVerticesM()[node]});
	}
	layout.reference_triangles={{{10,12,11}},{{10,13,12}}};
	layout.owned_reference_lumped_areas_m2={1./3.,1./6.,1./3.,1./6.};
	layout.reference_mesh_identity_sha256=iga::MaterialSurfacePatchMap::BuildReferenceIdentitySha256(
		material,8,vertices,layout.reference_triangles,{0,1},layout.owned_global_node_ids);
	layout.layout_identity_sha256=iga::BuildDistributedSurfaceLayoutIdentitySha256(layout);
	iga::DistributedSurfaceInterface interface;interface.id={"structure","membrane","bottom"};interface.subsystem_id="membrane";
	interface.boundary_labels={8};interface.reference_mesh_identity_sha256=layout.reference_mesh_identity_sha256;
	interface.provides={iga::SurfaceFieldQuantity::Displacement,iga::SurfaceFieldQuantity::Velocity};
	interface.requires={iga::SurfaceFieldQuantity::TractionOnStructure};
	return iga::MaterialSurfacePatchMap::Create(interface,layout,material,8,vertices,{0,1},layout.owned_global_node_ids);
}

void Run(MPI_Comm comm,const std::string& mode,bool background,bool aligned,bool deforming)
{
	int rank = 0,size = 1; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
	std::unique_ptr<iga::MovingCutGeometry> geometry,old_geometry;
	iga::ImmersedActiveLayout source_layout;
	std::unique_ptr<iga::ImmersedTransientFlowRuntime> reference;
	iga::ImmersedTransientFlowOptions options;
	std::vector<PetscScalar> committed,trial,expected_residual,expected_action;
	iga::ImmersedTransientFlowConservationDiagnostics conservation;
	iga::CollectiveLocalStage(comm,"transient operator reference",[&] {
		if (mode != "flow" && mode != "pressure" && mode != "traction" && mode != "closed" && mode != "inertial") throw std::invalid_argument("unknown transient operator mode");
		auto soup = Cube();for(auto& x:soup.vertices)for(double& coordinate:x)coordinate+=.1;
		auto moved=soup;
		for(auto& x:moved.vertices) {
			const double contraction=deforming?.03*(x[0]-.1):0.;
			x[0]+=(background ? (aligned ? .4 : .37) : .02)-contraction;
		}
		iga::PrescribedSurfaceMotion motion({{0,soup},{.125,moved}});
		iga::MovingCutGeometryOptions go; go.volume.max_depth = 2; go.volume_storage = iga::CutCellVolumeQuadratureStorageMode::Compact;
		const iga::CubicCartesianGridSpec grid{{{0,0,0}},{{background ? 2.4 : 1.2,1.2,1.2}},{{background ? 8u : 4u,3,3}}};
		old_geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(0,0,.125),go);
		geometry=iga::MovingCutGeometry::Build(grid,motion.Evaluate(.125,0,.125),go,old_geometry.get());
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
		source_layout=iga::ImmersedActiveLayout::Build(old_geometry->Domain(),old_geometry->Volume(),old_geometry->GeometryIdentitySha256(),layout.PortIds(),layout.HasGaugeRow());
		std::vector<std::array<double,4>> fields(source_layout.NodeIds().size());
		for(std::size_t i=0;i<fields.size();++i)for(int c=0;c<4;++c)fields[i][c]=Value(source_layout.NodeIds()[i],c,1.);
		const iga::ImmersedGlobalFlowState old_state(0,0,source_layout,fields,std::vector<double>(source_layout.PortIds().size()),source_layout.HasGaugeRow(),0.);
		const auto extension=iga::ImmersedVelocityExtension::Build(*old_geometry,source_layout,old_state,*geometry,layout,1);
		std::vector<double> pressure;for(const auto& value:fields)pressure.push_back(value[3]);
		const auto scalar=extension.ExtendScalar(pressure);
		std::vector<std::array<double,4>> seeded(layout.NodeIds().size());
		for(std::size_t i=0;i<seeded.size();++i) {
			for(int c=0;c<3;++c)seeded[i][c]=extension.TargetHistory().Velocities()[i][c];
			seeded[i][3]=scalar.target_values[i];
		}
		const iga::ImmersedGlobalFlowState seed(0,0,layout,seeded,std::vector<double>(layout.PortIds().size()),layout.HasGaugeRow(),0.);
		const auto map=iga::ImmersedMovingTrialMapIdentity::Create(old_state.HashSha256(),old_geometry->GeometryIdentitySha256(),geometry->GeometryIdentitySha256(),
			geometry->PublicationIdentitySha256(),source_layout.HashSha256(),layout.HashSha256(),extension.HashSha256(),extension.OperatorHashSha256(),
			extension.TargetHistory().HashSha256(),scalar.hash_sha256,layout.PortIds(),seed.PortMultipliers(),true,0.,0,.125,1,.125);
		committed.assign(source_layout.Rows(),0.);
		for(std::size_t i=0;i<fields.size();++i)for(int c=0;c<4;++c)committed[4*i+c]=fields[i][c];
		trial.resize(layout.Rows());for(std::size_t i=0;i<trial.size();++i)trial[i]=.01*std::sin(.13*(i+1));
		reference->BeginMovingTrial(.125,1,.125,extension.TargetHistory(),seed,map);
		reference->SetTrialState(trial);reference->Assemble();
		expected_residual = reference->AssembledNegativeResidual(); expected_action = reference->AssembledJacobianAction(trial);
		conservation = reference->ConservationDiagnostics();
		if (mode == "inertial") Require(reference->Diagnostics().wall_penalty.maximum_eta_t > 0,"reference inertial wall term is absent");
	});
	auto& g = *geometry;
	auto invalid = options;
	if (rank == size-1) invalid.wall_inertial_gamma0 = -1;
	Reject(comm,[&] { iga::ImmersedMovingTransientDistributedOperator rejected(comm,g,invalid); });
	if (size > 1) {
		auto different = options;
		if (rank == size-1) different.wall_inertial_gamma0 += .1;
		Reject(comm,[&] { iga::ImmersedMovingTransientDistributedOperator rejected(comm,g,different); });
	}
	const auto partition = background ? iga::ImmersedWorkPartition::FixedBackground : iga::ImmersedWorkPartition::CellCount;
	iga::ImmersedMovingTransientDistributedOperator op(comm,g,options,partition);
	if (background) {
		iga::ImmersedTransientDistributedOperator old_op(comm,old_geometry->Domain(),old_geometry->Volume(),
			old_geometry->Surface(),old_geometry->Ghost(),old_geometry->GeometryIdentitySha256(),options,partition);
		iga::CollectiveLocalStage(comm,"fixed background epoch ownership",[&] {
			Require(source_layout.NodeIds()!=op.Layout().NodeIds(),"moving fixture did not change active nodes");
			const auto count=g.Domain().Background().ElementCount();
			const auto begin=count/size*rank+std::min(count%size,static_cast<std::uint64_t>(rank));
			const auto end=count/size*(rank+1)+std::min(count%size,static_cast<std::uint64_t>(rank+1));
			for (auto cell:old_op.Inputs().OwnedCells()) Require(cell>=begin && cell<end,"source cell escaped fixed partition");
			for (auto cell:op.Inputs().OwnedCells()) Require(cell>=begin && cell<end,"target cell escaped fixed partition");
		});
		old_op.Close();
	}
	iga::DistributedImmersedVelocityExtension extension(comm,*old_geometry,source_layout,g,op.Layout(),1);
	Vec old = nullptr;
	iga::RequireCollectivePetscSuccess(comm,"transient test history vector",VecCreateMPI(comm,4*(source_layout.NodeIds().size()*(rank+1)/size-source_layout.NodeIds().size()*rank/size)+(rank==size-1?source_layout.Rows()-source_layout.NodeFieldRows():0),source_layout.Rows(),&old));
	try {
	const auto fill = [&](Vec vector,const std::vector<PetscScalar>& values) {
		iga::CollectiveLocalStage(comm,"transient test state insert",[&] {
			PetscInt first=0,last=0;Require(!VecGetOwnershipRange(vector,&first,&last),"vector ownership failed");
			for (PetscInt row = first; row < last; ++row)
				Require(!VecSetValue(vector,row,values[row],INSERT_VALUES),"transient vector insertion failed");
		});
		iga::RequireCollectivePetscSuccess(comm,"transient state begin",VecAssemblyBegin(vector));
		iga::RequireCollectivePetscSuccess(comm,"transient state end",VecAssemblyEnd(vector));
	};
	fill(old,committed); fill(op.Assembly().State(),trial);
	Reject(comm,[&] { op.Freeze(extension,old,source_layout,.01,0,.125,1,.115); });
	Reject(comm,[&] { op.Assemble(); });
	op.Freeze(extension,old,source_layout,0,0,.125,1,.125);
	if (!options.ports.empty()) Reject(comm,[&] { op.SetPortControlValue(options.ports[0].id,0); });
	const auto patch=BottomPatch(old_geometry->Evaluation());
	const auto captured=op.CaptureOwnedPatchState(patch);
	iga::CollectiveLocalStage(comm,"moving patch capture exact state",[&] {
		std::vector<std::uint64_t> expected;
		for (auto cell:op.Inputs().OwnedCells()) {
			if(g.Domain().Cells()[cell].classification!=iga::CellClassification::Cut)continue;
			const auto& points=g.Surface().UsableRule(g.Domain(),cell).Points();
			if(std::any_of(points.begin(),points.end(),[](const auto& p){return p.boundary_id==8;}))expected.push_back(cell);
		}
		Require(captured.size()==expected.size(),"patch capture cell coverage differs");
		for(std::size_t i=0;i<captured.size();++i) {
			Require(captured[i].cell_id==expected[i],"patch capture ownership differs");
			const auto element=g.Domain().Background().MaterializeElement(expected[i]);
			Require(captured[i].nodal_state.size()==element.connectivity.size(),"patch capture node count differs");
			for(std::size_t node=0;node<element.connectivity.size();++node)for(int field=0;field<4;++field)
				Require(captured[i].nodal_state[node][field]==PetscRealPart(trial[4*op.Layout().LocalNode(element.connectivity[node])+field]),"patch capture coefficient differs");
		}
	});

	std::array<double,4> errors{}; double physical_error = 0,geometric_rate_error=0;
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
					const auto& rule=g.Surface().UsableRule(g.Domain(),c);
					if(!std::any_of(rule.Points().begin(),rule.Points().end(),[&](const auto& point) { return point.boundary_id==port.boundary_label; }))continue;
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
		const auto material=op.MaterialConservationDiagnostics();
		iga::CollectiveLocalStage(comm,"moving material conservation comparison",[&] {
			Require(material.surface_flux_term_count==conservation.surface_flux_term_count,"moving flux term count differs");
			if(deforming) {
				// Unit cross section, affine contraction by 0.03 m in 0.125 s:
				// the endpoint material flux and exact volume rate are -0.24 m3/s.
				const double rate=(g.Diagnostics().closed_surface_physical_volume_m3
					-old_geometry->Diagnostics().closed_surface_physical_volume_m3)/.125;
				geometric_rate_error=std::max({geometric_rate_error,std::abs(rate+.24),
					std::abs(material.total_material_surface_outward_flow_m3_s+.24)});
				Require(geometric_rate_error<1e-12,"affine contraction violates analytic material volume rate");
			}
			const std::array<double iga::ImmersedTransientFlowConservationDiagnostics::*,13> fields{{
				&iga::ImmersedTransientFlowConservationDiagnostics::absolute_surface_flux_sum_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::endpoint_volume_divergence_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::total_surface_outward_flow_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::total_material_surface_outward_flow_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::total_material_wall_outward_flow_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::open_port_outward_flow_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::wall_outward_flow_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::wall_relative_leakage_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::divergence_theorem_defect_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::discrete_moving_wall_continuity_defect_m3_s,
				&iga::ImmersedTransientFlowConservationDiagnostics::normalized_open_balance,
				&iga::ImmersedTransientFlowConservationDiagnostics::normalized_wall_leakage,
				&iga::ImmersedTransientFlowConservationDiagnostics::normalized_discrete_moving_wall_continuity_defect}};
			std::size_t diagnostic_index=0;
			for(auto field:fields) {
				const double a=material.*field,b=conservation.*field;
				double tolerance=1e-11*std::max(1.,std::abs(b));
				if(diagnostic_index>=10) {
					const double actual_scale=diagnostic_index==12?material.discrete_moving_wall_continuity_normalization_scale_m3_s
						:std::max(options.flow_controller_reference_flow_m3_s,std::abs(material.open_port_outward_flow_m3_s));
					const double reference_scale=diagnostic_index==12?conservation.discrete_moving_wall_continuity_normalization_scale_m3_s
						:std::max(options.flow_controller_reference_flow_m3_s,std::abs(conservation.open_port_outward_flow_m3_s));
					// Propagate the unchanged raw-flux gate through the quotient;
					// near-zero denominators amplify reduction-order differences.
					tolerance=(1e-11+std::abs(b)*std::abs(actual_scale-reference_scale))/actual_scale
						+8*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(b));
				}
				if(!std::isfinite(a)||std::abs(a-b)>=tolerance) {
					std::ostringstream detail;detail<<std::setprecision(17)<<"moving material diagnostic "<<diagnostic_index<<" differs: actual="<<a<<" serial="<<b;
					throw std::runtime_error(detail.str());
				}
				++diagnostic_index;
			}
			for(const auto& item:conservation.material_surface_outward_flow_by_boundary_label_m3_s)
				Require(std::abs(material.material_surface_outward_flow_by_boundary_label_m3_s.at(item.first)-item.second)<1e-12,"moving material label flux differs");
			for(const auto& item:conservation.material_wall_outward_flow_by_boundary_label_m3_s)
				Require(std::abs(material.material_wall_outward_flow_by_boundary_label_m3_s.at(item.first)-item.second)<1e-12,"moving material wall label flux differs");
		});
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
		Reject(comm,[&] { op.Freeze(extension,old,source_layout,0,0,.125,1,.125); });
		op.SetPortControlValue("outlet",.0005); op.Freeze(extension,old,source_layout,0,0,.125,1,.125);
		expected_residual[op.Layout().ControllerRow(8)] += .0005; expected_residual[op.Layout().ControllerRow(9)] -= .0005;
		verify();
	}
	// A rank-local corrupt trial must reject collectively without closing the
	// operator; restoring the original owned field permits a healthy retry.
	iga::CollectiveLocalStage(comm,"moving diagnostic nonfinite state",[&] {
		if(rank==0) Require(!VecSetValue(op.Assembly().State(),op.Assembly().RowBegin(),
			std::numeric_limits<double>::quiet_NaN(),INSERT_VALUES),"nonfinite state insertion failed");
	});
	iga::RequireCollectivePetscSuccess(comm,"moving nonfinite state begin",VecAssemblyBegin(op.Assembly().State()));
	iga::RequireCollectivePetscSuccess(comm,"moving nonfinite state end",VecAssemblyEnd(op.Assembly().State()));
	Reject(comm,[&] { (void)op.MaterialConservationDiagnostics(); });
	fill(op.Assembly().State(),trial);verify();
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
	op.Close();extension.Close(); op.Close(); Reject(comm,[&] { op.Assemble(); });
	Reject(comm,[&] { (void)op.ConservationDiagnostics(); });
	Reject(comm,[&] { (void)op.MaterialConservationDiagnostics(); });
	iga::CollectiveLocalStage(comm,"transient operator close state",[&] {
		Require(!op.Inputs().History().Active() && op.Options().parameters.dt == 0,"closed transient input remains active");
	});
	std::cout << "immersed_moving_operator_mpi rank=" << rank << " ranks=" << size << " owned_rows=" << rows
		<< " halo_rows=" << halo << " owned_cells=" << work << " global_rows=" << op.Layout().Rows()
		<< " residual_relative_l2=" << residual_error << " action_relative_l2=" << action_error
		<< " deforming=" << deforming << " geometric_rate_error=" << geometric_rate_error << " mode=" << mode << " fixed_background=" << background << " translation_speed=" << (background ? (aligned ? 3.2 : 2.96) : .16) << " physical_error=" << global_physical << " passed\n";
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
		const std::string partition=argc > 3 ? argv[3] : "cell-count";
		Run(comm,argc > 1 ? argv[1] : "flow",partition == "background" || partition == "aligned" || partition == "deforming",partition == "aligned",partition == "deforming");
	} catch (const std::exception& error) { std::cerr << "rank " << rank << ": " << error.what() << '\n'; status = 1; }
	if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
	int global = 0; MPI_Allreduce(&status,&global,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);
	PetscFinalize(); return global;
}

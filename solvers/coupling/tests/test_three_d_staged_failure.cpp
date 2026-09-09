#include "BoundarySupport.hpp"
#include "ThreeDBodyFittedFlowTransportDomainAdapter.hpp"
#include <limits>
#include "../../cpu/tests/StringStreamFailure.hpp"

#define IGA_MULTIDOMAIN_FIXTURE_ONLY
#include "test_multidomain_flow_smoke.cpp"
#undef IGA_MULTIDOMAIN_FIXTURE_ONLY

namespace {

void Check(bool condition, const char* message)
{
	if (!condition) throw std::runtime_error(message);
}

template <class Work>
void Failure(MPI_Comm comm, const std::string& stage, Work&& work)
{
	std::string message;
	try { work(); } catch (const std::exception& error) { message = error.what(); }
	iga::CollectiveLocalStage(comm, "test diagnostic", [&] {
		if (message.find(stage) == std::string::npos) throw std::runtime_error("expected "+stage+", got: "+message);
	});
	iga::RequireCollectiveSameText(comm, "test common diagnostic", message);
}

void Near(double expected, double actual)
{
	if (!std::isfinite(expected) || !std::isfinite(actual)
		|| (expected == 0.0 ? std::abs(actual) > 1e-12 : std::abs((actual-expected)/expected) > 1e-6))
		throw std::runtime_error("staged quantity differs from reference");
}

std::vector<double> LocalState(Vec state)
{
	PetscInt count = 0;
	VecGetLocalSize(state, &count);
	iga::PetscReadArray view; view.Acquire(state);
	std::vector<double> result(count);
	for (PetscInt i = 0; i < count; ++i) result[i] = PetscRealPart(view.Data()[i]);
	return result;
}

std::vector<iga::CouplingPort> Ports()
{
	std::vector<iga::CouplingPort> ports;
	for (int label = 1; label <= 3; ++label) {
		iga::CouplingPort port;
		port.id = "port"+std::to_string(label); port.subsystem_id = "flow";
		port.locator_kind = "boundary_label"; port.locator = std::to_string(label);
		port.species = {"red", "blue"};
		port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
			iga::PortQuantity::MeanPressure, iga::PortQuantity::SpeciesConcentration, iga::PortQuantity::SpeciesFlux};
		port.requires = {iga::PortQuantity::SpeciesConcentration,
			label == 1 ? iga::PortQuantity::FlowRate : iga::PortQuantity::MeanPressure};
		ports.push_back(port);
	}
	return ports;
}

using Observations = std::map<std::string, std::vector<double>>;

void Compare(const Observations& expected, const Observations& actual)
{
	Check(expected.size() == actual.size(), "observation count differs");
	for (const auto& entry : expected) {
		const auto& values = actual.at(entry.first);
		Check(values.size() == entry.second.size(), "field size differs");
		double norm = 0.0, error = 0.0;
		for (std::size_t i = 0; i < values.size(); ++i) {
			Check(std::isfinite(values[i]) && std::isfinite(entry.second[i]), "nonfinite observation");
			norm = std::hypot(norm, entry.second[i]); error = std::hypot(error, values[i]-entry.second[i]);
		}
		Check(norm == 0.0 ? error <= 1e-12 : error/norm <= 1e-6, "staged field differs from reference");
	}
}

Observations Run(const fs::path& root, const fs::path& path, MPI_Comm comm, bool faults)
{
	int rank = 0, ranks = 1, cases = 0;
	MPI_Comm_rank(comm, &rank); MPI_Comm_size(comm, &ranks);
	iga::Database database(path.string());
	const auto mesh = iga::ReadLabeledHexMesh((root/"controlmesh.vtk").string(), 64, 1);
	const auto velocity = iga::ReadVelocity((root/"initial_velocityfield.txt").string(), 64);
	auto configuration = iga::ReadSimulationConfiguration((root/"simulation_config.json").string());
	const auto& definition = iga::FirstNavierStokesSystem(configuration);
	const auto boundaries = iga::ResolveFlowBoundaries(configuration, definition, mesh.labels, velocity);
	iga::TransientFlowRuntime flow(database, comm, true, true,
		{definition.density, definition.viscosity, kDt}, boundaries, mesh.labels,
		velocity, iga::WallTraceBasis(database, mesh, 0), {});
	flow.InitializeState(configuration);
	auto system = iga::CompileLinearSystem(configuration, "transport");
	system.terms.push_back({iga::TermKind::VolumeSource, 0, 0, 0.05, ""});
	system.terms.push_back({iga::TermKind::VolumeSource, 1, 1, 0.10, ""});
	// Construction copies the catalog. Faults must modify the runtime-owned
	// quadrature rather than the independent constructor input.
	std::map<std::uint64_t, iga::VolumeQuadratureRule> rules;
	iga::FullCell4x4x4VolumeQuadratureProvider provider(database.Load(0));
	rules.emplace(0, provider.Rule());
	iga::TransientTransportRuntime transport(database, comm, configuration, system, mesh.labels, std::move(rules));
	auto* weight = &const_cast<iga::VolumeQuadraturePoint&>(transport.VolumeRuleForTesting(0).Points().at(0)).weight;
	const auto original_weight = *weight;
	const double reference_flow = flow.ReferenceBoundaryFlow(1);
	Check(reference_flow < -1e-4, "fixture must have nonzero inlet flow");
	const auto ports = Ports();
	const std::map<std::string, std::string> bindings{{"red", "three_red"}, {"blue", "three_blue"}};
	const std::map<std::string, double> concentrations{{"red",1.25}, {"blue",3.5}};
	const std::map<std::string, double> reference{{"port1",reference_flow}};
	const iga::ThreeDFlowTransportDomainControls controls{{12,1e-8,1e-12,1e-6},1e-14};
	iga::ThreeDBodyFittedFlowTransportDomainAdapter adapter("flow",flow,transport,ports,configuration,root,bindings,reference,controls);
	const auto initial_flow = LocalState(flow.State());
	const auto initial_scalar = transport.GatherRequiredState();
	const auto check_initial = [&] {
		Check(flow.Phase() == iga::FlowStepPhase::Committed && transport.Phase() == iga::TransportStepPhase::Committed,
			"failed begin/abort left an open transaction");
		Check(LocalState(flow.State()) == initial_flow && transport.GatherRequiredState() == initial_scalar,
			"failed begin/abort changed initial snapshots");
		Check(transport.Steps() == 0, "failed trial advanced scalar clock");
	};
	const auto begin = [&] { adapter.BeginStep({0,0.0,kDt}); };
	const auto hydraulic = [&](double sign = 1.0) {
		for (int label = 1; label <= 3; ++label) {
			iga::PortBoundaryData input; input.time_s = kDt;
			if (label == 1) input.outward_flow_m3_s = sign*reference_flow;
			else input.mean_pressure_pa = 0.0;
			adapter.SetPortInput("port"+std::to_string(label),input);
		}
	};
	if (faults) {
		iga_test::ExpectStringStreamFailure(comm, "3d staged begin preparation", begin);
		check_initial(); cases += 3;
		Failure(comm,"checkpoint local read",[&] { transport.ReadState(root/"missing-initial.state"); });
		check_initial(); ++cases;
		// Invalid geometry is encountered after both transactions open. Begin
		// must coordinate the integral failure and abort both snapshots.
		if (rank == 0) *weight = -original_weight;
		Failure(comm,"transport mass integration",begin);
		if (rank == 0) *weight = original_weight;
		check_initial(); ++cases;
		if (ranks > 1) {
			for (int mode = 0; mode < 2; ++mode) {
				auto altered_ports = ports; auto altered_bindings = bindings;
				if (rank == ranks-1) {
					if (mode == 0) std::swap(altered_ports[0],altered_ports[1]);
					else std::swap(altered_bindings.at("red"),altered_bindings.at("blue"));
				}
				iga::ThreeDBodyFittedFlowTransportDomainAdapter altered("flow",flow,transport,
					altered_ports,configuration,root,altered_bindings,reference,controls);
				Failure(comm,"3d staged plan agreement",[&] { altered.BeginStep({0,0.0,kDt}); });
				check_initial(); ++cases;
			}
		}
	}
	begin();
	if (faults) {
		Failure(comm,"transport checkpoint write preparation",[&] { transport.WriteState(root/"trial-write.state"); });
		Check(!fs::exists(root/"trial-write.state") && transport.Steps() == 0
			&& transport.GatherRequiredState() == initial_scalar,"rejected trial write published or changed state"); ++cases;
		Failure(comm,"transport checkpoint preparation",[&] { transport.ReadState(root/"missing-trial.state"); });
		Check(transport.Phase() == iga::TransportStepPhase::TrialOpen && transport.Steps() == 0
			&& transport.GatherRequiredState() == initial_scalar,"rejected trial load changed transport state"); ++cases;
		Failure(comm,"3d accounting preparation",[&] { adapter.GetSpeciesStepAccounting(); }); ++cases;
		iga::PortBoundaryData mixed; mixed.time_s = kDt;
		mixed.outward_flow_m3_s = reference_flow; mixed.concentration = concentrations;
		iga_test::ExpectStringStreamFailure(comm, "3d staged input", [&] { adapter.SetPortInput("port1", mixed); });
		cases += 3;
		for (int mode = 0; mode < (ranks > 1 ? 3 : 2); ++mode) {
			auto invalid = mixed;
			if (rank == ranks-1) {
				if (mode == 0) invalid.time_s = 2*kDt;
				if (mode == 1) invalid.concentration.emplace("unknown",1.0);
				if (mode == 2) invalid.concentration.at("red") += 0.1;
			}
			Failure(comm,mode == 2 ? "3d staged input agreement" : "3d staged input",[&] { adapter.SetPortInput("port1",invalid); }); ++cases;
		}
		// This succeeds only if no successful peer published its rejected
		// mixed hydraulic/concentration candidate (duplicates are forbidden).
		adapter.SetPortInput("port1",mixed);
		adapter.AbortStep(); check_initial(); begin();
	}
	hydraulic(); adapter.SolveHydraulicTrial();
	const auto accepted_flow = LocalState(flow.State());
	if (faults) {
		Failure(comm,"3d staged species boundary",[&] { adapter.SolveTransportTrial(); }); ++cases;
		Check(LocalState(flow.State()) == accepted_flow && transport.Phase() == iga::TransportStepPhase::TrialOpen,
			"missing concentration changed hydraulic acceptance");
		for (int mode = 0; mode < (ranks > 1 ? 4 : 3); ++mode) {
			auto input = concentrations; double time = kDt;
			std::string id = "port1";
			if (rank == ranks-1) {
				if (mode == 0) input.erase("red");
				if (mode == 1) time *= 2;
				if (mode == 2) id = "unknown";
				if (mode == 3) input.at("red") += 0.1;
			}
			Failure(comm,mode == 3 ? "3d staged concentration agreement" : "3d staged concentration",
				[&] { adapter.SetTransportConcentration(id,time,input); }); ++cases;
		}
	}
	if (faults) {
		iga_test::ExpectStringStreamFailure(comm, "3d staged concentration", [&] {
			adapter.SetTransportConcentration("port1",kDt,concentrations);
		});
		cases += 3;
	}
	adapter.SetTransportConcentration("port1",kDt,concentrations);
	adapter.SolveTransportTrial();
	const auto scalar_trial = transport.GatherRequiredState();
	const auto accounting = adapter.GetSpeciesStepAccounting();
	Near(0.05*kDt,accounting.at("red").source_amount);
	Near(0.10*kDt,accounting.at("blue").source_amount);
	const auto masses = transport.TotalMass();
	const auto global_trial = transport.GatherState();
	if (faults) {
		iga_test::ExpectStringStreamFailure(comm, "transport global state preparation", [&] { transport.GatherState(); });
		iga_test::ExpectStringStreamFailure(comm, "transport mass preparation", [&] { transport.TotalMass(); });
		Check(transport.GatherState() == global_trial, "failed signature changed transport state");
		for (const auto& field : masses) Near(field.second, transport.TotalMass().at(field.first));
		cases += 6;
	}
	for (std::size_t node = 0; node < transport.RequiredNodes().size(); ++node)
		for (std::size_t field = 0; field < 2; ++field)
			Check(scalar_trial[node*2+field] == global_trial.at(transport.RequiredNodes()[node]*2+field),
				"global gather does not preserve node/field order");
	if (faults) {
		for (int mode = 0; mode < 4; ++mode) {
			if (mode == 3 && ranks == 1) continue;
			auto& fields = const_cast<iga::CompiledLinearSystem&>(transport.System()).fields;
			const auto saved = fields;
			if (rank == ranks-1) {
				if (mode == 0) fields.pop_back();
				if (mode == 1) fields[0] = fields[1];
				if (mode == 2) fields[0].clear();
				if (mode == 3) fields[0] = "another_field";
			}
			Failure(comm,mode == 3 ? "transport global state agreement" : "transport global state preparation",
				[&] { transport.GatherState(); });
			fields = saved;
			Check(transport.GatherState() == global_trial,"failed gather changed field or ordering");
			++cases;
		}
		// Illegal repeated calls must not clear accepted inputs or flags.
		Failure(comm,"3d staged hydraulic preparation",[&] { adapter.SolveHydraulicTrial(); });
		Failure(comm,"3d staged transport solve preparation",[&] { adapter.SolveTransportTrial(); });
		Failure(comm,"3d staged concentration",[&] { adapter.SetTransportConcentration("port1",kDt,concentrations); });
		cases += 3;
		for (int mode = 0; mode < 6; ++mode) {
			if (rank == 0 && mode < 2) *weight = -original_weight;
			auto& native = const_cast<iga::CompiledLinearSystem&>(transport.System());
			const auto saved_fields = native.fields; const auto saved_terms = native.terms;
			if (rank == ranks-1) {
				if (mode == 2) native.fields[0] = native.fields[1];
				if (mode == 3) native.terms.back().equation = native.fields.size();
				if (mode == 4) native.terms.back().coefficient += 0.1;
				if (mode == 5) native.fields[0] = "different";
			}
			if (mode >= 4 && ranks == 1) { native.fields = saved_fields; native.terms = saved_terms; continue; }
			const std::string stage = mode == 0 ? "transport mass integration" : mode == 1 ? "transport source integration"
				: mode == 2 ? "transport mass preparation" : mode == 3 ? "transport source preparation" : "transport integral agreement";
			Failure(comm,stage,[&] {
				if (mode == 1 || mode == 3 || mode == 4) transport.SourceIntegrals();
				else adapter.GetSpeciesStepAccounting();
			});
			*weight = original_weight; native.fields = saved_fields; native.terms = saved_terms; ++cases;
			for (const auto& field : masses) Near(field.second,transport.TotalMass().at(field.first));
		}
		auto global = transport.GatherState();
		if (rank == ranks-1) global.pop_back();
		Failure(comm,"transport mass preparation",[&] { transport.TotalMass(global); }); ++cases;
		global = transport.GatherState();
		if (rank == 0) global[0] = std::numeric_limits<double>::quiet_NaN();
		Failure(comm,"transport integral result",[&] { transport.TotalMass(global); }); ++cases;
		auto& native = const_cast<iga::CompiledLinearSystem&>(transport.System());
		const auto saved_terms = native.terms;
		for (int i = 0; i < 3; ++i)
			native.terms.push_back({iga::TermKind::VolumeSource,0,0,std::numeric_limits<double>::max(),""});
		Failure(comm,"transport integral result",[&] { transport.SourceIntegrals(); }); ++cases;
		native.terms = saved_terms;
		Check(transport.GatherRequiredState() == scalar_trial && LocalState(flow.State()) == accepted_flow,
			"accounting failure changed accepted fields");
	}
	adapter.RollbackTransportTrial();
	Check(LocalState(flow.State()) == accepted_flow && transport.GatherRequiredState() == initial_scalar,
		"scalar rollback changed hydraulic frame or snapshot");
	adapter.SetTransportConcentration("port1",kDt,concentrations); adapter.SolveTransportTrial();
	const auto retried = transport.GatherRequiredState();
	Observations result;
	for (std::size_t i = 0; i < retried.size(); ++i) {
		const auto key = i%2 == 0 ? "red_field" : "blue_field";
		result[key].push_back(retried[i]);
		Near(scalar_trial[i],retried[i]);
	}
	const auto final_accounting = adapter.GetSpeciesStepAccounting();
	for (const auto& entry : final_accounting) {
		Near(accounting.at(entry.first).initial_mass,entry.second.initial_mass);
		Near(accounting.at(entry.first).final_mass,entry.second.final_mass);
		result[entry.first+"_mass"] = {entry.second.initial_mass,entry.second.final_mass};
		result[entry.first+"_source"] = {entry.second.source_amount};
		for (const auto& amount : entry.second.outward_port_amount)
			result[entry.first+"_port_amount"].push_back(amount.second);
	}
	for (const auto& port : ports) {
		const auto value = adapter.GetTransportPortState(port.id);
		result["flow"].push_back(value.outward_flow_m3_s.value());
		result["pressure"].push_back(value.mean_pressure_pa.value());
	}
	adapter.PrepareCommitStep(); adapter.FinalizeCommitStep();
	const auto committed_flow = LocalState(flow.State());
	const auto committed_scalar = transport.GatherRequiredState();
	if (faults) {
		transport.WriteState(root/"readback.state");
		std::exception_ptr file_error;
		if (rank == 0) try {
			fs::copy_file(root/"readback.state",root/"truncated.state");
			fs::resize_file(root/"truncated.state",fs::file_size(root/"truncated.state")-1);
		} catch (...) { file_error = std::current_exception(); }
		iga::CollectiveLocalStage(comm,"transport checkpoint fixture",[&] {
			if (file_error) std::rethrow_exception(file_error);
		});
		// The native runtimes use valid configuration. Only this adapter replica
		// has a duplicate native species BC, encountered after its flow query.
		auto invalid = configuration;
		if (rank == ranks-1) invalid.boundaries[1].conditions.push_back(invalid.boundaries[1].conditions[1]);
		iga::ThreeDBodyFittedFlowTransportDomainAdapter bad("flow",flow,transport,ports,invalid,root,bindings,reference,controls);
		bad.BeginStep({1,kDt,kDt});
		for (int label = 1; label <= 3; ++label) {
			iga::PortBoundaryData input; input.time_s = 2*kDt;
			if (label == 1) input.outward_flow_m3_s = reference_flow; else input.mean_pressure_pa = 0.0;
			bad.SetPortInput("port"+std::to_string(label),input);
		}
		bad.SolveHydraulicTrial(); bad.SetTransportConcentration("port1",2*kDt,concentrations);
		Failure(comm,"3d staged species boundary",[&] { bad.SolveTransportTrial(); });
		bad.AbortStep(); ++cases;
		Check(LocalState(flow.State()) == committed_flow && transport.GatherRequiredState() == committed_scalar
			&& transport.Steps() == 1, "staged boundary failure did not restore committed state");
	}
	const auto reverse_hydraulics = [&] {
		adapter.BeginStep({1,kDt,kDt});
		for (int label = 1; label <= 3; ++label) {
			iga::PortBoundaryData input; input.time_s = 2*kDt;
			if (label == 1) input.outward_flow_m3_s = -reference_flow;
			else input.mean_pressure_pa = 0.0;
			adapter.SetPortInput("port"+std::to_string(label),input);
		}
		adapter.SolveHydraulicTrial();
		for (const auto& port : ports) {
			const auto value = adapter.GetHydraulicPortState(port.id).outward_flow_m3_s.value();
			Check(port.id == "port1" ? value > 1e-4 : value < -1e-4,"reversed port has wrong flow direction");
		}
	};
	if (faults) {
		reverse_hydraulics();
		adapter.SetTransportConcentration("port1",2*kDt,concentrations);
		Failure(comm,"3d staged species boundary",[&] { adapter.SolveTransportTrial(); });
		adapter.AbortStep(); ++cases;
		Check(LocalState(flow.State()) == committed_flow && transport.GatherRequiredState() == committed_scalar,
			"outward concentration rejection changed committed state");
	}
	reverse_hydraulics();
	for (const auto& id : {"port2","port3"}) adapter.SetTransportConcentration(id,2*kDt,concentrations);
	adapter.SolveTransportTrial();
	const auto reverse_state = transport.GatherRequiredState();
	for (std::size_t i = 0; i < reverse_state.size(); ++i)
		result[i%2 == 0 ? "reverse_red_field" : "reverse_blue_field"].push_back(reverse_state[i]);
	for (const auto& entry : adapter.GetSpeciesStepAccounting()) {
		result[entry.first+"_reverse_mass"] = {entry.second.initial_mass,entry.second.final_mass};
		result[entry.first+"_reverse_source"] = {entry.second.source_amount};
		for (const auto& amount : entry.second.outward_port_amount)
			result[entry.first+"_reverse_port_amount"].push_back(amount.second);
	}
	adapter.PrepareCommitStep(); adapter.FinalizeCommitStep();
	Check(transport.Steps() == 2,"accepted reverse step did not advance scalar clock");
	if (faults) {
		Failure(comm,"checkpoint write file preparation",[&] { transport.WriteState(root/"missing-parent"/"state"); });
		Check(transport.Steps() == 2 && transport.Phase() == iga::TransportStepPhase::Committed
			&& transport.GatherRequiredState() == reverse_state,"failed write changed transport state or clock"); ++cases;
		for (int mode = 0; mode < 3; ++mode) {
			fs::path checkpoint = root/(mode == 0 ? "missing-committed.state" : "truncated.state");
			if (mode == 2) {
				checkpoint = root/"readback.state";
				if (rank == ranks-1) checkpoint = checkpoint.string()+std::string(1,'\0');
			}
			Failure(comm,mode == 2 ? "checkpoint read preparation" : "checkpoint local read",
				[&] { transport.ReadState(checkpoint); }); ++cases;
			Check(transport.Steps() == 2 && transport.Phase() == iga::TransportStepPhase::Committed
				&& transport.GatherRequiredState() == reverse_state,"failed committed load changed state or clock");
		}
		Check(reverse_state != committed_scalar,"restart fixture needs different step states");
		transport.ReadState(root/"readback.state");
		Check(transport.Steps() == 1 && transport.GatherRequiredState() == committed_scalar,
			"successful load failed to restore field or warm-start clock");
		transport.BeginStep(); transport.AbortStep();
		Check(transport.Steps() == 1 && transport.GatherRequiredState() == committed_scalar,
			"abort after checkpoint load restored an obsolete snapshot");
	}
	if (rank == 0 && faults) std::cout << "three_d_staged_failure ranks=" << ranks << " cases=" << cases << " passed\n";
	return result;
}

} // namespace

int main(int argc,char** argv)
{
	PetscInitialize(&argc,&argv,nullptr,nullptr);
	int rank = 0, ranks = 1;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank); MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
	MPI_Comm group = MPI_COMM_NULL;
	try {
		Check(argc == 2 && ranks == 3,"usage: mpiexec -np 3 three_d_staged_failure_test FRESH_OUTPUT");
		for (int round = 0; round < 2; ++round) {
			const int color = round == 0 || rank == 0 ? 0 : 1;
			MPI_Comm_split(PETSC_COMM_WORLD,color,rank,&group);
			int local_rank = 0, local_size = 1;
			MPI_Comm_rank(group,&local_rank); MPI_Comm_size(group,&local_size);
			const auto root = fs::path(argv[1])/(std::to_string(round)+"-"+std::to_string(color));
			iga::CollectiveLocalStage(group,"fixture",[&] {
				if (local_rank == 0) {
					Check(fs::create_directories(root),"output already exists");
					WriteThreeDTransportCase(root); WriteDatabase(root/"group.ntiga",local_size); WriteDatabase(root/"serial.ntiga",1);
				}
			});
			Observations reference;
			iga::CollectiveLocalStage(group,"serial reference",[&] {
				if (local_rank == 0) reference = Run(root,root/"serial.ntiga",PETSC_COMM_SELF,false);
			});
			const auto actual = Run(root,root/"group.ntiga",group,true);
			iga::CollectiveLocalStage(group,"reference comparison",[&] { if (local_rank == 0) Compare(reference,actual); });
			MPI_Comm_free(&group);
		}
		if (rank == 0) std::cout << "three_d_staged_failure all comparisons passed\n";
	} catch (const std::exception& error) {
		std::cerr << "rank " << rank << ": " << error.what() << '\n';
		MPI_Abort(PETSC_COMM_WORLD,1); return 1;
	}
	PetscFinalize(); return 0;
}

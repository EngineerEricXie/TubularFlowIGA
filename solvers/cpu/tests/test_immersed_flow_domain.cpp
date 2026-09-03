#include "ThreeDImmersedFlowDomain.hpp"
#include "DomainRuntimeRegistry.hpp"
#include "SimulationGraph.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace {

iga::RawSurfaceTriangle Face(std::int64_t a, std::int64_t b, std::int64_t c, std::uint32_t label)
{ iga::RawSurfaceTriangle r; r.indices = {{a,b,c}}; r.boundary_id = label; return r; }

iga::CartesianDomainClassification Domain()
{
	iga::RawSurfaceSoup soup;
	soup.vertices = {{{{.1,.1,.1}},{{.9,.1,.1}},{{.9,.9,.1}},{{.1,.9,.1}},{{.1,.1,.9}},{{.9,.1,.9}},{{.9,.9,.9}},{{.1,.9,.9}}}};
	// The z caps are open ports and label 0 is the Nitsche wall.
	soup.triangles = {Face(0,2,1,1),Face(0,3,2,1),Face(4,5,6,2),Face(4,6,7,2), Face(0,1,5,0),Face(0,5,4,0),Face(1,2,6,0),Face(1,6,5,0),Face(2,3,7,0),Face(2,7,6,0),Face(3,0,4,0),Face(3,4,7,0)};
	const iga::CubicCartesianGridSpec spec{{{0,0,0}},{{1,1,1}},{{3,3,3}}};
	return iga::CartesianDomainClassification(iga::CubicCartesianBackground(spec), iga::SurfaceSpatialIndex(iga::ClosedTriangulatedSurface::Build(soup)));
}

template <class F> void Reject(F&& f)
{ bool rejected = false; try { f(); } catch (const std::exception&) { rejected = true; } assert(rejected); }

iga::CouplingPort Port(const char* id, int label, iga::PortQuantity input,
	std::set<iga::PortQuantity> provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate, iga::PortQuantity::MeanPressure, iga::PortQuantity::MeanNormalTraction})
{
	iga::CouplingPort p; p.id = id; p.subsystem_id = "immersed"; p.locator_kind = "boundary_label";
	p.locator = std::to_string(label); p.provides = std::move(provides); p.requires = {input}; return p;
}
iga::PortBoundaryData Flow(double time, double value)
{ iga::PortBoundaryData d; d.time_s = time; d.outward_flow_m3_s = value; return d; }
iga::PortBoundaryData Pressure(double time, double value = 0.0)
{ iga::PortBoundaryData d; d.time_s = time; d.mean_pressure_pa = value; return d; }
iga::PortBoundaryData Traction(double time, double value)
{ iga::PortBoundaryData d; d.time_s = time; d.mean_normal_traction_pa = value; return d; }

template <class T> void AppendBits(std::string& to, const T& value)
{
	static_assert(std::is_trivially_copyable<T>::value, "binary snapshot requires POD");
	to.append(reinterpret_cast<const char*>(&value), sizeof(value));
}
void AppendOptional(std::string& to, const std::optional<double>& value)
{ const bool has = value.has_value(); AppendBits(to, has); if (has) AppendBits(to, *value); }
std::string Serialize(const iga::PortState& s)
{
	std::string r; AppendBits(r, s.time_s); AppendOptional(r, s.area_m2); AppendOptional(r, s.outward_flow_m3_s);
	AppendOptional(r, s.mean_pressure_pa); AppendOptional(r, s.mean_normal_traction_pa); AppendOptional(r, s.total_pressure_pa);
	for (const auto& x : s.concentration) { const auto n = x.first.size(); AppendBits(r, n); r += x.first; AppendBits(r, x.second); }
	const auto concentration_end = s.concentration.size(); AppendBits(r, concentration_end);
	for (const auto& x : s.outward_species_flux) { const auto n = x.first.size(); AppendBits(r, n); r += x.first; AppendBits(r, x.second); }
	const auto flux_end = s.outward_species_flux.size(); AppendBits(r, flux_end); return r;
}
std::string Serialize(const std::map<std::string, iga::PortState>& states)
{
	std::string r; for (const auto& x : states) { const auto n = x.first.size(); AppendBits(r, n); r += x.first; r += Serialize(x.second); }
	const auto end = states.size(); AppendBits(r, end); return r;
}
bool SameBits(const std::vector<PetscScalar>& a, const std::vector<PetscScalar>& b)
{ return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size()*sizeof(PetscScalar)) == 0); }
std::string Work(const iga::ImmersedStaticFlowDiagnostics& d)
{
	std::string r; AppendBits(r, d.nonlinear_iterations); AppendBits(r, d.ksp_iterations); AppendBits(r, d.ksp_reason); AppendBits(r, d.residual_norm); AppendBits(r, d.damping); AppendBits(r, d.converged);
	for (const auto& s : d.newton_steps) { AppendBits(r, s.iteration); AppendBits(r, s.ksp_iterations); AppendBits(r, s.ksp_reason); AppendBits(r, s.residual_norm); AppendBits(r, s.update_norm); AppendBits(r, s.ksp_residual_norm); AppendBits(r, s.linear_residual_norm); AppendBits(r, s.linear_relative_residual); AppendBits(r, s.damping); AppendBits(r, s.candidate_residual_norm); }
	const auto end = d.newton_steps.size(); AppendBits(r, end); return r;
}
struct Counters { std::size_t ac, ap, ar, aa, bc, bp, bf, br; };
Counters Snapshot(const iga::ThreeDImmersedFlowDomain& a, const iga::ImmersedStaticFlowRuntime& r)
{ const auto& ad = a.Diagnostics(); const auto& rd = r.Diagnostics(); return {ad.committed_steps, ad.prepared_count, ad.rollback_count, ad.abort_count, rd.commit_count, rd.prepare_count, rd.finalize_count, rd.rollback_count}; }
bool SameCounters(const Counters& a, const Counters& b)
{ return a.ac == b.ac && a.ap == b.ap && a.ar == b.ar && a.aa == b.aa && a.bc == b.bc && a.bp == b.bp && a.bf == b.bf && a.br == b.br; }
void AssertCounterDelta(const Counters& before, const Counters& after,
	std::size_t ac = 0, std::size_t ap = 0, std::size_t ar = 0, std::size_t aa = 0,
	std::size_t bc = 0, std::size_t bp = 0, std::size_t bf = 0, std::size_t br = 0)
{
	assert(after.ac == before.ac+ac && after.ap == before.ap+ap && after.ar == before.ar+ar && after.aa == before.aa+aa);
	assert(after.bc == before.bc+bc && after.bp == before.bp+bp && after.bf == before.bf+bf && after.br == before.br+br);
}
iga::ImmersedStaticFlowOptions Options()
{
	iga::ImmersedStaticFlowOptions o; o.parameters = {1.0, 1.0, 0.0}; o.wall_labels = {0}; o.lu_pivot_shift = 1e-12;
	// A nonzero flow inlet plus pressure outlet makes this a real solve, not a vacuum fixture.
	o.ports = {{"inlet",1,iga::ImmersedFlowPortControlMode::FlowRate,-1e-3}, {"outlet",2,iga::ImmersedFlowPortControlMode::Pressure,0.0}}; return o;
}
void MainInputs(iga::ThreeDImmersedFlowDomain& a, double time, double flow, double pressure)
{ a.SetPortInput("inlet", Flow(time, flow)); a.SetPortInput("outlet", Pressure(time, pressure)); }

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr); int status = 0;
	try {
		const auto domain = Domain();
		// Transactional adapter semantics do not require deep cut-cell refinement;
		// dedicated aneurysm regressions retain that physical-refinement coverage.
		const iga::CutCellVolumeQuadratureCatalog volume(domain, {2,500000,500000,3000000});
		const iga::ImmersedSurfaceQuadratureCatalog surface(domain); const iga::CutCellGhostPenaltyCatalog ghost(domain, volume);
		const std::vector<iga::CouplingPort> metadata = {Port("inlet", 1, iga::PortQuantity::FlowRate), Port("outlet", 2, iga::PortQuantity::MeanPressure)};
		auto normal_metadata = Port("normal", 1, iga::PortQuantity::MeanNormalTraction);
		iga::ValidateThreeDImmersedFlowDomainMetadata("immersed", {normal_metadata});
		auto invalid_metadata = normal_metadata; invalid_metadata.requires = {iga::PortQuantity::TotalPressure};
		Reject([&] { iga::ValidateThreeDImmersedFlowDomainMetadata("immersed", {invalid_metadata}); });

		iga::ImmersedStaticFlowRuntime runtime(domain, volume, surface, ghost, Options());
		iga::ThreeDImmersedFlowDomain adapter("immersed", runtime, metadata);
		assert(adapter.Kind() == iga::DomainKind::ThreeDImmersedFlow);
		// Exact catalog binding: count, id, parsed label, canonical label, and mode all reject.
		Reject([&] { iga::ThreeDImmersedFlowDomain bad("immersed", runtime, {metadata[0]}); });
		auto wrong_id = metadata; wrong_id[0].id = "other"; Reject([&] { iga::ThreeDImmersedFlowDomain bad("immersed", runtime, wrong_id); });
		auto wrong_label = metadata; wrong_label[0].locator = "3"; Reject([&] { iga::ThreeDImmersedFlowDomain bad("immersed", runtime, wrong_label); });
		auto wrong_mode = metadata; wrong_mode[0].requires = {iga::PortQuantity::MeanPressure}; Reject([&] { iga::ThreeDImmersedFlowDomain bad("immersed", runtime, wrong_mode); });
		auto noncanonical = metadata; noncanonical[0].locator = "01"; Reject([&] { iga::ThreeDImmersedFlowDomain bad("immersed", runtime, noncanonical); });

		Reject([&] { adapter.SolveTrial(); }); Reject([&] { adapter.BeginStep({1,0.0,1.0}); });
		adapter.BeginStep({0,0.0,1.0}); const auto ready = Snapshot(adapter, runtime);
		Reject([&] { adapter.RollbackTrial(); }); // TrialReady cancellation is AbortStep, never rollback.
		assert(Snapshot(adapter, runtime).ar == ready.ar); adapter.AbortStep(); assert(adapter.Diagnostics().abort_count == ready.aa+1);

		adapter.BeginStep({0,0.0,1.0}); Reject([&] { adapter.BeginStep({1,1.0,1.0}); }); Reject([&] { adapter.SetPortInput("missing", Pressure(1.0)); }); Reject([&] { adapter.SetPortInput("inlet", Flow(2.0,-1e-3)); });
		iga::PortBoundaryData two_fields = Flow(1.0,-1e-3); two_fields.mean_pressure_pa = 0.0; Reject([&] { adapter.SetPortInput("inlet", two_fields); });
		const auto initial_state = adapter.CommittedBackendState(); const auto initial_ports = Serialize(adapter.CommittedPortStates());
		adapter.SetPortInput("inlet", Flow(1.0,-1e-3)); Reject([&] { adapter.SetPortInput("inlet", Flow(1.0,-2e-3)); }); // rejected input cannot overwrite accepted input
		adapter.SetPortInput("outlet", Pressure(1.0)); adapter.SolveTrial();
		const auto trial = runtime.TrialState(); assert(!SameBits(trial, initial_state));
		const auto inlet = adapter.GetPortState("inlet"), outlet = adapter.GetPortState("outlet");
		assert(inlet.time_s == 1.0 && inlet.area_m2 && inlet.outward_flow_m3_s && inlet.mean_pressure_pa && inlet.mean_normal_traction_pa);
		assert(outlet.time_s == 1.0 && outlet.area_m2 && outlet.outward_flow_m3_s && outlet.mean_pressure_pa && outlet.mean_normal_traction_pa);
		const auto port_trial = Serialize(inlet)+Serialize(outlet); const auto work = Work(runtime.Diagnostics());
		assert(runtime.Diagnostics().nonlinear_iterations > 0 || !runtime.Diagnostics().newton_steps.empty());
		const auto before_rollback = Snapshot(adapter, runtime); adapter.RollbackTrial();
		AssertCounterDelta(before_rollback, Snapshot(adapter, runtime), 0, 0, 1, 0, 0, 0, 0, 1);
		assert(SameBits(initial_state, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == initial_ports); Reject([&] { adapter.GetPortState("inlet"); });
		// A changed post-rollback input creates a different trial, while AbortStep
		// restores the precise committed image and publication counters.
		const auto before_changed_abort = Snapshot(adapter, runtime);
		MainInputs(adapter, 1.0, -2e-3, 0.5); adapter.SolveTrial(); assert(!SameBits(trial, runtime.TrialState())); adapter.AbortStep();
		const auto after_changed_abort = Snapshot(adapter, runtime);
		assert(SameBits(initial_state, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == initial_ports);
		AssertCounterDelta(before_changed_abort, after_changed_abort, 0, 0, 0, 1, 0, 0, 0, 1);
		assert(!runtime.Diagnostics().prepared && !runtime.Diagnostics().trial_active && runtime.Diagnostics().committed);
		adapter.BeginStep({0,0.0,1.0}); MainInputs(adapter, 1.0, -1e-3, 0.0); adapter.SolveTrial();
		assert(SameBits(trial, runtime.TrialState())); assert(port_trial == Serialize(adapter.GetPortState("inlet"))+Serialize(adapter.GetPortState("outlet"))); assert(work == Work(runtime.Diagnostics()));

		const auto before_failure_state = adapter.CommittedBackendState();
		const auto before_failure_ports = Serialize(adapter.CommittedPortStates());
		const auto before_failure = Snapshot(adapter, runtime);
		runtime.FailNextPrepareForTesting(); Reject([&] { adapter.PrepareCommitStep(); });
		// Failure preserves solved trial/prepared-port state, committed image/time, and all publication counters.
		assert(SameBits(trial, runtime.TrialState())); assert(port_trial == Serialize(adapter.GetPortState("inlet"))+Serialize(adapter.GetPortState("outlet"))); assert(work == Work(runtime.Diagnostics()));
		assert(SameBits(before_failure_state, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == before_failure_ports); assert(SameCounters(before_failure, Snapshot(adapter, runtime)));
		const auto before_prepare = Snapshot(adapter, runtime); adapter.PrepareCommitStep();
		const auto prepared = Snapshot(adapter, runtime); AssertCounterDelta(before_prepare, prepared, 0, 1, 0, 0, 0, 1);
		assert(SameBits(before_failure_state, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == before_failure_ports);
		const auto before_finalize = Snapshot(adapter, runtime); adapter.FinalizeCommitStep(); assert(SameBits(trial, adapter.CommittedBackendState())); assert(adapter.Diagnostics().committed_time_s == 1.0 && adapter.Diagnostics().committed_steps == 1);
		const auto finalized = Snapshot(adapter, runtime); AssertCounterDelta(before_finalize, finalized, 1, 0, 0, 0, 1, 0, 1);
		adapter.FinalizeCommitStep(); const auto finalized_twice = Snapshot(adapter, runtime); assert(SameCounters(finalized, finalized_twice));

		const auto committed = adapter.CommittedBackendState();
		const auto committed_ports = Serialize(adapter.CommittedPortStates());
		const auto committed_counters = Snapshot(adapter, runtime);
		adapter.BeginStep({1,1.0,1.0}); MainInputs(adapter, 2.0, -2e-3, 0.5); adapter.SolveTrial(); assert(!SameBits(committed, runtime.TrialState()));
		const auto before_solved_abort = Snapshot(adapter, runtime); adapter.AbortStep();
		AssertCounterDelta(before_solved_abort, Snapshot(adapter, runtime), 0, 0, 0, 1, 0, 0, 0, 1);
		assert(SameBits(committed, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == committed_ports); assert(adapter.Diagnostics().committed_time_s == 1.0 && Snapshot(adapter, runtime).ac == committed_counters.ac);
		assert(!runtime.Diagnostics().prepared && !runtime.Diagnostics().trial_active && runtime.Diagnostics().committed);
		adapter.BeginStep({1,1.0,1.0}); MainInputs(adapter, 2.0, -15e-4, 0.25); adapter.SolveTrial(); adapter.PrepareCommitStep(); Reject([&] { adapter.RollbackTrial(); });
		const auto before_prepared_abort = Snapshot(adapter, runtime); adapter.AbortStep();
		AssertCounterDelta(before_prepared_abort, Snapshot(adapter, runtime), 0, 0, 0, 1, 0, 0, 0, 1);
		assert(SameBits(committed, adapter.CommittedBackendState())); assert(Serialize(adapter.CommittedPortStates()) == committed_ports); assert(adapter.Diagnostics().committed_time_s == 1.0);
		assert(!runtime.Diagnostics().prepared && !runtime.Diagnostics().trial_active && runtime.Diagnostics().committed);

		// Normal-traction construction and input are executable, not metadata-only.
		auto traction_options = Options(); traction_options.ports = {{"normal",1,iga::ImmersedFlowPortControlMode::MeanNormalTraction,0.0}, {"outlet",2,iga::ImmersedFlowPortControlMode::Pressure,0.0}};
		iga::ImmersedStaticFlowRuntime traction_runtime(domain, volume, surface, ghost, traction_options);
		iga::ThreeDImmersedFlowDomain traction("immersed", traction_runtime, {Port("normal",1,iga::PortQuantity::MeanNormalTraction), Port("outlet",2,iga::PortQuantity::MeanPressure)});
		assert(!traction_runtime.HasGauge() && traction_runtime.GaugeDof() == -1);
		traction.BeginStep({0,0.0,1.0}); traction.SetPortInput("normal", Traction(1.0,.25)); traction.SetPortInput("outlet", Pressure(1.0)); traction.SolveTrial();
		const auto normal = traction.GetPortState("normal");
		const auto normal_diagnostic = std::find_if(traction_runtime.Diagnostics().ports.begin(), traction_runtime.Diagnostics().ports.end(), [](const auto& port) { return port.id == "normal"; });
		assert(normal_diagnostic != traction_runtime.Diagnostics().ports.end() && normal_diagnostic->target == .25);
		assert(normal.mean_normal_traction_pa && std::abs(*normal.mean_normal_traction_pa) > 1e-12);
		const auto traction_trial = traction_runtime.TrialState();
		assert(std::any_of(traction_trial.begin(), traction_trial.end(), [](PetscScalar value) { return PetscAbsScalar(value) > 1e-12; }));
		traction.AbortStep();

		// GetPortState exposes exactly declared outputs, never total pressure or species.
		iga::ImmersedStaticFlowRuntime subset_runtime(domain, volume, surface, ghost, Options());
		iga::ThreeDImmersedFlowDomain subset("immersed", subset_runtime, {Port("inlet",1,iga::PortQuantity::FlowRate,{iga::PortQuantity::Area}), Port("outlet",2,iga::PortQuantity::MeanPressure,{iga::PortQuantity::FlowRate})});
		subset.BeginStep({0,0.0,1.0}); MainInputs(subset, 1.0, -1e-3, 0.0); subset.SolveTrial(); const auto si = subset.GetPortState("inlet"), so = subset.GetPortState("outlet");
		assert(si.area_m2 && !si.outward_flow_m3_s && !si.mean_pressure_pa && !si.mean_normal_traction_pa && !si.total_pressure_pa && si.concentration.empty() && si.outward_species_flux.empty());
		assert(!so.area_m2 && so.outward_flow_m3_s && !so.mean_pressure_pa && !so.mean_normal_traction_pa && !so.total_pressure_pa && so.concentration.empty() && so.outward_species_flux.empty()); subset.AbortStep();

		const iga::SimulationGraph graph({{"immersed",iga::DomainKind::ThreeDImmersedFlow,metadata}}, {}); std::vector<std::unique_ptr<iga::CoupledDomainRuntime>> registered; registered.push_back(std::make_unique<iga::ThreeDImmersedFlowDomain>("immersed", runtime, metadata)); iga::DomainRuntimeRegistry registry(graph, std::move(registered)); assert(registry.Runtime("immersed").Kind() == iga::DomainKind::ThreeDImmersedFlow);
	} catch (const std::exception& error) { std::cerr << "immersed_flow_domain_test: " << error.what() << '\n'; status = 1; }
	PetscFinalize(); return status;
}

#include "ImmersedFlowCase.hpp"
#include "SurfaceReaders.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

iga::CouplingPort Port(const char* id, int label, iga::PortQuantity input)
{
	iga::CouplingPort port;
	port.id = id;
	port.subsystem_id = "immersed";
	port.locator_kind = "boundary_label";
	port.locator = std::to_string(label);
	port.provides = {iga::PortQuantity::Area, iga::PortQuantity::FlowRate,
		iga::PortQuantity::MeanPressure};
	port.requires = {input};
	return port;
}

double Norm(const std::vector<PetscScalar>& values, std::size_t block,
	std::size_t physical_dofs)
{
	double sum = 0.0;
	for (std::size_t row = 0; row < values.size(); ++row) {
		const std::size_t row_block = row >= physical_dofs ? 2 : ((row%4) == 3 ? 1 : 0);
		if (row_block == block) sum += PetscRealPart(values[row])*PetscRealPart(values[row]);
	}
	return std::sqrt(sum);
}

void CheckDirection(iga::ImmersedStaticFlowRuntime& runtime,
	const std::vector<PetscScalar>& zero, const std::vector<PetscScalar>& direction,
	const std::vector<PetscScalar>& action, double epsilon, const char* name)
{
	std::vector<PetscScalar> plus = zero, minus = zero;
	for (std::size_t i = 0; i < direction.size(); ++i) {
		plus[i] += epsilon*direction[i];
		minus[i] -= epsilon*direction[i];
	}
	runtime.SetCommittedState(plus); runtime.Assemble();
	const auto plus_rhs = runtime.AssembledNegativeResidual();
	runtime.SetCommittedState(minus); runtime.Assemble();
	const auto minus_rhs = runtime.AssembledNegativeResidual();
	std::vector<PetscScalar> defect(action.size());
	for (std::size_t row = 0; row < action.size(); ++row)
		defect[row] = action[row]+(plus_rhs[row]-minus_rhs[row])/(2.0*epsilon);
	const auto& diagnostics = runtime.Diagnostics();
	for (std::size_t block = 0; block < 3; ++block) {
		const double action_norm = Norm(action, block, diagnostics.physical_dofs);
		const double fd_norm = Norm(defect, block, diagnostics.physical_dofs);
		const double relative = action_norm > 0.0 ? fd_norm/action_norm : fd_norm;
		std::cout << "aneurysm zero-state centered-FD " << name << " output-block "
			<< block << " relative-defect=" << relative << '\n';
		if (action_norm == 0.0) {
			if (fd_norm > 1e-11)
				throw std::runtime_error(std::string(name)+" creates an unexpected block action");
			continue;
		}
		if (!(relative <= 1e-8))
			throw std::runtime_error(std::string(name)+" centered-FD Jacobian defect exceeds 1e-8");
	}
}

void CheckCapGeometry(const std::filesystem::path& root,
	const iga::ImmersedSurfaceQuadratureDiagnostics& surface_diagnostics)
{
	constexpr double expected_cap_area_m2 = 0.17677669529663687;
	for (const std::uint32_t label : {std::uint32_t{1}, std::uint32_t{2}}) {
		const auto area = surface_diagnostics.area_by_boundary_id.find(label);
		if (area == surface_diagnostics.area_by_boundary_id.end()
			|| std::abs(area->second-expected_cap_area_m2) > 2e-12)
			throw std::runtime_error("aneurysm cap area gate failed");
	}
	const auto source = iga::SurfaceReaders::ReadVtpPath((root/"surface.vtp").string());
	std::array<std::size_t, 3> cap_triangles{{0,0,0}};
	for (const auto& triangle : source.Triangles()) {
		if (triangle.boundary_id != 1 && triangle.boundary_id != 2) continue;
		const double expected_x = triangle.boundary_id == 1 ? -1.0 : 1.0;
		if (std::abs(triangle.outward_unit_normal[0]-expected_x) > 2e-14
			|| std::abs(triangle.outward_unit_normal[1]) > 2e-14
			|| std::abs(triangle.outward_unit_normal[2]) > 2e-14)
			throw std::runtime_error("aneurysm cap normal gate failed");
		++cap_triangles[triangle.boundary_id];
	}
	if (cap_triangles[1] != 8 || cap_triangles[2] != 8)
		throw std::runtime_error("aneurysm cap triangle-count gate failed");
}

void CheckNonnegativeFinite(double value, const char* name)
{
	if (!std::isfinite(value) || value < 0.0)
		throw std::runtime_error(std::string("aneurysm timing diagnostic is invalid: ")+name);
}

} // namespace

int main(int argc, char** argv)
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	int status = 0;
	try {
		const bool solve_only = argc == 3 && std::string(argv[1]) == "--solve-only";
		const bool supplied_case = argc == 2;
		if (argc != 1 && !solve_only && !supplied_case)
			throw std::invalid_argument("usage: immersed_aneurysm_jacobian_test [case-directory|--solve-only case-directory]");
		const auto root = solve_only ? std::filesystem::path(argv[2]) : supplied_case
			? std::filesystem::path(argv[1])
			: std::filesystem::path("../../examples/vascular_flow/immersed_aneurysm_chain/immersed");
		auto owner = iga::ImmersedFlowCase::Load(root, "immersed", {
			Port("inlet", 1, iga::PortQuantity::FlowRate),
			Port("outlet", 2, iga::PortQuantity::MeanPressure)}, 1);
		auto& runtime = owner->Runtime();
		runtime.SetPortControlValue("inlet", -1e-4);
		runtime.SetPortControlValue("outlet", 0.0);
		CheckCapGeometry(root, owner->SurfaceDiagnostics());
		const std::vector<PetscScalar> zero(runtime.Diagnostics().total_dofs, 0.0);
		runtime.SetCommittedState(zero); runtime.Assemble();

		if (!solve_only) {
			std::vector<PetscScalar> velocity(zero.size(), 0.0), pressure(zero.size(), 0.0), controller(zero.size(), 0.0);
			for (std::size_t node = 0; node < runtime.Diagnostics().active_nodes; ++node) {
				velocity[4*node] = 2e-4*(1.0+static_cast<double>(node%3));
				velocity[4*node+1] = -3e-4*(1.0+static_cast<double>(node%5));
				velocity[4*node+2] = 5e-4*(1.0+static_cast<double>(node%7));
				pressure[4*node+3] = 1.0+static_cast<double>(node%11);
			}
			controller[static_cast<std::size_t>(runtime.PortMultiplierDof("inlet"))] = 1.0;
			const auto velocity_action = runtime.AssembledJacobianAction(velocity);
			const auto pressure_action = runtime.AssembledJacobianAction(pressure);
			const auto controller_action = runtime.AssembledJacobianAction(controller);
			CheckDirection(runtime, zero, velocity, velocity_action, 1e-5, "velocity");
			CheckDirection(runtime, zero, pressure, pressure_action, 1e-5, "pressure");
			CheckDirection(runtime, zero, controller, controller_action, 1e-5, "controller");
		}

		runtime.SetCommittedState(zero);
		if (!runtime.SolveTrial()) throw std::runtime_error("aneurysm static solve did not converge");
		const auto& diagnostics = runtime.Diagnostics();
		if (diagnostics.newton_steps.empty()) throw std::runtime_error("aneurysm static solve accepted no Newton step");
		const auto& first = diagnostics.newton_steps.front();
		const bool newton_ksp_closed = diagnostics.converged
			&& diagnostics.nonlinear_iterations > 0
			&& first.ksp_reason > 0
			&& first.linear_relative_residual <= 1e-10
			&& first.candidate_residual_norm < first.residual_norm
			&& diagnostics.residual_norm <= 1e-11;
		const auto inlet = std::find_if(diagnostics.ports.begin(), diagnostics.ports.end(),
			[](const auto& port) { return port.id == "inlet"; });
		if (inlet == diagnostics.ports.end() || !std::isfinite(inlet->absolute_flow_residual_m3_s)
			|| inlet->absolute_flow_residual_m3_s > inlet->flow_tolerance_m3_s)
			throw std::runtime_error("aneurysm inlet flow target gate failed");
		const auto conservation = runtime.ConservationDiagnostics();
		if (!std::isfinite(conservation.open_port_outward_flow_m3_s)
			|| !std::isfinite(conservation.wall_outward_flow_m3_s)
			|| !std::isfinite(conservation.total_surface_outward_flow_m3_s)
			|| !std::isfinite(conservation.volume_divergence_integral_m3_s))
			throw std::runtime_error("aneurysm production-quadrature conservation diagnostic is not finite");
		if (conservation.surface_flow_by_boundary_label_m3_s.size() != 3
			|| conservation.surface_flow_by_boundary_label_m3_s.count(0) != 1
			|| conservation.surface_flow_by_boundary_label_m3_s.count(1) != 1
			|| conservation.surface_flow_by_boundary_label_m3_s.count(2) != 1)
			throw std::runtime_error("aneurysm production-quadrature label-flow map is incomplete");
		double port_sum = 0.0, label_sum = 0.0;
		for (const auto& entry : conservation.surface_flow_by_boundary_label_m3_s) label_sum += entry.second;
		for (const auto& port : diagnostics.ports) {
			const auto flow = conservation.surface_flow_by_boundary_label_m3_s.find(port.boundary_label);
			if (flow == conservation.surface_flow_by_boundary_label_m3_s.end()
				|| std::abs(flow->second-port.measurement.outward_flow_m3_s) > 1e-14)
				throw std::runtime_error("aneurysm production-quadrature port-flow diagnostic disagrees with port measurement");
			port_sum += flow->second;
		}
		if (std::abs(port_sum-conservation.open_port_outward_flow_m3_s) > 1e-14)
			throw std::runtime_error("aneurysm production-quadrature open-port sum is inconsistent");
		if (std::abs(conservation.surface_flow_by_boundary_label_m3_s.at(0)-conservation.wall_outward_flow_m3_s) > 1e-14
			|| std::abs(label_sum-conservation.total_surface_outward_flow_m3_s) > 1e-14
			|| std::abs(conservation.open_port_outward_flow_m3_s+conservation.wall_outward_flow_m3_s
				-conservation.total_surface_outward_flow_m3_s) > 1e-14)
			throw std::runtime_error("aneurysm production-quadrature surface-flow sums are inconsistent");
		const double throughflow_m3_s = std::abs(inlet->measurement.outward_flow_m3_s);
		if (!std::isfinite(throughflow_m3_s) || !(throughflow_m3_s > 0.0))
			throw std::runtime_error("aneurysm throughflow is invalid");
		const double open_normalized_balance = std::abs(conservation.open_port_outward_flow_m3_s)/throughflow_m3_s;
		const double wall_leakage_normalized = std::abs(conservation.wall_outward_flow_m3_s)/throughflow_m3_s;
		if (!(open_normalized_balance <= 1e-3) || !(wall_leakage_normalized <= 1e-3))
			throw std::runtime_error("aneurysm production-quadrature open/wall conservation gate failed");
		CheckNonnegativeFinite(diagnostics.last_assembly_seconds, "last assembly");
		CheckNonnegativeFinite(diagnostics.aggregate_assembly_seconds, "aggregate assembly");
		CheckNonnegativeFinite(diagnostics.last_linear_solve_seconds, "last linear solve");
		CheckNonnegativeFinite(diagnostics.aggregate_linear_solve_seconds, "aggregate linear solve");
		CheckNonnegativeFinite(diagnostics.last_line_search_candidate_assembly_seconds, "last line-search candidate assembly");
		CheckNonnegativeFinite(diagnostics.aggregate_line_search_candidate_assembly_seconds, "aggregate line-search candidate assembly");
		std::cout << "aneurysm static solve first-candidate=" << first.candidate_residual_norm
			<< " final-residual=" << diagnostics.residual_norm
			<< " first-linear-relative-residual=" << first.linear_relative_residual
			<< " inlet-target-error=" << inlet->absolute_flow_residual_m3_s << '\n';
		std::cout << "aneurysm conservation open-port-sum=" << conservation.open_port_outward_flow_m3_s
			<< " wall-surface-flow=" << conservation.wall_outward_flow_m3_s
			<< " total-surface-flow=" << conservation.total_surface_outward_flow_m3_s
			<< " volume-divergence=" << conservation.volume_divergence_integral_m3_s
			<< " open-normalized-balance=" << open_normalized_balance
			<< " wall-leakage-normalized=" << wall_leakage_normalized << '\n';
		std::cout << "aneurysm timings last-assembly=" << diagnostics.last_assembly_seconds
			<< " aggregate-assembly=" << diagnostics.aggregate_assembly_seconds
			<< " last-linear-solve=" << diagnostics.last_linear_solve_seconds
			<< " aggregate-linear-solve=" << diagnostics.aggregate_linear_solve_seconds
			<< " last-line-search-candidate-assembly=" << diagnostics.last_line_search_candidate_assembly_seconds
			<< " aggregate-line-search-candidate-assembly=" << diagnostics.aggregate_line_search_candidate_assembly_seconds << '\n';
		for (const auto& port : diagnostics.ports)
			std::cout << "aneurysm port " << port.id << " area=" << port.area_m2
				<< " mean-pressure=" << port.measurement.mean_pressure_pa
				<< " mean-traction=" << port.measurement.mean_normal_traction_pa << '\n';
		if (!newton_ksp_closed)
			throw std::runtime_error("aneurysm static Newton/KSP closure gate failed");
	} catch (const std::exception& error) {
		std::cerr << "immersed_aneurysm_jacobian_test: " << error.what() << '\n';
		status = 1;
	}
	PetscFinalize();
	return status;
}

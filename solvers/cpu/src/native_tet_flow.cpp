#include "NativeTetFem.hpp"
#include "CaseConfig.hpp"
#include "CheckedText.hpp"

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/resource.h>

namespace {

void Check(PetscErrorCode error, const char* context)
{
	if (error) throw std::runtime_error(std::string("PETSc failure in ")+context);
}

const iga::config_detail::JsonValue& Member(const iga::config_detail::JsonValue& value,
	const std::string& key, const std::string& context)
{
	const auto& object = iga::config_detail::RequireObject(value, context);
	const auto* member = iga::config_detail::Find(object, key);
	if (!member) throw std::runtime_error(context+" requires '"+key+"'");
	return *member;
}

struct Contract
{
	double density = 0.0, viscosity = 0.0, length = 0.0;
	double radius = 0.0, flow = 0.0, pressure_drop = 0.0;
	double maximum_mass_imbalance = 0.0, maximum_divergence_closure = 0.0;
	double maximum_outlet_backflow = 0.0;
};

Contract ReadContract(const std::string& path)
{
	std::ifstream input(path);
	if (!input) throw std::runtime_error("cannot open T2 contract: "+path);
	const auto root = iga::config_detail::JsonParser(iga::ReadCheckedText(input)).Parse();
	const auto& physics = Member(root, "physics", "contract");
	const auto& geometry = Member(root, "geometry", "contract");
	const auto& boundary = Member(root, "boundary_conditions", "contract");
	const auto& reference = Member(root, "analytic_reference", "contract");
	const auto& gates = Member(root, "qoi_gates", "contract");
	auto number = [](const auto& object, const char* key, const char* context) {
		return iga::config_detail::RequireNumber(Member(object, key, context),
			std::string(context)+"."+key);
	};
	Contract result;
	result.density = number(physics, "density_kg_m3", "physics");
	result.viscosity = number(physics, "dynamic_viscosity_pa_s", "physics");
	result.length = number(geometry, "length_m", "geometry");
	result.radius = number(geometry, "radius_m", "geometry");
	result.flow = number(boundary, "prescribed_flow_m3_s", "boundary_conditions");
	result.pressure_drop = number(reference, "pressure_drop_pa", "analytic_reference");
	result.maximum_mass_imbalance = number(gates,
		"maximum_relative_mass_imbalance", "qoi_gates");
	result.maximum_divergence_closure = number(gates,
		"maximum_surface_volume_divergence_relative_error", "qoi_gates");
	result.maximum_outlet_backflow = number(boundary,
		"maximum_outlet_backflow_area_fraction", "boundary_conditions");
	if (!(result.density > 0.0 && result.viscosity > 0.0 && result.length > 0.0
		&& result.radius > 0.0 && result.flow > 0.0 && result.pressure_drop > 0.0
		&& result.maximum_mass_imbalance > 0.0
		&& result.maximum_divergence_closure > 0.0
		&& result.maximum_outlet_backflow > 0.0))
		throw std::runtime_error("T2 contract physical values must be positive");
	return result;
}

std::array<double, 3> VelocityCoordinate(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology, std::uint32_t node)
{
	if (node < mesh.points.size()) return mesh.points[node];
	const auto& edge = topology.edges.at(node-mesh.points.size());
	std::array<double, 3> result{};
	for (int component = 0; component < 3; ++component)
		result[component] = 0.5*(mesh.points[edge[0]][component]+mesh.points[edge[1]][component]);
	return result;
}

std::vector<PetscInt> CellRows(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology, std::size_t cell)
{
	const auto velocity_nodes = mesh.points.size()+topology.edges.size();
	std::vector<PetscInt> rows;
	rows.reserve(iga::NativeTaylorHoodElementSystem::dofs);
	for (const auto node : topology.cell_velocity_nodes[cell])
		for (int component = 0; component < 3; ++component)
			rows.push_back(static_cast<PetscInt>(3*static_cast<std::uint64_t>(node)+component));
	for (const auto node : mesh.cells[cell].nodes)
		rows.push_back(static_cast<PetscInt>(3*velocity_nodes+node));
	return rows;
}

struct Metrics
{
	double inlet_flow = 0.0, outlet_flow = 0.0, wall_flow = 0.0;
	double pressure_drop = 0.0, velocity_relative_l2 = 0.0;
	double wall_shear_relative_l2 = 0.0, divergence_closure = 0.0;
	double mass_imbalance = 0.0, outlet_backflow_fraction = 0.0;
};

Metrics EvaluateMetrics(const iga::NativeTetMesh& mesh,
	const iga::NativeTaylorHoodTopology& topology, const PetscScalar* state,
	const Contract& contract)
{
	auto local_state = [&](std::size_t cell) {
		std::array<double, iga::NativeTaylorHoodElementSystem::dofs> result{};
		const auto rows = CellRows(mesh, topology, cell);
		for (std::size_t i = 0; i < rows.size(); ++i) result[i] = PetscRealPart(state[rows[i]]);
		return result;
	};
	double velocity_error = 0.0, velocity_reference = 0.0, volume_divergence = 0.0;
	const double pi = std::acos(-1.0);
	const double maximum_velocity = 2.0*contract.flow/(pi*contract.radius*contract.radius);
	for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
		const auto values = local_state(cell);
		const auto geometry = iga::EvaluateNativeTetGeometry(mesh, mesh.cells[cell]);
		for (const auto& point : iga::NativeTetDegreeFiveQuadrature()) {
			const auto basis = iga::EvaluateNativeTaylorHoodPhysicalBasis(geometry, point.reference);
			std::array<double,3> velocity{{0,0,0}}, coordinate{{0,0,0}};
			std::array<std::array<double,3>,3> gradient{};
			const std::array<double,4> lambda{{1.0-point.reference[0]-point.reference[1]
				-point.reference[2], point.reference[0], point.reference[1], point.reference[2]}};
			for (std::size_t a = 0; a < 4; ++a)
				for (int component = 0; component < 3; ++component)
					coordinate[component] += lambda[a]*mesh.points[mesh.cells[cell].nodes[a]][component];
			for (std::size_t a = 0; a < 10; ++a)
				for (int component = 0; component < 3; ++component) {
					velocity[component] += basis.velocity[a]*values[3*a+component];
					for (int derivative = 0; derivative < 3; ++derivative)
						gradient[component][derivative] +=
							basis.velocity_gradients[a][derivative]*values[3*a+component];
				}
			const double exact = maximum_velocity*(1.0-(coordinate[1]*coordinate[1]
				+coordinate[2]*coordinate[2])/(contract.radius*contract.radius));
			const double weight = point.weight*geometry.determinant;
			velocity_error += weight*((velocity[0]-exact)*(velocity[0]-exact)
				+velocity[1]*velocity[1]+velocity[2]*velocity[2]);
			velocity_reference += weight*exact*exact;
			volume_divergence += weight*(gradient[0][0]+gradient[1][1]+gradient[2][2]);
		}
	}
	using Face = std::array<std::uint32_t,3>;
	struct Owner { std::size_t cell; int opposite; };
	std::map<Face, Owner> owners;
	for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
		const auto& n = mesh.cells[cell].nodes;
		for (const auto& entry : {std::pair<Face,int>{{{n[1],n[2],n[3]}},0},
			std::pair<Face,int>{{{n[0],n[2],n[3]}},1},
			std::pair<Face,int>{{{n[0],n[1],n[3]}},2},
			std::pair<Face,int>{{{n[0],n[1],n[2]}},3}}) {
			auto face = entry.first;
			std::sort(face.begin(), face.end());
			if (!owners.emplace(face, Owner{cell, entry.second}).second) owners.erase(face);
		}
	}
	const std::array<double,4> gauss_x{{0.06943184420297371,0.33000947820757187,
		0.6699905217924281,0.9305681557970262}};
	const std::array<double,4> gauss_w{{0.17392742256872693,0.32607257743127307,
		0.32607257743127307,0.17392742256872693}};
	double inlet_area = 0.0, outlet_area = 0.0, inlet_pressure = 0.0, outlet_pressure = 0.0;
	double wall_shear_error = 0.0, wall_shear_reference = 0.0, backflow_area = 0.0;
	const double reference_shear = 4.0*contract.viscosity*contract.flow/
		(pi*contract.radius*contract.radius*contract.radius);
	Metrics result;
	for (const auto& triangle : mesh.boundary_triangles) {
		auto face = triangle.nodes;
		std::sort(face.begin(), face.end());
		const auto owner = owners.find(face);
		if (owner == owners.end()) throw std::runtime_error("cannot locate boundary triangle owner");
		const auto cell = owner->second.cell;
		const auto values = local_state(cell);
		const auto geometry = iga::EvaluateNativeTetGeometry(mesh, mesh.cells[cell]);
		const auto& p0 = mesh.points[triangle.nodes[0]];
		const auto& p1 = mesh.points[triangle.nodes[1]];
		const auto& p2 = mesh.points[triangle.nodes[2]];
		std::array<double,3> first{}, second{}, normal{};
		for (int c = 0; c < 3; ++c) { first[c] = p1[c]-p0[c]; second[c] = p2[c]-p0[c]; }
		normal = {{first[1]*second[2]-first[2]*second[1],
			first[2]*second[0]-first[0]*second[2], first[0]*second[1]-first[1]*second[0]}};
		const double surface_jacobian = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]
			+normal[2]*normal[2]);
		if (!(surface_jacobian > 0.0)) throw std::runtime_error("degenerate boundary triangle");
		const auto& opposite = mesh.points[mesh.cells[cell].nodes[owner->second.opposite]];
		double inward = 0.0;
		for (int c = 0; c < 3; ++c) inward += normal[c]*(opposite[c]-p0[c]);
		for (double& value : normal) value = (inward > 0.0 ? -value : value)/surface_jacobian;
		std::array<int,3> local_vertex{};
		for (int a = 0; a < 3; ++a) {
			local_vertex[a] = -1;
			for (int local = 0; local < 4; ++local)
				if (mesh.cells[cell].nodes[local] == triangle.nodes[a]) local_vertex[a] = local;
			if (local_vertex[a] < 0) throw std::runtime_error("boundary node is absent from owner cell");
		}
		for (std::size_t i = 0; i < 4; ++i)
			for (std::size_t j = 0; j < 4; ++j) {
				const double a = gauss_x[i], b = gauss_x[j];
				const std::array<double,3> triangle_lambda{{1.0-a,a*(1.0-b),a*b}};
				std::array<double,4> lambda{{0,0,0,0}};
				for (int vertex = 0; vertex < 3; ++vertex)
					lambda[local_vertex[vertex]] = triangle_lambda[vertex];
				const std::array<double,3> reference{{lambda[1],lambda[2],lambda[3]}};
				const auto basis = iga::EvaluateNativeTaylorHoodPhysicalBasis(geometry, reference);
				std::array<double,3> velocity{{0,0,0}};
				std::array<std::array<double,3>,3> gradient{};
				double pressure = 0.0;
				for (std::size_t node = 0; node < 10; ++node)
					for (int component = 0; component < 3; ++component) {
						velocity[component] += basis.velocity[node]*values[3*node+component];
						for (int derivative = 0; derivative < 3; ++derivative)
							gradient[component][derivative] +=
								basis.velocity_gradients[node][derivative]*values[3*node+component];
					}
				for (std::size_t node = 0; node < 4; ++node)
					pressure += basis.pressure[node]*values[30+node];
				const double weight = gauss_w[i]*gauss_w[j]*a*surface_jacobian;
				double flow = 0.0;
				for (int component = 0; component < 3; ++component)
					flow += velocity[component]*normal[component];
				if (triangle.boundary_label == 0) {
					result.wall_flow += weight*flow;
					std::array<double,3> traction{};
					for (int row = 0; row < 3; ++row)
						for (int column = 0; column < 3; ++column)
							traction[row] += ((row == column ? -pressure : 0.0)
								+contract.viscosity*(gradient[row][column]
								+gradient[column][row]))*normal[column];
					double normal_traction = 0.0;
					for (int component = 0; component < 3; ++component)
						normal_traction += traction[component]*normal[component];
					for (int component = 0; component < 3; ++component) {
						const double tangential = traction[component]-normal_traction*normal[component];
						const double exact = component == 0 ? -reference_shear : 0.0;
						wall_shear_error += weight*(tangential-exact)*(tangential-exact);
						wall_shear_reference += weight*exact*exact;
					}
				} else if (triangle.boundary_label == 1) {
					result.inlet_flow += weight*flow;
					inlet_area += weight;
					inlet_pressure += weight*pressure;
				} else if (triangle.boundary_label == 2) {
					result.outlet_flow += weight*flow;
					outlet_area += weight;
					outlet_pressure += weight*pressure;
					if (flow < 0.0) backflow_area += weight;
				} else throw std::runtime_error("unexpected boundary label in T2 metrics");
			}
	}
	result.pressure_drop = inlet_pressure/inlet_area-outlet_pressure/outlet_area;
	result.velocity_relative_l2 = std::sqrt(velocity_error/velocity_reference);
	result.wall_shear_relative_l2 = std::sqrt(wall_shear_error/wall_shear_reference);
	result.mass_imbalance = std::abs(result.inlet_flow+result.outlet_flow)/contract.flow;
	result.divergence_closure = std::abs(result.inlet_flow+result.outlet_flow+result.wall_flow
		-volume_divergence)/contract.flow;
	result.outlet_backflow_fraction = backflow_area/outlet_area;
	return result;
}

} // namespace

int main(int argc, char** argv)
{
	Check(PetscInitialize(&argc, &argv, nullptr, nullptr), "PetscInitialize");
	int exit_code = 0;
	try {
		int rank = 0, ranks = 1;
		Check(MPI_Comm_rank(PETSC_COMM_WORLD, &rank), "MPI_Comm_rank");
		Check(MPI_Comm_size(PETSC_COMM_WORLD, &ranks), "MPI_Comm_size");
		if (argc < 3) throw std::runtime_error(
			"usage: native_tet_flow volume.msh t2_fixed_flow_contract.json [PETSc options]");
		std::ifstream mesh_input(argv[1]);
		if (!mesh_input) throw std::runtime_error(std::string("cannot open mesh: ")+argv[1]);
		const auto mesh = iga::ReadNativeTetMeshGmsh41(mesh_input);
		const auto topology = iga::BuildNativeTaylorHoodTopology(mesh);
		const auto contract = ReadContract(argv[2]);
		PetscBool validation_mode = PETSC_FALSE;
		Check(PetscOptionsGetBool(nullptr, nullptr, "-native_tet_validation",
			&validation_mode, nullptr), "PetscOptionsGetBool validation");
		PetscReal nonlinear_absolute_tolerance = validation_mode ? 1.0e-14 : 1.0e-8;
		Check(PetscOptionsGetReal(nullptr, nullptr, "-native_tet_nonlinear_atol",
			&nonlinear_absolute_tolerance, nullptr), "PetscOptionsGetReal nonlinear atol");
		if (!(nonlinear_absolute_tolerance > 0.0)
			|| !std::isfinite(nonlinear_absolute_tolerance))
			throw std::runtime_error("native FEM nonlinear absolute tolerance must be positive");
		for (const int label : {0, 1, 2})
			if (topology.boundary_velocity_nodes.count(label) == 0)
				throw std::runtime_error("native FEM mesh requires wall/inlet/outlet labels 0/1/2");
		const std::uint64_t velocity_nodes = mesh.points.size()+topology.edges.size();
		const std::uint64_t dofs64 = 3*velocity_nodes+mesh.points.size();
		if (dofs64 > static_cast<std::uint64_t>(std::numeric_limits<PetscInt>::max()))
			throw std::runtime_error("native FEM system exceeds PETSc index range");
		const PetscInt dofs = static_cast<PetscInt>(dofs64);
		Mat jacobian = nullptr, schur_preconditioner = nullptr;
		Vec state = nullptr, residual = nullptr, negative_residual = nullptr, update = nullptr;
		Vec row_scale = nullptr, column_scale = nullptr;
		KSP solver = nullptr;
		Check(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, dofs, dofs,
			120, nullptr, 120, nullptr, &jacobian), "MatCreateAIJ");
		Check(MatSetOption(jacobian, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE), "MatSetOption");
		Check(MatCreateVecs(jacobian, &state, &residual), "MatCreateVecs");
		Check(VecDuplicate(residual, &negative_residual), "VecDuplicate residual");
		Check(VecDuplicate(state, &update), "VecDuplicate update");
		Check(VecDuplicate(state, &row_scale), "VecDuplicate row scale");
		Check(VecDuplicate(state, &column_scale), "VecDuplicate column scale");
		const double characteristic_velocity = contract.flow/
			(std::acos(-1.0)*contract.radius*contract.radius);
		const double characteristic_length = 2.0*contract.radius;
		const double pressure_scale = contract.viscosity*characteristic_velocity/
			characteristic_length;
		PetscInt scale_begin = 0, scale_end = 0;
		Check(VecGetOwnershipRange(row_scale, &scale_begin, &scale_end), "VecGetOwnershipRange scale");
		PetscScalar* row_values = nullptr;
		PetscScalar* column_values = nullptr;
		Check(VecGetArray(row_scale, &row_values), "VecGetArray row scale");
		Check(VecGetArray(column_scale, &column_values), "VecGetArray column scale");
		for (PetscInt global = scale_begin; global < scale_end; ++global) {
			const bool velocity_dof = global < static_cast<PetscInt>(3*velocity_nodes);
			row_values[global-scale_begin] = velocity_dof
				? 1.0/(contract.viscosity*characteristic_velocity)
				: 1.0/(characteristic_velocity*characteristic_length);
			column_values[global-scale_begin] = velocity_dof
				? characteristic_velocity : pressure_scale;
		}
		Check(VecRestoreArray(row_scale, &row_values), "VecRestoreArray row scale");
		Check(VecRestoreArray(column_scale, &column_values), "VecRestoreArray column scale");
		Check(KSPCreate(PETSC_COMM_WORLD, &solver), "KSPCreate");
		Check(KSPSetOptionsPrefix(solver, "native_tet_"), "KSPSetOptionsPrefix");
		Check(KSPSetType(solver, KSPPREONLY), "KSPSetType");
		PC pc = nullptr;
		Check(KSPGetPC(solver, &pc), "KSPGetPC");
		Check(PCSetType(pc, PCLU), "PCSetType");
		Check(KSPSetFromOptions(solver), "KSPSetFromOptions");
		PCType pc_type = nullptr;
		Check(PCGetType(pc, &pc_type), "PCGetType");
		if (pc_type && std::string(pc_type) == PCFIELDSPLIT) {
			IS velocity_is = nullptr, pressure_is = nullptr;
			PetscInt ownership_begin = 0, ownership_end = 0;
			Check(MatGetOwnershipRange(jacobian, &ownership_begin, &ownership_end),
				"MatGetOwnershipRange fieldsplit");
			const PetscInt velocity_end = static_cast<PetscInt>(3*velocity_nodes);
			std::vector<PetscInt> velocity_indices, pressure_indices;
			velocity_indices.reserve(static_cast<std::size_t>(ownership_end-ownership_begin));
			pressure_indices.reserve(static_cast<std::size_t>(ownership_end-ownership_begin));
			for (PetscInt row = ownership_begin; row < ownership_end; ++row)
				(row < velocity_end ? velocity_indices : pressure_indices).push_back(row);
			Check(ISCreateGeneral(PETSC_COMM_WORLD,
				static_cast<PetscInt>(velocity_indices.size()), velocity_indices.data(),
				PETSC_COPY_VALUES, &velocity_is), "ISCreateGeneral velocity");
			Check(ISCreateGeneral(PETSC_COMM_WORLD,
				static_cast<PetscInt>(pressure_indices.size()), pressure_indices.data(),
				PETSC_COPY_VALUES, &pressure_is), "ISCreateGeneral pressure");
			Check(PCFieldSplitSetIS(pc, "velocity", velocity_is), "PCFieldSplitSetIS velocity");
			Check(PCFieldSplitSetIS(pc, "pressure", pressure_is), "PCFieldSplitSetIS pressure");
			Check(PCSetUseAmat(pc, PETSC_TRUE), "PCSetUseAmat");
			const PetscInt pressure_dofs = static_cast<PetscInt>(mesh.points.size());
			const PetscInt local_pressure_dofs =
				static_cast<PetscInt>(pressure_indices.size());
			Check(MatCreateAIJ(PETSC_COMM_WORLD, local_pressure_dofs, local_pressure_dofs,
				pressure_dofs, pressure_dofs, 24, nullptr, 24, nullptr,
				&schur_preconditioner), "MatCreateAIJ Schur");
			Check(MatSetOption(schur_preconditioner, MAT_NEW_NONZERO_ALLOCATION_ERR,
				PETSC_FALSE), "MatSetOption Schur");
			for (std::size_t cell = static_cast<std::size_t>(rank); cell < mesh.cells.size();
					cell += static_cast<std::size_t>(ranks)) {
				const auto geometry = iga::EvaluateNativeTetGeometry(mesh, mesh.cells[cell]);
				std::array<PetscInt,4> rows{};
				std::array<PetscScalar,16> values{};
				for (std::size_t i = 0; i < 4; ++i) {
					rows[i] = static_cast<PetscInt>(mesh.cells[cell].nodes[i]);
					for (std::size_t j = 0; j < 4; ++j)
						values[4*i+j] = -geometry.determinant/
							(characteristic_length*characteristic_length
								*(i == j ? 60.0 : 120.0));
				}
				Check(MatSetValues(schur_preconditioner, 4, rows.data(), 4, rows.data(),
					values.data(), ADD_VALUES), "MatSetValues Schur");
			}
			Check(MatAssemblyBegin(schur_preconditioner, MAT_FINAL_ASSEMBLY),
				"MatAssemblyBegin Schur");
			Check(MatAssemblyEnd(schur_preconditioner, MAT_FINAL_ASSEMBLY),
				"MatAssemblyEnd Schur");
			Check(PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER,
				schur_preconditioner), "PCFieldSplitSetSchurPre");
			ISDestroy(&velocity_is);
			ISDestroy(&pressure_is);
		}

		std::set<PetscInt> wall_rows, inlet_rows;
		for (const auto node : topology.boundary_velocity_nodes.at(0))
			for (int component = 0; component < 3; ++component)
				wall_rows.insert(static_cast<PetscInt>(3*static_cast<std::uint64_t>(node)+component));
		for (const auto node : topology.boundary_velocity_nodes.at(1))
			for (int component = 0; component < 3; ++component)
				inlet_rows.insert(static_cast<PetscInt>(3*static_cast<std::uint64_t>(node)+component));
		std::set<PetscInt> all_rows = wall_rows;
		all_rows.insert(inlet_rows.begin(), inlet_rows.end());
		std::vector<PetscInt> boundary_rows;
		std::vector<PetscScalar> boundary_values;
		const double pi = std::acos(-1.0);
		const double maximum_velocity = 2.0*contract.flow/(pi*contract.radius*contract.radius);
		for (const auto row : all_rows) {
			const auto node = static_cast<std::uint32_t>(row/3);
			const int component = static_cast<int>(row%3);
			double value = 0.0;
			if (wall_rows.count(row) == 0 && component == 0) {
				const auto coordinate = VelocityCoordinate(mesh, topology, node);
				value = maximum_velocity*(1.0-(coordinate[1]*coordinate[1]
					+coordinate[2]*coordinate[2])/(contract.radius*contract.radius));
				if (value < 0.0 && value > -1.0e-12*maximum_velocity) value = 0.0;
				if (value < 0.0) throw std::runtime_error("inlet node lies outside the contract radius");
			}
			boundary_rows.push_back(row);
			boundary_values.push_back(value);
		}
		if (rank == 0) {
			for (std::uint32_t node = 0; node < velocity_nodes; ++node) {
				const auto coordinate = VelocityCoordinate(mesh, topology, node);
				double axial = maximum_velocity*(1.0-(coordinate[1]*coordinate[1]
					+coordinate[2]*coordinate[2])/(contract.radius*contract.radius));
				if (axial < 0.0) axial = 0.0;
				Check(VecSetValue(state, static_cast<PetscInt>(3*static_cast<std::uint64_t>(node)),
					axial, INSERT_VALUES), "VecSetValue initial velocity");
			}
			for (std::uint32_t node = 0; node < mesh.points.size(); ++node) {
				const PetscInt row = static_cast<PetscInt>(3*velocity_nodes+node);
				const double pressure = contract.pressure_drop*(1.0-mesh.points[node][0]/contract.length);
				Check(VecSetValue(state, row, pressure, INSERT_VALUES), "VecSetValue initial pressure");
			}
		}
		Check(VecAssemblyBegin(state), "VecAssemblyBegin state");
		Check(VecAssemblyEnd(state), "VecAssemblyEnd state");
		if (rank == 0)
			Check(VecSetValues(state, static_cast<PetscInt>(boundary_rows.size()),
				boundary_rows.data(), boundary_values.data(), INSERT_VALUES),
				"VecSetValues boundary state");
		Check(VecAssemblyBegin(state), "VecAssemblyBegin boundary state");
		Check(VecAssemblyEnd(state), "VecAssemblyEnd boundary state");
		PetscInt ownership_begin = 0, ownership_end = 0;
		Check(VecGetOwnershipRange(state, &ownership_begin, &ownership_end),
			"VecGetOwnershipRange boundary");
		std::vector<PetscInt> local_boundary_rows;
		std::vector<PetscScalar> local_boundary_values;
		for (std::size_t i = 0; i < boundary_rows.size(); ++i)
			if (boundary_rows[i] >= ownership_begin && boundary_rows[i] < ownership_end) {
				local_boundary_rows.push_back(boundary_rows[i]);
				local_boundary_values.push_back(boundary_values[i]);
			}

		VecScatter scatter = nullptr;
		Vec replicated = nullptr;
		Check(VecScatterCreateToAll(state, &scatter, &replicated), "VecScatterCreateToAll");
		const auto solve_start = std::chrono::steady_clock::now();
		double assembly_seconds = 0.0, linear_seconds = 0.0;
		double initial_norm = 0.0;
		PetscInt total_linear_iterations = 0;
		int completed_newton_iterations = 0;
		bool converged = false;
		KSPConvergedReason final_linear_reason = KSP_CONVERGED_ITERATING;
		for (int iteration = 0; iteration < 12; ++iteration) {
			const auto assembly_start = std::chrono::steady_clock::now();
			Check(VecScatterBegin(scatter, state, replicated, INSERT_VALUES, SCATTER_FORWARD),
				"VecScatterBegin state");
			Check(VecScatterEnd(scatter, state, replicated, INSERT_VALUES, SCATTER_FORWARD),
				"VecScatterEnd state");
			const PetscScalar* global_state = nullptr;
			Check(VecGetArrayRead(replicated, &global_state), "VecGetArrayRead state");
			Check(MatZeroEntries(jacobian), "MatZeroEntries");
			Check(VecSet(residual, 0.0), "VecSet residual");
			for (std::size_t cell = static_cast<std::size_t>(rank); cell < mesh.cells.size();
					cell += static_cast<std::size_t>(ranks)) {
				const auto rows = CellRows(mesh, topology, cell);
				std::array<double, iga::NativeTaylorHoodElementSystem::dofs> local_state{};
				for (std::size_t i = 0; i < rows.size(); ++i)
					local_state[i] = PetscRealPart(global_state[rows[i]]);
				const auto system = iga::BuildNativeTaylorHoodNavierStokesElement(mesh,
					mesh.cells[cell], local_state, {contract.density, contract.viscosity});
				Check(MatSetValues(jacobian, static_cast<PetscInt>(rows.size()), rows.data(),
					static_cast<PetscInt>(rows.size()), rows.data(), system.jacobian.data(), ADD_VALUES),
					"MatSetValues element");
				Check(VecSetValues(residual, static_cast<PetscInt>(rows.size()), rows.data(),
					system.residual.data(), ADD_VALUES), "VecSetValues element");
			}
			Check(VecRestoreArrayRead(replicated, &global_state), "VecRestoreArrayRead state");
			Check(MatAssemblyBegin(jacobian, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
			Check(MatAssemblyEnd(jacobian, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
			Check(VecAssemblyBegin(residual), "VecAssemblyBegin residual");
			Check(VecAssemblyEnd(residual), "VecAssemblyEnd residual");
			Check(MatZeroRowsColumns(jacobian,
				static_cast<PetscInt>(local_boundary_rows.size()), local_boundary_rows.data(),
				1.0, nullptr, nullptr), "MatZeroRowsColumns boundary");
			std::vector<PetscScalar> state_boundary(local_boundary_rows.size());
			Check(VecGetValues(state, static_cast<PetscInt>(local_boundary_rows.size()),
				local_boundary_rows.data(), state_boundary.data()), "VecGetValues boundary");
			for (std::size_t i = 0; i < state_boundary.size(); ++i)
				state_boundary[i] -= local_boundary_values[i];
			Check(VecSetValues(residual, static_cast<PetscInt>(local_boundary_rows.size()),
				local_boundary_rows.data(), state_boundary.data(), INSERT_VALUES),
				"VecSetValues boundary residual");
			Check(VecAssemblyBegin(residual), "VecAssemblyBegin constrained residual");
			Check(VecAssemblyEnd(residual), "VecAssemblyEnd constrained residual");
			PetscReal residual_norm = 0.0;
			Check(VecNorm(residual, NORM_2, &residual_norm), "VecNorm residual");
			assembly_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-assembly_start).count();
			if (iteration == 0) initial_norm = std::max(static_cast<double>(residual_norm), 1.0e-30);
			if (rank == 0) std::cout << "native_tet_newton iteration=" << iteration
				<< " residual=" << residual_norm << '\n';
			if (iteration > 0 && residual_norm <= std::max(
				static_cast<double>(nonlinear_absolute_tolerance), 1.0e-10*initial_norm)) {
				converged = true;
				completed_newton_iterations = iteration;
				break;
			}
			Check(MatDiagonalScale(jacobian, row_scale, column_scale), "MatDiagonalScale");
			Check(VecPointwiseMult(negative_residual, residual, row_scale),
				"VecPointwiseMult residual scale");
			Check(VecScale(negative_residual, -1.0), "VecScale residual");
			Check(KSPSetOperators(solver, jacobian, jacobian), "KSPSetOperators");
			const auto linear_start = std::chrono::steady_clock::now();
			Check(KSPSolve(solver, negative_residual, update), "KSPSolve");
			linear_seconds += std::chrono::duration<double>(
				std::chrono::steady_clock::now()-linear_start).count();
			KSPConvergedReason reason;
			PetscInt linear_iterations = 0;
			Check(KSPGetConvergedReason(solver, &reason), "KSPGetConvergedReason");
			final_linear_reason = reason;
			Check(KSPGetIterationNumber(solver, &linear_iterations), "KSPGetIterationNumber");
			if (reason <= 0) {
				PetscReal linear_residual = 0.0;
				Check(KSPGetResidualNorm(solver, &linear_residual), "KSPGetResidualNorm");
				throw std::runtime_error("native FEM linear solve did not converge: reason="
					+std::to_string(static_cast<int>(reason))+" iterations="
					+std::to_string(linear_iterations)+" residual="+std::to_string(linear_residual));
			}
			total_linear_iterations += linear_iterations;
			Check(VecPointwiseMult(update, update, column_scale),
				"VecPointwiseMult update scale");
			Check(VecAXPY(state, 1.0, update), "VecAXPY Newton update");
		}
		if (!converged) throw std::runtime_error("native FEM Newton solve did not converge");
		const double solve_seconds = std::chrono::duration<double>(
			std::chrono::steady_clock::now()-solve_start).count();
		Check(VecScatterBegin(scatter, state, replicated, INSERT_VALUES, SCATTER_FORWARD),
			"VecScatterBegin final state");
		Check(VecScatterEnd(scatter, state, replicated, INSERT_VALUES, SCATTER_FORWARD),
			"VecScatterEnd final state");
		struct rusage usage {};
		if (getrusage(RUSAGE_SELF, &usage))
			throw std::runtime_error("cannot query native FEM peak RSS");
		unsigned long long local_peak_rss =
			static_cast<unsigned long long>(usage.ru_maxrss)*1024ULL;
		unsigned long long peak_rss = 0;
		Check(MPI_Reduce(&local_peak_rss, &peak_rss, 1, MPI_UNSIGNED_LONG_LONG,
			MPI_MAX, 0, PETSC_COMM_WORLD), "MPI_Reduce peak RSS");
		if (rank == 0) {
			const PetscScalar* final_state = nullptr;
			Check(VecGetArrayRead(replicated, &final_state), "VecGetArrayRead final state");
			const auto metrics = EvaluateMetrics(mesh, topology, final_state, contract);
			Check(VecRestoreArrayRead(replicated, &final_state), "VecRestoreArrayRead final state");
			if (validation_mode && (metrics.mass_imbalance > contract.maximum_mass_imbalance
				|| metrics.divergence_closure > contract.maximum_divergence_closure
				|| metrics.outlet_backflow_fraction > contract.maximum_outlet_backflow))
				throw std::runtime_error("native FEM solution failed a conservation or backflow gate");
			char result_path[PETSC_MAX_PATH_LEN] = {};
			PetscBool write_result = PETSC_FALSE;
			Check(PetscOptionsGetString(nullptr, nullptr, "-native_tet_level_result",
				result_path, sizeof(result_path), &write_result), "PetscOptionsGetString result");
			if (write_result) {
				std::ofstream output(result_path);
				if (!output) throw std::runtime_error("cannot open native FEM level result");
				KSPType ksp_type = nullptr;
				PCType final_pc_type = nullptr;
				Check(KSPGetType(solver, &ksp_type), "KSPGetType result");
				Check(PCGetType(pc, &final_pc_type), "PCGetType result");
				output << std::setprecision(17)
					<< "{\n"
					<< "  \"schema_version\": 1,\n"
					<< "  \"result_classification\": \""
					<< (validation_mode ? "physical_validation_candidate" : "functional_smoke")
					<< "\",\n"
					<< "  \"validation_gates_enforced\": "
					<< (validation_mode ? "true" : "false") << ",\n"
					<< "  \"nonlinear_absolute_tolerance\": "
					<< nonlinear_absolute_tolerance << ",\n"
					<< "  \"backend\": \"native_cpp_petsc_tetrahedral_fem\",\n"
					<< "  \"backend_version\": \"PETSc " << PETSC_VERSION_MAJOR << '.'
					<< PETSC_VERSION_MINOR << '.' << PETSC_VERSION_SUBMINOR << "\",\n"
					<< "  \"velocity_space\": \"continuous tetrahedral Lagrange P2 vector\",\n"
					<< "  \"pressure_space\": \"continuous tetrahedral Lagrange P1 scalar\",\n"
					<< "  \"stabilization\": \"none (inf-sup-stable Taylor-Hood)\",\n"
					<< "  \"nonlinear_convergence_reason\": \"converged_newton_iterations_"
					<< completed_newton_iterations << "\",\n"
					<< "  \"linear_convergence_reason\": \"converged_petsc_reason_"
					<< static_cast<int>(final_linear_reason) << "\",\n"
					<< "  \"linear_solver\": \"" << (ksp_type ? ksp_type : "unknown")
					<< '+' << (final_pc_type ? final_pc_type : "unknown") << "\",\n"
					<< "  \"inlet_outward_flow_m3_s\": " << metrics.inlet_flow << ",\n"
					<< "  \"outlet_outward_flow_m3_s\": " << metrics.outlet_flow << ",\n"
					<< "  \"pressure_drop_pa\": " << metrics.pressure_drop << ",\n"
					<< "  \"velocity_relative_l2\": " << metrics.velocity_relative_l2 << ",\n"
					<< "  \"wall_shear_relative_l2\": " << metrics.wall_shear_relative_l2 << ",\n"
					<< "  \"outlet_backflow_area_fraction\": "
					<< metrics.outlet_backflow_fraction << ",\n"
					<< "  \"relative_mass_imbalance\": " << metrics.mass_imbalance << ",\n"
					<< "  \"surface_volume_divergence_relative_error\": "
					<< metrics.divergence_closure << ",\n"
					<< "  \"assembly_wall_s\": " << assembly_seconds << ",\n"
					<< "  \"solve_wall_s\": " << linear_seconds << ",\n"
					<< "  \"total_wall_s\": " << solve_seconds << ",\n"
					<< "  \"peak_rss_bytes\": " << peak_rss << ",\n"
					<< "  \"mpi_ranks\": " << ranks << ",\n"
					<< "  \"mesh_vertices\": " << mesh.points.size() << ",\n"
					<< "  \"velocity_edges\": " << topology.edges.size() << ",\n"
					<< "  \"tetrahedra\": " << mesh.cells.size() << ",\n"
					<< "  \"mixed_dofs\": " << dofs << ",\n"
					<< "  \"linear_iterations\": " << total_linear_iterations << "\n"
					<< "}\n";
				if (!output) throw std::runtime_error("cannot write native FEM level result");
			}
			std::cout << "native_tet_metrics inlet_flow=" << metrics.inlet_flow
				<< " outlet_flow=" << metrics.outlet_flow
				<< " pressure_drop=" << metrics.pressure_drop
				<< " velocity_relative_l2=" << metrics.velocity_relative_l2
				<< " wall_shear_relative_l2=" << metrics.wall_shear_relative_l2
				<< " mass_imbalance=" << metrics.mass_imbalance
				<< " divergence_closure=" << metrics.divergence_closure
				<< " outlet_backflow_fraction=" << metrics.outlet_backflow_fraction << '\n';
			std::cout << "native_tet_flow: "
				<< (validation_mode ? "VALIDATION_PASS" : "FUNCTIONAL_PASS")
				<< " nodes=" << mesh.points.size()
				<< " edges=" << topology.edges.size() << " cells=" << mesh.cells.size()
				<< " dofs=" << dofs << " newton=" << completed_newton_iterations
				<< " linear_iterations=" << total_linear_iterations
				<< " assembly_wall_s=" << assembly_seconds
				<< " linear_solve_wall_s=" << linear_seconds
				<< " total_wall_s=" << solve_seconds << '\n';
		}
		VecScatterDestroy(&scatter);
		VecDestroy(&replicated);
		KSPDestroy(&solver);
		VecDestroy(&update);
		VecDestroy(&negative_residual);
		VecDestroy(&residual);
		VecDestroy(&column_scale);
		VecDestroy(&row_scale);
		VecDestroy(&state);
		MatDestroy(&schur_preconditioner);
		MatDestroy(&jacobian);
	} catch (const std::exception& error) {
		int rank = 0;
		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
		if (rank == 0) std::cerr << "native_tet_flow: " << error.what() << '\n';
		exit_code = 1;
	}
	PetscFinalize();
	return exit_code;
}

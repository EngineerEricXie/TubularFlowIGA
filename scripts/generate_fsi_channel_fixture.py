#!/usr/bin/env python3
"""Generate a small labelled tetra fluid/wall case for the native FSI runner."""

import argparse
from collections import defaultdict
import json
from pathlib import Path


LENGTH, HEIGHT, WIDTH = 0.02, 0.004, 0.004


def node(i, j, k):
	return (i*2+j)*3+k


def determinant(points, nodes):
	a, b, c, d = (points[index] for index in nodes)
	x = [b[axis]-a[axis] for axis in range(3)]
	y = [c[axis]-a[axis] for axis in range(3)]
	z = [d[axis]-a[axis] for axis in range(3)]
	return (x[0]*(y[1]*z[2]-y[2]*z[1])
		-x[1]*(y[0]*z[2]-y[2]*z[0])
		+x[2]*(y[0]*z[1]-y[1]*z[0]))


def oriented(points, nodes):
	nodes = list(nodes)
	if determinant(points, nodes) < 0:
		nodes[1], nodes[2] = nodes[2], nodes[1]
	return tuple(nodes)


def fluid_mesh():
	points = [(LENGTH*i/2, HEIGHT*j, WIDTH*k/2)
		for i in range(3) for j in range(2) for k in range(3)]
	cells = []
	for i in range(2):
		for k in range(2):
			a, b, c, d = (node(i, 0, k), node(i+1, 0, k),
				node(i, 1, k), node(i+1, 1, k))
			e, f, g, h = (node(i, 0, k+1), node(i+1, 0, k+1),
				node(i, 1, k+1), node(i+1, 1, k+1))
			for vertices in ((a,b,d,h), (a,d,c,h), (a,c,g,h),
					(a,g,e,h), (a,e,f,h), (a,f,b,h)):
				cells.append(oriented(points, vertices))
	uses = defaultdict(list)
	for cell in cells:
		for opposite in range(4):
			face = tuple(cell[local] for local in range(4) if local != opposite)
			uses[tuple(sorted(face))].append(face)
	faces = []
	for occurrences in uses.values():
		if len(occurrences) != 1:
			continue
		face = occurrences[0]
		coordinates = [points[index] for index in face]
		if all(point[0] == 0 for point in coordinates):
			label = 1
		elif all(point[0] == LENGTH for point in coordinates):
			label = 2
		elif all(point[1] == HEIGHT for point in coordinates):
			label = 7
		else:
			label = 0
		faces.append((face, label))
	return points, cells, faces


def wall_mesh(fluid_points, fluid_faces):
	points = fluid_points.copy()
	cells, faces = [], []
	for interface, label in fluid_faces:
		if label != 7:
			continue
		apex = tuple(sum(fluid_points[index][axis] for index in interface)/3
			+(0.001 if axis == 1 else 0.0) for axis in range(3))
		apex_id = len(points)
		points.append(apex)
		cells.append(oriented(points, (*interface, apex_id)))
		faces.append((interface, 7))
		for first, second in ((0,1), (1,2), (2,0)):
			faces.append(((interface[first], interface[second], apex_id), 0))
	return points, cells, faces


def write_msh(path, points, cells, triangles):
	labels = sorted({label for _, label in triangles})
	lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$PhysicalNames",
		str(len(labels)+1)]
	lines += [f'2 {label+10} "boundary_label_{label}"' for label in labels]
	lines += ['3 1 "fluid"', '$EndPhysicalNames', '$Entities',
		f'0 0 {len(labels)} 1']
	for label in labels:
		selected = [points[index] for face, name in triangles if name == label
			for index in face]
		low = [min(point[axis] for point in selected) for axis in range(3)]
		high = [max(point[axis] for point in selected) for axis in range(3)]
		lines.append(' '.join(map(str, [label+10, *low, *high, 1, label+10, 0])))
	low = [min(point[axis] for point in points) for axis in range(3)]
	high = [max(point[axis] for point in points) for axis in range(3)]
	lines += [' '.join(map(str, [1, *low, *high, 1, 1, 0])), '$EndEntities',
		'$Nodes', f'1 {len(points)} 1 {len(points)}', f'3 1 0 {len(points)}']
	lines += [str(index+1) for index in range(len(points))]
	lines += [' '.join(map(str, point)) for point in points]
	lines += ['$EndNodes', '$Elements',
		f'{len(labels)+1} {len(triangles)+len(cells)} 1 {len(triangles)+len(cells)}']
	element = 1
	for label in labels:
		selected = [face for face, name in triangles if name == label]
		lines.append(f'2 {label+10} 2 {len(selected)}')
		for face in selected:
			lines.append(' '.join(map(str, [element, *(index+1 for index in face)])))
			element += 1
	lines.append(f'3 1 4 {len(cells)}')
	for cell in cells:
		lines.append(' '.join(map(str, [element, *(index+1 for index in cell)])))
		element += 1
	lines.append('$EndElements')
	path.write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('output_directory', type=Path)
	args = parser.parse_args()
	output = args.output_directory
	output.mkdir(parents=True, exist_ok=True)
	fluid_points, fluid_cells, fluid_faces = fluid_mesh()
	wall_points, wall_cells, wall_faces = wall_mesh(fluid_points, fluid_faces)
	write_msh(output/'fluid.msh', fluid_points, fluid_cells, fluid_faces)
	write_msh(output/'wall.msh', wall_points, wall_cells, wall_faces)
	interface = {index for face, label in wall_faces if label == 7 for index in face}
	fixed = [index for index, point in enumerate(wall_points)
		if index not in interface or point[0] in (0, LENGTH)
		or point[2] in (0, WIDTH)]
	case = {'schema_version': 1, 'fluid_mesh_file': 'fluid.msh',
		'solid_mesh_file': 'wall.msh', 'interface_label': 7,
		'fluid': {'density_kg_m3': 1000.0, 'dynamic_viscosity_pa_s': 0.004,
			'initial_edge_velocity_x_m_s': 0.0,
			'boundary_velocity_by_label_m_s': {'1': [0.2, 0.0, 0.0]},
			'natural_boundary_labels': [2]},
		'solid': {'young_modulus_pa': 20000.0, 'poisson_ratio': 0.3,
			'density_kg_m3': 1000.0, 'fixed_nodes': fixed},
		'time': {'dt_s': 0.01, 'steps': 2}}
	(output/'case.json').write_text(json.dumps(case, indent=2)+'\n', encoding='utf-8')
	fsi_species_case = {**case, 'species': [{
		'id': 'tracer', 'initial_concentration_mol_m3': 0.2,
		'diffusivity_m2_s': 0.01, 'source_mol_m3_s': 0.0,
		'first_order_decay_rate_s_inv': 0.1,
		'inflow_concentration_by_label_mol_m3': {
			'0': 1.0, '1': 1.0, '2': 1.0, '7': 1.0}},
		{'id': 'oxygen', 'initial_concentration_mol_m3': 0.4,
		'diffusivity_m2_s': 0.002, 'source_mol_m3_s': 0.0,
		'first_order_decay_rate_s_inv': 0.05,
		'inflow_concentration_by_label_mol_m3': {
			'0': 0.8, '1': 0.8, '2': 0.8, '7': 0.8}}]}
	(output/'fsi_species_case.json').write_text(
		json.dumps(fsi_species_case, indent=2)+'\n', encoding='utf-8')
	fsi_species_solver = {'schema_version': 'solver-v1',
		'time': {'mode': 'transient', 'dt_s': 0.01, 'steps': 2},
		'domains': [{'id': 'vessel', 'dimension': '3d', 'method': 'FEM',
			'device': 'CPU', 'backend': 'native_tet_fsi_steps',
			'physics': ['flow', 'solid', 'species'],
			'case_file': 'fsi_species_case.json'}],
		'connections': [], 'resources': {'mpi_ranks': 1},
		'output': {'directory': 'fsi_species_results', 'every_steps': 1}}
	(output/'fsi_species_solver.json').write_text(
		json.dumps(fsi_species_solver, indent=2)+'\n', encoding='utf-8')
	gpu_case = {'schema_version': 1, 'mesh_file': 'fluid.msh',
		'fluid': {'density_kg_m3': 1000.0, 'dynamic_viscosity_pa_s': 0.004,
			'initial_velocity_m_s': [0.0, 0.0, 0.0],
			'nonlinear_tolerance': 1e-8, 'maximum_iterations': 12},
		'boundaries': {'velocity_by_label_m_s': {
			'0': [0.0, 0.0, 0.0], '1': [0.2, 0.0, 0.0],
			'7': [0.0, 0.0, 0.0]},
			'pressure_by_label_pa': {'2': 0.0}},
		'time': {'dt_s': 0.01, 'steps': 2}}
	(output/'gpu_flow_case.json').write_text(
		json.dumps(gpu_case, indent=2)+'\n', encoding='utf-8')
	gpu_solver = {'schema_version': 'solver-v1',
		'time': {'mode': 'transient', 'dt_s': 0.01, 'steps': 2},
		'domains': [{'id': 'fluid', 'dimension': '3d', 'method': 'FEM',
			'device': 'GPU', 'backend': 'native_tet_flow_cuda',
			'physics': ['flow'], 'case_file': 'gpu_flow_case.json'}],
		'connections': [], 'resources': {'mpi_ranks': 1},
		'output': {'directory': 'gpu_flow_results'}}
	(output/'gpu_flow_solver.json').write_text(
		json.dumps(gpu_solver, indent=2)+'\n', encoding='utf-8')
	gpu_fsi_solver = {'schema_version': 'solver-v1',
		'time': {'mode': 'transient', 'dt_s': 0.01, 'steps': 2},
		'domains': [{'id': 'vessel', 'dimension': '3d', 'method': 'FEM',
			'device': 'GPU', 'backend': 'native_tet_fsi_cuda',
			'physics': ['flow', 'solid'], 'case_file': 'case.json'}],
		'connections': [], 'resources': {'mpi_ranks': 1},
		'output': {'directory': 'gpu_fsi_results', 'every_steps': 1}}
	(output/'gpu_fsi_solver.json').write_text(
		json.dumps(gpu_fsi_solver, indent=2)+'\n', encoding='utf-8')
	gpu_species_case = {'schema_version': 1, 'mesh_file': 'fluid.msh',
		'species_id': 'tracer', 'initial_concentration_mol_m3': 0.2,
		'diffusivity_m2_s': 0.01, 'source_mol_m3_s': 0.0,
		'first_order_decay_rate_s_inv': 0.1,
		'fluid_velocity_m_s': [0.0, 0.0, 0.0],
		'mesh_velocity_m_s': [0.0, 0.0, 0.0],
		'inflow_concentration_by_label_mol_m3': {
			'0': 1.0, '1': 1.0, '2': 1.0, '7': 1.0},
		'monotone': False, 'time': {'dt_s': 0.01, 'steps': 2}}
	(output/'gpu_species_case.json').write_text(
		json.dumps(gpu_species_case, indent=2)+'\n', encoding='utf-8')
	gpu_coupled = {'schema_version': 'solver-v1',
		'time': {'mode': 'transient', 'dt_s': 0.01, 'steps': 2},
		'domains': [
			{'id': 'fluid', 'dimension': '3d', 'method': 'FEM',
				'device': 'GPU', 'backend': 'native_tet_flow_cuda',
				'physics': ['flow'], 'case_file': 'gpu_flow_case.json',
				'ports': {'velocity': {}}},
			{'id': 'tracer', 'dimension': '3d', 'method': 'FEM',
				'device': 'GPU', 'backend': 'native_tet_species_cuda',
				'physics': ['species'], 'case_file': 'gpu_species_case.json',
				'ports': {'advection': {}}}],
		'connections': [{'from': 'fluid:velocity', 'to': 'tracer:advection',
			'exchange': ['velocity'], 'scheme': 'sequential'}],
		'resources': {'mpi_ranks': 1},
		'output': {'directory': 'gpu_flow_species_results', 'every_steps': 1}}
	(output/'gpu_flow_species_solver.json').write_text(
		json.dumps(gpu_coupled, indent=2)+'\n', encoding='utf-8')
	gpu_oxygen_case = {**gpu_species_case, 'species_id': 'oxygen',
		'initial_concentration_mol_m3': 0.4,
		'diffusivity_m2_s': 0.002,
		'first_order_decay_rate_s_inv': 0.05,
		'inflow_concentration_by_label_mol_m3': {
			'0': 0.8, '1': 0.8, '2': 0.8, '7': 0.8}}
	(output/'gpu_oxygen_case.json').write_text(
		json.dumps(gpu_oxygen_case, indent=2)+'\n', encoding='utf-8')
	gpu_multi = {**gpu_coupled,
		'domains': gpu_coupled['domains'] + [{
			'id': 'oxygen', 'dimension': '3d', 'method': 'FEM',
			'device': 'GPU', 'backend': 'native_tet_species_cuda',
			'physics': ['species'], 'case_file': 'gpu_oxygen_case.json',
			'ports': {'advection': {}}}],
		'connections': gpu_coupled['connections'] + [{
			'from': 'fluid:velocity', 'to': 'oxygen:advection',
			'exchange': ['velocity'], 'scheme': 'sequential'}],
		'output': {'directory': 'gpu_flow_multispecies_results', 'every_steps': 1}}
	(output/'gpu_flow_multispecies_solver.json').write_text(
		json.dumps(gpu_multi, indent=2)+'\n', encoding='utf-8')
	print(f'native FSI channel fixture: {len(fluid_cells)} fluid tetrahedra, '
		f'{len(wall_cells)} wall tetrahedra, case={output/"case.json"}')


if __name__ == '__main__':
	main()

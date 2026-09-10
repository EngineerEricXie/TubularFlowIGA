#!/usr/bin/env python3
"""Create the translating-cube graph regression without running a simulation."""
import argparse
import json
from pathlib import Path
import sys

from hpc_immersed_graph_regression import fixture
from hpc_inventory import digest


def create_case(repo, output, steps=4):
    if steps not in (2, 4):
        raise ValueError('steps must be 2 or 4')
    output.mkdir(parents=True, exist_ok=False)
    case = fixture(repo, output, 'aitken', transient=True, wall_inertial_gamma0=1.0)
    for relative in ('simulation_config.json', 'source/simulation_config.json',
                     'sink/simulation_config.json', 'immersed/simulation_config.json'):
        path = case/relative
        config = json.loads(path.read_text())
        config['time'].update(dt=.25/steps, steps=steps)
        for equation in config.get('equation_systems', []):
            equation['density'] = 1
            viscosity = 'viscosity' if relative.startswith('immersed/') else 'dynamic_viscosity'
            equation[viscosity] = .1
        path.write_text(json.dumps(config, indent=2)+'\n')

    path = case/'immersed/immersed_geometry.json'
    geometry = json.loads(path.read_text())
    geometry['grid'] = dict(lower_m=[0, 0, 0], upper_m=[2.4, 1.2, 1.2], cells=[8, 3, 3])
    # The validated moving fixture uses the runtime's default wall penalty.
    geometry['runtime'].pop('wall_gamma0', None)
    geometry['runtime'].update(flow_controller_reference_flow_m3_s=1,
                               minimum_damping=1/8192, lu_pivot_shift=0)
    geometry['prescribed_motion'] = dict(
        frames=[dict(time_s=0, surface='surface.vtp'), dict(time_s=.25, surface='motion.vtp')],
        extension_layers=2,
        conservation_limits=dict(divergence_theorem=.01, reynolds=1e-8, moving_mass=.01,
                                 wall_relative_leakage=.01, discrete_continuity=1e-8),
        volume_fitting=dict(support_expansion=3))
    path.write_text(json.dumps(geometry, indent=2)+'\n')

    # Retain the existing cube's directed triangles, labels and VTP encoding.
    template = (case/'immersed/surface.vtp').read_text().removeprefix('<?xml version="1.0"?>\n')
    begin = template.index('>', template.index('<DataArray', template.index('<Points>')))+1
    end = template.index('</DataArray>', begin)
    vertices = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for name, shift in [('surface.vtp', 0.), ('motion.vtp', .74)]:
        coordinates = ' '.join(format(value+.1+(shift if axis == 0 else 0.), '.17g')
                               for vertex in vertices for axis, value in enumerate(vertex))
        (case/'immersed'/name).write_text(template[:begin]+coordinates+template[end:])

    manifest = dict(schema_version=1, kind='moving_graph_fixture', steps=steps, dt_s=.25/steps,
                    case_files={str(path.relative_to(case)): digest(path)
                                for path in sorted(case.rglob('*')) if path.is_file()},
                    generator_sources={str(path.relative_to(repo)): digest(path) for path in
                                       (repo/'scripts/hpc_make_moving_graph_case.py',
                                        repo/'scripts/hpc_immersed_graph_regression.py',
                                        repo/'scripts/hpc_inventory.py')})
    (output/'fixture-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path, help='new directory; existing paths are rejected')
    parser.add_argument('--steps', type=int, choices=(2, 4), default=4,
                        help='4: accepted dt=.0625 fixture; 2: documented coarse-step nonconvergence probe')
    args = parser.parse_args()
    try:
        print(create_case(Path(__file__).resolve().parents[1], args.output_dir.resolve(), args.steps))
    except (OSError, ValueError) as error:
        print(f'moving graph fixture: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())

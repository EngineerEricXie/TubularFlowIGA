"""Synthetic regression tests; no solver or simulation data required."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from repair_network_vtkhdf_geometry import convert


class RepairNetworkGeometryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.source = self.directory / 'profile.vtkhdf'
        self.config = self.directory / 'simulation_config.json'
        # Bent inlet, bifurcation, and a closed cycle: endpoints must not shrink.
        self.points = np.array([[0, 0, 0], [1, 1, 0], [2, 1, 0],
                                [3, 2, 0], [3, 0, 0]], dtype=float) * .001
        self.edges = np.array([[1, 2], [0, 1], [2, 4], [3, 4], [2, 3]])
        self.times = np.array([0., .005, .01])
        obj = ['v ' + ' '.join(map(str, point / .001)) for point in self.points]
        obj += [f'l {a+1} {b+1}' for a, b in self.edges]
        (self.directory / 'network.obj').write_text('\n'.join(obj) + '\n')
        self.config.write_text(json.dumps({
            'geometry': {'file': 'network.obj', 'length_scale_to_m': .001},
            'equation_systems': [{'discretization': {'cells_per_segment': 1}}],
        }))

    def fixture(self, connectivity, native=False):
        n = len(connectivity)
        centers = self.points[connectivity].mean(axis=1)
        self.fields = {'pressure': np.arange(3*n, dtype=float) + .25,
                       'oxygen': np.arange(3*n, dtype=float) * .001,
                       'velocity': np.arange(3*n*3, dtype=float).reshape(3*n, 3)}
        with h5py.File(self.source, 'w') as file:
            root = file.create_group('VTKHDF')
            root.attrs['Type'] = np.bytes_('UnstructuredGrid' if native else 'PolyData')
            root.create_dataset('Points', data=centers if native else np.tile(centers, (3, 1)))
            root.create_dataset('NumberOfPoints', data=[n] if native else [n]*3)
            data = root.create_group('PointData')
            steps = root.create_group('Steps')
            steps.attrs['NSteps'] = 3
            steps.create_dataset('Values', data=self.times)
            steps.create_dataset('PointOffsets', data=[0]*3 if native else np.arange(3)*n)
            offsets = steps.create_group('PointDataOffsets')
            for name, values in self.fields.items():
                data.create_dataset(name, data=values)
                offsets.create_dataset(name, data=np.arange(3)*n)
            if native:
                data.create_dataset('parent_node_id', data=connectivity[:, 0]+1)
                data.create_dataset('child_node_id', data=connectivity[:, 1]+1)

    def check_conversion(self, connectivity, native=False):
        self.fixture(connectivity, native)
        original_hash = hashlib.sha256(self.source.read_bytes()).hexdigest()
        with contextlib.redirect_stdout(io.StringIO()):
            target = convert(self.source, self.config)
        self.assertEqual(original_hash, hashlib.sha256(self.source.read_bytes()).hexdigest())
        n = len(connectivity)
        with h5py.File(target, 'r') as file:
            root = file['VTKHDF']
            self.assertEqual(root.attrs['Type'], b'UnstructuredGrid')
            np.testing.assert_array_equal(root['Points'][:], self.points)
            np.testing.assert_array_equal(root['Connectivity'][:], connectivity.ravel())
            np.testing.assert_array_equal(root['Offsets'][:], np.arange(n+1)*2)
            self.assertEqual(root['Offsets'].dtype, root['Connectivity'].dtype)
            np.testing.assert_array_equal(root['Types'][:], np.full(n, 3))
            np.testing.assert_array_equal(root['Steps/Values'][:], self.times)
            self.assertEqual(len(root['PointData']), 0)
            for name, values in self.fields.items():
                np.testing.assert_array_equal(root['CellData'][name][:], values)
                np.testing.assert_array_equal(root['Steps/CellDataOffsets'][name][:],
                                              np.arange(3)*n)
            if native:
                np.testing.assert_array_equal(root['CellData/parent_node_id'][:],
                                              connectivity[:, 0]+1)
                np.testing.assert_array_equal(root['Steps/CellDataOffsets/parent_node_id'][:],
                                              np.zeros(3))
        report = json.loads(target.with_suffix('.validation.json').read_text())
        self.assertEqual(report['timesteps'], 3)
        self.assertEqual(report['lines'], n)
        with self.assertRaisesRegex(RuntimeError, 'overwrite'):
            convert(self.source, self.config)

    def test_obj_order_polydata(self):
        self.check_conversion(self.edges)

    def test_sorted_solver_order_polydata(self):
        self.check_conversion(np.array(sorted(map(tuple, self.edges))))

    def test_permuted_midpoint_lookup(self):
        self.check_conversion(self.edges[[3, 0, 4, 1, 2]])

    def test_native_node_ids_and_static_geometry(self):
        self.check_conversion(self.edges[[2, 4, 1, 3, 0]][:, ::-1], native=True)

    def test_subdivided_segments_rejected(self):
        config = json.loads(self.config.read_text())
        config['equation_systems'][0]['discretization']['cells_per_segment'] = 2
        self.config.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, 'one cell'):
            convert(self.source, self.config)

    def test_wrong_geometry_rejected(self):
        self.fixture(self.edges)
        with h5py.File(self.source, 'a') as file:
            file['VTKHDF/Points'][0] = [100, 100, 100]
        with self.assertRaisesRegex(ValueError, 'Midpoints'):
            convert(self.source, self.config)
        self.assertFalse(self.source.with_name('profile_geometry_corrected.vtkhdf').exists())


if __name__ == '__main__':
    unittest.main()

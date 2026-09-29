"""Restore OBJ edges in midpoint-based 0D/1D VTKHDF, without solving again."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.spatial import cKDTree

def require(condition, message):
    if not condition:
        raise ValueError(message)


def convert(source, config_path, target=None):
    source = Path(source).resolve()
    config_path = Path(config_path).resolve()
    config_dir = config_path.parent
    config = json.loads(config_path.read_text())
    require(config['equation_systems'][0]['discretization']['cells_per_segment'] == 1,
            'Only one cell per OBJ segment is supported')
    points, edges = [], []
    with (config_dir / config['geometry']['file']).open() as stream:
        for line in stream:
            fields = line.split('#', 1)[0].split()
            if fields and fields[0] == 'v':
                points.append([float(x) for x in fields[1:4]])
            elif fields and fields[0] == 'l':
                ids = [int(x.split('/')[0])-1 for x in fields[1:]]
                require(all(i >= 0 for i in ids), 'OBJ indices must be positive')
                edges.extend(zip(ids[:-1], ids[1:]))
    points = np.asarray(points) * config['geometry']['length_scale_to_m']
    edges = np.asarray(edges, dtype=np.int64)
    n = len(edges)
    require(points.ndim == 2 and points.shape[1] == 3 and np.isfinite(points).all(),
            'Expected finite 3D OBJ vertices')
    require(n > 0 and edges.max() < len(points), 'Empty or invalid OBJ connectivity')
    centers = (points[edges[:, 0]] + points[edges[:, 1]]) / 2
    target = (Path(target).resolve() if target is not None else
              source.with_name(source.stem + '_geometry_corrected.vtkhdf'))
    if target.exists() or target.with_suffix('.validation.json').exists():
        raise RuntimeError(f'Refusing to overwrite {target}')
    with h5py.File(source, 'r') as src:
        old = src['VTKHDF']
        steps = old['Steps']
        count = int(steps.attrs['NSteps'])
        times = steps['Values'][:count]
        require(count > 0 and len(times) == count and np.isfinite(times).all()
                and np.all(np.diff(times) > 0), 'Invalid time axis')
        old_points = old['Points'][:n]
        data = old['PointData']
        if 'parent_node_id' in data:
            conn = np.column_stack((data['parent_node_id'][:]-1, data['child_node_id'][:]-1))
            require(conn.shape == edges.shape, 'Node ID count does not match OBJ edges')
            require(sorted(map(tuple, np.sort(conn, axis=1))) ==
                    sorted(map(tuple, np.sort(edges, axis=1))), 'Node IDs do not match OBJ edges')
        elif np.allclose(old_points, centers, rtol=0, atol=1e-12):
            conn = edges
        else:
            ordered = np.sort(edges, axis=1)
            ordered = ordered[np.lexsort((ordered[:, 1], ordered[:, 0]))]
            if np.allclose(old_points, (points[ordered[:, 0]]+points[ordered[:, 1]])/2, rtol=0, atol=1e-12):
                conn = ordered
            else:
                distances, indices = cKDTree(centers).query(old_points, k=2)
                require(np.max(distances[:, 0]) < 1e-12, 'Midpoints do not match OBJ edges')
                require(np.min(distances[:, 1]) > 1e-12, 'Ambiguous midpoint mapping')
                require(len(np.unique(indices[:, 0])) == n, 'Non-bijective midpoint mapping')
                conn = edges[indices[:, 0]]
        conn = np.asarray(conn, dtype=np.int64)
        require(np.allclose(old_points, (points[conn[:, 0]]+points[conn[:, 1]])/2,
                            rtol=0, atol=1e-12), 'Connectivity does not match field locations')
        dynamic = list(steps['PointDataOffsets'])
        with h5py.File(target, 'x') as dst:
            root = dst.create_group('VTKHDF')
            root.attrs['Type'] = np.bytes_('UnstructuredGrid')
            root.attrs['Version'] = np.array([2, 1], dtype=np.int64)
            for name, value in {
                'Points': points, 'Connectivity': conn.ravel(),
                'Offsets': np.arange(n+1, dtype=np.int64)*2,
                'Types': np.full(n, 3, dtype=np.uint8),
                'NumberOfPoints': np.array([len(points)], dtype=np.int64),
                'NumberOfCells': np.array([n], dtype=np.int64),
                'NumberOfConnectivityIds': np.array([2*n], dtype=np.int64),
            }.items():
                root.create_dataset(name, data=value)
            root.create_group('PointData')
            cells = root.create_group('CellData')
            newsteps = root.create_group('Steps')
            newsteps.attrs['NSteps'] = np.int64(count)
            newsteps.create_dataset('Values', data=times)
            for name in ('PartOffsets', 'PointOffsets'):
                newsteps.create_dataset(name, data=np.zeros(count, dtype=np.int64))
            newsteps.create_dataset('NumberOfParts', data=np.ones(count, dtype=np.int64))
            for name in ('CellOffsets', 'ConnectivityIdOffsets'):
                newsteps.create_dataset(name, data=np.zeros((count, 1), dtype=np.int64))
            newsteps.create_group('PointDataOffsets')
            offsets = newsteps.create_group('CellDataOffsets')
            for name, field in data.items():
                if name not in dynamic:
                    require(len(field) == n, f'Invalid static field size: {name}')
                    cells.create_dataset(name, data=field[:])
                    offsets.create_dataset(name, data=np.zeros(count, dtype=np.int64))
                    continue
                output = cells.create_dataset(name, shape=(count*n,)+field.shape[1:], dtype=field.dtype)
                offsets.create_dataset(name, data=np.arange(count, dtype=np.int64)*n)
                for index in range(count):
                    start = int(steps['PointDataOffsets'][name][index])
                    values = field[start:start+n]
                    require(len(values) == n and np.isfinite(values).all(),
                            f'Invalid temporal field: {name}, step {index}')
                    output[index*n:(index+1)*n] = values
                print(f'{source.parent}: copied {name} ({count} states)', flush=True)
            # Confirm geometry stays fixed at every source time, not just t=0.
            for index in range(count):
                start = int(steps['PointOffsets'][index])
                require(np.array_equal(old['Points'][start:start+n], old_points),
                        'Moving geometry is not supported')
        with h5py.File(target, 'r') as dst:
            root = dst['VTKHDF']
            require(np.array_equal(root['Points'][:], points), 'Output point mismatch')
            require(np.array_equal(root['Connectivity'][:], conn.ravel()), 'Output edge mismatch')
            require(np.array_equal(root['Steps/Values'][:], times), 'Output time mismatch')
            for name in dynamic:
                for index in range(count):
                    start = int(steps['PointDataOffsets'][name][index])
                    require(np.array_equal(data[name][start:start+n],
                                           root['CellData'][name][index*n:(index+1)*n]),
                            f'Output field mismatch: {name}, step {index}')
        report = dict(source=str(source), output=str(target), timesteps=count,
                      first_time=float(times[0]), last_time=float(times[-1]),
                      points=len(points), lines=n, arrays=dynamic,
                      validation='Original OBJ connectivity restored; all temporal values exactly equal')
        with target.with_suffix('.validation.json').open('x') as stream:
            stream.write(json.dumps(report, indent=2)+'\n')
        print(json.dumps(report), flush=True)
        return target

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='Original midpoint-based VTKHDF')
    parser.add_argument('--config', type=Path, required=True,
                        help='simulation_config.json used for the original run')
    parser.add_argument('--output', type=Path, help='New VTKHDF path; never overwrites')
    args = parser.parse_args()
    convert(args.source, args.config, args.output)

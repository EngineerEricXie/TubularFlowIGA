from pathlib import Path
import tempfile

from paraview.simple import OpenDataFile


directory = Path(tempfile.gettempdir()) / "tubularflowiga-one-d-output-test"
path = directory / "vtkhdf/profile_1d.vtkhdf"
reader = OpenDataFile(str(path))
if reader is None:
    raise RuntimeError(f"ParaView could not open {path}")
reader.UpdatePipeline()
if list(reader.TimestepValues) != [0.0, 0.1]:
    raise RuntimeError(f"unexpected timesteps: {list(reader.TimestepValues)}")
reader.UpdatePipeline(time=0.1)
data = reader.GetClientSideObject().GetOutputDataObject(0)
if data.GetNumberOfPoints() != 5 or data.GetNumberOfCells() != 4:
    raise RuntimeError("unexpected 1D VTKHDF geometry dimensions")
if {data.GetCellType(cell) for cell in range(data.GetNumberOfCells())} != {3}:
    raise RuntimeError("1D VTKHDF cells are not VTK_LINE cells")
point_data = data.GetCellData()
for name in ("area", "flow_rate", "pressure", "velocity", "segment_id",
             "parent_node_id", "child_node_id", "segment_cell",
             "axial_position_m", "reference_radius_m"):
    if point_data.GetArray(name) is None:
        raise RuntimeError(f"missing cell array {name}")
if point_data.GetArray("pressure").GetValue(0) != 110.0:
    raise RuntimeError("ParaView did not select the final pressure timestep")
print(f"ParaView validated 1D temporal VTKHDF: {path}")

for subdivisions in (1, 3):
    for cycle in (False, True):
        stem = f"{'cycle' if cycle else 'branch'}_{subdivisions}"
        edges = 4 if cycle else 3
        for extension in ('vtkhdf', 'vtp'):
            reader = OpenDataFile(str(directory / f'{stem}.{extension}'))
            times = list(reader.TimestepValues) if extension == 'vtkhdf' else [0.0]
            if extension == 'vtkhdf' and times != [0.0, 0.1, 0.2]:
                raise RuntimeError(f'Unexpected resumed times for {stem}: {times}')
            for time in times:
                reader.UpdatePipeline(time=time)
                data = reader.GetClientSideObject().GetOutputDataObject(0)
                assert data.GetNumberOfPoints() == 4 + edges * (subdivisions - 1)
                assert data.GetNumberOfCells() == edges * subdivisions
                assert data.GetPoint(0) == (0.0, 0.0, 0.0)
                assert data.GetPoint(1) == (0.01, 0.01, 0.0)
                degree = [0] * data.GetNumberOfPoints()
                for cell in range(data.GetNumberOfCells()):
                    line = data.GetCell(cell)
                    assert line.GetCellType() == 3 and line.GetNumberOfPoints() == 2
                    for index in range(2):
                        degree[line.GetPointId(index)] += 1
                    resumed = 10.0 if extension == 'vtkhdf' and time > 0.0 else 0.0
                    assert data.GetCellData().GetArray('pressure').GetValue(cell) == 100.0 + cell + resumed
                    assert data.GetCellData().GetArray('oxygen').GetValue(cell) == 0.01 * cell
                assert degree[:4] == [1, 3, 2 if cycle else 1, 2 if cycle else 1]
                assert all(value == 2 for value in degree[4:])
                assert data.GetPointData().GetArray('pressure') is None
                if extension == 'vtkhdf':
                    for name in ('segment_id', 'parent_node_id', 'child_node_id',
                                 'segment_cell', 'axial_position_m', 'reference_radius_m'):
                        assert data.GetCellData().GetArray(name).GetNumberOfTuples() == edges * subdivisions
            print(f'ParaView validated native topology and cell values: {stem}.{extension}')

reader = OpenDataFile(str(directory / 'zero_d/profile_0d.vtkhdf'))
reader.UpdatePipeline(time=0.1)
data = reader.GetClientSideObject().GetOutputDataObject(0)
assert data.GetNumberOfPoints() == 5 and data.GetNumberOfCells() == 4
assert data.GetCellData().GetArray('pressure').GetValue(0) == 110.0
print('ParaView validated shared 0D native writer')

reader = OpenDataFile(str(directory / 'point_compatibility.vtkhdf'))
reader.UpdatePipeline(time=0.1)
data = reader.GetClientSideObject().GetOutputDataObject(0)
assert data.GetPointData().GetArray('pressure').GetValue(1) == 40.0
print('ParaView validated backward-compatible point writer')

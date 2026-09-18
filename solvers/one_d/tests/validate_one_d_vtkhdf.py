from pathlib import Path

from paraview.simple import OpenDataFile


path = Path("/tmp/tubularflowiga-one-d-output-test/vtkhdf/profile_1d.vtkhdf")
reader = OpenDataFile(str(path))
if reader is None:
    raise RuntimeError(f"ParaView could not open {path}")
reader.UpdatePipeline()
if list(reader.TimestepValues) != [0.0, 0.1]:
    raise RuntimeError(f"unexpected timesteps: {list(reader.TimestepValues)}")
reader.UpdatePipeline(time=0.1)
data = reader.GetClientSideObject().GetOutputDataObject(0)
if data.GetNumberOfPoints() != 4 or data.GetNumberOfCells() != 3:
    raise RuntimeError("unexpected 1D VTKHDF geometry dimensions")
if {data.GetCellType(cell) for cell in range(data.GetNumberOfCells())} != {3}:
    raise RuntimeError("1D VTKHDF cells are not VTK_LINE cells")
point_data = data.GetPointData()
for name in ("area", "flow_rate", "pressure", "velocity", "segment_id",
             "parent_node_id", "child_node_id", "segment_cell",
             "axial_position_m", "reference_radius_m"):
    if point_data.GetArray(name) is None:
        raise RuntimeError(f"missing point array {name}")
if point_data.GetArray("pressure").GetValue(0) != 110.0:
    raise RuntimeError("ParaView did not select the final pressure timestep")
print(f"ParaView validated 1D temporal VTKHDF: {path}")

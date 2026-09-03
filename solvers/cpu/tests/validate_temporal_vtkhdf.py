import argparse
from pathlib import Path

from paraview.simple import OpenDataFile


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


parser = argparse.ArgumentParser()
parser.add_argument("path", nargs="?", default=
                    "/tmp/tubularflowiga-temporal-vtkhdf-test/bezier.vtkhdf")
parser.add_argument("--points", type=int, default=64)
parser.add_argument("--cells", type=int, default=1)
parser.add_argument("--times", default="0,0.5,1")
parser.add_argument("--arrays", default="scalar,velocity")
parser.add_argument("--value", default="scalar:20")
args = parser.parse_args()

path = Path(args.path)
reader = OpenDataFile(str(path))
require(reader is not None, f"ParaView could not open {path}")
reader.UpdatePipeline()
times = list(reader.TimestepValues)
expected_times = [] if args.times == "none" else [
    float(value) for value in args.times.split(",") if value]
require(times == expected_times, f"unexpected timesteps: {times}")
if expected_times:
    reader.UpdatePipeline(time=expected_times[-1])
data = reader.GetClientSideObject().GetOutputDataObject(0)
require(data.GetNumberOfPoints() == args.points, "geometry point count changed")
require(data.GetNumberOfCells() == args.cells, "geometry cell count changed")
require(data.GetCellType(0) == 79, "cell is not VTK_BEZIER_HEXAHEDRON")
point_data = data.GetPointData()
require(point_data.GetArray("boundary_label") is not None,
        "missing static point boundary_label array")
cell_data = data.GetCellData()
require(cell_data.GetArray("boundary_label") is not None,
        "missing static cell boundary_label array")
face_labels = cell_data.GetArray("boundary_face_labels")
require(face_labels is not None, "missing boundary_face_labels array")
require(face_labels.GetNumberOfComponents() == 6,
        "boundary_face_labels does not have six components")
array_names = [] if args.arrays == "none" else [
    value for value in args.arrays.split(",") if value]
for name in array_names:
    require(point_data.GetArray(name) is not None, f"missing {name} array")
if args.value and args.value != "none":
    name, expected = args.value.split(":", 1)
    require(point_data.GetArray(name).GetValue(0) == float(expected),
            "time-dependent point data was not selected")
print(f"ParaView validated temporal VTKHDF: {path}")

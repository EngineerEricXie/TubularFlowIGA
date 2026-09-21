#!/usr/bin/env python3
"""Exercise fTetWild normalization, native-reader compatibility, and failures."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ftetwild_to_fem_volume import oriented_boundary


FAKE_FTETWILD = r'''#!/usr/bin/env python3
from pathlib import Path
import sys

arguments = sys.argv[1:]
def value(name): return arguments[arguments.index(name)+1]
vertices = []
for line in Path(value("--input")).read_text().splitlines():
	if line.startswith("v "):
		vertices.append(tuple(map(float, line.split()[1:4])))
if len(vertices) != 4:
	raise SystemExit("fake expects one tetrahedral surface")
output = Path(value("--output"))
with output.open("w") as stream:
	stream.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n4\n")
	for index, point in enumerate(vertices, 1):
		stream.write(f"{index} {point[0]} {point[1]} {point[2]}\n")
	stream.write("$EndNodes\n$Elements\n1\n1 4 0 1 2 3 4\n$EndElements\n")
'''


SURFACE = '''<?xml version="1.0"?>
<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian"><PolyData>
<Piece NumberOfPoints="4" NumberOfPolys="4">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
0 0 0 1 0 0 0 1 0 0 0 1
</DataArray></Points><Polys>
<DataArray type="Int32" Name="connectivity" format="ascii">0 2 1 0 1 3 0 3 2 1 2 3</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">3 6 9 12</DataArray>
</Polys><CellData><DataArray type="Int32" Name="boundary_id" format="ascii">0 1 2 0</DataArray>
</CellData></Piece></PolyData></VTKFile>
'''


def run(command):
	return subprocess.run(command, text=True, capture_output=True)


def main():
	root = Path(__file__).resolve().parents[2]
	nodes = {1: (0,0,0), 2: (1,0,0), 3: (0,1,0), 4: (0,0,1),
		5: (3,0,0), 6: (4,0,0), 7: (3,1,0), 8: (3,0,1)}
	for cells, expected in (([(1,2,3,4), (1,2,3,4)], "duplicate"),
			([(1,2,3,4), (5,6,7,8)], "disconnected")):
		try:
			oriented_boundary(nodes, cells)
		except RuntimeError as error:
			assert expected in str(error)
		else:
			raise AssertionError(f"{expected} tetrahedral output was accepted")
	adapter = root/"scripts"/"ftetwild_to_fem_volume.py"
	gmsh_adapter = root/"scripts"/"surface_to_fem_volume.py"
	preflight = root/"solvers"/"cpu"/"surface_fem_preflight"
	native_test = root/"solvers"/"cpu"/"native_tet_fem_test"
	with tempfile.TemporaryDirectory(prefix="tubularflow-ftetwild-test-") as directory:
		temporary = Path(directory)
		surface = temporary/"surface.vtp"
		surface.write_text(SURFACE, encoding="utf-8")
		fake = temporary/"FloatTetwild_bin"
		fake.write_text(FAKE_FTETWILD, encoding="utf-8"); fake.chmod(0o755)
		mesh, manifest = temporary/"volume.msh", temporary/"volume.json"
		command = [sys.executable, str(adapter), str(surface), str(mesh),
			"--manifest", str(manifest), "--ftetwild", str(fake),
			"--surface-preflight", str(preflight), "--target-size-m", "0.5",
			"--envelope-m", "1e-8", "--ftetwild-version", "test-double"]
		completed = run(command)
		assert completed.returncode == 0, completed.stdout+completed.stderr
		data = json.loads(manifest.read_text(encoding="utf-8"))
		assert data["route"] == "surface_to_fem_volume_ftetwild"
		assert data["mesher"]["name"] == "fTetWild"
		assert data["mesher"]["max_optimization_passes"] == 80
		assert data["resources"]["pipeline_wall_s"] > 0.0
		assert data["resources"]["peak_rss_bytes"] > 0
		assert data["volume_mesh"]["boundary_labels"] == {"0": 2, "1": 1, "2": 1}
		assert data["volume_mesh"]["boundary_matches_tetrahedra"] is True
		assert data["volume_mesh"]["maximum_source_surface_distance_m"] == 0.0
		assert data["volume_mesh"]["boundary_surface_area_relative_error"] < 1e-14
		assert data["volume_mesh"]["boundary_remeshed"] is True
		assert data["volume_mesh"]["implicit_surface_repair"] is False
		contract = data["geometry_contract"]
		assert contract["route"] == "surface_to_fem_volume_ftetwild"
		assert contract["source_geometry"]["identity_sha256"] == \
			data["surface_preflight"]["canonical_sha256"]
		assert contract["reference_geometry"]["identity_sha256"] == \
			data["volume_mesh"]["sha256"]
		assert contract["current_geometry"]["identity_sha256"] == \
			data["volume_mesh"]["sha256"]
		assert contract["stable_ids"]["volume_nodes"] == "output_msh_node_tag"
		change = data["geometry_change"]
		assert change["classification"] == "envelope_constrained_boundary_remesh_not_repair"
		assert change["implicit_repair"] is False
		assert change["maximum_source_surface_distance_m"] == 0.0
		assert change["output_shape_metrics"]["integrated_absolute_mean_curvature_m"] > 0.0
		parsed = run([str(native_test), str(mesh)])
		assert parsed.returncode == 0, parsed.stdout+parsed.stderr
		assert "external mesh passed" in parsed.stdout

		missing = run(command[:command.index("--ftetwild")]+[
			"--ftetwild", str(temporary/"missing"), *command[command.index("--surface-preflight"):]])
		assert missing.returncode == 2 and "missing fTetWild executable" in missing.stdout

		ambiguous = run(command+["--label-ambiguity-tolerance-m", "2"])
		assert ambiguous.returncode == 2 and "ambiguous" in ambiguous.stdout

		open_surface = temporary/"open-surface.vtp"
		open_surface.write_text(SURFACE.replace('NumberOfPolys="4"', 'NumberOfPolys="3"')
			.replace("0 2 1 0 1 3 0 3 2 1 2 3", "0 2 1 0 1 3 0 3 2")
			.replace("3 6 9 12", "3 6 9").replace(">0 1 2 0</", ">0 1 2</"),
			encoding="utf-8")
		gmsh_failure = run([sys.executable, str(gmsh_adapter), str(open_surface),
			str(temporary/"open-gmsh.msh"), "--surface-preflight", str(preflight),
			"--target-size-m", "0.5"])
		ftetwild_failure = run([sys.executable, str(adapter), str(open_surface),
			str(temporary/"open-ftetwild.msh"), "--ftetwild", str(fake),
			"--surface-preflight", str(preflight), "--target-size-m", "0.5",
			"--envelope-m", "1e-8"])
		for failure in (gmsh_failure, ftetwild_failure):
			assert failure.returncode == 2
			assert "surface edge is not used by exactly two faces" in failure.stdout+failure.stderr, \
				failure.stdout+failure.stderr
	print("ftetwild_to_fem_volume_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())

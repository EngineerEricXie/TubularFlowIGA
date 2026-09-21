#!/usr/bin/env python3
"""Verify that the body-fitted geometry identity is partition invariant."""

import json
from pathlib import Path
import struct
import subprocess
import tempfile


def write_database(path, ranks):
	header_size = 8+4+4+8+8+4+4+8+5*8
	index_size = 2*8+4
	record_offset = header_size+index_size
	connectivity = [0]
	record = bytearray()
	record += struct.pack("<QiiI", 17, 1, ranks-1, len(connectivity))
	record += struct.pack("<6i", 4, -1, -1, -1, -1, -1)
	record += struct.pack("<i", *connectivity)
	record += b"\0"
	record += struct.pack("<192d", *([0.0]*192))
	rank_index_offset = record_offset+len(record)
	contents = bytearray()
	contents += b"NTIGADB2"
	contents += struct.pack("<IIQQIIQ", 5, ranks, 1, 8, 64, 0, rank_index_offset)
	contents += struct.pack("<5d", 1.0, 2.0, 3.0, 4.0, 0.001)
	contents += struct.pack("<2Q", record_offset, rank_index_offset)
	contents += struct.pack("<i", ranks-1)
	contents += record
	contents += struct.pack(f"<{ranks+1}Q", *([0]*(ranks+1)))
	path.write_bytes(contents)


def inspect(executable, database, manifest, role="fluid"):
	return subprocess.run([str(executable), str(database), "--geometry-manifest",
		str(manifest), "--region-role", role], text=True, capture_output=True)


def main():
	root = Path(__file__).resolve().parents[2]
	executable = root/"solvers"/"cpu"/"iga_inspect"
	with tempfile.TemporaryDirectory(prefix="tubularflow-iga-manifest-") as directory:
		temporary = Path(directory)
		databases = [temporary/"rank2.ntiga", temporary/"rank3.ntiga"]
		manifests = [temporary/"rank2.json", temporary/"rank3.json"]
		for ranks, database, manifest in zip((2, 3), databases, manifests):
			write_database(database, ranks)
			completed = inspect(executable, database, manifest)
			assert completed.returncode == 0, completed.stderr
		data = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
		assert data[0]["reference_geometry"]["identity_sha256"] \
			== data[1]["reference_geometry"]["identity_sha256"]
		assert data[0]["reference_geometry"]["source_artifact_sha256"] \
			!= data[1]["reference_geometry"]["source_artifact_sha256"]
		for item in data:
			assert item["route"] == "centerline_to_iga_volume"
			assert item["current_geometry"]["identity_sha256"] \
				== item["reference_geometry"]["identity_sha256"]
			assert item["regions"] == [{"id": "volume", "role": "fluid", "dimension": 3}]
		bad = inspect(executable, databases[0], temporary/"bad.json", "organ")
		assert bad.returncode != 0
		assert not (temporary/"bad.json").exists()
	print("iga_geometry_manifest_test: PASS")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())

#!/usr/bin/env python3
"""Compare matching text fields without loading the full field into memory."""

import argparse
import itertools
import json
import math
from pathlib import Path
import sys

from hpc_inventory import digest


def compare(reference, candidate, relative_tolerance, absolute_tolerance, node_ids=False):
    if any(not math.isfinite(value) or value < 0 for value in (relative_tolerance, absolute_tolerance)):
        raise ValueError("tolerances must be finite and nonnegative")
    reference_norm = defect_norm = 0.0
    columns = None
    rows = 0
    previous_id = -1
    with reference.open() as first, candidate.open() as second:
        for line_number, pair in enumerate(itertools.zip_longest(first, second), 1):
            if None in pair:
                raise ValueError("field row counts differ")
            left, right = [line.split() for line in pair]
            if not left or len(left) != len(right):
                raise ValueError(f"empty or mismatched row {line_number}")
            if node_ids:
                first_id, second_id = int(left.pop(0)), int(right.pop(0))
                if first_id != second_id or first_id <= previous_id:
                    raise ValueError("node IDs must match and be strictly increasing")
                previous_id = first_id
            if not left:
                raise ValueError("field contains no values")
            if columns is None:
                columns = len(left)
            if len(left) != columns:
                raise ValueError("field width changes between rows")
            left, right = [[float(value) for value in values] for values in (left, right)]
            if not all(math.isfinite(value) for value in left + right):
                raise ValueError("field contains nonfinite values")
            reference_norm = math.hypot(reference_norm, *left)
            defect_norm = math.hypot(defect_norm, *(b-a for a, b in zip(left, right)))
            rows += 1
    if not rows or not all(math.isfinite(value) for value in (reference_norm, defect_norm)):
        raise ValueError("empty field or unrepresentable norm")
    relative_l2 = defect_norm / reference_norm if reference_norm else None
    passed = relative_l2 <= relative_tolerance if relative_l2 is not None else defect_norm <= absolute_tolerance
    return {"schema_version": 1, "kind": "hpc_field_comparison", "passed": passed,
            "rows": rows, "columns": columns, "node_ids": node_ids,
            "reference_sha256": digest(reference), "candidate_sha256": digest(candidate),
            "reference_l2": reference_norm, "absolute_l2": defect_norm,
            "relative_l2": relative_l2, "relative_tolerance": relative_tolerance,
            "zero_reference_absolute_tolerance": absolute_tolerance,
            "scope": "Nodal coefficient comparison; matching units, ordering, and pressure gauge required. Independent physical gates still required."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--rtol", type=float, required=True)
    parser.add_argument("--zero-atol", type=float, required=True)
    parser.add_argument("--node-ids", action="store_true")
    options = parser.parse_args()
    result = compare(options.reference, options.candidate, options.rtol, options.zero_atol, options.node_ids)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"hpc_compare_fields: {error}", file=sys.stderr)
        sys.exit(1)

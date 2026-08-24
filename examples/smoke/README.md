# Public smoke case

This three-point straight tube is intentionally small. It exercises mesh
generation, Bezier extraction, METIS partitioning, database packing, and the
boundary-condition validator without requiring PETSc. From the repository
root:

```bash
./scripts/prepare_smoke.sh
```

The script prints its temporary work directory. Keep it and run the optional
MPI/PETSc commands from [the quick start](../../docs/QUICKSTART.md) on a CPU
compute node.

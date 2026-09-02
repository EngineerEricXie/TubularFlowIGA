# Large NMO_06840 neuron transport

This source-only case transports the configured `N0` and `Nplus` fields through
a large 3D neurite bifurcation derived from NeuroMorpho.Org reconstruction
`NMO_06840`. The 13-node source subtree generates a 41,097-control-point,
36,540-element IGA mesh with one inlet, two outlets, and 8,660 external faces.

The same geometry was used to diagnose a Navier--Stokes solver regression, but
that does not make it a vascular example. The public configuration contains
only the `neuron_transport` system and uses the generated prescribed velocity
field; it does not run a blood-flow solve.

## Resource warning

This is not a quick-start example. Preparing the case creates a sparse spline
cache and packed database of about 363 MiB and 365 MiB, respectively. Allow at
least 1.5 GiB of free disk space, several GiB of RAM, and several minutes for
preprocessing, packing, solving, and validation. A first build or slower
laptop can take tens of minutes.

Use `neuron_transport/straight_neurite` first if the goal is only to check an
installation. Run this large case on a local workstation or an allocated
scheduler compute node, never on a shared cluster login node.

## Prepare and run with 12 local MPI ranks

Build the PETSc solver with the same MPI implementation that will launch it.
From the repository root:

```bash
./scripts/check_dependencies.sh preprocessing
./scripts/check_dependencies.sh cpu
make cpu-petsc PETSC_DIR="$PETSC_DIR" PETSC_ARCH="${PETSC_ARCH:-}"

NMO_SOURCE=examples/neuron_transport/nmo_06840_bifurcation
OMP_NUM_THREADS=12 ./scripts/generate_case.sh "$NMO_SOURCE" --ranks 12
NMO_WORK="$NMO_SOURCE/generated"
NMO_CASE="$NMO_WORK/preprocessing"
NMO_DB="$NMO_WORK/database/nmo_06840_bifurcation-12.ntiga"
NMO_RESULTS="$NMO_WORK/results"

mpiexec -np 12 ./solvers/cpu/iga_mesh_check "$NMO_DB"
mpiexec -np 12 ./solvers/cpu/iga_solve \
  "$NMO_DB" "$NMO_CASE" --system neuron_transport \
  --output "$NMO_RESULTS/neuron-cpu.txt" --output-every 1
```

The CPU MPI process count must equal `RANKS` used during preparation. Do not
mix an OpenMPI launcher with a PETSc executable built against MPICH. If the
launcher exposes fewer than 12 slots, choose a smaller `RANKS` value during
preparation and use it for both MPI commands. On a cluster, follow the local
scheduler policy and the [Bridges-2 guide](../../../docs/BRIDGES2.md).

The result has one row per database node: `node_id N0 Nplus`. The neighboring
`.fields` file records field order. With `--output-every 1`, the solver also
writes the initialized step, both solved VTU snapshots, and a PVD collection
for ParaView.

## CUDA command

A database packed for CPU MPI can also be read by the single-GPU backend,
which ignores ownership records:

```bash
make cuda CUDA_ARCHS=YOUR_GPU_ARCH
./solvers/cuda/iga_cuda device-info
./solvers/cuda/iga_cuda mesh-check "$NMO_DB"
./solvers/cuda/iga_cuda solve \
  "$NMO_DB" "$NMO_CASE" --system neuron_transport \
  --output "$NMO_RESULTS/neuron-cuda.txt" --output-every 1
```

Do not treat a successful CUDA build as a runtime result. Run `device-info`,
the packed geometry check, and the complete transport solve on a CUDA-capable
host. Consumer GPUs with weak FP64 throughput may be much slower than
scientific GPUs.

## Current geometry acceptance evidence

The schema-v4 source inputs reproduce 41,097 nodes and 36,540 elements. They
contain 8,120 wall faces and 180 faces on each of the inlet and two outlet
caps. Control-mesh generation reports `min_scaled_J=0.288586`, zero invalid
elements, and zero boundary-surface intersections. Two-rank packed geometry
validation reports `minimum_detJ=7.24399e-11`, zero bad elements, and zero bad
quadrature samples. Geometry preflight also reports zero hard errors.

## Historical two-step transport evidence

The pre-schema-v4 29,238-node, 25,920-element mesh was run on 2026-08-26 under
WSL on an Intel Core i9-14900KF with PETSc 3.15.5, OpenMPI 4.1.2, and 12 MPI
ranks. These solver values are historical and must not be presented as results
from the current adaptive mesh:

| Measurement | Result |
|---|---:|
| Transport assembly | 26.6397 s |
| Two-step time-loop solve | 109.36 s |
| Total Krylov iterations | 10,520 |
| Final coefficient L2 | `330.035` |
| Output rows and fields | 29,238; `N0`, `Nplus` |
| Final `N0` range | `0.0403333` to `1.0001424` |
| Final `Nplus` range | `1.7740167` to `2.0012732` |

Every output row contained exactly the node ID and two finite field values.
The PVD collection contained the initialized state and both solved states. The
short horizon is an execution regression, not a claim that the transported
field has reached a biological steady state. Increase `time.steps` only when
the additional runtime and output volume are intentional. A CUDA runtime result
has not yet been recorded for this public input.

## Source, attribution, and scientific scope

The committed SWC is a deterministic subtree of NeuroMorpho.Org `NMO_06840`,
[Stevens archive](https://neuromorpho.org/NeuroMorpho_ArchiveLinkout.jsp?ARCHIVE=Stevens&DATE=2011-11-08)
neuron `20061101z174r1c1`, accessed 2026-08-26. It selects the original
node-100-to-node-108 trunk and node-108 descendants ending at leaves 110 and
112. The original SWC SHA-256 is
`45fe0342082bb667626b89b82e0345c315be2761434335e1f1c9d6f409b3cbeb`;
the committed derived SWC SHA-256 is
`c9e2b1258f75d1cc359bfbfb4dd08a2fa0a5d9d1a3ce513babaee1b2cbbb309b`.

NeuroMorpho.Org material is available under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Follow the
[NeuroMorpho.Org Terms of Use](https://neuromorpho.org/useterm.jsp), which
require attribution to the reconstruction's original publication,
NeuroMorpho.Org (RRID:SCR_002145), and the repository publication:

1. Lee S, Stevens CF. General design principle for scalable neural circuits in
   a vertebrate retina. *PNAS*. 2007;104(31):12931--12935.
   [doi:10.1073/pnas.0705469104](https://doi.org/10.1073/pnas.0705469104).
2. NeuroMorpho.Org (RRID:SCR_002145), reconstruction `NMO_06840`, Stevens
   archive, neuron `20061101z174r1c1`.
3. Tecuatl C, Ljungquist B, Ascoli GA. Accelerating the continuous community
   sharing of digital neuromorphology data. *FASEB BioAdvances*.
   2024;6(7):207--221.
   [doi:10.1096/fba.2024-00048](https://doi.org/10.1096/fba.2024-00048).

This is material transport on a morphology-derived numerical domain, not
membrane voltage, action potentials, synapses, calibrated neuronal biology, or
vascular hemodynamics.

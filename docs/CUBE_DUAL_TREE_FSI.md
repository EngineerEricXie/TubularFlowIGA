# Idealized cube dual-tree FSI case

This functional case couples two synthetic vessel trees through a fixed Darcy
region and applies quasi-steady vessel-wall feedback. The geometry and all
parameters are defined in
[`cases/idealized_cube_dual_tree.json`](../cases/idealized_cube_dual_tree.json).

The baseline tree contains seven cylindrical segments and four terminals. The arterial
tree enters at `x=0`; the venous tree exits at `x=0.03 m`. Gmsh OCC Boolean
operations create five disjoint, conforming tetrahedral regions:

1. arterial lumen;
2. arterial wall;
3. venous lumen;
4. venous wall; and
5. the surrounding Darcy region.

Gmsh supplies geometry and meshing only. CoupledFlow implements the FEM weak
forms, elements, assembly, transfer, and acceptance checks; PETSc supplies
sparse algebra and MPI.

The case is designed to verify geometry, coupling, conservation, and
cross-rank reproducibility. It is not an anatomical or physiological liver
model. The companion [passive tracer case](CUBE_DUAL_TREE_OXYGEN.md) advances
concentration on the converged, frozen post-FSI flow field.

## Run

The geometry scripts require the Gmsh Python package and NumPy. The solver
requires MPI and PETSc. PNG rendering additionally requires VTK Python
bindings and Matplotlib; ParaView can read the PVTU output without them.

```bash
make -C solvers/cpu native_tet_cube_dual_tree_fixed_flow \
  PETSC_DIR=/usr/lib/petscdir/petsc3.15/x86_64-linux-gnu-real
python3 -m unittest scripts.tests.test_idealized_cube_dual_tree

python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py \
  /tmp/idealized-cube-r1 --ranks 1
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py \
  /tmp/idealized-cube-r2 --ranks 2
python3 scripts/compare_idealized_cube_dual_tree.py \
  /tmp/idealized-cube-r1 /tmp/idealized-cube-r2

python3 scripts/visualize_idealized_cube_dual_tree.py \
  /tmp/idealized-cube-r1
python3 scripts/export_idealized_cube_paraview.py \
  /tmp/idealized-cube-r1
```

Output directories must not exist. The runner records case, solver, geometry,
partition, convergence, conservation, and FSI identities in `summary.json`.
Rank-owned VTU/PVTU fields are written under `fields/`. WSL and container
launches must permit local MPI loopback sockets.

## Numerical model

- The arterial and venous lumens use the project P2-velocity/P1-pressure
  steady tetrahedral flow formulation.
- The tissue region uses P1 Darcy pressure and conservative RT0 face fluxes.
- Arterial terminal face fluxes are mapped to adjacent Darcy cells by matching
  geometry, unique ownership, and opposing normals.
- Darcy flux at the four venous interfaces supplies the venous inlet controls.
- Each wall uses small-strain P1 tetrahedral linear elasticity. The outer wall
  and root cap are fixed; the lumen surface receives fluid Cauchy traction.
- Elastic ALE extension moves each lumen mesh. Flow, Darcy pressure, and wall
  displacement are then resolved until the displacement fixed-point criterion
  and conservation checks pass.

The hydraulic transfer is conservative and one-way: arterial flux drives
Darcy sources, and Darcy flux drives venous inlet flow. The case does not
enforce vessel--tissue pressure continuity. The surrounding Darcy mesh remains
fixed, and the procedure is quasi-steady rather than time-accurate transient
FSI.

The baseline functional parameters are:

| Parameter | Value |
|---|---:|
| Inlet velocity | `0.002 m/s` |
| Fluid density | `1000 kg/m³` |
| Dynamic viscosity | `0.004 Pa·s` |
| Darcy mobility | `10⁻⁸ m²/(Pa·s)` |
| Wall thickness | `0.0006 m` |
| Young's modulus | `2×10⁴ Pa` |
| Poisson ratio | `0.3` |

The four arterial exchange caps, venous outlet, and four Darcy-to-vein caps
use the pressure conditions specified by the case file. Remaining tissue
boundaries are no-flux.

## Recorded validation

A 2026-09-20 WSL validation compared independent one- and two-rank runs of the
baseline case:

| Quantity | Result |
|---|---:|
| Tetrahedra | `70,791` |
| Minimum scaled Jacobian | `0.0127652` |
| Artery to Darcy flow | `6.544721406613861×10⁻⁹ m³/s` |
| Darcy to vein flow | `6.544721406613780×10⁻⁹ m³/s` |
| Venous outlet flow | `6.544721406613783×10⁻⁹ m³/s` |
| Maximum arterial/venous wall displacement | `2.8535×10⁻⁸ / 2.2939×10⁻⁸ m` |
| FSI displacement residual | `2.3430×10⁻⁵` |
| Minimum deformed-fluid Jacobian ratio | `0.9999969` |
| Maximum cross-rank absolute field difference | `3.2685×10⁻¹²` |

The geometry and partition hashes matched across runs, and all 36 rank-owned
field arrays were compared by global ID.

## Refined thin-wall variant

[`cases/idealized_cube_dual_tree_refined.json`](../cases/idealized_cube_dual_tree_refined.json)
keeps the same centerlines, lumen radii, materials, and boundary conditions but
changes the wall thickness from `0.60 mm` to `0.10 mm`. Curvature-based sizing
increases the approximate circumference resolution from 8 to 32 nodes and
uses a minimum mesh size of `0.04 mm`.

```bash
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py \
  /tmp/idealized-cube-refined \
  --case cases/idealized_cube_dual_tree_refined.json --ranks 8
python3 scripts/visualize_idealized_cube_dual_tree.py \
  /tmp/idealized-cube-refined \
  --case cases/idealized_cube_dual_tree_refined.json --magnification 5000
python3 scripts/export_idealized_cube_paraview.py \
  /tmp/idealized-cube-refined --magnification 5000
```

The refined mesh contains `271,870` tetrahedra with minimum scaled Jacobian
`0.00330819`. Five regions, 12 matching interfaces, native P2 topology, and
positive Jacobians pass the mesh checks. Independent two- and eight-rank runs
recorded:

| Quantity | Result |
|---|---:|
| Artery to Darcy flow | `8.133254228498714×10⁻⁹ m³/s` |
| Darcy to vein flow | `8.133254228498833×10⁻⁹ m³/s` |
| Venous outlet flow | `8.133254228498835×10⁻⁹ m³/s` |
| Maximum Darcy cell-balance defect | `8.17×10⁻²³ m³/s` |
| FSI displacement residual | `2.80×10⁻⁵` |
| Minimum deformed-fluid Jacobian ratio | `0.999979` |
| Maximum arterial/venous wall displacement | `6.92×10⁻⁹ / 5.15×10⁻⁹ m` |
| Maximum cross-rank absolute field difference | `6.89×10⁻¹¹` |

The two runs used identical case, solver, geometry, and partition hashes. Their
wall-clock times are not treated as a scaling benchmark because they partly
overlapped on the same machine.

## Template-free smooth-tree variant

[`cases/idealized_cube_dual_tree_smooth.json`](../cases/idealized_cube_dual_tree_smooth.json)
keeps the refined thin wall and replaces each centerline edge with a cubic
Hermite path sampled into overlapping capsules. Short Boolean seam edges are
collapsed, then a windowed-sinc pass rounds the vessel walls and junctions
without moving inlet or terminal boundary rings. A conforming discrete-surface
remesh creates the tetrahedra. Node tangents combine the
incoming direction with the mean outgoing direction, so the construction does
not depend on bifurcation templates and accepts branch nodes with two, three,
or more children. The fixed-flow reference executable still requires four
arterial and four venous terminal caps.

```bash
python3 scripts/generate_idealized_cube_dual_tree.py \
  /tmp/idealized-cube-smooth \
  --case cases/idealized_cube_dual_tree_smooth.json
python3 scripts/run_idealized_cube_dual_tree_fixed_flow.py \
  /tmp/idealized-cube-smooth-run \
  --case cases/idealized_cube_dual_tree_smooth.json --ranks 8
```

The geometry-only check produces five conforming regions and all 12 matching
interfaces while retaining the native FEM minimum scaled-Jacobian gate of
`0.001`. Gmsh 4.8.4 produced `272,751` first-order tetrahedra with a measured
minimum of `0.00499229`. This variant has not replaced the recorded flow and
FSI evidence above; use the refined case when
reproducing those numerical values.

## ParaView output

Open `paraview/arterial_pipe_x1000.pvd` and
`paraview/venous_pipe_x1000.pvd`, apply both as surfaces, and switch between
states `0` and `1`. These states mean reference and converged quasi-steady FSI;
they are not seconds in a transient simulation. Files ending in `_actual.pvd`
contain the physical displacement. Magnified files change display geometry
only; the `displacement_m` array remains in SI units.

Use `tissue_darcy_velocity.vtu` with Clip or Slice to inspect the interior
Darcy mesh. Cell arrays `darcy_velocity_m_s` and `darcy_speed_m_s` carry the
RT0-reconstructed velocity. `tissue_darcy_direction_arrows.vtp` uses fixed
arrow lengths for visibility, so color the cell field when comparing
magnitudes.

The lumen and wall deformation files contain `deformation_jacobian` and
`relative_cell_volume_change`. `volume_audit.json` records integrated reference
and deformed volumes. Native PVTU coordinates are already displaced; applying
Warp By Vector to them would duplicate the motion.

The exporter rejects magnification factors that invert first-order tetrahedra.
This geometric display check does not replace a high-order surface
self-intersection analysis.

# Explicit vessel-to-tissue source map: generic functional contract

This is a project-owned **data/volume-source mapping layer**, not a completed
vessel–liver coupling. It has no inferred anatomy, physiology, geometry
intersection, or external FEM dependency. `NativeTetVesselTissueSourceMap.hpp`
takes, for each explicitly named vessel port, a signed outward flow `Q_v`
in `m³/s` and a list of target tissue tetrahedron IDs with positive weights
`w_i` summing to one. It returns Darcy cell sources
`s_i = Σ_ports Q_v w_i / V_i` in `1/s`, where `V_i` is the current tetrahedron
volume in `m³`. The tissue receives `Q_t=Σ_i s_i V_i=Q_v`; reverse flow is
negative. Equivalently, if tissue exchange is defined outward-positive,
`Q_v + Q_t,out = 0` for each port.

The API rejects absent/duplicate cell IDs, duplicate/empty port names,
missing targets, nonpositive/unnormalized weights, nonfinite flows and
invalid tetrahedra. A port's target cells are supplied by the user; **cell
IDs alone are not evidence that a vessel and a tissue region touch**. The
caller must verify that geometric relation, tissue material, and measured
vessel-port flow before any organ interpretation. Multiple distinct ports
may contribute to one tissue cell, and their signed sources add; port IDs
must remain unique so that the same vessel flux is not silently counted
twice.

The local acceptance case uses a manufactured unit cube: port A transfers
`+0.6 m³/s` to cells 1 and 2 with weights `0.25/0.75`; port B returns
`0.2 m³/s` from cell 3 to its vessel. Each tetra has volume `1/48 m³`, so
the resulting sources are `7.2`, `21.6`, `−9.6 s⁻¹`, and the net Darcy
volume source is `+0.4 m³/s`. Native Darcy solves this source field and
its recovered RT0 face flux passes per-cell and global balance on 1/2/4
MPI ranks. Bad weights, absent cells and repeated port names are rejected.
Run `make t6-native-tet-darcy-test PETSC_DIR=...`.

The second project-owned entrypoint,
`MapNativeTetVesselStateToTissueSource`, takes a P2/P1 vessel state and an
explicit `name → boundary_label → tissue cell IDs/weights` table. It
integrates the **actual supplied FEM field** over each labelled boundary
using `EvaluateNativeTetBoundaryFlows`, then applies the same conservative
source map. Missing or duplicate boundary labels, zero-area labels and
nonfinite states fail closed. A manufactured constant `u_x=0.6 m/s` on a
unit cube yields `Q(x=1)=+0.6 m³/s` and `Q(x=0)=−0.6 m³/s`; the mapped
Darcy source and recovered RT0 flux balance pass on 1/2/4 MPI ranks.

This still does **not** prove that any labelled vessel boundary touches the
declared tissue cells, nor has it connected a solved vascular state to a
same-case liver mesh or supplied liver-specific BC/material values. Those
remain outside the supported scope described in the [technology matrix](../TECHNOLOGY_MATRIX.md).
The constant-field cube is a functional/manufactured test, not physical or
physiological validation.

## Exact matching-face geometry gate

`NativeTetMatchingInterfaceSource.hpp` adds an explicitly labelled,
**matching-facet-only** route. Each port names one vessel boundary label and
one tissue boundary label. The routine requires an exact one-to-one match of
the selected triangular facets' Cartesian coordinates (metres), unique
boundary tetrahedron ownership, and opposed outward normals. It integrates
the project-owned P2 vessel velocity **on every matched facet**, then puts
that signed `Q_f [m³/s]` into the adjacent tissue tetrahedron as
`s_i = Σ_{f→i} Q_f/V_i [1/s]`; no user-entered cell weights or proximity
guess are involved. P2 mid-edge basis integration is exact for quadratic
velocity traces. Missing, displaced, duplicate, same-side, or wrongly labelled
facets, nonfinite state, repeated port labels, and reusing one geometric facet
under two different port labels fail closed. The result also
returns an ordered per-facet ledger with the two adjacent global cell IDs,
sorted triangle coordinates, port name, and signed vessel-outward `Q_f`; the
same values feed the cell-source accumulation, so reverse flow remains visible
instead of disappearing into a net port total.
The tissue interface face uses its zero natural-flux condition in this
**adjacent-cell volumetric-source model**; the transfer is conservative per
facet/cell in the integrated sense, not a claim of pointwise velocity or
traction continuity across the two PDEs.

A manufactured pair of adjacent tetrahedra gives `Q=+0.3 m³/s` and
`s=+1.8 s⁻¹` for constant `u_x=0.6 m/s`; reversing velocity gives the
negative values. Native P1 Darcy with `K=1 m²/(Pa s)`, one `0 Pa` exit face,
and zero natural flux elsewhere recovers `±0.3 m³/s` conservative outflow
and per-cell balance on 1/2/4 ranks. A variable quadratic P2 trace agrees
with the existing boundary quadrature. A two-facet case sends different
flows `1/3` and `8/15 m³/s` to their respective neighbouring tissue cells,
rather than distributing one global total by area. These are dimensioned
manufactured tests, not physiological parameters.

This route does not construct a matching mesh, validate that entire 3D
regions are disjoint, or support nonconforming/cut interfaces. The early,
independently generated TCIA liver/vessel trial meshes cross and are **not**
eligible. A later same-SEG, voxel-faithful **local ROI** does provide disjoint,
matching Portal/tissue tetrahedra and passes this functional map with an
explicitly artificial solved fluid field. Its 158-facet ledger is checked
against both mesh interfaces and every tissue cell source, including reverse
flows; see `docs/LIVER_GEOMETRY_CANDIDATES.md`. This does not complete the
full-organ L1–L4 gates or establish physiological liver perfusion.

The idealized Y runtime now also tests this route with its **solved** native
P2/P1 blood velocity, not only a manufactured constant field. The test
constructs an explicitly artificial, external tetrahedral cone behind each
outlet cap (25 exact matching facets per cap, 50 tissue tetrahedra total),
without reusing any vessel tetrahedron volume. Each matched outlet face feeds
only its adjacent Darcy cell; the tissue's external faces use `0 Pa`, the
interface faces zero natural flux, and the cell sources model the transfer.
On 1/2/4 ranks, the summed source and conservative Darcy external outflow are
both about `3.82492727e−7 m³/s`, maximum Darcy cell-balance defect is about
`1.2e−23 m³/s` or less, and its PETSc converged reason is `3`.
The source-to-outflow test tolerance is `1e−8` relative to the exchange flow.
The cones and `K=1 m²/(Pa s)` are dimensioned **functional fixtures**, not
liver tissue morphology or physiological parameters. This strengthens the
generic field→matching-face→adjacent-cell→Darcy test; it does not justify
applying this interface map to the independent overlapping TCIA meshes.

An additional idealized Y-bifurcation regression now takes a **solved**
project-owned P2/P1 vessel state, integrates its two outlet-labelled flows,
and maps both into one explicitly named manufactured tissue tetra, then
solves that tetra with the project-owned P1 Darcy FEM and RT0 conservative
flux recovery (`K=1 m²/(Pa s)`, zero-pressure artificial tissue exit).
These values and the disconnected tissue shape are solely functional inputs.
The vessel pressure BC is `0 Pa` on both outlet labels; a point pressure gauge was
insufficient for this regression because it yielded about `2.47%` global
boundary-flow imbalance despite a small algebraic residual. With the
explicit pressure outlets the 1/2/4-rank steady solve reports inlet
`−3.67717e−7 m³/s`, outlet flows about `1.83464e−7` and
`1.84253e−7 m³/s`, global imbalance near roundoff, and per-port
source-map balance at the declared `1e−12` relative tolerance. Darcy's
conservative outlet flow equals the mapped `3.677167777e−7 m³/s` source,
its cell balance defect is zero in the test, and PETSc converged reason is
positive (`3`) on 1/2/4 ranks. The tissue
tet is **not geometrically connected** to the Y; this tests the numerical
field→flux→source interface, not a vessel–tissue anatomy or physiology.

The WSL smoke can be regenerated without committing patient data or mesh
artifacts. The generator creates an idealized Y surface with labels 0–3;
the converter runs boundary/quality gates and writes the mesh/source hashes
to `y.json`. The runtime writes one result per rank count; its JSON includes
convergence reasons and source/flux balances. Replace `PETSC_DIR` if the
local installation differs. The `0.003 m` surface target and `0.006 m`
volume target below passed the 2-rank test locally; the coarser
`prepare_t2_bifurcation_meshes.py --maximum-level 1` mesh failed the
test's near-even outlet-split gate (0.562), despite passing conservation,
and must not be substituted as a passing fixture.

```bash
make -C solvers/cpu surface_fem_preflight
mkdir -p /tmp/native-y-source-map
python3 scripts/generate_y_pipe_surface.py /tmp/native-y-source-map/y.vtp \
  --target-size-m 0.003
python3 scripts/surface_to_fem_volume.py /tmp/native-y-source-map/y.vtp \
  /tmp/native-y-source-map/y.msh --manifest /tmp/native-y-source-map/y.json \
  --target-size-m 0.006
make -C solvers/cpu native_tet_bifurcation_runtime_test PETSC_DIR=/usr/lib/petsc
for ranks in 1 2 4; do
  mpiexec -np "$ranks" solvers/cpu/native_tet_bifurcation_runtime_test \
    /tmp/native-y-source-map/y.msh \
    "/tmp/native-y-source-map/solved-coupling-r${ranks}.json" \
    -bifurcation_steady true -bifurcation_steps 1
done
```

The mesh is an analytic-functionality fixture, not a liver input. Solver
results can vary slightly across rank counts; compare declared conservation
gates and numerical fields, not bitwise hashes of the result JSON.
On the regenerated 2,536-tetra fixture, 1/2/4 ranks all passed: the 2-rank
mapped source and conservative Darcy outflow were both
`3.8249272699962092e−7 m³/s`, its tissue cell balance defect was zero,
the vessel relative mass imbalance was `6.92e−17`, and Darcy's PETSc
converged reason was `3`.

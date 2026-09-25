# Supported input, unit, and label contract

The authoritative machine-readable contract is
[`benchmarks/supported_input_contract.json`](../benchmarks/supported_input_contract.json).
It records the locally implemented input routes, exact minimum case assets,
source-to-SI conversions, boundary-label rules, generated solver interfaces,
and explicit rejection policy. `make t0-audit` checks both its schema and a
small set of high-risk claims directly against the production readers.

The short operational rules are:

- Solver-facing geometry is in metres; time, pressure, velocity, and flow are
  seconds, pascals, metres per second, and cubic metres per second. SWC/line
  OBJ and surface preprocessing accept only an explicit positive scale to
  metres. Production immersed geometry is already in metres.
- A centerline case uses strict seven-column SWC or the documented
  radius-annotated **line** OBJ. Surface OBJ is not accepted by that route.
- The common surface preflight accepts closed triangular VTP or STL. VTP
  carries a nonnegative integer `boundary_id` per triangle. STL has no labels,
  so it is usable only with one explicit default label.
- OFF/OBJ/STL/PLY are capabilities of the separately installed fTetWild CLI.
  The adapter still begins from the common VTP/STL preflight; therefore OFF,
  OBJ, and PLY are not advertised as project surface inputs.
- Production immersed cases accept VTP, `simulation_config.json`, and
  `immersed_geometry.json`. Wall and port labels must be disjoint and exactly
  cover every surface label. Moving frames retain topology and material IDs.
- Native tetrahedral FEM consumes ASCII Gmsh 4.1 with first-order triangles
  and tetrahedra, physical surface names `boundary_label_<integer>`, and one
  volume named `fluid`. The repository owns P2/P1 basis construction,
  quadrature, assembly, boundary conditions, and QoIs; PETSc is algebra/MPI.
- No route or physical boundary role is guessed from organ name, axis,
  numeric label order, file order, or MPI ownership. A case configuration
  assigns wall, inlet, outlet, coupling, or constraint semantics explicitly.

The generated `.ntiga`, `.msh`, cache, partition, profile, and manifest files
remain case data. They preserve the documented file interfaces and are
not repository source artifacts.

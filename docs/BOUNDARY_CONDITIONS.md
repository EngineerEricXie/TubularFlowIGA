# Boundary-condition configuration

New cases should define field-specific conditions in the canonical
[`simulation_config.json`](PDE_CONFIGURATION.md). That schema supports arbitrary
field names, fixed or profile-scaled velocity Dirichlet data, and scalar
Dirichlet/no-flux/outflow/flux/Robin conditions. The CPU and CUDA
configured-transport solvers assemble matching scalar `flux` and `robin`
surface integrals. Named temporal waveforms can now be parsed and
evaluated. CPU and CUDA configured transport update waveform-backed Dirichlet
values at each step. CPU backward-Euler Navier–Stokes also updates velocity
Dirichlet and pressure-traction values at each physical step; steady flow
rejects temporal values. Validation evidence is recorded in the
[CPU](../solvers/cpu/VALIDATION.md) and
[CUDA](../solvers/cuda/VALIDATION.md) validation reports.

For configured transport, `flux` specifies the outward diffusive flux
`D grad(c) dot n = q` and contributes `dt q N_a` to the right-hand side.
`robin` uses `D grad(c) dot n = h(c_ext - c)`, contributing
`dt h N_a N_b` to the left-hand side and `dt h c_ext N_a` to the
right-hand side. Surface execution requires a version 4 `.ntiga` database
containing element-face labels. A configured surface label with no packed face
is rejected with a request to rerun `iga_pack`.

## Physiological pressure outlets

An open flow boundary should normally use `pressure_traction`:

```json
{"field": "pressure", "type": "pressure_traction", "value": 80.0}
```

It applies the natural Cauchy traction `t = -p n` to the momentum residual and
keeps the continuity equation on every outlet pressure degree of freedom. It may
name a temporal `waveform`; the scalar value is multiplied by that function at
each backward-Euler step. A pressure `dirichlet` condition is still accepted for
legacy cases or an explicit gauge, but constraining an entire outlet pressure
field replaces those nodes' continuity rows and is not equivalent to a pressure
traction.

CPU and CUDA Navier–Stokes accept three pressure outlet types. R, RC, and RCR
all passed the CUDA fixed-point path, including RCR capacitor-state restart.
All flow rates use the outward normal, so a healthy outlet normally has
positive `Q`.

- `resistance`: `p = reference_pressure + resistance * Q`.
- `windkessel_rc`: a parallel RC model advanced with backward Euler.
- `windkessel_rcr`: a proximal resistance followed by a parallel distal
  resistance/capacitance model.

For example:

```json
{
  "field": "pressure",
  "type": "windkessel_rcr",
  "proximal_resistance": 1200.0,
  "distal_resistance": 8800.0,
  "capacitance": 0.00012,
  "reference_pressure": 0.0,
  "initial_pressure": 80.0
}
```

The solver integrates `u dot n` from version 4 packed boundary faces and
fixed-point couples each outlet pressure to the 3D solve. R/RC/RCR values are
lowered to the same natural `-p n` traction, not pressure Dirichlet nodes.
RC/RCR require `backward_euler`; a pure resistance also works in steady flow. Checkpoints store
flow, applied pressure, and capacitor pressure for consistent restart. CUDA
uses the same host-side surface quadrature and outlet equations around its GPU
3D solve, and stores the same outlet metadata in raw-state checkpoints.


## Transitional schema v1

Older cases may contain an optional `case_config.json`. CPU and CUDA read the
same file and resolve it against the integer point labels in `controlmesh.vtk`.
The mesh generator uses this label contract:

- `0`: tube wall
- `1`: root/inlet cap
- `2` and above: terminal/outlet caps
- `-1`: unconstrained interior point

If the file is absent, the historic behavior is preserved: label 0 is a
no-slip wall, label 1 uses `vplus * initial_velocityfield.txt` and the
`N0bc`/`Nplusbc` transport values, and every label of 2 or greater has zero
pressure.

An explicit case can override all boundaries:

```json
{
  "schema_version": 1,
  "boundaries": {
    "inherit_legacy": false,
    "conditions": [
      {"label": 0, "name": "wall", "type": "wall"},
      {
        "label": 1,
        "name": "inlet",
        "type": "inlet",
        "velocity_scale": 1.0,
        "transport": {"N0": 1.0, "Nplus": 2.0}
      },
      {"label": 2, "name": "outlet A", "type": "outlet", "pressure": 0.0},
      {"label": 3, "name": "outlet B", "type": "outlet", "pressure": 0.25}
    ]
  }
}
```

For an inlet, use either `velocity_scale` to multiply the generated reference
field or `velocity: [vx, vy, vz]` for a fixed vector. A wall defaults to
`[0, 0, 0]` but may also specify a fixed velocity. An outlet may specify its
pressure. `transport` is valid on inlet rules and sets the two transported
Dirichlet values.

With `inherit_legacy: true`, listed labels override the historic rules and
unlisted labels retain them. With `false`, every non-negative label present in
the mesh must have exactly one rule. Unknown keys, duplicate labels, invalid
type/value combinations, mesh labels missing from the configuration, and
configured labels absent from the mesh are rejected.

Validate a prepared case without launching PETSc or CUDA:

```bash
./solvers/cpu/iga_case_check DATABASE.ntiga CASE_DIR
```

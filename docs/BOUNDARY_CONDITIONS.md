# Boundary-condition configuration

New cases should define field-specific conditions in the canonical
[`simulation_config.json`](PDE_CONFIGURATION.md). That schema supports arbitrary
field names, fixed or profile-scaled velocity Dirichlet data, and scalar
Dirichlet/no-flux/outflow/flux/Robin conditions. The CPU and CUDA
configured-transport solvers assemble matching scalar `flux` and `robin`
surface integrals. Named temporal waveforms can now be parsed and
evaluated. CPU configured transport updates waveform-backed Dirichlet values at
each step. Steady Navier–Stokes explicitly rejects waveforms pending transient
time stepping, so cardiovascular inlet profile scaling remains constant in
actual flow runs. Track these increments in the [development roadmap](ROADMAP.md).

For configured transport, `flux` specifies the outward diffusive flux
`D grad(c) dot n = q` and contributes `dt q N_a` to the right-hand side.
`robin` uses `D grad(c) dot n = h(c_ext - c)`, contributing
`dt h N_a N_b` to the left-hand side and `dt h c_ext N_a` to the
right-hand side. Surface execution requires a version 4 `.ntiga` database
containing element-face labels. A configured surface label with no packed face
is rejected with a request to rerun `iga_pack`.


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

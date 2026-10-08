# Prescribed LV numerical QoI diagnostic

The native P2/P1 tetrahedral ALE solver runs the idealized 16-step ED–ES–ED
motion with density `1050 kg/m³`, viscosity `0.012 Pa·s`, outward inlet target
`-1e-6 m³/s`, zero-gauge pressure outlet, and optional outlet backflow term
`beta=0.5`. The source surface, labels, wall motion, and port conditions are
the same for the 96-, 384-, and 672-tetra interior meshes. FEM integration and
assembly remain project-owned; PETSc supplies algebra and MPI.

At each accepted step the diagnostic reports chamber volume, volume-mean P1
pressure (relative to the zero-pressure outlet gauge), kinetic energy
`0.5*rho*integral(|u|² dV)`, and outward port flows. It also sums
`0.5*(p_previous+p_current)*(V_current-V_previous)` over the 16 steps, taking
the initial mean pressure as zero. This pressure–volume integral is a
**gauge-relative numerical proxy**, not myocardial stroke work: the wall
motion is prescribed, there is no myocardial model, and a full mechanical
energy balance would additionally require boundary traction work, viscous
dissipation, and flux terms.

| Native ALE volume mesh | Peak absolute mean pressure | Step-9 mean pressure | Final mean pressure | Pressure–volume integral | Initial → final kinetic energy |
|---|---:|---:|---:|---:|---:|
| 96 tetra | 23.0454 Pa | −23.0454 Pa | −0.869730 Pa | −8.36161e−5 J | 7.69526e−7 → 1.59267e−5 J |
| 384 tetra | 20.6076 Pa | −20.6076 Pa | −0.665956 Pa | −7.59315e−5 J | 7.69526e−7 → 1.38402e−5 J |
| 672 tetra | 20.0952 Pa | −20.0952 Pa | −0.632867 Pa | −7.42476e−5 J | 7.69526e−7 → 1.35102e−5 J |

All three cases pass the inlet-flow, boundary-flow, and GCL gates, and geometric
volume returns to its initial value within `3.2e-14` relatively. The
pressure–volume integral changes by 10.12% from 96 to 384 tetrahedra (384
denominator), then by 2.27% from 384 to 672 (672 denominator). The largest
mean-pressure difference between the first pair is 2.438 Pa at the first
expansion step; between the second pair it is 0.512 Pa at the same step.
The fluid kinetic energy does **not** return to its
initial value, so one prescribed-motion cycle is not evidence of a periodic
fluid state. The radial grading changes between meshes, and only the first
two levels have the two-cycle diagnostic below. Neither asymptotic mesh
convergence nor physiological pressure validation has been established.

An explicit second prescribed-motion cycle confirms that distinction. Both
96- and 384-tetra meshes accept all 32 steps while closing inlet flow, boundary
flow, and GCL.
At the identical end-diastolic geometry, the relative P2 velocity L2
difference between steps 16 and 32 (normalized by the step-32 field) is
`27.40%` on 96 tetrahedra and `37.22%` on 384 tetrahedra. The respective
kinetic energies change from `1.59267e-5` to `1.86007e-5 J`, and from
`1.38402e-5` to `1.63290e-5 J`. Thus two geometry cycles still do not
establish a periodic flow solution. More cycles or a dedicated periodic-state
solve would be needed before interpreting a pressure–volume loop as a
steady-cycle quantity.

The high-precision per-step series and definitions are in
`benchmarks/t6_native_lv_qoi_evidence.json`; `make t6-native-lv-qoi-audit`
recomputes the pressure–volume sums and checks the scope of the claims. The
three native one-cycle runs are reproduced with
`native_tet_ale_idealized_lv_backflow_flow_test` and
`native_tet_ale_idealized_lv_radial_backflow_flow_test` and
`native_tet_ale_idealized_lv_two_shell_backflow_flow_test` after a PETSc build on
an allocated compute resource. The two-cycle diagnostics use the corresponding
`native_tet_ale_idealized_lv_two_cycle_backflow_flow_test` and
`native_tet_ale_idealized_lv_radial_two_cycle_backflow_flow_test` binaries.

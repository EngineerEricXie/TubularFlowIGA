# Native zero-dimensional circuit networks

The directory separates lumped circuits from spatially distributed 1D PDE
examples:

| Case | Formulation | Purpose |
|---|---|---|
| `steady_resistive_straight` | `steady_r` | Poiseuille R network and analytic pressure drop |
| `rc_rcr_bifurcation` | `transient_rc` | vessel compliance with RC and RCR terminal beds |

The `rc_rcr_bifurcation` case is a transient lumped-parameter circuit on a
centerline tree. Vessel resistances and compliances are selected explicitly by
the child node ID of each segment. One terminal uses a two-element RC model and
the other uses a three-element RCR Windkessel.

```bash
make zero-d-petsc
./solvers/one_d/iga_0d examples/zero_d/steady_resistive_straight --check
./solvers/one_d/iga_0d examples/zero_d/steady_resistive_straight
./solvers/one_d/iga_0d examples/zero_d/rc_rcr_bifurcation --check
./solvers/one_d/iga_0d examples/zero_d/rc_rcr_bifurcation
```

The backward-Euler circuit solve is MPI-distributed through PETSc. Inspect
`outlet_timeseries.csv` for terminal pressure, capacitor pressure, distal
flow, and capacitor storage rate. `flow_timeseries.csv` and
`branch_timeseries.csv` contain the network balance and segment states.

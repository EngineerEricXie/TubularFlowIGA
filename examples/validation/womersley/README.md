# Womersley reference generation

`iga_womersley_reference` creates analytical point samples for a straight,
rigid, circular tube in `controlmesh.vtk` node ordering. These are useful for
visualization and inlet-profile preparation, but they are not generally the IGA
control coefficients of the analytical field.

From the repository root:

```bash
make cpu
./solvers/cpu/iga_womersley_reference \
  CASE.ntiga CASE_DIR/controlmesh.vtk \
  examples/validation/womersley/womersley_reference.json REFERENCE_DIR

./solvers/cpu/iga_flow_validate CASE.ntiga \
  --womersley examples/validation/womersley/womersley_reference.json \
  NUMERICAL_DIR velocity_series.csv
```

The pressure-gradient convention is positive driving force along
`axis_direction`:

```text
G(t) = -dp/ds
     = mean_pressure_gradient
       + sum(a_n cos(n omega t) + b_n sin(n omega t)).
```

For a tube of length `L` with zero outlet pressure, prescribe inlet pressure
coefficients equal to `L` times the gradient coefficients. Density, dynamic
viscosity, coordinates, radius, time, and pressure gradient must use one
consistent unit system. Snapshot times in the configuration select which
numerical manifest times are validated; the numerical manifest may also contain
earlier startup snapshots. `--womersley` evaluates the numerical spline and
analytical solution at element volume quadrature points and integrates the
physical relative L2 norm. Use `--compare-manifests` only for CPU/CUDA
coefficient parity or an independently projected coefficient reference. Edit
the supplied values for the actual tube; this file is a schema and workflow
example, not a physiological parameter set.

The analytical implementation includes the steady Poiseuille component and
the complex-Bessel Womersley response for every supplied Fourier harmonic.
Points farther outside `radius` than the legacy six-decimal VTK tolerance are
rejected.

`axis_origin`, `radius`, and all other geometric values are expressed in the
normalized database coordinate system. The spline stage subtracts each input
axis minimum and divides all coordinates by the shortest input-axis extent.
The generator applies that same transformation to `controlmesh.vtk`; for the
example tube spanning `y,z = [-0.5,0.5]`, its database center is therefore
`[0,0.5,0.5]`.

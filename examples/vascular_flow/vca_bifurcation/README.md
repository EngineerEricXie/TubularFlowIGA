# Native 3D VCA bifurcation smoke case

This source-only case exercises the CPU 3D explicit-staggered VCA path with a
flow-controlled pump, one reservoir species, and two venous outlets. Generate
the mesh/database on a suitable local or allocated resource; generated files
are intentionally not committed.

With the PETSc-enabled solver built, run the complete one-versus-two-rank
validation below. `prepare_example.sh` requires an empty work directory and
generates the two-rank partition; the `awk` command derives the all-owned
single-rank partition from it. Use a scheduler allocation for large production
cases.

The recommended form is the automated validator; set
`VCA_RELATIVE_TOLERANCE`, `VCA_PRESSURE_RELATIVE_TOLERANCE`, or
`VCA_VOLUME_BALANCE_TOLERANCE` if the site-specific solver tolerance requires
it. The pressure default is `1e-4`: near-zero pressure fields make a tighter
relative-only comparison misleading.

```bash
VCA_RELATIVE_TOLERANCE=1e-6 ./scripts/validate_vca_bifurcation.sh /tmp/vca-bifurcation
```

The commands below are the equivalent expanded procedure.

```bash
WORK=/tmp/vca-bifurcation
RANKS=2 ./scripts/prepare_example.sh vascular_flow/vca_bifurcation "$WORK"

awk '{print 0}' "$WORK/bzmeshinfo.txt.epart.2" > "$WORK/bzmeshinfo.txt.epart.1"
./solvers/cpu/iga_pack "$WORK" 1 "$WORK/vca_bifurcation-1.ntiga"

mkdir -p "$WORK"/{one,two}/{full,split,resumed}

mpiexec -np 1 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-1.ntiga" "$WORK" \
  --checkpoint "$WORK/one/full/checkpoint" --output "$WORK/one/full/flow.txt"
mpiexec -np 1 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-1.ntiga" "$WORK" \
  --stop-after-step 2 --checkpoint "$WORK/one/split/checkpoint" --output "$WORK/one/split/flow.txt"
mpiexec -np 1 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-1.ntiga" "$WORK" \
  --restart "$WORK/one/split/checkpoint" --checkpoint "$WORK/one/resumed/checkpoint" \
  --output "$WORK/one/resumed/flow.txt"

mpiexec -np 2 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-2.ntiga" "$WORK" \
  --checkpoint "$WORK/two/full/checkpoint" --output "$WORK/two/full/flow.txt"
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-2.ntiga" "$WORK" \
  --stop-after-step 2 --checkpoint "$WORK/two/split/checkpoint" --output "$WORK/two/split/flow.txt"
mpiexec -np 2 ./solvers/cpu/iga_navier_stokes "$WORK/vca_bifurcation-2.ntiga" "$WORK" \
  --restart "$WORK/two/split/checkpoint" --checkpoint "$WORK/two/resumed/checkpoint" \
  --output "$WORK/two/resumed/flow.txt"

cmp "$WORK/one/full/checkpoint.state" "$WORK/one/resumed/checkpoint.state"
cmp "$WORK/one/full/checkpoint.vca_transport.state" "$WORK/one/resumed/checkpoint.vca_transport.state"
cmp "$WORK/two/full/checkpoint.state" "$WORK/two/resumed/checkpoint.state"
cmp "$WORK/two/full/checkpoint.vca_transport.state" "$WORK/two/resumed/checkpoint.vca_transport.state"

relative_l2() {
  awk 'NR == FNR { for (i = 1; i <= NF; ++i) reference[FNR, i] = $i; columns[FNR] = NF; next }
       NF != columns[FNR] { exit 2 }
       { for (i = 1; i <= NF; ++i) { delta = $i-reference[FNR, i]; error += delta*delta; norm += reference[FNR, i]*reference[FNR, i] } }
       END { printf "relative_l2=%.17g\\n", sqrt(error)/(sqrt(norm) > 0 ? sqrt(norm) : 1) }' "$1" "$2"
}
relative_l2 "$WORK/one/full/flow.txt" "$WORK/two/full/flow.txt"
relative_l2 "$WORK/one/full/flow.txt.pressure" "$WORK/two/full/flow.txt.pressure"
```

For the one-versus-two-rank field check, compute the relative L2 difference
for both `flow.txt` and its `.pressure` companion; retain the command output
and the six `coupling_manifest.json` files with the run record. Each full-run
manifest must contain four arterial and venous records, while split and resumed
manifests contain two; all must report outlet labels `2` and `3`, oxygen species
flux, and reservoir volume changes close to zero for an incompressible closed
loop. Compare their final reservoir state and port history at the chosen
numerical tolerance. Do not commit any generated case, database, partition,
checkpoint, or result file.

# CoupledFlow technology matrix

Updated 2026-09-21.

This matrix maps physical problems to the numerical routes available in the
repository. It is an implementation index, not a development roadmap or a
general survey of CFD and biomechanics.

## Status labels

| Label | Meaning |
|---|---|
| **Supported** | Implemented, documented, and exercised by repository tests or cases |
| **Functional** | A bounded workflow or reference implementation exists; use its linked validation scope |
| **Planned** | No complete repository workflow currently exists |

`FEM`, `IGA`, and `FVM` describe spatial discretizations. `Body-fitted` means
that the computational boundary follows the physical geometry. `ALE` uses a
moving volume mesh, while an immersed method represents the physical boundary
inside a background mesh. A shell or membrane is a structural surface model;
1D models discretize a path or network, and 0D models contain no spatial mesh.

## Three-dimensional models

| Physical domain | Primary route | Required geometry/data | Status |
|---|---|---|---|
| Rigid vascular flow | Body-fitted IGA Navier--Stokes or tetrahedral P2/P1 FEM | IGA control mesh and Bezier extraction, or labeled tetrahedral lumen mesh | **Supported** IGA CPU/CUDA workflow; **Functional** native FEM workflow with manufactured and idealized tube/Y evidence |
| Deforming vascular flow | Tetrahedral ALE coupled to a solid, or moving immersed flow | Fluid volume mesh plus wall/surface model and motion or material data | **Functional** ALE, moving immersed, and matching ALE--solid routes |
| Thin vessel wall | NURBS Kirchhoff--Love shell or pre-tensioned membrane | Analysis-suitable midsurface, thickness, material, supports | **Functional** single-patch shell and membrane reference runtimes |
| Thick vessel wall or local stress | Nonlinear tetrahedral solid FEM | Wall volume mesh, material, loads, and constraints | **Functional** compressible neo-Hookean and stabilized mixed P1/P1 formulations |
| Prescribed chamber flow | ALE or immersed flow with supplied wall motion | Chamber volume/surface, motion law, and labeled ports | **Functional** idealized LV comparison cases |
| Myocardium and active tissue | Nonlinear anisotropic solid with activation | Tissue volume, fibers, material, prestress, support, activation | **Planned**; current solids cover isotropic reference problems only |
| Valves and contact | Shell/solid structure coupled to flow and contact | Leaflet geometry, material, contact, and support data | **Planned** beyond the shell/membrane foundations |
| Porous tissue perfusion | P1 Darcy pressure with RT0 conservative flux | Labeled tissue tetrahedra, mobility/permeability, sources, and boundaries | **Functional** single-compartment steady workflow, including matching-face source transfer |
| Tissue deformation with fluid interaction | Poroelastic or biphasic FEM | Tissue solid/fluid material and boundary data | **Planned** |
| Species transport | Advection--diffusion--reaction in the corresponding domain | Concentration fields, velocity/flux, diffusion, reactions, sources, interfaces | **Supported** body-fitted IGA multispecies transport; **Functional** native tetrahedral P1 ALE transport and wall/reservoir exchange |

The [idealized dual-tree case](docs/CUBE_DUAL_TREE_FSI.md) combines arterial
and venous lumens, fixed Darcy tissue, and quasi-steady wall feedback. The
[geometry audit](docs/LIVER_GEOMETRY_CANDIDATES.md) documents a separate,
patient-derived local ROI workflow. Neither case supplies a calibrated
whole-organ physiology model.

## One-dimensional models

| Purpose | Formulation | Data | Status |
|---|---|---|---|
| Rigid network flow | `steady_poiseuille` | Rooted SWC or radius-annotated line OBJ | **Supported** |
| Rigid network with inertia | `rigid_inertance`, backward Euler | Fixed-area network | **Supported** |
| Compliant pulse propagation | Explicit finite-volume A/Q with Rusanov flux | Multi-cell segments and linear or Olufsen wall law | **Supported** |
| Implicit network flow | `pressure_network`, `linearized_aq`, `nonlinear_aq`, `implicit_1d_pde` | Lumped nodes/branches or multi-cell network, depending on formulation | **Supported** PETSc route |
| Species transport | Conservative `A*C` advection, diffusion, reaction, source, and wall flux | Same cell network as flow | **Supported** multispecies route |

The executable name `iga_1d` does not imply that every 1D formulation uses an
IGA basis. Inputs must include radii and a supported rooted topology. See the
[1D guide](docs/ONE_D.md) and [examples](examples/one_d/README.md).

## Zero-dimensional models

| Purpose | Model | Status |
|---|---|---|
| Terminal vascular bed | R, RC, RLC, or RCR/Windkessel with dynamic capacitor pressure | **Supported** |
| Flow source or reservoir | Storage, pump inflow, and resistance relation | **Supported** |
| VCA/extracorporeal circuit | Well-mixed reservoir, pump, and selected oxygenator/dialyzer/infusion components | **Supported** by dedicated 1D and CPU 3D paths |
| Generic species compartment | Variable-volume, well-mixed multispecies balance with staged graph adapters | **Functional** |
| Complete closed-loop heart/circulation | Time-varying elastance chambers, valves, and circulation graph | **Planned** as a general workflow |

See the [0D guide](docs/ZERO_D.md),
[`ZeroDFlowDomain`](include/ZeroDFlowDomain.hpp), and the
[generic species contract](docs/T7_GENERIC_ZERO_D_SPECIES.md).

## Geometry and mesh routes

| Input | Target | Status |
|---|---|---|
| Radius-annotated centerline | 1D network | **Supported** SWC and line-OBJ validation |
| Radius-annotated centerline | Body-fitted 3D IGA | **Supported** smoothing, hexahedral control mesh, spline/Bezier extraction, METIS, and `.ntiga` packing |
| Closed triangulated surface | Immersed IGA | **Supported** Cartesian cubic B-spline background with cut volume/surface quadrature |
| Closed triangulated surface | Body-fitted tetrahedral FEM | **Functional** Gmsh exact-boundary and fTetWild envelope-remesh adapters producing labeled Gmsh 4.1 meshes |
| Surface or spline midsurface | IGA shell | **Functional** single-patch analysis-suitable input; no automatic triangle-to-NURBS fitting |
| Organ exterior surface | Complete tissue/vessel model | **Planned**; an exterior alone does not define internal regions, ports, fibers, or materials |
| Arbitrary surface | Body-fitted volume IGA | **Planned** general patch decomposition and positive-Jacobian parameterization |

The centerline pipeline targets tubular rooted trees and supported
bifurcations. Surface-to-FEM adapters require labeled, closed geometry and do
not infer missing caps or boundary conditions.

## Coupling and execution

| Layer | Available route | Status |
|---|---|---|
| 0D/1D/3D pressure and flow | Named outward-positive SI ports; explicit or partitioned strong coupling | **Supported** for validated acyclic graph topologies |
| 1D/3D species | Flow-direction-aware concentration and conservative mass transfer | **Supported** for native 1D/body-fitted 3D; **Functional** native tetra moving chains |
| Fluid/structure | Dirichlet--Neumann iteration with Aitken relaxation; matching ALE--solid transfer | **Functional** immersed membrane and native matching reference workflows |
| Vessel/tissue | Matching-face or named-port source transfer | **Functional** conservative flow and selected species exchange workflows |
| Geometry motion | Fixed, ALE, or immersed | **Supported/Functional** by route; remeshing and general history transfer are not yet common services |
| Time integration | Backward Euler plus formulation-specific 1D explicit/implicit schemes | **Supported**; no framework-wide multirate integrator |
| CPU | C++, MPI/PETSc, and selected OpenMP assembly | **Supported** |
| GPU | Single-GPU CUDA for the body-fitted IGA subset and native FEM workflows exposed by the shared entry | **Functional** by solver/case; no general multi-GPU graph runtime |
| Checkpoint/restart | Standalone and graph-specific accepted-state formats | **Supported** where listed in [coupled restart](docs/COUPLED_RESTART.md) |
| Output | Text fields, VTU/PVTU/PVD, VTKHDF, QoI, and conservation ledgers | **Supported** by backend |

## Choosing a route

| Goal | Recommended starting point |
|---|---|
| Rigid tubular flow with a smooth spline geometry | Body-fitted 3D IGA |
| Surface-defined fixed flow | Tetrahedral FEM or immersed IGA |
| Pulse propagation over a large vascular network | Compliant 1D A/Q |
| Terminal impedance or storage | 0D R/RC/RLC/RCR |
| Prescribed moving chamber | ALE and immersed comparison workflows |
| Thin spline structure | Single-patch IGA shell foundation |
| Conservative tissue source/sink flow | Tetrahedral Darcy with RT0 recovery |
| Mixed-dimensional pressure/species coupling | Native graph runtime with named ports |

## Evidence index

- [Project overview](README.md)
- [Numerical benchmarks](docs/BENCHMARKS.md)
- [CPU solver](solvers/cpu/README.md)
- [CUDA solver](solvers/cuda/README.md)
- [0D guide](docs/ZERO_D.md)
- [1D guide](docs/ONE_D.md)
- [Moving-domain architecture](docs/architecture/MOVING_DOMAIN_ARCHITECTURE.md)
- [FSI architecture](docs/architecture/FSI_ARCHITECTURE.md)
- [Checkpoint/restart](docs/COUPLED_RESTART.md)

Update this matrix when a public entry point, supported topology, file format,
or validation status changes. Detailed measurements belong in the linked
validation documents rather than in this summary.

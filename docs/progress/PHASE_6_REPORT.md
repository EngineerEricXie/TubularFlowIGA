# Phase 6 Progress

## PR6.1 — Backend-neutral immersed 3D domain kind

Status: complete.

The coupling graph now distinguishes `three_d_immersed_flow` from the existing
body-fitted 3D backend while treating both as semantically three-dimensional.
Schema-v5 accepts immersed domains without a body-fitted `.ntiga` database;
their case directory and boundary-label ports remain the narrow configuration
contract for a later runtime factory.

Immersed ports use canonical nonnegative integer `boundary_label` locators,
outward orientation `+1`, and the current pressure/flow interface quantities.
Total-pressure and species extensions are not part of PR6.1.

The approved first Phase 6 integration target is quasi-static immersed 3D.
This milestone does not claim backward-Euler equivalence.

Validated by `simulation_graph_test`, `multidomain_config_test`, and
`pressure_flow_executor_test`: schema parsing/canonical serialization,
dimension-aware graph validation, exact runtime-kind matching, and
backend-neutral pressure/flow executor behavior.

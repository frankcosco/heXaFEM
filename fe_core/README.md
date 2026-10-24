# fe_core — generic hexahedral FE solver

Generic finite-element building blocks: mesher, Hexa8/20/27/32/64 shape
functions, Gauss integration schemes (1/8/14/27/64 points), stiffness/mass
assembly (`assmK_MFem.m`), plotting, and a node-matching utility
(`OutputRequests.m`). Nothing in this folder is specific to the aluminium
cantilever beam or to the Scientific Reports paper this repository was
assembled for.

## Why there are more element types / integration rules than this paper uses

This is a shared library, used across several of the group's publications on
higher-order finite elements, not a paper-specific solver. For the
**Scientific Reports** paper this repository accompanies ("Influence of
Interpolation Order on Finite Element Model Updating: Accuracy, Robustness,
Computational Efficiency, and Parameter Compensation Effects"), only four
element/integration combinations are ever used (see Table 2 of the
manuscript and `model_updating/execution/Modal_updating_clean_LM.m`,
`EXA_types`/`EXA_orders`):

| `hexa_type` (nodes) | `hexa_order` (Gauss points) | Used as |
|---|---|---|
| 8  | 8  | Hexa8-Int8  |
| 20 | 27 | Hexa20-Int27 |
| 27 | 27 | Hexa27-Int27 |
| 64 | 64 | Hexa64-Int64 |

`HEXA_shapes_32.m` (a 32-node element) and `HEXA_integration_14.m` (a
14-point Gauss rule) are reachable from `assmK_MFem.m`'s dispatch switches
but are **never selected by any script in `model_updating/`** for this
paper — confirmed by checking every value `EXA_types`/`EXA_orders` take in
`Modal_updating_clean_LM.m` and in every script under
`model_updating/sensitivity_runs/`. They are not dead code in an absolute
sense, though: this same `fe_core` library has been used in the group's
other higher-order finite-element work, where these configurations are
relevant:

- Karpik, A., Cosco, F. & Mundo, D. "Higher-Order Hexahedral Finite Elements
  for Structural Dynamics: A Comparative Review." *Machines* **11**(3), 326
  (2023). <https://doi.org/10.3390/machines11030326>
- Karpik, A., Cosco, F. & Mundo, D. "Dynamic Error Estimation in
  Higher-Order Finite Elements." *Designs* **8**(4), 79 (2024).
  <https://doi.org/10.3390/designs8040079>

They are kept here, undocumented-per-configuration beyond this note, rather
than removed, so that this repository stays usable as the one shared
`fe_core` across the group's related publications instead of forking a
paper-specific copy. If you are only trying to reproduce the Scientific
Reports results, you only need the four combinations in the table above.

## Files

- `assmK_MFem.m` — global stiffness/mass assembly; dispatches on `hexa_type`
  (element shape functions) and `hexa_order` (Gauss rule).
- `HEXA_shapes_8/20/27/32/64.m` — isoparametric shape functions and
  derivatives for each element.
- `HEXA_integration_1/8/14/27/64.m` — Gauss point coordinates and weights
  for each integration rule.
- `stiffnessMatrix_HEXA_new.m` — alternative/earlier stiffness-matrix
  routine, kept alongside `assmK_MFem.m` for reference (see its own header
  comments for which one is live).
- `brick_hexa_mesher.m` / `brick_hexa8_mesher_advanced.m` — structured brick
  mesh generators for the cantilever-beam geometry, parametrized by element
  type and per-axis element counts.
- `OutputRequests.m` — nearest-node matching between two node sets (used to
  snap a design mesh's boundary nodes onto experimentally-measured DIC grid
  coordinates).
- `FEplot.m` — mesh/result plotting utility.

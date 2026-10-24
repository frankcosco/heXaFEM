# heXaFEM

**heXa**hedral finite elements, e**X**tended beyond plain FE: a small,
paper-specific toolbox for higher-order hexahedral elements plus a family
of modules built on top of it. The first module is **heXaFEM-MU**, a
finite-element model-updating (FEMU) pipeline.

> **Scope note**: this is not a general-purpose FE framework (for that,
> see established libraries such as [GetFEM](https://getfem.org/),
> [deal.II](https://www.dealii.org/), or [Elmer](https://www.elmerfem.org/)).
> It is the reproducibility companion code for the three publications
> listed below: it implements only the hexahedral element families those
> papers use, and the model-updating module only the calibration pipeline
> they describe.

## What's in here

- **Core library** (`fe_core/`): a generic hexahedral FE solver — mesher,
  Hexa8/20/27/32/64 shape functions, Gauss integration rules,
  stiffness/mass assembly, plotting.
- **heXaFEM-MU** (`model_updating/`): the model-updating pipeline that
  calibrates material parameters (density, Young's modulus, Poisson's
  ratio) of an FE model against experimental modal data (frequencies +
  full-field DIC mode shapes), plus a set of scripts quantifying the
  identifiability and robustness of that calibration across element
  types/orders and mesh densities.
- **Experimental-data readers** (`experimental_data/`) for the test case:
  UFF modal-test exports and DIC-based mode-shape re-alignment.

See `docs/VALIDATION_NOTES.md` for the detailed build/validation log, and
each subfolder's own `README.md` for more on that part specifically
(`fe_core/README.md` in particular explains which element/integration
combinations are used where across the related publications below).

## Related publications

This repository is the companion code for three publications on higher-order
hexahedral finite elements and their use in model updating:

1. **A. Karpik, F. Cosco, D. Mundo (2023).** *Higher-Order Hexahedral Finite
   Elements for Structural Dynamics: A Comparative Review.* Machines, 11(3).
   https://doi.org/10.3390/machines11030326
   Systematically compares the common hexahedral element formulations
   (Hexa8, serendipity/complete Hexa20/27, Hexa64), including reduced
   integration schemes, and benchmarks their convergence rate and
   computational efficiency on the structural eigenvalue problem.

2. **A. Karpik, F. Cosco, D. Mundo (2024).** *Dynamic Error Estimation in
   Higher-Order Finite Elements.* Designs, 8(4), 79. https://doi.org/10.3390/designs8040079
   Extends the comparison to transient dynamic solutions, and shows that
   the eigenenergy norm computed from the eigensolution is a reliable,
   much cheaper predictor of the time- and frequency-domain accuracy
   (TRAC/FRAC) that would otherwise require a full transient run.

3. **F. Cosco, A. Karpik, M. Perrelli, D. Mundo.** *Influence of
   Interpolation Order on Finite Element Model Updating: Accuracy,
   Robustness, Computational Efficiency, and Parameter Compensation
   Effects.* Submitted to Scientific Reports (under review).
   Uses the element library from (1)-(2) inside a nonlinear least-squares
   model-updating pipeline (`heXaFEM-MU`) to investigate how finite-element
   interpolation order affects the accuracy, robustness and physical
   plausibility of the identified material parameters. The paper proposes
   normalized compensation indices (PCI₀, PCI∞) to quantify the extent to
   which a reduced modal residual is accompanied by a genuine improvement
   in parameter identification, and evaluates them on both numerical and
   DIC-based experimental modal data from an aluminium cantilever beam.
   Manuscript currently under peer review; details may change before
   publication.

## Code authors

This codebase originated from, and grew out of, the work of **Anna
Karpik** and **Francesco Cosco** during Anna Karpik's PhD research,
carried out under the continuous supervision of Francesco Cosco. The
library is currently maintained by Francesco Cosco.

This is the authorship of the *code*, not of the papers: the full author
lists of the associated publications are given under "Related
publications" above and are unaffected by this.

## Reproducing the Scientific Reports results

1. MATLAB, Optimization Toolbox required.
2. Current Folder = `model_updating/execution/`, run
   `Modal_updating_clean_LM.m` for the full mesh-convergence campaign
   (Table 2 of the paper), or:
3. Current Folder = `model_updating/sensitivity_runs/`, run
   `compute_and_cache_experimental_data.m` once, then any of the
   identifiability/robustness scripts in that folder independently — see
   `docs/VALIDATION_NOTES.md` for what each one tests.

`model_updating/reference_results/` holds the exact log files behind the
paper's published Figures 7–11, for direct comparison.

## Data

The raw DIC calibration images and raw modal-test export used by the code
are not included in this public release (see `.gitignore`); the
already-extracted modal data (frequencies, mode shapes) are published in
tabular form instead, in
[`experimental_data/tabular_modal_data/`](experimental_data/tabular_modal_data/README.md).


# Validation & build notes

Companion code for "Influence of Interpolation Order on Finite Element
Model Updating: Accuracy, Robustness, Computational Efficiency, and
Parameter Compensation Effects" (submitted to Scientific Reports).

> **Naming note**: the repository is named `heXaFEM`, with `heXaFEM-MU`
> as the model-updating module built on top of it — chosen to leave room
> for future modules beyond model updating.
>
> **Margin note for later**: it may be worth evaluating, post-publication,
> whether to migrate this codebase onto an established, general-purpose FE
> library (e.g. GetFEM, deal.II, or Elmer) instead of maintaining
> `fe_core/` as bespoke MATLAB code. Those libraries already provide
> higher-order hexahedral elements, generic assembly, and (GetFEM in
> particular) a generic weak-form language, which could remove a lot of
> the maintenance burden `fe_core/` currently carries. This is a
> non-trivial undertaking (different language/ecosystem — C++/Python
> rather than MATLAB — for most of these) that would need to be scoped
> properly before starting.

## Repository structure

- `fe_core/` — generic FE solver: mesher, Hexa8/20/27/32/64 shape
  functions, Gauss integration schemes, stiffness/mass assembly, plotting,
  node-matching utility. Shared with other publications — see
  `fe_core/README.md` for which element/integration combinations this
  paper actually uses, and the related citations.
- `experimental_data/` — the single experimental test case: UFF/.unv
  readers and DIC node re-alignment (`riallineaDICtests.m`, with its
  homography helpers). See "Data" below for what is and isn't included.
- `model_updating/execution/` — runs the paper's single-case calibration:
  `Modal_updating_clean_LM.m` (entry point), `objectiveFunction.m`,
  `outFun.m` (lsqnonlin output function / logger), the MAC-based error
  metric.
- `model_updating/post_processing/` — reads already-saved `.txt` result
  logs and produces the paper's plots/Tikz coordinates — does NOT run any
  FE analysis. `parseModelUpdatingReport.m`, `parseModelUpdatingFolder.m`,
  `parseFinelPaper.m`, `buildTikzCoordinatesFromTable.m`.
- `model_updating/reference_results/` — the 13 `.txt` log files that are
  the exact outputs behind the paper's published Figures 7–11, i.e. the
  mesh-convergence DoE of Table 2 (5 Hexa8-Int8 densities, 3 Hexa20-Int27,
  3 Hexa27-Int27, 2 Hexa64-Int64).
- `model_updating/sensitivity_runs/` — the numerical identifiability and
  robustness evidence discussed in the paper (see below for the full list
  of scripts).
- `model_updating/logs/` — run outputs only (no code): `.txt` diary logs
  from `execution/`, and `.csv`/`.mat` results from `sensitivity_runs/`.

Any script under `model_updating/*` that opens a data file by bare
relative filename resolves it against MATLAB's Current Folder, so each
subfolder's script can be run standalone with Current Folder set to that
subfolder.

## Data

The raw DIC calibration images and the raw modal-test export (`.unv`) are
not included in this public release. The already-extracted modal data
(frequencies, mode shapes) used by the paper are published instead, in
tabular form, in `experimental_data/tabular_modal_data/`.

## Identifiability / robustness evidence (`sensitivity_runs/`)

Numerical campaign behind the paper's discussion of whether jointly
fitting (E, ρ, ν) is methodologically sound, and of the
discretization-dependent "Parameter Compensation Effect" the paper
reports. Run `compute_and_cache_experimental_data.m` once first (produces
`aligned_experimental_data.mat`, shared by every other script below);
then any of the others can be run independently, each with Current Folder
set to this subfolder.

| Script | Question it answers | Results saved to |
|---|---|---|
| `compute_and_cache_experimental_data.m` | Caches the DIC-aligned experimental data once, so every sensitivity run below uses an identical reference (isolates optimizer sensitivity from DIC-fit run-to-run noise). | `aligned_experimental_data.mat` (this folder) |
| `multi_start_hexa8.m` | Does Hexa8-Int8 (26x234x2) converge to the same compensated (ρ,E,ν) region regardless of starting point, or is it an artifact of always starting from the nominal x0? 5 trials (nominal, solver default, Hexa64 asymptotic reference, both bound corners). | `../logs/multi_start_hexa8_results.csv` / `.mat` |
| `jacobian_svd_hexa8.m` | SVD + correlation matrix of the residual Jacobian at the Hexa8-Int8 (26x234x2) optimum — which parameter combination is poorly constrained by the modal data. | `../logs/jacobian_svd_hexa8_results.mat` |
| `jacobian_svd_hexa64_small.m` | Same SVD/correlation analysis on a different element type/order (Hexa64-Int64, 13x117x1) — checks the flat direction is a property of the parametrization, not one specific mesh. | `../logs/jacobian_svd_hexa64_small_results.mat` |
| `verify_rho_E_gauge_invariance.m` | Direct empirical check of the exact ρ–E gauge invariance on the real FE+experimental pipeline: rescaling a converged optimum to any reference ρ reproduces the identical residual (no re-optimization). | `../logs/gauge_invariance_check_results.csv` |
| `jacobian_svd_hexa27_fixedRho_freeL.m` | Follow-up: with ρ fixed (not fitted) and beam length freed instead, is there a different/new flat direction (E vs L), and does fixing ρ make E, ν well identified? (Hexa27-Int27, 13x117x1.) | `../logs/jacobian_svd_hexa27_fixedRho_freeL_results.mat` |
| `jacobian_svd_three_meshes_fixedRho_fixedL.m` | Diagnostic only (not the paper's methodology, which leaves ρ free per mesh): with both ρ and length fixed at nominal, repeats the SVD/correlation check on Hexa8-Int8, Hexa27-Int27 and Hexa64-Int64, to quantify how much fixing both degenerate directions improves conditioning, and whether an externally-fixed ρ distorts the fit differently per mesh/order. | `../logs/jacobian_svd_three_meshes_fixedRho_fixedL_summary.csv` / `_results.mat` |
| `jacobian_svd_hexa27_fixedRho_fixedL_retry.m` | Standalone retry of the Hexa27-Int27 case above with the ν bound raised to 0.45 (the three-mesh run clipped ν at its 0.40 upper bound). | `../logs/jacobian_svd_hexa27_fixedRho_fixedL_retry_results.mat` |
| `jacobian_svd_along_trajectory_hexa8.m` | Is the ρ–E ill-conditioning present throughout the optimization trajectory, or only near the final converged point? Captures the Hexa8-Int8 iterate trajectory via `OutputFcn`, recomputes an independent central-difference Jacobian at the initial guess, an intermediate iterate, and the final optimum. | `../logs/jacobian_svd_along_trajectory_hexa8_summary.csv` |

**Headline findings**: multi-start converges to the same region from 5
different starting points (ρ/E spread ~1%, ν spread ~0.5%); the ρ–E gauge
invariance holds to within ~10⁻⁸ on the real pipeline; SVD condition
numbers are 170–420 depending on mesh/order with corr(ρ,E) ≈
0.9999–1.0000; the flattest Jacobian direction matches the theoretically
predicted exact-invariance direction to within ~0.06° at every point
along the optimization trajectory, not only at convergence.

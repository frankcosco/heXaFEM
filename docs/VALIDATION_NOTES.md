# Internal validation & build notes

> This is the internal working log kept during the 23-24/10/2026 cleanup and numerical-verification pass. The public-facing landing page is the top-level `README.md`; this file has the detailed, session-by-session validation history, bug fixes, and status notes that don't belong on a public landing page but are useful to keep for provenance.

> **Naming note (24/10/2026)**: the repository is named `heXaFEM`, with
> `heXaFEM-MU` as the model-updating module built on top of it (chosen
> over `HexaFEMU` to leave room for future modules beyond model updating).
>
> **Margin note for later (post-publication, not for this revision)**:
> once the paper is through review, it would be worth spending some time
> evaluating whether to migrate this codebase onto an established,
> general-purpose FE library (e.g. GetFEM, deal.II, or Elmer) instead of
> maintaining `fe_core/` as bespoke MATLAB code. Those libraries already
> provide higher-order hexahedral elements, generic assembly, and (GetFEM
> in particular) a generic weak-form language, which could remove a lot
> of the maintenance burden `fe_core/` currently carries and make the
> model-updating pipeline easier for others to extend. This is explicitly
> deferred — not relevant to the current Scientific Reports revision, and
> a non-trivial undertaking (different language/ecosystem for most of
> these: C++/Python rather than MATLAB) that should be scoped properly
> before starting, not squeezed in alongside the revision.

# Clean reference copy — FE model updating (Scientific Reports revision)

Companion code for: "Influence of Interpolation Order on Finite Element
Model Updating: Accuracy, Robustness, Computational Efficiency, and
Parameter Compensation Effects" (Scientific Reports, submission
95ce1094-6dfb-4062-8b70-8e68b17806a0).

Built 23/10/2026 from the original working copy in
`dev2\__final__ModelUpdatingBeam__AnnaPHD` (and a few shared FE-core files
living directly in `dev2\`), which contains many redundant/historical
variants. This folder keeps only the files actually needed to run the
pipeline and reproduce the paper's results, organized by role.

## Repository structure

- `fe_core/` — generic FE solver: mesher, Hexa8/20/27/32/64 shape
  functions, Gauss integration schemes, stiffness/mass assembly, plotting,
  node-matching utility. Shared with other publications of the group — see
  `fe_core/README.md` for which element/integration combinations this paper
  actually uses, and the citations to the other work.
- `experimental_data/` — the single experimental test case: UFF/.unv
  readers, DIC node re-alignment (`riallineaDICtests.m`, with its
  homography helpers), and the raw data files (`real_modeset.unv`, two
  reference `.tif` images — see the Data note below).
- `model_updating/execution/` — runs the paper's single-case calibration:
  `Modal_updating_clean_LM.m` (entry point), `objectiveFunction.m`,
  `outFun.m` (lsqnonlin output function / logger), the MAC-based error
  metric.
- `model_updating/post_processing/` — reads already-saved `.txt` result
  logs and produces the paper's plots/Tikz coordinates — does NOT run any
  FE analysis. `parseModelUpdatingReport.m`, `parseModelUpdatingFolder.m`,
  `parseFinelPaper.m`, `buildTikzCoordinatesFromTable.m`.
- `model_updating/reference_results/` — the 13 `.txt` log files (copied
  verbatim, not regenerated) that are the exact outputs behind the paper's
  published Figures 7–11, i.e. the mesh-convergence DoE of Table 2 (5
  Hexa8-Int8 densities, 3 Hexa20-Int27, 3 Hexa27-Int27, 2 Hexa64-Int64).
  See that folder's own `README.md` for a provenance note on where these
  came from and a discrepancy that was investigated and resolved on
  23/10/2026.
- `model_updating/sensitivity_runs/` — the numerical identifiability /
  robustness evidence behind the response to Reviewers 3, 4 and 5 (see
  below for the full list of scripts).
- `model_updating/logs/` — run OUTPUTS only (no code): the `.txt` diary
  logs from `execution/`, and the `.csv`/`.mat` results from
  `sensitivity_runs/`. Kept separate so the code folders stay clean.

Any script under `model_updating/*` that opens a data file by bare relative
filename (`real_modeset.unv`, the two `.tif`s — fopen/imread resolve
against MATLAB's Current Folder, not its path) keeps its own working copy
of those 3 files alongside it, so each subfolder's script can be run
standalone with Current Folder set to that subfolder.

## Data note — what is / is not published alongside the code

The raw DIC calibration images (`referenceDIC-FE-UPDATING.tif`,
`referenceDic-grid.tif`) and the raw modal test export
(`real_modeset.unv`) are present in this working copy (in
`experimental_data/`, `model_updating/execution/`, and
`model_updating/sensitivity_runs/`) because the code needs them to run
end-to-end. **Decision (23/10/2026): these raw files are not meant to go
into the public GitHub/Zenodo release.** Re-deriving results from them
would require re-running the DIC homography fit from scratch, which is not
necessary for reproducing the paper's reported numbers and is not what a
reader of the paper needs. Instead, the already-extracted modal data
(frequencies, mode shapes) should be published in tabular form, the same
way `model_updating/reference_results/` already does for the FEMU outputs.
**Action item before pushing a public release**: exclude the three raw
files above (in all three locations) from what gets published — e.g. via
`.gitignore` or by not including them in the release archive — and confirm
where the already-extracted tabular experimental data should live in this
repository (not yet decided — see the project tracker).

## Status — numerically validated (23/10/2026)

`model_updating/execution/Modal_updating_clean_LM.m` has ONE correction
applied relative to the live script in `dev2`: `var2optimize` was
`[true true true true false false]` there (thickness free), but every
saved result in `revisedFrequencies/materialOnly/` reports the optimized
thickness unchanged from its initial 4.0 mm across all cases — proof the
historical "material only" runs (consistent with the paper's Fig. 8, which
reports only E, rho, nu) were run with thickness FIXED. Corrected here to
`[true true true false false false]`.

An initial validation run (23/10/2026, Hexa8-Int8 26x234x2) did not
exactly reproduce the historical `materialOnly` log byte-for-byte: the
INITIAL error (before any optimizer iteration, same x0, same mesh) already
differed. Root cause identified: `riallineaDICtests.m` runs its own
internal `lsqnonlin` (Levenberg-Marquardt, in `fineTuneHomography`/
`residui_vincolati`) to refine the DIC pixel→physical homography — a
separate nonlinear fit from the main model-updating optimization, and
apparently not perfectly reproducible run-to-run (likely
MATLAB/Optimization-Toolbox version drift between April and October 2026).
This does not affect the qualitative conclusions below, but means a fresh
run is not expected to be byte-for-byte identical to the paper's original
logs. `sensitivity_runs/compute_and_cache_experimental_data.m` works
around this by computing the DIC alignment ONCE and caching it
(`aligned_experimental_data.mat`), so that all `sensitivity_runs/` scripts
use the identical experimental reference data and are directly comparable
to each other.

## Identifiability / robustness evidence (`sensitivity_runs/`)

Numerical campaign behind the response to Reviewers 3, 4 and 5 on whether
jointly fitting (E, ρ, ν) is methodologically sound, and on the
discretization-dependent "Parameter Compensation Effect" the paper
reports. Run `compute_and_cache_experimental_data.m` once first (produces
`aligned_experimental_data.mat`, shared by every other script below); then
any of the others can be run independently, each with Current Folder set
to this subfolder.

| Script | Question it answers | Results saved to |
|---|---|---|
| `compute_and_cache_experimental_data.m` | Caches the DIC-aligned experimental data once, so every sensitivity run below uses an identical reference (isolates optimizer sensitivity from DIC-fit run-to-run noise). | `aligned_experimental_data.mat` (this folder) |
| `multi_start_hexa8.m` | R3.1: does Hexa8-Int8 (26x234x2) converge to the same compensated (ρ,E,ν) region regardless of starting point, or is it an artifact of always starting from the nominal x0? 5 trials (paper nominal, script nominal, Hexa64 asymptotic reference, both bound corners). | `../logs/multi_start_hexa8_results.csv` / `.mat` |
| `jacobian_svd_hexa8.m` | R3.3/R4.2/R5.3: SVD + correlation matrix of the residual Jacobian at the Hexa8-Int8 (26x234x2) optimum — which parameter combination is poorly constrained by the modal data. | `../logs/jacobian_svd_hexa8_results.mat` |
| `jacobian_svd_hexa64_small.m` | Same SVD/correlation analysis on a different element type/order (Hexa64-Int64, 13x117x1) — checks the flat direction is a property of the parametrization, not one specific mesh. | `../logs/jacobian_svd_hexa64_small_results.mat` |
| `verify_rho_E_gauge_invariance.m` | Direct empirical check of the exact ρ–E gauge invariance on the real FE+experimental pipeline: rescaling a converged optimum to any reference ρ reproduces the identical residual (no re-optimization). | `../logs/gauge_invariance_check_results.csv` |
| `jacobian_svd_hexa27_fixedRho_freeL.m` | Follow-up: with ρ fixed (not fitted) and beam length freed instead, is there a different/new flat direction (E vs L), and does fixing ρ make E, ν well identified? (Hexa27-Int27, 13x117x1.) | `../logs/jacobian_svd_hexa27_fixedRho_freeL_results.mat` |
| `jacobian_svd_three_meshes_fixedRho_fixedL.m` | With BOTH ρ and length fixed at nominal (only E, ν free), repeats the SVD/correlation check on Hexa8-Int8, Hexa27-Int27 and Hexa64-Int64, to quantify how much fixing both known degenerate directions improves conditioning — and whether forcing an externally-fixed ρ distorts the fit differently per mesh/order (diagnostic test, not a recommended methodology: the paper's own per-mesh free-ρ fit is the correct choice — see `verify_rho_E_gauge_invariance.m` and the project tracker for why). | `../logs/jacobian_svd_three_meshes_fixedRho_fixedL_summary.csv` / `_results.mat` |
| `jacobian_svd_hexa27_fixedRho_fixedL_retry.m` | Standalone retry of the Hexa27-Int27 case above: the three-mesh run clipped ν at its 0.40 upper bound (not an interior optimum), so this reruns it alone with the ν bound raised to 0.45. | `../logs/jacobian_svd_hexa27_fixedRho_fixedL_retry_results.mat` |
| `jacobian_svd_along_trajectory_hexa8.m` | Is the ρ–E ill-conditioning present throughout the optimization trajectory (as the exact invariance argument predicts), or only near the final converged point? Captures the Hexa8-Int8 iterate trajectory via `OutputFcn`, recomputes an independent central-difference Jacobian at the initial guess, an intermediate iterate, and the final optimum. | `../logs/jacobian_svd_along_trajectory_hexa8_summary.csv` |

**Headline findings** (full narrative and numbers in the response to
Reviewers 3/4/5 and in the project tracker, not duplicated here):
multi-start converges to the same region from 5 different starting points
(ρ/E spread ~1%, ν spread ~0.5%); the ρ–E gauge invariance holds to within
~10⁻⁸ on the real pipeline; SVD condition numbers are 170–420 depending on
mesh/order with corr(ρ,E) ≈ 0.9999–1.0000; the flattest Jacobian direction
matches the theoretically predicted exact-invariance direction to within
~0.06° at every point along the optimization trajectory, not only at
convergence.

## Not yet done

- Decide and implement where the tabular (non-image) experimental data
  should live for the public release, and exclude the raw DIC/UFF files
  from the published archive (see Data note above).
- A LICENSE file.
- Final check that every `sensitivity_runs/` script, as currently saved,
  reproduces the numbers cited in the reviewer response if re-run from
  scratch (one inconsistency — `jacobian_svd_hexa64_small.m` had `ρ`
  incorrectly fixed instead of free — was found and corrected on
  23–24/10/2026; the rest were checked and are consistent with their own
  saved results).

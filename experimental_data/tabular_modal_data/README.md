# Tabular experimental modal data

Already-extracted experimental modal data for the aluminium cantilever
beam test case, exported in plain tabular form — generated 24/10/2026 from
`model_updating/sensitivity_runs/aligned_experimental_data.mat` (the cached
output of the UFF read + DIC re-alignment pipeline, see
`model_updating/sensitivity_runs/compute_and_cache_experimental_data.m`).

**This is the public substitute for the raw DIC calibration images and raw
UFF export** (decision, 23/10/2026): re-deriving these numbers from the raw
data would require re-running the DIC homography fit from scratch, which
is not needed to reproduce the paper's reported model-updating results —
only this already-extracted modal data is.

## Files

- `experimental_modal_frequencies.csv` — the 4 identified natural
  frequencies (Hz) used as model-updating targets, in the order they are
  stored throughout the pipeline (`ExpModal_freqs` / `ExpLambda`): 1st
  bending, 2nd bending, torsional, 3rd bending.
- `experimental_mode_shapes.csv` — one row per DIC-measured grid node
  (1115 nodes): its DIC-realigned physical coordinates (`x_mm`, `y_mm`,
  `z_mm`), followed by the 3D displacement components (`ux`, `uy`, `uz`,
  dimensionless mode-shape amplitude) for each of the 4 modes above
  (`mode1_*`matches frequency row 1, etc.).

## Provenance notes

- Node coordinates are the DIC-realigned grid (`ExpMNodes` in the cached
  `.mat`), i.e. already matched onto the FE mesh's coordinate system by
  `riallineaDICtests.m` — not the raw pixel/camera coordinates.
- The mode-shape sign convention here already includes the sign flips
  applied in `compute_and_cache_experimental_data.m` (modes 2 and 3 are
  negated relative to the raw UFF read) — use these values as-is, do not
  re-flip.
- This cache was computed once (23/10/2026) and reused across the whole
  identifiability/robustness study in `model_updating/sensitivity_runs/`,
  so these numbers are exactly what every script in that folder was run
  against.

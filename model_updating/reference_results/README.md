# Reference results (materialOnly)

These 13 files are copied verbatim from the historical `dev2` repository
(`__final__ModelUpdatingBeam__AnnaPHD/revisedFrequencies/materialOnly/`),
copied on 23/10/2026 and **not regenerated** -- they are the exact outputs
behind the manuscript's published Figures 7, 8, 9, 10 and 11.

They are read by `post_processing/parseFinelPaper.m`
(`parseModelUpdatingFolder("revisedFrequencies\materialOnly", ...)` in the
original script -- the path there assumes this folder is still named
`revisedFrequencies\materialOnly`; update the call if you keep the new
`reference_results` name).

"materialOnly" = density (rho), Young's modulus (E) and Poisson's ratio (nu)
free; thickness/width/length FIXED at nominal (4, 40, 360 mm) -- this is the
configuration actually used for the paper's headline results (NOT the
`allParams`/`allParams2` variants, which also freed geometry and were
superseded; those are not copied here, see `dev2` if needed).

## Provenance note (added 23/10/2026): Hexa8-Int8_65x585x5

Table 2 lists 5 Hexa8-Int8 densities (26x234x2, 39x351x3, 52x468x4,
65x585x5, 78x702x6). The `materialOnly/` folder on `dev2` only had 4 of
these active at the top level; `Hexa8-Int8_65x585x5.txt` was sitting in a
`materialOnly/discarded/` subfolder instead, which `parseModelUpdatingFolder`
does not recurse into.

Two hypotheses were checked:

1. **"The root `revisedFrequencies/` folder (not `materialOnly/`) is the
   real source of the figures."** Ruled out: diffed byte-for-byte, root's
   Hexa8 files and `materialOnly`'s are genuinely different optimization
   runs (same x0, but different initial residuals `Err_PI` -- e.g.
   `[0.1234 0.1372 0.2390 0.1678]` in root vs `[0.1243 0.1367 0.1839
   0.1885]` in materialOnly for 26x234x2), and root's files all predate
   materialOnly's by several hours on the same day (21/04/2026, root ~01:05-
   03:46, materialOnly ~11:10-19:00). Root is an earlier/superseded run;
   `materialOnly` is the later, corrected one that the live code actually
   reads from (confirmed both by the `parseFinelPaper.m` source and by
   these timestamps).

2. **"`Hexa8-Int8_65x585x5.txt` in `discarded/` was actually part of the
   same materialOnly batch and belongs in the 5-point series."**
   Confirmed: its timestamp (21/04 15:15) sits right inside the window of
   the other active `materialOnly` files (11:10-19:00), and its raw
   optimum (rho=2.6644e-06, E=6.5341e7 kPa, nu=0.3593) only LOOKS like an
   outlier (rho and E both dip below their neighbours) because it landed
   on a different point of the same rho-E degenerate manifold documented
   elsewhere in this revision. Rescaling it to the common reference
   rho_ref=2.7747e-06 (E_new = E*(rho_ref/rho), the same gauge-fixing used
   in the identifiability analysis) gives E=6.805e7 kPa -- which falls
   exactly in sequence with the other four points' rescaled E (6.348,
   6.608, 6.717, 6.805, 6.860 x1e7 kPa, monotonically increasing with mesh
   refinement). So the file is numerically sound and consistent; there is
   no technical reason found so far to exclude it, and it has been added
   back here as the 13th file to complete Table 2's DoE.

Do not edit these files -- if the pipeline is rerun, save new results
elsewhere (e.g. `model_updating/logs/`) and compare against these as the
reference, rather than overwriting them.

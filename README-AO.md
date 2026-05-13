# Analysis Extension — Branch `aoliveram/analysis-extensions`

This branch extends the e-cigarette threshold analysis on the USC ADVANCE
dataset documented in `netdiffuseR_past6mo_report.docx`. Yuchan's pipeline
(`codes/260427_netdiffuseR_past6mo 4.17.35 PM.R`) is left untouched; new work
lives in `codes/NN-*-AO.R` and outputs in `outputs_AO/`.

## TL;DR

1. **Replication** of Tables 2–3 from the report on the refreshed
   (post-2026-04-28) wave data succeeds. Coefficients shift slightly
   because C2 schools gained ~47 adopters in the data refresh; substantive
   conclusions are preserved.
2. **Diagnostic finding**: median ego out-degree at TOA is **2**, and
   **26.6% of adopters are network isolates** at their adoption time
   (degree = 0 in the consenting friendship network). The threshold
   variable is therefore effectively discrete with mass on rational
   fractions {0, ½, ⅓, ⅔, 1}, and Yuchan's quasibinomial GLM treats this
   discrete support as a continuous proportion.
3. **Methodological extension**: refit the threshold regression as a
   **grouped binomial** `glm(cbind(k_users, n_alters − k_users) ~ X,
   family = binomial)`, dropping isolates. This properly weights each
   adopter by their degree and ignores the uninformative isolate stratum.
4. **Robustness check**: the Hispanic threshold-paradox effect
   (`OR = 1.64 → 1.48`, still significant) and the perceived-friend-use
   effect (`OR = 1.14 → 1.24`, now `p < .001`) both **survive the
   methodological correction**. GAD/MDD remain non-significant. Sex-min
   coefficient changes direction (becomes marginally negative,
   `OR = 0.71, p = .093`), worth a follow-up.

Includes a small codebook-aligned fix in the `female` coding
(`wN_DEM_GENDER == 0` per the ADVANCE codebook).

## File map

```
codes/
  01-rebuild-threshold-data-AO.R         # rebuild thr_sub from raw waves
  02-threshold-empirical-plots-AO.R      # boxplot+violin by covariate
  03-threshold-model-predictions-AO.R    # marginal predictions, Yuchan-spec
  04-degree-diagnostics-AO.R             # ego degree-at-TOA distribution
  05-grouped-binomial-AO.R               # 3-model coef comparison
  06-grouped-binomial-predictions-AO.R   # marginal predictions, Model C

outputs_AO/
  intermediate/
    thr_data-AO.rds                      # thr_sub + diffnet objects
    thr_data_with_kn-AO.rds              # + n_alters and k_users per adopter
  empirical/
    threshold_empirical_by_covariate-AO.pdf
    threshold_empirical_facetted_by_toa-AO.pdf
    threshold_empirical_summary-AO.csv
  diagnostics/
    degree_diagnostics-AO.pdf
    degree_summary-AO.csv
  model/
    threshold_model_predictions_{GAD,MDD}-AO.pdf      # Yuchan-spec marginals
    threshold_model_predictions-AO.csv
    threshold_model_summaries-AO.txt
    grouped_binomial_coefficients-AO.csv
    grouped_binomial_forestplot-AO.pdf                # A vs B vs C
    grouped_binomial_predictions_{GAD,MDD}-AO.pdf
    grouped_binomial_predictions-AO.csv
    grouped_binomial_summaries-AO.txt
```

## Threshold regression — comparison across specifications

Reference: Yuchan baseline = `glm(threshold ~ X, family = quasibinomial)`
including isolates with threshold = 0 (per `netdiffuseR::threshold()`).

| Parameter | A. Yuchan baseline (n = 481) | C. Grouped binomial (n = 404) |
|-----------|------------------------------|-------------------------------|
| **Hispanic** | OR = 1.64 [1.14, 2.37], *p* = .008 | **OR = 1.48 [1.09, 2.03], *p* = .013** |
| **Friends_ecig** (PFU) | OR = 1.14 [1.03, 1.25], *p* = .009 | **OR = 1.24 [1.15, 1.35], *p* < .001** |
| GAD | OR = 0.98 [0.80, 1.20], NS | OR = 1.00 [0.84, 1.19], NS |
| MDD | OR = 0.86 [0.69, 1.08], NS | OR = 0.90 [0.72, 1.13], NS |
| Asian | OR = 0.96, NS | OR = 0.99, NS |
| Female (post-fix) | OR = 1.11, NS | OR = 1.06, NS |
| **Sex-min** | OR = 1.10, NS | **OR = 0.71 [0.47, 1.06], *p* = .093** |
| Cohort C2 | OR = 0.71, *p* = .10 | OR = 0.78, NS |

See `outputs_AO/model/grouped_binomial_forestplot-AO.pdf` for the visual.

## Degree diagnostics — why we re-specified

| | Overall | Asian-maj. | Hispanic-maj. | Other |
|---|---:|---:|---:|---:|
| Adopters | 828 | 297 | 363 | 168 |
| Median ego degree at TOA | 2 | 2 | 2 | 2 |
| % with degree = 0 (isolates) | **26.6** | 25.3 | 26.4 | 29.2 |
| % with degree ≤ 1 | 39.2 | 35.4 | 41.1 | 42.3 |
| % with degree ≤ 3 | 69.8 | 62.6 | 73.3 | 75.0 |

See `outputs_AO/diagnostics/degree_diagnostics-AO.pdf` for the visual.
Report §9 already notes the consenting-participants restriction as a
limitation; the 26.6% degree-0 rate quantifies its impact.

## How to reproduce

Requires the ADVANCE raw waves locally. The default `data_path` in
`01-rebuild-threshold-data-AO.R` points to the companion disadoption
project's `data/advance/Cleaned-Data/`; adjust if running elsewhere.

```bash
Rscript codes/01-rebuild-threshold-data-AO.R
Rscript codes/02-threshold-empirical-plots-AO.R
Rscript codes/03-threshold-model-predictions-AO.R
Rscript codes/04-degree-diagnostics-AO.R
Rscript codes/05-grouped-binomial-AO.R
Rscript codes/06-grouped-binomial-predictions-AO.R
```

R packages used beyond `netdiffuseR`: `ggplot2`, `ggpubr`, `dplyr`,
`tidyr`, `patchwork`, `marginaleffects`, `openxlsx`.

## Open issues / next steps

- **School-level clustering** (report §9, limitation #4): refit Model C
  with random effects for school (`glmer` / `glmmTMB`).
- **Cross-study comparison** with the companion disadoption analysis
  (Olivera et al., 2026): same predictors, opposite directionality —
  joint table strengthening the bidirectional-mechanism claim.
- **Sex-min direction flip**: the change of sign under the corrected
  specification is unexpected; warrants a focused diagnostic.

---
Author: Anibal Olivera. Branch opened 2026-05-12.

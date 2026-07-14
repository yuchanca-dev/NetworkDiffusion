# 05-grouped-binomial-AO.R
# Compare three specifications of the threshold regression to test whether
# Yuchan's quasibinomial is biased by (a) treating degree=0 isolates as
# threshold=0 and (b) ignoring degree-based information weights.
#
# Models (estimated for both GAD and MDD variants of mental health covariate):
#   (A) Yuchan baseline: glm(threshold ~ X, quasibinomial)
#                        — includes isolates (n=0) with threshold=0
#   (B) Quasibinomial without isolates: same model but n_alters > 0 only
#                        — isolates the effect of dropping isolates
#   (C) Grouped binomial: glm(cbind(k, n-k) ~ X, binomial), n>0 only
#                        — weights by degree, treats k/n as binomial proportion
#
# Outputs:
#   outputs_AO/model/grouped_binomial_coefficients-AO.csv
#   outputs_AO/model/grouped_binomial_forestplot-AO.pdf
#   outputs_AO/model/grouped_binomial_summaries-AO.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
in_rds       <- file.path(project_root, "outputs_AO", "intermediate", "thr_data_with_kn-AO.rds")
out_dir      <- file.path(project_root, "outputs_AO", "model")

thr <- readRDS(in_rds)
cat(sprintf("Loaded thr_data_with_kn-AO.rds: %d rows\n", nrow(thr)))
cat(sprintf("With degree>0: %d   isolates (degree=0): %d\n",
            sum(thr$n_alters > 0, na.rm = TRUE),
            sum(thr$n_alters == 0, na.rm = TRUE)))

# ---- Fit helpers -----------------------------------------------------------
covs <- c("factor(toa)","cohort","female","hispanic","asian","sex_min",
          "<mh>","friends_ecig")

fit_quasi <- function(data, mh) {
  req <- c("threshold","toa","cohort","female","hispanic","asian","sex_min",
           mh,"friends_ecig")
  d   <- data[complete.cases(data[, req]), ]
  fml <- as.formula(paste("threshold ~ factor(toa) + cohort + female + hispanic +",
                          "asian + sex_min +", mh, "+ friends_ecig"))
  list(fit = glm(fml, data = d, family = quasibinomial(link = "logit")),
       n   = nrow(d))
}

fit_grouped <- function(data, mh) {
  data <- data[data$n_alters > 0, ]
  req <- c("k_users","n_alters","toa","cohort","female","hispanic","asian",
           "sex_min", mh,"friends_ecig")
  d   <- data[complete.cases(data[, req]), ]
  fml <- as.formula(paste("cbind(k_users, n_alters - k_users) ~ factor(toa) + ",
                          "cohort + female + hispanic + asian + sex_min +",
                          mh, "+ friends_ecig"))
  list(fit = glm(fml, data = d, family = binomial(link = "logit")),
       n   = nrow(d))
}

# Models for each MH variant
specs <- list(
  list(model = "A_quasi_all",    mh = "gad",
       title = "(A) Quasibinomial — all (Yuchan baseline) — GAD"),
  list(model = "A_quasi_all",    mh = "mdd",
       title = "(A) Quasibinomial — all (Yuchan baseline) — MDD"),
  list(model = "B_quasi_noiso",  mh = "gad",
       title = "(B) Quasibinomial — non-isolates only — GAD"),
  list(model = "B_quasi_noiso",  mh = "mdd",
       title = "(B) Quasibinomial — non-isolates only — MDD"),
  list(model = "C_grouped",      mh = "gad",
       title = "(C) Grouped binomial cbind(k, n-k) — GAD"),
  list(model = "C_grouped",      mh = "mdd",
       title = "(C) Grouped binomial cbind(k, n-k) — MDD")
)

run_one <- function(s) {
  if (s$model == "A_quasi_all") {
    res <- fit_quasi(thr, s$mh)
  } else if (s$model == "B_quasi_noiso") {
    res <- fit_quasi(thr[thr$n_alters > 0, ], s$mh)
  } else if (s$model == "C_grouped") {
    res <- fit_grouped(thr, s$mh)
  }
  ct <- summary(res$fit)$coefficients
  parm <- rownames(ct)
  parm <- gsub("factor\\(toa\\)", "GP", parm)
  data.frame(model     = s$model,
             mh        = s$mh,
             title     = s$title,
             parameter = parm,
             estimate  = ct[, 1],
             se        = ct[, 2],
             zval      = ct[, 3],
             pval      = ct[, 4],
             or        = exp(ct[, 1]),
             ci_lo     = ct[, 1] - 1.96 * ct[, 2],
             ci_hi     = ct[, 1] + 1.96 * ct[, 2],
             n_obs     = res$n,
             stringsAsFactors = FALSE)
}

all_coefs <- bind_rows(lapply(specs, run_one)) |>
  mutate(across(c(estimate, se, zval, pval, or, ci_lo, ci_hi),
                ~ round(.x, 4)))

write.csv(all_coefs,
          file.path(out_dir, "grouped_binomial_coefficients-AO.csv"),
          row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "grouped_binomial_coefficients-AO.csv"), "\n")

# ---- Forest plot: compare estimates across models -------------------------
key_params <- c("hispanic","asian","female","sex_min","cohortC2",
                "friends_ecig","gad","mdd")
plot_df <- all_coefs |>
  filter(parameter %in% key_params) |>
  mutate(
    parameter = factor(parameter,
                       levels = c("hispanic","asian","female","sex_min",
                                  "cohortC2","friends_ecig","gad","mdd")),
    model_lbl = recode(model,
                       "A_quasi_all"   = "A. Quasi-bin (all, Yuchan)",
                       "B_quasi_noiso" = "B. Quasi-bin (no isolates)",
                       "C_grouped"     = "C. Grouped bin (cbind k,n-k)"),
    mh_lbl    = recode(mh, "gad" = "GAD model", "mdd" = "MDD model")
  )

p_forest <- ggplot(plot_df,
                   aes(x = estimate, y = parameter, color = model_lbl,
                       shape = model_lbl)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.55),
                  size = 0.6, linewidth = 0.9) +
  scale_color_manual(values = c("A. Quasi-bin (all, Yuchan)" = "#888888",
                                "B. Quasi-bin (no isolates)" = "#3b82c4",
                                "C. Grouped bin (cbind k,n-k)" = "#e76f51")) +
  scale_shape_manual(values = c("A. Quasi-bin (all, Yuchan)" = 16,
                                "B. Quasi-bin (no isolates)" = 17,
                                "C. Grouped bin (cbind k,n-k)" = 15)) +
  facet_wrap(~ mh_lbl, ncol = 2) +
  labs(title = "Threshold regression — coefficient comparison across model specs",
       subtitle = "Logit-scale estimates with 95% CI. Reference: Yuchan baseline (grey).",
       x = "Estimate (log-OR)", y = NULL,
       color = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "grey92"))

pdf(file.path(out_dir, "grouped_binomial_forestplot-AO.pdf"),
    width = 11, height = 6.5)
print(p_forest)
dev.off()
cat("Saved:", file.path(out_dir, "grouped_binomial_forestplot-AO.pdf"), "\n")

# ---- Full summaries to a text file ---------------------------------------
sink(file.path(out_dir, "grouped_binomial_summaries-AO.txt"))
for (s in specs) {
  cat(sprintf("\n\n=== %s ===\n", s$title))
  if (s$model == "A_quasi_all") {
    res <- fit_quasi(thr, s$mh)
  } else if (s$model == "B_quasi_noiso") {
    res <- fit_quasi(thr[thr$n_alters > 0, ], s$mh)
  } else if (s$model == "C_grouped") {
    res <- fit_grouped(thr, s$mh)
  }
  cat(sprintf("n used = %d\n", res$n))
  print(summary(res$fit))
}
sink()
cat("Saved:", file.path(out_dir, "grouped_binomial_summaries-AO.txt"), "\n")

# ---- Concise comparison table for the report ------------------------------
key_rows <- all_coefs |>
  filter(parameter %in% c("hispanic","friends_ecig","gad","mdd","cohortC2","asian","female","sex_min")) |>
  mutate(or_ci = sprintf("%.2f (%.2f, %.2f)", or, exp(ci_lo), exp(ci_hi)),
         pstar = ifelse(pval < .001, "***",
                  ifelse(pval < .01,  "**",
                  ifelse(pval < .05,  "*",
                  ifelse(pval < .10,  ".", ""))))) |>
  select(model, mh, parameter, or_ci, pval, pstar, n_obs)
print(key_rows, n = Inf)

cat("\nDone.\n")

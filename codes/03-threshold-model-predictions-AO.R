# 03-threshold-model-predictions-AO.R
# Refit Yuchan's quasibinomial threshold regression and compute
# model-predicted threshold distributions across each covariate.
# Side-by-side comparison: empirical (boxplot+violin) vs predicted (point + 95% CI).
#
# Model = exactly Yuchan's specification:
#   threshold ~ factor(toa) + cohort + female + hispanic + asian + sex_min + <mh> + friends_ecig
# where <mh> is gad (primary, per the report) or mdd.
#
# Two flavors of model-based predictions:
#   1. AME-style "average prediction by group" (averages over all other covariates
#      in the sample at each level of the focal covariate) via marginaleffects::avg_predictions().
#   2. Typical-individual prediction at modal/mean covariates via marginaleffects::predictions(datagrid()).
# We use (1) as the primary contrast vs the empirical distributions.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(marginaleffects)
  library(patchwork)
})

project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
in_rds       <- file.path(project_root, "outputs_AO", "intermediate", "thr_data-AO.rds")
out_dir      <- file.path(project_root, "outputs_AO", "model")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

d <- readRDS(in_rds)
thr_sub  <- d$thr_sub
gp_labels <- d$gp_labels
all_gps   <- d$all_gps

# Helpful labels mirroring _02
thr_sub <- thr_sub |>
  mutate(
    gp_label = factor(gp_labels[toa - min(all_gps) + 1L], levels = gp_labels),
    hispanic_f = factor(ifelse(hispanic == 1, "Hispanic", "Non-Hispanic"),
                        levels = c("Non-Hispanic", "Hispanic")),
    asian_f    = factor(ifelse(asian == 1, "Asian", "Non-Asian"),
                        levels = c("Non-Asian", "Asian")),
    female_f   = factor(ifelse(female == 1, "Female", "Male"),
                        levels = c("Male", "Female")),
    sex_min_f  = factor(ifelse(sex_min == 1, "Sexual minority", "Non-sex-minority"),
                        levels = c("Non-sex-minority", "Sexual minority")),
    cohort_f   = factor(cohort, levels = c("C1", "C2"))
  )

# ---- Fit Yuchan's threshold regression (both MH variants) -----------------
req_g <- c("threshold","toa","cohort","female","hispanic","asian","sex_min","gad","friends_ecig")
req_m <- c("threshold","toa","cohort","female","hispanic","asian","sex_min","mdd","friends_ecig")
dat_g <- thr_sub[complete.cases(thr_sub[, req_g]), ]
dat_m <- thr_sub[complete.cases(thr_sub[, req_m]), ]
cat(sprintf("GAD-model n: %d   MDD-model n: %d\n", nrow(dat_g), nrow(dat_m)))

fit_gad <- glm(threshold ~ factor(toa) + cohort + female + hispanic + asian + sex_min +
                          gad + friends_ecig,
               data = dat_g, family = quasibinomial(link = "logit"))
fit_mdd <- glm(threshold ~ factor(toa) + cohort + female + hispanic + asian + sex_min +
                          mdd + friends_ecig,
               data = dat_m, family = quasibinomial(link = "logit"))

# ---- Predicted thresholds by group (model = GAD primary) ------------------
# avg_predictions averages predicted threshold over all *other* covariates' empirical
# distributions, at each level of the focal covariate (i.e., g-computation / marginal mean).
get_pred <- function(model, by, transform_x = identity, focal_levels = NULL) {
  p <- avg_predictions(model, by = by)
  p <- as.data.frame(p)
  if (!is.null(focal_levels)) p[[by]] <- factor(p[[by]], levels = focal_levels)
  p$x_label <- transform_x(p[[by]])
  p
}

binary_label <- function(var, lab0, lab1) {
  function(x) factor(ifelse(x == 1, lab1, lab0), levels = c(lab0, lab1))
}

bin_specs <- list(
  list(focal = "hispanic", focal_levels = c(0, 1),
       relabel = binary_label("hispanic", "Non-Hispanic", "Hispanic"),
       title   = "Hispanic ethnicity"),
  list(focal = "asian",    focal_levels = c(0, 1),
       relabel = binary_label("asian",    "Non-Asian",    "Asian"),
       title   = "Asian race"),
  list(focal = "female",   focal_levels = c(0, 1),
       relabel = binary_label("female",   "Male",         "Female"),
       title   = "Gender (sex at birth, 0=Female fixed)"),
  list(focal = "sex_min",  focal_levels = c(0, 1),
       relabel = binary_label("sex_min",  "Non-sex-minority","Sexual minority"),
       title   = "Sexual minority status")
)

cat_specs <- list(
  list(focal = "cohort", focal_levels = c("C1","C2"),
       relabel = function(x) factor(x, levels = c("C1","C2")),
       title = "Cohort"),
  list(focal = "toa", focal_levels = sort(unique(dat_g$toa)),
       relabel = function(x) factor(gp_labels[as.integer(as.character(x)) - min(all_gps) + 1L],
                                    levels = gp_labels),
       title = "Grade period of adoption")
)

# Continuous: predict at a fine grid for friends_ecig, gad, mdd.
cont_grid_pred <- function(model, focal, grid_vals) {
  args <- setNames(list(grid_vals), focal)
  nd   <- do.call(datagrid, c(list(model = model), args))
  p    <- predictions(model, newdata = nd)
  as.data.frame(p)
}

# ---- Build the empirical-vs-predicted twin panel for one categorical/binary covariate
twin_panel <- function(emp_df, focal_factor_var, pred_df, pred_x_var, title, xlab) {
  emp_df <- emp_df[!is.na(emp_df[[focal_factor_var]]) & !is.na(emp_df$threshold), ]
  p_emp <- ggplot(emp_df,
                  aes(x = .data[[focal_factor_var]], y = threshold,
                      fill = .data[[focal_factor_var]])) +
    geom_violin(alpha = 0.55, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.07, height = 0, alpha = 0.18, size = 0.4) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Empirical distribution",
         x = xlab, y = "Threshold") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))

  p_mod <- ggplot(pred_df,
                  aes(x = .data[[pred_x_var]], y = estimate)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high,
                        color = .data[[pred_x_var]]), size = 0.9, linewidth = 0.9) +
    scale_y_continuous(limits = c(0, max(0.3, max(pred_df$conf.high, na.rm = TRUE) * 1.1))) +
    labs(title = "Model-predicted threshold (95% CI)\n(avg over other covariates)",
         x = xlab, y = "Predicted threshold") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))

  (p_emp | p_mod) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
}

# Continuous twin panel: scatter (empirical) + line (predicted with CI ribbon)
twin_panel_cont <- function(emp_df, focal_x_var, pred_df, focal_x_var_pred,
                             title, xlab) {
  emp_df <- emp_df[!is.na(emp_df[[focal_x_var]]) & !is.na(emp_df$threshold), ]
  p_emp <- ggplot(emp_df, aes(x = .data[[focal_x_var]], y = threshold)) +
    geom_jitter(width = 0.04, height = 0, alpha = 0.25, size = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "steelblue",
                alpha = 0.18, formula = y ~ x) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Empirical (jitter + LOESS)", x = xlab, y = "Threshold") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))

  p_mod <- ggplot(pred_df, aes(x = .data[[focal_x_var_pred]], y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "tomato", alpha = 0.20) +
    geom_line(color = "tomato", linewidth = 1.1) +
    scale_y_continuous(limits = c(0, max(0.5, max(pred_df$conf.high, na.rm = TRUE) * 1.1))) +
    labs(title = "Model-predicted (95% CI)", x = xlab, y = "Predicted threshold") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))

  (p_emp | p_mod) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
}

# ---- Generate plots for GAD model (primary) -------------------------------
pdf(file.path(out_dir, "threshold_model_predictions_GAD-AO.pdf"),
    width = 11, height = 5.5)

# Binary covariates
for (s in bin_specs) {
  pred <- get_pred(fit_gad, by = s$focal, transform_x = s$relabel)
  pred$x_label <- s$relabel(pred[[s$focal]])
  print(twin_panel(emp_df = thr_sub, focal_factor_var = paste0(s$focal, "_f"),
                   pred_df = pred, pred_x_var = "x_label",
                   title = paste0("Threshold by ", s$title, " (GAD model)"),
                   xlab = s$title))
}

# Categorical covariates
for (s in cat_specs) {
  pred <- get_pred(fit_gad, by = s$focal, transform_x = s$relabel)
  pred$x_label <- s$relabel(pred[[s$focal]])
  emp_var <- if (s$focal == "toa") "gp_label" else paste0(s$focal, "_f")
  print(twin_panel(emp_df = thr_sub, focal_factor_var = emp_var,
                   pred_df = pred, pred_x_var = "x_label",
                   title = paste0("Threshold by ", s$title, " (GAD model)"),
                   xlab = s$title))
}

# Continuous: friends_ecig (1-5)
fe_vals  <- sort(unique(dat_g$friends_ecig))
pred_fe  <- cont_grid_pred(fit_gad, "friends_ecig", fe_vals)
print(twin_panel_cont(thr_sub, "friends_ecig", pred_fe, "friends_ecig",
                      "Threshold by perceived friend e-cig use (GAD model)",
                      "Friends using e-cig (self-report, 1-5)"))

# Continuous: GAD
gad_grid <- seq(0, max(dat_g$gad, na.rm = TRUE), length.out = 25)
pred_gad <- cont_grid_pred(fit_gad, "gad", gad_grid)
print(twin_panel_cont(thr_sub, "gad", pred_gad, "gad",
                      "Threshold by GAD (anxiety) at TOA",
                      "GAD (RCADS anxiety subscale mean)"))

dev.off()
cat("Saved:", file.path(out_dir, "threshold_model_predictions_GAD-AO.pdf"), "\n")

# ---- Generate plots for MDD model (parallel) ------------------------------
pdf(file.path(out_dir, "threshold_model_predictions_MDD-AO.pdf"),
    width = 11, height = 5.5)

for (s in bin_specs) {
  pred <- get_pred(fit_mdd, by = s$focal, transform_x = s$relabel)
  pred$x_label <- s$relabel(pred[[s$focal]])
  print(twin_panel(emp_df = thr_sub, focal_factor_var = paste0(s$focal, "_f"),
                   pred_df = pred, pred_x_var = "x_label",
                   title = paste0("Threshold by ", s$title, " (MDD model)"),
                   xlab = s$title))
}

for (s in cat_specs) {
  pred <- get_pred(fit_mdd, by = s$focal, transform_x = s$relabel)
  pred$x_label <- s$relabel(pred[[s$focal]])
  emp_var <- if (s$focal == "toa") "gp_label" else paste0(s$focal, "_f")
  print(twin_panel(emp_df = thr_sub, focal_factor_var = emp_var,
                   pred_df = pred, pred_x_var = "x_label",
                   title = paste0("Threshold by ", s$title, " (MDD model)"),
                   xlab = s$title))
}

fe_vals_m  <- sort(unique(dat_m$friends_ecig))
pred_fe_m  <- cont_grid_pred(fit_mdd, "friends_ecig", fe_vals_m)
print(twin_panel_cont(thr_sub, "friends_ecig", pred_fe_m, "friends_ecig",
                      "Threshold by perceived friend e-cig use (MDD model)",
                      "Friends using e-cig (self-report, 1-5)"))

mdd_grid <- seq(0, max(dat_m$mdd, na.rm = TRUE), length.out = 25)
pred_mdd <- cont_grid_pred(fit_mdd, "mdd", mdd_grid)
print(twin_panel_cont(thr_sub, "mdd", pred_mdd, "mdd",
                      "Threshold by MDD (depression) at TOA",
                      "MDD (RCADS depression subscale mean)"))

dev.off()
cat("Saved:", file.path(out_dir, "threshold_model_predictions_MDD-AO.pdf"), "\n")

# ---- Save tidy predictions CSV -------------------------------------------
collect_preds <- function(model, model_label) {
  rows <- list()
  for (s in c(bin_specs, cat_specs)) {
    p <- avg_predictions(model, by = s$focal)
    p <- as.data.frame(p)
    p$model    <- model_label
    p$covariate <- s$focal
    rows[[length(rows) + 1]] <- p
  }
  bind_rows(rows)
}
pred_table <- bind_rows(collect_preds(fit_gad, "GAD"),
                        collect_preds(fit_mdd, "MDD"))
write.csv(pred_table |>
            mutate(across(where(is.numeric), ~ round(.x, 4))),
          file.path(out_dir, "threshold_model_predictions-AO.csv"),
          row.names = FALSE)
cat("Saved:", file.path(out_dir, "threshold_model_predictions-AO.csv"), "\n")

# ---- Save model summaries -------------------------------------------------
sink(file.path(out_dir, "threshold_model_summaries-AO.txt"))
cat("===== Threshold regression — GAD model =====\n")
print(summary(fit_gad))
cat("\n===== Threshold regression — MDD model =====\n")
print(summary(fit_mdd))
sink()
cat("Saved:", file.path(out_dir, "threshold_model_summaries-AO.txt"), "\n")
cat("Done.\n")

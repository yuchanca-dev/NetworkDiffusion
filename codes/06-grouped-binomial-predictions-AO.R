# 06-grouped-binomial-predictions-AO.R
# Marginal threshold predictions using Model C (grouped binomial,
# cbind(k, n-k) ~ X, non-isolates only). Mirrors script 03 but with
# the corrected specification.
#
# Outputs:
#   outputs_AO/model/grouped_binomial_predictions_GAD-AO.pdf
#   outputs_AO/model/grouped_binomial_predictions_MDD-AO.pdf
#   outputs_AO/model/grouped_binomial_predictions-AO.csv

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(marginaleffects)
  library(patchwork)
})

project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
in_rds       <- file.path(project_root, "outputs_AO", "intermediate", "thr_data_with_kn-AO.rds")
out_dir      <- file.path(project_root, "outputs_AO", "model")

thr <- readRDS(in_rds)
gp_labels <- c("10Fa","10Sp","11Fa","11Sp","12Fa","12Sp")
all_gps   <- 3:8

# Add factor labels (mirrors script 03)
thr <- thr |>
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

# Fit Model C: grouped binomial, non-isolates only
fit_grouped <- function(data, mh) {
  data <- data[data$n_alters > 0, ]
  req <- c("k_users","n_alters","toa","cohort","female","hispanic","asian",
           "sex_min", mh, "friends_ecig")
  d   <- data[complete.cases(data[, req]), ]
  fml <- as.formula(paste("cbind(k_users, n_alters - k_users) ~ factor(toa) + ",
                          "cohort + female + hispanic + asian + sex_min +",
                          mh, "+ friends_ecig"))
  list(fit = glm(fml, data = d, family = binomial(link = "logit")),
       data = d)
}

g <- fit_grouped(thr, "gad")
m <- fit_grouped(thr, "mdd")
cat(sprintf("GAD model C: n=%d   MDD model C: n=%d\n",
            nrow(g$data), nrow(m$data)))

# ---- Twin-panel plotting (mirrors script 03) ------------------------------
twin_panel <- function(emp_df, focal_factor_var, pred_df, pred_x_var, title, xlab) {
  emp_df <- emp_df[!is.na(emp_df[[focal_factor_var]]) & !is.na(emp_df$threshold), ]
  p_emp <- ggplot(emp_df,
                  aes(x = .data[[focal_factor_var]], y = threshold,
                      fill = .data[[focal_factor_var]])) +
    geom_violin(alpha = 0.55, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.07, height = 0, alpha = 0.18, size = 0.4) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Empirical (non-isolates only)", x = xlab, y = "Threshold") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))
  p_mod <- ggplot(pred_df,
                  aes(x = .data[[pred_x_var]], y = estimate)) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high,
                        color = .data[[pred_x_var]]), size = 0.9, linewidth = 0.9) +
    scale_y_continuous(limits = c(0, max(0.3, max(pred_df$conf.high, na.rm = TRUE) * 1.1))) +
    labs(title = "Grouped binomial — predicted prob alter is using (95% CI)",
         x = xlab, y = "Predicted P(alter using)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"))
  (p_emp | p_mod) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
}

twin_panel_cont <- function(emp_df, focal_x_var, pred_df, focal_x_var_pred, title, xlab) {
  emp_df <- emp_df[!is.na(emp_df[[focal_x_var]]) & !is.na(emp_df$threshold), ]
  emp_df <- emp_df[emp_df$n_alters > 0, ]
  p_emp <- ggplot(emp_df, aes(x = .data[[focal_x_var]], y = threshold)) +
    geom_jitter(width = 0.04, height = 0, alpha = 0.25, size = 0.6) +
    geom_smooth(method = "loess", se = TRUE, color = "steelblue", fill = "steelblue",
                alpha = 0.18, formula = y ~ x) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Empirical (non-isolates) + LOESS", x = xlab, y = "Threshold") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))
  p_mod <- ggplot(pred_df, aes(x = .data[[focal_x_var_pred]], y = estimate)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "tomato", alpha = 0.20) +
    geom_line(color = "tomato", linewidth = 1.1) +
    scale_y_continuous(limits = c(0, max(0.5, max(pred_df$conf.high, na.rm = TRUE) * 1.1))) +
    labs(title = "Grouped binomial — predicted (95% CI)", x = xlab,
         y = "Predicted P(alter using)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold"))
  (p_emp | p_mod) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
}

binary_label <- function(var, lab0, lab1) {
  function(x) factor(ifelse(x == 1, lab1, lab0), levels = c(lab0, lab1))
}

bin_specs <- list(
  list(focal = "hispanic", relabel = binary_label("hispanic", "Non-Hispanic", "Hispanic"),
       title = "Hispanic ethnicity"),
  list(focal = "asian",    relabel = binary_label("asian", "Non-Asian", "Asian"),
       title = "Asian race"),
  list(focal = "female",   relabel = binary_label("female", "Male", "Female"),
       title = "Gender (sex at birth)"),
  list(focal = "sex_min",  relabel = binary_label("sex_min", "Non-sex-minority", "Sexual minority"),
       title = "Sexual minority status")
)
cat_specs <- list(
  list(focal = "cohort",
       relabel = function(x) factor(x, levels = c("C1","C2")),
       title = "Cohort"),
  list(focal = "toa",
       relabel = function(x) factor(gp_labels[as.integer(as.character(x)) - min(all_gps) + 1L],
                                    levels = gp_labels),
       title = "Grade period of adoption")
)
cont_grid_pred <- function(model, focal, grid_vals) {
  args <- setNames(list(grid_vals), focal)
  nd   <- do.call(datagrid, c(list(model = model), args))
  p    <- predictions(model, newdata = nd)
  as.data.frame(p)
}
get_pred <- function(model, by, relabel) {
  p <- as.data.frame(avg_predictions(model, by = by))
  p$x_label <- relabel(p[[by]])
  p
}

render_pdf <- function(fit, dat, mh_name, out_path) {
  pdf(out_path, width = 11, height = 5.5)
  for (s in bin_specs) {
    pred <- get_pred(fit, by = s$focal, relabel = s$relabel)
    print(twin_panel(thr[thr$n_alters > 0, ], paste0(s$focal, "_f"),
                     pred, "x_label",
                     paste0("Threshold by ", s$title,
                            " (Grouped binomial, ", toupper(mh_name), ")"),
                     s$title))
  }
  for (s in cat_specs) {
    pred <- get_pred(fit, by = s$focal, relabel = s$relabel)
    emp_var <- if (s$focal == "toa") "gp_label" else paste0(s$focal, "_f")
    print(twin_panel(thr[thr$n_alters > 0, ], emp_var,
                     pred, "x_label",
                     paste0("Threshold by ", s$title,
                            " (Grouped binomial, ", toupper(mh_name), ")"),
                     s$title))
  }
  # Continuous: friends_ecig
  fe_vals <- sort(unique(dat$friends_ecig))
  pred_fe <- cont_grid_pred(fit, "friends_ecig", fe_vals)
  print(twin_panel_cont(thr, "friends_ecig", pred_fe, "friends_ecig",
                        paste0("Threshold by perceived friend use (Grouped binomial, ",
                               toupper(mh_name), ")"),
                        "Friends using e-cig (1-5)"))
  # Continuous: mental health
  mh_grid <- seq(0, max(dat[[mh_name]], na.rm = TRUE), length.out = 25)
  pred_mh <- cont_grid_pred(fit, mh_name, mh_grid)
  print(twin_panel_cont(thr, mh_name, pred_mh, mh_name,
                        paste0("Threshold by ", toupper(mh_name),
                               " at TOA (Grouped binomial)"),
                        paste0(toupper(mh_name), " (RCADS subscale mean)")))
  dev.off()
}

render_pdf(g$fit, g$data, "gad",
           file.path(out_dir, "grouped_binomial_predictions_GAD-AO.pdf"))
cat("Saved:", file.path(out_dir, "grouped_binomial_predictions_GAD-AO.pdf"), "\n")
render_pdf(m$fit, m$data, "mdd",
           file.path(out_dir, "grouped_binomial_predictions_MDD-AO.pdf"))
cat("Saved:", file.path(out_dir, "grouped_binomial_predictions_MDD-AO.pdf"), "\n")

# ---- Tidy predictions CSV -------------------------------------------------
collect_preds <- function(model, model_label) {
  rows <- list()
  for (s in c(bin_specs, cat_specs)) {
    p <- as.data.frame(avg_predictions(model, by = s$focal))
    p$model <- model_label
    p$covariate <- s$focal
    rows[[length(rows) + 1]] <- p
  }
  bind_rows(rows)
}
pred_tbl <- bind_rows(collect_preds(g$fit, "GAD"),
                      collect_preds(m$fit, "MDD")) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))
write.csv(pred_tbl,
          file.path(out_dir, "grouped_binomial_predictions-AO.csv"),
          row.names = FALSE)
cat("Saved:", file.path(out_dir, "grouped_binomial_predictions-AO.csv"), "\n")
cat("Done.\n")

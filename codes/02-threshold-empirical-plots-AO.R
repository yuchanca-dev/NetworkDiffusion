# 02-threshold-empirical-plots-AO.R
# Empirical adoption-threshold distributions stratified by covariates.
# Reads outputs_AO/intermediate/thr_data-AO.rds (produced by 01-rebuild-...).
# Writes:
#   outputs_AO/empirical/threshold_empirical_by_covariate-AO.pdf  (multi-page)
#   outputs_AO/empirical/threshold_empirical_summary-AO.csv       (per-group N, median, IQR, test p)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
})

project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
in_rds       <- file.path(project_root, "outputs_AO", "intermediate", "thr_data-AO.rds")
out_dir      <- file.path(project_root, "outputs_AO", "empirical")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

d <- readRDS(in_rds)
thr_sub  <- d$thr_sub
gp_labels <- d$gp_labels
all_gps   <- d$all_gps

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
    cohort_f   = factor(cohort, levels = c("C1", "C2")),
    sch_type_f = factor(sch_type,
                        levels = c("Asian-majority", "Hispanic-majority", "Other")),
    school_f   = factor(school)
  )

# Helper: tertile bins for continuous covariates
tertile_bin <- function(x, label_prefix) {
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  cuts <- cut(x, breaks = q, include.lowest = TRUE, labels = c(
    sprintf("%s low (%.2f-%.2f)",  label_prefix, q[1], q[2]),
    sprintf("%s mid (%.2f-%.2f)",  label_prefix, q[2], q[3]),
    sprintf("%s high (%.2f-%.2f)", label_prefix, q[3], q[4])
  ))
  cuts
}
thr_sub$mdd_bin <- tertile_bin(thr_sub$mdd, "MDD")
thr_sub$gad_bin <- tertile_bin(thr_sub$gad, "GAD")

# friends_ecig is already a 1-5 Likert; treat as ordered factor
thr_sub$friends_ecig_f <- factor(thr_sub$friends_ecig, ordered = TRUE,
                                 levels = sort(unique(thr_sub$friends_ecig[!is.na(thr_sub$friends_ecig)])))

# ---- Plot helper -----------------------------------------------------------
plot_threshold_by_group <- function(df, group_var, title, xlab,
                                    test = c("auto", "wilcox", "kruskal", "none")) {
  test <- match.arg(test)
  df  <- df[!is.na(df[[group_var]]) & !is.na(df$threshold), ]
  if (nrow(df) == 0) return(NULL)

  ng <- length(unique(df[[group_var]]))
  if (test == "auto") test <- if (ng == 2) "wilcox" else if (ng > 2) "kruskal" else "none"

  n_per <- df |>
    group_by(.data[[group_var]]) |>
    summarise(n = n(), med = median(threshold), .groups = "drop")
  x_labels <- sprintf("%s\nn=%d\nmed=%.2f",
                      as.character(n_per[[group_var]]), n_per$n, n_per$med)

  p <- ggplot(df, aes(x = .data[[group_var]], y = threshold, fill = .data[[group_var]])) +
    geom_violin(alpha = 0.55, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.07, height = 0, alpha = 0.25, size = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_discrete(labels = x_labels) +
    labs(title = title, x = xlab, y = "Adoption threshold\n(fraction of friends already using)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12),
          axis.text.x = element_text(size = 9))

  if (test == "wilcox") {
    p <- p + stat_compare_means(method = "wilcox.test",
                                label = "p.format", label.y = 0.97, size = 3.4)
  } else if (test == "kruskal") {
    p <- p + stat_compare_means(method = "kruskal.test",
                                label = "p.format", label.y = 0.97, size = 3.4)
  }
  p
}

# ---- Spec covariates to plot ----------------------------------------------
specs <- list(
  list(var = "hispanic_f",   title = "Threshold by Hispanic ethnicity",          xlab = "Hispanic"),
  list(var = "asian_f",      title = "Threshold by Asian race",                  xlab = "Asian"),
  list(var = "female_f",     title = "Threshold by gender",                      xlab = "Gender"),
  list(var = "sex_min_f",    title = "Threshold by sexual minority status",      xlab = "Sexual minority"),
  list(var = "cohort_f",     title = "Threshold by cohort",                      xlab = "Cohort"),
  list(var = "sch_type_f",   title = "Threshold by school type",                 xlab = "School type"),
  list(var = "gp_label",     title = "Threshold by grade period of adoption",    xlab = "Grade period (TOA)"),
  list(var = "school_f",     title = "Threshold by school (14 schools)",         xlab = "School"),
  list(var = "mdd_bin",      title = "Threshold by MDD (tertiles at TOA)",       xlab = "MDD tertile"),
  list(var = "gad_bin",      title = "Threshold by GAD (tertiles at TOA)",       xlab = "GAD tertile"),
  list(var = "friends_ecig_f", title = "Threshold by perceived friend e-cig use", xlab = "Friends using e-cig (self-report, 1-5)")
)

plots <- lapply(specs, function(s) {
  plot_threshold_by_group(thr_sub, s$var, s$title, s$xlab)
})

# Save multi-page PDF
pdf(file.path(out_dir, "threshold_empirical_by_covariate-AO.pdf"),
    width = 8.5, height = 5.5)
for (p in plots) if (!is.null(p)) print(p)
dev.off()
cat("Saved:", file.path(out_dir, "threshold_empirical_by_covariate-AO.pdf"), "\n")

# ---- Summary CSV: N, median, IQR, test p per covariate × group ------------
summarise_by <- function(df, group_var) {
  df  <- df[!is.na(df[[group_var]]) & !is.na(df$threshold), ]
  if (nrow(df) == 0) return(NULL)
  s <- df |>
    group_by(.data[[group_var]]) |>
    summarise(n = n(),
              median  = median(threshold),
              q25     = quantile(threshold, 0.25),
              q75     = quantile(threshold, 0.75),
              mean    = mean(threshold),
              sd      = sd(threshold),
              .groups = "drop") |>
    rename(group = 1) |>
    mutate(covariate = group_var,
           group     = as.character(group)) |>
    select(covariate, group, n, median, q25, q75, mean, sd)
  # add test p
  ng <- length(unique(df[[group_var]]))
  if (ng == 2) {
    s$test     <- "wilcox"
    s$test_p   <- tryCatch(
      wilcox.test(threshold ~ df[[group_var]], data = df)$p.value,
      error = function(e) NA_real_)
  } else if (ng > 2) {
    s$test   <- "kruskal"
    s$test_p <- tryCatch(
      kruskal.test(threshold ~ df[[group_var]], data = df)$p.value,
      error = function(e) NA_real_)
  } else {
    s$test <- NA_character_; s$test_p <- NA_real_
  }
  s
}

summary_tbl <- bind_rows(lapply(specs, function(s) summarise_by(thr_sub, s$var)))
summary_tbl <- summary_tbl |>
  mutate(across(c(median, q25, q75, mean, sd, test_p),
                ~ round(.x, 4)))

write.csv(summary_tbl,
          file.path(out_dir, "threshold_empirical_summary-AO.csv"),
          row.names = FALSE)
cat("Saved:", file.path(out_dir, "threshold_empirical_summary-AO.csv"), "\n")

# ---- Bonus: faceted views (toa × key binary covariates) -------------------
faceted_plot <- function(df, group_var, title) {
  df <- df[!is.na(df[[group_var]]) & !is.na(df$threshold), ]
  ggplot(df, aes(x = .data[[group_var]], y = threshold, fill = .data[[group_var]])) +
    geom_violin(alpha = 0.55, color = NA, scale = "width") +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.shape = NA, color = "black") +
    facet_wrap(~ gp_label, nrow = 1) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = title, x = NULL,
         y = "Adoption threshold") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(face = "bold"))
}

pdf(file.path(out_dir, "threshold_empirical_facetted_by_toa-AO.pdf"),
    width = 11, height = 5)
for (s in list(
  list(var = "hispanic_f", title = "Threshold by Hispanic × grade period of adoption"),
  list(var = "asian_f",    title = "Threshold by Asian × grade period of adoption"),
  list(var = "sex_min_f",  title = "Threshold by Sex-min × grade period of adoption"),
  list(var = "cohort_f",   title = "Threshold by Cohort × grade period of adoption"),
  list(var = "sch_type_f", title = "Threshold by School-type × grade period of adoption")
)) {
  print(faceted_plot(thr_sub, s$var, s$title))
}
dev.off()
cat("Saved:", file.path(out_dir, "threshold_empirical_facetted_by_toa-AO.pdf"), "\n")

cat("Done.\n")

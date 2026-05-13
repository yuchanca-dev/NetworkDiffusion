# 04-degree-diagnostics-AO.R
# Diagnostic plots & tables to assess whether the heavy mass at low-fraction
# threshold values in Yuchan's distribution is driven by low ego out-degree.
# Inputs:  outputs_AO/intermediate/thr_data-AO.rds
# Outputs: outputs_AO/diagnostics/degree_diagnostics-AO.pdf
#          outputs_AO/diagnostics/degree_summary-AO.csv

suppressPackageStartupMessages({
  library(netdiffuseR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
in_rds       <- file.path(project_root, "outputs_AO", "intermediate", "thr_data-AO.rds")
out_dir      <- file.path(project_root, "outputs_AO", "diagnostics")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

d        <- readRDS(in_rds)
thr_sub  <- d$thr_sub
dn       <- d$diffnet_all
gp_labels <- d$gp_labels
all_gps   <- d$all_gps

# ---- Recover degree-at-TOA and # alters using at TOA ----------------------
dgr_mat   <- dgr(dn, cmode = "outdegree")         # [n × T]
cumadopt  <- dn$cumadopt                          # [n × T]
toa_vec   <- dn$toa
ego_ids   <- names(toa_vec)

# Column indices for dgr/cumadopt by GP
gp_col_idx <- function(gp) match(as.character(gp), colnames(dgr_mat))

# numerator k = # alters who had adopted by TOA[i]
# = (graph_t %*% cumadopt[, t])[i]  computed per time period
n_alters_at_toa <- integer(length(toa_vec))
k_users_at_toa  <- integer(length(toa_vec))
for (gp in all_gps) {
  idx <- which(toa_vec == gp)
  if (length(idx) == 0) next
  g   <- as.matrix(dn$graph[[as.character(gp)]])
  c_t <- cumadopt[, gp_col_idx(gp)]
  k_t <- as.integer(g %*% c_t)              # # alters using at gp
  n_t <- as.integer(rowSums(g))             # out-degree at gp
  n_alters_at_toa[idx] <- n_t[idx]
  k_users_at_toa[idx]  <- k_t[idx]
}
names(n_alters_at_toa) <- ego_ids
names(k_users_at_toa)  <- ego_ids

# Attach to thr_sub via the diffnet id (column `id` in thr_sub = vertex id)
adopt_diag <- data.frame(
  id        = ego_ids,
  n_alters  = n_alters_at_toa,
  k_users   = k_users_at_toa,
  stringsAsFactors = FALSE
)
adopt_diag <- adopt_diag[!is.na(toa_vec), ]

thr_diag <- merge(thr_sub, adopt_diag, by = "id", all.x = TRUE) |>
  mutate(
    threshold_implied = ifelse(n_alters > 0, k_users / n_alters, NA_real_),
    degree_bin = cut(n_alters, breaks = c(-Inf, 1, 3, 5, Inf),
                     labels = c("1 alter", "2-3 alters", "4-5 alters", "6+ alters"))
  )

# Sanity: compare implied threshold against netdiffuseR's threshold()
agree <- with(thr_diag,
              mean(abs(threshold - threshold_implied) < 1e-9, na.rm = TRUE))
cat(sprintf("Sanity check: %.1f%% of egos have threshold == k/n exactly.\n", 100 * agree))

# ---- Summary table --------------------------------------------------------
summary_overall <- thr_diag |>
  summarise(
    n_adopters   = n(),
    median_dgr   = median(n_alters, na.rm = TRUE),
    iqr_low      = quantile(n_alters, 0.25, na.rm = TRUE),
    iqr_high     = quantile(n_alters, 0.75, na.rm = TRUE),
    mean_dgr     = mean(n_alters, na.rm = TRUE),
    pct_dgr_le_1 = 100 * mean(n_alters <= 1, na.rm = TRUE),
    pct_dgr_le_2 = 100 * mean(n_alters <= 2, na.rm = TRUE),
    pct_dgr_le_3 = 100 * mean(n_alters <= 3, na.rm = TRUE),
    pct_dgr_eq_0 = 100 * mean(n_alters == 0, na.rm = TRUE)
  ) |>
  mutate(scope = "OVERALL")

summary_by_school_type <- thr_diag |>
  group_by(sch_type) |>
  summarise(
    n_adopters   = n(),
    median_dgr   = median(n_alters, na.rm = TRUE),
    iqr_low      = quantile(n_alters, 0.25, na.rm = TRUE),
    iqr_high     = quantile(n_alters, 0.75, na.rm = TRUE),
    mean_dgr     = mean(n_alters, na.rm = TRUE),
    pct_dgr_le_1 = 100 * mean(n_alters <= 1, na.rm = TRUE),
    pct_dgr_le_2 = 100 * mean(n_alters <= 2, na.rm = TRUE),
    pct_dgr_le_3 = 100 * mean(n_alters <= 3, na.rm = TRUE),
    pct_dgr_eq_0 = 100 * mean(n_alters == 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  rename(scope = sch_type)

summary_tbl <- bind_rows(summary_overall, summary_by_school_type) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) |>
  select(scope, everything())
print(summary_tbl)

write.csv(summary_tbl, file.path(out_dir, "degree_summary-AO.csv"), row.names = FALSE)

# ---- Plots ----------------------------------------------------------------
pdf(file.path(out_dir, "degree_diagnostics-AO.pdf"), width = 10, height = 6.5)

# (1) Histogram of degree-at-TOA (overall)
p1 <- ggplot(thr_diag, aes(x = n_alters)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.85) +
  geom_vline(xintercept = median(thr_diag$n_alters, na.rm = TRUE),
             linetype = "dashed", color = "tomato", linewidth = 0.9) +
  scale_x_continuous(breaks = 0:max(thr_diag$n_alters, na.rm = TRUE)) +
  labs(title = "Distribution of ego out-degree at time of adoption (TOA)",
       subtitle = sprintf("n_adopters = %d   median = %d   mean = %.2f",
                          nrow(thr_diag),
                          median(thr_diag$n_alters, na.rm = TRUE),
                          mean(thr_diag$n_alters, na.rm = TRUE)),
       x = "Number of alters at TOA (out-degree)",
       y = "Adopters") +
  theme_minimal(base_size = 12)
print(p1)

# (2) Histogram faceted by school type
p2 <- ggplot(thr_diag, aes(x = n_alters, fill = sch_type)) +
  geom_histogram(binwidth = 1, color = "white", alpha = 0.85) +
  facet_wrap(~ sch_type, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = 0:max(thr_diag$n_alters, na.rm = TRUE)) +
  labs(title = "Out-degree at TOA — by school type",
       x = "Number of alters at TOA", y = "Adopters") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
print(p2)

# (3) Scatter threshold vs degree, with rational guidelines k/n overlaid
# Build the rational lines for n in 1..max
max_n <- max(thr_diag$n_alters, na.rm = TRUE)
ratlines <- expand.grid(n = 1:max_n, k = 0:max_n) |>
  filter(k <= n) |>
  mutate(y = k / n)
p3 <- ggplot(thr_diag |> filter(n_alters > 0),
             aes(x = n_alters, y = threshold)) +
  geom_segment(data = ratlines,
               aes(x = n - 0.4, xend = n + 0.4, y = y, yend = y),
               color = "grey75", linewidth = 0.3, inherit.aes = FALSE) +
  geom_jitter(aes(color = hispanic == 1), width = 0.15, height = 0,
              alpha = 0.45, size = 0.9) +
  scale_color_manual(values = c("FALSE" = "#3b82c4", "TRUE" = "#e76f51"),
                     labels = c("Non-Hispanic", "Hispanic"), name = NULL) +
  scale_x_continuous(breaks = 1:max_n) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(title = "Threshold = k/n is constrained to rational fractions",
       subtitle = "Grey ticks mark allowed values; coloured by Hispanic",
       x = "Out-degree n at TOA",
       y = "Threshold (k users / n alters)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
print(p3)

# (4) Threshold distribution by degree bin
p4 <- ggplot(thr_diag |> filter(!is.na(degree_bin)),
             aes(x = degree_bin, y = threshold, fill = degree_bin)) +
  geom_violin(alpha = 0.55, color = NA, scale = "width") +
  geom_boxplot(width = 0.18, alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.07, height = 0, alpha = 0.25, size = 0.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(title = "Threshold distribution by degree bin at TOA",
       subtitle = "Low-degree bins show sharp modes (rational-fraction artefact)",
       x = NULL, y = "Threshold") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
print(p4)

# (5) Counts per degree value with median threshold annotation
counts_by_n <- thr_diag |>
  filter(!is.na(n_alters)) |>
  group_by(n_alters) |>
  summarise(n = n(),
            median_thr = median(threshold, na.rm = TRUE),
            mean_thr   = mean(threshold, na.rm = TRUE),
            .groups = "drop")
p5 <- ggplot(counts_by_n, aes(x = n_alters, y = n)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  geom_text(aes(label = sprintf("n=%d\nmed=%.2f", n, median_thr)),
            vjust = -0.3, size = 3) +
  scale_x_continuous(breaks = 0:max(counts_by_n$n_alters)) +
  labs(title = "Adopters and median threshold, by out-degree at TOA",
       x = "Out-degree n at TOA", y = "Adopters") +
  theme_minimal(base_size = 12) +
  coord_cartesian(ylim = c(0, max(counts_by_n$n) * 1.15))
print(p5)

# (6) Empirical threshold density vs predicted-from-model density (overlay)
# Just for the overall sample, no covariate adjustment — to show the
# multimodality more starkly.
p6 <- ggplot(thr_diag |> filter(!is.na(threshold)),
             aes(x = threshold)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "white", alpha = 0.85) +
  scale_x_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.2)) +
  labs(title = "Threshold distribution (all adopters, no covariates)",
       subtitle = "Spikes at 0, 1/4, 1/3, 1/2, 2/3, 1 = denominator artefacts",
       x = "Threshold (k/n)", y = "Adopters") +
  theme_minimal(base_size = 12)
print(p6)

dev.off()
cat("\nSaved:", file.path(out_dir, "degree_diagnostics-AO.pdf"), "\n")
cat("Saved:", file.path(out_dir, "degree_summary-AO.csv"), "\n")

# ---- Save k and n back to a richer intermediate ---------------------------
saveRDS(thr_diag,
        file.path(project_root, "outputs_AO", "intermediate", "thr_data_with_kn-AO.rds"))
cat("Saved: outputs_AO/intermediate/thr_data_with_kn-AO.rds (thr_sub + n_alters + k_users)\n")
cat("Done.\n")

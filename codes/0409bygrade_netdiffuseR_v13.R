rm(list = ls())
library(netdiffuseR)
library(openxlsx)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

# ---------------------------------------------------------------
# School definitions
# ---------------------------------------------------------------
c1_early_schools <- c(101, 102, 103, 104, 105)
c1_late_schools  <- c(106, 107, 112, 113, 114)
c2_schools       <- c(201, 212, 213, 214)
exclude_schools  <- c(108)

all_schools      <- setdiff(c(c1_early_schools, c1_late_schools, c2_schools), exclude_schools)
asian_schools    <- c(103, 105, 112, 113, 212, 213)
hispanic_schools <- c(102, 106, 107, 114, 214)

# ---------------------------------------------------------------
# Grade-wave mapping (per Dr. Valente, Apr 7 2026)
# Analysis window: Grade 10 Fall to Grade 12 Spring
#   C1 (101-114): W3-W8  (W3=GP3=10Fa ... W8=GP8=12Sp)
#   C2 (201-214): W5-W10 (W5=GP3=10Fa ... W10=GP8=12Sp)
# Set INCLUDE_GRADE9 <- TRUE to extend to Grade 9 Fall
# ---------------------------------------------------------------
INCLUDE_GRADE9 <- FALSE
all_gps   <- 3:8
gp_labels <- c("10Fa","10Sp","11Fa","11Sp","12Fa","12Sp")

school_waves <- function(sch) {
  if (sch %in% c1_early_schools) return(if (INCLUDE_GRADE9) 1:8 else 3:8)
  if (sch %in% c1_late_schools)  return(3:8)
  if (sch %in% c2_schools)       return(if (INCLUDE_GRADE9) 3:10 else 5:10)
  return(integer(0))
}

w_to_gp <- function(sch, w) {
  if (sch %in% c(c1_early_schools, c1_late_schools)) return(as.integer(w))
  if (sch %in% c2_schools) return(as.integer(w) - 2L)
  return(NA_integer_)
}

get_cohort <- function(sch) {
  if (sch %in% c1_early_schools) return("C1-early")
  if (sch %in% c1_late_schools)  return("C1-late")
  return("C2")
}

get_schtype <- function(sch) {
  if (sch %in% asian_schools)    return("Asian-majority")
  if (sch %in% hispanic_schools) return("Hispanic-majority")
  return("Other")
}

# ---------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------
get_filename <- function(w, type = "data") {
  if (type == "data") return(file.path(data_path, sprintf("w%d_adv_data.csv", w)))
  if (type == "edge") return(file.path(data_path, sprintf("w%dedges_clean.csv", w)))
}

build_school_adj <- function(edgelist, raw_ids, vertex_ids) {
  edges <- edgelist[edgelist$ego %in% raw_ids & edgelist$alter %in% raw_ids, ]
  n     <- length(raw_ids)
  mat   <- matrix(0L, nrow = n, ncol = n,
                  dimnames = list(vertex_ids, vertex_ids))
  if (nrow(edges) > 0) {
    for (k in seq_len(nrow(edges))) {
      i <- match(edges$ego[k],   raw_ids)
      j <- match(edges$alter[k], raw_ids)
      if (!is.na(i) & !is.na(j)) mat[i, j] <- 1L
    }
  }
  diag(mat) <- 0L
  mat
}

get_ecig_ever <- function(df, w) {
  v_bin <- paste0("w", w, "_life_use_sub_3a")
  v_ord <- paste0("w", w, "_life_use_sub_3")
  if (v_bin %in% names(df)) {
    x <- suppressWarnings(as.numeric(df[[v_bin]]))
    return(as.integer(x == 1))
  }
  if (v_ord %in% names(df)) {
    x <- df[[v_ord]]; x[x == ""] <- NA
    return(as.integer(suppressWarnings(as.numeric(x)) > 0))
  }
  rep(NA_integer_, nrow(df))
}

# Case-insensitive column lookup
find_col <- function(df, colname) {
  hits <- grep(paste0("^", colname, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) return(NULL)
  df[[hits[1]]]
}

# ---------------------------------------------------------------
# Load all waves
# ---------------------------------------------------------------
needed_waves <- sort(unique(unlist(lapply(all_schools, school_waves))))
cat("Loading data files...\n")
loaded_data <- list()
loaded_edge <- list()
for (w in needed_waves) {
  f_data <- get_filename(w, "data")
  f_edge <- get_filename(w, "edge")
  if (file.exists(f_data)) {
    loaded_data[[as.character(w)]] <- read.csv(f_data)
    cat(sprintf("  W%d data loaded\n", w))
  } else cat(sprintf("  W%d data NOT FOUND\n", w))
  if (file.exists(f_edge)) {
    loaded_edge[[as.character(w)]] <- read.csv(f_edge)
    cat(sprintf("  W%d edges loaded\n", w))
  } else cat(sprintf("  W%d edges NOT FOUND\n", w))
}

# ---------------------------------------------------------------
# Build diffnet per school
# vertex_ids format: "{school}_{record_id}" e.g. "101_101003"
# covariate_list uses the same format — consistent throughout
# ---------------------------------------------------------------
diffnet_list   <- list()
results_table  <- data.frame()
covariate_list <- list()   # id = "{school}_{record_id}", covariates from baseline wave

for (sch in all_schools) {

  sch_waves  <- school_waves(sch)
  sch_gps    <- sapply(sch_waves, function(w) w_to_gp(sch, w))
  sch_cohort <- get_cohort(sch)
  sch_type   <- get_schtype(sch)

  raw_ids_per_wave <- lapply(sch_waves, function(w) {
    wk <- as.character(w)
    if (!wk %in% names(loaded_data)) return(character(0))
    d  <- loaded_data[[wk]]
    as.character(d[[id_var]][d[[schoolid_var]] == sch & !is.na(d[[schoolid_var]])])
  })
  raw_ids <- sort(Reduce(union, raw_ids_per_wave))

  if (length(raw_ids) < 20) {
    cat(sprintf("School %d: n=%d < 20, skipped\n", sch, length(raw_ids)))
    next
  }
  n_students <- length(raw_ids)
  vertex_ids <- paste0(sch, "_", raw_ids)   # e.g. "101_101003"

  # Behavior matrix (rows=students, cols=grade periods)
  beh_mat <- matrix(NA_integer_, nrow = n_students, ncol = length(all_gps),
                    dimnames = list(vertex_ids, as.character(all_gps)))
  for (wi in seq_along(sch_waves)) {
    w  <- sch_waves[wi]
    gp <- sch_gps[wi]
    wk <- as.character(w)
    if (!wk %in% names(loaded_data)) next
    d   <- loaded_data[[wk]]
    idx <- match(raw_ids, as.character(d[[id_var]]))
    r   <- d[idx, ]
    beh_mat[, as.character(gp)] <- get_ecig_ever(r, w)
  }

  # Enforce monotonicity (ever-use cannot revert)
  for (i in seq_len(nrow(beh_mat))) {
    adopted <- FALSE
    for (gp_col in as.character(all_gps)) {
      if (!is.na(beh_mat[i, gp_col]) && beh_mat[i, gp_col] == 1) adopted <- TRUE
      if (adopted && !is.na(beh_mat[i, gp_col])) beh_mat[i, gp_col] <- 1L
    }
  }

  # toa = actual grade period number (3-8); NA = never adopted
  toa <- apply(beh_mat, 1, function(x) {
    first <- which(x == 1)
    if (length(first) == 0) return(NA_integer_)
    as.integer(all_gps[first[1]])
  })
  names(toa) <- vertex_ids

  # Network array (n x n x 6 grade periods)
  net_array <- array(0L, dim = c(n_students, n_students, length(all_gps)),
                     dimnames = list(vertex_ids, vertex_ids, as.character(all_gps)))
  for (wi in seq_along(sch_waves)) {
    w  <- sch_waves[wi]
    gp <- sch_gps[wi]
    wk <- as.character(w)
    if (!wk %in% names(loaded_edge)) next
    net_array[,, as.character(gp)] <- build_school_adj(loaded_edge[[wk]], raw_ids, vertex_ids)
  }

  # Create diffnet object
  dn <- tryCatch(
    as_diffnet(graph = net_array, toa = toa,
               t0 = min(all_gps), t1 = max(all_gps)),
    error = function(e) {
      cat(sprintf("School %d diffnet error: %s\n", sch, e$message)); NULL
    }
  )
  if (is.null(dn)) next
  diffnet_list[[as.character(sch)]] <- dn

  # Extract baseline covariates (from first wave of this school)
  baseline_w  <- sch_waves[1]
  baseline_wk <- as.character(baseline_w)
  if (baseline_wk %in% names(loaded_data)) {
    d_base  <- loaded_data[[baseline_wk]]
    idx_b   <- match(raw_ids, as.character(d_base[[id_var]]))
    r_base  <- d_base[idx_b, ]

    gen_r   <- find_col(d_base, paste0("w", baseline_w, "_dem_gender"))[idx_b]
    eth_r   <- find_col(d_base, paste0("w", baseline_w, "_eth"))[idx_b]
    race_r  <- find_col(d_base, paste0("w", baseline_w, "_race"))[idx_b]
    sex_r   <- find_col(d_base, paste0("w", baseline_w, "_dem_sexuality"))[idx_b]

    female   <- as.integer(!is.na(gen_r)  & gen_r  == 2)
    hispanic <- as.integer(!is.na(eth_r)  & eth_r  == 1)
    asian    <- as.integer(!is.na(race_r) & suppressWarnings(as.numeric(race_r)) == 2)
    sex_raw2 <- suppressWarnings(as.numeric(sex_r))
    sex_min  <- as.integer(!is.na(sex_raw2) & sex_raw2 != 1)

    covariate_list[[as.character(sch)]] <- data.frame(
      id       = vertex_ids,   # "{school}_{record_id}" format
      female   = female,
      hispanic = hispanic,
      asian    = asian,
      sex_min  = sex_min,
      stringsAsFactors = FALSE
    )
  }

  n_adopters <- sum(!is.na(toa))
  results_table <- rbind(results_table, data.frame(
    school     = sch,
    cohort     = sch_cohort,
    sch_type   = sch_type,
    n_students = n_students,
    n_adopters = n_adopters,
    final_prev = round(n_adopters / n_students * 100, 1),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("School %d (%s, %s): n=%d, adopters=%d (%.1f%%)\n",
              sch, sch_cohort, sch_type, n_students, n_adopters,
              n_adopters / n_students * 100))
}

# ---------------------------------------------------------------
# Combine all schools into pooled diffnet
# NOTE: do.call(c, diffnet_list) adds a "{key}." prefix to vertex IDs
# e.g. "101_101003" becomes "101.101_101003"
# We strip this prefix wherever needed using sub("^[^.]+\\.", "", id)
# ---------------------------------------------------------------
cat("\nCombining schools...\n")
diffnet_all <- tryCatch(
  do.call(c, diffnet_list),
  error = function(e) { cat("Combine error:", e$message, "\n"); NULL }
)
if (is.null(diffnet_all)) stop("diffnet_all is NULL — check errors above.")

cat("\n===== POOLED DIFFUSION SUMMARY =====\n")
print(summary(diffnet_all))

# Build covariate lookup using ORIGINAL (unprefixed) IDs
cov_df <- do.call(rbind, covariate_list)   # id = "101_101003"

# Build school/cohort/type lookup using ORIGINAL IDs
attrs_df <- do.call(rbind, lapply(names(diffnet_list), function(k) {
  sch_int <- as.integer(k)
  data.frame(
    id       = names(diffnet_list[[k]]$toa),   # original "101_101003"
    school   = sch_int,
    cohort   = get_cohort(sch_int),
    sch_type = get_schtype(sch_int),
    stringsAsFactors = FALSE
  )
}))
attrs_df <- merge(attrs_df, cov_df, by = "id", all.x = TRUE)

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
pdf(file.path(data_path, "0409_netdiffuseR_adopters.pdf"), width = 9, height = 5)
plot_adopters(diffnet_all,
              main = "E-cigarette Ever-Use: Cumulative Adoption by Grade Period\n(Cohorts 1 & 2)",
              xlab = "Grade Period", ylab = "Proportion")
axis(1, at = all_gps, labels = gp_labels, tick = FALSE, line = 1.5, cex.axis = 0.75)
dev.off()
cat("Saved: adopters.pdf\n")

pdf(file.path(data_path, "0409_netdiffuseR_hazard.pdf"), width = 9, height = 5)
hazard_rate(diffnet_all,
            main = "Hazard Rate of E-cigarette Initiation by Grade Period\n(Cohorts 1 & 2)",
            xlab = "Grade Period")
axis(1, at = all_gps, labels = gp_labels, tick = FALSE, line = 1.5, cex.axis = 0.75)
dev.off()
cat("Saved: hazard.pdf\n")

thr_vals <- threshold(diffnet_all)
toa_vals <- diffnet_all$toa
thr_data <- data.frame(toa = as.integer(toa_vals), threshold = as.numeric(thr_vals))
thr_data <- thr_data[!is.na(thr_data$toa) & !is.na(thr_data$threshold), ]
present_gps    <- sort(unique(thr_data$toa))
present_labels <- gp_labels[present_gps - min(all_gps) + 1]
box_cols <- c("10Fa"="#f4a261","10Sp"="#f4a261","11Fa"="#e76f51",
              "11Sp"="#e76f51","12Fa"="#c0392b","12Sp"="#c0392b")

pdf(file.path(data_path, "0409_netdiffuseR_threshold_boxplot.pdf"), width = 9, height = 5)
par(mar = c(5, 4, 4, 2))
boxplot(threshold ~ toa, data = thr_data,
        names  = present_labels,
        xlab   = "Grade Period of Adoption",
        ylab   = "Threshold (fraction of friends already using)",
        main   = "Adoption Threshold by Grade Period\n(Cohorts 1 & 2)",
        col    = box_cols[present_labels],
        ylim   = c(0, 1))
abline(h = 0.5, lty = 2, col = "red", lwd = 1.5)
n_per_gp <- table(thr_data$toa)
for (i in seq_along(n_per_gp))
  mtext(paste0("n=", n_per_gp[i]), side = 1, at = i, line = 3.5, cex = 0.7, col = "gray40")
dev.off()
cat("Saved: threshold_boxplot.pdf\n")

all_colors <- c(
  "101"="#1f77b4","102"="#d62728","103"="#2ca02c","104"="#9467bd","105"="#8c564b",
  "106"="#e377c2","107"="#ff7f0e","112"="#17becf","113"="#bcbd22","114"="#7f7f7f",
  "201"="#aec7e8","212"="#ffbb78","213"="#98df8a","214"="#ff9896"
)
pdf(file.path(data_path, "0409_netdiffuseR_byschool.pdf"), width = 11, height = 6)
par(mar = c(5, 4, 4, 10))
plot(NULL, xlim = c(min(all_gps), max(all_gps)), ylim = c(0, 0.65),
     xlab = "Grade Period", ylab = "Cumulative Adoption Proportion",
     main = "E-cigarette Cumulative Adoption by School\n(solid=C1, dashed=C2)", xaxt = "n")
axis(1, at = all_gps, labels = gp_labels, cex.axis = 0.85)
for (sch_key in names(diffnet_list)) {
  dn_sch   <- diffnet_list[[sch_key]]
  toa_sch  <- dn_sch$toa
  n_sch    <- length(toa_sch)
  cum_prev <- sapply(all_gps, function(gp) sum(!is.na(toa_sch) & toa_sch <= gp) / n_sch)
  col_val  <- ifelse(is.na(all_colors[sch_key]), "black", all_colors[sch_key])
  lty_val  <- ifelse(as.integer(sch_key) >= 200, 2, 1)
  lines(all_gps, cum_prev, col = col_val, lwd = 2, lty = lty_val)
  points(all_gps, cum_prev, col = col_val, pch = 16, cex = 0.7)
}
legend("topright", inset = c(-0.32, 0), xpd = TRUE,
       legend = paste0(names(diffnet_list),
                       ifelse(as.integer(names(diffnet_list)) >= 200, " (C2)", " (C1)")),
       col = ifelse(is.na(all_colors[names(diffnet_list)]), "black",
                    all_colors[names(diffnet_list)]),
       lwd = 2, pch = 16, cex = 0.65, bty = "n")
dev.off()
cat("Saved: byschool.pdf\n")

# ---------------------------------------------------------------
# Node-level outcomes table
# Uses ORIGINAL (unprefixed) IDs throughout
# ---------------------------------------------------------------
node_df <- data.frame(
  id        = names(toa_vals),                    # prefixed "101.101_101003"
  id_orig   = sub("^[^.]+\\.", "", names(toa_vals)),  # original "101_101003"
  toa       = as.integer(toa_vals),
  threshold = as.numeric(thr_vals),
  stringsAsFactors = FALSE
)
node_df <- merge(node_df, attrs_df, by.x = "id_orig", by.y = "id", all.x = TRUE)

adopters <- node_df[!is.na(node_df$toa), ]
cat(sprintf("\nTotal adopters: %d\n", nrow(adopters)))
cat(sprintf("Adopters with covariate data: %d (female non-NA)\n",
            sum(!is.na(adopters$female))))

# ---------------------------------------------------------------
# Task (3a): Threshold regression with covariates
# ---------------------------------------------------------------
cat("\n===== THRESHOLD REGRESSION =====\n")
thr_sub <- adopters[!is.na(adopters$threshold), ]
cat(sprintf("N with threshold data: %d\n", nrow(thr_sub)))

# Build formula: include sex_min only if enough non-NA observations
thr_covs <- c("factor(toa)", "female", "hispanic", "asian")
if (sum(!is.na(thr_sub$sex_min)) > 30 &&
    length(unique(na.omit(thr_sub$sex_min))) > 1) {
  thr_covs <- c(thr_covs, "sex_min")
  cat("sex_min included in threshold model\n")
} else {
  cat("sex_min excluded (insufficient data)\n")
}
fml_thr <- as.formula(paste("threshold ~", paste(thr_covs, collapse = " + ")))
cat("Formula:", deparse(fml_thr), "\n")

fit_thr <- tryCatch(
  glm(fml_thr, data = thr_sub[complete.cases(thr_sub[, c("threshold","female","hispanic","asian")]),],
      family = quasibinomial(link = "logit")),
  error = function(e) { cat("Threshold GLM error:", e$message, "\n"); NULL }
)
thr_result <- NULL
if (!is.null(fit_thr)) {
  ct_thr <- summary(fit_thr)$coefficients
  rownames(ct_thr) <- gsub("factor\\(toa\\)", "GP", rownames(ct_thr))
  print(round(ct_thr, 4))
  thr_result <- data.frame(
    parameter = rownames(ct_thr),
    estimate  = round(ct_thr[,1], 4),
    se        = round(ct_thr[,2], 4),
    t_value   = round(ct_thr[,3], 4),
    pval      = round(ct_thr[,4], 6),
    stringsAsFactors = FALSE
  )
  write.csv(thr_result,
            file.path(data_path, "0409_netdiffuseR_threshold_regression.csv"),
            row.names = FALSE)
  cat("Saved: threshold_regression.csv\n")
}

# Save adopters-only node outcomes
write.csv(adopters,
          file.path(data_path, "0409_netdiffuseR_node_outcomes_adopters.csv"),
          row.names = FALSE)
cat("Saved: node_outcomes_adopters.csv\n")

# ---------------------------------------------------------------
# Exposure regression
# KEY: strip the do.call(c,...) prefix from toa_vals names
# so IDs match the unprefixed covariate_list IDs
# ---------------------------------------------------------------
cat("\nCalculating network exposure...\n")
expo      <- exposure(diffnet_all)
adopt_mat <- diffnet_all$cumadopt
n_period  <- ncol(adopt_mat)

# Strip prefix: "101.101_101003" -> "101_101003"
vertex_ids_clean <- sub("^[^.]+\\.", "", names(toa_vals))

reg_rows <- list()
for (gi in 1:(n_period - 1)) {
  not_yet <- which(adopt_mat[, gi] == 0)
  reg_rows[[gi]] <- data.frame(
    id           = vertex_ids_clean[not_yet],
    adopt_next   = as.integer(adopt_mat[not_yet, gi + 1]),
    exposure     = as.numeric(expo[not_yet, gi]),
    grade_period = all_gps[gi],
    stringsAsFactors = FALSE
  )
}
reg_data <- do.call(rbind, reg_rows)

# Merge with covariates using original IDs (both sides now match)
reg_data <- merge(reg_data, cov_df, by = "id", all.x = TRUE)

# Also add school type and cohort
reg_data <- merge(reg_data,
                  attrs_df[, c("id","cohort","sch_type")],
                  by = "id", all.x = TRUE)

cat(sprintf("Rows before filtering: %d\n", nrow(reg_data)))
cat(sprintf("female NAs: %d\n",   sum(is.na(reg_data$female))))
cat(sprintf("hispanic NAs: %d\n", sum(is.na(reg_data$hispanic))))
cat(sprintf("sex_min NAs: %d\n",  sum(is.na(reg_data$sex_min))))

# Build model: include sex_min if enough data
exp_covs <- c("female", "hispanic", "asian")
if (sum(!is.na(reg_data$sex_min)) > 100 &&
    length(unique(na.omit(reg_data$sex_min))) > 1) {
  exp_covs <- c(exp_covs, "sex_min")
  cat("sex_min included in exposure model\n")
} else {
  cat("sex_min excluded (insufficient data)\n")
}

# Filter to complete cases on required variables
req_vars <- c("adopt_next", "exposure", "grade_period", "female", "hispanic", "asian")
reg_data_final <- reg_data[complete.cases(reg_data[, req_vars]), ]
cat(sprintf("Exposure regression: %d person-period obs\n", nrow(reg_data_final)))

fml_exp <- as.formula(paste(
  "adopt_next ~ exposure + factor(grade_period) +",
  paste(exp_covs, collapse = " + ")
))
cat("Formula:", deparse(fml_exp), "\n")

fit_exp <- glm(fml_exp, data = reg_data_final, family = binomial)
cat("\n===== EXPOSURE LOGISTIC REGRESSION =====\n")
ct <- summary(fit_exp)$coefficients
rownames(ct) <- gsub("factor\\(grade_period\\)", "GP", rownames(ct))
print(round(ct, 4))

or_exp <- exp(coef(fit_exp)["exposure"])
ci_exp <- exp(confint(fit_exp, "exposure"))
cat(sprintf("\nExposure OR = %.3f (95%% CI: %.3f - %.3f), p = %.4f\n",
            or_exp, ci_exp[1], ci_exp[2], ct["exposure", 4]))

exp_result <- data.frame(
  parameter = rownames(ct),
  estimate  = round(ct[,1], 4),
  se        = round(ct[,2], 4),
  z         = round(ct[,3], 4),
  pval      = round(ct[,4], 6),
  OR        = round(exp(ct[,1]), 4),
  stringsAsFactors = FALSE
)
write.csv(exp_result,
          file.path(data_path, "0409_netdiffuseR_exposure_regression.csv"),
          row.names = FALSE)
cat("Saved: exposure_regression.csv\n")

# ---------------------------------------------------------------
# Adopter classification
# ---------------------------------------------------------------
cat("\n===== ADOPTER CLASSIFICATION =====\n")
cls <- classify(diffnet_all, include_censored = TRUE)
print(ftable(cls))
cls_df <- as.data.frame(ftable(cls))
write.csv(cls_df,
          file.path(data_path, "0409_netdiffuseR_classification.csv"),
          row.names = FALSE)
cat("Saved: classification.csv\n")

# ---------------------------------------------------------------
# Per-school summary
# ---------------------------------------------------------------
write.csv(results_table,
          file.path(data_path, "0409_netdiffuseR_school_summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------
# Task (4): Combine all CSVs into one Excel workbook
# ---------------------------------------------------------------
cat("\nBuilding combined Excel workbook...\n")
wb <- createWorkbook()

read_if_exists <- function(path) {
  if (file.exists(path)) read.csv(path) else data.frame(note = "File not found")
}

sheets <- list(
  school_summary       = "0409_netdiffuseR_school_summary.csv",
  exposure_regression  = "0409_netdiffuseR_exposure_regression.csv",
  threshold_regression = "0409_netdiffuseR_threshold_regression.csv",
  classification       = "0409_netdiffuseR_classification.csv",
  node_outcomes        = "0409_netdiffuseR_node_outcomes_adopters.csv"
)
for (sheet_name in names(sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name,
            read_if_exists(file.path(data_path, sheets[[sheet_name]])))
}
saveWorkbook(wb,
             file.path(data_path, "0409_netdiffuseR_ALL_RESULTS.xlsx"),
             overwrite = TRUE)
cat("Saved: 0409_netdiffuseR_ALL_RESULTS.xlsx\n")


# ---------------------------------------------------------------
# Task (5): plot_diffnet() — network graphs colored by adoption status
# One PDF per school, each page = one grade period
# Nodes: grey = non-adopter, red = adopter
# ---------------------------------------------------------------
cat("\nGenerating plot_diffnet() network graphs...\n")

plot_dir <- file.path(data_path, "0409_netdiffuseR_network_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

for (sch_key in names(diffnet_list)) {
  dn_sch   <- diffnet_list[[sch_key]]
  sch_int  <- as.integer(sch_key)
  sch_type <- get_schtype(sch_int)
  cohort   <- get_cohort(sch_int)
  n_nodes  <- length(dn_sch$toa)

  # Skip very large networks (plot would be unreadable)
  if (n_nodes > 500) {
    cat(sprintf("  School %s: n=%d too large for network plot, skipping\n",
                sch_key, n_nodes))
    next
  }

  out_file <- file.path(plot_dir,
                         sprintf("0409_diffnet_school%s.pdf", sch_key))
  pdf(out_file, width = 12, height = 10)

  plot_diffnet(dn_sch, slices = seq_along(all_gps), vertex.color = c("grey70", "tomato"))

  dev.off()
  cat(sprintf("  Saved: 0409_diffnet_school%s.pdf\n", sch_key))
}

# Also plot one representative school per type for the paper
# Pick smallest school in each type (easier to see patterns)
representative <- list(
  "Asian-majority"    = "103",   # smallest Asian school
  "Hispanic-majority" = "102",   # smallest Hispanic school
  "Other"             = "201"    # smallest Other school
)

pdf(file.path(data_path, "0409_netdiffuseR_diffnet_representative.pdf"),
    width = 14, height = 10)

for (label in names(representative)) {
  sch_key <- representative[[label]]
  if (!sch_key %in% names(diffnet_list)) next
  dn_sch <- diffnet_list[[sch_key]]

  plot_diffnet(dn_sch, slices = seq_along(all_gps), vertex.color = c("grey70", "tomato"))
}

dev.off()
cat("Saved: 0409_netdiffuseR_diffnet_representative.pdf\n")


# ---------------------------------------------------------------
# Save R objects
# ---------------------------------------------------------------
save(diffnet_list, diffnet_all, results_table, node_df,
     file = file.path(data_path, "0409_netdiffuseR_everuse.RData"))

cat("\n===== DONE =====\n")
cat("Output folder:", data_path, "\n")
cat("PDFs  : adopters / hazard / threshold_boxplot / byschool\n")
cat("        diffnet_representative / network_plots/ (per school)\n")
cat("CSVs  : school_summary / exposure_regression / threshold_regression\n")
cat("        classification / node_outcomes_adopters\n")
cat("Excel : 0409_netdiffuseR_ALL_RESULTS.xlsx\n")
cat("RData : 0409_netdiffuseR_everuse.RData\n")

rm(list = ls())
library(netdiffuseR)
library(openxlsx)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned_Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

# School definitions
c1_early_schools <- c(101, 102, 103, 104, 105)
c1_late_schools  <- c(106, 107, 112, 113, 114)
c2_schools       <- c(201, 212, 213, 214)
exclude_schools  <- c(108)

all_schools      <- setdiff(c(c1_early_schools, c1_late_schools, c2_schools), exclude_schools)
asian_schools    <- c(103, 105, 112, 113, 212, 213)
hispanic_schools <- c(102, 106, 107, 114, 214)

# Grade-wave mapping
# C1 (101-114): W3-W8  (GP3=10Fa ... GP8=12Sp)
# C2 (201-214): W5-W10 (GP3=10Fa ... GP8=12Sp)
#
# WHY NOT GP1 OR GP2?
#   GP1 = Grade 9 Fall, GP2 = Grade 9 Spring.
#   C1 late schools (106-114) only entered the study at W3 = Grade 10 Fall.
#   C2 schools (201-214) entered at W5 = Grade 10 Fall (mapped to GP3).
#   To align BOTH cohorts across the same six grade periods,
#   we must start at Grade 10 Fall (GP3). Including Grade 9 would mean
#   C1 late schools and all C2 schools have no data for those periods,
#   making diffusion analysis impossible across the full sample.
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
  if (sch %in% c(c1_early_schools, c1_late_schools)) return("C1")
  return("C2")
}

get_schtype <- function(sch) {
  if (sch %in% asian_schools)    return("Asian-majority")
  if (sch %in% hispanic_schools) return("Hispanic-majority")
  return("Other")
}

get_filename <- function(w, type = "data") {
  if (type == "data") return(file.path(data_path, sprintf("w%d_adv_data.csv", w)))
  if (type == "edge") return(file.path(data_path, sprintf("w%dedges_clean.csv", w)))
}

build_school_adj <- function(edgelist, raw_ids, vertex_ids) {
  edges <- edgelist[edgelist$ego %in% raw_ids & edgelist$alter %in% raw_ids, ]
  n     <- length(raw_ids)
  mat   <- matrix(0L, nrow = n, ncol = n, dimnames = list(vertex_ids, vertex_ids))
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

get_ecig_past6mo <- function(df, w) {
  v <- paste0("w", w, "_past_6mo_use_3")
  if (v %in% names(df)) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    return(as.integer(x == 1))
  }
  rep(NA_integer_, nrow(df))
}

get_gad <- function(df, w) {
  v <- paste0("w", w, "_rcads_gad_mean")
  hits <- grep(paste0("^", v, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(suppressWarnings(as.numeric(df[[hits[1]]])))
  rep(NA_real_, nrow(df))
}

get_mdd <- function(df, w) {
  v <- paste0("w", w, "_rcads_mdd_mean")
  hits <- grep(paste0("^", v, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(suppressWarnings(as.numeric(df[[hits[1]]])))
  rep(NA_real_, nrow(df))
}

get_friends_ecig <- function(df, w) {
  v <- paste0("w", w, "_friends_use_ecig")
  hits <- grep(paste0("^", v, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) {
    x <- suppressWarnings(as.numeric(df[[hits[1]]]))
    x[x == 6] <- NA
    return(x)
  }
  rep(NA_real_, nrow(df))
}

find_col <- function(df, colname) {
  hits <- grep(paste0("^", colname, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) return(NULL)
  df[[hits[1]]]
}

# Load all waves
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

# Build diffnet per school
diffnet_list   <- list()
results_table  <- data.frame()
covariate_list <- list()
timevar_list   <- list()

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
  vertex_ids <- paste0(sch, "_", raw_ids)

  # Behavior matrix
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
    beh_mat[, as.character(gp)] <- get_ecig_past6mo(r, w)
  }

  # Time-varying covariates
  tv_rows <- list()
  for (wi in seq_along(sch_waves)) {
    w  <- sch_waves[wi]
    gp <- sch_gps[wi]
    wk <- as.character(w)
    if (!wk %in% names(loaded_data)) next
    d   <- loaded_data[[wk]]
    idx <- match(raw_ids, as.character(d[[id_var]]))
    r   <- d[idx, ]
    tv_rows[[wi]] <- data.frame(
      id           = vertex_ids,
      grade_period = gp,
      gad          = get_gad(r, w),
      mdd          = get_mdd(r, w),
      friends_ecig = get_friends_ecig(r, w),
      stringsAsFactors = FALSE
    )
  }
  timevar_list[[as.character(sch)]] <- do.call(rbind, tv_rows)

  # Enforce monotonicity
  for (i in seq_len(nrow(beh_mat))) {
    adopted <- FALSE
    for (gp_col in as.character(all_gps)) {
      if (!is.na(beh_mat[i, gp_col]) && beh_mat[i, gp_col] == 1) adopted <- TRUE
      if (adopted && !is.na(beh_mat[i, gp_col])) beh_mat[i, gp_col] <- 1L
    }
  }

  toa <- apply(beh_mat, 1, function(x) {
    first <- which(x == 1)
    if (length(first) == 0) return(NA_integer_)
    as.integer(all_gps[first[1]])
  })
  names(toa) <- vertex_ids

  net_array <- array(0L, dim = c(n_students, n_students, length(all_gps)),
                     dimnames = list(vertex_ids, vertex_ids, as.character(all_gps)))
  for (wi in seq_along(sch_waves)) {
    w  <- sch_waves[wi]
    gp <- sch_gps[wi]
    wk <- as.character(w)
    if (!wk %in% names(loaded_edge)) next
    net_array[,, as.character(gp)] <- build_school_adj(loaded_edge[[wk]], raw_ids, vertex_ids)
  }

  dn <- tryCatch(
    as_diffnet(graph = net_array, toa = toa, t0 = min(all_gps), t1 = max(all_gps)),
    error = function(e) { cat(sprintf("School %d diffnet error: %s\n", sch, e$message)); NULL }
  )
  if (is.null(dn)) next
  diffnet_list[[as.character(sch)]] <- dn

  # Baseline covariates
  # dem_gender: 0=Female, 1=Male, 3=Prefer not to disclose
  #
  # par_edu — three-step harmonisation (aligned with Anibal's disadoption study):
  #
  #   Step 1 — Read raw value at every wave
  #
  #   Step 2 — Harmonise W7-W10 9-level scale → legacy 7-level scale
  #     New scale (W7-W10): 1  2  3  4  5  6  7  8  9
  #     Legacy scale (W1-W6):1  2  3  4  4  4  5  6  NA
  #     (new 5/6 = vocational/associate → legacy 4 "some college";
  #      new 7 = bachelor's → legacy 5 "college graduate";
  #      new 8 = master's/doctoral → legacy 6 "advanced degree";
  #      new 9 = don't know → NA; legacy 7 = don't know → NA)
  #     Codes 1-4 are identical across both scales.
  #
  #   Step 3 — LOCF (Last Observation Carried Forward) per student
  #     Treat par_edu as slowly-varying within student.
  #     Carry most recent valid value forward to fill NAs.
  #     Observed values are never replaced by earlier observations.
  #     After LOCF, any remaining NA → school modal value (last resort).

  remap_par_edu <- function(x, w) {
    x <- suppressWarnings(as.numeric(x))
    if (w >= 7) {
      new_to_legacy <- c(1, 2, 3, 4, 4, 4, 5, 6, NA)
      x_out <- ifelse(!is.na(x) & x >= 1 & x <= 9, new_to_legacy[x], NA_real_)
    } else {
      x_out <- ifelse(!is.na(x) & x == 7, NA_real_, x)
    }
    x_out
  }

  baseline_w  <- sch_waves[1]
  baseline_wk <- as.character(baseline_w)
  if (baseline_wk %in% names(loaded_data)) {
    d_base <- loaded_data[[baseline_wk]]
    idx_b  <- match(raw_ids, as.character(d_base[[id_var]]))
    gen_r  <- find_col(d_base, paste0("w", baseline_w, "_dem_gender"))[idx_b]
    eth_r  <- find_col(d_base, paste0("w", baseline_w, "_eth"))[idx_b]
    race_r <- find_col(d_base, paste0("w", baseline_w, "_race"))[idx_b]
    sex_r  <- find_col(d_base, paste0("w", baseline_w, "_dem_sexuality"))[idx_b]

    female   <- as.integer(!is.na(gen_r) & suppressWarnings(as.numeric(gen_r)) == 0)
    hispanic <- as.integer(!is.na(eth_r)  & eth_r  == 1)
    asian    <- as.integer(!is.na(race_r) & suppressWarnings(as.numeric(race_r)) == 2)
    sex_raw2 <- suppressWarnings(as.numeric(sex_r))
    sex_min  <- as.integer(!is.na(sex_raw2) & sex_raw2 != 1)

    # Step 1+2: build harmonised par_edu matrix across all school waves
    par_matrix <- matrix(NA_real_, nrow = n_students, ncol = length(sch_waves))
    for (wi in seq_along(sch_waves)) {
      w  <- sch_waves[wi]
      wk <- as.character(w)
      if (!wk %in% names(loaded_data)) next
      d_w   <- loaded_data[[wk]]
      idx_w <- match(raw_ids, as.character(d_w[[id_var]]))
      raw_w <- find_col(d_w, paste0("w", w, "_dem_high_par_edu"))[idx_w]
      par_matrix[, wi] <- remap_par_edu(raw_w, w)
    }

    # Step 3: LOCF — carry forward most recent valid value per student
    par_locf <- par_matrix
    for (wi in 2:ncol(par_locf)) {
      still_na <- is.na(par_locf[, wi])
      par_locf[still_na, wi] <- par_locf[still_na, wi - 1]
    }

    # Take baseline value (first wave column) after LOCF
    par_edu <- par_locf[, 1]

    # Last resort: school modal value for any remaining NA
    n_still_na <- sum(is.na(par_edu))
    if (n_still_na > 0) {
      valid_vals <- par_edu[!is.na(par_edu)]
      if (length(valid_vals) > 0) {
        modal_val <- as.numeric(names(sort(table(valid_vals), decreasing = TRUE))[1])
        par_edu[is.na(par_edu)] <- modal_val
        cat(sprintf("  School %d: %d par_edu still NA after LOCF → school modal = %g\n",
                    sch, n_still_na, modal_val))
      }
    }

    cat(sprintf("  School %d: par_edu range [%g, %g], NAs remaining = %d\n",
                sch, min(par_edu, na.rm = TRUE), max(par_edu, na.rm = TRUE),
                sum(is.na(par_edu))))

    covariate_list[[as.character(sch)]] <- data.frame(
      id = vertex_ids, female = female, hispanic = hispanic,
      asian = asian, sex_min = sex_min, par_edu = par_edu,
      stringsAsFactors = FALSE
    )
  }

  n_adopters <- sum(!is.na(toa))
  results_table <- rbind(results_table, data.frame(
    school = sch, cohort = sch_cohort, sch_type = sch_type,
    n_students = n_students, n_adopters = n_adopters,
    final_prev = round(n_adopters / n_students * 100, 1),
    stringsAsFactors = FALSE
  ))
  cat(sprintf("School %d (%s, %s): n=%d, adopters=%d (%.1f%%)\n",
              sch, sch_cohort, sch_type, n_students, n_adopters,
              n_adopters / n_students * 100))
}

cat("\nCombining schools...\n")
diffnet_all <- tryCatch(
  do.call(c, diffnet_list),
  error = function(e) { cat("Combine error:", e$message, "\n"); NULL }
)
if (is.null(diffnet_all)) stop("diffnet_all is NULL — check errors above.")
print(summary(diffnet_all))

cov_df <- do.call(rbind, covariate_list)

attrs_df <- do.call(rbind, lapply(names(diffnet_list), function(k) {
  sch_int <- as.integer(k)
  data.frame(id = names(diffnet_list[[k]]$toa), school = sch_int,
             cohort = get_cohort(sch_int), sch_type = get_schtype(sch_int),
             stringsAsFactors = FALSE)
}))
attrs_df <- merge(attrs_df, cov_df, by = "id", all.x = TRUE)

timevar_df <- do.call(rbind, timevar_list)

# Summary plots
pdf(file.path(data_path, "0511_netdiffuseR_past6mo_adopters.pdf"), width = 9, height = 5)
plot_adopters(diffnet_all,
              main = "E-cigarette Past-6-Month Use: Cumulative Adoption by Grade Period\n(Cohorts 1 & 2)",
              xlab = "Grade Period", ylab = "Proportion")
axis(1, at = all_gps, labels = gp_labels, tick = FALSE, line = 1.5, cex.axis = 0.75)
dev.off()

pdf(file.path(data_path, "0511_netdiffuseR_past6mo_hazard.pdf"), width = 9, height = 5)
hazard_rate(diffnet_all,
            main = "Hazard Rate of E-cigarette Past-6-Month Use by Grade Period\n(Cohorts 1 & 2)",
            xlab = "Grade Period")
axis(1, at = all_gps, labels = gp_labels, tick = FALSE, line = 1.5, cex.axis = 0.75)
dev.off()

thr_vals <- threshold(diffnet_all)
toa_vals <- diffnet_all$toa
thr_data <- data.frame(toa = as.integer(toa_vals), threshold = as.numeric(thr_vals))
thr_data <- thr_data[!is.na(thr_data$toa) & !is.na(thr_data$threshold), ]
present_gps    <- sort(unique(thr_data$toa))
present_labels <- gp_labels[present_gps - min(all_gps) + 1]
box_cols <- c("10Fa"="#f4a261","10Sp"="#f4a261","11Fa"="#e76f51",
              "11Sp"="#e76f51","12Fa"="#c0392b","12Sp"="#c0392b")

pdf(file.path(data_path, "0511_netdiffuseR_past6mo_threshold_boxplot.pdf"), width = 9, height = 5)
par(mar = c(5, 4, 4, 2))
boxplot(threshold ~ toa, data = thr_data, names = present_labels,
        xlab = "Grade Period of Adoption",
        ylab = "Threshold (fraction of friends already using)",
        main = "Adoption Threshold by Grade Period\n(Cohorts 1 & 2)",
        col = box_cols[present_labels], ylim = c(0, 1))
abline(h = 0.5, lty = 2, col = "red", lwd = 1.5)
n_per_gp <- table(thr_data$toa)
for (i in seq_along(n_per_gp))
  mtext(paste0("n=", n_per_gp[i]), side = 1, at = i, line = 3.5, cex = 0.7, col = "gray40")
dev.off()

all_colors <- c(
  "101"="#1f77b4","102"="#d62728","103"="#2ca02c","104"="#9467bd","105"="#8c564b",
  "106"="#e377c2","107"="#ff7f0e","112"="#17becf","113"="#bcbd22","114"="#7f7f7f",
  "201"="#aec7e8","212"="#ffbb78","213"="#98df8a","214"="#ff9896"
)
pdf(file.path(data_path, "0511_netdiffuseR_past6mo_byschool.pdf"), width = 11, height = 6)
par(mar = c(5, 4, 4, 10))
plot(NULL, xlim = c(min(all_gps), max(all_gps)), ylim = c(0, 0.65),
     xlab = "Grade Period", ylab = "Cumulative Adoption Proportion",
     main = "E-cigarette Past-6-Month Use by School\n(solid=C1, dashed=C2)", xaxt = "n")
axis(1, at = all_gps, labels = gp_labels, cex.axis = 0.85)
for (sch_key in names(diffnet_list)) {
  dn_sch  <- diffnet_list[[sch_key]]
  toa_sch <- dn_sch$toa
  n_sch   <- length(toa_sch)
  cum_prev <- sapply(all_gps, function(gp) sum(!is.na(toa_sch) & toa_sch <= gp) / n_sch)
  col_val <- ifelse(is.na(all_colors[sch_key]), "black", all_colors[sch_key])
  lty_val <- ifelse(as.integer(sch_key) >= 200, 2, 1)
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

# Node outcomes
node_df <- data.frame(
  id      = names(toa_vals),
  id_orig = sub("^[^.]+\\.", "", names(toa_vals)),
  toa     = as.integer(toa_vals),
  threshold = as.numeric(thr_vals),
  stringsAsFactors = FALSE
)
node_df  <- merge(node_df, attrs_df, by.x = "id_orig", by.y = "id", all.x = TRUE)
adopters <- node_df[!is.na(node_df$toa), ]
cat(sprintf("\nTotal adopters: %d\n", nrow(adopters)))

# Bivariate threshold analyses (per 学姐's suggestion)
# Run before full model to check if individual predictors are significant
# If bivariate significant but full model not → possible multicollinearity
run_threshold_bivariate <- function(thr_sub, mh_var) {
  cat("\n===== THRESHOLD BIVARIATE ANALYSES =====\n")
  predictors <- c("friends_ecig", "hispanic", "asian", "sex_min",
                  "cohort", "female", "par_edu", mh_var)
  results <- list()
  for (pred in predictors) {
    if (!pred %in% names(thr_sub)) next
    if (sum(!is.na(thr_sub[[pred]])) < 10) next
    fml <- as.formula(paste("threshold ~", pred))
    fit <- tryCatch(
      glm(fml, data = thr_sub[!is.na(thr_sub[[pred]]) & !is.na(thr_sub$threshold), ],
          family = quasibinomial(link = "logit")),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    ct <- summary(fit)$coefficients
    if (nrow(ct) < 2) next
    cat(sprintf("  %-15s  estimate=%.3f  OR=%.3f  p=%.4f\n",
                pred, ct[2,1], exp(ct[2,1]), ct[2,4]))
    results[[pred]] <- data.frame(
      predictor = pred,
      estimate  = round(ct[2,1], 4),
      OR        = round(exp(ct[2,1]), 4),
      pval      = round(ct[2,4], 6)
    )
  }
  do.call(rbind, results)
}

# Threshold regression (full model)
# Includes: female, par_edu, factor(school) as fixed effect (per Dr. Valente)
run_threshold_reg <- function(thr_sub, mh_var, label) {
  cat(sprintf("\n===== THRESHOLD REGRESSION — FULL MODEL (%s) =====\n", label))
  thr_covs <- c("factor(toa)", "cohort", "female", "hispanic", "asian",
                "par_edu", "factor(school)", mh_var, "friends_ecig")
  if (sum(!is.na(thr_sub$sex_min)) > 30 && length(unique(na.omit(thr_sub$sex_min))) > 1)
    thr_covs <- c(thr_covs, "sex_min")
  fml <- as.formula(paste("threshold ~", paste(thr_covs, collapse = " + ")))
  cat("Formula:", deparse(fml), "\n")
  req <- c("threshold", "cohort", "female", "hispanic", "asian", "par_edu", "school", mh_var, "friends_ecig")
  fit_data <- thr_sub[complete.cases(thr_sub[, req[req %in% names(thr_sub)]]), ]
  fit <- tryCatch(
    glm(fml, data = fit_data, family = quasibinomial(link = "logit")),
    error = function(e) { cat("Error:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coefficients
  rownames(ct) <- gsub("factor\\(toa\\)", "GP", rownames(ct))
  rownames(ct) <- gsub("factor\\(school\\)", "School", rownames(ct))
  out <- data.frame(
    parameter = rownames(ct),
    estimate  = round(ct[,1], 4),
    se        = round(ct[,2], 4),
    t_value   = round(ct[,3], 4),
    pval      = round(ct[,4], 6),
    OR        = round(exp(ct[,1]), 4),
    stringsAsFactors = FALSE
  )
  print(out)

  # VIF check for multicollinearity
  cat("\n--- VIF (multicollinearity check) ---\n")
  tryCatch({
    library(car)
    v <- vif(fit)
    print(round(v, 3))
    high_vif <- names(v[v > 5])
    if (length(high_vif) > 0)
      cat(sprintf("WARNING: VIF > 5 for: %s\n", paste(high_vif, collapse = ", ")))
    else
      cat("All VIF < 5 — no multicollinearity concern\n")
  }, error = function(e) cat("VIF could not be computed:", e$message, "\n"))

  out
}

# Merge time-varying covariates at TOA for threshold
thr_sub <- adopters[!is.na(adopters$threshold), ]
thr_sub <- merge(thr_sub, timevar_df,
                 by.x = c("id_orig", "toa"), by.y = c("id", "grade_period"), all.x = TRUE)

thr_gad <- run_threshold_reg(thr_sub, "gad", "GAD")
thr_biv_gad <- run_threshold_bivariate(thr_sub, "gad")
if (!is.null(thr_gad)) {
  write.csv(thr_gad,
            file.path(data_path, "0512_netdiffuseR_past6mo_threshold_GAD.csv"),
            row.names = FALSE)
  cat("Saved: threshold_GAD.csv\n")
}
if (!is.null(thr_biv_gad)) {
  write.csv(thr_biv_gad,
            file.path(data_path, "0512_netdiffuseR_past6mo_threshold_bivariate_GAD.csv"),
            row.names = FALSE)
  cat("Saved: threshold_bivariate_GAD.csv\n")
}

thr_mdd <- run_threshold_reg(thr_sub, "mdd", "MDD")
thr_biv_mdd <- run_threshold_bivariate(thr_sub, "mdd")
if (!is.null(thr_mdd)) {
  write.csv(thr_mdd,
            file.path(data_path, "0512_netdiffuseR_past6mo_threshold_MDD.csv"),
            row.names = FALSE)
  cat("Saved: threshold_MDD.csv\n")
}
if (!is.null(thr_biv_mdd)) {
  write.csv(thr_biv_mdd,
            file.path(data_path, "0512_netdiffuseR_past6mo_threshold_bivariate_MDD.csv"),
            row.names = FALSE)
  cat("Saved: threshold_bivariate_MDD.csv\n")
}

write.csv(adopters,
          file.path(data_path, "0511_netdiffuseR_past6mo_node_outcomes_adopters.csv"),
          row.names = FALSE)

# Exposure regression
cat("\nCalculating network exposure...\n")
expo      <- exposure(diffnet_all)
adopt_mat <- diffnet_all$cumadopt
n_period  <- ncol(adopt_mat)
vertex_ids_clean <- sub("^[^.]+\\.", "", names(toa_vals))

reg_rows <- list()
for (gi in 1:n_period) {
  not_yet <- which(adopt_mat[, gi] == 0)
  adopt_next_vals <- if (gi < n_period) as.integer(adopt_mat[not_yet, gi + 1]) else rep(0L, length(not_yet))
  reg_rows[[gi]] <- data.frame(
    id = vertex_ids_clean[not_yet], adopt_next = adopt_next_vals,
    exposure = as.numeric(expo[not_yet, gi]), grade_period = all_gps[gi],
    stringsAsFactors = FALSE
  )
}
reg_data <- do.call(rbind, reg_rows)
reg_data <- merge(reg_data, cov_df, by = "id", all.x = TRUE)
reg_data <- merge(reg_data, attrs_df[, c("id","cohort","sch_type")], by = "id", all.x = TRUE)
reg_data <- merge(reg_data, timevar_df,
                  by.x = c("id", "grade_period"), by.y = c("id", "grade_period"), all.x = TRUE)

cat(sprintf("gad NAs: %d | mdd NAs: %d | friends_ecig NAs: %d out of %d rows\n",
            sum(is.na(reg_data$gad)), sum(is.na(reg_data$mdd)),
            sum(is.na(reg_data$friends_ecig)), nrow(reg_data)))

# Exposure regression
# Includes: female, par_edu, factor(school) as fixed effect (per Dr. Valente)
run_exposure_reg <- function(reg_data, mh_var, label) {
  cat(sprintf("\n===== EXPOSURE REGRESSION (%s) =====\n", label))
  exp_covs <- c("cohort", "female", "hispanic", "asian", "par_edu",
                "factor(school)", mh_var, "friends_ecig")
  if (sum(!is.na(reg_data$sex_min)) > 100 && length(unique(na.omit(reg_data$sex_min))) > 1)
    exp_covs <- c(exp_covs, "sex_min")
  req_vars <- c("adopt_next", "exposure", "grade_period", "cohort",
                "female", "hispanic", "asian", "par_edu", "school", mh_var, "friends_ecig")
  reg_final <- reg_data[complete.cases(reg_data[, req_vars[req_vars %in% names(reg_data)]]), ]
  cat(sprintf("Person-period obs: %d\n", nrow(reg_final)))
  fml <- as.formula(paste("adopt_next ~ exposure + factor(grade_period) +",
                           paste(exp_covs, collapse = " + ")))
  cat("Formula:", deparse(fml), "\n")
  fit <- glm(fml, data = reg_final, family = binomial)
  ct  <- summary(fit)$coefficients
  rownames(ct) <- gsub("factor\\(grade_period\\)", "GP", rownames(ct))
  rownames(ct) <- gsub("factor\\(school\\)", "School", rownames(ct))
  out <- data.frame(
    parameter = rownames(ct),
    estimate  = round(ct[,1], 4),
    se        = round(ct[,2], 4),
    z         = round(ct[,3], 4),
    pval      = round(ct[,4], 6),
    OR        = round(exp(ct[,1]), 4),
    stringsAsFactors = FALSE
  )
  print(out)
  or_exp <- exp(coef(fit)["exposure"])
  ci_exp <- exp(confint(fit, "exposure"))
  cat(sprintf("Exposure OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
              or_exp, ci_exp[1], ci_exp[2], ct["exposure", 4]))

  # VIF check
  cat("\n--- VIF (multicollinearity check) ---\n")
  tryCatch({
    library(car)
    v <- vif(fit)
    print(round(v, 3))
    high_vif <- names(v[v > 5])
    if (length(high_vif) > 0)
      cat(sprintf("WARNING: VIF > 5 for: %s\n", paste(high_vif, collapse = ", ")))
    else
      cat("All VIF < 5 — no multicollinearity concern\n")
  }, error = function(e) cat("VIF could not be computed:", e$message, "\n"))

  out
}

exp_gad <- run_exposure_reg(reg_data, "gad", "GAD")
write.csv(exp_gad,
          file.path(data_path, "0512_netdiffuseR_past6mo_exposure_GAD.csv"),
          row.names = FALSE)
cat("Saved: exposure_GAD.csv\n")

exp_mdd <- run_exposure_reg(reg_data, "mdd", "MDD")
write.csv(exp_mdd,
          file.path(data_path, "0512_netdiffuseR_past6mo_exposure_MDD.csv"),
          row.names = FALSE)
cat("Saved: exposure_MDD.csv\n")

# School bivariate analyses (per Dr. Valente: bivariate first, then fixed effect)
cat("\n===== SCHOOL BIVARIATE ANALYSES =====\n")
outcomes <- list(
  exposure  = reg_data[, c("adopt_next","exposure","school","sch_type")],
  threshold = thr_sub[, c("threshold","school","sch_type")[
    c("threshold","school","sch_type") %in% names(thr_sub)]]
)
for (out_name in names(outcomes)) {
  df_out <- outcomes[[out_name]]
  y_var  <- names(df_out)[1]
  if (!"sch_type" %in% names(df_out)) next
  cat(sprintf("\n%s ~ school type:\n", y_var))
  fml <- as.formula(paste(y_var, "~ sch_type"))
  fit_biv <- tryCatch(lm(fml, data = df_out), error = function(e) NULL)
  if (!is.null(fit_biv)) print(summary(fit_biv)$coefficients)
}

# Classification
cat("\n===== ADOPTER CLASSIFICATION =====\n")
cls <- classify(diffnet_all, include_censored = TRUE)
print(ftable(cls))
write.csv(as.data.frame(ftable(cls)),
          file.path(data_path, "0512_netdiffuseR_past6mo_classification.csv"),
          row.names = FALSE)

write.csv(adopters,
          file.path(data_path, "0512_netdiffuseR_past6mo_node_outcomes_adopters.csv"),
          row.names = FALSE)

write.csv(results_table,
          file.path(data_path, "0512_netdiffuseR_past6mo_school_summary.csv"),
          row.names = FALSE)

cat("\nBuilding Excel workbook...\n")
wb <- createWorkbook()
read_if_exists <- function(path) {
  if (file.exists(path)) read.csv(path) else data.frame(note = "File not found")
}
sheets <- list(
  school_summary         = "0512_netdiffuseR_past6mo_school_summary.csv",
  exposure_GAD           = "0512_netdiffuseR_past6mo_exposure_GAD.csv",
  exposure_MDD           = "0512_netdiffuseR_past6mo_exposure_MDD.csv",
  threshold_GAD          = "0512_netdiffuseR_past6mo_threshold_GAD.csv",
  threshold_MDD          = "0512_netdiffuseR_past6mo_threshold_MDD.csv",
  threshold_biv_GAD      = "0512_netdiffuseR_past6mo_threshold_bivariate_GAD.csv",
  threshold_biv_MDD      = "0512_netdiffuseR_past6mo_threshold_bivariate_MDD.csv",
  classification         = "0512_netdiffuseR_past6mo_classification.csv",
  node_outcomes          = "0512_netdiffuseR_past6mo_node_outcomes_adopters.csv"
)
for (sheet_name in names(sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, read_if_exists(file.path(data_path, sheets[[sheet_name]])))
}
saveWorkbook(wb,
             file.path(data_path, "0512_netdiffuseR_past6mo_ALL_RESULTS.xlsx"),
             overwrite = TRUE)
cat("Saved: 0512_netdiffuseR_past6mo_ALL_RESULTS.xlsx\n")

pdf(file.path(data_path, "0512_netdiffuseR_past6mo_diffnet_allschools.pdf"),
    width = 14, height = 10)
for (sch_key in names(diffnet_list)) {
  dn_sch  <- diffnet_list[[sch_key]]
  n_nodes <- length(dn_sch$toa)
  sch_int <- as.integer(sch_key)
  label   <- sprintf("School %s (%s, %s, n=%d)",
                     sch_key, get_cohort(sch_int), get_schtype(sch_int), n_nodes)
  if (n_nodes > 500) {
    cat(sprintf("  School %s: n=%d too large, skipping\n", sch_key, n_nodes))
    next
  }
  plot_diffnet(dn_sch, slices = seq_along(all_gps),
               vertex.color = c("grey70", "tomato"),
               main = label)
  cat(sprintf("  Plotted school %s\n", sch_key))
}
dev.off()
cat("Saved: 0512_netdiffuseR_past6mo_diffnet_allschools.pdf\n")

save(diffnet_list, diffnet_all, results_table, node_df,
     file = file.path(data_path, "0512_netdiffuseR_past6mo.RData"))

cat("\n===== DONE =====\n")
cat("Outputs in:", data_path, "\n")
cat("Excel: 0512_netdiffuseR_past6mo_ALL_RESULTS.xlsx\n")
cat("  Sheets: school_summary / exposure_GAD / exposure_MDD\n")
cat("          threshold_GAD / threshold_MDD / threshold_biv_GAD / threshold_biv_MDD\n")
cat("          classification / node_outcomes\n")
cat("Network PDF: 0512_netdiffuseR_past6mo_diffnet_allschools.pdf\n")

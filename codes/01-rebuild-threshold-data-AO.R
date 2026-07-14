# 01-rebuild-threshold-data-AO.R
# Rebuild Yuchan's threshold dataframe (`thr_sub`) with full covariate set,
# including time-varying covariates (mdd, gad, friends_ecig) merged at TOA.
# Output: outputs_AO/intermediate/thr_data-AO.rds
#
# This script is a parametrized derivative of:
#   codes/260427_netdiffuseR_past6mo 4.17.35 PM.R
# Yuchan's script is left untouched. All differences are confined to:
#   - configurable data_path
#   - tolerant filename matching for w9/w10 quirks
#   - additional sch_type column carried through to thr_sub
#   - serialization of intermediate objects for downstream plotting scripts

suppressPackageStartupMessages({
  library(netdiffuseR)
})

# ---- Config -----------------------------------------------------------------
project_root <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/Z-Network-Diffusion-Yuchan"
data_path    <- "/Users/anibaloliveramorales/Desktop/Doctorado/-Projects-/A - Network-Disadoption/data/advance/Cleaned-Data"
out_dir      <- file.path(project_root, "outputs_AO", "intermediate")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

id_var       <- "record_id"
schoolid_var <- "schoolid"

c1_early_schools <- c(101, 102, 103, 104, 105)
c1_late_schools  <- c(106, 107, 112, 113, 114)
c2_schools       <- c(201, 212, 213, 214)
exclude_schools  <- c(108)
all_schools      <- setdiff(c(c1_early_schools, c1_late_schools, c2_schools), exclude_schools)
asian_schools    <- c(103, 105, 112, 113, 212, 213)
hispanic_schools <- c(102, 106, 107, 114, 214)

INCLUDE_GRADE9 <- FALSE
all_gps        <- 3:8
gp_labels      <- c("10Fa","10Sp","11Fa","11Sp","12Fa","12Sp")

school_waves <- function(sch) {
  if (sch %in% c1_early_schools) return(if (INCLUDE_GRADE9) 1:8 else 3:8)
  if (sch %in% c1_late_schools)  return(3:8)
  if (sch %in% c2_schools)       return(if (INCLUDE_GRADE9) 3:10 else 5:10)
  integer(0)
}
w_to_gp <- function(sch, w) {
  if (sch %in% c(c1_early_schools, c1_late_schools)) return(as.integer(w))
  if (sch %in% c2_schools) return(as.integer(w) - 2L)
  NA_integer_
}
get_cohort <- function(sch) if (sch %in% c(c1_early_schools, c1_late_schools)) "C1" else "C2"
get_schtype <- function(sch) {
  if (sch %in% asian_schools)    return("Asian-majority")
  if (sch %in% hispanic_schools) return("Hispanic-majority")
  "Other"
}

# Tolerant filename resolver for w9/w10 quirks
resolve_data_file <- function(w) {
  candidates <- c(
    sprintf("w%d_adv_data.csv", w),
    sprintf("w%dadv_data.csv", w),
    sprintf("w%dadv_data_clean.csv", w),
    sprintf("w%d_adv_data_clean.csv", w)
  )
  for (c in candidates) {
    p <- file.path(data_path, c)
    if (file.exists(p)) return(p)
  }
  NA_character_
}
resolve_edge_file <- function(w) {
  candidates <- c(
    sprintf("w%dedges_clean.csv", w),
    sprintf("w%d_edges_clean.csv", w)
  )
  for (c in candidates) {
    p <- file.path(data_path, c)
    if (file.exists(p)) return(p)
  }
  NA_character_
}

# ---- Column readers (verbatim from Yuchan, plus light helpers) --------------
get_ecig_past6mo <- function(df, w) {
  v <- paste0("w", w, "_past_6mo_use_3")
  if (v %in% names(df)) {
    x <- suppressWarnings(as.numeric(df[[v]]))
    return(as.integer(x == 1))
  }
  rep(NA_integer_, nrow(df))
}
get_field <- function(df, pattern) {
  hits <- grep(paste0("^", pattern, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) return(suppressWarnings(as.numeric(df[[hits[1]]])))
  rep(NA_real_, nrow(df))
}
get_gad          <- function(df, w) get_field(df, paste0("w", w, "_rcads_gad_mean"))
get_mdd          <- function(df, w) get_field(df, paste0("w", w, "_rcads_mdd_mean"))
get_friends_ecig <- function(df, w) {
  x <- get_field(df, paste0("w", w, "_friends_use_ecig"))
  x[x == 6] <- NA
  x
}
find_col <- function(df, colname) {
  hits <- grep(paste0("^", colname, "$"), names(df), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) return(NULL)
  df[[hits[1]]]
}
build_school_adj <- function(edgelist, raw_ids, vertex_ids) {
  edges <- edgelist[edgelist$ego %in% raw_ids & edgelist$alter %in% raw_ids, ]
  n     <- length(raw_ids)
  mat   <- matrix(0L, nrow = n, ncol = n, dimnames = list(vertex_ids, vertex_ids))
  if (nrow(edges) > 0) {
    for (k in seq_len(nrow(edges))) {
      i <- match(edges$ego[k],   raw_ids)
      j <- match(edges$alter[k], raw_ids)
      if (!is.na(i) && !is.na(j)) mat[i, j] <- 1L
    }
  }
  diag(mat) <- 0L
  mat
}

# ---- Load all waves ---------------------------------------------------------
needed_waves <- sort(unique(unlist(lapply(all_schools, school_waves))))
cat("Loading data files from:", data_path, "\n")
loaded_data <- list(); loaded_edge <- list()
for (w in needed_waves) {
  f_data <- resolve_data_file(w)
  f_edge <- resolve_edge_file(w)
  if (!is.na(f_data)) {
    loaded_data[[as.character(w)]] <- read.csv(f_data)
    cat(sprintf("  W%d data: %s (%d rows)\n", w, basename(f_data), nrow(loaded_data[[as.character(w)]])))
  } else cat(sprintf("  W%d data NOT FOUND\n", w))
  if (!is.na(f_edge)) {
    loaded_edge[[as.character(w)]] <- read.csv(f_edge)
    cat(sprintf("  W%d edges: %s (%d rows)\n", w, basename(f_edge), nrow(loaded_edge[[as.character(w)]])))
  } else cat(sprintf("  W%d edges NOT FOUND\n", w))
}

# ---- Build diffnet per school (verbatim logic, Yuchan's) --------------------
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

  beh_mat <- matrix(NA_integer_, nrow = n_students, ncol = length(all_gps),
                    dimnames = list(vertex_ids, as.character(all_gps)))
  for (wi in seq_along(sch_waves)) {
    w <- sch_waves[wi]; gp <- sch_gps[wi]; wk <- as.character(w)
    if (!wk %in% names(loaded_data)) next
    d   <- loaded_data[[wk]]
    idx <- match(raw_ids, as.character(d[[id_var]]))
    r   <- d[idx, ]
    beh_mat[, as.character(gp)] <- get_ecig_past6mo(r, w)
  }

  tv_rows <- list()
  for (wi in seq_along(sch_waves)) {
    w <- sch_waves[wi]; gp <- sch_gps[wi]; wk <- as.character(w)
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
    w <- sch_waves[wi]; gp <- sch_gps[wi]; wk <- as.character(w)
    if (!wk %in% names(loaded_edge)) next
    net_array[,, as.character(gp)] <- build_school_adj(loaded_edge[[wk]], raw_ids, vertex_ids)
  }

  dn <- tryCatch(
    as_diffnet(graph = net_array, toa = toa, t0 = min(all_gps), t1 = max(all_gps)),
    error = function(e) { cat(sprintf("School %d diffnet error: %s\n", sch, e$message)); NULL }
  )
  if (is.null(dn)) next
  diffnet_list[[as.character(sch)]] <- dn

  baseline_w  <- sch_waves[1]
  baseline_wk <- as.character(baseline_w)
  if (baseline_wk %in% names(loaded_data)) {
    d_base <- loaded_data[[baseline_wk]]
    idx_b  <- match(raw_ids, as.character(d_base[[id_var]]))
    gen_r  <- find_col(d_base, paste0("w", baseline_w, "_dem_gender"))[idx_b]
    eth_r  <- find_col(d_base, paste0("w", baseline_w, "_eth"))[idx_b]
    race_r <- find_col(d_base, paste0("w", baseline_w, "_race"))[idx_b]
    sex_r  <- find_col(d_base, paste0("w", baseline_w, "_dem_sexuality"))[idx_b]
    # Bug fix vs Yuchan's script: codebook says wN_DEM_GENDER is 0=Female, 1=Male, 3=Decline.
    # Yuchan compared against ==2 (which never appears), making `female` always 0 and
    # silently dropped by R from the regressions. Real coding: female = (gen_r == 0).
    female   <- as.integer(!is.na(gen_r)  & gen_r  == 0)
    hispanic <- as.integer(!is.na(eth_r)  & eth_r  == 1)
    asian    <- as.integer(!is.na(race_r) & suppressWarnings(as.numeric(race_r)) == 2)
    sex_raw2 <- suppressWarnings(as.numeric(sex_r))
    sex_min  <- as.integer(!is.na(sex_raw2) & sex_raw2 != 1)
    covariate_list[[as.character(sch)]] <- data.frame(
      id = vertex_ids, female = female, hispanic = hispanic,
      asian = asian, sex_min = sex_min, stringsAsFactors = FALSE
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
diffnet_all <- do.call(c, diffnet_list)

cov_df <- do.call(rbind, covariate_list)
attrs_df <- do.call(rbind, lapply(names(diffnet_list), function(k) {
  sch_int <- as.integer(k)
  data.frame(id = names(diffnet_list[[k]]$toa), school = sch_int,
             cohort = get_cohort(sch_int), sch_type = get_schtype(sch_int),
             stringsAsFactors = FALSE)
}))
attrs_df   <- merge(attrs_df, cov_df, by = "id", all.x = TRUE)
timevar_df <- do.call(rbind, timevar_list)

# Thresholds, then build thr_sub with time-varying covariates joined at TOA
thr_vals <- threshold(diffnet_all)
toa_vals <- diffnet_all$toa
node_df  <- data.frame(
  id        = names(toa_vals),
  id_orig   = sub("^[^.]+\\.", "", names(toa_vals)),
  toa       = as.integer(toa_vals),
  threshold = as.numeric(thr_vals),
  stringsAsFactors = FALSE
)
node_df  <- merge(node_df, attrs_df, by.x = "id_orig", by.y = "id", all.x = TRUE)
adopters <- node_df[!is.na(node_df$toa), ]

thr_sub <- adopters[!is.na(adopters$threshold), ]
thr_sub <- merge(thr_sub, timevar_df,
                 by.x = c("id_orig", "toa"), by.y = c("id", "grade_period"),
                 all.x = TRUE)

cat(sprintf("\nthr_sub: %d adopters, %d with mdd, %d with gad, %d with friends_ecig\n",
            nrow(thr_sub),
            sum(!is.na(thr_sub$mdd)),
            sum(!is.na(thr_sub$gad)),
            sum(!is.na(thr_sub$friends_ecig))))

# ---- Save -------------------------------------------------------------------
saveRDS(list(
  thr_sub      = thr_sub,
  diffnet_all  = diffnet_all,
  diffnet_list = diffnet_list,
  attrs_df     = attrs_df,
  timevar_df   = timevar_df,
  cov_df       = cov_df,
  results_table = results_table,
  gp_labels    = gp_labels,
  all_gps      = all_gps
), file = file.path(out_dir, "thr_data-AO.rds"))

cat("\nSaved: ", file.path(out_dir, "thr_data-AO.rds"), "\n", sep = "")
cat("Done.\n")

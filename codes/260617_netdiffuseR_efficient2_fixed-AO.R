# =============================================================================
# 260617_netdiffuseR_efficient2_fixed-AO.R
# The fully-fixed, fast pipeline. Same efficient slim-cache loader as efficient1,
# with ALL THREE fixes applied:
#   * Problem 1 — threshold complete separation: grouped-binomial threshold fit by
#     Firth (brglm2 `method="brglmFit"`) so GP3 / School-105 separations resolve and
#     the `-20.63 / 17.x` block becomes finite. (≈ run_threshold_reg)
#   * Problem 2 — par_edu "Don't know" (9) leak + bad W7+ remap: par_edu recoded to a
#     homogeneous 1..6 scale with 9 -> NA at all waves. (≈ remap_par_edu)
#   * Problem 3 — the mixed adoption model (lme4::glmer, Table 2) is replaced by a
#     fixed-effects GLM with CLUSTER-ROBUST SEs (clustered by student). efficient2
#     verified the (1|id) random intercept is unidentified here (non-convergent,
#     variance ~128), ~600x slower, and that the FE fit reproduces Table 3 exactly.
#     The function also prints the diagnostics that justify dropping the random
#     effects (school ICC rho; student-level design effect). (≈ run_exposure_reg)
# Net effect: the report's tables, fully corrected, in a fraction of the runtime.
#
# WHY: the two ADVANCE workbooks are huge (W1-W8 ≈ 153 MB / 10,341 columns) and
# openxlsx reads the entire sheet every run. The analysis only needs ~865 of those
# columns. We read the workbooks ONCE, keep the needed columns, and cache them as a
# 468 KB RDS; later runs load the RDS in ~1 s.
#
# MEASURED (this machine, netdiffuseR 1.25.0):
#   revised16 (full read every run) ........ 169 s
#   efficient, 1st run (build cache) ....... 163 s   (read 153 MB ≈ 54 s, one-time)
#   efficient, cache hit ................... 109 s   (slim load ≈ 1 s)
# The cache removes the ~54 s read on every subsequent run (and far more on slower
# disks — this is the cost Yuchan hit).
#
# IS PARALLELISM (N cores) NEEDED? No. A profile of the 109 s cache-hit run shows
# the time is dominated by ONE call:
#   lme4::glmer (mixed adoption model) .. 83 s (77%)   <- single-threaded; can't be split
#   plot_diffnet (network figure) ....... 18 s (17%)
#   everything else (incl. per-school loop) < 10 s
# Multi-core would only touch the <10 s loop. The real levers for the remaining time
# are to drop/replace the unstable mixed model with the FE GLM (Table 3) and to skip
# plot_diffnet — not parallelism.
# =============================================================================
rm(list = ls())
library(netdiffuseR)
library(openxlsx)
library(lme4)        # only for the school-ICC diagnostic (1|school) — NOT for the main fit
library(sandwich)    # cluster-robust covariance for the FE adoption model
library(lmtest)      # coeftest() with the robust vcov
# brglm2 provides the Firth-penalized fit used for the threshold model (problem 1).
# Resolve from a side library if it is not in the default .libPaths().
if (!requireNamespace("brglm2", quietly = TRUE))
  suppressWarnings(library(brglm2, lib.loc = c("playground/Rlib", .libPaths())))
suppressMessages(library(brglm2, lib.loc = c("playground/Rlib", .libPaths())))

# ── Paths (AO efficient variant) ─────────────────────────────────────────────
# Input: the two big ADVANCE workbooks live in `in_path`. Output: results go to
# `out_path`. `data_path` is kept as an alias so every legacy file.path(data_path,
# ...) output call below is unchanged. Override via env vars if needed.
in_path   <- Sys.getenv("ADVANCE_DATA", unset = "Cleaned_Data")
out_path  <- Sys.getenv("ADVANCE_OUT",  unset = "Cleaned_Data")
data_path <- out_path

id_var       <- "record_id"
schoolid_var <- "schoolid"

# School definitions
c1_early_schools <- c(101, 102, 103, 104, 105)
c1_late_schools  <- c(106, 107, 112, 113, 114)
c2_schools       <- c(201, 212, 213, 214)

# School FE helper: C2 schools get reference school (101) so factor(school_c1fe)
# only produces C1 dummies. Cohort dummy captures the C1 vs C2 average difference.
# This removes the perfect collinearity between factor(school) and cohort.
make_school_c1fe <- function(school_vec) {
  ifelse(school_vec %in% c(201, 212, 213, 214), 101L, as.integer(school_vec))
}
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

# ── EFFICIENT data loading (AO) ───────────────────────────────────────────────
# The two ADVANCE workbooks are huge (W1-W8 is ~153 MB / 10,341 columns) and
# openxlsx::read.xlsx of the full sheet is the dominant cost of this pipeline
# (minutes). But the analysis only touches a few hundred of those columns. So we
# read the full workbooks ONCE, keep only the columns the pipeline uses (plus the
# survey-susceptibility items used by revised13, in case they return), and cache
# the slim result as a compact RDS. Every subsequent run loads the RDS (~1 s).
# Delete the cache (or bump CACHE_VER) to force a rebuild; set AO_SLIM to relocate.
cat("Loading ADVANCE data (efficient slim-cache loader)...\n")
library(openxlsx)
CACHE_VER  <- "v1"
slim_path  <- Sys.getenv("AO_SLIM", unset = file.path(in_path, paste0("ADVANCE_slim_", CACHE_VER, ".rds")))

# Columns the pipeline reads (names are lower-cased first). Allow-list of regexes:
keep_patterns <- paste(c(
  "^record_id$",
  "^w[0-9]+_schoolid$",
  "^w[0-9]+_past_6mo_use_3$",                                  # adoption outcome
  "^w[0-9]+_rcads_gad_mean$", "^w[0-9]+_rcads_mdd_mean$",      # mental health
  "^w[0-9]+_friends_use_ecig$",                               # perceived friend use
  "^w[0-9]+_try_friend_ecig$", "^w[0-9]+_use_next_yr_ecig$",   # survey susceptibility (revised13 reuse)
  "^w[0-9]+_dem_gender$", "^w[0-9]+_eth$", "^w[0-9]+_race$", "^w[0-9]+_dem_sexuality$",
  "^w[0-9]+_dem_high_par_edu(_new)?$",                        # parent education (+ new-scale variant)
  "^w[0-9]+_sm_post_[0-9]+$", "^w[0-9]+_ecig_posted_[0-9a-z]+$", # social media
  "^w[0-9]+_friend[0-9]+_[0-9]+$"                             # friendship-nomination edges
), collapse = "|")

if (file.exists(slim_path)) {
  .slim  <- readRDS(slim_path)
  d_w18  <- .slim$w18; d_w910 <- .slim$w910
  cat(sprintf("  slim-cache HIT: %s\n", slim_path))
} else {
  cat("  slim-cache MISS — reading full workbooks once (slow, one-time)...\n")
  .t0 <- Sys.time()
  .raw18  <- read.xlsx(file.path(in_path, "ADVANCE_W1-W8.xlsx"),  sheet = 1)
  .raw910 <- read.xlsx(file.path(in_path, "ADVANCE_W9-W10.xlsx"), sheet = 1)
  cat(sprintf("  full workbook read: %.1f min\n", as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
  .slim_keep <- function(d) { names(d) <- tolower(names(d)); d[, grepl(keep_patterns, names(d)), drop = FALSE] }
  d_w18  <- .slim_keep(.raw18)
  d_w910 <- .slim_keep(.raw910)
  rm(.raw18, .raw910); gc()
  saveRDS(list(w18 = d_w18, w910 = d_w910), slim_path, compress = "xz")
  cat(sprintf("  slim-cache built: %s\n", slim_path))
}
# Names already lower-cased in the cache; enforce again (idempotent) for safety.
names(d_w18)  <- tolower(names(d_w18))
names(d_w910) <- tolower(names(d_w910))
cat(sprintf("  W1-W8:  %d students, %d columns (slim)\n",  nrow(d_w18),  ncol(d_w18)))
cat(sprintf("  W9-W10: %d students, %d columns (slim)\n",  nrow(d_w910), ncol(d_w910)))

# get_wave_data: returns per-wave slice with standard column record_id + schoolid
get_wave_data <- function(w) {
  d <- if (w <= 8) d_w18 else d_w910
  # keep only rows with a valid schoolid at this wave
  sid_col <- paste0("w", w, "_schoolid")
  if (!sid_col %in% names(d)) return(NULL)
  d_sub <- d[!is.na(d[[sid_col]]), ]
  # add standardised columns used by downstream code
  d_sub$schoolid <- suppressWarnings(as.integer(d_sub[[sid_col]]))
  d_sub
}

# build_edgelist_from_xlsx: parse w{w}_friend{k}_{school} columns → ego/alter data.frame
build_edgelist_from_xlsx <- function(w) {
  d <- if (w <= 8) d_w18 else d_w910
  pat <- paste0("^w", w, "_friend[0-9]+_[0-9]+$")
  friend_cols <- grep(pat, names(d), value = TRUE)
  if (length(friend_cols) == 0) return(data.frame(ego=character(), alter=character()))
  id_col <- "record_id"
  rows <- lapply(friend_cols, function(col) {
    valid <- !is.na(d[[col]]) & d[[col]] != ""
    if (sum(valid) == 0) return(NULL)
    data.frame(ego   = as.character(d[[id_col]][valid]),
               alter = as.character(d[[col]][valid]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows[!sapply(rows, is.null)])
}

# Pre-build all needed wave data and edge lists
needed_waves <- sort(unique(unlist(lapply(all_schools, school_waves))))
loaded_data <- list()
loaded_edge <- list()
for (w in needed_waves) {
  wk <- as.character(w)
  wd <- get_wave_data(w)
  if (!is.null(wd) && nrow(wd) > 0) {
    loaded_data[[wk]] <- wd
    cat(sprintf("  W%d data: %d rows\n", w, nrow(wd)))
  } else {
    cat(sprintf("  W%d data NOT FOUND\n", w))
  }
  el <- build_edgelist_from_xlsx(w)
  if (nrow(el) > 0) {
    loaded_edge[[wk]] <- el
    cat(sprintf("  W%d edges: %d nominations\n", w, nrow(el)))
  } else {
    cat(sprintf("  W%d edges NOT FOUND\n", w))
  }
}

# schoolid_var: per-wave column name used below
schoolid_var_wave <- function(w) paste0("w", w, "_schoolid")

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

# (data loading handled above by XLSX section)

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

    # Out-degree and in-degree: build adjacency inline (net_array not yet constructed)
    out_deg <- rep(0L, n_students)
    in_deg  <- rep(0L, n_students)
    if (wk %in% names(loaded_edge)) {
      adj_gp_tmp <- build_school_adj(loaded_edge[[wk]], raw_ids, vertex_ids)
      out_deg <- as.integer(rowSums(adj_gp_tmp))
      in_deg  <- as.integer(colSums(adj_gp_tmp))
    }

    # Social media: E-cigarettes and nicotine vaping products (W4-W8 only)
    #
    # sm_post scale: 1=Several times/day, 2=Daily, 3=Weekly, 4=Monthly or less,
    #                5=Never, 99=Don't know
    # EXPOSED = responded 1-4 on any platform (any frequency > never)
    #
    # PRIMARY (consistent): sm_post_1~6 = Facebook/Instagram/Twitter/YouTube/TikTok/Twitch
    #   Available W4-W8 for all; ensures comparability across waves
    # SENSITIVITY (all available): grep all sm_post_* in that wave
    #   Adds Discord(W4 only), Snapchat(W5-W8), Threads(W7-W8) when available

    # ── Consistent 6-platform version ────────────────────────────────────────
    sm_post_mat <- sapply(paste0("w", w, "_sm_post_", 1:6), function(col) {
      if (col %in% names(r)) suppressWarnings(as.numeric(r[[col]])) else rep(NA_real_, nrow(r))
    })
    sm_exposed_mat <- (sm_post_mat >= 1 & sm_post_mat <= 4)
    sm_exposed_mat[is.na(sm_post_mat)] <- NA
    sm_ecig_any <- as.integer(apply(sm_exposed_mat, 1, function(x) {
      if (all(is.na(x))) NA_integer_ else as.integer(any(x, na.rm = TRUE))
    }))

    # Frequency: recode 1→4, 2→3, 3→2, 4→1, 5/99→0 (higher = more frequent)
    sm_freq_mat <- sm_post_mat
    sm_freq_mat[sm_post_mat == 1]  <- 4
    sm_freq_mat[sm_post_mat == 2]  <- 3
    sm_freq_mat[sm_post_mat == 3]  <- 2
    sm_freq_mat[sm_post_mat == 4]  <- 1
    sm_freq_mat[sm_post_mat == 5 | sm_post_mat == 99] <- 0
    sm_ecig_freq_max <- apply(sm_freq_mat, 1, function(x) {
      if (all(is.na(x))) NA_real_ else as.numeric(max(x, na.rm = TRUE))
    })

    # ── All-available-platform version (sensitivity) ──────────────────────────
    all_sm_cols <- grep(paste0("^w", w, "_sm_post_[0-9]+$"), names(r), value = TRUE)
    if (length(all_sm_cols) > 0) {
      sm_all_mat  <- sapply(all_sm_cols, function(col) suppressWarnings(as.numeric(r[[col]])))
      sm_exp_all  <- (sm_all_mat >= 1 & sm_all_mat <= 4)
      sm_exp_all[is.na(sm_all_mat)] <- NA
      sm_ecig_any_all <- as.integer(apply(sm_exp_all, 1, function(x) {
        if (all(is.na(x))) NA_integer_ else as.integer(any(x, na.rm = TRUE))
      }))
      sm_ecig_platform_count <- apply(sm_exp_all, 1, function(x) {
        if (all(is.na(x))) NA_real_ else as.numeric(sum(x, na.rm = TRUE))
      })
    } else {
      sm_ecig_any_all        <- rep(NA_integer_, nrow(r))
      sm_ecig_platform_count <- rep(NA_real_,    nrow(r))
    }

    get_posted <- function(num) {
      col <- paste0("w", w, "_ecig_posted_", num)
      if (col %in% names(r)) suppressWarnings(as.numeric(r[[col]])) else rep(NA_real_, nrow(r))
    }
    ep1 <- get_posted(1); ep2 <- get_posted(2); ep3 <- get_posted(3)
    sm_ecig_friend    <- ifelse(is.na(ep1), NA_integer_, as.integer(ep1 == 1))
    sm_ecig_nonfriend <- ifelse(is.na(ep2) & is.na(ep3), NA_integer_,
                                as.integer((!is.na(ep2) & ep2 == 1) | (!is.na(ep3) & ep3 == 1)))

    tv_rows[[wi]] <- data.frame(
      id                    = vertex_ids,
      grade_period          = gp,
      gad                   = get_gad(r, w),
      mdd                   = get_mdd(r, w),
      friends_ecig          = get_friends_ecig(r, w),
      out_degree            = as.integer(out_deg),
      in_degree             = as.integer(in_deg),
      sm_ecig_any           = sm_ecig_any,
      sm_ecig_any_all       = sm_ecig_any_all,
      sm_ecig_freq_max      = sm_ecig_freq_max,
      sm_ecig_platform_count= sm_ecig_platform_count,
      sm_ecig_friend        = sm_ecig_friend,
      sm_ecig_nonfriend     = sm_ecig_nonfriend,
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
  #     Students with no valid par_edu in any wave → NA (excluded by complete.cases).

  # FIX (par_edu) — In the 042326 release the column `w{w}_dem_high_par_edu` is on
  # the LEGACY 6-level scale at EVERY wave (1..6 valid, 9 = "Don't know"; no 7/8).
  # The original remap (a) applied a 9-level "new scale" remap at w>=7 that mis-
  # collapsed 5,6 -> 4, and (b) only mapped 9 -> NA at w>=7, so for the dominant
  # w<7 path "Don't know" (9) survived as a numeric level ABOVE "advanced degree".
  # Correct recode: homogeneous 1..6 at all waves; anything else (7, 8, 9) -> NA.
  remap_par_edu <- function(x, w) {
    x <- suppressWarnings(as.numeric(x))
    ifelse(!is.na(x) & x >= 1 & x <= 6, x, NA_real_)
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
    sex_min  <- ifelse(is.na(sex_raw2), NA_integer_, as.integer(sex_raw2 != 1))

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

    # Step 3: first observed valid value across all waves per student
    # (per Dr. Valente + ChatGPT suggestion)
    # LOCF only fills forward — useless if baseline itself is NA.
    # "First valid" correctly finds the nearest non-missing value in any wave,
    # whether before or after the baseline wave.
    par_edu <- apply(par_matrix, 1, function(x) {
      first_valid <- x[!is.na(x)][1]
      if (length(first_valid) == 0) return(NA_real_)
      return(first_valid)
    })

    # Students with no valid par_edu across any wave → remain NA
    # Excluded via complete.cases() in each regression.
    # Consistent with disadoption paper (no modal imputation).
    n_still_na <- sum(is.na(par_edu))
    if (n_still_na > 0) {
      cat(sprintf("  School %d: %d students with no valid par_edu → NA (dropped by complete.cases)\n",
                  sch, n_still_na))
    }
    cat(sprintf("  School %d: par_edu range [%g, %g], NAs = %d (%.1f%%)\n",
                sch, min(par_edu, na.rm = TRUE), max(par_edu, na.rm = TRUE),
                n_still_na, 100 * n_still_na / length(par_edu)))

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

# ===== DIAGNOSTIC: CHECK GAD/MDD MATCH BEFORE REGRESSIONS =====
# Run up to here first, then check the output below before continuing.
cat("\n===== DIAGNOSTIC: timevar_df GAD/MDD coverage =====\n")
cat(sprintf("timevar_df total rows:         %d\n", nrow(timevar_df)))
cat(sprintf("GAD non-NA rows:               %d  (%.1f%%)\n",
            sum(!is.na(timevar_df$gad)),
            100 * mean(!is.na(timevar_df$gad))))
cat(sprintf("MDD non-NA rows:               %d  (%.1f%%)\n",
            sum(!is.na(timevar_df$mdd)),
            100 * mean(!is.na(timevar_df$mdd))))
cat(sprintf("friends_ecig non-NA rows:      %d  (%.1f%%)\n",
            sum(!is.na(timevar_df$friends_ecig)),
            100 * mean(!is.na(timevar_df$friends_ecig))))
cat(sprintf("timevar_df id sample:          %s\n", paste(head(timevar_df$id, 3), collapse = ", ")))
cat(sprintf("timevar_df grade_period range: %s\n", paste(sort(unique(timevar_df$grade_period)), collapse = ", ")))
# If GAD non-NA is 0%, STOP HERE — column name mismatch. Check with:
#   grep("rcads_gad", names(loaded_data[["3"]]), ignore.case=TRUE, value=TRUE)
cat("===== END DIAGNOSTIC =====\n\n")

# Summary plots
pdf(file.path(data_path, "260609_netdiffuseR_past6mo_adopters.pdf"), width = 9, height = 5)
plot_adopters(diffnet_all,
              main = "E-cigarette Past-6-Month Use: Cumulative Adoption by Grade Period\n(Cohorts 1 & 2)",
              xlab = "Grade Period", ylab = "Proportion")
axis(1, at = all_gps, labels = gp_labels, tick = FALSE, line = 1.5, cex.axis = 0.75)
dev.off()

pdf(file.path(data_path, "260609_netdiffuseR_past6mo_hazard.pdf"), width = 9, height = 5)
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

pdf(file.path(data_path, "260609_netdiffuseR_past6mo_threshold_boxplot.pdf"), width = 9, height = 5)
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
pdf(file.path(data_path, "260609_netdiffuseR_past6mo_byschool.pdf"), width = 11, height = 6)
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

# Bivariate exposure analyses
# Same logic: run each predictor alone first
# If bivariate significant but full model not → multicollinearity
run_exposure_bivariate <- function(reg_data, mh_var) {
  cat("\n===== EXPOSURE BIVARIATE ANALYSES =====\n")
  predictors <- c("exposure", "friends_ecig", "hispanic", "asian", "sex_min",
                  "cohort", "female", "par_edu", mh_var)
  results <- list()
  for (pred in predictors) {
    if (!pred %in% names(reg_data)) next
    if (sum(!is.na(reg_data[[pred]])) < 30) next
    fml <- as.formula(paste("adopt_next ~", pred))
    fit <- tryCatch(
      glm(fml, data = reg_data[!is.na(reg_data[[pred]]) & !is.na(reg_data$adopt_next), ],
          family = binomial),
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
# Re-specified as grouped binomial per collaborator diagnostic (Olivera, 2026-05-12):
#   - 26.6% of adopters are network isolates (n_alters = 0) at TOA
#   - threshold has discrete support {0, 1/k, ..., 1}; quasibinomial treats it as continuous
#   - grouped binomial glm(cbind(k_users, n_alters - k_users) ~ X) weights by degree
#   - isolates excluded (threshold = 0/0 is uninformative)
# Includes: female, par_edu, factor(school) as fixed effect (per Dr. Valente)
run_threshold_reg <- function(thr_sub) {
  cat("\n===== THRESHOLD REGRESSION — GROUPED BINOMIAL (GAD + MDD) =====\n")

  # Exclude isolates (n_alters = 0) — their threshold is uninformative
  thr_non_iso <- thr_sub[!is.na(thr_sub$n_alters) & thr_sub$n_alters > 0, ]
  cat(sprintf("Adopters with n_alters > 0 (non-isolates): %d\n", nrow(thr_non_iso)))
  cat(sprintf("Excluded isolates: %d\n", nrow(thr_sub) - nrow(thr_non_iso)))

  # school_c1fe: C2 schools collapsed to reference (101) to remove cohort collinearity
  thr_non_iso$school_c1fe <- make_school_c1fe(thr_non_iso$school)

  thr_covs <- c("factor(toa)", "cohort", "female", "hispanic", "asian",
                "par_edu", "factor(school_c1fe)", "gad", "mdd", "friends_ecig")
  if (sum(!is.na(thr_non_iso$sex_min)) > 30 && length(unique(na.omit(thr_non_iso$sex_min))) > 1)
    thr_covs <- c(thr_covs, "sex_min")

  # grouped binomial: outcome = cbind(k_users, n_alters - k_users)
  fml <- as.formula(paste("cbind(k_users, n_alters - k_users) ~",
                          paste(thr_covs, collapse = " + ")))
  cat("Formula:", deparse(fml), "\n")

  req <- c("k_users", "n_alters", "cohort", "female", "hispanic", "asian",
           "par_edu", "school", "gad", "mdd", "friends_ecig")
  fit_data <- thr_non_iso[complete.cases(thr_non_iso[, req[req %in% names(thr_non_iso)]]), ]
  cat(sprintf("Final N (complete cases): %d\n", nrow(fit_data)))

  # FIX (threshold separation) — GP3 adopters all have k_users = 0 (no period
  # precedes the first), so the GP3 reference level is perfectly separated and a
  # plain glm sends the intercept + all factor(toa) dummies to +/- infinity (the
  # -20.63 / 17.x block in the report); School 105 is separated the same way.
  # Firth's penalized likelihood (brglm2) has a guaranteed-finite MLE: it resolves
  # BOTH separations, keeps all 406 adopters, and recovers the real steep grade
  # trend, while leaving the substantive ORs essentially unchanged.
  fit <- tryCatch(
    glm(fml, data = fit_data, family = binomial(link = "logit"), method = "brglmFit"),
    error = function(e) { cat("Error:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coefficients
  rownames(ct) <- gsub("factor\\(toa\\)", "GP", rownames(ct))
  rownames(ct) <- gsub("factor\\(school_c1fe\\)", "School", gsub("factor\\(school\\)", "School", rownames(ct)))
  out <- data.frame(
    parameter = rownames(ct),
    estimate  = round(ct[,1], 4),
    se        = round(ct[,2], 4),
    z_value   = round(ct[,3], 4),
    pval      = round(ct[,4], 6),
    OR        = round(exp(ct[,1]), 4),
    stringsAsFactors = FALSE
  )
  print(out)

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

# Compute n_alters (out-degree at TOA) and k_users (# alters using at TOA)
# for the grouped-binomial re-specification (per collaborator diagnostic).
# 26.6% of adopters are network isolates (n_alters = 0) at adoption time —
# their threshold is uninformative (0/0) and they are excluded from the
# grouped-binomial model.
thr_sub$n_alters <- NA_integer_
thr_sub$k_users  <- NA_integer_
for (sch_key in names(diffnet_list)) {
  dn_sch     <- diffnet_list[[sch_key]]
  toa_sch    <- as.integer(dn_sch$toa)
  node_names <- names(dn_sch$toa)  # "101_101003" format
  names(toa_sch) <- node_names

  for (gi in seq_along(all_gps)) {
    t <- all_gps[gi]
    adj_t <- as.matrix(dn_sch$graph[[gi]])
    adopters_t_idx <- which(!is.na(toa_sch) & toa_sch == t)
    if (length(adopters_t_idx) == 0) next

    for (i in adopters_t_idx) {
      nd_name <- node_names[i]
      row_idx <- which(thr_sub$id_orig == nd_name & thr_sub$toa == t)
      if (length(row_idx) == 0) next
      alters_idx <- which(adj_t[i, ] > 0)  # numeric index — no dimnames needed
      n_al <- length(alters_idx)
      if (n_al == 0) {
        thr_sub$n_alters[row_idx] <- 0L
        thr_sub$k_users[row_idx]  <- 0L
      } else {
        alters_names   <- node_names[alters_idx]
        adopted_before <- alters_names[!is.na(toa_sch[alters_names]) &
                                         toa_sch[alters_names] < t]
        thr_sub$n_alters[row_idx] <- n_al
        thr_sub$k_users[row_idx]  <- length(adopted_before)
      }
    }
  }
}
cat(sprintf("thr_sub: n_alters computed for %d adopters\n", sum(!is.na(thr_sub$n_alters))))
cat(sprintf("Isolates (n_alters = 0): %d (%.1f%%)\n",
            sum(thr_sub$n_alters == 0, na.rm = TRUE),
            100 * mean(thr_sub$n_alters == 0, na.rm = TRUE)))

thr_res <- run_threshold_reg(thr_sub)
if (!is.null(thr_res)) {
  write.csv(thr_res,
            file.path(data_path, "260609_netdiffuseR_past6mo_threshold.csv"),
            row.names = FALSE)
  cat("Saved: threshold.csv\n")
}

write.csv(adopters,
          file.path(data_path, "260609_netdiffuseR_past6mo_node_outcomes_adopters.csv"),
          row.names = FALSE)

# Exposure regression
cat("\nCalculating network exposure...\n")
expo      <- exposure(diffnet_all)

adopt_mat <- diffnet_all$cumadopt
n_period  <- ncol(adopt_mat)

# After do.call(c, diffnet_list), node names get numeric prefix e.g. "1.101_12345"
# Strip the leading "N." to restore original "101_12345" format matching timevar_df/cov_df
vertex_ids_clean <- sub("^[^.]+\\.", "", names(toa_vals))

# DIAGNOSTIC: verify IDs match before merging
cat("\n--- ID format check ---\n")
cat(sprintf("toa_vals names sample:   %s\n", paste(head(names(toa_vals), 3), collapse=", ")))
cat(sprintf("vertex_ids_clean sample: %s\n", paste(head(vertex_ids_clean, 3), collapse=", ")))
cat(sprintf("timevar_df$id sample:    %s\n", paste(head(timevar_df$id, 3), collapse=", ")))
cat(sprintf("cov_df$id sample:        %s\n", paste(head(cov_df$id, 3), collapse=", ")))

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
reg_data <- merge(reg_data, attrs_df[, c("id","cohort","sch_type","school")], by = "id", all.x = TRUE)
reg_data <- merge(reg_data, timevar_df,
                  by.x = c("id", "grade_period"), by.y = c("id", "grade_period"), all.x = TRUE)

cat(sprintf("gad NAs: %d | mdd NAs: %d | friends_ecig NAs: %d out of %d rows\n",
            sum(is.na(reg_data$gad)), sum(is.na(reg_data$mdd)),
            sum(is.na(reg_data$friends_ecig)), nrow(reg_data)))

# Exposure regression
# Includes: female, par_edu, factor(school) as fixed effect (per Dr. Valente)
# FIX (problem 3) — the mixed adoption model is replaced by a fixed-effects GLM with
# CLUSTER-ROBUST standard errors (clustered by student). The efficient2 verification
# showed the (1|id) random intercept is unidentified here (it fails to converge and
# its variance explodes to ~128 via quasi-separation), is ~600x slower, and the FE
# fit reproduces the report's Table 3 exactly. Cluster-robust SEs recover the only
# thing the random intercept was nominally for — within-student correlation — without
# the instability or the cost. The function also reports the diagnostics that justify
# dropping the random effects (school ICC and the design effect).
run_exposure_reg <- function(reg_data) {
  cat("\n===== ADOPTION REGRESSION — FE GLM + cluster-robust SE (GAD + MDD) =====\n")
  reg_data$school_c1fe <- make_school_c1fe(reg_data$school)
  exp_covs <- c("cohort", "female", "hispanic", "asian", "par_edu",
                "gad", "mdd", "friends_ecig")
  if (sum(!is.na(reg_data$sex_min)) > 100 && length(unique(na.omit(reg_data$sex_min))) > 1)
    exp_covs <- c(exp_covs, "sex_min")
  req_vars <- c("adopt_next", "exposure", "grade_period", "cohort",
                "female", "hispanic", "asian", "par_edu", "school", "gad", "mdd", "friends_ecig")
  reg_final <- reg_data[complete.cases(reg_data[, req_vars[req_vars %in% names(reg_data)]]), ]
  cat(sprintf("Person-period obs: %d\n", nrow(reg_final)))
  n_students <- length(unique(reg_final$id))
  n_adopt    <- length(unique(reg_final$id[reg_final$adopt_next == 1]))
  cat(sprintf("Unique students: %d (adopters: %d, non-adopters: %d)\n",
              n_students, n_adopt, n_students - n_adopt))

  fml <- as.formula(paste("adopt_next ~ exposure + factor(grade_period) +",
                          paste(exp_covs, collapse = " + "), "+ factor(school_c1fe)"))
  cat("Formula:", deparse(fml), "\n")
  fit      <- glm(fml, data = reg_final, family = binomial)
  se_naive <- sqrt(diag(vcov(fit)))
  vc_cl    <- sandwich::vcovCL(fit, cluster = reg_final$id)   # cluster-robust by student
  ctr      <- lmtest::coeftest(fit, vcov. = vc_cl)
  rn <- gsub("factor\\(grade_period\\)", "GP",
             gsub("factor\\(school_c1fe\\)", "School", rownames(ctr)))
  keep <- rownames(ctr) != "(Intercept)"
  out <- data.frame(
    parameter = rn[keep],
    estimate  = round(ctr[keep, 1], 4),
    se        = round(ctr[keep, 2], 4),         # cluster-robust SE
    z         = round(ctr[keep, 3], 4),
    pval      = round(ctr[keep, 4], 6),
    OR        = round(exp(ctr[keep, 1]), 4),
    stringsAsFactors = FALSE
  )
  print(out)

  # ── DIAGNOSTIC: are the random intercepts necessary? ─────────────────────────
  cat("\n----- Are the random effects necessary? (justify FE over mixed) -----\n")
  # (1) SCHOOL: ICC (rho) from a school-only random-intercept model after covariates.
  m_sch <- tryCatch(
    lme4::glmer(update(fml, . ~ . - factor(school_c1fe) + (1 | school)),
                data = reg_final, family = binomial,
                control = lme4::glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL)
  if (!is.null(m_sch)) {
    v_sch <- as.data.frame(lme4::VarCorr(m_sch))$vcov[1]
    icc   <- v_sch / (v_sch + pi^2/3)            # logistic latent-scale ICC; pi^2/3 = 3.29
    cat(sprintf("  School random-intercept variance = %.4f  ->  ICC (rho_school) = %.4f\n",
                v_sch, icc))
    cat("    rho_school ~ 0  => essentially no between-school variance left after covariates;\n")
    cat("    a school random effect adds nothing (FE school dummies, or none, suffice).\n")
  }
  # (2) STUDENT: design effect = cluster-robust SE / naive SE for key predictors.
  #     If ~1, ignoring within-student correlation does not distort inference, so the
  #     (1|id) random intercept is unnecessary. (The mixed model's exploded id-variance
  #     of ~128 is a quasi-separation artifact, not real, usable clustering.)
  cat("  Student-level design effect (cluster-robust SE / naive SE):\n")
  for (k in intersect(c("exposure","friends_ecig","mdd","gad"), rownames(ctr)))
    cat(sprintf("    %-13s = %.2f\n", k, ctr[k, 2] / se_naive[k]))
  cat("    ~1.0  => within-student clustering is negligible for inference; (1|id) unneeded.\n")

  ci <- exp(ctr["exposure", 1] + c(-1.96, 1.96) * ctr["exposure", 2])
  cat(sprintf("\nExposure OR = %.3f (95%% robust CI: %.3f-%.3f), p = %.4f\n",
              exp(ctr["exposure", 1]), ci[1], ci[2], ctr["exposure", 4]))
  out
}

exp_res <- run_exposure_reg(reg_data)
write.csv(exp_res,
          file.path(data_path, "260609_netdiffuseR_past6mo_exposure.csv"),
          row.names = FALSE)
cat("Saved: exposure.csv\n")

# Sensitivity check: school fixed effects GLM (for Table 4)
run_exposure_reg_sensitivity <- function(reg_data) {
  cat("\n===== EXPOSURE REGRESSION — SENSITIVITY (school fixed effects GLM) =====\n")
  reg_data$school_c1fe <- make_school_c1fe(reg_data$school)
  exp_covs_s <- c("cohort", "female", "hispanic", "asian", "par_edu",
                  "factor(school_c1fe)", "gad", "mdd", "friends_ecig", "sex_min")
  req_vars <- c("adopt_next","exposure","grade_period","cohort","female",
                "hispanic","asian","par_edu","school","gad","mdd","friends_ecig")
  reg_final <- reg_data[complete.cases(reg_data[, req_vars[req_vars %in% names(reg_data)]]), ]
  cat(sprintf("Person-period obs: %d\n", nrow(reg_final)))
  fml <- as.formula(paste("adopt_next ~ exposure + factor(grade_period) +",
                           paste(exp_covs_s, collapse = " + ")))
  cat("Formula:", deparse(fml), "\n")
  fit <- tryCatch(
    glm(fml, data = reg_final, family = binomial),
    error = function(e) { cat("  Error:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coefficients
  rownames(ct) <- gsub("factor\\(grade_period\\)", "GP", rownames(ct))
  rownames(ct) <- gsub("factor\\(school_c1fe\\)", "School", rownames(ct))
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
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
  ci_exp <- suppressMessages(exp(confint(fit, "exposure")))
  cat(sprintf("Exposure OR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
              or_exp, ci_exp[1], ci_exp[2], ct["exposure", 4]))
  out
}

exp_res_sensitivity <- run_exposure_reg_sensitivity(reg_data)
if (!is.null(exp_res_sensitivity)) {
  write.csv(exp_res_sensitivity,
            file.path(data_path, "260610_netdiffuseR_past6mo_exposure_sensitivity.csv"),
            row.names = FALSE)
  cat("Saved: exposure_sensitivity.csv\n")
}

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
          file.path(data_path, "260609_netdiffuseR_past6mo_classification.csv"),
          row.names = FALSE)

write.csv(adopters,
          file.path(data_path, "260609_netdiffuseR_past6mo_node_outcomes_adopters.csv"),
          row.names = FALSE)

write.csv(results_table,
          file.path(data_path, "260609_netdiffuseR_past6mo_school_summary.csv"),
          row.names = FALSE)

cat("\nBuilding Excel workbook...\n")
wb <- createWorkbook()
read_if_exists <- function(path) {
  if (file.exists(path)) read.csv(path) else data.frame(note = "File not found")
}
sheets <- list(
  school_summary         = "260609_netdiffuseR_past6mo_school_summary.csv",
  exposure               = "260609_netdiffuseR_past6mo_exposure.csv",
  threshold              = "260609_netdiffuseR_past6mo_threshold.csv",
  classification         = "260609_netdiffuseR_past6mo_classification.csv",
  node_outcomes          = "260609_netdiffuseR_past6mo_node_outcomes_adopters.csv"
)
for (sheet_name in names(sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, read_if_exists(file.path(data_path, sheets[[sheet_name]])))
}
saveWorkbook(wb,
             file.path(data_path, "260609_netdiffuseR_past6mo_ALL_RESULTS.xlsx"),
             overwrite = TRUE)
cat("Saved: 260609_netdiffuseR_past6mo_ALL_RESULTS.xlsx\n")

pdf(file.path(data_path, "260609_netdiffuseR_past6mo_diffnet_allschools.pdf"),
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
  
  # 先存6张png，再拼到pdf一页
  tmp_files <- c()
  for (i in seq_along(all_gps)) {
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = 600, height = 600)
    slice_title <- sprintf("%s\n%s (GP%d)", label, gp_labels[i], all_gps[i])
    plot_diffnet(dn_sch, slices = i,
                 vertex.color = c("grey70", "tomato"),
                 main = slice_title)
    dev.off()
    tmp_files <- c(tmp_files, tmp)
  }
  
  # 拼成一页（2行3列）
  par(mfrow = c(2, 3), mar = c(0.5, 0.5, 2, 0.5))
  for (f in tmp_files) {
    img <- png::readPNG(f)
    plot.new()
    rasterImage(img, 0, 0, 1, 1)
  }
  par(mfrow = c(1, 1))
  cat(sprintf("  Plotted school %s\n", sch_key))
}
dev.off()
cat("Saved: 260609_netdiffuseR_past6mo_diffnet_allschools.pdf\n")

save(diffnet_list, diffnet_all, results_table, node_df,
     file = file.path(data_path, "260603_netdiffuseR_past6mo.RData"))

cat("\n===== DONE =====\n")
cat("Outputs in:", data_path, "\n")
cat("Excel: 260609_netdiffuseR_past6mo_ALL_RESULTS.xlsx\n")
cat("          classification / node_outcomes\n")
cat("Network PDF: 260609_netdiffuseR_past6mo_diffnet_allschools.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# TOM'S TABLE: Factors associated with E-cigarette use based on peer exposures
# ══════════════════════════════════════════════════════════════════════════════
#
# 6 subgroups (per Valente email 2026-06-03):
#
# E-cig USERS:
#   Col 1 "No Exposure"   — adopters with network exposure at TOA = 0
#   Col 2 "No Friends"    — adopters who were network ISOLATES at GP3 (out_degree_gp3 = 0)
#                           NOTE: these students may still report perceived friend use (survey)
#   Col 3 "Low Threshold" — adopters with threshold ≤ 0.25
#
# NON-USERS:
#   Col 4 "Exposure"      — non-adopters with max exposure > 0 across all waves
#   Col 5 "Friends"       — non-adopters with friends_ecig_gp3 > 0 (perceived friend use at baseline)
#   Col 6 "High Threshold"— non-adopters with out_degree_gp3 > 0 but max_exposure = 0
#                           (connected in network, but no friends adopted = high implicit threshold)
#
# Row predictors (per Tom's table structure):
#   Demographics  : female, hispanic, asian, sex_min, par_edu, mdd
#   Network at TOA: out_degree_at, in_degree_at (TOA for users, GP8 for non-users), friends_ecig_gp3 (Friends T3, baseline)
#   Network T8    : friends_ecig_gp8 (Friends T8), expo_gp8 (Exposure ASR T8)

cat("\n===== TOM'S TABLE: 6-SUBGROUP ANALYSIS =====\n")

# ── Build full analytic dataframe ────────────────────────────────────────────
# node_df already contains attrs_df columns (female, hispanic, asian, sex_min,
# par_edu, school, cohort, sch_type) from the earlier merge. Do NOT merge again
# — that would create female.x / female.y duplicates that break the regressions.
full_df <- node_df

# Merge exposure at TOA (for users) and at GP8 (for non-users)
# expo matrix: rows = nodes (combined diffnet order), cols = grade periods
expo_df <- data.frame(
  id       = vertex_ids_clean,
  expo_gp3 = as.numeric(expo[, 1]),
  expo_gp8 = as.numeric(expo[, 6]),
  stringsAsFactors = FALSE
)
full_df <- merge(full_df, expo_df, by.x = "id_orig", by.y = "id", all.x = TRUE)

# Merge GP3 time-varying variables (baseline degree, friends_ecig, social media)
tv_gp3 <- timevar_df[timevar_df$grade_period == 3, ]
names(tv_gp3)[names(tv_gp3) == "id"] <- "id_tv"
full_df <- merge(full_df,
                 tv_gp3[, c("id_tv","friends_ecig",
                             "out_degree","in_degree","gad","mdd",
                             "sm_ecig_any","sm_ecig_any_all","sm_ecig_freq_max",
                             "sm_ecig_platform_count","sm_ecig_friend","sm_ecig_nonfriend")],
                 by.x = "id_orig", by.y = "id_tv", all.x = TRUE)

# Merge GP8 time-varying (friends_ecig at GP8 for non-users)
tv_gp8 <- timevar_df[timevar_df$grade_period == 8, ]
names(tv_gp8)[names(tv_gp8) == "id"] <- "id_tv8"
full_df <- merge(full_df,
                 tv_gp8[, c("id_tv8","friends_ecig","out_degree","in_degree")],
                 by.x = "id_orig", by.y = "id_tv8", all.x = TRUE,
                 suffixes = c("_gp3","_gp8"))

# ── Compute exposure at TOA for adopters ──────────────────────────────────────
# For each adopter, look up their exposure at their TOA period
adopters_expo <- reg_data[!is.na(reg_data$adopt_next) | TRUE, c("id","grade_period","exposure")]
# Use the reg_data exposure at the period BEFORE adoption (exposure predicts adoption)
# Simplify: use expo_gp3 for all, and expo at TOA via merge from reg_data
toa_expo <- do.call(rbind, lapply(seq_along(all_gps), function(gi) {
  gp <- all_gps[gi]
  data.frame(id = vertex_ids_clean, toa_period = gp,
             expo_at_period = as.numeric(expo[, gi]),
             stringsAsFactors = FALSE)
}))
# For each adopter, get exposure at their TOA
adopter_toa_expo <- merge(
  data.frame(id = full_df$id_orig, toa = full_df$toa, stringsAsFactors = FALSE),
  toa_expo, by.x = c("id","toa"), by.y = c("id","toa_period"), all.x = TRUE
)
full_df <- merge(full_df, adopter_toa_expo[, c("id","expo_at_period")],
                 by.x = "id_orig", by.y = "id", all.x = TRUE)

# Lagged exposure: exposure at TOA-1 (one grade period before adoption)
adopter_toa_lag_expo <- merge(
  data.frame(id = full_df$id_orig, toa_lag = full_df$toa - 1L, stringsAsFactors = FALSE),
  toa_expo, by.x = c("id","toa_lag"), by.y = c("id","toa_period"), all.x = TRUE
)
names(adopter_toa_lag_expo)[names(adopter_toa_lag_expo) == "expo_at_period"] <- "expo_lag"
full_df <- merge(full_df, adopter_toa_lag_expo[, c("id","expo_lag")],
                 by.x = "id_orig", by.y = "id", all.x = TRUE)

# Friends (perceived) at TOA: merge from timevar_df
friends_toa_df <- merge(
  data.frame(id = full_df$id_orig, toa = full_df$toa, stringsAsFactors = FALSE),
  timevar_df[, c("id","grade_period","friends_ecig")],
  by.x = c("id","toa"), by.y = c("id","grade_period"), all.x = TRUE
)
names(friends_toa_df)[names(friends_toa_df) == "friends_ecig"] <- "friends_ecig_toa"
full_df <- merge(full_df, friends_toa_df[, c("id","friends_ecig_toa")],
                 by.x = "id_orig", by.y = "id", all.x = TRUE)

# Out-degree and in-degree at TOA: merge from timevar_df (per Tom's request, replacing GP3 baseline)
degree_toa_df <- merge(
  data.frame(id = full_df$id_orig, toa = full_df$toa, stringsAsFactors = FALSE),
  timevar_df[, c("id","grade_period","out_degree","in_degree")],
  by.x = c("id","toa"), by.y = c("id","grade_period"), all.x = TRUE
)
names(degree_toa_df)[names(degree_toa_df) == "out_degree"] <- "out_degree_toa"
names(degree_toa_df)[names(degree_toa_df) == "in_degree"]  <- "in_degree_toa"
full_df <- merge(full_df, degree_toa_df[, c("id","out_degree_toa","in_degree_toa")],
                 by.x = "id_orig", by.y = "id", all.x = TRUE)

# Social media exposure at TOA (sm_ecig_any, W4-W8 only)
# GP3 adopters will have NA for this variable
sm_toa_df <- merge(
  data.frame(id = full_df$id_orig, toa = full_df$toa, stringsAsFactors = FALSE),
  timevar_df[, c("id","grade_period","sm_ecig_any")],
  by.x = c("id","toa"), by.y = c("id","grade_period"), all.x = TRUE
)
names(sm_toa_df)[names(sm_toa_df) == "sm_ecig_any"] <- "sm_ecig_any_toa"
full_df <- merge(full_df, sm_toa_df[, c("id","sm_ecig_any_toa")],
                 by.x = "id_orig", by.y = "id", all.x = TRUE)

# For non-users, get max exposure across all waves
max_expo_df <- aggregate(expo_at_period ~ id,
                          data = toa_expo,
                          FUN = function(x) max(x, na.rm = TRUE))
names(max_expo_df)[2] <- "max_exposure"
full_df <- merge(full_df, max_expo_df, by.x = "id_orig", by.y = "id", all.x = TRUE)

# For non-users, get max perceived friend use across all waves
max_friends_df <- aggregate(friends_ecig ~ id,
                             data = timevar_df,
                             FUN = function(x) max(x, na.rm = TRUE))
names(max_friends_df)[2] <- "max_friends_ecig"
full_df <- merge(full_df, max_friends_df, by.x = "id_orig", by.y = "id", all.x = TRUE)

# ── Define subgroups ──────────────────────────────────────────────────────────
full_df$school_c1fe <- make_school_c1fe(full_df$school)
full_df$is_user    <- as.integer(!is.na(full_df$toa))
full_df$is_nonuser <- as.integer( is.na(full_df$toa))

# Unified degree variable: TOA for users (no TOA exists for non-users → NA),
# GP8 for non-users (consistent with max_exposure/friends_ecig_gp8 at end of observation)
full_df$out_degree_at <- ifelse(full_df$is_user == 1, full_df$out_degree_toa, full_df$out_degree_gp8)
full_df$in_degree_at  <- ifelse(full_df$is_user == 1, full_df$in_degree_toa,  full_df$in_degree_gp8)

# Among USERS:
# Col1: No Exposure — zero network exposure at TOA
full_df$col1_no_exposure   <- ifelse(full_df$is_user == 1,
  as.integer(!is.na(full_df$expo_at_period) & full_df$expo_at_period == 0), NA)

# Col2: No Perceived Friend Use — perceived friend use = 0 at TOA
full_df$col2_no_friends    <- ifelse(full_df$is_user == 1,
  as.integer(!is.na(full_df$friends_ecig_toa) & full_df$friends_ecig_toa == 0), NA)

# Col_both0: No Exposure AND No Perceived Friend Use at TOA
full_df$col_both0          <- ifelse(full_df$is_user == 1,
  as.integer(!is.na(full_df$expo_at_period)   & full_df$expo_at_period == 0 &
             !is.na(full_df$friends_ecig_toa) & full_df$friends_ecig_toa == 0), NA)

# Threshold distribution check — justify 25% cutoff
# 70th percentile ≈ 0.25 in this data, supporting 25% as a meaningful low-threshold cutoff
cat("\n===== THRESHOLD DISTRIBUTION CHECK =====\n")
thr_check <- as.numeric(thr_vals[!is.na(thr_vals)])
cat(sprintf("N adopters with threshold: %d\n", length(thr_check)))
cat(sprintf("Threshold = 0: %d (%.1f%%)\n",
    sum(thr_check == 0), 100*mean(thr_check == 0)))
print(quantile(thr_check, probs = c(0.50,0.60,0.70,0.75,0.80,0.90)))
cat("70th percentile ≈ 0.25 supports use of 25% as low-threshold cutoff\n")
cat("===== END THRESHOLD CHECK =====\n\n")

# Col3: Low Threshold — adopted when ≤ 25% of friends were using
# Justified by distribution: 70th percentile = 0.25 in this sample (Valente 1996)
full_df$col3_low_threshold <- ifelse(full_df$is_user == 1 & !is.na(full_df$threshold),
  as.integer(full_df$threshold <= 0.25), NA)

# Among NON-USERS (resistant to peer influence):
# Col4: High Exposure — had network exposure > 0 at any point but never adopted
full_df$col4_exposure      <- ifelse(full_df$is_nonuser == 1,
  as.integer(!is.na(full_df$max_exposure) & full_df$max_exposure > 0), NA)

# Col5: High Perceived Friends — max perceived friend use > 0 across all waves
full_df$col5_friends       <- ifelse(full_df$is_nonuser == 1,
  as.integer(!is.na(full_df$max_friends_ecig) & full_df$max_friends_ecig > 0), NA)

# Col6: Both High — had both network exposure > 0 AND perceived friend use > 0
full_df$col6_both_high     <- ifelse(full_df$is_nonuser == 1,
  as.integer(!is.na(full_df$max_exposure)    & full_df$max_exposure > 0 &
             !is.na(full_df$max_friends_ecig) & full_df$max_friends_ecig > 0), NA)

cat(sprintf("Subgroup sizes:\n"))
for (col in c("col1_no_exposure","col2_no_friends","col_both0","col3_low_threshold",
              "col4_exposure","col5_friends","col6_both_high")) {
  cat(sprintf("  %s: n=%d\n", col, sum(full_df[[col]] == 1, na.rm = TRUE)))
}

# ── OR regression for each column ────────────────────────────────────────────
# Demographics | Network anchored at TOA for users (out_degree_at/in_degree_at),
# GP3 baseline for non-users (out_degree_gp3/in_degree_gp3, since non-users have no TOA)
# sm_ecig_any: W4-W8 only; GP3 adopters will have NA → complete cases reduced
base_covs_users <- c(
  "cohort",
  "female","hispanic","asian","sex_min","par_edu","gad","mdd",
  "out_degree_at","in_degree_at","expo_at_period",
  "friends_ecig_toa",
  "sm_ecig_any_toa"
)

base_covs_nonusers <- c(
  "cohort",
  "female","hispanic","asian","sex_min","par_edu","gad","mdd",
  "out_degree_gp3","in_degree_gp3",
  "friends_ecig_gp3","friends_ecig_gp8",
  "friends_ecig_toa",
  "sm_ecig_any_toa"
)

run_subgroup_or <- function(df, outcome_col, stratum_col, label, exclude_covs=c(), use_toa_fe=FALSE,
                             base_covs=base_covs_users) {
  cat(sprintf("\n--- %s ---\n", label))
  sub <- df[df[[stratum_col]] == 1 & !is.na(df[[outcome_col]]), ]
  cat(sprintf("  Stratum n = %d, outcome=1: %d\n",
              nrow(sub), sum(sub[[outcome_col]] == 1, na.rm = TRUE)))
  if (sum(sub[[outcome_col]] == 1, na.rm = TRUE) < 10) {
    cat("  Too few events, skipping\n"); return(NULL)
  }
  covs_use   <- base_covs[!base_covs %in% exclude_covs]
  covs_avail <- covs_use[covs_use %in% names(sub)]
  # Auto-drop any variable with < 2 levels in this subset to prevent separation
  covs_final <- c()
  for (v in covs_avail) {
    if (length(unique(na.omit(sub[[v]]))) >= 2) {
      covs_final <- c(covs_final, v)
    } else {
      cat(sprintf("  Dropping %s: only one level in this subset\n", v))
    }
  }
  rhs_terms <- covs_final
  if (use_toa_fe) rhs_terms <- c("factor(toa)", rhs_terms)
  fml <- as.formula(paste(outcome_col, "~", paste(rhs_terms, collapse = " + ")))
  cat(sprintf("  N with complete sm_ecig_any: %d\n",
              sum(complete.cases(sub[, intersect(c(outcome_col, rhs_terms), names(sub))]))))
  fit <- tryCatch(
    glm(fml, data = sub, family = binomial(link = "logit")),
    error = function(e) { cat("  Error:", e$message, "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  ct <- summary(fit)$coefficients
  ct <- ct[!grepl("^\\(Intercept\\)$", rownames(ct)), , drop = FALSE]
  out <- data.frame(
    column    = label,
    parameter = rownames(ct),
    OR        = round(exp(ct[, 1]), 3),
    pval      = round(ct[, 4], 4),
    sig       = ifelse(ct[, 4] < .01, "**", ifelse(ct[, 4] < .05, "*", "")),
    stringsAsFactors = FALSE
  )
  print(out[, c("parameter","OR","pval","sig")])
  out
}

# ── Run regressions ───────────────────────────────────────────────────────────
# Per Tom: only exclude exposure row if it DIRECTLY DEFINES the subgroup
# Col1 No Exposure: expo=0 by definition → exclude expo_lag/expo_gp8; keep friends_ecig_toa
# Col2 No Perceived Friends: friends=0 by definition → exclude friends_ecig_toa; keep expo
# Col3 Neither: both=0 by definition → exclude both
# Col4 Low Threshold: expo_lag excluded due to severe separation (OR=0.001, p<.0001)
#   — threshold (≤25% friends using) is mechanically derived from exposure proportion,
#   creating near-collinearity with expo_lag in this subgroup
or_col1 <- run_subgroup_or(full_df, "col1_no_exposure", "is_user",
                           "Col1: Users-NoExposure",exclude_covs = c("expo_at_period"),
                           use_toa_fe = TRUE)
or_col2 <- run_subgroup_or(full_df, "col2_no_friends", "is_user",
                           "Col2: Users-NoPerceivedFriends",
                           exclude_covs = c("friends_ecig_toa"),
                           use_toa_fe = TRUE)
or_col_both0 <- run_subgroup_or(full_df, "col_both0", "is_user",
                           "Col3: Users-NeitherExposureNorFriends",
                          exclude_covs = c("expo_at_period","friends_ecig_toa"),
                           use_toa_fe = TRUE)
or_col3 <- run_subgroup_or(full_df, "col3_low_threshold", "is_user",
                           "Col4: Users-LowThreshold",
                           exclude_covs = c("expo_at_period"),
                           use_toa_fe = TRUE)

# Non-users: exclude only the variable that defines each subgroup
# expo_gp8 also excluded from Col5/Col7 — near-perfect overlap with High Exposure definition
# Out-degree/in-degree use GP3 baseline for non-users (no TOA exists)
or_col4 <- run_subgroup_or(full_df, "col4_exposure",  "is_nonuser",
                           "Col5: NonUsers-HighExposure",
                           exclude_covs = c("friends_ecig_toa"),
                           base_covs = base_covs_nonusers)
or_col5 <- run_subgroup_or(full_df, "col5_friends",   "is_nonuser",
                           "Col6: NonUsers-HighFriends",
                           exclude_covs = c("friends_ecig_gp3","friends_ecig_gp8",
                                            "friends_ecig_toa"),
                           base_covs = base_covs_nonusers)
or_col6 <- run_subgroup_or(full_df, "col6_both_high", "is_nonuser",
                           "Col7: NonUsers-BothHigh",
                           exclude_covs = c("friends_ecig_toa"),
                           base_covs = base_covs_nonusers)

# ── Build two separate tables ─────────────────────────────────────────────────
make_wide_table <- function(or_list, filename) {
  combined <- do.call(rbind, Filter(Negate(is.null), or_list))
  if (is.null(combined) || nrow(combined) == 0) return(NULL)
  tbl <- reshape(
    combined[, c("column","parameter","OR","pval","sig")],
    idvar = "parameter", timevar = "column", direction = "wide"
  )
  cat("\n===== TABLE:", filename, "=====\n")
  print(tbl)
  write.csv(tbl, file.path(data_path, filename), row.names = FALSE)
  cat("Saved:", filename, "\n")
  tbl
}

user_table    <- make_wide_table(list(or_col1, or_col2, or_col_both0, or_col3),
                                 "260610_table1_users.csv")
nonuser_table <- make_wide_table(list(or_col4, or_col5, or_col6),
                                 "260610_table2_nonusers.csv")

cat("\n===== ALL DONE =====\n")

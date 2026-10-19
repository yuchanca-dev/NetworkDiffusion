# 261015 MINIMAL MISSING-OUTCOME FIX to original 260718 script.
# Only new sample criterion in run_tom_sm_event_reg(): observed 0/1 at outcome wave.
# All other original analysis choices remain unchanged; outputs isolated.
# =============================================================================
# 260704_FULL_CODE_WITH_STATA_STYLE_LOG.R
# Self-logging R script: runs the full analysis and saves a Stata-style .txt log
# =============================================================================
# How to use:
#   1. Put this file in your codes/ folder.
#   2. In RStudio, open this file and click Source, OR run:
#        source("260704_FULL_CODE_WITH_STATA_STYLE_LOG.R")
#   3. The full console log will be saved to ../Cleaned_Data/ as:
#        260704_R_full_run_log_YYYYMMDD_HHMMSS.txt
#
# The log captures:
#   - the R code that was executed (echo = TRUE)
#   - model output printed to console
#   - cat()/print() output
#   - messages, warnings, and errors
# =============================================================================

.get_this_script_path <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && nzchar(ofile)) {
      return(normalizePath(ofile, mustWork = TRUE))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(p) && nzchar(p)) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  stop("Could not determine this script path. Please run with source('260704_FULL_CODE_WITH_STATA_STYLE_LOG.R').")
}

# Parent run: start log, then source this same file in child mode with echo = TRUE.
# Child run: execute the actual analysis code in the else block below.
if (Sys.getenv("AO_SELF_LOG_CHILD") != "1") {
  script_path <- .get_this_script_path()
  out_dir <- "../Cleaned_Data"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_root <- file.path(out_dir, "261015_provisional_original_degree_results")
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  log_dir <- tempfile(paste0("run_",format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir=result_root)
  if (!dir.create(log_dir, recursive=TRUE)) stop("Could not create unique output directory")
  Sys.setenv(AO_FINAL_RUN_DIR=normalizePath(log_dir,mustWork=TRUE))
  log_file <- file.path(log_dir, "261015_FINAL_full_run_log.txt")
  log_con <- file(log_file, open = "wt")

  sink(log_con, split = TRUE)
  sink(log_con, type = "message")

  cat("===== R FULL RUN LOG STARTED =====\n")
  cat("Start time:", as.character(Sys.time()), "\n")
  cat("Working directory:", getwd(), "\n")
  cat("Script:", script_path, "\n")
  cat("Log file:", log_file, "\n\n")

  Sys.setenv(AO_SELF_LOG_CHILD = "1")
  run_error <- NULL
  tryCatch(
    {
      # Run the actual analysis in a separate environment so the child script's
      # rm(list = ls()) does not delete the parent logging objects (log_con, log_file).
      child_env <- new.env(parent = globalenv())
      source(script_path, local = child_env, echo = TRUE,
             max.deparse.length = Inf, keep.source = TRUE)
    },
    error = function(e) {
      run_error <<- e
      cat("\n===== ERROR OCCURRED =====\n")
      cat(conditionMessage(e), "\n")
      cat("\nCall:\n")
      print(conditionCall(e))
    }
  )
  Sys.unsetenv("AO_SELF_LOG_CHILD")

  cat("\n===== WARNINGS =====\n")
  print(warnings())

  cat("\n===== R FULL RUN LOG ENDED =====\n")
  cat("End time:", as.character(Sys.time()), "\n")

  sink(type = "message")
  sink()
  close(log_con)

  cat("Full R log saved to:", log_file, "\n")
  if (!is.null(run_error)) stop(run_error)
} else {

  # =============================================================================
  # 260617_netdiffuseR_efficient2_fixed-AO.R
  # =============================================================================
  rm(list = ls())
  library(netdiffuseR)
  library(openxlsx)
  library(lme4)        # only for the school-ICC diagnostic (1|school) — NOT for the main fit
  library(sandwich)    # cluster-robust covariance for the FE adoption model
  library(lmtest)      # coeftest() with the robust vcov
  # brglm2 provides the Firth-penalized fit used for the threshold model (problem 1).
  # (Also requires: lme4, sandwich, lmtest.) A local side-library is added to the
  # search path only if it exists; otherwise the default .libPaths() is used.
  .extra_lib <- if (dir.exists("playground/Rlib")) "playground/Rlib" else character(0)
  suppressMessages(library(brglm2, lib.loc = c(.extra_lib, .libPaths())))

  # ── Paths (simple relative path) ─────────────────────────────────────────────
  # Just put this script in the codes/ folder and run it directly.
  # It will read from ../Cleaned_Data/ and write results there.
  in_path   <- "../Cleaned_Data"
  out_path <- Sys.getenv("AO_FINAL_RUN_DIR", unset="")
  if (!nzchar(out_path)) stop("Missing unique run directory; source this complete file normally")
  if (!dir.exists(out_path)) stop("Run directory does not exist")
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

  # ── Data loading (AO) — SLIM RDS ONLY ─────────────────────────────────────────
  # This script NEVER touches the 153 MB ADVANCE workbooks. It loads a pre-built slim
  # cache (~468 KB) that holds only the columns the pipeline needs. Build it ONCE with
  # the companion prep script 260617_netdiffuseR_prepare_slim-AO.R; then ship just that
  # RDS + this script to reproduce every result. Override the location with AO_SLIM.
  cat("Loading ADVANCE data (slim RDS only)...\n")
  slim_path <- Sys.getenv("AO_SLIM", unset = file.path(in_path, "ADVANCE_slim_v1.rds"))
  if (!file.exists(slim_path))
    stop("Slim cache not found:\n  ", slim_path,
         "\nBuild it once with:\n  Rscript codes/260617_netdiffuseR_prepare_slim-AO.R",
         "\n(or set AO_SLIM to an existing ADVANCE_slim_v1.rds).")
  .slim  <- readRDS(slim_path)
  d_w18  <- .slim$w18
  d_w910 <- .slim$w910
  names(d_w18)  <- tolower(names(d_w18))
  names(d_w910) <- tolower(names(d_w910))
  cat(sprintf("  slim cache: %s\n", slim_path))
  cat(sprintf("  W1-W8:  %d students, %d columns\n",  nrow(d_w18),  ncol(d_w18)))
  cat(sprintf("  W9-W10: %d students, %d columns\n",  nrow(d_w910), ncol(d_w910)))

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
      # ── SOCIAL MEDIA: TOM'S METHOD (average + reverse-coded) ──────────────────
      # Tom Valente's method (from Stata analysis):
      # 1. Find all sm_post_* columns
      # 2. Recode 99 to NA (missing, not 0)
      # 3. Calculate mean across all platforms
      # 4. Reverse-code: 6 - mean (so higher = more frequent exposure to e-cig content)
    
      all_sm_cols <- grep(paste0("^w", w, "_sm_post_[0-9]+$"), names(r), value = TRUE)
      if (length(all_sm_cols) > 0) {
        sm_all_mat <- sapply(all_sm_cols, function(col) suppressWarnings(as.numeric(r[[col]])))
        sm_all_mat[sm_all_mat == 99] <- NA
      
        sm_post_mean <- apply(sm_all_mat, 1, function(x) {
          if (all(is.na(x))) return(NA_real_)
          mean(x, na.rm = TRUE)
        })
      
        sm_post_avg_rev <- 6 - sm_post_mean
      } else {
        sm_post_avg_rev <- rep(NA_real_, nrow(r))
      }

      # Instagram-only social media exposure: W4-W8 sm_post_2
      # Same direction as sm_post_avg_rev: higher = more frequent exposure.
      # Kept as a separate variable to compare with Tom's Stata Instagram model.
      insta_col <- paste0("w", w, "_sm_post_2")
      if (insta_col %in% names(r)) {
        insta_raw <- suppressWarnings(as.numeric(r[[insta_col]]))
        insta_raw[insta_raw == 99] <- NA
        insta_rev <- 6 - insta_raw
      } else {
        insta_rev <- rep(NA_real_, nrow(r))
      }

      # Snapchat-only social media exposure: sm_post_8
      # Same direction as the other social media variables:
      # higher values indicate more frequent exposure.
      snap_col <- paste0("w", w, "_sm_post_8")
      if (snap_col %in% names(r)) {
        snap_raw <- suppressWarnings(as.numeric(r[[snap_col]]))
        snap_raw[snap_raw == 99] <- NA
        snap_rev <- 6 - snap_raw
      } else {
        snap_rev <- rep(NA_real_, nrow(r))
      }

      tv_rows[[wi]] <- data.frame(
        id                    = vertex_ids,
        grade_period          = gp,
        gad                   = get_gad(r, w),
        mdd                   = get_mdd(r, w),
        friends_ecig          = get_friends_ecig(r, w),
        out_degree            = as.integer(out_deg),
        in_degree             = as.integer(in_deg),
        sm_post_avg_rev       = sm_post_avg_rev,
        insta_rev             = insta_rev,
        snap_rev              = snap_rev,
        sm_ecig_any           = NA_integer_,
        sm_ecig_any_all       = NA_integer_,
        sm_ecig_freq_max      = NA_real_,
        sm_ecig_platform_count= NA_real_,
        sm_ecig_friend        = NA_integer_,
        sm_ecig_nonfriend     = NA_integer_,
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

  # ── Lagged social media variables for Tom-requested lagged models ────────────
  # These represent the previous grade period's social media exposure. They are
  # used for sensitivity models where the outcome is at TOA but the social media
  # predictor is lagged by one time step.
  timevar_df <- timevar_df[order(timevar_df$id, timevar_df$grade_period), ]

  timevar_df$sm_post_avg_rev_lag <- ave(
    timevar_df$sm_post_avg_rev,
    timevar_df$id,
    FUN = function(x) c(NA, x[-length(x)])
  )

  timevar_df$insta_rev_lag <- ave(
    timevar_df$insta_rev,
    timevar_df$id,
    FUN = function(x) c(NA, x[-length(x)])
  )

  timevar_df$snap_rev_lag <- ave(
    timevar_df$snap_rev,
    timevar_df$id,
    FUN = function(x) c(NA, x[-length(x)])
  )

  cat("\n===== DIAGNOSTIC: lagged social media coverage =====\n")
  cat(sprintf("sm_post_avg_rev_lag non-NA rows: %d\n",
              sum(!is.na(timevar_df$sm_post_avg_rev_lag))))
  cat(sprintf("insta_rev_lag non-NA rows:       %d\n",
              sum(!is.na(timevar_df$insta_rev_lag))))
  cat(sprintf("snap_rev_lag non-NA rows:        %d\n",
              sum(!is.na(timevar_df$snap_rev_lag))))
  cat("===== END DIAGNOSTIC =====\n\n")

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
    predictors <- c("hispanic", "asian", "sex_min",
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
    predictors <- c("exposure", "hispanic", "asian", "sex_min",
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
    cat("\n===== THRESHOLD REGRESSION — GROUPED BINOMIAL (MDD ONLY; GAD + PERCEIVED FRIENDS EXCLUDED) =====\n")

    # Exclude isolates (n_alters = 0) — their threshold is uninformative
    thr_non_iso <- thr_sub[!is.na(thr_sub$n_alters) & thr_sub$n_alters > 0, ]
    cat(sprintf("Adopters with n_alters > 0 (non-isolates): %d\n", nrow(thr_non_iso)))
    cat(sprintf("Excluded isolates: %d\n", nrow(thr_sub) - nrow(thr_non_iso)))

    # school_c1fe: C2 schools collapsed to reference (101) to remove cohort collinearity
    thr_non_iso$school_c1fe <- make_school_c1fe(thr_non_iso$school)

    thr_covs <- c("factor(toa)", "cohort", "female", "hispanic", "asian",
                  "par_edu", "factor(school_c1fe)", "mdd")
    if (sum(!is.na(thr_non_iso$sex_min)) > 30 && length(unique(na.omit(thr_non_iso$sex_min))) > 1)
      thr_covs <- c(thr_covs, "sex_min")

    # grouped binomial: outcome = cbind(k_users, n_alters - k_users)
    fml <- as.formula(paste("cbind(k_users, n_alters - k_users) ~",
                            paste(thr_covs, collapse = " + ")))
    cat("Formula:", deparse(fml), "\n")

    req <- c("k_users", "n_alters", "cohort", "female", "hispanic", "asian",
             "par_edu", "school", "mdd")
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
  #exposure regression--count
  expo      <- exposure(diffnet_all, normalized = FALSE)

  # Proportion of nominated friends who had already adopted at each grade period.
  # Used only for the prior-semester (lagged) network low-threshold outcomes.
  expo_prop <- exposure(diffnet_all, normalized = TRUE)
  expo_prop[!is.finite(expo_prop)] <- NA_real_

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

  cat(sprintf("mdd NAs: %d out of %d rows\n",
              sum(is.na(reg_data$mdd)), nrow(reg_data)))

  # ══════════════════════════════════════════════════════════════════════════════
  # TOM-REQUESTED SOCIAL MEDIA EVENT-HISTORY MODELS
  # ══════════════════════════════════════════════════════════════════════════════
  # Tom requested 4 current-TOA models and the same 4 models with social media
  # lagged by one grade period:
  #   1) No-network-exposure users × average social media
  #   2) Low-threshold users (threshold <= 0.25) × average social media
  #   3) No-network-exposure users × Instagram-only
  #   4) Low-threshold users (threshold <= 0.25) × Instagram-only
  #
  # Covariates included per Tom + follow-up confirmation:
  #   time, cohort, female, sexual minority, parent education, Asian, Hispanic,
  #   major depressive symptoms, out-degree, in-degree, and school fixed effects.
  #
  # GAD is intentionally excluded in this Tom-requested set.
  # School FE note: original data include 14 schools. Because cohort is included,
  # C2 schools are collapsed to the reference school via make_school_c1fe(), giving
  # 10 effective school FE levels; one reference level is omitted in regression output.
  #
  # Model type: logistic regression with student-clustered robust SEs.
  # ══════════════════════════════════════════════════════════════════════════════

  run_tom_sm_event_reg <- function(expo, toa_vals, thr_vals, vertex_ids_clean,
                                   cov_df, attrs_df, timevar_df,
                                   outcome_type = c("no_exposure", "low_threshold"),
                                   sm_var = "sm_post_avg_rev",
                                   out_file = "tom_sm_model.csv",
                                   low_threshold_cut = 0.25,
                                   network_timing = c("current", "lagged"),
                                   expo_prop = NULL) {
    outcome_type <- match.arg(outcome_type)
    network_timing <- match.arg(network_timing)

    if (network_timing == "lagged" && outcome_type == "low_threshold" && is.null(expo_prop)) {
      stop("expo_prop is required for a lagged-network low-threshold model.")
    }

    cat(sprintf(
      "\n===== TOM SM MODEL: %s | network=%s | social media=%s =====\n",
      outcome_type, network_timing, sm_var
    ))

    if (!sm_var %in% names(timevar_df)) {
      stop(sprintf("%s not found in timevar_df. Check social media coding / lag creation section.", sm_var))
    }

    # Make TOA and current-TOA threshold vectors match cleaned IDs.
    toa_vec <- as.integer(toa_vals)
    names(toa_vec) <- vertex_ids_clean
    thr_vec <- as.numeric(thr_vals)
    names(thr_vec) <- vertex_ids_clean

    # Build person-period rows up to and including TOA. For lagged-network
    # specifications, the network measure attached to period t comes from t-1.
    tom_rows <- list()
    for (gi in seq_along(all_gps)) {
      gp <- all_gps[gi]

      current_count <- as.numeric(expo[, gi])
      lag_count <- if (gi > 1) as.numeric(expo[, gi - 1]) else rep(NA_real_, nrow(expo))
      lag_prop <- if (gi > 1 && !is.null(expo_prop)) {
        as.numeric(expo_prop[, gi - 1])
      } else {
        rep(NA_real_, nrow(expo))
      }

      tom_rows[[gi]] <- data.frame(
        id = vertex_ids_clean,
        grade_period = gp,
        exposure = current_count,
        network_exposure_lag = lag_count,
        network_threshold_lag = lag_prop,
        toa = as.integer(toa_vec[vertex_ids_clean]),
        threshold = as.numeric(thr_vec[vertex_ids_clean]),
        stringsAsFactors = FALSE
      )
    }

    tom_pp <- do.call(rbind, tom_rows)
    tom_pp <- tom_pp[is.na(tom_pp$toa) | tom_pp$grade_period <= tom_pp$toa, ]

    # Define the outcome at TOA using either same-semester or prior-semester
    # network exposure. Rows without the required prior network measure are set
    # to NA and removed by complete-case filtering rather than coded as controls.
    at_toa <- !is.na(tom_pp$toa) & tom_pp$grade_period == tom_pp$toa

    if (outcome_type == "no_exposure") {
      net_count <- if (network_timing == "current") {
        tom_pp$exposure
      } else {
        tom_pp$network_exposure_lag
      }

      network_available <- !is.na(net_count)
      tom_pp$outcome <- ifelse(
        network_available,
        as.integer(at_toa & net_count == 0),
        NA_integer_
      )
      outcome_label <- sprintf(
        "No-network-exposure adoption (%s network exposure)",
        network_timing
      )
    } else {
      if (network_timing == "current") {
        # Preserve the original current-TOA definition and results.
        tom_pp$outcome <- as.integer(
          at_toa &
            !is.na(tom_pp$threshold) &
            tom_pp$threshold <= low_threshold_cut
        )
      } else {
        network_available <- !is.na(tom_pp$network_threshold_lag)
        tom_pp$outcome <- ifelse(
          network_available,
          as.integer(at_toa & tom_pp$network_threshold_lag <= low_threshold_cut),
          NA_integer_
        )
      }
      outcome_label <- sprintf(
        "Low-threshold adoption (threshold <= %.2f; %s network exposure)",
        low_threshold_cut, network_timing
      )
    }

    cat(sprintf("Outcome: %s\n", outcome_label))
    cat(sprintf("Person-period rows before covariate merge: %d\n", nrow(tom_pp)))
    cat(sprintf("Rows with defined outcome: %d\n", sum(!is.na(tom_pp$outcome))))
    cat(sprintf("Events before covariate filtering: %d\n",
                sum(tom_pp$outcome == 1, na.rm = TRUE)))
    cat(sprintf("Unique students before filtering: %d\n",
                length(unique(tom_pp$id))))

    # Merge baseline/static and time-varying covariates.
    tom_pp <- merge(tom_pp, cov_df, by = "id", all.x = TRUE)
    tom_pp <- merge(tom_pp,
                    attrs_df[, c("id", "cohort", "sch_type", "school")],
                    by = "id", all.x = TRUE)
    tom_pp <- merge(tom_pp, timevar_df,
                    by.x = c("id", "grade_period"),
                    by.y = c("id", "grade_period"),
                    all.x = TRUE)

    tom_pp$school_c1fe <- make_school_c1fe(tom_pp$school)

    # 261015 MINIMAL CHANGE: require an ACTUALLY OBSERVED outcome at the event wave.
    # All original event-type definitions, earlier-history window, network degree
    # convention, covariates/timing, school FE and student-clustered inference stay
    # UNCHANGED. Wave = GP for C1; Wave = GP + 2 for C2.
    raw_id_for_check <- sub("^[^_]+_", "", as.character(tom_pp$id))
    outcome_wave_for_check <- ifelse(tom_pp$school %in% c2_schools,
                                     tom_pp$grade_period + 2L,
                                     tom_pp$grade_period)
    observed_outcome_for_check <- rep(NA_real_, nrow(tom_pp))
    for (ww in sort(unique(outcome_wave_for_check[!is.na(outcome_wave_for_check)]))) {
      selected <- which(outcome_wave_for_check == ww)
      src <- if (ww <= 8L) d_w18 else d_w910
      raw_col <- paste0("w", ww, "_past_6mo_use_3")
      if (!raw_col %in% names(src)) stop("Missing raw outcome field: ", raw_col)
      jj <- match(raw_id_for_check[selected], as.character(src$record_id))
      observed_outcome_for_check[selected] <- suppressWarnings(
        as.numeric(as.character(src[[raw_col]][jj]))
      )
    }
    valid_next_outcome <- !is.na(observed_outcome_for_check) &
                          observed_outcome_for_check %in% c(0, 1)
    old_defined <- !is.na(tom_pp$outcome)
    excluded <- sum(old_defined & !valid_next_outcome)
    excluded_controls <- sum(old_defined & !valid_next_outcome & tom_pp$outcome == 0)
    excluded_events <- sum(old_defined & !valid_next_outcome & tom_pp$outcome == 1)
    cat("MINIMAL MISSING-OUTCOME FIX: originally defined rows with no valid event-wave response: ",
        excluded, "\n", sep = "")
    cat("  Among those, old coded non-events: ",
        excluded_controls, "\n", sep = "")
    cat("  Among those, old coded events: ",
        excluded_events, "\n", sep = "")
    # Only additional sample restriction: unobserved event-wave outcome is unknown.
    tom_pp$outcome[!valid_next_outcome] <- NA_integer_

    covs <- c(
      "cohort",
      "female", "hispanic", "asian",
      "par_edu",
      "mdd",
      "out_degree", "in_degree",
      sm_var
    )
    if (sum(!is.na(tom_pp$sex_min)) > 100 &&
        length(unique(na.omit(tom_pp$sex_min))) > 1) {
      covs <- c(covs, "sex_min")
    }

    req_vars <- c(
      "outcome", "grade_period", "cohort",
      "female", "hispanic", "asian",
      "par_edu", "school",
      "mdd", "out_degree", "in_degree",
      sm_var
    )
    if ("sex_min" %in% covs) req_vars <- c(req_vars, "sex_min")

    reg_final <- tom_pp[
      complete.cases(tom_pp[, req_vars[req_vars %in% names(tom_pp)], drop = FALSE]),
    ]

    cat(sprintf("Person-period obs used: %d\n", nrow(reg_final)))
    cat(sprintf("Unique students used: %d\n", length(unique(reg_final$id))))
    cat(sprintf("Events used: %d\n", sum(reg_final$outcome == 1, na.rm = TRUE)))
    cat(sprintf("Effective school FE levels used: %d\n",
                length(unique(reg_final$school_c1fe))))
    cat("School FE levels: ", paste(sort(unique(reg_final$school_c1fe)), collapse = ", "), "\n", sep = "")

    if (sum(reg_final$outcome == 1, na.rm = TRUE) == 0) {
      stop("No events after filtering. Check outcome definition or missing social media/covariates.")
    }

    fml <- as.formula(paste(
      "outcome ~ factor(grade_period) +",
      paste(covs, collapse = " + "),
      "+ factor(school_c1fe)"
    ))
    cat("Formula:", deparse(fml), "\n")

    fit <- glm(fml, data = reg_final, family = binomial(link = "logit"))
    vc_cl <- sandwich::vcovCL(fit, cluster = reg_final$id)
    ctr <- lmtest::coeftest(fit, vcov. = vc_cl)

    rn <- gsub("factor\\(grade_period\\)", "GP",
               gsub("factor\\(school_c1fe\\)", "School", rownames(ctr)))
    keep <- rownames(ctr) != "(Intercept)"

    out <- data.frame(
      outcome_type = outcome_type,
      network_timing = network_timing,
      sm_variable = sm_var,
      parameter = rn[keep],
      estimate = round(ctr[keep, 1], 4),
      se       = round(ctr[keep, 2], 4),
      z        = round(ctr[keep, 3], 4),
      pval     = round(ctr[keep, 4], 6),
      OR       = round(exp(ctr[keep, 1]), 4),
      stringsAsFactors = FALSE
    )
    print(out)

    write.csv(out, file.path(data_path, out_file), row.names = FALSE)
    cat("Saved:", out_file, "\n")
    minimal_audit <- data.frame(
      outcome_type = outcome_type, network_timing = network_timing,
      sm_variable = sm_var, N_eligible_before_new_rule = sum(old_defined),
      N_excluded_missing_event_wave = excluded,
      N_excluded_old_non_events = excluded_controls,
      N_excluded_old_events = excluded_events,
      N_model = nrow(reg_final), N_students = length(unique(reg_final$id)),
      N_events = sum(reg_final$outcome == 1),
      SM_beta = if (sm_var %in% rownames(ctr)) unname(ctr[sm_var, 1]) else NA_real_,
      SM_SE = if (sm_var %in% rownames(ctr)) unname(ctr[sm_var, 2]) else NA_real_,
      SM_OR = if (sm_var %in% rownames(ctr)) exp(ctr[sm_var, 1]) else NA_real_,
      SM_CI_low = if (sm_var %in% rownames(ctr)) exp(ctr[sm_var, 1] - 1.96 * ctr[sm_var, 2]) else NA_real_,
      SM_CI_high = if (sm_var %in% rownames(ctr)) exp(ctr[sm_var, 1] + 1.96 * ctr[sm_var, 2]) else NA_real_,
      SM_p = if (sm_var %in% rownames(ctr)) unname(ctr[sm_var, 4]) else NA_real_
    )
    diagnostic_filename <- paste0("261015_missing_only_diagnostic_", sub("\\.csv$", "", out_file), ".csv")
    write.csv(minimal_audit, file.path(data_path, diagnostic_filename), row.names = FALSE)

    invisible(list(
      results = out,
      analytic_data = reg_final,
      model = fit,
      network_timing = network_timing
    ))
  }
  # ══════════════════════════════════════════════════════════════════════════════
  # TOM-REQUESTED CHECK: OVERALL E-CIGARETTE ADOPTION
  # Outcome = any first adoption at TOA, regardless of network exposure status.
  # Current overall social media and current Instagram exposure only.
  # Uses the same covariates, fixed effects, and student-clustered SEs as the
  # primary current-TOA models.
  # ══════════════════════════════════════════════════════════════════════════════
  
  run_overall_adoption_sm_reg <- function(
    sm_var,
    out_file
  ) {
    
    cat(
      sprintf(
        "\n===== OVERALL E-CIGARETTE ADOPTION MODEL: %s =====\n",
        sm_var
      )
    )
    
    if (!sm_var %in% names(timevar_df)) {
      stop(
        sprintf(
          "%s not found in timevar_df.",
          sm_var
        )
      )
    }
    
    # Clean TOA IDs so they match the covariate and time-varying datasets.
    toa_vec <- as.integer(toa_vals)
    names(toa_vec) <- vertex_ids_clean
    
    # Construct person-grade-semester rows up to and including adoption.
    overall_rows <- list()
    
    for (gi in seq_along(all_gps)) {
      
      gs <- all_gps[gi]
      
      overall_rows[[gi]] <- data.frame(
        id = vertex_ids_clean,
        grade_period = gs,
        toa = as.integer(toa_vec[vertex_ids_clean]),
        stringsAsFactors = FALSE
      )
    }
    
    overall_pp <- do.call(rbind, overall_rows)
    
    # Keep students at risk through their first adoption semester.
    overall_pp <- overall_pp[
      is.na(overall_pp$toa) |
        overall_pp$grade_period <= overall_pp$toa,
    ]
    
    # Event = first observed e-cigarette adoption, regardless of peer exposure.
    overall_pp$outcome <- as.integer(
      !is.na(overall_pp$toa) &
        overall_pp$grade_period == overall_pp$toa
    )
    
    cat(
      sprintf(
        "Person-period rows before covariate merge: %d\n",
        nrow(overall_pp)
      )
    )
    cat(
      sprintf(
        "Events before covariate filtering: %d\n",
        sum(overall_pp$outcome == 1, na.rm = TRUE)
      )
    )
    
    # Merge the same baseline and time-varying covariates used previously.
    overall_pp <- merge(
      overall_pp,
      cov_df,
      by = "id",
      all.x = TRUE
    )
    
    overall_pp <- merge(
      overall_pp,
      attrs_df[, c("id", "cohort", "sch_type", "school")],
      by = "id",
      all.x = TRUE
    )
    
    overall_pp <- merge(
      overall_pp,
      timevar_df,
      by = c("id", "grade_period"),
      all.x = TRUE
    )
    
    overall_pp$school_c1fe <- make_school_c1fe(overall_pp$school)
    
    covs <- c(
      "cohort",
      "female",
      "hispanic",
      "asian",
      "par_edu",
      "mdd",
      "out_degree",
      "in_degree",
      sm_var
    )
    
    if (
      sum(!is.na(overall_pp$sex_min)) > 100 &&
      length(unique(na.omit(overall_pp$sex_min))) > 1
    ) {
      covs <- c(covs, "sex_min")
    }
    
    req_vars <- c(
      "outcome",
      "grade_period",
      "cohort",
      "female",
      "hispanic",
      "asian",
      "par_edu",
      "school",
      "mdd",
      "out_degree",
      "in_degree",
      sm_var
    )
    
    if ("sex_min" %in% covs) {
      req_vars <- c(req_vars, "sex_min")
    }
    
    reg_final <- overall_pp[
      complete.cases(
        overall_pp[, req_vars, drop = FALSE]
      ),
    ]
    
    cat(sprintf("Person-period observations used: %d\n", nrow(reg_final)))
    cat(sprintf("Unique students used: %d\n", length(unique(reg_final$id))))
    cat(sprintf("Overall adoption events used: %d\n",
                sum(reg_final$outcome == 1, na.rm = TRUE)))
    
    if (sum(reg_final$outcome == 1, na.rm = TRUE) == 0) {
      stop("No adoption events remained after complete-case filtering.")
    }
    
    fml <- as.formula(
      paste(
        "outcome ~ factor(grade_period) +",
        paste(covs, collapse = " + "),
        "+ factor(school_c1fe)"
      )
    )
    
    cat("Formula:", deparse(fml), "\n")
    
    fit <- glm(
      fml,
      data = reg_final,
      family = binomial(link = "logit")
    )
    
    vc_cl <- sandwich::vcovCL(
      fit,
      cluster = reg_final$id
    )
    
    ctr <- lmtest::coeftest(
      fit,
      vcov. = vc_cl
    )
    
    rn <- gsub(
      "factor\\(grade_period\\)",
      "GS",
      gsub(
        "factor\\(school_c1fe\\)",
        "School",
        rownames(ctr)
      )
    )
    
    keep <- rownames(ctr) != "(Intercept)"
    
    out <- data.frame(
      outcome_type = "overall_adoption",
      sm_variable = sm_var,
      parameter = rn[keep],
      estimate = round(ctr[keep, 1], 4),
      se = round(ctr[keep, 2], 4),
      z = round(ctr[keep, 3], 4),
      pval = round(ctr[keep, 4], 6),
      OR = round(exp(ctr[keep, 1]), 4),
      stringsAsFactors = FALSE
    )
    
    print(out)
    
    write.csv(
      out,
      file.path(data_path, out_file),
      row.names = FALSE
    )
    
    cat("Saved:", out_file, "\n")
    
    invisible(
      list(
        results = out,
        analytic_data = reg_final,
        model = fit
      )
    )
  }
  
  
  # Current overall social media exposure
  tom_overall_adoption_allSM_res <- run_overall_adoption_sm_reg(
    sm_var = "sm_post_avg_rev",
    out_file = "260716_overall_adoption_allSM.csv"
  )
  
  # Current Instagram-specific exposure
  tom_overall_adoption_insta_res <- run_overall_adoption_sm_reg(
    sm_var = "insta_rev",
    out_file = "260716_overall_adoption_insta.csv"
  )
  
  
  # Extract only the two social media estimates for easy review.
  extract_overall_sm_row <- function(model_object, parameter_name, label) {
    
    r <- model_object$results
    row <- r[r$parameter == parameter_name, , drop = FALSE]
    d <- model_object$analytic_data
    
    data.frame(
      model = label,
      person_period_observations = nrow(d),
      unique_students = length(unique(d$id)),
      adoption_events = sum(d$outcome == 1, na.rm = TRUE),
      AOR = if (nrow(row) == 1) row$OR else NA_real_,
      SE = if (nrow(row) == 1) row$se else NA_real_,
      p_value = if (nrow(row) == 1) row$pval else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  tom_overall_adoption_summary <- rbind(
    extract_overall_sm_row(
      tom_overall_adoption_allSM_res,
      "sm_post_avg_rev",
      "Overall social media exposure"
    ),
    extract_overall_sm_row(
      tom_overall_adoption_insta_res,
      "insta_rev",
      "Instagram-specific exposure"
    )
  )
  
  cat("\n===== TOM OVERALL ADOPTION SUMMARY =====\n")
  print(tom_overall_adoption_summary)
  
  write.csv(
    tom_overall_adoption_summary,
    file.path(
      data_path,
      "260716_overall_adoption_SM_summary.csv"
    ),
    row.names = FALSE
  )
  
  cat("Saved: 260716_overall_adoption_SM_summary.csv\n")
  
  
  # Continue to the subgroup, lagged, and Snapchat-only models below.
  # stop("Overall-adoption models completed successfully. Intentional stop.")
  

  # ---- Current-TOA social media versions ----
  tom_current_noexp_allSM_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "sm_post_avg_rev",
    out_file = "260704_tom_current_noexposure_allSM.csv"
  )

  tom_current_lowthr_allSM_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "sm_post_avg_rev",
    out_file = "260704_tom_current_lowthreshold_allSM.csv"
  )

  tom_current_noexp_insta_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "insta_rev",
    out_file = "260704_tom_current_noexposure_insta.csv"
  )

  tom_current_lowthr_insta_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "insta_rev",
    out_file = "260704_tom_current_lowthreshold_insta.csv"
  )

  # ---- Current Snapchat-only versions ----
  tom_current_noexp_snap_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "snap_rev",
    out_file = "260718_tom_current_noexposure_snapchat.csv"
  )

  tom_current_lowthr_snap_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "snap_rev",
    out_file = "260718_tom_current_lowthreshold_snapchat.csv"
  )

  # ---- One-period lagged social media versions ----
  tom_lag_noexp_allSM_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "sm_post_avg_rev_lag",
    out_file = "260704_tom_lag_noexposure_allSM.csv"
  )

  tom_lag_lowthr_allSM_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "sm_post_avg_rev_lag",
    out_file = "260704_tom_lag_lowthreshold_allSM.csv"
  )

  tom_lag_noexp_insta_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "insta_rev_lag",
    out_file = "260704_tom_lag_noexposure_insta.csv"
  )

  tom_lag_lowthr_insta_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "insta_rev_lag",
    out_file = "260704_tom_lag_lowthreshold_insta.csv"
  )

  # ---- One-period lagged Snapchat-only versions ----
  tom_lag_noexp_snap_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure",
    sm_var = "snap_rev_lag",
    out_file = "260718_tom_lag_noexposure_snapchat.csv"
  )

  tom_lag_lowthr_snap_res <- run_tom_sm_event_reg(
    expo = expo, toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold",
    sm_var = "snap_rev_lag",
    out_file = "260718_tom_lag_lowthreshold_snapchat.csv"
  )

  # ══════════════════════════════════════════════════════════════════════════════
  # PRIOR-SEMESTER NETWORK EXPOSURE MODELS
  # Cross prior-semester network exposure with current and lagged social media.
  # This adds the two timing combinations not covered by the original models.
  # ══════════════════════════════════════════════════════════════════════════════

  # Lagged network + current social media
  tom_lagnet_current_noexp_allSM_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "sm_post_avg_rev",
    out_file = "260718_lagnet_currentSM_noexposure_allSM.csv"
  )

  tom_lagnet_current_lowthr_allSM_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "sm_post_avg_rev",
    out_file = "260718_lagnet_currentSM_lowthreshold_allSM.csv"
  )

  tom_lagnet_current_noexp_insta_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "insta_rev",
    out_file = "260718_lagnet_currentSM_noexposure_insta.csv"
  )

  tom_lagnet_current_lowthr_insta_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "insta_rev",
    out_file = "260718_lagnet_currentSM_lowthreshold_insta.csv"
  )

  tom_lagnet_current_noexp_snap_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "snap_rev",
    out_file = "260718_lagnet_currentSM_noexposure_snapchat.csv"
  )

  tom_lagnet_current_lowthr_snap_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "snap_rev",
    out_file = "260718_lagnet_currentSM_lowthreshold_snapchat.csv"
  )

  # Lagged network + lagged social media
  tom_lagnet_lagSM_noexp_allSM_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "sm_post_avg_rev_lag",
    out_file = "260718_lagnet_lagSM_noexposure_allSM.csv"
  )

  tom_lagnet_lagSM_lowthr_allSM_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "sm_post_avg_rev_lag",
    out_file = "260718_lagnet_lagSM_lowthreshold_allSM.csv"
  )

  tom_lagnet_lagSM_noexp_insta_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "insta_rev_lag",
    out_file = "260718_lagnet_lagSM_noexposure_insta.csv"
  )

  tom_lagnet_lagSM_lowthr_insta_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "insta_rev_lag",
    out_file = "260718_lagnet_lagSM_lowthreshold_insta.csv"
  )

  tom_lagnet_lagSM_noexp_snap_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "no_exposure", network_timing = "lagged",
    sm_var = "snap_rev_lag",
    out_file = "260718_lagnet_lagSM_noexposure_snapchat.csv"
  )

  tom_lagnet_lagSM_lowthr_snap_res <- run_tom_sm_event_reg(
    expo = expo, expo_prop = expo_prop,
    toa_vals = toa_vals, thr_vals = thr_vals,
    vertex_ids_clean = vertex_ids_clean,
    cov_df = cov_df, attrs_df = attrs_df, timevar_df = timevar_df,
    outcome_type = "low_threshold", network_timing = "lagged",
    sm_var = "snap_rev_lag",
    out_file = "260718_lagnet_lagSM_lowthreshold_snapchat.csv"
  )

  # Combined result table for the 12 subgroup social-media models.
  tom_sm_results_all <- do.call(rbind, list(
    tom_current_noexp_allSM_res$results,
    tom_current_lowthr_allSM_res$results,
    tom_current_noexp_insta_res$results,
    tom_current_lowthr_insta_res$results,
    tom_current_noexp_snap_res$results,
    tom_current_lowthr_snap_res$results,
    tom_lag_noexp_allSM_res$results,
    tom_lag_lowthr_allSM_res$results,
    tom_lag_noexp_insta_res$results,
    tom_lag_lowthr_insta_res$results,
    tom_lag_noexp_snap_res$results,
    tom_lag_lowthr_snap_res$results
  ))
  write.csv(tom_sm_results_all,
            file.path(data_path, "260718_tom_requested_12_SM_models_combined.csv"),
            row.names = FALSE)
  cat("Saved: 260718_tom_requested_12_SM_models_combined.csv\n")

  lagged_network_results_all <- do.call(rbind, list(
    tom_lagnet_current_noexp_allSM_res$results,
    tom_lagnet_current_lowthr_allSM_res$results,
    tom_lagnet_current_noexp_insta_res$results,
    tom_lagnet_current_lowthr_insta_res$results,
    tom_lagnet_current_noexp_snap_res$results,
    tom_lagnet_current_lowthr_snap_res$results,
    tom_lagnet_lagSM_noexp_allSM_res$results,
    tom_lagnet_lagSM_lowthr_allSM_res$results,
    tom_lagnet_lagSM_noexp_insta_res$results,
    tom_lagnet_lagSM_lowthr_insta_res$results,
    tom_lagnet_lagSM_noexp_snap_res$results,
    tom_lagnet_lagSM_lowthr_snap_res$results
  ))

  write.csv(
    lagged_network_results_all,
    file.path(data_path, "260718_lagged_network_12_models_combined.csv"),
    row.names = FALSE
  )
  cat("Saved: 260718_lagged_network_12_models_combined.csv\n")

  extract_timing_sm_row <- function(model_object, parameter_name, label) {
    r <- model_object$results
    row <- r[r$parameter == parameter_name, , drop = FALSE]
    d <- model_object$analytic_data

    data.frame(
      model = label,
      network_timing = unique(r$network_timing)[1],
      social_media_variable = parameter_name,
      person_period_observations = nrow(d),
      unique_students = length(unique(d$id)),
      events = sum(d$outcome == 1, na.rm = TRUE),
      AOR = if (nrow(row) == 1) row$OR else NA_real_,
      SE = if (nrow(row) == 1) row$se else NA_real_,
      p_value = if (nrow(row) == 1) row$pval else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  network_timing_summary <- rbind(
    extract_timing_sm_row(tom_lagnet_current_noexp_allSM_res, "sm_post_avg_rev",
                          "Lagged network + current overall SM: no exposure"),
    extract_timing_sm_row(tom_lagnet_current_lowthr_allSM_res, "sm_post_avg_rev",
                          "Lagged network + current overall SM: low threshold"),
    extract_timing_sm_row(tom_lagnet_current_noexp_insta_res, "insta_rev",
                          "Lagged network + current Instagram: no exposure"),
    extract_timing_sm_row(tom_lagnet_current_lowthr_insta_res, "insta_rev",
                          "Lagged network + current Instagram: low threshold"),
    extract_timing_sm_row(tom_lagnet_current_noexp_snap_res, "snap_rev",
                          "Lagged network + current Snapchat: no exposure"),
    extract_timing_sm_row(tom_lagnet_current_lowthr_snap_res, "snap_rev",
                          "Lagged network + current Snapchat: low threshold"),
    extract_timing_sm_row(tom_lagnet_lagSM_noexp_allSM_res, "sm_post_avg_rev_lag",
                          "Lagged network + lagged overall SM: no exposure"),
    extract_timing_sm_row(tom_lagnet_lagSM_lowthr_allSM_res, "sm_post_avg_rev_lag",
                          "Lagged network + lagged overall SM: low threshold"),
    extract_timing_sm_row(tom_lagnet_lagSM_noexp_insta_res, "insta_rev_lag",
                          "Lagged network + lagged Instagram: no exposure"),
    extract_timing_sm_row(tom_lagnet_lagSM_lowthr_insta_res, "insta_rev_lag",
                          "Lagged network + lagged Instagram: low threshold"),
    extract_timing_sm_row(tom_lagnet_lagSM_noexp_snap_res, "snap_rev_lag",
                          "Lagged network + lagged Snapchat: no exposure"),
    extract_timing_sm_row(tom_lagnet_lagSM_lowthr_snap_res, "snap_rev_lag",
                          "Lagged network + lagged Snapchat: low threshold")
  )

  cat("\n===== PRIOR-SEMESTER NETWORK TIMING SUMMARY =====\n")
  print(network_timing_summary)
  write.csv(
    network_timing_summary,
    file.path(data_path, "260718_prior_semester_network_timing_summary.csv"),
    row.names = FALSE
  )
  cat("Saved: 260718_prior_semester_network_timing_summary.csv\n")

  extract_snap_row <- function(model_object, parameter_name, label) {
    r <- model_object$results
    row <- r[r$parameter == parameter_name, , drop = FALSE]
    d <- model_object$analytic_data

    data.frame(
      model = label,
      person_period_observations = nrow(d),
      unique_students = length(unique(d$id)),
      events = sum(d$outcome == 1, na.rm = TRUE),
      AOR = if (nrow(row) == 1) row$OR else NA_real_,
      SE = if (nrow(row) == 1) row$se else NA_real_,
      p_value = if (nrow(row) == 1) row$pval else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  snapchat_summary <- rbind(
    extract_snap_row(tom_current_noexp_snap_res, "snap_rev",
                     "Current Snapchat: no network exposure"),
    extract_snap_row(tom_current_lowthr_snap_res, "snap_rev",
                     "Current Snapchat: low threshold"),
    extract_snap_row(tom_lag_noexp_snap_res, "snap_rev_lag",
                     "Lagged Snapchat: no network exposure"),
    extract_snap_row(tom_lag_lowthr_snap_res, "snap_rev_lag",
                     "Lagged Snapchat: low threshold")
  )

  cat("\n===== SNAPCHAT SUPPLEMENTAL SUMMARY =====\n")
  print(snapchat_summary)
  write.csv(
    snapchat_summary,
    file.path(data_path, "260718_snapchat_supplemental_summary.csv"),
    row.names = FALSE
  )
  cat("Saved: 260718_snapchat_supplemental_summary.csv\n")
  # ══════════════════════════════════════════════════════════════════════════════
  # ENHANCED TABLE 1: Characteristics of adopters in Tom-requested current-TOA models
  # ══════════════════════════════════════════════════════════════════════════════
  
  make_table1_for_tom_adopters <- function(noexp_data, lowthr_data, out_file) {
    
    # Keep only event/adopter rows from each analytic dataset
    noexp_adopters <- noexp_data[noexp_data$outcome == 1, ]
    lowthr_adopters <- lowthr_data[lowthr_data$outcome == 1, ]
    
    pct <- function(x) round(mean(x == 1, na.rm = TRUE) * 100, 1)
    mn  <- function(x) round(mean(x, na.rm = TRUE), 2)
    sdv <- function(x) round(sd(x, na.rm = TRUE), 2)
    
    make_col <- function(d) {
      c(
        "N adopters/events" = nrow(d),
        "N unique adopters" = length(unique(d$id)),
        "N unique students in analytic model" = NA,
        "N person-period observations in analytic model" = NA,
        
        "Cohort 2, %" = pct(d$cohort == "C2"),
        "Female, %" = pct(d$female),
        "Sexual minority, %" = pct(d$sex_min),
        "Asian, %" = pct(d$asian),
        "Hispanic/Latine, %" = pct(d$hispanic),
        
        "Parent education, mean" = mn(d$par_edu),
        "Parent education, SD" = sdv(d$par_edu),
        
        "Major depressive symptoms, mean" = mn(d$mdd),
        "Major depressive symptoms, SD" = sdv(d$mdd),
        
        "Out-degree, mean" = mn(d$out_degree),
        "Out-degree, SD" = sdv(d$out_degree),
        
        "In-degree, mean" = mn(d$in_degree),
        "In-degree, SD" = sdv(d$in_degree),
        
        "Overall social media exposure, mean" = mn(d$sm_post_avg_rev),
        "Overall social media exposure, SD" = sdv(d$sm_post_avg_rev),
        
        "Instagram-only exposure, mean" = mn(d$insta_rev),
        "Instagram-only exposure, SD" = sdv(d$insta_rev),
        
        "Adopted in GP3 / 10Fa, %" = pct(d$grade_period == 3),
        "Adopted in GP4 / 10Sp, %" = pct(d$grade_period == 4),
        "Adopted in GP5 / 11Fa, %" = pct(d$grade_period == 5),
        "Adopted in GP6 / 11Sp, %" = pct(d$grade_period == 6),
        "Adopted in GP7 / 12Fa, %" = pct(d$grade_period == 7),
        "Adopted in GP8 / 12Sp, %" = pct(d$grade_period == 8)
      )
    }
    
    table1 <- data.frame(
      Characteristic = names(make_col(noexp_adopters)),
      `No-network-exposure adopters` = as.numeric(make_col(noexp_adopters)),
      `Low-threshold adopters` = as.numeric(make_col(lowthr_adopters)),
      check.names = FALSE
    )
    
    # Add analytic sample N from the full model datasets
    table1[table1$Characteristic == "N unique students in analytic model",
           "No-network-exposure adopters"] <- length(unique(noexp_data$id))
    
    table1[table1$Characteristic == "N unique students in analytic model",
           "Low-threshold adopters"] <- length(unique(lowthr_data$id))
    
    table1[table1$Characteristic == "N person-period observations in analytic model",
           "No-network-exposure adopters"] <- nrow(noexp_data)
    
    table1[table1$Characteristic == "N person-period observations in analytic model",
           "Low-threshold adopters"] <- nrow(lowthr_data)
    
    write.csv(table1, file.path(data_path, out_file), row.names = FALSE)
    
    cat("\n===== ENHANCED TABLE 1 SAVED =====\n")
    cat("Saved:", out_file, "\n")
    print(table1)
    
    return(table1)
  }
  
  table1_tom_enhanced <- make_table1_for_tom_adopters(
    noexp_data = tom_current_noexp_allSM_res$analytic_data,
    lowthr_data = tom_current_lowthr_allSM_res$analytic_data,
    out_file = "260704_table1_tom_adopter_characteristics_enhanced.csv"
  )
  
  # ══════════════════════════════════════════════════════════════════════════════
  # TABLE 2: Current-TOA social media models
  # ══════════════════════════════════════════════════════════════════════════════
  
  format_or <- function(results_df, param_name) {
    row <- results_df[results_df$parameter == param_name, ]
    if (nrow(row) == 0) return("")
    
    or <- row$OR[1]
    p  <- row$pval[1]
    
    stars <- ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                           ifelse(p < 0.05, "*", "")))
    
    paste0(sprintf("%.2f", or), stars)
  }
  
  make_table2_current_models <- function(
    noexp_allSM_res,
    lowthr_allSM_res,
    noexp_insta_res,
    lowthr_insta_res,
    out_file = "260704_table2_current_TOA_social_media_models.csv"
  ) {
    
    # Results tables
    r1 <- noexp_allSM_res$results
    r2 <- lowthr_allSM_res$results
    r3 <- noexp_insta_res$results
    r4 <- lowthr_insta_res$results
    
    # Analytic data
    d1 <- noexp_allSM_res$analytic_data
    d2 <- lowthr_allSM_res$analytic_data
    d3 <- noexp_insta_res$analytic_data
    d4 <- lowthr_insta_res$analytic_data
    
    table2 <- data.frame(
      Characteristic = c(
        "N person-period observations",
        "N unique students",
        "N events",
        "Grade-period fixed effects",
        "School fixed effects",
        "Cohort 2",
        "Female",
        "Sexual minority",
        "Parent education",
        "Asian",
        "Hispanic/Latine",
        "Major depressive symptoms",
        "Out-degree",
        "In-degree",
        "Overall social media exposure",
        "Instagram-only exposure"
      ),
      
      `No-network-exposure adopters: Overall social media` = c(
        nrow(d1),
        length(unique(d1$id)),
        sum(d1$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r1, "cohortC2"),
        format_or(r1, "female"),
        format_or(r1, "sex_min"),
        format_or(r1, "par_edu"),
        format_or(r1, "asian"),
        format_or(r1, "hispanic"),
        format_or(r1, "mdd"),
        format_or(r1, "out_degree"),
        format_or(r1, "in_degree"),
        format_or(r1, "sm_post_avg_rev"),
        ""
      ),
      
      `Low-threshold adopters: Overall social media` = c(
        nrow(d2),
        length(unique(d2$id)),
        sum(d2$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r2, "cohortC2"),
        format_or(r2, "female"),
        format_or(r2, "sex_min"),
        format_or(r2, "par_edu"),
        format_or(r2, "asian"),
        format_or(r2, "hispanic"),
        format_or(r2, "mdd"),
        format_or(r2, "out_degree"),
        format_or(r2, "in_degree"),
        format_or(r2, "sm_post_avg_rev"),
        ""
      ),
      
      `No-network-exposure adopters: Instagram only` = c(
        nrow(d3),
        length(unique(d3$id)),
        sum(d3$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r3, "cohortC2"),
        format_or(r3, "female"),
        format_or(r3, "sex_min"),
        format_or(r3, "par_edu"),
        format_or(r3, "asian"),
        format_or(r3, "hispanic"),
        format_or(r3, "mdd"),
        format_or(r3, "out_degree"),
        format_or(r3, "in_degree"),
        "",
        format_or(r3, "insta_rev")
      ),
      
      `Low-threshold adopters: Instagram only` = c(
        nrow(d4),
        length(unique(d4$id)),
        sum(d4$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r4, "cohortC2"),
        format_or(r4, "female"),
        format_or(r4, "sex_min"),
        format_or(r4, "par_edu"),
        format_or(r4, "asian"),
        format_or(r4, "hispanic"),
        format_or(r4, "mdd"),
        format_or(r4, "out_degree"),
        format_or(r4, "in_degree"),
        "",
        format_or(r4, "insta_rev")
      ),
      
      check.names = FALSE
    )
    
    write.csv(table2, file.path(data_path, out_file), row.names = FALSE)
    
    cat("\n===== TABLE 2 CURRENT-TOA MODELS SAVED =====\n")
    cat("Saved:", out_file, "\n")
    print(table2)
    
    return(table2)
  }
  
  table2_current_TOA <- make_table2_current_models(
    noexp_allSM_res = tom_current_noexp_allSM_res,
    lowthr_allSM_res = tom_current_lowthr_allSM_res,
    noexp_insta_res = tom_current_noexp_insta_res,
    lowthr_insta_res = tom_current_lowthr_insta_res,
    out_file = "260704_table2_current_TOA_social_media_models.csv"
  )
  
  # ══════════════════════════════════════════════════════════════════════════════
  # TABLE 3: One-period lagged social media models
  # ══════════════════════════════════════════════════════════════════════════════
  
  format_or <- function(results_df, param_name) {
    row <- results_df[results_df$parameter == param_name, ]
    if (nrow(row) == 0) return("")
    
    or <- row$OR[1]
    p  <- row$pval[1]
    
    stars <- ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                           ifelse(p < 0.05, "*", "")))
    
    paste0(sprintf("%.2f", or), stars)
  }
  
  make_table3_lagged_models <- function(
    noexp_allSM_res,
    lowthr_allSM_res,
    noexp_insta_res,
    lowthr_insta_res,
    out_file = "260704_table3_lagged_social_media_models.csv"
  ) {
    
    # Results tables
    r1 <- noexp_allSM_res$results
    r2 <- lowthr_allSM_res$results
    r3 <- noexp_insta_res$results
    r4 <- lowthr_insta_res$results
    
    # Analytic data
    d1 <- noexp_allSM_res$analytic_data
    d2 <- lowthr_allSM_res$analytic_data
    d3 <- noexp_insta_res$analytic_data
    d4 <- lowthr_insta_res$analytic_data
    
    table3 <- data.frame(
      Characteristic = c(
        "N person-period observations",
        "N unique students",
        "N events",
        "Grade-period fixed effects",
        "School fixed effects",
        "Cohort 2",
        "Female",
        "Sexual minority",
        "Parent education",
        "Asian",
        "Hispanic/Latine",
        "Major depressive symptoms",
        "Out-degree",
        "In-degree",
        "Lagged overall social media exposure",
        "Lagged Instagram-only exposure"
      ),
      
      `No-network-exposure adopters: Lagged overall social media` = c(
        nrow(d1),
        length(unique(d1$id)),
        sum(d1$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r1, "cohortC2"),
        format_or(r1, "female"),
        format_or(r1, "sex_min"),
        format_or(r1, "par_edu"),
        format_or(r1, "asian"),
        format_or(r1, "hispanic"),
        format_or(r1, "mdd"),
        format_or(r1, "out_degree"),
        format_or(r1, "in_degree"),
        format_or(r1, "sm_post_avg_rev_lag"),
        ""
      ),
      
      `Low-threshold adopters: Lagged overall social media` = c(
        nrow(d2),
        length(unique(d2$id)),
        sum(d2$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r2, "cohortC2"),
        format_or(r2, "female"),
        format_or(r2, "sex_min"),
        format_or(r2, "par_edu"),
        format_or(r2, "asian"),
        format_or(r2, "hispanic"),
        format_or(r2, "mdd"),
        format_or(r2, "out_degree"),
        format_or(r2, "in_degree"),
        format_or(r2, "sm_post_avg_rev_lag"),
        ""
      ),
      
      `No-network-exposure adopters: Lagged Instagram only` = c(
        nrow(d3),
        length(unique(d3$id)),
        sum(d3$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r3, "cohortC2"),
        format_or(r3, "female"),
        format_or(r3, "sex_min"),
        format_or(r3, "par_edu"),
        format_or(r3, "asian"),
        format_or(r3, "hispanic"),
        format_or(r3, "mdd"),
        format_or(r3, "out_degree"),
        format_or(r3, "in_degree"),
        "",
        format_or(r3, "insta_rev_lag")
      ),
      
      `Low-threshold adopters: Lagged Instagram only` = c(
        nrow(d4),
        length(unique(d4$id)),
        sum(d4$outcome == 1, na.rm = TRUE),
        "Yes",
        "Yes",
        format_or(r4, "cohortC2"),
        format_or(r4, "female"),
        format_or(r4, "sex_min"),
        format_or(r4, "par_edu"),
        format_or(r4, "asian"),
        format_or(r4, "hispanic"),
        format_or(r4, "mdd"),
        format_or(r4, "out_degree"),
        format_or(r4, "in_degree"),
        "",
        format_or(r4, "insta_rev_lag")
      ),
      
      check.names = FALSE
    )
    
    write.csv(table3, file.path(data_path, out_file), row.names = FALSE)
    
    cat("\n===== TABLE 3 LAGGED MODELS SAVED =====\n")
    cat("Saved:", out_file, "\n")
    print(table3)
    
    return(table3)
  }
  
  table3_lagged <- make_table3_lagged_models(
    noexp_allSM_res = tom_lag_noexp_allSM_res,
    lowthr_allSM_res = tom_lag_lowthr_allSM_res,
    noexp_insta_res = tom_lag_noexp_insta_res,
    lowthr_insta_res = tom_lag_lowthr_insta_res,
    out_file = "260704_table3_lagged_social_media_models.csv"
  )
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
    cat("\n===== ADOPTION REGRESSION — FE GLM + cluster-robust SE (MDD ONLY; GAD AND PERCEIVED FRIENDS EXCLUDED) =====\n")
    reg_data$school_c1fe <- make_school_c1fe(reg_data$school)
    exp_covs <- c("cohort", "female", "hispanic", "asian", "par_edu",
                  "mdd")
    if (sum(!is.na(reg_data$sex_min)) > 100 && length(unique(na.omit(reg_data$sex_min))) > 1)
      exp_covs <- c(exp_covs, "sex_min")
    req_vars <- c("adopt_next", "exposure", "sm_post_avg_rev", "grade_period", "cohort",
                  "female", "hispanic", "asian", "par_edu", "school", "mdd")
    reg_final <- reg_data[complete.cases(reg_data[, req_vars[req_vars %in% names(reg_data)]]), ]
    cat(sprintf("Person-period obs: %d\n", nrow(reg_final)))
    n_students <- length(unique(reg_final$id))
    n_adopt    <- length(unique(reg_final$id[reg_final$adopt_next == 1]))
    cat(sprintf("Unique students: %d (adopters: %d, non-adopters: %d)\n",
                n_students, n_adopt, n_students - n_adopt))

    fml <- as.formula(paste("adopt_next ~ exposure + sm_post_avg_rev + factor(grade_period) +",
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
    for (k in intersect(c("exposure","mdd"), rownames(ctr)))
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
                    "factor(school_c1fe)", "mdd", "sex_min")
    req_vars <- c("adopt_next","exposure","sm_post_avg_rev","grade_period","cohort","female",
                  "hispanic","asian","par_edu","school","mdd")
    reg_final <- reg_data[complete.cases(reg_data[, req_vars[req_vars %in% names(reg_data)]]), ]
    cat(sprintf("Person-period obs: %d\n", nrow(reg_final)))
    fml <- as.formula(paste("adopt_next ~ exposure + sm_post_avg_rev + factor(grade_period) +",
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
    tom_nonexp_allSM       = "260702_tom_nonexp_TOA_event_history_allSM.csv",
    tom_nonexp_insta       = "260702_tom_nonexp_TOA_event_history_insta.csv",
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

  # Main decision table: original-model specification plus ONLY observed outcome fix.
  focus_models <- list(
    overall_sm_no_exposure = tom_lagnet_lagSM_noexp_allSM_res,
    overall_sm_low_threshold = tom_lagnet_lagSM_lowthr_allSM_res,
    instagram_no_exposure = tom_lagnet_lagSM_noexp_insta_res,
    instagram_low_threshold = tom_lagnet_lagSM_lowthr_insta_res
  )
  focus_benchmarks <- c(.046428, .040639, .029441, .030354)
  focused <- do.call(rbind, lapply(seq_along(focus_models), function(i) {
    nm <- names(focus_models)[i]; ob <- focus_models[[i]]
    v <- if (grepl("instagram", nm)) "insta_rev_lag" else "sm_post_avg_rev_lag"
    r <- ob$results[ob$results$parameter == v, , drop = FALSE]
    if (nrow(r) != 1L) stop("Cannot identify social media coefficient for ", nm)
    data.frame(model = nm, N = nrow(ob$analytic_data),
               students = length(unique(ob$analytic_data$id)),
               events = sum(ob$analytic_data$outcome == 1),
               OR = r$OR, p = r$pval,
               original_p_benchmark = focus_benchmarks[i],
               stringsAsFactors = FALSE)
  }))
  write.csv(focused, file.path(data_path, "261015_missing_only_FOUR_PRIMARY_RESULTS.csv"),
            row.names = FALSE)
  cat("\n===== FOUR PRIMARY MODELS: MISSING OUTCOME FIX ONLY =====\n")
  print(focused, row.names = FALSE)

  cat("\n===== DONE =====\n")
  cat("Outputs in:", data_path, "\n")
  cat("Excel: 260609_netdiffuseR_past6mo_ALL_RESULTS.xlsx\n")
  cat("          classification / node_outcomes\n")
  cat("Network PDF: 260609_netdiffuseR_past6mo_diffnet_allschools.pdf\n")
  cat("Tom-requested 8 SM models: 260704_tom_requested_8_SM_models_combined.csv\n")

  # Old TOM'S TABLE / 6-subgroup analysis removed in this clean version.
  # Current deliverables are the 8 person-period social-media models above.



  # =====================================================================
  # PROVISIONAL TWO-STAGE REPORT -- ORIGINAL DEGREE CONVENTION UNCHANGED
  # Stage 1 original <=25% event-type model; Stage 2 continuous peer
  # proportion with netdiffuseR's original zero-degree-as-zero convention.
  # IMPORTANT: 0/0 is mathematically undefined. Treat results as provisional
  # until collaborator clarifies definition of zero-degree ties.
  # =====================================================================
  cat("\n===== PROVISIONAL STAGE 1: ORIGINAL <=25% THRESHOLD =====\n")
  .safe_coef <- function(fit, d, term, label) {
    vv <- sandwich::vcovCL(fit, cluster = d$id)
    if (!term %in% names(coef(fit))) stop("Missing coefficient ", term, ": ", label)
    bb <- unname(coef(fit)[term]); ss <- sqrt(vv[term, term])
    data.frame(model = label, term = term, N = nrow(d),
      students_school_prefixed = length(unique(d$id)),
      events = sum(d$outcome == 1), beta = bb, robust_se = ss,
      AOR = exp(bb), CI_lower = exp(bb - 1.96*ss),
      CI_upper = exp(bb + 1.96*ss), p_value = 2*pnorm(-abs(bb/ss)),
      converged = fit$converged, stringsAsFactors = FALSE)
  }
  .gp_degree <- timevar_df[, c("id", "grade_period", "out_degree")]
  .gp_degree$grade_period <- .gp_degree$grade_period + 1L
  names(.gp_degree)[3] <- "matched_degree_previous_GP"
  if (anyDuplicated(.gp_degree[c("id", "grade_period")]))
    stop("Previous-GP degree lookup has duplicate student/GP pairs")

  threshold_results <- list()
  threshold_flow <- list()
  for (.name in names(focus_models)) {
    .o <- focus_models[[.name]]
    .term <- if (grepl("instagram", .name)) "insta_rev_lag" else "sm_post_avg_rev_lag"
    # DO NOT REFIT OR FILTER by matched degree. Use the original fitted model.
    threshold_results[[length(threshold_results)+1L]] <-
      .safe_coef(.o$model, .o$analytic_data, .term, .name)
    .d <- merge(.o$analytic_data, .gp_degree,
                by=c("id", "grade_period"), all.x=TRUE, sort=FALSE)
    if (nrow(.d) != nrow(.o$analytic_data)) stop("Degree audit join changed row count")
    .zero <- !is.na(.d$matched_degree_previous_GP) & .d$matched_degree_previous_GP == 0
    threshold_flow[[length(threshold_flow)+1L]] <- data.frame(
      model=.name, N=nrow(.d), events=sum(.d$outcome),
      zero_degree_periods=sum(.zero), zero_degree_events=sum(.d$outcome[.zero]),
      unknown_degree_periods=sum(is.na(.d$matched_degree_previous_GP)),
      stringsAsFactors=FALSE)
  }
  write.csv(do.call(rbind, threshold_results),
            file.path(data_path,"261015_PROVISIONAL_01_ORIGINAL_THRESHOLD.csv"),row.names=FALSE)
  write.csv(do.call(rbind, threshold_flow),
            file.path(data_path,"261015_PROVISIONAL_02_DEGREE_AUDIT_ONLY.csv"),row.names=FALSE)

  # Stage 2: Original risk set and netdiffuseR normalized peer proportion.
  # No positive matched-degree eligibility restriction is applied.
  # Keep the original zero-degree-as-zero value, but flag those rows.
  cat("\n===== PROVISIONAL STAGE 2: CONTINUOUS PEER (ORIGINAL DEGREE RULE) =====\n")
  .toa <- as.integer(toa_vals)
  names(.toa) <- vertex_ids_clean
  cont_rows <- list(); cont_flow <- list()
  for (.sm_name in c("overall_sm_low_threshold", "instagram_low_threshold")) {
    .var <- if (grepl("instagram", .sm_name)) "insta_rev_lag" else "sm_post_avg_rev_lag"
    .raw <- focus_models[[.sm_name]]$analytic_data
    .d <- merge(.raw, .gp_degree, by=c("id","grade_period"), all.x=TRUE, sort=FALSE)
    if (nrow(.d) != nrow(.raw)) stop("Continuous join changed row count")
    .d$first_adoption_gp <- unname(.toa[match(.d$id,names(.toa))])
    if (any(!.d$id %in% names(.toa))) stop("Unmatched student IDs in toa")
    .d$outcome <- as.integer(!is.na(.d$first_adoption_gp) &
                              .d$grade_period == .d$first_adoption_gp)
    if (!"network_threshold_lag" %in% names(.d)) stop("Lagged peer proportion missing")
    .d$peer_prop <- as.numeric(.d$network_threshold_lag)
    .d$zero_degree_flag <- !is.na(.d$matched_degree_previous_GP) &
                             .d$matched_degree_previous_GP == 0
    # Explicitly retain original netdiffuseR zero-degree=0 coding, if present.
    # Do not overwrite other missing normalized exposure values.
    .d$peer_prop[.d$zero_degree_flag] <- 0
    .kept <- is.finite(.d$peer_prop) & .d$peer_prop >= 0 & .d$peer_prop <= 1
    .d <- .d[.kept, , drop=FALSE]
    if (nrow(.d) < 100 || sum(.d$outcome) < 10)
      stop("Too few analyzable provisional continuous peer observations")
    .d$peer25 <- .d$peer_prop / 0.25
    .d$sm_center <- .d[[.var]] - mean(.d[[.var]])
    .sex <- "sex_min" %in% attr(terms(focus_models[[.sm_name]]$model),"term.labels")
    .terms <- c("factor(grade_period)","cohort","female","hispanic","asian",
                "par_edu","mdd","out_degree","in_degree",
                if (.sex) "sex_min", "factor(school_c1fe)")
    .f0 <- as.formula(paste("outcome ~ sm_center + peer25 +",paste(.terms,collapse=" + ")))
    .f1 <- as.formula(paste("outcome ~ sm_center * peer25 +",paste(.terms,collapse=" + ")))
    .m0 <- glm(.f0, data=.d, family=binomial(link="logit"),na.action=na.fail)
    .m1 <- glm(.f1, data=.d, family=binomial(link="logit"),na.action=na.fail)
    .sm_label <- if (grepl("instagram",.sm_name)) "instagram" else "overall_sm"
    for (.mod in c("additive","interaction")) {
      .fit <- if (.mod=="additive") .m0 else .m1
      for (.coef in c("sm_center","peer25",if (.mod=="interaction") "sm_center:peer25")) {
        cont_rows[[length(cont_rows)+1L]] <-
          .safe_coef(.fit,.d,.coef,paste0(.sm_label,"__",.mod))
      }
    }
    cont_flow[[length(cont_flow)+1L]] <- data.frame(
      social_media=.sm_label, original_analytic_person_periods=nrow(.raw),
      original_low_threshold_events=sum(.raw$outcome),
      provisional_continuous_person_periods=nrow(.d),
      students_school_prefixed=length(unique(.d$id)),
      overall_adoption_events=sum(.d$outcome),
      zero_degree_person_periods=sum(.d$zero_degree_flag),
      zero_degree_overall_events=sum(.d$outcome[.d$zero_degree_flag]),
      peer_zero_share=mean(.d$peer_prop==0),
      peer_leq_25_share=mean(.d$peer_prop<=0.25),
      additive_converged=.m0$converged,
      interaction_converged=.m1$converged)
  }
  write.csv(do.call(rbind,cont_rows),
            file.path(data_path,"261015_PROVISIONAL_03_CONTINUOUS.csv"),row.names=FALSE)
  write.csv(do.call(rbind,cont_flow),
            file.path(data_path,"261015_PROVISIONAL_04_CONTINUOUS_FLOW.csv"),row.names=FALSE)
  writeLines(c(
    "PROVISIONAL ANALYSIS -- KEEP ORIGINAL MATCHED-DEGREE CONVENTION UNTIL ANIBAL REVIEW.",
    "Stage 1 uses original fully lagged <=25% threshold event-type model, no degree exclusion.",
    "Stage 2 uses continuous lagged normalized matched-peer exposure and overall first-observed adoption.",
    "Stage 2 records with matched degree = 0 are provisionally assigned peer proportion = 0.",
    "This is the ORIGINAL computational convention, NOT a valid mathematical fraction (0/0).",
    "The zero-degree classification includes unmatched nominations or unknown survey completion.",
    "Compare with alternative degree policy later; do not treat provisional results as finalized.",
    "Original covariate timing, pooled C2 school FE, and school-prefixed ID clustering retained."
  ),file.path(data_path,"261015_PROVISIONAL_README.txt"))
  cat("\n===== PROVISIONAL TABLES SAVED =====\n")

}

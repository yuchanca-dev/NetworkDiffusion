rm(list = ls())
library(RSiena)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

cohort1_schools  <- c(101, 102, 103, 104, 105, 106, 107, 112, 113, 114)
cohort2_schools  <- c(201, 212, 213, 214)
asian_schools    <- c(103, 105, 112, 113)
hispanic_schools <- c(102, 106, 107, 114)
exclude_schools  <- c(108)

# Grade transition -> wave pair mapping
# C1: 9th=W1/W2, 10th=W3/W4, 11th=W5/W6, 12th=W7/W8
# C2: 9th=W3/W4, 10th=W5/W6, 11th=W7/W8
grade_wave_map <- list(
  list(grade = "9to10_Fall",    c1 = c(1,3), c2 = c(3,5)),
  list(grade = "10to11_Fall",   c1 = c(3,5), c2 = c(5,7)),
  list(grade = "11to12_Fall",   c1 = c(5,7), c2 = NULL),
  list(grade = "9to10_Spring",  c1 = c(2,4), c2 = c(4,6)),
  list(grade = "10to11_Spring", c1 = c(4,6), c2 = c(6,8)),
  list(grade = "11to12_Spring", c1 = c(6,8), c2 = NULL)
)

get_filename <- function(w, type = "data") {
  if (type == "data") {
    if (w == 9)  return(file.path(data_path, "w9adv_data_clean.csv"))
    if (w == 10) return(file.path(data_path, "w10adv_data.csv"))
    return(file.path(data_path, sprintf("w%d_adv_data.csv", w)))
  }
  if (type == "edge") {
    return(file.path(data_path, sprintf("w%dedges_clean.csv", w)))
  }
}

build_school_adj <- function(edgelist, student_ids) {
  edges <- edgelist[edgelist$ego %in% student_ids & edgelist$alter %in% student_ids, ]
  n   <- length(student_ids)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(student_ids, student_ids))
  if (nrow(edges) > 0) {
    for (k in seq_len(nrow(edges))) {
      i <- match(edges$ego[k], student_ids)
      j <- match(edges$alter[k], student_ids)
      if (!is.na(i) & !is.na(j)) mat[i, j] <- 1
    }
  }
  diag(mat) <- 0
  return(mat)
}

# Get ever-use of e-cigarettes (binary 0/1).
# Tries the binary _3a variable first (available W1-W4);
# falls back to the ordinal life_use_sub_3 variable (W1-W9), binarized as > 0.
get_ecig_ever <- function(df, w) {
  v_bin <- paste0("w", w, "_life_use_sub_3a")
  v_ord <- paste0("w", w, "_life_use_sub_3")
  if (v_bin %in% names(df)) {
    x <- suppressWarnings(as.numeric(df[[v_bin]]))
    return(as.integer(x == 1))
  }
  if (v_ord %in% names(df)) {
    x <- df[[v_ord]]
    x[x == ""] <- NA
    x_num <- suppressWarnings(as.numeric(x))
    return(as.integer(x_num > 0))
  }
  return(rep(NA_integer_, nrow(df)))
}

# Enforce monotonicity: ever-use cannot revert from 1 to 0.
enforce_monotone <- function(beh_mat) {
  revert_rows <- which(beh_mat[, 1] == 1 & beh_mat[, 2] == 0)
  if (length(revert_rows) > 0) beh_mat[revert_rows, 2] <- 1L
  return(beh_mat)
}

# Check if a covariate has meaningful variance (more than one unique value).
# FIX: asian_cov is all-zero in Hispanic-majority schools, which causes a
# coCovar() warning and makes the sameX effect inestimable.
# We skip the asian covariate entirely for such schools.
has_variance <- function(v) {
  vals <- v[!is.na(v)]
  length(unique(vals)) > 1
}

# ---------------------------------------------------------------
# Distribution table
# ---------------------------------------------------------------

distribution_table <- data.frame()

for (gm in grade_wave_map) {
  grade_label <- gm$grade

  school_wave_list <- list()
  for (sch in cohort1_schools) {
    if (sch %in% exclude_schools) next
    if (!is.null(gm$c1)) school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c1)
  }
  for (sch in cohort2_schools) {
    if (sch %in% exclude_schools) next
    if (!is.null(gm$c2)) school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c2)
  }

  loaded <- list()

  for (sw in school_wave_list) {
    sch <- sw$sch
    wt1 <- sw$wp[1]; wt2 <- sw$wp[2]

    f1 <- get_filename(wt1, "data"); f2 <- get_filename(wt2, "data")
    e1 <- get_filename(wt1, "edge"); e2 <- get_filename(wt2, "edge")

    cohort_label   <- ifelse(sch %in% cohort1_schools, "C1", "C2")
    sch_type_label <- ifelse(sch %in% asian_schools,    "Asian-majority",
                       ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))

    if (!all(file.exists(f1, f2, e1, e2))) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade = grade_label, school = sch, cohort = cohort_label,
        school_type = sch_type_label,
        n_students = NA, n_everuse_t1 = NA, n_everuse_t2 = NA,
        prev_t1_pct = NA, prev_t2_pct = NA,
        n_new_onset = NA, n_reversals_raw = NA,
        included = FALSE, skip_reason = "Data files missing",
        stringsAsFactors = FALSE
      ))
      next
    }

    cache_key1 <- paste0("d", wt1); cache_key2 <- paste0("d", wt2)
    if (is.null(loaded[[cache_key1]])) loaded[[cache_key1]] <- read.csv(f1)
    if (is.null(loaded[[cache_key2]])) loaded[[cache_key2]] <- read.csv(f2)

    s1 <- loaded[[cache_key1]]; s2 <- loaded[[cache_key2]]

    ids1 <- s1[[id_var]][s1[[schoolid_var]] == sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]] == sch & !is.na(s2[[schoolid_var]])]
    pids <- sort(intersect(ids1, ids2))

    if (length(pids) < 10) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade = grade_label, school = sch, cohort = cohort_label,
        school_type = sch_type_label,
        n_students = length(pids), n_everuse_t1 = NA, n_everuse_t2 = NA,
        prev_t1_pct = NA, prev_t2_pct = NA,
        n_new_onset = NA, n_reversals_raw = NA,
        included = FALSE,
        skip_reason = paste0("Too few students (n=", length(pids), ")"),
        stringsAsFactors = FALSE
      ))
      next
    }

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]

    beh_raw     <- cbind(get_ecig_ever(r1, wt1), get_ecig_ever(r2, wt2))
    n_reversals <- sum(beh_raw[,1] == 1 & beh_raw[,2] == 0, na.rm = TRUE)
    beh_mat     <- enforce_monotone(beh_raw)

    n_total <- nrow(beh_mat)
    n_t1    <- sum(beh_mat[,1] == 1, na.rm = TRUE)
    n_t2    <- sum(beh_mat[,2] == 1, na.rm = TRUE)
    prev_t1 <- n_t1 / n_total
    prev_t2 <- n_t2 / n_total
    n01     <- sum(beh_mat[,1] == 0 & beh_mat[,2] == 1, na.rm = TRUE)

    pass_n     <- n_t2 >= 10
    pass_prev  <- prev_t2 >= 0.10
    pass_trans <- n01 >= 10
    included   <- pass_n & pass_prev & pass_trans

    skip_reason <- ifelse(!pass_n,
      paste0("Too few ever-users at T2 (n=", n_t2, ")"),
      ifelse(!pass_prev,
        paste0("Prevalence too low at T2 (", round(prev_t2*100,1), "%)"),
        ifelse(!pass_trans,
          paste0("Too few new onsets (n01=", n01, ")"),
          "Included")))

    distribution_table <- rbind(distribution_table, data.frame(
      grade = grade_label, school = sch, cohort = cohort_label,
      school_type = sch_type_label,
      n_students = n_total,
      n_everuse_t1 = n_t1,
      n_everuse_t2 = n_t2,
      prev_t1_pct = round(prev_t1 * 100, 1),
      prev_t2_pct = round(prev_t2 * 100, 1),
      n_new_onset = n01,
      n_reversals_raw = n_reversals,
      included = included, skip_reason = skip_reason,
      stringsAsFactors = FALSE
    ))
  }
}

distribution_table <- distribution_table[order(distribution_table$grade,
                                                distribution_table$school), ]

write.csv(distribution_table,
          file.path(data_path, "siena_distribution_by_school_everuse.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------
# Build sienaData objects
# ---------------------------------------------------------------

school_data_list <- list()

for (gm in grade_wave_map) {
  grade_label <- gm$grade

  school_wave_list <- list()
  for (sch in cohort1_schools) {
    if (sch %in% exclude_schools) next
    if (!is.null(gm$c1)) school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c1)
  }
  for (sch in cohort2_schools) {
    if (sch %in% exclude_schools) next
    if (!is.null(gm$c2)) school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c2)
  }

  loaded <- list()

  for (sw in school_wave_list) {
    sch <- sw$sch
    wt1 <- sw$wp[1]; wt2 <- sw$wp[2]
    key <- paste(grade_label, sch, sep = "_")

    dist_row <- distribution_table[distribution_table$grade == grade_label &
                                   distribution_table$school == sch, ]
    if (nrow(dist_row) == 0 || !dist_row$included) next

    f1 <- get_filename(wt1, "data"); f2 <- get_filename(wt2, "data")
    e1 <- get_filename(wt1, "edge"); e2 <- get_filename(wt2, "edge")
    if (!all(file.exists(f1, f2, e1, e2))) next

    cache_key1 <- paste0("d", wt1); cache_key2 <- paste0("d", wt2)
    cache_e1   <- paste0("e", wt1); cache_e2   <- paste0("e", wt2)
    if (is.null(loaded[[cache_key1]])) loaded[[cache_key1]] <- read.csv(f1)
    if (is.null(loaded[[cache_key2]])) loaded[[cache_key2]] <- read.csv(f2)
    if (is.null(loaded[[cache_e1]]))   loaded[[cache_e1]]   <- read.csv(e1)
    if (is.null(loaded[[cache_e2]]))   loaded[[cache_e2]]   <- read.csv(e2)

    s1  <- loaded[[cache_key1]]; s2  <- loaded[[cache_key2]]
    ed1 <- loaded[[cache_e1]];   ed2 <- loaded[[cache_e2]]

    ids1 <- s1[[id_var]][s1[[schoolid_var]] == sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]] == sch & !is.na(s2[[schoolid_var]])]
    pids <- sort(intersect(ids1, ids2))

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]

    beh_mat <- enforce_monotone(
      cbind(get_ecig_ever(r1, wt1), get_ecig_ever(r2, wt2))
    )

    adj1      <- build_school_adj(ed1, pids)
    adj2      <- build_school_adj(ed2, pids)
    net_array <- array(c(adj1, adj2), dim = c(length(pids), length(pids), 2))

    sex_v  <- paste0("w", wt1, "_dem_gender")
    eth_v  <- paste0("w", wt1, "_eth")
    race_v <- paste0("w", wt1, "_race")
    sex_cov   <- as.integer(r1[[sex_v]])
    hisp_cov  <- as.integer(r1[[eth_v]]); hisp_cov[is.na(hisp_cov)] <- 0L
    asian_cov <- as.integer(as.integer(r1[[race_v]]) == 2); asian_cov[is.na(asian_cov)] <- 0L

    # FIX: only include asian_cov if it has actual variance in this school.
    # In Hispanic-majority schools asian_cov is all-zero, causing a coCovar()
    # warning and making the sameX(asian) effect inestimable.
    use_asian <- has_variance(asian_cov)

    sdata <- tryCatch({
      if (use_asian) {
        sienaDataCreate(
          siena_net   = sienaDependent(net_array),
          siena_nic   = sienaDependent(beh_mat, type = "behavior"),
          siena_sex   = coCovar(sex_cov),
          siena_hisp  = coCovar(hisp_cov),
          siena_asian = coCovar(asian_cov)
        )
      } else {
        sienaDataCreate(
          siena_net  = sienaDependent(net_array),
          siena_nic  = sienaDependent(beh_mat, type = "behavior"),
          siena_sex  = coCovar(sex_cov),
          siena_hisp = coCovar(hisp_cov)
        )
      }
    }, error = function(e) { NULL })

    if (!is.null(sdata)) {
      # Store whether asian was included so the effects loop knows
      attr(sdata, "use_asian") <- use_asian
      school_data_list[[key]] <- sdata
    }
  }
}

# ---------------------------------------------------------------
# SIENA algorithms
# FIX: removed thetaBound -- it is not a valid argument in this
# version of RSiena. Convergence is handled instead by the
# three-attempt retry logic and the avAlt fix-to-zero fallback.
# ---------------------------------------------------------------

# Primary algorithm
algo <- sienaAlgorithmCreate(
  projname = "advance_bygrade_everuse",
  n3       = 2000,   # manual Section 10: 2000-4000 for publication
  seed     = 42
)

# Retry algorithm: standard initial values (manual Section 11.2)
# useStdInits = TRUE sets density/reciprocity to sensible non-zero
# starting values and all other parameters to zero, which is more
# stable than continuing from a diverged parameter vector.
algo_std <- sienaAlgorithmCreate(
  projname    = "advance_bygrade_everuse_retry",
  n3          = 2000,
  seed        = 42,
  useStdInits = TRUE
)

results_list     <- list()
fixed_avalt_list <- list()

for (key in names(school_data_list)) {
  sdata     <- school_data_list[[key]]
  use_asian <- isTRUE(attr(sdata, "use_asian"))
  eff       <- getEffects(sdata)

  # Network evolution effects
  eff <- includeEffects(eff, transTrip)
  eff <- includeEffects(eff, sameX, interaction1 = "siena_sex")
  eff <- includeEffects(eff, sameX, interaction1 = "siena_hisp")
  # FIX: only add asian homophily when the variable has variance in this school
  if (use_asian) {
    eff <- includeEffects(eff, sameX, interaction1 = "siena_asian")
  }

  # Behavior evolution: avAlt = average behavior of alters (peer exposure).
  #eff <- includeEffects(eff, avAlt, name = "siena_nic", interaction1 = "siena_net")

  # --- Attempt 1: primary algorithm ---
  res  <- tryCatch(
    siena07(algo, data = sdata, effects = eff, batch = TRUE, verbose = FALSE,
            returnDeps = FALSE, initC = TRUE, useCluster = TRUE, nbrNodes = 4),
    error = function(e) { NULL }
  )
  conv <- if (!is.null(res)) max(abs(res$tconv.max), na.rm = TRUE) else Inf
  cat(sprintf("%s | attempt 1 | conv = %.3f\n", key, conv))

  # --- Attempt 2: restart from standard initial values (manual Section 11.2) ---
  if (conv >= 0.25) {
    res2 <- tryCatch(
      siena07(algo_std, data = sdata, effects = eff, batch = TRUE, verbose = FALSE,
              returnDeps = FALSE, initC = TRUE, useCluster = TRUE, nbrNodes = 4),
      error = function(e) { NULL }
    )
    conv2 <- if (!is.null(res2)) max(abs(res2$tconv.max), na.rm = TRUE) else Inf
    cat(sprintf("%s | attempt 2 (std inits) | conv = %.3f\n", key, conv2))
    if (conv2 < conv) { res <- res2; conv <- conv2 }
  }

  # --- Attempt 3: continue from best result so far using prevAns (manual Section 2.8) ---
  if (conv >= 0.25 && !is.null(res)) {
    res3 <- tryCatch(
      siena07(algo, data = sdata, effects = eff, batch = TRUE, verbose = FALSE,
              returnDeps = FALSE, initC = TRUE, useCluster = TRUE, nbrNodes = 4,
              prevAns = res),
      error = function(e) { NULL }
    )
    conv3 <- if (!is.null(res3)) max(abs(res3$tconv.max), na.rm = TRUE) else Inf
    cat(sprintf("%s | attempt 3 (prevAns) | conv = %.3f\n", key, conv3))
    if (conv3 < conv) { res <- res3; conv <- conv3 }
  }

  # --- Attempt 4: fix avAlt to 0 if its estimate is pathologically large ---
  # Manual Section 6.2.1: when both estimate and SE are large (>5), fix the parameter.
  # avAlt explodes under low prevalence because nearly all actors share state=0,
  # leaving too little information to estimate peer exposure.
  # The fixed model is stored separately and NOT used in the main meta-analysis
  # (consistent with Piombo et al. 2025: 3 schools excluded for inestimable avAlt).
  if (conv >= 0.25 && !is.null(res)) {
    avalt_row <- which(eff$shortName == "avAlt" & eff$include)
    if (length(avalt_row) > 0) {
      avalt_est <- tryCatch(abs(res$theta[avalt_row]), error = function(e) 0)
      if (!is.na(avalt_est) && avalt_est > 5) {
        eff_fixed <- eff
        eff_fixed$fix[avalt_row]          <- TRUE
        eff_fixed$initialValue[avalt_row] <- 0
        res_fx <- tryCatch(
          siena07(algo_std, data = sdata, effects = eff_fixed, batch = TRUE, verbose = FALSE,
                  returnDeps = FALSE, initC = TRUE, useCluster = TRUE, nbrNodes = 4),
          error = function(e) { NULL }
        )
        conv_fx <- if (!is.null(res_fx)) max(abs(res_fx$tconv.max), na.rm = TRUE) else Inf
        cat(sprintf("%s | attempt 4 (avAlt fixed=0) | conv = %.3f\n", key, conv_fx))
        if (!is.null(res_fx)) fixed_avalt_list[[key]] <- res_fx
      }
    }
  }

  if (!is.null(res)) {
    results_list[[key]] <- res
    cat(sprintf("%s | FINAL conv = %.3f | %s\n\n", key, conv,
                ifelse(conv < 0.25, "CONVERGED", "did not converge")))
  }
}

# ---------------------------------------------------------------
# Per-school results table
# ---------------------------------------------------------------

if (length(results_list) > 0) {
  per_school_table <- do.call(rbind, lapply(names(results_list), function(key) {
    parts    <- strsplit(key, "_")[[1]]
    sch      <- as.integer(tail(parts, 1))
    grade    <- paste(parts[1:(length(parts)-1)], collapse = "_")
    res      <- results_list[[key]]
    cohort   <- ifelse(sch %in% cohort1_schools, "C1", "C2")
    sch_type <- ifelse(sch %in% asian_schools,    "Asian-majority",
                  ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))
    conv     <- max(abs(res$tconv.max), na.rm = TRUE)
    effs     <- res$effects[res$effects$include, ]
    data.frame(
      grade       = grade,
      school      = sch,
      cohort      = cohort,
      school_type = sch_type,
      convergence = round(conv, 4),
      converged   = conv < 0.25,
      effectName  = paste(effs$shortName, effs$interaction1),
      estimate    = round(res$theta, 4),
      se          = round(res$se, 4),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(per_school_table,
            file.path(data_path, "siena_bygrade_results_per_school_everuse.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------------
# Meta-analysis by grade (>= 2 converged schools, avAlt not fixed)
# ---------------------------------------------------------------

converged_keys <- names(results_list)[sapply(results_list, function(r) {
  max(abs(r$tconv.max), na.rm = TRUE) < 0.25
})]

# Exclude schools where avAlt was fixed (no peer exposure estimate)
converged_keys_full <- converged_keys[!sapply(converged_keys, function(k) {
  any(results_list[[k]]$effects$shortName == "avAlt" &
      results_list[[k]]$effects$include &
      results_list[[k]]$effects$fix)
})]

grade_labels <- unique(sapply(strsplit(converged_keys_full, "_(?=[^_]+$)", perl=TRUE), `[`, 1))

meta_tables <- list()

for (gl in grade_labels) {
  keys_gl <- converged_keys_full[startsWith(converged_keys_full, gl)]
  if (length(keys_gl) < 2) next

  meta_res     <- siena08(results_list[keys_gl])
  m            <- meta_res
  ref_res      <- results_list[[keys_gl[1]]]
  effect_names <- ref_res$effects$effectName[ref_res$effects$include]

  meta_tables[[gl]] <- data.frame(
    grade      = gl,
    effectName = effect_names,
    estimate   = round(as.numeric(m$muhat),    4),
    se         = round(as.numeric(m$se.muhat), 4),
    ci.lb      = round(as.numeric(m$muhat) - 1.96 * as.numeric(m$se.muhat), 4),
    ci.ub      = round(as.numeric(m$muhat) + 1.96 * as.numeric(m$se.muhat), 4),
    pval       = round(2 * pnorm(-abs(as.numeric(m$muhat) / as.numeric(m$se.muhat))), 6),
    n_schools  = length(keys_gl),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

if (length(meta_tables) > 0) {
  meta_all <- do.call(rbind, meta_tables)
  write.csv(meta_all,
            file.path(data_path, "siena_bygrade_meta_output_everuse.csv"),
            row.names = FALSE)
}

# ---------------------------------------------------------------
# Goodness-of-fit for converged models
# ---------------------------------------------------------------

gof_results <- list()

for (key in converged_keys_full) {
  res <- results_list[[key]]
  gof_net <- tryCatch(
    sienaGOF(res, verbose = FALSE, varName = "siena_net",
             IndegreeDistribution, join = TRUE, cumulative = FALSE),
    error = function(e) { NULL }
  )
  gof_beh <- tryCatch(
    sienaGOF(res, verbose = FALSE, varName = "siena_nic",
             BehaviorDistribution, join = TRUE, cumulative = FALSE),
    error = function(e) { NULL }
  )
  gof_results[[key]] <- list(net = gof_net, beh = gof_beh)
}

# ---------------------------------------------------------------
# Convergence summary
# ---------------------------------------------------------------

cat("\n===== CONVERGENCE SUMMARY =====\n")
for (key in names(results_list)) {
  conv        <- max(abs(results_list[[key]]$tconv.max), na.rm = TRUE)
  avalt_fixed <- any(results_list[[key]]$effects$shortName == "avAlt" &
                     results_list[[key]]$effects$include &
                     results_list[[key]]$effects$fix)
  status <- if (conv < 0.25) "CONVERGED" else "NOT CONVERGED"
  note   <- if (avalt_fixed) " [avAlt fixed]" else ""
  cat(sprintf("  %-35s | conv = %.4f | %s%s\n", key, conv, status, note))
}
cat(sprintf("\nConverged (full model): %d / %d\n",
            length(converged_keys_full), length(results_list)))
if (length(fixed_avalt_list) > 0) {
  cat(sprintf("avAlt fixed models (excluded from meta): %d\n", length(fixed_avalt_list)))
  cat("  Keys:", paste(names(fixed_avalt_list), collapse = ", "), "\n")
}
cat("================================\n\n")

# ---------------------------------------------------------------
# Save all results
# ---------------------------------------------------------------

save(results_list, fixed_avalt_list, school_data_list,
     distribution_table, meta_tables, gof_results,
     file = file.path(data_path, "advance_bygrade_results_everuse.RData"))

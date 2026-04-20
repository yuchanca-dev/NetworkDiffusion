# rm(list = ls())
library(RSiena)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

# School metadata (for labeling results only, not used inside SIENA models)
cohort1_schools  <- c(101, 102, 103, 104, 105, 106, 107, 112, 113, 114)
cohort2_schools  <- c(201, 212, 213, 214)
asian_schools    <- c(103, 105, 112, 113)
hispanic_schools <- c(102, 106, 107, 114)

# Wave pairs: Fall 11→12th (W5→W7) and Spring 11→12th (W6→W8)
# Per Dr. Valente: highest prevalence at 11→12th grade transition
wave_pairs <- list(c(5, 7), c(6, 8))

# Helper: get correct filename (W9 and W10 have different naming)
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

# Helper: build adjacency matrix for one school
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

# Helper: get past30 variable 
get_nic <- function(df, v_any, v_freq) {
  if (v_any %in% names(df)) return(as.integer(df[[v_any]]))
  if (v_freq %in% names(df)) {
    x <- df[[v_freq]]
    x[x == ""] <- NA
    x_num <- suppressWarnings(as.numeric(x))
    return(as.integer(x_num > 0))
  }
  return(rep(NA_integer_, nrow(df)))
}

# ---- STEP 1: BUILD sienaData FOR EACH SCHOOL x WAVE PAIR ----

school_data_list <- list()

for (wp in wave_pairs) {
  wt1 <- wp[1]; wt2 <- wp[2]

  f1 <- get_filename(wt1, "data"); f2 <- get_filename(wt2, "data")
  e1 <- get_filename(wt1, "edge"); e2 <- get_filename(wt2, "edge")
  if (!all(file.exists(f1, f2, e1, e2))) {
    cat(sprintf("W%d-W%d: file missing, skipping\n", wt1, wt2)); next
  }

  s1 <- read.csv(f1); s2 <- read.csv(f2)
  ed1 <- read.csv(e1); ed2 <- read.csv(e2)

  stopifnot(!anyDuplicated(s1[[id_var]]))
  stopifnot(!anyDuplicated(s2[[id_var]]))

  # Past 30-day NIC/ecig use: try _any first, fallback to frequency dichotomized
  # Per Dr. Valente: use ever use to increase prevalence and reduce convergence issues
  ev_v1  <- paste0("w", wt1, "_past30_sub_3_any")
  ev_v2  <- paste0("w", wt2, "_past30_sub_3_any")
  sex_v  <- paste0("w", wt1, "_dem_gender")
  eth_v  <- paste0("w", wt1, "_eth")
  race_v <- paste0("w", wt1, "_race")

  schools <- sort(intersect(
    unique(s1[[schoolid_var]][!is.na(s1[[schoolid_var]])]),
    unique(s2[[schoolid_var]][!is.na(s2[[schoolid_var]])])
  ))

  for (sch in schools) {
    key <- paste(sch, wt1, wt2, sep = "_")

    ids1 <- s1[[id_var]][s1[[schoolid_var]] == sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]] == sch & !is.na(s2[[schoolid_var]])]
    pids <- sort(intersect(ids1, ids2))
    if (length(pids) < 10) next

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]

    freq_v1 <- paste0("w", wt1, "_past30_sub_3")
    freq_v2 <- paste0("w", wt2, "_past30_sub_3")
    beh_mat <- cbind(get_nic(r1, ev_v1, freq_v1), get_nic(r2, ev_v2, freq_v2))

    # Filter 1: minimum 10 users at each wave
    n_t1 <- sum(beh_mat[,1] == 1, na.rm = TRUE)
    n_t2 <- sum(beh_mat[,2] == 1, na.rm = TRUE)
    if (n_t1 < 10) { cat(sprintf("SKIPPED %s: only %d ever-users at T1\n", key, n_t1)); next }
    if (n_t2 < 10) { cat(sprintf("SKIPPED %s: only %d ever-users at T2\n", key, n_t2)); next }

    # Filter 2: minimum 10% prevalence at each wave
    prev_t1 <- n_t1 / nrow(beh_mat)
    prev_t2 <- n_t2 / nrow(beh_mat)
    if (prev_t1 < 0.10 | prev_t2 < 0.10) {
      cat(sprintf("SKIPPED %s: prevalence too low (%.1f%%, %.1f%%)\n",
                  key, prev_t1*100, prev_t2*100)); next
    }

    # Filter 3: minimum 10 behavior transitions (0→1 or 1→0)
    n01 <- sum(beh_mat[,1] == 0 & beh_mat[,2] == 1, na.rm = TRUE)
    n10 <- sum(beh_mat[,1] == 1 & beh_mat[,2] == 0, na.rm = TRUE)
    if ((n01 + n10) < 10) {
      cat(sprintf("SKIPPED %s: too few transitions (n01=%d, n10=%d)\n", key, n01, n10)); next
    }

    # Print distribution for records
    cat(sprintf("\n--- %s ever-use distribution ---\n", key))
    cat(sprintf("T1: %d ever-users / %d total (%.1f%%)\n", n_t1, nrow(beh_mat), prev_t1*100))
    cat(sprintf("T2: %d ever-users / %d total (%.1f%%)\n", n_t2, nrow(beh_mat), prev_t2*100))
    cat(sprintf("Transitions: 0→1=%d, 1→0=%d\n", n01, n10))

    adj1 <- build_school_adj(ed1, pids)
    adj2 <- build_school_adj(ed2, pids)
    net_array <- array(c(adj1, adj2), dim = c(length(pids), length(pids), 2))

    # Covariates: replace NA with 0 (per Dr. Valente)
    sex_cov   <- as.integer(r1[[sex_v]])
    hisp_cov  <- as.integer(r1[[eth_v]]);   hisp_cov[is.na(hisp_cov)]   <- 0L
    asian_cov <- as.integer(as.integer(r1[[race_v]]) == 2); asian_cov[is.na(asian_cov)] <- 0L

    sdata <- tryCatch({
      sienaDataCreate(
        siena_net   = sienaDependent(net_array),
        siena_beh   = sienaDependent(beh_mat, type = "behavior"),
        siena_sex   = coCovar(sex_cov),
        siena_hisp  = coCovar(hisp_cov),
        siena_asian = coCovar(asian_cov)
      )
    }, error = function(e) { cat(sprintf("sienaDataCreate failed for %s: %s\n", key, e$message)); NULL })

    if (!is.null(sdata))
      school_data_list[[key]] <- sdata
  }
}

cat(sprintf("\nTotal school-wave pairs to model: %d\n", length(school_data_list)))

# ---- STEP 2: RUN MODEL FOR EACH SCHOOL x WAVE PAIR ----

algo <- sienaAlgorithmCreate(projname = "advance_5768_everuse", n3 = 1000, seed = 42)

results_list <- list()

for (key in names(school_data_list)) {
  sdata <- school_data_list[[key]]
  eff   <- getEffects(sdata)

  # Network effects
  eff <- includeEffects(eff, transTrip)
  eff <- includeEffects(eff, sameX, interaction1 = "siena_sex")
  eff <- includeEffects(eff, sameX, interaction1 = "siena_hisp")
  eff <- includeEffects(eff, sameX, interaction1 = "siena_asian")

  # Behavior effects (peer influence = avSim, network exposure = indeg)
  # Note: effFrom removed due to instability in prior runs
  eff <- includeEffects(eff, avSim, name = "siena_beh", interaction1 = "siena_net")
  #eff <- includeEffects(eff, indeg, name = "siena_beh", interaction1 = "siena_net")

  res <- tryCatch(
    siena07(algo, data = sdata, effects = eff, batch = TRUE, verbose = FALSE,
            returnDeps = TRUE, initC = TRUE, useCluster = TRUE, nbrNodes = 4),
    error = function(e) { cat(sprintf("siena07 failed for %s: %s\n", key, e$message)); NULL }
  )

  if (!is.null(res)) {
    results_list[[key]] <- res
    parts  <- strsplit(key, "_")[[1]]
    sch    <- as.integer(parts[1])
    cohort <- ifelse(sch %in% cohort1_schools, "C1", "C2")
    conv   <- max(abs(res$tconv.max), na.rm = TRUE)
    cat(sprintf("[%s] School %s W%s-%s [%s]: conv=%.3f %s\n",
                Sys.time(), parts[1], parts[2], parts[3], cohort, conv,
                ifelse(conv < 0.25, "OK", "POOR")))
  }
}

# ---- STEP 3: EXTRACT PER-SCHOOL RESULTS TABLE ----

if (length(results_list) > 0) {
  per_school_table <- do.call(rbind, lapply(names(results_list), function(key) {
    parts    <- strsplit(key, "_")[[1]]
    sch      <- as.integer(parts[1])
    wt1      <- as.integer(parts[2])
    wt2      <- as.integer(parts[3])
    res      <- results_list[[key]]
    cohort   <- ifelse(sch %in% cohort1_schools, "C1", "C2")
    sch_type <- ifelse(sch %in% asian_schools,    "Asian-majority",
                 ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))
    conv     <- max(abs(res$tconv.max), na.rm = TRUE)
    effs     <- res$effects[res$effects$include, ]
    data.frame(
      school      = sch,
      wave_from   = wt1,
      wave_to     = wt2,
      cohort      = cohort,
      school_type = sch_type,
      convergence = round(conv, 4),
      effectName  = paste(effs$shortName, effs$interaction1),
      estimate    = round(res$theta, 4),
      se          = round(res$se, 4),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(per_school_table,
            file.path(data_path, "siena_5768_results_per_school.csv"),
            row.names = FALSE)
  cat("\nPer-school results saved.\n")
}

# ---- STEP 4: META-ANALYSIS (converged models only) ----

converged_keys <- names(results_list)[sapply(results_list, function(r) {
  max(abs(r$tconv.max), na.rm = TRUE) < 0.25
})]
cat(sprintf("\nConverged: %d / %d models\n", length(converged_keys), length(results_list)))

if (length(converged_keys) >= 2) {
  meta_results <- siena08(results_list[converged_keys])
  m <- meta_results

  # Get effect names from first converged model
  ref_res      <- results_list[[converged_keys[1]]]
  effect_names <- ref_res$effects$effectName[ref_res$effects$include]

  meta_table <- data.frame(
    effectName = effect_names,
    estimate   = round(as.numeric(m$muhat),    4),
    se         = round(as.numeric(m$se.muhat), 4),
    ci.lb      = round(as.numeric(m$muhat) - 1.96 * as.numeric(m$se.muhat), 4),
    ci.ub      = round(as.numeric(m$muhat) + 1.96 * as.numeric(m$se.muhat), 4),
    pval       = round(2 * pnorm(-abs(as.numeric(m$muhat) / as.numeric(m$se.muhat))), 6),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  write.csv(meta_table,
            file.path(data_path, "siena_5768_meta_output.csv"),
            row.names = FALSE)
  print(meta_table)
}

# ---- STEP 4.5: GOF FOR CONVERGED MODELS ----

gof_results <- list()

for (key in converged_keys) {
  res <- results_list[[key]]
  cat(sprintf("\nRunning GOF for: %s\n", key))

  gof_net <- tryCatch(
    sienaGOF(res, verbose = FALSE, varName = "siena_net",
             IndegreeDistribution, join = TRUE, cumulative = FALSE),
    error = function(e) { cat("GOF net failed:", e$message, "\n"); NULL }
  )

  gof_beh <- tryCatch(
    sienaGOF(res, verbose = FALSE, varName = "siena_beh",
             BehaviorDistribution, join = TRUE, cumulative = FALSE),
    error = function(e) { cat("GOF beh failed:", e$message, "\n"); NULL }
  )

  gof_results[[key]] <- list(net = gof_net, beh = gof_beh)

  if (!is.null(gof_net)) cat(sprintf("  Network GOF p = %.3f\n", gof_net$Joint$p.value))
  if (!is.null(gof_beh)) cat(sprintf("  Behavior GOF p = %.3f\n", gof_beh$Joint$p.value))
}

# ---- STEP 5: SAVE ALL OBJECTS ----

save(results_list, school_data_list, per_school_table, gof_results,
     file = file.path(data_path, "advance_5768_everuse_results.RData"))

cat("\nDone. Output files:\n")
cat(" - siena_5768_results_per_school.csv\n")
cat(" - siena_5768_meta_output.csv\n")
cat(" - advance_5768_everuse_results.RData\n")

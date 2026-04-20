rm(list = ls())
library(RSiena)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

# C1 only (W3-W8 = Grade 10 Fall to Grade 12 Spring, per Dr. Valente)
# C2 excluded: no Grade 12 data
cohort1_schools  <- c(101, 102, 103, 104, 105, 106, 107, 112, 113, 114)
asian_schools    <- c(103, 105, 112, 113)
hispanic_schools <- c(102, 106, 107, 114)
exclude_schools  <- c(108)

# 4 overlapping 3-wave windows across W3-W8 (C1 only)
# C1 grade mapping: W3=10Fall, W4=10Spr, W5=11Fall, W6=11Spr, W7=12Fall, W8=12Spr
grade_wave_map <- list(
  list(grade = "10to11_early", c1 = c(3,4,5)),  # 10Fall -> 10Spr -> 11Fall
  list(grade = "10to11_late",  c1 = c(4,5,6)),  # 10Spr  -> 11Fall -> 11Spr
  list(grade = "11to12_early", c1 = c(5,6,7)),  # 11Fall -> 11Spr  -> 12Fall
  list(grade = "11to12_late",  c1 = c(6,7,8))   # 11Spr  -> 12Fall -> 12Spr
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
      i <- match(edges$ego[k],   student_ids)
      j <- match(edges$alter[k], student_ids)
      if (!is.na(i) & !is.na(j)) mat[i, j] <- 1
    }
  }
  diag(mat) <- 0
  return(mat)
}

# Get ever-use of e-cigarettes (binary 0/1).
# Tries the binary _3a variable first (available W1-W4);
# falls back to the ordinal life_use_sub_3 variable, binarized as > 0.
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

# Enforce monotonicity across 3 waves: ever-use cannot revert from 1 to 0.
enforce_monotone_3w <- function(beh_mat) {
  rev12 <- which(beh_mat[,1] == 1 & beh_mat[,2] == 0)
  if (length(rev12) > 0) beh_mat[rev12, 2] <- 1L
  rev23 <- which(beh_mat[,2] == 1 & beh_mat[,3] == 0)
  if (length(rev23) > 0) beh_mat[rev23, 3] <- 1L
  return(beh_mat)
}

# Skip covariate if all-zero (e.g. asian_cov in Hispanic-majority schools)
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
    school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c1)
  }

  loaded <- list()

  for (sw in school_wave_list) {
    sch <- sw$sch
    wt1 <- sw$wp[1]; wt2 <- sw$wp[2]; wt3 <- sw$wp[3]

    f1 <- get_filename(wt1,"data"); f2 <- get_filename(wt2,"data"); f3 <- get_filename(wt3,"data")
    e1 <- get_filename(wt1,"edge"); e2 <- get_filename(wt2,"edge"); e3 <- get_filename(wt3,"edge")

    sch_type_label <- ifelse(sch %in% asian_schools,    "Asian-majority",
                       ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))

    if (!all(file.exists(f1,f2,f3,e1,e2,e3))) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade=grade_label, school=sch, school_type=sch_type_label,
        n_students=NA, n_everuse_t1=NA, n_everuse_t2=NA, n_everuse_t3=NA,
        prev_t1_pct=NA, prev_t2_pct=NA, prev_t3_pct=NA,
        n_new_onset=NA, n_reversals_raw=NA,
        included=FALSE, skip_reason="Data files missing",
        stringsAsFactors=FALSE))
      next
    }

    for (wt in c(wt1,wt2,wt3)) {
      ck <- paste0("d",wt)
      if (is.null(loaded[[ck]])) loaded[[ck]] <- read.csv(get_filename(wt,"data"))
    }

    s1 <- loaded[[paste0("d",wt1)]]
    s2 <- loaded[[paste0("d",wt2)]]
    s3 <- loaded[[paste0("d",wt3)]]

    ids1 <- s1[[id_var]][s1[[schoolid_var]] == sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]] == sch & !is.na(s2[[schoolid_var]])]
    ids3 <- s3[[id_var]][s3[[schoolid_var]] == sch & !is.na(s3[[schoolid_var]])]
    pids <- sort(Reduce(intersect, list(ids1, ids2, ids3)))

    if (length(pids) < 10) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade=grade_label, school=sch, school_type=sch_type_label,
        n_students=length(pids), n_everuse_t1=NA, n_everuse_t2=NA, n_everuse_t3=NA,
        prev_t1_pct=NA, prev_t2_pct=NA, prev_t3_pct=NA,
        n_new_onset=NA, n_reversals_raw=NA,
        included=FALSE,
        skip_reason=paste0("Too few students (n=",length(pids),")"),
        stringsAsFactors=FALSE))
      next
    }

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]
    r3 <- s3[match(pids, s3[[id_var]]), ]

    beh_raw <- cbind(get_ecig_ever(r1,wt1), get_ecig_ever(r2,wt2), get_ecig_ever(r3,wt3))
    n_reversals <- sum((beh_raw[,1]==1 & beh_raw[,2]==0) |
                       (beh_raw[,2]==1 & beh_raw[,3]==0), na.rm=TRUE)
    beh_mat <- enforce_monotone_3w(beh_raw)

    n_total <- nrow(beh_mat)
    n_t1    <- sum(beh_mat[,1]==1, na.rm=TRUE)
    n_t2    <- sum(beh_mat[,2]==1, na.rm=TRUE)
    n_t3    <- sum(beh_mat[,3]==1, na.rm=TRUE)
    prev_t3 <- n_t3 / n_total
    # New onsets: 0 at t1 -> 1 at t3
    n01     <- sum(beh_mat[,1]==0 & beh_mat[,3]==1, na.rm=TRUE)

    pass_n     <- n_t3 >= 10
    pass_prev  <- prev_t3 >= 0.10
    pass_trans <- n01  >= 10
    included   <- pass_n & pass_prev & pass_trans

    skip_reason <- ifelse(!pass_n,
      paste0("Too few ever-users at T3 (n=",n_t3,")"),
      ifelse(!pass_prev,
        paste0("Prevalence too low at T3 (",round(prev_t3*100,1),"%)"),
        ifelse(!pass_trans,
          paste0("Too few new onsets (n01=",n01,")"),
          "Included")))

    distribution_table <- rbind(distribution_table, data.frame(
      grade=grade_label, school=sch, school_type=sch_type_label,
      n_students=n_total,
      n_everuse_t1=n_t1, n_everuse_t2=n_t2, n_everuse_t3=n_t3,
      prev_t1_pct=round(n_t1/n_total*100,1),
      prev_t2_pct=round(n_t2/n_total*100,1),
      prev_t3_pct=round(prev_t3*100,1),
      n_new_onset=n01, n_reversals_raw=n_reversals,
      included=included, skip_reason=skip_reason,
      stringsAsFactors=FALSE))
  }
}

distribution_table <- distribution_table[order(distribution_table$grade,
                                                distribution_table$school), ]

write.csv(distribution_table,
          file.path(data_path, "0326bygrade_distribution.csv"),
          row.names=FALSE)

# ---------------------------------------------------------------
# Build sienaData objects (3-wave)
# ---------------------------------------------------------------

school_data_list <- list()

for (gm in grade_wave_map) {
  grade_label <- gm$grade

  school_wave_list <- list()
  for (sch in cohort1_schools) {
    if (sch %in% exclude_schools) next
    school_wave_list[[length(school_wave_list)+1]] <- list(sch=sch, wp=gm$c1)
  }

  loaded <- list()

  for (sw in school_wave_list) {
    sch <- sw$sch
    wt1 <- sw$wp[1]; wt2 <- sw$wp[2]; wt3 <- sw$wp[3]
    key <- paste(grade_label, sch, sep="_")

    dist_row <- distribution_table[distribution_table$grade==grade_label &
                                   distribution_table$school==sch, ]
    if (nrow(dist_row)==0 || !dist_row$included) next

    f1 <- get_filename(wt1,"data"); f2 <- get_filename(wt2,"data"); f3 <- get_filename(wt3,"data")
    e1 <- get_filename(wt1,"edge"); e2 <- get_filename(wt2,"edge"); e3 <- get_filename(wt3,"edge")
    if (!all(file.exists(f1,f2,f3,e1,e2,e3))) next

    for (wt in c(wt1,wt2,wt3)) {
      ck <- paste0("d",wt); ce <- paste0("e",wt)
      if (is.null(loaded[[ck]])) loaded[[ck]] <- read.csv(get_filename(wt,"data"))
      if (is.null(loaded[[ce]])) loaded[[ce]] <- read.csv(get_filename(wt,"edge"))
    }

    s1  <- loaded[[paste0("d",wt1)]]; s2  <- loaded[[paste0("d",wt2)]]; s3  <- loaded[[paste0("d",wt3)]]
    ed1 <- loaded[[paste0("e",wt1)]]; ed2 <- loaded[[paste0("e",wt2)]]; ed3 <- loaded[[paste0("e",wt3)]]

    ids1 <- s1[[id_var]][s1[[schoolid_var]]==sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]]==sch & !is.na(s2[[schoolid_var]])]
    ids3 <- s3[[id_var]][s3[[schoolid_var]]==sch & !is.na(s3[[schoolid_var]])]
    pids <- sort(Reduce(intersect, list(ids1, ids2, ids3)))

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]
    r3 <- s3[match(pids, s3[[id_var]]), ]

    # 3-column behavior matrix with monotonicity enforced
    beh_mat <- enforce_monotone_3w(
      cbind(get_ecig_ever(r1,wt1), get_ecig_ever(r2,wt2), get_ecig_ever(r3,wt3))
    )

    # 3-wave network array
    adj1 <- build_school_adj(ed1, pids)
    adj2 <- build_school_adj(ed2, pids)
    adj3 <- build_school_adj(ed3, pids)
    net_array <- array(c(adj1,adj2,adj3), dim=c(length(pids),length(pids),3))

    sex_v  <- paste0("w",wt1,"_dem_gender")
    eth_v  <- paste0("w",wt1,"_eth")
    race_v <- paste0("w",wt1,"_race")
    sex_cov   <- as.integer(r1[[sex_v]])
    hisp_cov  <- as.integer(r1[[eth_v]]); hisp_cov[is.na(hisp_cov)]   <- 0L
    asian_cov <- as.integer(as.integer(r1[[race_v]])==2); asian_cov[is.na(asian_cov)] <- 0L
    use_asian <- has_variance(asian_cov)

    sdata <- tryCatch({
      if (use_asian) {
        sienaDataCreate(
          siena_net   = sienaDependent(net_array),
          siena_nic   = sienaDependent(beh_mat, type="behavior"),
          siena_sex   = coCovar(sex_cov),
          siena_hisp  = coCovar(hisp_cov),
          siena_asian = coCovar(asian_cov)
        )
      } else {
        sienaDataCreate(
          siena_net  = sienaDependent(net_array),
          siena_nic  = sienaDependent(beh_mat, type="behavior"),
          siena_sex  = coCovar(sex_cov),
          siena_hisp = coCovar(hisp_cov)
        )
      }
    }, error = function(e) { NULL })

    if (!is.null(sdata)) {
      attr(sdata, "use_asian") <- use_asian
      school_data_list[[key]] <- sdata
    }
  }
}

# ---------------------------------------------------------------
# SIENA algorithms
# ---------------------------------------------------------------

algo <- sienaAlgorithmCreate(
  projname = "0326bygrade_everuse",
  n3       = 2000,
  seed     = 42
)

algo_std <- sienaAlgorithmCreate(
  projname    = "0326bygrade_everuse_retry",
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
  eff <- includeEffects(eff, sameX, interaction1="siena_sex")
  eff <- includeEffects(eff, sameX, interaction1="siena_hisp")
  if (use_asian) {
    eff <- includeEffects(eff, sameX, interaction1="siena_asian")
  }

  # Behavior evolution: avAlt = average peer exposure (key peer influence effect)
  eff <- includeEffects(eff, avAlt, name="siena_nic", interaction1="siena_net")
  eff <- includeEffects(eff, effFrom, name="siena_nic", interaction1="siena_sex")
  eff <- includeEffects(eff, effFrom, name="siena_nic", interaction1="siena_hisp")

  # --- Attempt 1: primary algorithm ---
  res  <- tryCatch(
    siena07(algo, data=sdata, effects=eff, batch=TRUE, verbose=FALSE,
            returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4),
    error = function(e) { NULL }
  )
  conv <- if (!is.null(res)) max(abs(res$tconv.max), na.rm=TRUE) else Inf
  cat(sprintf("%s | attempt 1 | conv = %.3f\n", key, conv))

  # --- Attempt 2: standard initial values ---
  if (conv >= 0.25) {
    res2 <- tryCatch(
      siena07(algo_std, data=sdata, effects=eff, batch=TRUE, verbose=FALSE,
              returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4),
      error = function(e) { NULL }
    )
    conv2 <- if (!is.null(res2)) max(abs(res2$tconv.max), na.rm=TRUE) else Inf
    cat(sprintf("%s | attempt 2 (std inits) | conv = %.3f\n", key, conv2))
    if (conv2 < conv) { res <- res2; conv <- conv2 }
  }

  # --- Attempt 3: continue from best result (prevAns) ---
  if (conv >= 0.25 && !is.null(res)) {
    res3 <- tryCatch(
      siena07(algo, data=sdata, effects=eff, batch=TRUE, verbose=FALSE,
              returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4,
              prevAns=res),
      error = function(e) { NULL }
    )
    conv3 <- if (!is.null(res3)) max(abs(res3$tconv.max), na.rm=TRUE) else Inf
    cat(sprintf("%s | attempt 3 (prevAns) | conv = %.3f\n", key, conv3))
    if (conv3 < conv) { res <- res3; conv <- conv3 }
  }

  # --- Attempt 4: fix avAlt=0 if estimate is pathologically large (>5) ---
  # Stored separately, excluded from main meta-analysis (per Piombo et al. 2025)
  if (conv >= 0.25 && !is.null(res)) {
    avalt_row <- which(eff$shortName == "avAlt" & eff$include)
    if (length(avalt_row) > 0) {
      avalt_est <- tryCatch(abs(res$theta[avalt_row]), error=function(e) 0)
      if (!is.na(avalt_est) && avalt_est > 5) {
        eff_fixed <- eff
        eff_fixed$fix[avalt_row]          <- TRUE
        eff_fixed$initialValue[avalt_row] <- 0
        res_fx <- tryCatch(
          siena07(algo_std, data=sdata, effects=eff_fixed, batch=TRUE, verbose=FALSE,
                  returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4),
          error = function(e) { NULL }
        )
        conv_fx <- if (!is.null(res_fx)) max(abs(res_fx$tconv.max), na.rm=TRUE) else Inf
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
    grade    <- paste(parts[1:(length(parts)-1)], collapse="_")
    res      <- results_list[[key]]
    sch_type <- ifelse(sch %in% asian_schools,    "Asian-majority",
                  ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))
    conv     <- max(abs(res$tconv.max), na.rm=TRUE)
    effs     <- res$effects[res$effects$include, ]
    data.frame(
      grade       = grade,
      school      = sch,
      school_type = sch_type,
      convergence = round(conv, 4),
      converged   = conv < 0.25,
      effectName  = paste(effs$shortName, effs$interaction1),
      estimate    = round(res$theta, 4),
      se          = round(res$se,    4),
      stringsAsFactors = FALSE
    )
  }))

  write.csv(per_school_table,
            file.path(data_path, "0326bygrade_results_per_school.csv"),
            row.names=FALSE)
}

# ---------------------------------------------------------------
# Meta-analysis by grade (>= 2 converged schools, avAlt not fixed)
# ---------------------------------------------------------------

converged_keys <- names(results_list)[sapply(results_list, function(r) {
  max(abs(r$tconv.max), na.rm=TRUE) < 0.25
})]

# Exclude schools where avAlt was fixed (no valid peer exposure estimate)
converged_keys_full <- converged_keys[!sapply(converged_keys, function(k) {
  any(results_list[[k]]$effects$shortName == "avAlt" &
      results_list[[k]]$effects$include    &
      results_list[[k]]$effects$fix)
})]

grade_labels <- unique(sapply(strsplit(converged_keys_full, "_(?=[^_]+$)", perl=TRUE), `[`, 1))
meta_tables  <- list()

for (gl in grade_labels) {
  keys_gl <- converged_keys_full[startsWith(converged_keys_full, gl)]
  if (length(keys_gl) < 2) next

  meta_res     <- siena08(results_list[keys_gl])
  ref_res      <- results_list[[keys_gl[1]]]
  effect_names <- ref_res$effects$effectName[ref_res$effects$include]

  meta_tables[[gl]] <- data.frame(
    grade      = gl,
    effectName = effect_names,
    estimate   = round(as.numeric(meta_res$muhat),    4),
    se         = round(as.numeric(meta_res$se.muhat), 4),
    ci.lb      = round(as.numeric(meta_res$muhat) - 1.96*as.numeric(meta_res$se.muhat), 4),
    ci.ub      = round(as.numeric(meta_res$muhat) + 1.96*as.numeric(meta_res$se.muhat), 4),
    pval       = round(2*pnorm(-abs(as.numeric(meta_res$muhat)/as.numeric(meta_res$se.muhat))), 6),
    n_schools  = length(keys_gl),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

if (length(meta_tables) > 0) {
  meta_all <- do.call(rbind, meta_tables)
  write.csv(meta_all,
            file.path(data_path, "0326bygrade_meta_output.csv"),
            row.names=FALSE)
}

# ---------------------------------------------------------------
# Goodness-of-fit for converged models
# ---------------------------------------------------------------

gof_results <- list()

for (key in converged_keys_full) {
  res <- results_list[[key]]
  gof_net <- tryCatch(
    sienaGOF(res, verbose=FALSE, varName="siena_net",
             IndegreeDistribution, join=TRUE, cumulative=FALSE),
    error = function(e) { NULL }
  )
  gof_beh <- tryCatch(
    sienaGOF(res, verbose=FALSE, varName="siena_nic",
             BehaviorDistribution, join=TRUE, cumulative=FALSE),
    error = function(e) { NULL }
  )
  gof_results[[key]] <- list(net=gof_net, beh=gof_beh)
}

# ---------------------------------------------------------------
# Convergence summary
# ---------------------------------------------------------------

cat("\n===== CONVERGENCE SUMMARY =====\n")
for (key in names(results_list)) {
  conv        <- max(abs(results_list[[key]]$tconv.max), na.rm=TRUE)
  avalt_fixed <- any(results_list[[key]]$effects$shortName == "avAlt" &
                     results_list[[key]]$effects$include              &
                     results_list[[key]]$effects$fix)
  status <- if (conv < 0.25) "CONVERGED" else "NOT CONVERGED"
  note   <- if (avalt_fixed) " [avAlt fixed]" else ""
  cat(sprintf("  %-40s | conv = %.4f | %s%s\n", key, conv, status, note))
}
cat(sprintf("\nConverged (full model): %d / %d\n",
            length(converged_keys_full), length(results_list)))
if (length(fixed_avalt_list) > 0) {
  cat(sprintf("avAlt fixed models (excluded from meta): %d\n", length(fixed_avalt_list)))
  cat("  Keys:", paste(names(fixed_avalt_list), collapse=", "), "\n")
}
cat("================================\n\n")

# ---------------------------------------------------------------
# Save all results
# ---------------------------------------------------------------

save(results_list, fixed_avalt_list, school_data_list,
     distribution_table, meta_tables, gof_results,
     file = file.path(data_path, "0326bygrade_everuse.RData"))

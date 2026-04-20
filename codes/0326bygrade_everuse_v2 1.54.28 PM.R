rm(list = ls())
library(RSiena)

data_path <- "/Users/yuchancao/Downloads/25fallRA/advance/Cleaned Data"

id_var       <- "record_id"
schoolid_var <- "schoolid"

# C1 only (W3-W8 = Grade 10 Fall to Grade 12 Spring, per Dr. Valente)
cohort1_schools  <- c(101, 102, 103, 104, 105, 106, 107, 112, 113, 114)
asian_schools    <- c(103, 105, 112, 113)
hispanic_schools <- c(102, 106, 107, 114)
exclude_schools  <- c(108)

# 4 overlapping 3-wave windows (C1 only, W3-W8)
grade_wave_map <- list(
  list(grade = "10to11_early", c1 = c(3,4,5)),
  list(grade = "10to11_late",  c1 = c(4,5,6)),
  list(grade = "11to12_early", c1 = c(5,6,7)),
  list(grade = "11to12_late",  c1 = c(6,7,8))
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

enforce_monotone_3w <- function(beh_mat) {
  rev12 <- which(beh_mat[,1] == 1 & beh_mat[,2] == 0)
  if (length(rev12) > 0) beh_mat[rev12, 2] <- 1L
  rev23 <- which(beh_mat[,2] == 1 & beh_mat[,3] == 0)
  if (length(rev23) > 0) beh_mat[rev23, 3] <- 1L
  return(beh_mat)
}

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
  loaded <- list()

  for (sch in cohort1_schools) {
    if (sch %in% exclude_schools) next

    wt1 <- gm$c1[1]; wt2 <- gm$c1[2]; wt3 <- gm$c1[3]
    f1 <- get_filename(wt1,"data"); f2 <- get_filename(wt2,"data"); f3 <- get_filename(wt3,"data")
    e1 <- get_filename(wt1,"edge"); e2 <- get_filename(wt2,"edge"); e3 <- get_filename(wt3,"edge")

    sch_type <- ifelse(sch %in% asian_schools, "Asian-majority",
                  ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))

    if (!all(file.exists(f1,f2,f3,e1,e2,e3))) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade=grade_label, school=sch, school_type=sch_type,
        n_students=NA, n_t1=NA, n_t2=NA, n_t3=NA,
        prev_t1=NA, prev_t2=NA, prev_t3=NA, n_new=NA,
        included=FALSE, skip_reason="Files missing", stringsAsFactors=FALSE))
      next
    }

    for (wt in c(wt1,wt2,wt3)) {
      ck <- paste0("d",wt)
      if (is.null(loaded[[ck]])) loaded[[ck]] <- read.csv(get_filename(wt,"data"))
    }

    s1 <- loaded[[paste0("d",wt1)]]
    s2 <- loaded[[paste0("d",wt2)]]
    s3 <- loaded[[paste0("d",wt3)]]

    ids1 <- s1[[id_var]][s1[[schoolid_var]]==sch & !is.na(s1[[schoolid_var]])]
    ids2 <- s2[[id_var]][s2[[schoolid_var]]==sch & !is.na(s2[[schoolid_var]])]
    ids3 <- s3[[id_var]][s3[[schoolid_var]]==sch & !is.na(s3[[schoolid_var]])]
    pids <- sort(Reduce(intersect, list(ids1,ids2,ids3)))

    if (length(pids) < 10) {
      distribution_table <- rbind(distribution_table, data.frame(
        grade=grade_label, school=sch, school_type=sch_type,
        n_students=length(pids), n_t1=NA, n_t2=NA, n_t3=NA,
        prev_t1=NA, prev_t2=NA, prev_t3=NA, n_new=NA,
        included=FALSE,
        skip_reason=paste0("Too few students (n=",length(pids),")"),
        stringsAsFactors=FALSE))
      next
    }

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]
    r3 <- s3[match(pids, s3[[id_var]]), ]

    beh <- enforce_monotone_3w(
      cbind(get_ecig_ever(r1,wt1), get_ecig_ever(r2,wt2), get_ecig_ever(r3,wt3))
    )

    n    <- nrow(beh)
    n_t1 <- sum(beh[,1]==1, na.rm=TRUE)
    n_t2 <- sum(beh[,2]==1, na.rm=TRUE)
    n_t3 <- sum(beh[,3]==1, na.rm=TRUE)
    n_new <- sum(beh[,1]==0 & beh[,3]==1, na.rm=TRUE)
    prev_t3 <- n_t3/n

    included <- (n_t3>=10) & (prev_t3>=0.10) & (n_new>=10)
    skip_reason <- ifelse(!included,
      paste0("t3: n=",n_t3,", prev=",round(prev_t3*100,1),"%, new=",n_new),
      "Included")

    distribution_table <- rbind(distribution_table, data.frame(
      grade=grade_label, school=sch, school_type=sch_type,
      n_students=n, n_t1=n_t1, n_t2=n_t2, n_t3=n_t3,
      prev_t1=round(n_t1/n*100,1),
      prev_t2=round(n_t2/n*100,1),
      prev_t3=round(prev_t3*100,1),
      n_new=n_new,
      included=included, skip_reason=skip_reason,
      stringsAsFactors=FALSE))
  }
}

distribution_table <- distribution_table[order(distribution_table$grade,
                                                distribution_table$school), ]
write.csv(distribution_table,
          file.path(data_path, "0326bygrade_v2_distribution.csv"),
          row.names=FALSE)

# ---------------------------------------------------------------
# Build sienaData objects
# ---------------------------------------------------------------

school_data_list <- list()

for (gm in grade_wave_map) {
  grade_label <- gm$grade
  loaded <- list()

  for (sch in cohort1_schools) {
    if (sch %in% exclude_schools) next
    key <- paste(grade_label, sch, sep="_")

    dist_row <- distribution_table[distribution_table$grade==grade_label &
                                   distribution_table$school==sch, ]
    if (nrow(dist_row)==0 || !dist_row$included) next

    wt1 <- gm$c1[1]; wt2 <- gm$c1[2]; wt3 <- gm$c1[3]
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
    pids <- sort(Reduce(intersect, list(ids1,ids2,ids3)))

    r1 <- s1[match(pids, s1[[id_var]]), ]
    r2 <- s2[match(pids, s2[[id_var]]), ]
    r3 <- s3[match(pids, s3[[id_var]]), ]

    beh_mat <- enforce_monotone_3w(
      cbind(get_ecig_ever(r1,wt1), get_ecig_ever(r2,wt2), get_ecig_ever(r3,wt3))
    )

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
    use_hisp  <- has_variance(hisp_cov)   # FALSE in Hispanic-majority schools

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
      attr(sdata, "use_hisp")  <- use_hisp
      school_data_list[[key]]  <- sdata
    }
  }
}

# ---------------------------------------------------------------
# SIENA algorithms
# thetaBound=100: prevents interactive pause when parameters explode.
# The 4-attempt retry logic below handles non-convergence automatically.
# ---------------------------------------------------------------

algo <- sienaAlgorithmCreate(
  projname = "0326bygrade_v2",
  n3       = 2000,
  seed     = 42
)

algo_std <- sienaAlgorithmCreate(
  projname    = "0326bygrade_v2_retry",
  n3          = 2000,
  seed        = 42,
  useStdInits = TRUE
)

results_list     <- list()
fixed_list       <- list()   # models with avAlt and/or effFrom(hisp) fixed=0

for (key in names(school_data_list)) {
  sdata     <- school_data_list[[key]]
  use_asian <- isTRUE(attr(sdata, "use_asian"))
  use_hisp  <- isTRUE(attr(sdata, "use_hisp"))
  eff       <- getEffects(sdata)

  # Network evolution
  eff <- includeEffects(eff, transTrip)
  eff <- includeEffects(eff, sameX, interaction1="siena_sex")
  eff <- includeEffects(eff, sameX, interaction1="siena_hisp")
  if (use_asian) {
    eff <- includeEffects(eff, sameX, interaction1="siena_asian")
  }

  # Behavior evolution
  # avAlt: peer exposure (core effect of interest)
  eff <- includeEffects(eff, avAlt, name="siena_nic", interaction1="siena_net")

  # effFrom(sex): almost always stable, always include
  eff <- includeEffects(eff, effFrom, name="siena_nic", interaction1="siena_sex")

  # effFrom(hisp): skip in Hispanic-majority schools (no variance -> explosion)
  if (use_hisp) {
    eff <- includeEffects(eff, effFrom, name="siena_nic", interaction1="siena_hisp")
  }

  run_model <- function(algo_obj, prev=NULL) {
    tryCatch(
      siena07(algo_obj, data=sdata, effects=eff, batch=TRUE, verbose=FALSE,
              returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4,
              thetaBound=100, prevAns=prev),
      error = function(e) NULL
    )
  }

  get_conv <- function(res) {
    if (is.null(res)) return(Inf)
    max(abs(res$tconv.max), na.rm=TRUE)
  }

  # Attempt 1: primary
  res  <- run_model(algo)
  conv <- get_conv(res)
  cat(sprintf("%s | attempt 1 | conv = %.3f\n", key, conv))

  # Attempt 2: standard initial values
  if (conv >= 0.25) {
    res2 <- run_model(algo_std)
    conv2 <- get_conv(res2)
    cat(sprintf("%s | attempt 2 (stdInits) | conv = %.3f\n", key, conv2))
    if (conv2 < conv) { res <- res2; conv <- conv2 }
  }

  # Attempt 3: continue from best so far
  if (conv >= 0.25 && !is.null(res)) {
    res3 <- run_model(algo, prev=res)
    conv3 <- get_conv(res3)
    cat(sprintf("%s | attempt 3 (prevAns) | conv = %.3f\n", key, conv3))
    if (conv3 < conv) { res <- res3; conv <- conv3 }
  }

  # Attempt 4: if avAlt is still explosive (>5), fix it to 0 and rerun.
  # Also fix effFrom(hisp) if it's explosive.
  # Per Piombo et al. 2025: fixed models stored separately, excluded from meta.
  if (conv >= 0.25 && !is.null(res)) {
    eff_fixed <- eff

    avalt_row <- which(eff$shortName=="avAlt" & eff$include)
    if (length(avalt_row) > 0 &&
        !is.na(res$theta[avalt_row]) && abs(res$theta[avalt_row]) > 5) {
      eff_fixed$fix[avalt_row]          <- TRUE
      eff_fixed$initialValue[avalt_row] <- 0
      cat(sprintf("%s | fixing avAlt=0\n", key))
    }

    hisp_row <- which(eff$shortName=="effFrom" &
                      eff$interaction1=="siena_hisp" & eff$include)
    if (length(hisp_row) > 0 &&
        !is.na(res$theta[hisp_row]) && abs(res$theta[hisp_row]) > 5) {
      eff_fixed$fix[hisp_row]          <- TRUE
      eff_fixed$initialValue[hisp_row] <- 0
      cat(sprintf("%s | fixing effFrom(hisp)=0\n", key))
    }

    # Only rerun if something was actually fixed
    any_fixed <- any(eff_fixed$fix & eff_fixed$include &
                     !(eff$fix & eff$include))
    if (any_fixed) {
      res_fx   <- tryCatch(
        siena07(algo_std, data=sdata, effects=eff_fixed, batch=TRUE, verbose=FALSE,
                returnDeps=FALSE, initC=TRUE, useCluster=TRUE, nbrNodes=4,
                thetaBound=100),
        error = function(e) NULL
      )
      conv_fx <- get_conv(res_fx)
      cat(sprintf("%s | attempt 4 (fixed) | conv = %.3f\n", key, conv_fx))
      if (!is.null(res_fx)) fixed_list[[key]] <- res_fx
      # Use fixed model as best if it converged and full model didn't
      if (conv_fx < conv) { res <- res_fx; conv <- conv_fx }
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
  per_school_rows <- lapply(names(results_list), function(key) {
    parts    <- strsplit(key, "_")[[1]]
    sch      <- as.integer(tail(parts, 1))
    grade    <- paste(parts[1:(length(parts)-1)], collapse="_")
    res      <- results_list[[key]]
    sch_type <- ifelse(sch %in% asian_schools, "Asian-majority",
                  ifelse(sch %in% hispanic_schools, "Hispanic-majority", "Other"))
    conv     <- max(abs(res$tconv.max), na.rm=TRUE)
    effs     <- res$effects[res$effects$include, ]
    data.frame(
      grade=grade, school=sch, school_type=sch_type,
      convergence=round(conv,4), converged=(conv<0.25),
      effectName=paste(effs$shortName, effs$interaction1),
      estimate=round(res$theta, 4),
      se=round(res$se, 4),
      fixed=effs$fix,
      stringsAsFactors=FALSE
    )
  })
  per_school_table <- do.call(rbind, per_school_rows)
  write.csv(per_school_table,
            file.path(data_path, "0326bygrade_v2_per_school.csv"),
            row.names=FALSE)
}

# ---------------------------------------------------------------
# Meta-analysis: only converged schools where avAlt was NOT fixed
# ---------------------------------------------------------------

converged_keys <- names(results_list)[sapply(results_list, function(r) {
  max(abs(r$tconv.max), na.rm=TRUE) < 0.25
})]

# Exclude any model where avAlt was fixed (inestimable peer influence)
converged_full <- converged_keys[!sapply(converged_keys, function(k) {
  any(results_list[[k]]$effects$shortName=="avAlt" &
      results_list[[k]]$effects$include &
      results_list[[k]]$effects$fix)
})]

grade_labels <- unique(sapply(
  strsplit(converged_full, "_(?=[^_]+$)", perl=TRUE), `[`, 1
))
meta_tables <- list()

for (gl in grade_labels) {
  keys_gl <- converged_full[startsWith(converged_full, gl)]
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
    pval       = round(2*pnorm(-abs(as.numeric(meta_res$muhat)/
                                    as.numeric(meta_res$se.muhat))), 6),
    n_schools  = length(keys_gl),
    row.names  = NULL, stringsAsFactors=FALSE
  )
}

if (length(meta_tables) > 0) {
  meta_all <- do.call(rbind, meta_tables)
  write.csv(meta_all,
            file.path(data_path, "0326bygrade_v2_meta.csv"),
            row.names=FALSE)
}

# ---------------------------------------------------------------
# Convergence summary
# ---------------------------------------------------------------

cat("\n===== CONVERGENCE SUMMARY =====\n")
for (key in names(results_list)) {
  res   <- results_list[[key]]
  conv  <- max(abs(res$tconv.max), na.rm=TRUE)
  fixed_params <- paste(
    res$effects$shortName[res$effects$include & res$effects$fix],
    collapse=", "
  )
  note   <- if (nchar(fixed_params)>0) paste0(" [fixed: ",fixed_params,"]") else ""
  status <- if (conv < 0.25) "CONVERGED" else "NOT CONVERGED"
  cat(sprintf("  %-40s | conv=%.4f | %s%s\n", key, conv, status, note))
}
cat(sprintf("\nTotal attempted: %d\n", length(results_list)))
cat(sprintf("Converged (full avAlt model): %d\n", length(converged_full)))
cat(sprintf("Fixed models (excluded from meta): %d\n", length(fixed_list)))
cat("================================\n\n")

# ---------------------------------------------------------------
# Save
# ---------------------------------------------------------------

save(results_list, fixed_list, school_data_list,
     distribution_table, meta_tables,
     file = file.path(data_path, "0326bygrade_v2_everuse.RData"))

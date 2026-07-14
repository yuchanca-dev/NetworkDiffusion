# 09-fully-corrected-tables-AO.R
# ─────────────────────────────────────────────────────────────────────────────
# Fully-corrected versions of the 260612 report tables (1a, 1b, 2, 3), applying
# EVERY fix identified in diffusion-study.md §2.4.1 and §5:
#
#   (i)  par_edu recoded — homogeneous legacy scale 1..6; "Don't know" (9) -> NA at
#        ALL waves; the spurious W7+ "new-scale" remap (which mis-collapsed 5,6 -> 4)
#        is dropped because in the 042326 release the column is legacy-scaled at every
#        wave. par_edu_corrected = first valid (1..6) value across the student's waves.
#   (ii) THRESHOLD + SUBGROUP/SOCIAL-MEDIA models fit by FIRTH penalized likelihood
#        (brglm2), which gives guaranteed-finite estimates under the GP3 / School-105
#        separation and corrects small-sample bias. Threshold excludes isolates only.
#   (iii)ADOPTION (school-FE GLM, Table 3) refit by Firth too (its GP8 dummy is
#        separated); corrected par_edu throughout. The mixed model (Table 2) is not
#        re-estimated — Firth does not extend to glmer, and §4 shows that model is
#        numerically unstable; the Firth-FE adoption fit is the recommended stable form.
#
# INPUT : cached frames thr_sub / full_df / reg_data (Yuchan's 260612 pipeline) +
#         the two ADVANCE XLSX (for the par_edu recompute).
# OUTPUT: outputs_AO/model/fully_corrected_{threshold,adoption,table1a,table1b}-AO.csv
# ─────────────────────────────────────────────────────────────────────────────

suppressWarnings(suppressMessages(library(openxlsx)))
cache_dir <- Sys.getenv("AO_CACHE", unset = "playground/outputs_260612")
in_path   <- Sys.getenv("AO_IN",    unset = "playground/data")
firth_lib <- "playground/Rlib"
stopifnot(requireNamespace("brglm2", quietly = TRUE, lib.loc = c(firth_lib, .libPaths())))
suppressMessages(library(brglm2, lib.loc = c(firth_lib, .libPaths())))

rd <- function(n) readRDS(file.path(cache_dir, paste0("AO_", n, ".rds")))
thr_sub <- rd("thr_sub"); full_df <- rd("full_df"); reg_data <- rd("reg_data")

# ── (i) corrected par_edu ─────────────────────────────────────────────────────
d18  <- read.xlsx(file.path(in_path, "ADVANCE_W1-W8.xlsx"),  sheet = 1); names(d18)  <- tolower(names(d18))
d910 <- read.xlsx(file.path(in_path, "ADVANCE_W9-W10.xlsx"), sheet = 1); names(d910) <- tolower(names(d910))
recode_par <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(!is.na(x) & x >= 1 & x <= 6, x, NA_real_) }
par_mat <- sapply(1:10, function(w) {
  d <- if (w <= 8) d18 else d910
  col <- paste0("w", w, "_dem_high_par_edu")
  if (col %in% names(d)) recode_par(d[[col]])[match(d18$record_id, d$record_id)] else rep(NA_real_, nrow(d18))
})
par_corr <- apply(par_mat, 1, function(r) { v <- r[!is.na(r)]; if (length(v)) v[1] else NA_real_ })
par_map  <- data.frame(record_id = as.character(d18$record_id), par_edu_c = par_corr, stringsAsFactors = FALSE)

attach_par <- function(df, idcol) {
  rid <- sub("^[^_]+_", "", as.character(df[[idcol]]))
  df$par_edu <- par_map$par_edu_c[match(rid, par_map$record_id)]   # OVERWRITE with corrected
  df
}
thr_sub  <- attach_par(thr_sub, "id_orig")
full_df  <- attach_par(full_df, "id_orig")
reg_data <- attach_par(reg_data, "id")
cat(sprintf("par_edu recoded. Old range had 9 (Don't know); new range [%g, %g], NAs=%d.\n",
            min(full_df$par_edu, na.rm = TRUE), max(full_df$par_edu, na.rm = TRUE),
            sum(is.na(full_df$par_edu))))

make_school_c1fe <- function(s) ifelse(s %in% c(201,212,213,214), 101L, as.integer(s))
firth <- function(f, d) glm(f, data = d, family = binomial, method = "brglmFit")
ORtab <- function(fit, keep = NULL) {
  ct <- summary(fit)$coefficients
  rownames(ct) <- gsub("factor\\(toa\\)","GP", gsub("factor\\(grade_period\\)","GP",
                   gsub("factor\\(school_c1fe\\)","School", rownames(ct))))
  out <- data.frame(parameter = rownames(ct), OR = round(exp(ct[,1]),3),
                    p = round(ct[,4],4), sig = ifelse(ct[,4]<.01,"**",ifelse(ct[,4]<.05,"*","")),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (!is.null(keep)) out <- out[out$parameter %in% keep, ]
  out
}

# ── (ii) THRESHOLD — Firth grouped binomial, corrected par_edu ────────────────
thr <- thr_sub[!is.na(thr_sub$n_alters) & thr_sub$n_alters > 0, ]
thr$school_c1fe <- make_school_c1fe(thr$school)
thr_covs <- c("factor(toa)","cohort","female","hispanic","asian","par_edu",
              "factor(school_c1fe)","gad","mdd","friends_ecig","sex_min")
req <- c("k_users","n_alters","cohort","female","hispanic","asian","par_edu","school","gad","mdd","friends_ecig")
thr <- thr[complete.cases(thr[, req]), ]
fthr <- glm(as.formula(paste("cbind(k_users, n_alters - k_users) ~", paste(thr_covs, collapse=" + "))),
            data = thr, family = binomial, method = "brglmFit")
thr_out <- ORtab(fthr); write.csv(thr_out, "outputs_AO/model/fully_corrected_threshold-AO.csv", row.names = FALSE)

# ── (iii) ADOPTION — Firth school-FE GLM, corrected par_edu ───────────────────
reg_data$school_c1fe <- make_school_c1fe(reg_data$school)
adopt_covs <- c("exposure","factor(grade_period)","cohort","female","hispanic","asian",
                "par_edu","factor(school_c1fe)","gad","mdd","friends_ecig","sex_min")
reqa <- c("adopt_next","exposure","grade_period","cohort","female","hispanic","asian",
          "par_edu","school","gad","mdd","friends_ecig")
rega <- reg_data[complete.cases(reg_data[, reqa]), ]
fadopt <- glm(as.formula(paste("adopt_next ~", paste(adopt_covs, collapse=" + "))),
              data = rega, family = binomial, method = "brglmFit")
adopt_out <- ORtab(fadopt); write.csv(adopt_out, "outputs_AO/model/fully_corrected_adoption-AO.csv", row.names = FALSE)

# ── (ii) SUBGROUPS 1a / 1b — Firth logistic, corrected par_edu ────────────────
base_covs <- c("cohort","female","hispanic","asian","sex_min","par_edu","gad","mdd",
               "out_degree_gp3","in_degree_gp3","friends_ecig_gp3","expo_lag",
               "friends_ecig_toa","friends_ecig_gp8","expo_gp8","sm_ecig_any_toa")
sub_fit <- function(out, stratum, excl, use_toa) {
  covs <- base_covs[!base_covs %in% excl]
  d <- full_df[full_df[[stratum]] == 1 & !is.na(full_df[[out]]), ]
  covs <- covs[vapply(covs, function(v) v %in% names(d) && length(unique(na.omit(d[[v]]))) >= 2, logical(1))]
  rhs <- c(if (use_toa) "factor(toa)", covs)
  d <- d[complete.cases(d[, c(out, intersect(rhs, names(d)))]), ]
  if (sum(d[[out]] == 1) < 10) return(NULL)
  fit <- tryCatch(firth(as.formula(paste(out, "~", paste(rhs, collapse=" + "))), d),
                  error = function(e) { message(out, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  o <- ORtab(fit); o$parameter[o$parameter == "(Intercept)"] <- NA; o[!is.na(o$parameter), ]
}
users <- list(
  c1 = sub_fit("col1_no_exposure",   "is_user", c("expo_lag","expo_gp8"), TRUE),
  c2 = sub_fit("col2_no_friends",    "is_user", c("friends_ecig_toa"), TRUE),
  c3 = sub_fit("col_both0",          "is_user", c("expo_lag","expo_gp8","friends_ecig_toa"), TRUE),
  c4 = sub_fit("col3_low_threshold", "is_user", c("expo_lag","expo_gp8"), TRUE))
nonusers <- list(
  c5 = sub_fit("col4_exposure",  "is_nonuser", c("friends_ecig_toa","expo_gp8"), FALSE),
  c6 = sub_fit("col5_friends",   "is_nonuser", c("friends_ecig_gp3","friends_ecig_gp8","friends_ecig_toa","expo_lag","expo_gp8"), FALSE),
  c7 = sub_fit("col6_both_high", "is_nonuser", c("friends_ecig_toa","expo_gp8"), FALSE))
merge_tbl <- function(lst, labs) {
  Reduce(function(a, b) merge(a, b, by = "parameter", all = TRUE),
         Map(function(x, l) { if (is.null(x)) return(NULL); names(x)[2:4] <- paste0(c("OR.","p.","sig."), l); x },
             lst, labs)[!sapply(lst, is.null)])
}
t1a <- merge_tbl(users,    c("NoExp","NoFriends","Neither","LowThr"))
t1b <- merge_tbl(nonusers, c("HighExp","HighFriends","BothHigh"))
write.csv(t1a, "outputs_AO/model/fully_corrected_table1a-AO.csv", row.names = FALSE)
write.csv(t1b, "outputs_AO/model/fully_corrected_table1b-AO.csv", row.names = FALSE)

# ── console ───────────────────────────────────────────────────────────────────
cat("\n===== FULLY CORRECTED — THRESHOLD (Firth + corrected par_edu) =====\n")
print(thr_out[thr_out$parameter %in% c("(Intercept)","GP4","GP5","GP8","cohortC2","female",
      "hispanic","asian","par_edu","gad","mdd","friends_ecig","sex_min","School104","School114"), ], row.names = FALSE)
cat("\n===== FULLY CORRECTED — ADOPTION (Firth FE + corrected par_edu) =====\n")
print(adopt_out[adopt_out$parameter %in% c("exposure","GP4","GP8","cohortC2","female","hispanic",
      "asian","par_edu","gad","mdd","friends_ecig","sex_min"), ], row.names = FALSE)
cat("\n===== FULLY CORRECTED — TABLE 1a (users) =====\n"); print(t1a, row.names = FALSE)
cat("\n===== FULLY CORRECTED — TABLE 1b (non-users) =====\n"); print(t1b, row.names = FALSE)
cat("\nSaved: outputs_AO/model/fully_corrected_{threshold,adoption,table1a,table1b}-AO.csv\n")

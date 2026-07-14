# 07-threshold-corrected-AO.R
# ─────────────────────────────────────────────────────────────────────────────
# CORRECTED version of Yuchan's grouped-binomial threshold model
# (the "Threshold" column of Tables 2–3 in 260612_threshold_report_final).
#
# ERROR CORRECTED — complete separation in factor(toa).
#   The diffnet starts at GP3, so a GP3 adopter can have NO out-alter who adopted
#   in an earlier period: k_users = 0 for 100% of GP3 adopters (structural).
#   The reference time level therefore predicts zero successes perfectly, so the
#   intercept and EVERY GP dummy diverge. This is exactly the degenerate block in
#   the report (Constant = -20.63; GP4..GP8 = 17.23, 17.88, 18.11, 18.46, 18.84).
#   Those rows are numerical artifacts and must not be interpreted.
#
# FIX (two options; we RECOMMEND the Firth fit, model C) —
#   B. Exclude GP3 adopters (toa > 3): the susceptibility model already does this. GP3
#      adopters are all at the threshold floor (k=0) and carry no information about what
#      predicts a HIGHER threshold, so dropping them is information-lossless for the slopes
#      (proof: B_intercept - A_intercept = A's toa4 coef; B's GP dummies = A's GP differences).
#      CAVEAT: a SECOND separation remains in this fit — factor(school_c1fe)105 (9 rows, all
#      k=0) — which excluding GP3 does NOT resolve.
#   C. Firth-penalized GLM (brglm2) on the FULL data (GP3 included): RECOMMENDED primary.
#      It keeps all 406 adopters, has a guaranteed-finite MLE, fixes BOTH separations (GP and
#      school105), and recovers the genuine steep grade-period trend. Substantive ORs are
#      *substantively similar* to A/B (Firth's mild shrinkage moves cohortC2/mdd/hispanic at
#      the 2nd-3rd decimal; PFU and Asian are unchanged), not bit-identical.
#   Note on n: 149 GP3 adopters are non-isolates; 44 survive complete-case deletion into the
#   fitted frame — all with k=0, so the separation holds either way.
#
# INPUT  : analytic frame `thr_sub` from Yuchan's 260612 pipeline. If not already
#          in the session, it is read from a cached RDS (set AO_THR_SUB to override).
# OUTPUT : outputs_AO/model/threshold_corrected-AO.csv  (+ console comparison)
# ─────────────────────────────────────────────────────────────────────────────

suppressWarnings(suppressMessages({}))

thr_rds <- Sys.getenv("AO_THR_SUB",
  unset = "playground/outputs_260612/AO_thr_sub.rds")
if (!exists("thr_sub")) {
  if (!file.exists(thr_rds))
    stop("thr_sub not in session and cache not found: ", thr_rds,
         "\n  Run Yuchan's 260612 pipeline first (it saves AO_thr_sub.rds).")
  thr_sub <- readRDS(thr_rds)
  message("Loaded cached thr_sub: ", thr_rds, " (", nrow(thr_sub), " adopters)")
}

make_school_c1fe <- function(s) ifelse(s %in% c(201, 212, 213, 214), 101L, as.integer(s))

# Yuchan's exact covariate block and complete-case requirement
covs <- c("factor(toa)", "cohort", "female", "hispanic", "asian", "par_edu",
          "factor(school_c1fe)", "gad", "mdd", "friends_ecig", "sex_min")
req  <- c("k_users", "n_alters", "cohort", "female", "hispanic", "asian",
          "par_edu", "school", "gad", "mdd", "friends_ecig")

prep <- function(d, drop_gp3) {
  d <- d[!is.na(d$n_alters) & d$n_alters > 0, ]        # drop isolates (as Yuchan does)
  d$school_c1fe <- make_school_c1fe(d$school)
  if (drop_gp3) d <- d[d$toa > 3, ]                    # THE FIX
  d[complete.cases(d[, req]), ]
}
fml <- as.formula(paste("cbind(k_users, n_alters - k_users) ~", paste(covs, collapse = " + ")))

tidy <- function(fit, model, n) {
  ct <- summary(fit)$coefficients
  data.frame(model = model, n = n,
             parameter = rownames(ct),
             estimate  = round(ct[, 1], 4),
             se        = round(ct[, 2], 4),
             OR        = round(exp(ct[, 1]), 4),
             p         = round(ct[, 4], 5),
             separated = abs(ct[, 1]) > 10,
             row.names = NULL, stringsAsFactors = FALSE)
}

# ── A. As-is (Yuchan): GP3 included → separated ──────────────────────────────
dA  <- prep(thr_sub, drop_gp3 = FALSE)
fA  <- glm(fml, data = dA, family = binomial)
tA  <- tidy(fA, "A_yuchan_asis", nrow(dA))

# ── B. Corrected: GP3 excluded → finite GP/intercept, same substantive ORs ───
dB  <- prep(thr_sub, drop_gp3 = TRUE)
fB  <- glm(fml, data = dB, family = binomial)
tB  <- tidy(fB, "B_AO_excl_GP3", nrow(dB))

# ── C. Robustness: Firth penalty on FULL set (brglm2), GP3 included ───────────
tC <- NULL
firth_lib <- "playground/Rlib"
if (requireNamespace("brglm2", quietly = TRUE,
                     lib.loc = c(firth_lib, .libPaths()))) {
  library(brglm2, lib.loc = c(firth_lib, .libPaths()))
  fC <- tryCatch(glm(fml, data = dA, family = binomial, method = "brglmFit"),
                 error = function(e) { message("Firth fit failed: ", conditionMessage(e)); NULL })
  if (!is.null(fC)) tC <- tidy(fC, "C_firth_full", nrow(dA))
} else message("brglm2 not available — skipping Firth robustness fit.")

allres <- do.call(rbind, Filter(Negate(is.null), list(tA, tB, tC)))
dir.create("outputs_AO/model", recursive = TRUE, showWarnings = FALSE)
write.csv(allres, "outputs_AO/model/threshold_corrected-AO.csv", row.names = FALSE)

# ── Focused console comparison ───────────────────────────────────────────────
keep <- c("(Intercept)", "GP4", "GP5", "GP8", "asian", "friends_ecig",
          "hispanic", "gad", "mdd", "cohortC2", "sex_min")
pretty <- function(t) {
  t$parameter <- gsub("factor\\(toa\\)", "GP", t$parameter)
  t[t$parameter %in% keep, c("parameter", "OR", "p", "separated")]
}
cat("\n================ THRESHOLD: as-is vs corrected ================\n")
cat(sprintf("\n--- A. Yuchan as-is (n=%d, GP3 included) ---\n", nrow(dA)))
print(pretty(tA), row.names = FALSE)
cat(sprintf("\n--- B. AO corrected (n=%d, GP3 excluded) ---\n", nrow(dB)))
print(pretty(tB), row.names = FALSE)
if (!is.null(tC)) {
  cat(sprintf("\n--- C. Firth penalized (n=%d, GP3 included, finite) ---\n", nrow(dA)))
  print(pretty(tC), row.names = FALSE)
}
cat("\nSeparated parameters (|beta|>10):\n")
cat("  A as-is     :", paste(gsub("factor\\(toa\\)","GP",tA$parameter[tA$separated]), collapse=", "), "\n")
cat("  B corrected :", paste(gsub("factor\\(toa\\)","GP",tB$parameter[tB$separated]), collapse=", "), "\n")
cat("\nSaved: outputs_AO/model/threshold_corrected-AO.csv\n")

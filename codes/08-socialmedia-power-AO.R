# 08-socialmedia-power-AO.R
# ─────────────────────────────────────────────────────────────────────────────
# CORRECTED version of Yuchan's social-media subgroup analysis
# (the `sm_ecig_any (W4–W8)` row of Table 1a in 260612_threshold_report_final).
#
# ERROR ADDRESSED — under-power / structural missingness.
#   Social-media items exist ONLY W4–W8. The 224 GP3 (= W3) adopters — the single
#   largest adoption cohort — have NO social-media measurement at adoption and are
#   silently dropped (the report note: "GP3 adopters excluded from social media
#   models"). Of 841 adopters only ~456 are measurable, and each subgroup's
#   sm model is fit on ~210–243 complete cases while carrying ~14 covariates,
#   across 4 subgroups, with no multiplicity control. The two "significant" cells
#   (No Exposure OR 2.31*, No Perceived Friends OR 0.40*) point in OPPOSITE
#   directions — a fragility signature, not a coherent channel effect.
#
# WHAT THIS SCRIPT DOES
#   1. Reproduces Yuchan's Table-1a sm_ecig_any OR per user subgroup (baseline).
#   2. Quantifies the realized N, sm prevalence, events, and the MINIMUM DETECTABLE
#      OR at 80% power given each model's own precision (are the nulls informative?).
#   3. Robustness: (a) parsimonious covariate block (restores events-per-variable);
#      (b) a cumulative "ever saw e-cig content W4..TOA" operationalization that
#      recovers late-adopter signal; (c) Benjamini–Hochberg + Bonferroni across the
#      four user-subgroup tests.
#
# INPUT  : full_df (+ timevar_df for the cumulative variant) from Yuchan's 260612
#          pipeline; read from cache if not already in the session.
# OUTPUT : outputs_AO/model/socialmedia_corrected-AO.csv
# ─────────────────────────────────────────────────────────────────────────────

cache_dir <- Sys.getenv("AO_CACHE", unset = "playground/outputs_260612")
get_obj <- function(name) {
  if (exists(name, inherits = TRUE)) return(get(name, inherits = TRUE))
  p <- file.path(cache_dir, paste0("AO_", name, ".rds"))
  if (!file.exists(p)) stop(name, " not in session and cache missing: ", p)
  readRDS(p)
}
full_df    <- get_obj("full_df")
timevar_df <- tryCatch(get_obj("timevar_df"), error = function(e) NULL)
if (!"id_orig" %in% names(full_df))
  full_df$id_orig <- rownames(full_df)            # id_orig is the merge key

z <- qnorm(c(.975, .80))                           # alpha=.05 two-sided, 80% power
min_detectable_OR <- function(se) exp(sum(z) * se) # given a coefficient SE

# user subgroups exactly as in Table 1a, with Yuchan's definitional exclusions
subgroups <- list(
  list(out = "col1_no_exposure",   label = "No Exposure",         excl = c("expo_lag","expo_gp8")),
  list(out = "col2_no_friends",    label = "No Perceived Friends", excl = c("friends_ecig_toa")),
  list(out = "col_both0",          label = "Neither",             excl = c("expo_lag","expo_gp8","friends_ecig_toa")),
  list(out = "col3_low_threshold", label = "Low Threshold",       excl = c("expo_lag","expo_gp8"))
)
base_covs <- c("cohort","female","hispanic","asian","sex_min","par_edu","gad","mdd",
               "out_degree_gp3","in_degree_gp3","friends_ecig_gp3","expo_lag",
               "friends_ecig_toa","friends_ecig_gp8","expo_gp8")  # + factor(toa) + sm var
parsimony <- c("cohort","female","par_edu")                       # too-sparse demo block
# Temporally-valid PRE-ADOPTION (baseline GP3) confounders only. The fair middle ground:
# avoids the post-treatment collider friends_ecig_toa AND the under-adjustment of `parsimony`.
baseline  <- c("cohort","female","par_edu","friends_ecig_gp3","out_degree_gp3","in_degree_gp3")

# brglm2 (Firth) from the side lib, for finite small-sample estimates
has_firth <- requireNamespace("brglm2", quietly = TRUE,
                              lib.loc = c("playground/Rlib", .libPaths()))
if (has_firth) suppressMessages(library(brglm2, lib.loc = c("playground/Rlib", .libPaths())))

fit_sm <- function(df, out, covs, sm_var, use_toa = TRUE, firth = FALSE) {
  rhs <- covs[covs %in% names(df)]
  rhs <- rhs[vapply(rhs, function(v) length(unique(na.omit(df[[v]]))) >= 2, logical(1))]
  if (use_toa) rhs <- c("factor(toa)", rhs)
  rhs <- c(rhs, sm_var)
  d <- df[df$is_user == 1 & !is.na(df[[out]]) & !is.na(df[[sm_var]]), ]
  d <- d[complete.cases(d[, c(out, intersect(rhs, names(d)))]), ]
  if (sum(d[[out]] == 1) < 10 || length(unique(na.omit(d[[sm_var]]))) < 2) return(NULL)
  f <- as.formula(paste(out, "~", paste(rhs, collapse = " + ")))
  fit <- tryCatch(
    if (firth && has_firth) glm(f, data = d, family = binomial, method = "brglmFit")
    else glm(f, data = d, family = binomial),
    error = function(e) NULL)
  if (is.null(fit) || !sm_var %in% rownames(summary(fit)$coefficients)) return(NULL)
  ct <- summary(fit)$coefficients[sm_var, ]
  list(n = nrow(d), events = sum(d[[out]] == 1),
       sm_prev = mean(d[[sm_var]] > 0, na.rm = TRUE),
       epv = sum(d[[out]] == 1) / (length(rhs) + 1),   # events-per-variable
       OR = exp(ct[1]), se = ct[2], p = ct[4],
       ci_lo = exp(ct[1] - 1.96*ct[2]), ci_hi = exp(ct[1] + 1.96*ct[2]),
       mdor = min_detectable_OR(ct[2]))
}

# ── cumulative "ever exposed by TOA" indicator (recovers late-adopter signal) ──
sm_ever <- NULL
if (!is.null(timevar_df) && all(c("id","grade_period","sm_ecig_any") %in% names(timevar_df))) {
  tv <- timevar_df[!is.na(timevar_df$sm_ecig_any), c("id","grade_period","sm_ecig_any")]
  key <- full_df[full_df$is_user == 1, c("id_orig","toa")]
  m <- merge(tv, key, by.x = "id", by.y = "id_orig")
  m <- m[m$grade_period <= m$toa, ]                # waves up to & including adoption
  sm_ever <- aggregate(sm_ecig_any ~ id, data = m, FUN = function(x) as.integer(any(x == 1)))
  names(sm_ever) <- c("id_orig", "sm_ecig_ever_toa")
  full_df <- merge(full_df, sm_ever, by = "id_orig", all.x = TRUE)
}

rows <- list()
for (g in subgroups) {
  covs_full <- base_covs[!base_covs %in% g$excl]
  full  <- fit_sm(full_df, g$out, covs_full, "sm_ecig_any_toa")               # reproduce Table 1a
  firth <- fit_sm(full_df, g$out, covs_full, "sm_ecig_any_toa", firth = TRUE) # Firth (finite, small-sample)
  base  <- fit_sm(full_df, g$out, baseline,  "sm_ecig_any_toa")               # temporally-valid baseline
  pars  <- fit_sm(full_df, g$out, parsimony, "sm_ecig_any_toa")               # too-sparse demo
  cum   <- if ("sm_ecig_ever_toa" %in% names(full_df))
             fit_sm(full_df, g$out, baseline, "sm_ecig_ever_toa") else NULL    # cumulative (note: ≠ snapshot)
  mk <- function(r, spec, sm) if (is.null(r)) NULL else data.frame(
    subgroup = g$label, spec = spec, sm_var = sm,
    n = r$n, events = r$events, sm_prev = round(r$sm_prev,3), epv = round(r$epv,1),
    OR = round(r$OR,3), ci_lo = round(r$ci_lo,3), ci_hi = round(r$ci_hi,3),
    p = round(r$p,4), min_detectable_OR_80 = round(r$mdor,2))
  rows <- c(rows, list(mk(full, "1_full_yuchan",      "sm_ecig_any_toa"),
                       mk(firth,"2_firth_full",       "sm_ecig_any_toa"),
                       mk(base, "3_baseline_GP3",     "sm_ecig_any_toa"),
                       mk(pars, "4_parsimonious",     "sm_ecig_any_toa"),
                       mk(cum,  "5_cumulative_W4toTOA","sm_ecig_ever_toa")))
}
res <- do.call(rbind, Filter(Negate(is.null), rows))

# ── multiplicity correction across the 4 user-subgroup primary tests ──────────
prim <- res[res$spec == "1_full_yuchan", ]
prim$p_BH    <- round(p.adjust(prim$p, "BH"), 4)
prim$p_bonf  <- round(p.adjust(prim$p, "bonferroni"), 4)
res <- merge(res, prim[, c("subgroup","p_BH","p_bonf")], by = "subgroup", all.x = TRUE)
res$p_BH[res$spec != "1_full_yuchan"]   <- NA
res$p_bonf[res$spec != "1_full_yuchan"] <- NA
res <- res[order(match(res$subgroup, sapply(subgroups, `[[`, "label")), res$spec), ]

dir.create("outputs_AO/model", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "outputs_AO/model/socialmedia_corrected-AO.csv", row.names = FALSE)

cat("\n=========== SOCIAL MEDIA (sm_ecig_any) — power & robustness ===========\n")
nU <- sum(full_df$is_user == 1, na.rm = TRUE)
nSnap <- sum(!is.na(full_df$sm_ecig_any_toa[full_df$is_user == 1]))
cat(sprintf("Adopters: %d. Measurable at TOA (snapshot): %d (%.0f%%).\n",
            nU, nSnap, 100*nSnap/nU))
if ("sm_ecig_ever_toa" %in% names(full_df)) {
  nCum <- sum(!is.na(full_df$sm_ecig_ever_toa[full_df$is_user == 1]))
  cat(sprintf("Measurable cumulatively (any wave W4..TOA): %d (%.0f%%) -- recovers %d adopters.\n",
              nCum, 100*nCum/nU, nCum - nSnap))
}
print(res[, c("subgroup","spec","n","events","epv","sm_prev","OR","ci_lo","ci_hi",
              "p","p_BH","p_bonf","min_detectable_OR_80")], row.names = FALSE)
cat("\nReading (see diffusion-study.md §5.2):\n")
cat(" 1 full_yuchan  : reproduces Table 1a (15-cov block + factor(toa)). p_BH/p_bonf over the 4 tests.\n")
cat(" 2 firth_full   : same model, Firth penalty -> finite small-sample estimate (over-fitting check).\n")
cat(" 3 baseline_GP3 : temporally-VALID pre-adoption confounders only (the fair middle ground).\n")
cat(" 4 parsimonious : deliberately too-sparse (cohort,female,par_edu) -> shows under-adjustment.\n")
cat(" 5 cumulative   : 'ever exposed W4..TOA'. NOTE: mechanically tied to TOA, not apples-to-apples.\n")
cat(" * min_detectable_OR_80 = post-hoc Wald MDE at the observed SE (informativeness proxy, not a-priori power).\n")
cat(" * epv = events-per-variable; <10 signals over-fitting risk.\n")
cat("Saved: outputs_AO/model/socialmedia_corrected-AO.csv\n")

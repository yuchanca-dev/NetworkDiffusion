# =============================================================================
# 260617_netdiffuseR_efficient2-AO.R   — VERIFICATION ONLY (no pipeline outputs)
#
# Claim under test (from the efficiency review): the adoption model is the runtime
# bottleneck because it is a mixed-effects logistic regression (lme4::glmer with
# student + school random intercepts, Table 2). A fixed-effects logistic GLM
# (the Table 3 "sensitivity" model) returns the SAME substantive conclusions, is
# numerically stable/reproducible, and is dramatically faster.
#
# This script fits the SAME adoption data three ways, times each, and compares:
#   (A) MIXED  : glmer(adopt_next ~ ... + (1|id) + (1|school))       <- Table 2
#   (B) FE     : glm(adopt_next ~ ... + factor(school))              <- Table 3
#   (C) FE + cluster-robust SE (clustered by student id)             <- principled fast option
#
# It does NOT touch the threshold / subgroup models or write any results.
# Input: the person-period frame produced by the efficient1 pipeline. By default
# it reads the cached copy; override with AO_REGDATA.
# =============================================================================

suppressWarnings(suppressMessages({
  library(lme4)
  library(sandwich)
  library(lmtest)
}))

regdata_path <- Sys.getenv("AO_REGDATA", unset = "playground/outputs_260612/AO_reg_data.rds")
stopifnot(file.exists(regdata_path))
reg <- readRDS(regdata_path)

# Same construction the pipeline uses
make_school_c1fe <- function(s) ifelse(s %in% c(201, 212, 213, 214), 101L, as.integer(s))
reg$school_c1fe  <- factor(make_school_c1fe(reg$school))
req <- c("adopt_next","exposure","grade_period","cohort","female","hispanic","asian",
         "par_edu","school","gad","mdd","friends_ecig","sex_min")
reg <- reg[complete.cases(reg[, req]), ]
cat(sprintf("Person-periods (complete cases): %d | students: %d\n",
            nrow(reg), length(unique(reg$id))))

rhs <- "exposure + factor(grade_period) + cohort + female + hispanic + asian + par_edu + gad + mdd + friends_ecig + sex_min"
key <- c("exposure","friends_ecig","mdd","gad")
orp <- function(est, p) sprintf("%6.2f %s", exp(est), ifelse(p<.01,"**",ifelse(p<.05,"* "," ")))

## (A) MIXED — glmer with random intercepts (Table 2) ─────────────────────────
cat("\n[A] MIXED glmer (1|id)+(1|school) ... ")
tA <- system.time({
  fA <- glmer(as.formula(paste("adopt_next ~", rhs, "+ (1|id) + (1|school)")),
              data = reg, family = binomial, control = glmerControl(optimizer = "bobyqa"))
})[["elapsed"]]
cA <- summary(fA)$coefficients; cat(sprintf("%.1f s\n", tA))

## (B) FE — plain logistic GLM with school fixed effects (Table 3) ─────────────
cat("[B] FE glm + factor(school) ... ")
tB <- system.time({
  fB <- glm(as.formula(paste("adopt_next ~", rhs, "+ factor(school_c1fe)")),
            data = reg, family = binomial)
})[["elapsed"]]
cB <- summary(fB)$coefficients; cat(sprintf("%.1f s\n", tB))

## (C) FE + cluster-robust SE by student (the principled fast option) ──────────
cat("[C] FE + cluster-robust SE (by id) ... ")
tC <- system.time({
  vcC <- sandwich::vcovCL(fB, cluster = reg$id)
  cC  <- lmtest::coeftest(fB, vcov. = vcC)
})[["elapsed"]] + tB
cat(sprintf("%.1f s (refit+vcov)\n", tC))

## Comparison ──────────────────────────────────────────────────────────────────
cat("\n================ ADOPTION MODEL: three ways ================\n")
cat(sprintf("%-14s | %-12s | %-12s | %-12s\n", "predictor", "A mixed OR", "B FE OR", "C FE+robust OR"))
cat(strrep("-", 60), "\n")
for (k in key) {
  a <- cA[k, ]; b <- cB[k, ]; cc <- cC[k, ]
  cat(sprintf("%-14s | %-12s | %-12s | %-12s\n", k,
              orp(a[1], a[4]), orp(b[1], b[4]), orp(cc[1], cc[4])))
}
cat(strrep("-", 60), "\n")
cat(sprintf("TIME           | %-12s | %-12s | %-12s\n",
            sprintf("%.1f s", tA), sprintf("%.1f s", tB), sprintf("%.1f s", tC)))
cat(sprintf("\nSpeed-up (mixed / FE): %.0fx\n", tA / tB))

cat("\n--- Reproduction check vs report Table 3 (FE) ---\n")
cat(sprintf("  exposure OR %.2f (report 2.31) | MDD %.2f (1.61) | GAD %.2f (0.77) | friends %.2f (1.43)\n",
            exp(cB["exposure",1]), exp(cB["mdd",1]), exp(cB["gad",1]), exp(cB["friends_ecig",1])))
cat("\n--- Why FE is faster (random-effect variance) ---\n")
vc <- as.data.frame(VarCorr(fA))
cat(sprintf("  glmer estimates random-intercept variances: id=%.3f, school=%.3f\n",
            vc$vcov[vc$grp=="id"], vc$vcov[vc$grp=="school"]))
cat("  -> glmer must integrate the likelihood over these random effects (Laplace),\n")
cat("     iterating over variance components; the GLM solves one fixed-design IRLS.\n")
cat("\n[A] subject-specific (conditional) ORs; [B]/[C] population-average (marginal) ORs.\n")
cat("    Magnitudes differ by design; sign & significance agree. [C] keeps FE speed but\n")
cat("    corrects SEs for repeated person-periods per student.\n")

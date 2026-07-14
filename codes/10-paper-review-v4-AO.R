# =============================================================================
# 10-paper-review-v4-AO.R
# Two deliverables for the review of the "Low Threshold Adopters V4" draft — the
# items referenced in the 2026-07-13 message to Y. Cao:
#   (1) PFU sensitivity: the four current-TOA social-media models re-fit with
#       perceived friends' e-cig use (friends_ecig, at TOA) added as a covariate.
#       PFU is the strongest predictor in all four and the social-media OR
#       attenuates. (NB: PFU is contemporaneous, so it may be a mediator rather
#       than a confounder — reported as a sensitivity analysis, not a refutation.)
#   (2) In-/out-degree histograms at TOA for the low-threshold (n=231) and
#       no-network-exposure (n=173) adopters, on a shared 0-10 axis.
#
# INPUT: the workspace snapshot of Yuchan's 260704 pipeline (its objects: expo,
#   toa_vals, thr_vals, vertex_ids_clean, cov_df, attrs_df, timevar_df, all_gps,
#   make_school_c1fe). It is produced by running that pipeline with a save.image()
#   at the end. Override the path with AO_WORKSPACE.
# OUTPUT: outputs_AO/paper_review_v4/
# =============================================================================
suppressMessages({ library(netdiffuseR); library(sandwich); library(lmtest) })

ws <- Sys.getenv("AO_WORKSPACE", unset = "playground/outputs_260704/AO_workspace_260704.RData")
if (!file.exists(ws)) stop("Workspace snapshot not found: ", ws,
  "\n  Produce it by running Yuchan's 260704 pipeline with a save.image() at the end.")
load(ws)
out_dir <- "outputs_AO/paper_review_v4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Rebuild the person-period analytic frame (as run_tom_sm_event_reg does) ────
toa_vec <- as.integer(toa_vals);  names(toa_vec) <- vertex_ids_clean
thr_vec <- as.numeric(thr_vals);  names(thr_vec) <- vertex_ids_clean
pp <- do.call(rbind, lapply(seq_along(all_gps), function(gi) data.frame(
  id = vertex_ids_clean, grade_period = all_gps[gi],
  exposure = as.numeric(expo[, gi]),
  toa = as.integer(toa_vec[vertex_ids_clean]),
  threshold = as.numeric(thr_vec[vertex_ids_clean]),
  stringsAsFactors = FALSE)))
pp <- pp[is.na(pp$toa) | pp$grade_period <= pp$toa, ]
pp <- merge(pp, cov_df, by = "id", all.x = TRUE)
pp <- merge(pp, attrs_df[, c("id","cohort","sch_type","school")], by = "id", all.x = TRUE)
pp <- merge(pp, timevar_df, by = c("id","grade_period"), all.x = TRUE)
pp$school_c1fe <- make_school_c1fe(pp$school)

base_covs <- c("cohort","female","hispanic","asian","par_edu","mdd",
               "out_degree","in_degree","sex_min")
def_noexp  <- function(d) as.integer(!is.na(d$toa) & d$grade_period == d$toa &
                                     !is.na(d$exposure)  & d$exposure == 0)
def_lowthr <- function(d) as.integer(!is.na(d$toa) & d$grade_period == d$toa &
                                     !is.na(d$threshold) & d$threshold <= 0.25)

# ── (1) PFU sensitivity ───────────────────────────────────────────────────────
fit_sm <- function(outcome_def, sm_var, add_pfu) {
  d <- pp; d$outcome <- outcome_def(d)
  covs <- c(base_covs, sm_var, if (add_pfu) "friends_ecig")
  req  <- c("outcome","grade_period","school", covs)
  d <- d[complete.cases(d[, req]), ]
  fit <- glm(as.formula(paste("outcome ~ factor(grade_period) +",
              paste(covs, collapse=" + "), "+ factor(school_c1fe)")),
             data = d, family = binomial)
  ctr <- lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = d$id))
  g <- function(v) if (v %in% rownames(ctr)) ctr[v, ] else rep(NA, 4)
  s <- g(sm_var); p <- if (add_pfu) g("friends_ecig") else rep(NA, 4)
  data.frame(
    n = nrow(d), events = sum(d$outcome == 1),
    sm_OR = round(exp(s[1]), 3),
    sm_CI = if (is.na(s[1])) NA else sprintf("[%.3f, %.3f]", exp(s[1]-1.96*s[2]), exp(s[1]+1.96*s[2])),
    sm_p  = round(s[4], 4),
    pfu_OR = round(exp(p[1]), 3), pfu_p = round(p[4], 4),
    row.names = NULL, stringsAsFactors = FALSE)
}
specs <- list(
  list("No-network-exposure", "overall (sm_post_avg_rev)", def_noexp,  "sm_post_avg_rev"),
  list("Low-threshold",       "overall (sm_post_avg_rev)", def_lowthr, "sm_post_avg_rev"),
  list("No-network-exposure", "Instagram (insta_rev)",     def_noexp,  "insta_rev"),
  list("Low-threshold",       "Instagram (insta_rev)",     def_lowthr, "insta_rev"))
pfu_tab <- do.call(rbind, lapply(specs, function(s) {
  a <- fit_sm(s[[3]], s[[4]], FALSE); b <- fit_sm(s[[3]], s[[4]], TRUE)
  rbind(cbind(outcome = s[[1]], sm = s[[2]], spec = "paper (no PFU)", a),
        cbind(outcome = s[[1]], sm = s[[2]], spec = "+ PFU",          b))
}))
cat("===== (1) PFU SENSITIVITY (current-TOA models) =====\n")
print(pfu_tab, row.names = FALSE)
write.csv(pfu_tab, file.path(out_dir, "pfu_sensitivity-AO.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "pfu_sensitivity-AO.csv"), "\n\n")

# ── (2) degree histograms (shared 0-10 axis, integer ticks) ───────────────────
req0 <- c("grade_period","cohort","female","hispanic","asian","par_edu","school",
          "mdd","out_degree","in_degree","sm_post_avg_rev","sex_min")
an <- pp[complete.cases(pp[, req0]), ]
lt <- an[def_lowthr(an) == 1, ]
ne <- an[def_noexp(an)  == 1, ]

XMAX <- 10
draw_pair <- function(d, group_lab) {
  od <- pmin(d$out_degree, XMAX); id <- pmin(d$in_degree, XMAX)
  ymax <- max(table(factor(od, 0:XMAX)), table(factor(id, 0:XMAX)))
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 3, 0))
  for (side in 1:2) {
    v   <- if (side == 1) od else id
    col <- if (side == 1) "#e76f51" else "#3b82c4"
    ttl <- if (side == 1) "Out-degree at time of adoption" else "In-degree at time of adoption"
    xlb <- if (side == 1) "Out-degree (friends nominated)" else "In-degree (nominations received)"
    hist(v, breaks = seq(-0.5, XMAX + 0.5, 1), col = col, border = "white",
         main = ttl, xlab = xlb, ylab = "Number of adopters",
         xaxt = "n", xlim = c(-0.5, XMAX + 0.5), ylim = c(0, ymax * 1.05))
    axis(1, at = 0:XMAX, labels = c(0:(XMAX - 1), "10+"), cex.axis = 0.85, gap.axis = -1)
    abline(v = mean(v), lty = 2)
  }
  mtext(group_lab, side = 3, line = 1, outer = TRUE, cex = 1.05, font = 2)
}
emit <- function(base, d, lab) {
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(file.path(out_dir, paste0(base, ".pdf")), width = 9, height = 4.6)
    else png(file.path(out_dir, paste0(base, ".png")), width = 9, height = 4.6, units = "in", res = 200)
    draw_pair(d, lab); dev.off()
  }
  cat("Saved:", file.path(out_dir, paste0(base, ".{pdf,png}")), "\n")
}
emit("lowthreshold_degree_histograms-AO", lt,
     sprintf("Low-threshold adopters (threshold <= 0.25) in the analytic sample, n = %d", nrow(lt)))
emit("noexposure_degree_histograms-AO", ne,
     sprintf("No-network-exposure adopters in the analytic sample, n = %d", nrow(ne)))

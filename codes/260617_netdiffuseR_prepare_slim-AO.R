# =============================================================================
# 260617_netdiffuseR_prepare_slim-AO.R   — ONE-TIME DATA PREP
#
# Reads the two big ADVANCE workbooks (W1-W8 ≈ 153 MB / 10,341 cols; W9-W10) ONCE,
# keeps only the ~865 columns the threshold pipeline actually uses, and writes a
# compact slim cache (~468 KB RDS). Run this once wherever the workbooks live.
# Thereafter 260617_netdiffuseR_efficient2_fixed-AO.R reads ONLY this RDS, so you
# can hand off just the RDS + that script to reproduce every result without the
# 153 MB workbooks.
#
#   Rscript codes/260617_netdiffuseR_prepare_slim-AO.R
#
# Paths (override via env vars):
#   ADVANCE_DATA : folder holding ADVANCE_W1-W8.xlsx + ADVANCE_W9-W10.xlsx (default Cleaned_Data)
#   AO_SLIM      : output RDS path (default <ADVANCE_DATA>/ADVANCE_slim_v1.rds)
# =============================================================================
library(openxlsx)

in_path   <- Sys.getenv("ADVANCE_DATA", unset = "Cleaned_Data")
slim_path <- Sys.getenv("AO_SLIM",      unset = file.path(in_path, "ADVANCE_slim_v1.rds"))

# Allow-list of columns the pipeline reads (names are lower-cased first). MUST stay
# in sync with the variables used by 260617_netdiffuseR_efficient2_fixed-AO.R; it
# also keeps the survey-susceptibility items from revised13 in case they return.
keep_patterns <- paste(c(
  "^record_id$",
  "^w[0-9]+_schoolid$",
  "^w[0-9]+_past_6mo_use_3$",                                   # adoption outcome
  "^w[0-9]+_rcads_gad_mean$", "^w[0-9]+_rcads_mdd_mean$",       # mental health
  "^w[0-9]+_friends_use_ecig$",                                # perceived friend use
  "^w[0-9]+_try_friend_ecig$", "^w[0-9]+_use_next_yr_ecig$",    # survey susceptibility (revised13 reuse)
  "^w[0-9]+_dem_gender$", "^w[0-9]+_eth$", "^w[0-9]+_race$", "^w[0-9]+_dem_sexuality$",
  "^w[0-9]+_dem_high_par_edu(_new)?$",                         # parent education (+ new-scale variant)
  "^w[0-9]+_sm_post_[0-9]+$", "^w[0-9]+_ecig_posted_[0-9a-z]+$", # social media
  "^w[0-9]+_friend[0-9]+_[0-9]+$"                              # friendship-nomination edges
), collapse = "|")

f18  <- file.path(in_path, "ADVANCE_W1-W8.xlsx")
f910 <- file.path(in_path, "ADVANCE_W9-W10.xlsx")
stopifnot(file.exists(f18), file.exists(f910))

cat("Reading full ADVANCE workbooks (slow, one-time)...\n")
t0     <- Sys.time()
raw18  <- read.xlsx(f18,  sheet = 1)
raw910 <- read.xlsx(f910, sheet = 1)
cat(sprintf("  read in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

slim_keep <- function(d) {
  names(d) <- tolower(names(d))
  d[, grepl(keep_patterns, names(d)), drop = FALSE]
}
d_w18  <- slim_keep(raw18)
d_w910 <- slim_keep(raw910)
rm(raw18, raw910); gc()

saveRDS(list(w18 = d_w18, w910 = d_w910), slim_path, compress = "xz")
cat(sprintf("\nSaved slim cache: %s  (%.0f KB)\n", slim_path, file.size(slim_path) / 1024))
cat(sprintf("  W1-W8:  %d students x %d columns\n",  nrow(d_w18),  ncol(d_w18)))
cat(sprintf("  W9-W10: %d students x %d columns\n",  nrow(d_w910), ncol(d_w910)))
cat("\nNow run:  Rscript codes/260617_netdiffuseR_efficient2_fixed-AO.R\n")

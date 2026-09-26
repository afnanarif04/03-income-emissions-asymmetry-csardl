# =====================================================================
# 08_temporal_stability_check.R
# NOT a manuscript table. This script reproduces a diagnostic robustness
# check discussed only in the Conclusion's Limitations paragraph: whether
# the symmetric long-run income and energy-use elasticities are stable
# across two long sub-windows (1990-2006 vs 2007-2022).
#
# An asymmetric or Fourier-augmented estimator could not be estimated on
# either sub-window (T = 16-17 years each): every primary estimator
# tested (CS-ARDL, CS-NARDL, PMG-ARDL) produced economically implausible
# coefficients (population elasticities beyond -14, negative-shock
# coefficients beyond -47, error-correction speeds exceeding 1 in
# absolute value, and CS-ARDL failed to identify at all for 2007-2022).
# Only the symmetric CS-DL specification, which has no lagged dependent
# variable and one long-run coefficient rather than two, remains
# estimable at this sample length. This script is retained for
# transparency and reproducibility; its output is not reported as a
# manuscript table.
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")

stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                        ifelse(p < 0.10, "*", ""))))
}

split_estimate <- function(yr_lo, yr_hi, label) {
  p <- filter(panel, year >= yr_lo & year <= yr_hi)
  p <- p %>% group_by(id) %>% filter(n() >= 7) %>% ungroup()
  res <- cs_dl(p, KEY, CTRL, pT = 0, asym = FALSE, fourier = FALSE)
  if (nrow(res) == 0) return(tibble())
  res %>% filter(term %in% c("lgdp","lenergy")) %>%
    mutate(sample = label, stars = stars(p))
}

split_res <- bind_rows(
  split_estimate(1990, 2006, "1990-2006"),
  split_estimate(2007, 2022, "2007-2022")
)
cat("\n===== TEMPORAL STABILITY (diagnostic only, not a manuscript table) =====\n")
print(as.data.frame(
  split_res %>% mutate(cell = sprintf("%.3f%s", coef, stars)) %>%
    select(sample, term, cell) %>%
    pivot_wider(names_from = term, values_from = cell)), row.names = FALSE)

write_csv(split_res, "output/supplementary_temporal_stability_not_a_table.csv")

cat("\n08_temporal_stability_check.R complete\n\n")

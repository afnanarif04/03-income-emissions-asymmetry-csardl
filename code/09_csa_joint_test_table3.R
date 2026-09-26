# =====================================================================
# 09_csa_joint_test_table3.R
# Produces output/table3e_csa_joint_f_main_panels.csv (feeds the "CSA
# terms" row of manuscript Table 3, Panels A and B, for the full panel
# and the three carbon lock-in sub-panels).
#
# Joint significance (F-test) of the cross-sectional average (CSA) terms
# in each unit-level regression, mean-grouped across units. Reports
# whether the CSA augmentation is jointly significant -- i.e. whether
# common-factor proxying is actually doing work in the regression -- for
# each primary/supplementary estimator and each sample.
# H0: all CSA coefficients = 0.
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")
KEY <- "lgdp"; CTRL <- c("lenergy", "lpop")

# ---------------------------------------------------------------------
# Mean-group the unit F-statistics for one (sample, estimator) cell.
# ---------------------------------------------------------------------
mg_csa_f <- function(panel, sub, method) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  # CS-ARDL and the CS-DL family are LINEAR in the key regressor, so the
  # level-based CSA is the correct nonlinear-transform-consistent proxy
  # for them (Kapetanios, Mitchell & Shin 2014 only bind for a nonlinear
  # conditional mean). CS-NARDL/FCS-NARDL enter the regression through
  # key_pos/key_neg, so their CSA must be built from those partial sums
  # -- see the corrected cs_nardl() in 03_estimators.R for the full
  # citation and reasoning. pT is 1 for the linear methods (unchanged)
  # and 0 for the nonlinear methods (verified in 12_csnardl_verified.R:
  # the extra CSA variable from the partial-sum split pushes |ECT| above
  # 1 at pT=1 on this panel's T=33).
  csa_level <- make_csa(p, c("lco2", KEY, CTRL))
  p_ps <- p %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[KEY]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[KEY]])), 0))) %>%
    ungroup()
  csa_ps <- make_csa(p_ps, c("lco2", "key_pos", "key_neg", CTRL))
  ids <- unique(p$id)
  fs <- map_dbl(ids, function(i) {
    d <- filter(p, id == i)
    d_ps <- filter(p_ps, id == i)
    switch(method,
      "CS-ARDL"   = csa_joint_f_ardl(d, csa_level, KEY, CTRL, pT = 1),
      "CS-NARDL"  = csa_joint_f_nardl(d_ps, csa_ps, KEY, CTRL, pT = 0, fourier = FALSE),
      "FCS-NARDL" = csa_joint_f_nardl(d_ps, csa_ps, KEY, CTRL, pT = 0, fourier = TRUE),
      "CS-DL"     = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = FALSE, fourier = FALSE),
      "CS-NDL"    = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = TRUE,  fourier = FALSE),
      "FCS-NDL"   = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = TRUE,  fourier = TRUE),
      NA_real_)
  })
  mean(fs, na.rm = TRUE)
}

samples <- list(All = NULL, coal = "coal", gas = "gas", oil = "oil")
methods <- c("CS-ARDL", "CS-NARDL", "FCS-NARDL", "CS-DL", "CS-NDL", "FCS-NDL")

csa_f_table <- map_dfr(names(samples), function(sn) {
  map_dfr(methods, function(m) {
    tibble(sample = if (is.null(samples[[sn]])) "All countries" else
                    paste0(str_to_title(samples[[sn]]), "-dominant"),
           method = m,
           csa_joint_F = mg_csa_f(panel, samples[[sn]], m))
  })
})

cat("Mean-group CSA joint-significance F-test (H0: all CSA coefficients = 0):\n")
print(as.data.frame(csa_f_table))

write_csv(csa_f_table, "output/table3e_csa_joint_f_main_panels.csv")
cat("\n09_csa_joint_test_table3.R complete\n\n")

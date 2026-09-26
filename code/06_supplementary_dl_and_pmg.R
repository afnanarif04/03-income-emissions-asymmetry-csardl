# =====================================================================
# 06_supplementary_dl_and_pmg.R
# Produces output/table3_panelB_supplementary_dl.csv and
# output/table3_panelC_pmg.csv (manuscript Table 3, Panels B and C, ALL
# FOUR samples: All countries, Coal-dominant, Gas-dominant, Oil-dominant).
#
# FIX (this version): Panel B previously ran cs_dl() on the full,
# unfiltered panel only, leaving the Coal-/Gas-/Oil-dominant long-run
# columns of Table 3 Panel B blank, while the "CSA terms" row for those
# same columns (fed separately by 09_csa_joint_test_table3.R, which does
# loop over sub-panels) was populated -- an internal inconsistency
# spotted on manuscript review. This version loops over the same four
# samples as Panel A (05_main_estimation_table3.R) and Panel C's PMG
# block below, so every cell in Panel B is now populated from the same
# reproducible pipeline. Verified in sandbox: cs_dl() runs cleanly on
# all three carbon lock-in sub-panels with no estimation failures --
# the missing cells were a script gap, not a genuine non-identification
# case.
#
# None of the estimators here is presented as a contribution of this
# study:
#   (a) distributed-lag counterparts CS-DL / CS-NDL / FCS-NDL that drop
#       the lagged dependent variable (a specification robustness check);
#   (b) PMG pooled long-run estimates (Pesaran, Shin & Smith 1999
#       concentrated maximum likelihood; a slope-homogeneity benchmark).
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")

KEY <- "lgdp"; CTRL <- c("lenergy", "lpop")

stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                        ifelse(p < 0.10, "*", ""))))
}

samples <- list(All = NULL, coal = "coal", gas = "gas", oil = "oil")

# ---- long-run symmetry Wald test for the CS-NDL/FCS-NDL group ---------
# Mirrors symmetry_wald() in 05_main_estimation_table3.R exactly, but
# calls cs_dl_unit(asym=TRUE) instead of cs_nardl_unit2(). IMPORTANT: the
# CS-DL family's cross-sectional average is built from the LEVEL of the
# key regressor, not from its partial sums, because CS-ARDL and CS-DL
# are both linear in the key regressor (the Kapetanios, Mitchell & Shin
# 2014 nonlinear-transform-consistent CSA only binds for the nonlinear
# CS-NARDL conditional mean) -- see the note at the top of
# 09_csa_joint_test_table3.R for the same point applied to the CSA
# joint-F test. Using the partial-sum CSA here (copying
# symmetry_wald()'s CS-NARDL convention verbatim) would silently proxy
# for the wrong common factor and was avoided.
symmetry_wald_dl <- function(panel, sub = NULL, fourier = FALSE) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  csa <- make_csa(p, c("lco2", KEY, CTRL))
  ids <- unique(p$id)
  diffs <- map_dbl(ids, function(i) {
    d <- filter(p, id == i)
    u <- cs_dl_unit(d, csa, KEY, CTRL, pT = 1, asym = TRUE, fourier = fourier)
    if (is.null(u)) return(NA_real_)
    as.numeric(u["key_pos"] - u["key_neg"])
  })
  diffs <- diffs[!is.na(diffs)]
  N <- length(diffs); mdiff <- mean(diffs); se <- sd(diffs) / sqrt(N)
  z <- mdiff / se
  tibble(diff_pos_neg = mdiff, se = se, z = z, wald_F = z^2,
         p_value = 2 * (1 - pnorm(abs(z))), N = N)
}

# ---- (a) supplementary distributed-lag checks, all four samples ------
run_dl_block <- function(panel, sub) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  lbl <- if (is.null(sub)) "All countries" else paste0(str_to_title(sub), "-dominant")
  bind_rows(
    cs_dl(p, KEY, CTRL, asym = FALSE, fourier = FALSE) %>% mutate(method = "CS-DL"),
    cs_dl(p, KEY, CTRL, asym = TRUE,  fourier = FALSE) %>% mutate(method = "CS-NDL"),
    cs_dl(p, KEY, CTRL, asym = TRUE,  fourier = TRUE)  %>% mutate(method = "FCS-NDL")
  ) %>% mutate(sample = lbl, stars = stars(p))
}

dl_res <- map_dfr(names(samples), function(sn) run_dl_block(panel, samples[[sn]]))

cat("\n===== SUPPLEMENTARY distributed-lag checks (long-run), all four samples =====\n")
print(as.data.frame(
  dl_res %>% filter(term %in% c("lgdp","key_pos","key_neg","lenergy","lpop")) %>%
    mutate(cell = sprintf("%.3f%s", coef, stars)) %>%
    select(sample, method, term, cell) %>%
    pivot_wider(names_from = term, values_from = cell)), row.names = FALSE)
write_csv(dl_res, "output/table3_panelB_supplementary_dl.csv")

# ---- Wald symmetry test for CS-NDL and FCS-NDL, all four samples -----
dl_wald <- map_dfr(names(samples), function(sn) {
  w  <- symmetry_wald_dl(panel, samples[[sn]], fourier = FALSE)
  wf <- symmetry_wald_dl(panel, samples[[sn]], fourier = TRUE)
  tibble(sample = if (is.null(samples[[sn]])) "All countries" else
                  paste0(str_to_title(samples[[sn]]), "-dominant"),
         wald_diff_CSNDL = w$diff_pos_neg, wald_F_CSNDL = w$wald_F, wald_p_CSNDL = w$p_value,
         wald_diff_FCSNDL = wf$diff_pos_neg, wald_F_FCSNDL = wf$wald_F, wald_p_FCSNDL = wf$p_value)
})
cat("\n===== SUPPLEMENTARY distributed-lag checks: LONG-RUN SYMMETRY WALD TEST =====\n")
print(as.data.frame(dl_wald), row.names = FALSE)
write_csv(dl_wald, "output/table3_panelB_wald_symmetry.csv")

# ---- (b) PMG homogeneity benchmark, all four samples (unchanged from
# the prior version -- this block already looped over samples) --------
pmg_by_panel <- map_dfr(names(samples), function(sn) {
  p <- if (is.null(samples[[sn]])) panel else filter(panel, lockin == samples[[sn]])
  lbl <- if (is.null(samples[[sn]])) "All countries" else
         paste0(str_to_title(samples[[sn]]), "-dominant")
  bind_rows(
    pmg_est(p, KEY, CTRL, asym = FALSE) %>% mutate(method = "PMG-ARDL"),
    pmg_est(p, KEY, CTRL, asym = TRUE)  %>% mutate(method = "PMG-NARDL"),
    pmg_est(p, KEY, CTRL, asym = TRUE, fourier = TRUE) %>% mutate(method = "FPMG-NARDL")
  ) %>% mutate(sample = lbl)
})
cat("\n===== PMG pooled long-run benchmark (concentrated ML), all four samples =====\n")
print(as.data.frame(pmg_by_panel), row.names = FALSE)
write_csv(pmg_by_panel, "output/table3_panelC_pmg.csv")

cat("\n06_supplementary_dl_and_pmg.R complete\n\n")

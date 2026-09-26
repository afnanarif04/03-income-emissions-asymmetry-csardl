# =====================================================================
# 07_subsample_sensitivity_table4.R
# Produces output/table4_income_terciles.csv (manuscript Table 4).
#
# STRATEGY A (verified in sandbox trial, session of 2026-07-05): the
# OECD sub-sample used in the previous version of this script has been
# replaced with three income-level terciles -- low, middle, and
# high-income economies -- formed by splitting all 153 countries into
# three equal-sized groups (51/51/51) on time-averaged GDP per capita
# (current USD). This directly tests whether the paper's headline
# finding (long-run symmetry not rejected; energy use dominates) holds
# across the development spectrum, mirroring the income-level framing
# of Sim & Sek (2025, Energy & Environment) and Namahoro et al. (2021,
# Energy 237), the paper's two locked interpretation benchmarks.
#
# NOTE ON SCOPE: the OECD sub-sample and the "coal-dependent economies"
# sub-sample used in earlier drafts are both dropped. OECD membership
# is an institutional/policy grouping, not a development-level grouping,
# and duplicated most of the same high-income countries already
# implicitly covered by the new high-income tercile. The coal-dependent
# sub-sample was already dropped in the prior version because it
# duplicates Table 3's coal-dominant carbon lock-in panel exactly (same
# 10 countries; see data/README_data.md for the classification).
#
# Runs all nine estimators (primary CS-ARDL group, supplementary CS-DL
# group, and PMG group) on each income tercile, together with the
# cross-sectional-average joint-significance F-test AND the long-run
# symmetry Wald test (the latter is new to Table 4 in this version --
# it was previously Table-3-only, but the income-tercile split makes it
# directly relevant here: the trial run found symmetry rejected for the
# middle-income tercile at CS-NARDL, p = 0.026, and marginally at
# FCS-NARDL, p = 0.064, while low- and high-income terciles do not
# reject symmetry. This is a genuine, reportable heterogeneity result).
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")

KEY  <- "lgdp"
CTRL <- c("lenergy", "lpop")

stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                        ifelse(p < 0.10, "*", ""))))
}

# ---- build income terciles from time-averaged GDP per capita ----------
# gdp_per_capita (current USD, level) is already in Panel_Data; not used
# by any estimator directly (only its natural log, lgdp, is), so this
# grouping variable is new but requires no new data collection.
avg_gdp <- panel %>%
  group_by(iso3) %>%
  summarise(avg_gdp_pc = mean(gdp_per_capita, na.rm = TRUE), .groups = "drop")

income_cuts <- quantile(avg_gdp$avg_gdp_pc, probs = c(1/3, 2/3), na.rm = TRUE)

avg_gdp <- avg_gdp %>%
  mutate(income_group = dplyr::case_when(
    avg_gdp_pc <= income_cuts[1] ~ "low",
    avg_gdp_pc <= income_cuts[2] ~ "middle",
    TRUE ~ "high"
  ))

cat("Income tercile cutoffs (avg GDP per capita, current USD):\n")
print(income_cuts)
cat("\nIncome group counts:\n")
print(table(avg_gdp$income_group))

low_iso3    <- avg_gdp %>% filter(income_group == "low")    %>% pull(iso3)
middle_iso3 <- avg_gdp %>% filter(income_group == "middle") %>% pull(iso3)
high_iso3   <- avg_gdp %>% filter(income_group == "high")   %>% pull(iso3)

# ---- run every estimator (primary + supplementary DL + PMG) on one
# income-tercile sub-sample, tidied to the same term set used throughout
# the study (unchanged from the prior OECD version of this function) ---
sub_estimate_full <- function(iso_set, label) {
  p <- filter(panel, iso3 %in% iso_set)
  if (n_distinct(p$id) < 5) return(NULL)

  primary <- bind_rows(
    cs_ardl(p, KEY, CTRL, pT = 1) %>% mutate(method = "CS-ARDL"),
    cs_nardl(p, KEY, CTRL, pT = 0, fourier = FALSE) %>% mutate(method = "CS-NARDL"),
    cs_nardl(p, KEY, CTRL, pT = 0, fourier = TRUE)  %>% mutate(method = "FCS-NARDL")
  ) %>% mutate(stars = stars(p))
  suppl <- bind_rows(
    cs_dl(p, KEY, CTRL, pT = 1, asym = FALSE, fourier = FALSE) %>% mutate(method = "CS-DL"),
    cs_dl(p, KEY, CTRL, pT = 1, asym = TRUE,  fourier = FALSE) %>% mutate(method = "CS-NDL"),
    cs_dl(p, KEY, CTRL, pT = 1, asym = TRUE,  fourier = TRUE)  %>% mutate(method = "FCS-NDL")
  ) %>% mutate(stars = stars(p))
  # PMG's stars column is already computed inside pmg_est() from its own
  # Hessian-based t-statistics (see 02_estimator_functions.R); it must NOT
  # be overwritten here by stars(p), since PMG's tibble has no p column of
  # its own (mutate(stars = stars(p)) on PMG rows would silently null out
  # every PMG significance star, as happened before this fix -- confirmed
  # by direct inspection: PMG's se/coef were correct but its stars were
  # entirely NaN in the written CSV until this line was corrected).
  pmg_sym  <- pmg_est(p, KEY, CTRL, asym = FALSE) %>% mutate(method = "PMG-ARDL")
  pmg_asym <- pmg_est(p, KEY, CTRL, asym = TRUE)  %>% mutate(method = "PMG-NARDL")
  pmg_fourier <- pmg_est(p, KEY, CTRL, asym = TRUE, fourier = TRUE) %>%
    mutate(method = "FPMG-NARDL")

  bind_rows(primary, suppl, pmg_sym, pmg_asym, pmg_fourier) %>%
    mutate(sample = label)
}

# ---- CSA joint-significance F-test for the same sub-sample (unchanged
# from the prior OECD version of this function) -------------------------
sub_csa_f <- function(iso_set, label) {
  p <- filter(panel, iso3 %in% iso_set)
  if (n_distinct(p$id) < 5) return(tibble())
  csa_level <- make_csa(p, c("lco2", KEY, CTRL))
  p_ps <- p %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[KEY]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[KEY]])), 0))) %>%
    ungroup()
  csa_ps <- make_csa(p_ps, c("lco2", "key_pos", "key_neg", CTRL))
  ids <- unique(p$id)
  methods6 <- c("CS-ARDL", "CS-NARDL", "FCS-NARDL", "CS-DL", "CS-NDL", "FCS-NDL")
  map_dfr(methods6, function(m) {
    fs <- map_dbl(ids, function(i) {
      d <- filter(p, id == i)
      d_ps <- filter(p_ps, id == i)
      switch(m,
        "CS-ARDL"   = csa_joint_f_ardl(d, csa_level, KEY, CTRL, pT = 1),
        "CS-NARDL"  = csa_joint_f_nardl(d_ps, csa_ps, KEY, CTRL, pT = 0, fourier = FALSE),
        "FCS-NARDL" = csa_joint_f_nardl(d_ps, csa_ps, KEY, CTRL, pT = 0, fourier = TRUE),
        "CS-DL"     = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = FALSE, fourier = FALSE),
        "CS-NDL"    = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = TRUE,  fourier = FALSE),
        "FCS-NDL"   = csa_joint_f_dl(d, csa_level, KEY, CTRL, pT = 1, asym = TRUE,  fourier = TRUE),
        NA_real_)
    })
    tibble(sample = label, method = m, csa_joint_F = mean(fs, na.rm = TRUE))
  })
}

# ---- long-run symmetry Wald test, adapted from symmetry_wald() in
# 05_main_estimation_table3.R -- NEW to Table 4 in this version, since
# income level (unlike OECD membership) is the exact dimension the
# paper's interpretation benchmarks (Sim & Sek 2025; Namahoro et al.
# 2021) test symmetry against. Identical logic to Table 3's version,
# filtered by iso3 set instead of the lockin label. ----------------------
symmetry_wald_iso <- function(iso_set, fourier = FALSE) {
  p <- filter(panel, iso3 %in% iso_set)
  p_ps <- p %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[KEY]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[KEY]])), 0))) %>%
    ungroup()
  csa <- make_csa(p_ps, c("lco2", "key_pos", "key_neg", CTRL))
  ids <- unique(p_ps$id)
  diffs <- map_dbl(ids, function(i) {
    d <- filter(p_ps, id == i); if (fourier) d <- attach_fourier(d)
    u <- cs_nardl_unit2(d, csa, KEY, CTRL, pT = 0, fourier = fourier)
    if (is.null(u)) return(NA_real_)
    as.numeric(u["key_pos"] - u["key_neg"])
  })
  diffs <- diffs[!is.na(diffs)]
  N <- length(diffs); mdiff <- mean(diffs); se <- sd(diffs) / sqrt(N)
  z <- mdiff / se
  tibble(diff_pos_neg = mdiff, se = se, z = z, wald_F = z^2,
         p_value = 2 * (1 - pnorm(abs(z))), N = N)
}

groups <- list(
  list(iso = low_iso3,    lbl = "Low-income economies"),
  list(iso = middle_iso3, lbl = "Middle-income economies"),
  list(iso = high_iso3,   lbl = "High-income economies")
)

# ---- long-run symmetry Wald test for the CS-NDL/FCS-NDL group ---------
# Mirrors symmetry_wald_dl() in 06_supplementary_dl_and_pmg.R, filtered
# by iso3 set instead of the lockin label. As in that function, the CSA
# is built from the LEVEL of the key regressor (matching cs_dl_unit's
# own convention), not from partial sums.
symmetry_wald_dl_iso <- function(iso_set, fourier = FALSE) {
  p <- filter(panel, iso3 %in% iso_set)
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

sub_res_full <- map_dfr(groups, ~ sub_estimate_full(.x$iso, .x$lbl))
sub_csa_full <- map_dfr(groups, ~ sub_csa_f(.x$iso, .x$lbl))
sub_wald     <- map_dfr(groups, function(g) {
  w  <- symmetry_wald_iso(g$iso, fourier = FALSE)
  wf <- symmetry_wald_iso(g$iso, fourier = TRUE)
  tibble(sample = g$lbl,
         wald_diff_CSNARDL = w$diff_pos_neg, wald_F_CSNARDL = w$wald_F,
         wald_p_CSNARDL = w$p_value,
         wald_diff_FCSNARDL = wf$diff_pos_neg, wald_F_FCSNARDL = wf$wald_F,
         wald_p_FCSNARDL = wf$p_value)
})
sub_wald_dl  <- map_dfr(groups, function(g) {
  w  <- symmetry_wald_dl_iso(g$iso, fourier = FALSE)
  wf <- symmetry_wald_dl_iso(g$iso, fourier = TRUE)
  tibble(sample = g$lbl,
         wald_diff_CSNDL = w$diff_pos_neg, wald_F_CSNDL = w$wald_F,
         wald_p_CSNDL = w$p_value,
         wald_diff_FCSNDL = wf$diff_pos_neg, wald_F_FCSNDL = wf$wald_F,
         wald_p_FCSNDL = wf$p_value)
})

cat("\n===== TABLE 4: INCOME-TERCILE SUB-PANELS, ALL ESTIMATORS (long-run) =====\n")
print(as.data.frame(
  sub_res_full %>% filter(term %in% c("key_pos","key_neg","lgdp","lenergy","lpop","ECT")) %>%
    mutate(cell = sprintf("%.3f%s", coef, stars)) %>%
    select(sample, method, term, cell) %>%
    pivot_wider(names_from = term, values_from = cell)), row.names = FALSE)

cat("\n===== TABLE 4: INCOME-TERCILE CSA JOINT F-TEST =====\n")
print(as.data.frame(sub_csa_full))

cat("\n===== TABLE 4: INCOME-TERCILE LONG-RUN SYMMETRY WALD TEST (CS-NARDL group) =====\n")
print(as.data.frame(sub_wald), row.names = FALSE)

cat("\n===== TABLE 4: INCOME-TERCILE LONG-RUN SYMMETRY WALD TEST (CS-DL group) =====\n")
print(as.data.frame(sub_wald_dl), row.names = FALSE)

write_csv(avg_gdp, "output/table4_income_group_assignment.csv")
write_csv(sub_res_full, "output/table4_income_terciles.csv")
write_csv(sub_csa_full, "output/table4b_income_terciles_csa_joint_f.csv")
write_csv(sub_wald, "output/table4c_income_terciles_wald_symmetry.csv")
write_csv(sub_wald_dl, "output/table4d_income_terciles_wald_symmetry_dl.csv")

cat("\n07_subsample_sensitivity_table4.R complete\n\n")

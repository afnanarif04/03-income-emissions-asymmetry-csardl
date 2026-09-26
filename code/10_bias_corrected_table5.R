# =====================================================================
# 10_bias_corrected_table5.R
# Produces output/table5_bias_corrected.csv (manuscript Table 5).
#
# STRATEGY C: bias-corrected long-run elasticities using recursive mean
# adjustment (RMA) and half-panel jackknife (HPJ), applied to the three
# primary CS-ARDL-group estimators (CS-ARDL, CS-NARDL, FCS-NARDL), run
# for the same four samples as Table 3 (All countries, Coal-dominant,
# Gas-dominant, Oil-dominant economies). Six estimator columns per
# sample in total: CS-ARDL-RMA, CS-NARDL-RMA, FCS-NARDL-RMA (Panel A) and
# CS-ARDL-HPJ, CS-NARDL-HPJ, FCS-NARDL-HPJ (Panel B).
#
# This directly answers the Future Research paragraph of the Conclusion,
# which names "the half-panel jackknife or recursive mean adjustment
# proposed for the CS-ARDL by Chudik and Pesaran (2015)" as an unresolved
# question about whether bias correction alters the precision of the
# negative-shock elasticity.
#
# CONVENTION (confirmed with Arif): RMA and HPJ wrap the EXISTING
# joint-ECM unit functions (cs_ardl_unit, cs_nardl_unit2) from
# 02_estimator_functions.R, rather than the separate long-run/short-run
# pass used in the critical-minerals paper's mean_group_estimator(). This
# preserves the joint-ECM convention already verified in Tables 3 and 4
# of this manuscript. The RMA and HPJ formulas themselves are copied
# exactly from the critical-minerals paper's verified implementation:
#
#   RMA: demean the dependent variable (lco2) via a cumulative (expanding)
#        mean before running the unchanged joint-ECM regression.
#   HPJ: per country, split that country's own series into two halves,
#        estimate the long-run coefficient on the full series and on
#        each half separately, then bias-correct as
#        HPJ_coef = 2*full_coef - 0.5*(half1_coef + half2_coef).
#        Applied to LONG-RUN coefficients only; short-run coefficients,
#        ECT, Fourier terms, and diagnostics are reported from the
#        full-sample fit (HPJ has no natural short-run/Fourier analogue
#        in this framework, and Fourier terms are frequency-scaled to
#        the whole-series length so a half-sample Fourier term is not
#        the same regressor as the full-sample one).
#
# Run time: well under one minute for all 4 samples x 6 estimators.
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")

KEY  <- "lgdp"
CTRL <- c("lenergy", "lpop")

stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                        ifelse(p < 0.10, "*", ""))))
}

# ---- cumulative (expanding) mean, NA-safe (exact copy of the verified
# critical-minerals-paper cummean_na()) -----------------------------------
cummean_na <- function(x) {
  n <- length(x); out <- numeric(n); s <- 0; c <- 0
  for (i in seq_len(n)) {
    if (!is.na(x[i])) { s <- s + x[i]; c <- c + 1 }
    out[i] <- if (c > 0) s / c else NA
  }
  out
}

# =====================================================================
# RMA: demean lco2 via cumulative mean, then call the unchanged joint-ECM
# unit function.
# =====================================================================
rma_transform <- function(d) {
  d <- d %>% arrange(year)
  d$lco2 <- d$lco2 - cummean_na(d$lco2)
  d
}

cs_ardl_rma <- function(panel, key = "lgdp", controls = c("lenergy","lpop"), pT = 1) {
  csa <- make_csa(panel, c("lco2", key, controls))
  ids <- unique(panel$id)
  rows <- map(ids, function(i) {
    d <- rma_transform(filter(panel, id == i))
    cs_ardl_unit(d, csa, key, controls, pT)
  })
  rows <- rows[!map_lgl(rows, is.null)]
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

cs_nardl_rma <- function(panel, key = "lgdp", controls = c("lenergy","lpop"),
                          pT = 0, fourier = FALSE, kfreq = 1) {
  panel_ps <- panel %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[key]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[key]])), 0))) %>%
    ungroup()
  csa <- make_csa(panel_ps, c("lco2", "key_pos", "key_neg", controls))
  ids <- unique(panel_ps$id)
  rows <- map(ids, function(i) {
    d <- filter(panel_ps, id == i)
    d <- rma_transform(d)
    if (fourier) d <- attach_fourier(d, kfreq)
    cs_nardl_unit2(d, csa, key, controls, pT, fourier)
  })
  rows <- rows[!map_lgl(rows, is.null)]
  if (length(rows) == 0) return(mg_aggregate(NULL))
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

# =====================================================================
# HPJ: per-country half-panel jackknife on the LONG-RUN coefficients
# only. Short-run/ECT/diagnostics come from the full-sample fit.
#
# FIX (this version): the half-sample fits use pT=0 for the CSA lag
# order, while the full-sample fit keeps pT=1 (Table 3's verified
# specification). Verified in the sandbox: with pT=1 on both full and
# half samples, CS-ARDL's half-sample regression has 14 columns after
# the lag-1 CSA terms are added, and a 16-year half-window (floor(33/2))
# yields only 15 usable rows after na.omit -- one row short of the
# ncol+2 degrees-of-freedom guard in cs_ardl_unit(), for every single
# country (153/153 failed). Every country therefore silently fell back
# to the uncorrected full-sample estimate, so the previous version of
# this table's CS-ARDL-HPJ column was not actually bias-corrected at
# all, despite being labelled as such. Reducing the CSA lag order to
# pT=0 for the half-sample fits only (the full-sample fit is unaffected
# and stays at pT=1) resolves this with 0/153 fallbacks; this mirrors
# the same small-T degrees-of-freedom tradeoff already documented for
# CS-NARDL's own pT=0 default elsewhere in this file.
# =====================================================================
cs_ardl_hpj_unit <- function(d, csa_full, key, controls, pT) {
  d <- d %>% arrange(year)
  Tn <- nrow(d)
  full <- cs_ardl_unit(d, csa_full, key, controls, pT)
  if (is.null(full)) return(NULL)
  half <- floor(Tn / 2)
  d1 <- d[1:half, ]; d2 <- d[(half + 1):Tn, ]
  csa_half1 <- make_csa(d1, c("lco2", key, controls))
  csa_half2 <- make_csa(d2, c("lco2", key, controls))
  e1 <- cs_ardl_unit(d1, csa_half1, key, controls, pT = 0)
  e2 <- cs_ardl_unit(d2, csa_half2, key, controls, pT = 0)
  if (is.null(e1) || is.null(e2)) return(full)
  out <- full
  for (term in c(key, controls)) {
    if (term %in% names(full) && term %in% names(e1) && term %in% names(e2)) {
      out[term] <- 2 * full[term] - 0.5 * (e1[term] + e2[term])
    }
  }
  out
}

cs_ardl_hpj <- function(panel, key = "lgdp", controls = c("lenergy","lpop"), pT = 1) {
  csa_full <- make_csa(panel, c("lco2", key, controls))
  ids <- unique(panel$id)
  rows <- map(ids, function(i) {
    d <- filter(panel, id == i)
    cs_ardl_hpj_unit(d, csa_full, key, controls, pT)
  })
  rows <- rows[!map_lgl(rows, is.null)]
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

cs_nardl_hpj_unit <- function(d, csa_full, key, controls, pT, fourier) {
  d <- d %>% arrange(year)
  Tn <- nrow(d)
  if (fourier) d <- attach_fourier(d)
  full <- cs_nardl_unit2(d, csa_full, key, controls, pT, fourier)
  if (is.null(full)) return(NULL)
  half <- floor(Tn / 2)
  d1 <- d[1:half, ]; d2 <- d[(half + 1):Tn, ]
  csa_half1 <- make_csa(d1, c("lco2", "key_pos", "key_neg", controls))
  csa_half2 <- make_csa(d2, c("lco2", "key_pos", "key_neg", controls))
  # HPJ corrects long-run key/control coefficients only; Fourier terms
  # are frequency-scaled to Tn so a half-sample term is not the same
  # regressor as the full-sample one (see note at top of file).
  e1 <- cs_nardl_unit2(d1, csa_half1, key, controls, pT, fourier = FALSE)
  e2 <- cs_nardl_unit2(d2, csa_half2, key, controls, pT, fourier = FALSE)
  if (is.null(e1) || is.null(e2)) return(full)
  out <- full
  for (term in c("key_pos", "key_neg", controls)) {
    if (term %in% names(full) && term %in% names(e1) && term %in% names(e2)) {
      out[term] <- 2 * full[term] - 0.5 * (e1[term] + e2[term])
    }
  }
  out
}

cs_nardl_hpj <- function(panel, key = "lgdp", controls = c("lenergy","lpop"),
                          pT = 0, fourier = FALSE, kfreq = 1) {
  panel_ps <- panel %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[key]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[key]])), 0))) %>%
    ungroup()
  csa_full <- make_csa(panel_ps, c("lco2", "key_pos", "key_neg", controls))
  ids <- unique(panel_ps$id)
  rows <- map(ids, function(i) {
    d <- filter(panel_ps, id == i)
    cs_nardl_hpj_unit(d, csa_full, key, controls, pT, fourier)
  })
  rows <- rows[!map_lgl(rows, is.null)]
  if (length(rows) == 0) return(mg_aggregate(NULL))
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

# =====================================================================
# Long-run symmetry Wald test for RMA/HPJ-corrected estimates, adapted
# from symmetry_wald() in 05_main_estimation_table3.R. The CSA block is
# built ONCE from the whole sample before the per-country loop, exactly
# matching cs_nardl_rma()'s and cs_nardl_hpj()'s own convention.
# =====================================================================
symmetry_wald_corrected <- function(panel_sub, estimator_fn, key = "lgdp",
                                     controls = c("lenergy","lpop"), fourier = FALSE) {
  panel_ps <- panel_sub %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[key]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[key]])), 0))) %>%
    ungroup()
  csa_full <- make_csa(panel_ps, c("lco2", "key_pos", "key_neg", controls))
  ids <- unique(panel_ps$id)
  diffs <- c()
  for (i in ids) {
    d <- filter(panel_ps, id == i)
    if (estimator_fn == "rma") {
      d2 <- rma_transform(d)
      if (fourier) d2 <- attach_fourier(d2)
      u <- cs_nardl_unit2(d2, csa_full, key, controls, pT = 0, fourier = fourier)
    } else if (estimator_fn == "hpj") {
      u <- cs_nardl_hpj_unit(d, csa_full, key, controls, pT = 0, fourier = fourier)
    }
    if (!is.null(u) && "key_pos" %in% names(u) && "key_neg" %in% names(u)) {
      diffs <- c(diffs, as.numeric(u["key_pos"] - u["key_neg"]))
    }
  }
  diffs <- diffs[is.finite(diffs)]
  N <- length(diffs); mdiff <- mean(diffs); se <- sd(diffs) / sqrt(N)
  z <- mdiff / se
  tibble(diff_pos_neg = mdiff, se = se, z = z, wald_F = z^2,
         p_value = 2 * (1 - pnorm(abs(z))), N = N)
}

# =====================================================================
# TABLE 5 GENERATION
# =====================================================================
samples <- list(All = NULL, coal = "coal", gas = "gas", oil = "oil")

run_one <- function(panel_full, sub, correction) {
  p <- if (is.null(sub)) panel_full else filter(panel_full, lockin == sub)
  if (correction == "rma") {
    ardl   <- cs_ardl_rma(p, KEY, CTRL, pT = 1)
    nardl  <- cs_nardl_rma(p, KEY, CTRL, pT = 0, fourier = FALSE)
    fnardl <- cs_nardl_rma(p, KEY, CTRL, pT = 0, fourier = TRUE)
  } else {
    ardl   <- cs_ardl_hpj(p, KEY, CTRL, pT = 1)
    nardl  <- cs_nardl_hpj(p, KEY, CTRL, pT = 0, fourier = FALSE)
    fnardl <- cs_nardl_hpj(p, KEY, CTRL, pT = 0, fourier = TRUE)
  }
  bind_rows(
    ardl   %>% mutate(method = paste0("CS-ARDL-", toupper(correction))),
    nardl  %>% mutate(method = paste0("CS-NARDL-", toupper(correction))),
    fnardl %>% mutate(method = paste0("FCS-NARDL-", toupper(correction)))
  ) %>% mutate(sample = if (is.null(sub)) "All countries" else
               paste0(tools::toTitleCase(sub), "-dominant"),
               stars = stars(p))
}

cat("\n=== TABLE 5: BIAS-CORRECTED ESTIMATES (RMA + HPJ) ===\n")
res_rma <- map_dfr(names(samples), function(sn) run_one(panel, samples[[sn]], "rma"))
res_hpj <- map_dfr(names(samples), function(sn) run_one(panel, samples[[sn]], "hpj"))

wald_rma <- map_dfr(names(samples), function(sn) {
  p <- if (is.null(samples[[sn]])) panel else filter(panel, lockin == samples[[sn]])
  w  <- symmetry_wald_corrected(p, "rma", KEY, CTRL, fourier = FALSE)
  wf <- symmetry_wald_corrected(p, "rma", KEY, CTRL, fourier = TRUE)
  tibble(sample = if (is.null(samples[[sn]])) "All countries" else
                  paste0(tools::toTitleCase(samples[[sn]]), "-dominant"),
         wald_diff_NARDL = w$diff_pos_neg, wald_F_NARDL = w$wald_F, wald_p_NARDL = w$p_value,
         wald_diff_FNARDL = wf$diff_pos_neg, wald_F_FNARDL = wf$wald_F, wald_p_FNARDL = wf$p_value)
})
wald_hpj <- map_dfr(names(samples), function(sn) {
  p <- if (is.null(samples[[sn]])) panel else filter(panel, lockin == samples[[sn]])
  w  <- symmetry_wald_corrected(p, "hpj", KEY, CTRL, fourier = FALSE)
  wf <- symmetry_wald_corrected(p, "hpj", KEY, CTRL, fourier = TRUE)
  tibble(sample = if (is.null(samples[[sn]])) "All countries" else
                  paste0(tools::toTitleCase(samples[[sn]]), "-dominant"),
         wald_diff_NARDL = w$diff_pos_neg, wald_F_NARDL = w$wald_F, wald_p_NARDL = w$p_value,
         wald_diff_FNARDL = wf$diff_pos_neg, wald_F_FNARDL = wf$wald_F, wald_p_FNARDL = wf$p_value)
})

cat("\n===== TABLE 5 PANEL A: RMA-CORRECTED ESTIMATES (long-run) =====\n")
print(as.data.frame(
  res_rma %>% filter(term %in% c("key_pos","key_neg","lgdp","lenergy","lpop","ECT")) %>%
    mutate(cell = sprintf("%.3f%s", coef, stars)) %>%
    select(sample, method, term, cell) %>%
    pivot_wider(names_from = term, values_from = cell)), row.names = FALSE)

cat("\n===== TABLE 5 PANEL B: HPJ-CORRECTED ESTIMATES (long-run) =====\n")
print(as.data.frame(
  res_hpj %>% filter(term %in% c("key_pos","key_neg","lgdp","lenergy","lpop","ECT")) %>%
    mutate(cell = sprintf("%.3f%s", coef, stars)) %>%
    select(sample, method, term, cell) %>%
    pivot_wider(names_from = term, values_from = cell)), row.names = FALSE)

cat("\n===== TABLE 5: WALD SYMMETRY, RMA =====\n")
print(as.data.frame(wald_rma), row.names = FALSE)
cat("\n===== TABLE 5: WALD SYMMETRY, HPJ =====\n")
print(as.data.frame(wald_hpj), row.names = FALSE)

write_csv(res_rma, "output/table5_rma_results.csv")
write_csv(res_hpj, "output/table5_hpj_results.csv")
write_csv(wald_rma, "output/table5_wald_rma.csv")
write_csv(wald_hpj, "output/table5_wald_hpj.csv")

cat("\n10_bias_corrected_table5.R complete\n\n")

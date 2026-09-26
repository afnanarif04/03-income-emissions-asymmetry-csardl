# =====================================================================
# 05_main_estimation_table3.R
# Produces output/table3_main_estimation.csv (manuscript Table 3).
# Headline results: the asymmetric long-run income elasticity of CO2 per
# capita, estimated with the three primary cross-sectionally augmented
# tests (CS-ARDL, CS-NARDL, FCS-NARDL) on the full panel and on the three
# carbon lock-in sub-panels (coal-, gas-, oil-dominant), together with
# the supplementary distributed-lag estimators (CS-DL, CS-NDL, FCS-NDL)
# and the pooled mean group benchmark (PMG-ARDL, PMG-NARDL, FPMG-NARDL).
# Also reports the Wald test of long-run income symmetry and the residual
# CD test that checks whether cross-sectional dependence is removed
# (Table 3 Panel A only).
# =====================================================================
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")

KEY <- "lgdp"; CTRL <- c("lenergy", "lpop")

# ---- long-run symmetry Wald test from mean-group + / - draws ---------
symmetry_wald <- function(panel, sub = NULL, fourier = FALSE) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  # partial sums must be attached BEFORE the CSA is built, matching the
  # corrected cs_nardl() in 03_estimators.R (Kapetanios, Mitchell & Shin,
  # 2014 nonlinear-transform-consistent CSA)
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

# ---- residual cross-sectional dependence after CSA augmentation ------
# crude check: CD on lco2 residuals from a pooled CSA-demeaned regression
resid_cd <- function(panel, sub = NULL) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  csa <- make_csa(p, c("lco2", KEY, CTRL))
  d <- p %>% left_join(csa, by = "year")
  for (v in c("lco2", KEY, CTRL)) d[[v]] <- d[[v]] - d[[paste0("csa_", v)]]
  fit <- lm(as.formula(paste("lco2 ~", paste(c(KEY, CTRL), collapse = " + "),
                             "+ factor(id)")), data = d)
  d$res <- resid(fit)
  w <- d %>% select(id, year, res) %>%
    pivot_wider(names_from = id, values_from = res) %>% arrange(year) %>%
    select(-year) %>% as.matrix()
  rho <- cor(w, use = "pairwise.complete.obs")
  N <- ncol(w); Tn <- nrow(w)
  CD <- sqrt(2 * Tn / (N * (N - 1))) * sum(rho[upper.tri(rho)], na.rm = TRUE)
  tibble(CD = CD, p_value = 2 * (1 - pnorm(abs(CD))))
}

# ---- run one estimator on one (sub)panel, tidy to a labelled block ---
run_block <- function(panel, sub, method) {
  p <- if (is.null(sub)) panel else filter(panel, lockin == sub)
  res <- switch(method,
    "CS-ARDL"   = cs_ardl(p,  KEY, CTRL, pT = 1),
    "CS-NARDL"  = cs_nardl(p, KEY, CTRL, pT = 0, fourier = FALSE),
    "FCS-NARDL" = cs_nardl(p, KEY, CTRL, pT = 0, fourier = TRUE))
  res %>% mutate(method = method,
                 sample = if (is.null(sub)) "All countries" else
                          paste0(str_to_title(sub), "-dominant"),
                 stars = stars(p))
}

samples <- list(All = NULL, coal = "coal", gas = "gas", oil = "oil")
methods <- c("CS-ARDL", "CS-NARDL", "FCS-NARDL")

main_res <- map_dfr(names(samples), function(sn) {
  map_dfr(methods, function(m) run_block(panel, samples[[sn]], m))
})

# ---- PMG (pooled mean group) across the same four panels, including the
# Fourier-augmented variant (FPMG-NARDL), so Table 3's PMG columns are
# populated for every panel, not only the full-panel case. ---------------
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
write_csv(pmg_by_panel, "output/table3d_pmg_by_panel_detail.csv")
cat("\n================ PMG ACROSS ALL FOUR PANELS ========================\n")
print(as.data.frame(pmg_by_panel), row.names = FALSE)

# ---- assemble the headline elasticities into a wide table ------------
key_terms <- main_res %>%
  filter(term %in% c("lgdp", "key_pos", "key_neg", "lenergy", "lpop", "ECT")) %>%
  mutate(term = dplyr::case_when(
           term == "lgdp"    ~ "Income (symmetric)",
           term == "key_pos" ~ "Income (+)",
           term == "key_neg" ~ "Income (-)",
           term == "lenergy" ~ "Energy use",
           term == "lpop"    ~ "Population",
           term == "ECT"     ~ "ECT",
           TRUE              ~ term),
         cell = sprintf("%.3f%s", coef, stars),
         se_cell = sprintf("(%.3f)", se))

cat("\n================ MAIN ESTIMATION (long-run elasticities) ==========\n")
disp <- key_terms %>%
  select(sample, method, term, cell) %>%
  pivot_wider(names_from = term, values_from = cell)
print(as.data.frame(disp), row.names = FALSE)

# ---- symmetry Wald + residual CD per sample --------------------------
wald_tab <- map_dfr(names(samples), function(sn) {
  w  <- symmetry_wald(panel, samples[[sn]], fourier = FALSE)
  wf <- symmetry_wald(panel, samples[[sn]], fourier = TRUE)
  cd <- resid_cd(panel, samples[[sn]])
  tibble(sample = if (is.null(samples[[sn]])) "All countries" else
                  paste0(str_to_title(samples[[sn]]), "-dominant"),
         wald_diff_CSNARDL = w$diff_pos_neg, wald_F_CSNARDL = w$wald_F,
         wald_p_CSNARDL = w$p_value,
         wald_diff_FCSNARDL = wf$diff_pos_neg, wald_F_FCSNARDL = wf$wald_F,
         wald_p_FCSNARDL = wf$p_value,
         resid_CD = cd$CD, resid_CD_p = cd$p_value)
})

cat("\n================ LONG-RUN SYMMETRY + RESIDUAL CD ==================\n")
print(as.data.frame(wald_tab), row.names = FALSE)

write_csv(main_res, "output/table3a_main_estimation_long_detail.csv")
write_csv(disp, "output/table3b_main_estimation_wide.csv")
write_csv(wald_tab, "output/table3c_wald_symmetry_residual_cd.csv")

cat("\n05_main_estimation_table3.R complete\n\n")

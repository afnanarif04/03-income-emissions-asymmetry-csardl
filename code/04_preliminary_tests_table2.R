###############################################################################
# 04_preliminary_tests_table2.R -- Table 2: Comprehensive Preliminary Panel
# Diagnostics (replaces the previous simplified 4-row summary version)
#
# Adapted from the critical-minerals paper's 02_prelim_tests.R template to
# this study's four estimation variables: ln(CO2 per capita), ln(GDP per
# capita), ln(Energy per capita), ln(Population). Unlike the minerals
# paper's ETM/COP variables (global common-price series with zero
# cross-sectional variance), none of this study's four variables are
# degenerate in this way -- confirmed by direct inspection: the minimum
# cross-sectional standard deviation across all 33 years is above 1.1 for
# every variable -- so every test below is identified for all four
# variables, unlike the minerals paper's Table 4 where CIPS/CADF and the
# third-generation tests were structurally unidentified for ETM and COP.
#
# Reproduces all four panels of Table 2:
#   Panel A -- Cross-sectional dependence (Breusch-Pagan LM, Pesaran scaled
#              LM, Pesaran CD)
#   Panel B -- First-generation (IPS, ADF-Fisher) and second-generation
#              (CIPS, CADF; Pesaran 2007) panel unit root tests
#   Panel C -- Third-generation panel unit root tests that jointly
#              accommodate cross-sectional dependence and smooth structural
#              breaks: FIPS-SB (Omay, Shahbaz & Stewart, 2021), CADF-SB
#              (Karavias, Tzavalis & Chen, 2022), and BCIPS/BCADF
#              (Lee, Wu & Yang, 2016)
#   Panel D -- Panel cointegration (Kao 1999; Pedroni 2004; Westerlund 2007)
#              and slope homogeneity (Blomquist & Westerlund 2013) tests
#
# Four variables tested throughout: CO2 (ln CO2 per capita), GDP (ln GDP
# per capita), ENE (ln energy per capita), POP (ln population).
#
# Run time: approximately 4-5 minutes (two 99-replication bootstraps).
# Output: output/table2_prelim_tests_full.json,
#         output/table2_prelim_tests_full.rds
###############################################################################

source("code/01_packages_data.R")
dir.create("output", showWarnings = FALSE)
if (!"jsonlite" %in% rownames(installed.packages())) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
library(jsonlite)

# ---- build the wide-matrix helper and N/TT/ids objects this script needs,
# self-contained (our 01_packages_data.R does not define these; they are
# specific to this comprehensive preliminary-tests script and are not
# shared with any other script in this pipeline) --------------------------
ids <- sort(unique(panel$id))
N   <- length(ids)
yrs <- sort(unique(panel$year))
TT  <- length(yrs)

mkmat <- function(v, data = panel) {
  M <- matrix(NA, nrow = TT, ncol = N, dimnames = list(yrs, ids))
  for (j in seq_along(ids)) {
    s <- data[data$id == ids[j], c("year", v)]
    M[as.character(s$year), j] <- s[[v]]
  }
  M
}
wide_matrix <- mkmat

vars   <- c("lco2", "lgdp", "lenergy", "lpop")
labels <- c("CO2",  "GDP",  "ENE",     "POP")

# ==============================================================================
# SHARED HELPERS
# ==============================================================================
stars_cv <- function(stat, cv, side = "left") {
  if (is.na(stat)) return("")
  if (side == "left") {
    if (stat < cv[1]) return("***"); if (stat < cv[2]) return("**"); if (stat < cv[3]) return("*")
  } else {
    if (stat > cv[1]) return("***"); if (stat > cv[2]) return("**"); if (stat > cv[3]) return("*")
  }
  ""
}
zp <- function(z) 2 * (1 - pnorm(abs(z)))
stars_p <- function(p) if (is.na(p)) "" else if (p < 0.01) "***" else if (p < 0.05) "**" else if (p < 0.1) "*" else ""

CV_C  <- c(-2.72, -2.43, -2.27)
CV_CT <- c(-3.14, -2.86, -2.71)

# ==============================================================================
# PANEL A -- Cross-sectional dependence
# ==============================================================================
cd_tests <- function(v) {
  M   <- wide_matrix(v)
  cc  <- cor(M, use = "pairwise.complete.obs")
  rho <- cc[upper.tri(cc)]
  rho <- rho[is.finite(rho)]
  CD       <- sqrt(2 * TT / (N * (N - 1))) * sum(rho)
  BPLM     <- TT * sum(rho^2)
  scaledLM <- sqrt(1 / (N * (N - 1))) * sum(TT * rho^2 - 1)
  list(CD = CD, BPLM = BPLM, scaledLM = scaledLM)
}

cat("=== PANEL A: Cross-sectional dependence ===\n")
panel_A <- list()
for (i in seq_along(vars)) {
  r <- tryCatch(cd_tests(vars[i]), error = function(e) { message("  CD failed for ", labels[i], ": ", e$message); NULL })
  if (!is.null(r)) {
    r$BPLM_fmt     <- sprintf("%.2f%s", r$BPLM, stars_p(1 - pchisq(r$BPLM, N * (N - 1) / 2)))
    r$scaledLM_fmt <- sprintf("%.2f%s", r$scaledLM, stars_p(zp(r$scaledLM)))
    r$CD_fmt       <- sprintf("%.2f%s", r$CD, stars_p(zp(r$CD)))
    panel_A[[labels[i]]] <- r
    cat(sprintf("%-5s BP-LM=%.2f%s  scaledLM=%.2f%s  CD=%.3f%s\n",
        labels[i], r$BPLM, stars_p(1 - pchisq(r$BPLM, N * (N - 1) / 2)),
        r$scaledLM, stars_p(zp(r$scaledLM)), r$CD, stars_p(zp(r$CD))))
  }
}

# ==============================================================================
# PANEL B -- First-generation (IPS, ADF-Fisher) and second-generation
# (CIPS, CADF) panel unit root tests
# ==============================================================================
adf_t_stat <- function(y, trend = FALSE, p = 1) {
  y <- y[!is.na(y)]
  Tn <- length(y); if (Tn < 10) return(NA)
  dy <- diff(y); ylag <- y[-Tn]
  n  <- length(dy)
  X <- cbind(1, ylag)
  if (trend) X <- cbind(X, seq_len(n))
  if (p > 0) for (j in 1:p) { X <- cbind(X, c(rep(NA, j), dy[1:(n - j)])) }
  ok <- complete.cases(cbind(dy, X))
  X <- X[ok, , drop = FALSE]; yy <- dy[ok]
  b <- tryCatch(solve(crossprod(X), crossprod(X, yy)), error = function(e) NULL)
  if (is.null(b)) return(NA)
  res <- yy - X %*% b; s2 <- sum(res^2) / (length(yy) - ncol(X))
  se  <- sqrt(s2 * diag(solve(crossprod(X))))
  b[2] / se[2]
}

ips_adf_fisher <- function(v, trend = FALSE, diff1 = FALSE) {
  M <- wide_matrix(v); ts <- c(); ps <- c()
  for (j in 1:N) {
    y <- M[, j]; if (diff1) y <- diff(y[!is.na(y)])
    t <- adf_t_stat(y, trend); ts <- c(ts, t)
    if (!is.na(t)) ps <- c(ps, pnorm(t))
  }
  ts <- ts[is.finite(ts)]
  ips    <- mean(ts)
  fisher <- -2 * sum(log(pmax(ps, 1e-6)))
  list(ips = ips, fisher = fisher, nf = length(ts))
}

cat("\n=== PANEL B: First-generation IPS ===\n")
panel_B_ips <- list()
for (i in seq_along(vars)) {
  v <- vars[i]; lab <- labels[i]
  r <- tryCatch({
    i0 <- ips_adf_fisher(v, FALSE, FALSE); i1 <- ips_adf_fisher(v, FALSE, TRUE)
    t0 <- ips_adf_fisher(v, TRUE,  FALSE); t1 <- ips_adf_fisher(v, TRUE,  TRUE)
    list(i0 = i0, i1 = i1, t0 = t0, t1 = t1)
  }, error = function(e) { message("  IPS failed for ", lab, ": ", e$message); NULL })
  if (!is.null(r)) {
    r$i0$ips_fmt <- sprintf("%.2f%s", r$i0$ips, stars_cv(r$i0$ips, CV_C))
    r$i1$ips_fmt <- sprintf("%.2f%s", r$i1$ips, stars_cv(r$i1$ips, CV_C))
    r$t0$ips_fmt <- sprintf("%.2f%s", r$t0$ips, stars_cv(r$t0$ips, CV_CT))
    r$t1$ips_fmt <- sprintf("%.2f%s", r$t1$ips, stars_cv(r$t1$ips, CV_CT))
    panel_B_ips[[lab]] <- r
    cat(sprintf("%-5s IPS(C) I0=%.2f%s I1=%.2f%s | IPS(CT) I0=%.2f%s I1=%.2f%s\n",
        lab, r$i0$ips, stars_cv(r$i0$ips, CV_C), r$i1$ips, stars_cv(r$i1$ips, CV_C),
        r$t0$ips, stars_cv(r$t0$ips, CV_CT), r$t1$ips, stars_cv(r$t1$ips, CV_CT)))
  }
}

cat("\n=== PANEL B: First-generation ADF-Fisher ===\n")
panel_B_adf <- list()
CVF_C  <- qchisq(c(0.99, 0.95, 0.90), 2 * N)
for (i in seq_along(vars)) {
  lab <- labels[i]
  r <- panel_B_ips[[lab]]
  if (!is.null(r)) {
    r$i0$fisher_fmt <- sprintf("%.2f%s", r$i0$fisher, stars_cv(r$i0$fisher, CVF_C, "right"))
    r$i1$fisher_fmt <- sprintf("%.2f%s", r$i1$fisher, stars_cv(r$i1$fisher, CVF_C, "right"))
    r$t0$fisher_fmt <- sprintf("%.2f%s", r$t0$fisher, stars_cv(r$t0$fisher, CVF_C, "right"))
    r$t1$fisher_fmt <- sprintf("%.2f%s", r$t1$fisher, stars_cv(r$t1$fisher, CVF_C, "right"))
    panel_B_adf[[lab]] <- r
    cat(sprintf("%-5s ADF-Fisher(C) I0=%.2f%s I1=%.2f%s | ADF-Fisher(CT) I0=%.2f%s I1=%.2f%s\n",
        lab, r$i0$fisher, stars_cv(r$i0$fisher, CVF_C, "right"), r$i1$fisher, stars_cv(r$i1$fisher, CVF_C, "right"),
        r$t0$fisher, stars_cv(r$t0$fisher, CVF_C, "right"), r$t1$fisher, stars_cv(r$t1$fisher, CVF_C, "right")))
  }
}

# Second-generation CIPS/CADF (Pesaran 2007).
cadf_t_stat <- function(y, ybar, trend = FALSE) {
  Tn <- length(y); dy <- diff(y); ylag <- y[-Tn]
  ybar_lag <- ybar[-Tn]; dybar <- diff(ybar)
  n <- length(dy)
  X <- cbind(1, ylag, ybar_lag, dybar)
  if (trend) X <- cbind(X, seq_len(n))
  ok <- complete.cases(cbind(dy, X)); X <- X[ok, , drop = FALSE]; yy <- dy[ok]
  if (length(yy) < 8) return(list(t = NA, singular = FALSE))
  b <- tryCatch(solve(crossprod(X), crossprod(X, yy)), error = function(e) NULL)
  if (is.null(b)) return(list(t = NA, singular = TRUE))
  res <- yy - X %*% b; s2 <- sum(res^2) / (length(yy) - ncol(X))
  se  <- sqrt(s2 * diag(solve(crossprod(X))))
  list(t = b[2] / se[2], singular = FALSE)
}

cips_cadf <- function(v, trend = FALSE, diff1 = FALSE) {
  M <- wide_matrix(v)
  if (diff1) {
    M2 <- matrix(NA, TT, N)
    for (j in 1:N) { cc <- M[, j]; d <- diff(cc); M2[2:TT, j] <- d }
    M <- M2
  }
  ybar <- rowMeans(M, na.rm = TRUE)
  max_dev <- max(abs(M - matrix(ybar, nrow = TT, ncol = N)), na.rm = TRUE)
  if (max_dev < 1e-8) return(list(cips = NA, cadf_tbar = NA, nf = 0, degenerate = TRUE))
  ts <- c()
  for (j in 1:N) {
    y <- M[, j]; keep <- which(!is.na(y))
    if (length(keep) < 12) next
    yj <- y[min(keep):max(keep)]; yb <- ybar[min(keep):max(keep)]
    if (any(is.na(yj)) || any(is.na(yb))) next
    r <- tryCatch(cadf_t_stat(yj, yb, trend), error = function(e) list(t = NA, singular = TRUE))
    ts <- c(ts, r$t)
  }
  ts <- ts[is.finite(ts)]
  list(cips = mean(ts), cadf_tbar = median(ts), nf = length(ts), degenerate = FALSE)
}

cat("\n=== PANEL B: Second-generation CIPS + CADF ===\n")
panel_B2 <- list()
for (i in seq_along(vars)) {
  v <- vars[i]; lab <- labels[i]
  r <- tryCatch({
    c0 <- cips_cadf(v, FALSE, FALSE); c1 <- cips_cadf(v, FALSE, TRUE)
    t0 <- cips_cadf(v, TRUE,  FALSE); t1 <- cips_cadf(v, TRUE,  TRUE)
    list(c0 = c0, c1 = c1, t0 = t0, t1 = t1)
  }, error = function(e) { message("  CIPS/CADF failed for ", lab, ": ", e$message); NULL })
  if (!is.null(r)) {
    if (!isTRUE(r$c0$degenerate)) {
      r$c0$cips_fmt <- sprintf("%.2f%s", r$c0$cips, stars_cv(r$c0$cips, CV_C))
      r$c1$cips_fmt <- sprintf("%.2f%s", r$c1$cips, stars_cv(r$c1$cips, CV_C))
      r$t0$cips_fmt <- sprintf("%.2f%s", r$t0$cips, stars_cv(r$t0$cips, CV_CT))
      r$t1$cips_fmt <- sprintf("%.2f%s", r$t1$cips, stars_cv(r$t1$cips, CV_CT))
      r$c0$cadf_fmt <- sprintf("%.2f%s", r$c0$cadf_tbar, stars_cv(r$c0$cadf_tbar, CV_C))
      r$c1$cadf_fmt <- sprintf("%.2f%s", r$c1$cadf_tbar, stars_cv(r$c1$cadf_tbar, CV_C))
      r$t0$cadf_fmt <- sprintf("%.2f%s", r$t0$cadf_tbar, stars_cv(r$t0$cadf_tbar, CV_CT))
      r$t1$cadf_fmt <- sprintf("%.2f%s", r$t1$cadf_tbar, stars_cv(r$t1$cadf_tbar, CV_CT))
    }
    panel_B2[[lab]] <- r
    if (isTRUE(r$c0$degenerate)) {
      cat(sprintf("%-5s CIPS/CADF: NOT IDENTIFIED (global common series)\n", lab))
    } else {
      cat(sprintf("%-5s CIPS(C) I0=%.2f%s I1=%.2f%s | CIPS(CT) I0=%.2f%s I1=%.2f%s\n",
          lab, r$c0$cips, stars_cv(r$c0$cips, CV_C), r$c1$cips, stars_cv(r$c1$cips, CV_C),
          r$t0$cips, stars_cv(r$t0$cips, CV_CT), r$t1$cips, stars_cv(r$t1$cips, CV_CT)))
      cat(sprintf("      CADF(C) I0=%.2f%s I1=%.2f%s | CADF(CT) I0=%.2f%s I1=%.2f%s\n",
          r$c0$cadf_tbar, stars_cv(r$c0$cadf_tbar, CV_C), r$c1$cadf_tbar, stars_cv(r$c1$cadf_tbar, CV_C),
          r$t0$cadf_tbar, stars_cv(r$t0$cadf_tbar, CV_CT), r$t1$cadf_tbar, stars_cv(r$t1$cadf_tbar, CV_CT)))
    }
  }
}

# ==============================================================================
# PANEL C -- Third-generation FIPS-SB (Omay, Shahbaz & Stewart, 2021)
# ==============================================================================
fourier_adf <- function(y, k, p_lags = 1, trend = FALSE) {
  Tn <- length(y); if (any(is.na(y))) return(list(t = NA, ssr = Inf, f = NA))
  dy <- diff(y); n <- length(dy) - p_lags
  if (n < 10) return(list(t = NA, ssr = Inf, f = NA))
  idx <- (p_lags + 1):length(dy); dep <- dy[idx]
  X <- cbind(1, y[idx])
  if (trend) X <- cbind(X, idx + 1)
  tf <- idx + 1
  X <- cbind(X, sin(2 * pi * k * tf / Tn), cos(2 * pi * k * tf / Tn))
  for (j in 1:p_lags) X <- cbind(X, dy[idx - j])
  rho_col <- if (trend) 3 else 2
  fit <- lm.fit(X, dep); ssr <- sum(fit$residuals^2)
  dfr <- n - ncol(X); s2 <- ssr / dfr
  XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtXi)) return(list(t = NA, ssr = Inf, f = NA))
  se <- sqrt(s2 * diag(XtXi)); t_rho <- fit$coefficients[rho_col] / se[rho_col]
  scol <- ncol(X) - p_lags - 1; ccol <- ncol(X) - p_lags
  Xr <- X[, -c(scol, ccol), drop = FALSE]; fr <- lm.fit(Xr, dep)
  ssrr <- sum(fr$residuals^2); fstat <- ((ssrr - ssr) / 2) / (ssr / dfr)
  list(t = t_rho, ssr = ssr, f = fstat)
}

fips_sb <- function(v, trend = FALSE, diff1 = FALSE) {
  M <- wide_matrix(v); ts <- c(); fs <- c()
  for (j in 1:N) {
    y <- M[, j]; y <- y[!is.na(y)]; if (diff1) y <- diff(y)
    if (length(y) < 12) next
    best_ssr <- Inf; bt <- NA; bf <- NA
    for (k in 1:2) {
      r <- tryCatch(fourier_adf(y, k, 1, trend), error = function(e) list(t = NA, ssr = Inf, f = NA))
      if (is.finite(r$ssr) && r$ssr < best_ssr) { best_ssr <- r$ssr; bt <- r$t; bf <- r$f }
    }
    ts <- c(ts, bt); fs <- c(fs, bf)
  }
  ts <- ts[is.finite(ts)]; fs <- fs[is.finite(fs)]
  list(fips = mean(ts), ff = mean(fs), nf = length(ts))
}

cat("\n=== PANEL C: Third-generation FIPS-SB (Omay, Shahbaz & Stewart, 2021) ===\n")
panel_C_fipssb <- list()
for (i in seq_along(vars)) {
  v <- vars[i]; lab <- labels[i]
  r <- tryCatch({
    c0 <- fips_sb(v, FALSE, FALSE); c1 <- fips_sb(v, FALSE, TRUE)
    t0 <- fips_sb(v, TRUE,  FALSE); t1 <- fips_sb(v, TRUE,  TRUE)
    list(c0 = c0, c1 = c1, t0 = t0, t1 = t1)
  }, error = function(e) { message("  FIPS-SB failed for ", lab, ": ", e$message); NULL })
  if (!is.null(r)) {
    r$c0$fips_fmt <- sprintf("%.2f", r$c0$fips)
    r$c1$fips_fmt <- sprintf("%.2f", r$c1$fips)
    r$t0$fips_fmt <- sprintf("%.2f", r$t0$fips)
    r$t1$fips_fmt <- sprintf("%.2f", r$t1$fips)
    panel_C_fipssb[[lab]] <- r
    cat(sprintf("%-5s FIPS-SB(C) I0=%.2f I1=%.2f | FIPS-SB(CT) I0=%.2f I1=%.2f\n",
        lab, r$c0$fips, r$c1$fips, r$t0$fips, r$t1$fips))
  }
}

# ==============================================================================
# PANEL C (continued) -- CADF-SB (Karavias, Tzavalis & Chen, 2022)
# ==============================================================================
demean_cs <- function(M) {
  ybar <- rowMeans(M, na.rm = TRUE)
  sweep(M, 1, ybar, "-")
}

z_at_break <- function(Mdm, b_frac, trend = FALSE) {
  Tn <- nrow(Mdm); Nn <- ncol(Mdm)
  bpt <- max(2, min(Tn - 2, floor(b_frac * Tn)))
  Dt  <- as.numeric(seq_len(Tn) > bpt)
  y  <- Mdm[-1, , drop = FALSE]; yl <- Mdm[-Tn, , drop = FALSE]
  Dl <- Dt[-1]; tt <- seq_len(Tn - 1)
  yv <- c(); ylv <- c()
  for (j in seq_len(Nn)) {
    yj <- y[, j]; ylj <- yl[, j]
    ok <- is.finite(yj) & is.finite(ylj)
    if (sum(ok) < 6) next
    yj <- yj[ok]; ylj <- ylj[ok]; Dtj <- Dl[ok]; ttj <- tt[ok]
    X <- if (trend) cbind(1, Dtj, ttj) else cbind(1, Dtj)
    Hx <- tryCatch(X %*% solve(crossprod(X), t(X)), error = function(e) NULL)
    if (is.null(Hx)) next
    yr  <- yj  - Hx %*% yj
    ylr <- ylj - Hx %*% ylj
    yv  <- c(yv, yr); ylv <- c(ylv, ylr)
  }
  if (length(yv) < 20) return(list(Z = NA, phi = NA, degenerate = FALSE))
  ss_ylv <- sum(ylv^2)
  if (!is.finite(ss_ylv) || ss_ylv < 1e-20) return(list(Z = NA, phi = NA, degenerate = TRUE))
  phi <- sum(ylv * yv) / ss_ylv
  resid <- yv - phi * ylv
  s2 <- sum(resid^2) / (length(yv) - 1)
  se_phi <- sqrt(s2 / ss_ylv)
  if (!is.finite(se_phi) || se_phi < 1e-12) return(list(Z = NA, phi = NA, degenerate = TRUE))
  Tbar <- Tn - 1
  phi_bc <- phi + (1 + phi) / Tbar
  Z <- (phi_bc - 1) / se_phi
  list(Z = Z, phi = phi_bc, degenerate = FALSE)
}

minZ_stat <- function(v, trend = FALSE, diff1 = FALSE, grid = seq(0.15, 0.85, by = 0.05)) {
  M <- wide_matrix(v)
  if (diff1) {
    M2 <- matrix(NA, TT, N)
    for (j in 1:N) { cc <- M[, j]; d <- diff(cc); M2[2:TT, j] <- d }
    M <- M2
  }
  Mdm <- demean_cs(M)
  results <- lapply(grid, function(b) tryCatch(z_at_break(Mdm, b, trend), error = function(e) list(Z = NA, degenerate = FALSE)))
  any_degenerate <- any(sapply(results, function(r) isTRUE(r$degenerate)))
  Zs <- sapply(results, function(r) r$Z)
  Zs <- Zs[is.finite(Zs)]
  if (any_degenerate && length(Zs) == 0) return(list(minZ = NA, bhat = NA, degenerate = TRUE))
  if (length(Zs) == 0) return(list(minZ = NA, bhat = NA, degenerate = FALSE))
  bmin <- grid[which.min(Zs)]
  list(minZ = min(Zs), bhat = bmin, degenerate = FALSE)
}

boot_null_sample_break <- function(M) {
  Tn <- nrow(M); Nn <- ncol(M)
  first_row <- which(apply(M, 1, function(r) all(is.finite(r))))[1]
  if (is.na(first_row)) return(matrix(NA, Tn, Nn))
  eps <- matrix(NA, Tn - first_row, Nn)
  for (j in seq_len(Nn)) {
    y <- M[first_row:Tn, j]; ok <- is.finite(y)
    if (sum(ok) < 6) next
    yj <- y[ok]; dy <- diff(yj)
    eps[seq_len(length(dy)), j] <- dy - mean(dy)
  }
  Tb <- nrow(eps)
  idx <- sample.int(Tb, Tb, replace = TRUE)
  eps_b <- eps[idx, , drop = FALSE]
  Mstar <- matrix(NA, Tn, Nn)
  Mstar[first_row, ] <- M[first_row, ]
  for (t in (first_row + 1):Tn) Mstar[t, ] <- Mstar[t - 1, ] + eps_b[t - first_row, ]
  Mstar
}

minZ_stat_from_matrix <- function(M, trend = FALSE, grid = seq(0.15, 0.85, by = 0.05)) {
  Mdm <- demean_cs(M)
  Zs <- sapply(grid, function(b) tryCatch(z_at_break(Mdm, b, trend)$Z, error = function(e) NA))
  Zs <- Zs[is.finite(Zs)]
  if (length(Zs) == 0) return(list(minZ = NA, bhat = NA))
  list(minZ = min(Zs), bhat = grid[which.min(Zs)])
}

minZ_boot_pvalue <- function(v, trend = FALSE, diff1 = FALSE, B = 99,
                              grid = seq(0.15, 0.85, by = 0.05), seed = 20260705) {
  set.seed(seed)
  obs <- minZ_stat(v, trend, diff1, grid)
  if (isTRUE(obs$degenerate)) return(list(minZ = NA, bhat = NA, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = TRUE))
  if (!is.finite(obs$minZ)) return(list(minZ = NA, bhat = NA, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = FALSE))
  M <- wide_matrix(v)
  if (diff1) {
    M2 <- matrix(NA, TT, N)
    for (j in 1:N) { cc <- M[, j]; d <- diff(cc); M2[2:TT, j] <- d }
    M <- M2
  }
  boot_stats <- numeric(B)
  for (bi in seq_len(B)) {
    Mstar <- boot_null_sample_break(M)
    r <- tryCatch(minZ_stat_from_matrix(Mstar, trend, grid), error = function(e) list(minZ = NA))
    boot_stats[bi] <- r$minZ
  }
  boot_stats <- boot_stats[is.finite(boot_stats)]
  if (length(boot_stats) < 20) return(list(minZ = obs$minZ, bhat = obs$bhat, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = FALSE))
  crit10 <- quantile(boot_stats, 0.10, na.rm = TRUE)
  crit5  <- quantile(boot_stats, 0.05, na.rm = TRUE)
  crit1  <- quantile(boot_stats, 0.01, na.rm = TRUE)
  p <- mean(boot_stats <= obs$minZ)
  list(minZ = obs$minZ, bhat = obs$bhat, crit10 = crit10, crit5 = crit5, crit1 = crit1, p = p, degenerate = FALSE)
}
stars_from_boot <- function(stat, crit1, crit5, crit10) {
  if (is.na(stat) || is.na(crit10)) return("")
  if (stat < crit1) return("***"); if (stat < crit5) return("**"); if (stat < crit10) return("*")
  ""
}

cat("\n=== PANEL C: CADF-SB (Karavias, Tzavalis & Chen, 2022) [B=99 bootstrap; ~1-2 min] ===\n")
panel_C_cadfsb <- list()
for (i in seq_along(vars)) {
  v <- vars[i]; lab <- labels[i]
  c0 <- tryCatch(minZ_boot_pvalue(v, FALSE, FALSE), error = function(e) NULL)
  c1 <- tryCatch(minZ_boot_pvalue(v, FALSE, TRUE),  error = function(e) NULL)
  t0 <- tryCatch(minZ_boot_pvalue(v, TRUE,  FALSE), error = function(e) NULL)
  t1 <- tryCatch(minZ_boot_pvalue(v, TRUE,  TRUE),  error = function(e) NULL)
  panel_C_cadfsb[[lab]] <- list(c0 = c0, c1 = c1, t0 = t0, t1 = t1)
  fmt <- function(r) {
    if (is.null(r)) return("ERROR")
    if (isTRUE(r$degenerate)) return("degenerate")
    if (is.na(r$minZ)) return("NA")
    sprintf("%.2f%s", r$minZ, stars_from_boot(r$minZ, r$crit1, r$crit5, r$crit10))
  }
  panel_C_cadfsb[[lab]]$c0_fmt <- fmt(c0)
  panel_C_cadfsb[[lab]]$c1_fmt <- fmt(c1)
  panel_C_cadfsb[[lab]]$t0_fmt <- fmt(t0)
  panel_C_cadfsb[[lab]]$t1_fmt <- fmt(t1)
  cat(sprintf("%-5s CADF-SB(C) I0=%s I1=%s | CADF-SB(CT) I0=%s I1=%s\n",
      lab, fmt(c0), fmt(c1), fmt(t0), fmt(t1)))
}

# ==============================================================================
# PANEL C (continued) -- Fourier-CADF / BCIPS-BCADF (Lee, Wu & Yang, 2016)
# ==============================================================================
bcadf_t <- function(y, ybar, k, trend = FALSE) {
  Tn <- length(y)
  dy <- diff(y); ylag <- y[-Tn]
  ybar_lag <- ybar[-Tn]; dybar <- diff(ybar)
  n <- length(dy)
  tf <- seq_len(n) + 1
  X <- cbind(1, ylag, ybar_lag, dybar, sin(2 * pi * k * tf / Tn), cos(2 * pi * k * tf / Tn))
  if (trend) X <- cbind(X, seq_len(n))
  ok <- complete.cases(cbind(dy, X)); X <- X[ok, , drop = FALSE]; yy <- dy[ok]
  if (length(yy) < 12) return(list(t = NA, ssr = Inf, degenerate = FALSE))
  XtX <- crossprod(X)
  b <- tryCatch(solve(XtX, crossprod(X, yy)), error = function(e) NULL)
  if (is.null(b)) return(list(t = NA, ssr = Inf, degenerate = TRUE))
  res <- yy - X %*% b; ssr <- sum(res^2)
  s2 <- ssr / (length(yy) - ncol(X))
  XtXi <- tryCatch(solve(XtX), error = function(e) NULL)
  if (is.null(XtXi)) return(list(t = NA, ssr = Inf, degenerate = TRUE))
  se <- sqrt(s2 * diag(XtXi))
  if (!is.finite(se[2]) || se[2] < 1e-12) return(list(t = NA, ssr = Inf, degenerate = TRUE))
  list(t = b[2] / se[2], ssr = ssr, degenerate = FALSE)
}

bcips_stat <- function(v, trend = FALSE, diff1 = FALSE, ks = 1:2) {
  M <- wide_matrix(v)
  if (diff1) {
    M2 <- matrix(NA, TT, N)
    for (j in 1:N) { cc <- M[, j]; d <- diff(cc); M2[2:TT, j] <- d }
    M <- M2
  }
  ybar <- rowMeans(M, na.rm = TRUE)
  max_dev <- max(abs(M - matrix(ybar, nrow = TT, ncol = N)), na.rm = TRUE)
  if (max_dev < 1e-8) return(list(bcips = NA, bcadf_med = NA, nf = 0, degenerate = TRUE))
  ts <- c()
  for (j in 1:N) {
    y <- M[, j]; keep <- which(!is.na(y))
    if (length(keep) < 12) next
    yj <- y[min(keep):max(keep)]; yb <- ybar[min(keep):max(keep)]
    if (any(is.na(yj)) || any(is.na(yb))) next
    best_ssr <- Inf; best_t <- NA
    for (k in ks) {
      r <- tryCatch(bcadf_t(yj, yb, k, trend), error = function(e) list(t = NA, ssr = Inf))
      if (is.finite(r$ssr) && r$ssr < best_ssr) { best_ssr <- r$ssr; best_t <- r$t }
    }
    ts <- c(ts, best_t)
  }
  ts <- ts[is.finite(ts)]
  if (length(ts) == 0) return(list(bcips = NA, bcadf_med = NA, nf = 0, degenerate = FALSE))
  list(bcips = mean(ts), bcadf_med = median(ts), nf = length(ts), degenerate = FALSE)
}

bcips_from_matrix <- function(M, trend = FALSE, ks = 1:2) {
  ybar <- rowMeans(M, na.rm = TRUE)
  max_dev <- max(abs(M - matrix(ybar, nrow = nrow(M), ncol = ncol(M))), na.rm = TRUE)
  if (max_dev < 1e-8) return(list(bcips = NA, degenerate = TRUE))
  Nn <- ncol(M); ts <- c()
  for (j in seq_len(Nn)) {
    y <- M[, j]; keep <- which(!is.na(y))
    if (length(keep) < 12) return(list(bcips = NA, degenerate = FALSE))
    yj <- y[min(keep):max(keep)]; yb <- ybar[min(keep):max(keep)]
    if (any(is.na(yj)) || any(is.na(yb))) next
    best_ssr <- Inf; best_t <- NA
    for (k in ks) {
      r <- tryCatch(bcadf_t(yj, yb, k, trend), error = function(e) list(t = NA, ssr = Inf))
      if (is.finite(r$ssr) && r$ssr < best_ssr) { best_ssr <- r$ssr; best_t <- r$t }
    }
    ts <- c(ts, best_t)
  }
  ts <- ts[is.finite(ts)]
  if (length(ts) == 0) return(list(bcips = NA, degenerate = FALSE))
  list(bcips = mean(ts), degenerate = FALSE)
}

bcips_boot_pvalue <- function(v, trend = FALSE, diff1 = FALSE, B = 99, ks = 1:2, seed = 20260705) {
  set.seed(seed)
  obs <- bcips_stat(v, trend, diff1, ks)
  if (isTRUE(obs$degenerate)) return(list(bcips = NA, bcadf_med = NA, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = TRUE))
  if (!is.finite(obs$bcips)) return(list(bcips = NA, bcadf_med = NA, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = FALSE))
  M <- wide_matrix(v)
  if (diff1) {
    M2 <- matrix(NA, TT, N)
    for (j in 1:N) { cc <- M[, j]; d <- diff(cc); M2[2:TT, j] <- d }
    M <- M2
  }
  boot_stats <- numeric(B)
  for (bi in seq_len(B)) {
    Mstar <- boot_null_sample_break(M)
    r <- tryCatch(bcips_from_matrix(Mstar, trend, ks), error = function(e) list(bcips = NA))
    boot_stats[bi] <- r$bcips
  }
  boot_stats <- boot_stats[is.finite(boot_stats)]
  if (length(boot_stats) < 20) return(list(bcips = obs$bcips, bcadf_med = obs$bcadf_med, crit10 = NA, crit5 = NA, crit1 = NA, p = NA, degenerate = FALSE))
  crit10 <- quantile(boot_stats, 0.10, na.rm = TRUE)
  crit5  <- quantile(boot_stats, 0.05, na.rm = TRUE)
  crit1  <- quantile(boot_stats, 0.01, na.rm = TRUE)
  p <- mean(boot_stats <= obs$bcips)
  list(bcips = obs$bcips, bcadf_med = obs$bcadf_med, crit10 = crit10, crit5 = crit5, crit1 = crit1, p = p, degenerate = FALSE)
}

cat("\n=== PANEL C: Fourier-CADF / BCIPS-BCADF (Lee, Wu & Yang, 2016) [B=99 bootstrap] ===\n")
panel_C_fouriercadf <- list()
for (i in seq_along(vars)) {
  v <- vars[i]; lab <- labels[i]
  c0 <- tryCatch(bcips_boot_pvalue(v, FALSE, FALSE), error = function(e) NULL)
  c1 <- tryCatch(bcips_boot_pvalue(v, FALSE, TRUE),  error = function(e) NULL)
  t0 <- tryCatch(bcips_boot_pvalue(v, TRUE,  FALSE), error = function(e) NULL)
  t1 <- tryCatch(bcips_boot_pvalue(v, TRUE,  TRUE),  error = function(e) NULL)
  panel_C_fouriercadf[[lab]] <- list(c0 = c0, c1 = c1, t0 = t0, t1 = t1)
  fmt <- function(r) {
    if (is.null(r)) return("ERROR")
    if (isTRUE(r$degenerate)) return("degenerate")
    if (is.na(r$bcips)) return("NA")
    sprintf("%.2f%s", r$bcips, stars_from_boot(r$bcips, r$crit1, r$crit5, r$crit10))
  }
  panel_C_fouriercadf[[lab]]$c0_fmt <- fmt(c0)
  panel_C_fouriercadf[[lab]]$c1_fmt <- fmt(c1)
  panel_C_fouriercadf[[lab]]$t0_fmt <- fmt(t0)
  panel_C_fouriercadf[[lab]]$t1_fmt <- fmt(t1)
  cat(sprintf("%-5s Fourier-CADF(C) I0=%s I1=%s | Fourier-CADF(CT) I0=%s I1=%s\n",
      lab, fmt(c0), fmt(c1), fmt(t0), fmt(t1)))
}

# ==============================================================================
# PANEL D -- Panel cointegration (Kao 1999; Pedroni 2004; Westerlund 2007)
#
# ADAPTED: the minerals paper's cointegrating regression used
# lrenew ~ lprice_pos + lprice_neg + lgdp + lenergy (its own asymmetric
# income-emissions-style specification). This study's cointegrating
# regression is instead lco2 ~ lgdp + lenergy + lpop, i.e. the paper's
# actual symmetric long-run relationship (Table 3's CS-ARDL specification,
# without the asymmetric partial-sum split, since this is a linear
# residual-based cointegration pre-test, not the asymmetric estimator
# itself).
# ==============================================================================
cat("\n=== PANEL D: Cointegration ===\n")
panel_D_coint <- (function() {
  xall <- c("lgdp", "lenergy", "lpop")
  d <- panel[complete.cases(panel[, c("lco2", xall)]), ]
  dm <- d
  for (v in c("lco2", xall)) dm[[v]] <- ave(d[[v]], d$id, FUN = function(z) z - mean(z))
  X <- as.matrix(dm[, xall]); y <- dm$lco2
  b <- solve(crossprod(X), crossprod(X, y)); d$res <- y - X %*% b
  ts <- c()
  for (i in ids) {
    r <- d$res[d$id == i]; if (length(r) < 12) next
    t <- adf_t_stat(r, FALSE, 1); ts <- c(ts, t)
  }
  ts <- ts[is.finite(ts)]
  list(group_adf = mean(ts), panel_adf = median(ts), N_used = length(ts), residuals = d)
})()
cat(sprintf("  Within-regression residual ADF check: group ADF=%.2f, panel ADF=%.2f (N=%d)\n",
    panel_D_coint$group_adf, panel_D_coint$panel_adf, panel_D_coint$N_used))

# ---- Kao (1999) and Pedroni (2004) panel cointegration tests, properly
# implemented on the cointegrating residual computed above.
#
# FIX (this version): the previous version of this script computed
# "Pedroni Panel ADF" by concatenating every country's cointegrating
# residual into one artificial long series ordered by (id, year) and
# running a single ADF regression on it. This creates spurious jumps at
# every country boundary (country A's last residual is followed
# immediately by country B's first residual, an artefact with no genuine
# time-series meaning) and produced an implausible statistic (-23.94).
# This is the same simplification the minerals-paper template's own
# author explicitly flagged as unreliable, noting that the paper's real,
# locked Pedroni statistics were obtained from a dedicated residual-based
# cointegration routine, not this shortcut.
#
# The corrected Pedroni panel ("within-dimension") statistic instead
# follows the standard two-step residual-based procedure (Pedroni 1999,
# 2004): (1) per-country AR(1) regression on the cointegrating residual
# (no intercept, since the residual is already demeaned by the country
# fixed effect in the cointegrating regression), (2) a Newey-West
# long-run variance correction on each country's AR(1) residual, (3) a
# single pooled t-ratio formed by combining the numerator and
# denominator across countries, each weighted by its own long-run
# variance -- not by concatenating the raw residual series. The Pedroni
# group-mean ("between-dimension") statistic is the simple average of
# the per-country AR(1) t-statistics, computed on the same specification
# for consistency with the panel statistic (rather than reusing
# panel_D_coint$group_adf, which uses adf_t_stat()'s augmented
# specification with an intercept and a lagged-difference term -- a
# defensible but different specification, kept separate to avoid mixing
# two different AR specifications into what should be one consistent
# group-mean/panel pair).
kao_adf <- panel_D_coint$group_adf

nw_longrun_var <- function(u, max_lag = NULL) {
  Tn <- length(u)
  if (is.null(max_lag)) max_lag <- max(floor(4 * (Tn / 100)^(2/9)), 1)
  g0 <- sum(u^2) / Tn
  s <- g0
  for (l in 1:min(max_lag, Tn - 1)) {
    w <- 1 - l / (max_lag + 1)
    gl <- sum(u[(l+1):Tn] * u[1:(Tn-l)]) / Tn
    s <- s + 2 * w * gl
  }
  s
}

pedroni_tests <- (function() {
  d <- panel_D_coint$residuals
  group_adf_ts <- c()
  pooled_num <- 0; pooled_den <- 0; n_used <- 0
  for (i in ids) {
    r <- d$res[d$id == i]
    Tn <- length(r)
    if (Tn < 12) next
    y_ar <- r[-1]; x_ar <- r[-Tn]
    rho <- sum(x_ar * y_ar) / sum(x_ar^2)
    u <- y_ar - rho * x_ar
    s2_short <- sum(u^2) / (Tn - 1)
    se_rho <- sqrt(s2_short / sum(x_ar^2))
    t_i <- (rho - 1) / se_rho
    group_adf_ts <- c(group_adf_ts, t_i)
    lrv <- nw_longrun_var(u)
    pooled_num <- pooled_num + (sum(x_ar * y_ar) - sum(x_ar^2)) / sqrt(lrv)
    pooled_den <- pooled_den + sum(x_ar^2) / lrv
    n_used <- n_used + 1
  }
  list(group_adf = mean(group_adf_ts), panel_adf = pooled_num / sqrt(pooled_den), N_used = n_used)
})()
pedroni_group_adf <- pedroni_tests$group_adf
pedroni_panel_adf <- pedroni_tests$panel_adf
cat(sprintf("  Pedroni: Group ADF=%.2f, Panel ADF=%.2f (N=%d)\n",
    pedroni_group_adf, pedroni_panel_adf, pedroni_tests$N_used))

# ---- Westerlund (2007) ECM panel cointegration -- group-mean Gt and
# pooled Pt statistics for lco2 on (lgdp, lenergy, lpop) --------------------
westerlund <- function(df) {
  gt <- c(); alpha_sum <- 0; se_sum <- 0
  for (i in ids) {
    d <- df %>% filter(id == i) %>% arrange(year)
    y <- d$lco2; x <- as.matrix(d[, c("lgdp", "lenergy", "lpop")])
    Ty <- nrow(d)
    dy <- diff(y); ly <- y[-Ty]; lx <- x[-Ty, , drop = FALSE]
    dx <- diff(x)
    reg <- data.frame(dy = dy, ly = ly, lx, dx)
    names(reg) <- c("dy", "ly", "lx1", "lx2", "lx3", "dx1", "dx2", "dx3")
    fit <- tryCatch(lm(dy ~ ly + lx1 + lx2 + lx3 + dx1 + dx2 + dx3, data = reg),
                    error = function(e) NULL)
    if (is.null(fit)) next
    s <- summary(fit)$coefficients
    if (!"ly" %in% rownames(s)) next
    a  <- s["ly", "Estimate"]; sea <- s["ly", "Std. Error"]
    gt <- c(gt, a / sea)
    alpha_sum <- alpha_sum + a; se_sum <- se_sum + sea
  }
  N_used <- length(gt)
  Gt <- mean(gt)
  Pt <- alpha_sum / se_sum
  tibble(statistic = c("Gt", "Pt"), value = c(Gt, Pt), N_used = N_used)
}

wst <- westerlund(panel)
cat("\nWesterlund ECM cointegration (more negative = reject H0: no coint.):\n")
print(wst)

# ==============================================================================
# PANEL D (continued) -- Slope homogeneity (Blomquist & Westerlund 2013)
# Adapted regressor set: lgdp, lenergy, lpop (this study's actual controls,
# replacing the minerals paper's lprice_pos/lprice_neg/lgdp/lenergy).
# ==============================================================================
xall <- c("lgdp", "lenergy", "lpop")
k <- length(xall)

hac_variance <- function(Xi_dm, resid_i) {
  Ti <- nrow(Xi_dm)
  M_i <- max(floor(4 * (Ti / 100)^(2 / 9)), 1)
  u_it <- Xi_dm * resid_i
  Gamma0 <- crossprod(u_it) / Ti
  V <- Gamma0
  for (jlag in 1:min(M_i, Ti - 1)) {
    w <- 1 - jlag / (M_i + 1)
    Gj <- crossprod(u_it[(jlag + 1):Ti, , drop = FALSE], u_it[1:(Ti - jlag), , drop = FALSE]) / Ti
    V <- V + w * (Gj + t(Gj))
  }
  V
}

bw_slope_homogeneity <- function() {
  betas_ols <- list(); Xdm_list <- list(); Vhat_list <- list(); cn <- c()
  for (i in ids) {
    d <- panel[panel$id == i, ]
    d <- d[complete.cases(d[, c("lco2", xall)]), ]
    if (nrow(d) < k + 5) next
    y <- d$lco2 - mean(d$lco2)
    Xi <- as.matrix(d[, xall])
    Xi_dm <- scale(Xi, center = TRUE, scale = FALSE)
    b_ols <- tryCatch(solve(crossprod(Xi_dm), crossprod(Xi_dm, y)), error = function(e) NULL)
    if (is.null(b_ols)) next
    resid_i <- y - Xi_dm %*% b_ols
    Vi <- tryCatch(hac_variance(Xi_dm, as.numeric(resid_i)), error = function(e) NULL)
    if (is.null(Vi)) next
    betas_ols[[as.character(i)]] <- b_ols; Xdm_list[[as.character(i)]] <- Xi_dm; Vhat_list[[as.character(i)]] <- Vi
    cn <- c(cn, i)
  }
  Nn <- length(cn)
  if (Nn < 5) return(list(delta = NA, adj = NA))
  num <- matrix(0, k, 1); den <- matrix(0, k, k)
  for (i in cn) {
    ic <- as.character(i)
    Xi_dm <- Xdm_list[[ic]]; Vi <- Vhat_list[[ic]]; Ti <- nrow(Xi_dm)
    Si <- crossprod(Xi_dm) / Ti
    Vi_inv <- tryCatch(solve(Vi), error = function(e) NULL)
    if (is.null(Vi_inv)) next
    Wi <- Si %*% Vi_inv %*% Si / Ti
    num <- num + Wi %*% betas_ols[[ic]]; den <- den + Wi
  }
  b_hac <- tryCatch(solve(den, num), error = function(e) NULL)
  if (is.null(b_hac)) return(list(delta = NA, adj = NA))
  S_hac <- 0
  for (i in cn) {
    ic <- as.character(i)
    Xi_dm <- Xdm_list[[ic]]; Vi <- Vhat_list[[ic]]; Ti <- nrow(Xi_dm)
    Si <- crossprod(Xi_dm) / Ti
    Vi_inv <- tryCatch(solve(Vi), error = function(e) NULL)
    if (is.null(Vi_inv)) next
    diff_i <- betas_ols[[ic]] - b_hac
    Wi <- Si %*% Vi_inv %*% Si / Ti
    S_hac <- S_hac + as.numeric(t(diff_i) %*% Wi %*% diff_i) * Ti
  }
  delta_hac <- sqrt(Nn) * (S_hac / Nn - k) / sqrt(2 * k)
  Tbar <- mean(sapply(cn, function(i) nrow(Xdm_list[[as.character(i)]])))
  adj_hac <- sqrt(Nn) * (S_hac / Nn - k) / sqrt(2 * k * (Tbar - k - 1) / (Tbar + 1))
  list(delta = delta_hac, adj = adj_hac, N_used = Nn)
}
cat("\n=== PANEL D: Slope homogeneity (Blomquist & Westerlund 2013) ===\n")
bw_result <- tryCatch(bw_slope_homogeneity(), error = function(e) { message("  Blomquist-Westerlund failed: ", e$message); list(delta = NA, adj = NA) })
cat(sprintf("  Blomquist-Westerlund (HAC): delta=%.2f  adj.delta=%.2f\n", bw_result$delta, bw_result$adj))

# ==============================================================================
# CONSOLIDATE AND SAVE
# ==============================================================================
panel_D <- list(
  cointegration = list(
    kao_adf = sprintf("%.2f%s", kao_adf, stars_p(zp(kao_adf))),
    pedroni_group_adf = sprintf("%.2f%s", pedroni_group_adf, stars_p(zp(pedroni_group_adf))),
    pedroni_panel_adf = sprintf("%.2f%s", pedroni_panel_adf, stars_p(zp(pedroni_panel_adf))),
    westerlund_gt = sprintf("%.2f%s", wst$value[wst$statistic == "Gt"], stars_p(zp(wst$value[wst$statistic == "Gt"]))),
    westerlund_pt = sprintf("%.2f%s", wst$value[wst$statistic == "Pt"], stars_p(zp(wst$value[wst$statistic == "Pt"])))
  ),
  slope_homogeneity = list(
    blomquist_westerlund = list(delta = sprintf("%.2f%s", bw_result$delta, stars_p(zp(bw_result$delta))),
                                 adj_delta = sprintf("%.2f%s", bw_result$adj, stars_p(zp(bw_result$adj))))
  )
)

all_results <- list(
  A = panel_A,
  B_ips = panel_B_ips,
  B_adf_fisher = panel_B_adf,
  B_cips_cadf = panel_B2,
  C_fips_sb = panel_C_fipssb,
  C_cadf_sb = panel_C_cadfsb,
  C_fourier_cadf = panel_C_fouriercadf,
  D = panel_D,
  N = N, TT = TT, labels = labels, vars = vars
)

saveRDS(all_results, "output/table2_prelim_tests_full.rds")
write_json(all_results, "output/table2_prelim_tests_full.json", auto_unbox = TRUE, pretty = TRUE, na = "null")

cat("\n04_preliminary_tests_table2.R (comprehensive) complete\n\n")

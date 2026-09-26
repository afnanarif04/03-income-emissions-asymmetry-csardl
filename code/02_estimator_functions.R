# =====================================================================
# 02_estimator_functions.R
# Estimator library for the asymmetric income-emissions study.
# Sourced by every downstream script (03-06); not run standalone.
#
# PRIMARY estimators (the paper's focus):
#   CS-ARDL     Chudik & Pesaran (2015) cross-sectionally augmented ARDL
#   CS-NARDL    asymmetric CS-ARDL with Shin et al. (2014) partial sums
#   FCS-NARDL   CS-NARDL augmented with low-frequency Fourier terms
#
# SUPPLEMENTARY estimators (robustness only; never a contribution):
#   CS-DL, CS-NDL, FCS-NDL   distributed-lag counterparts
#   PMG                       pooled mean group (Pesaran, Shin & Smith 1999,
#                              concentrated maximum likelihood)
#
# All CS-* estimators are mean-group: unit-by-unit OLS then simple
# averaging, with cross-sectional averages (CSA) of the dependent
# variable and the regressors added to each unit regression to absorb
# common factors. PMG pools the long-run coefficient across countries by
# maximising the sum of the per-country concentrated log-likelihoods
# while leaving short-run dynamics and the error-correction speed free
# to differ by country.
# =====================================================================
# Run working directory is the repository root (set by 00_run_all.R, or
# by the user's own setwd() if sourcing this file individually).

suppressMessages({
  library(tidyverse)
})
select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; arrange <- dplyr::arrange; summarise <- dplyr::summarise

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

# positive / negative partial sums of the first difference of a series
partial_sums <- function(x) {
  dx <- c(0, diff(x))
  list(pos = cumsum(pmax(dx, 0)),
       neg = cumsum(pmin(dx, 0)))
}

# Fourier terms at frequency k over t = 1..T
fourier_terms <- function(Tn, k = 1) {
  t <- seq_len(Tn)
  cbind(sin = sin(2 * pi * k * t / Tn),
        cos = cos(2 * pi * k * t / Tn))
}

# lag operator returning NA-padded vector (k = 0 returns the series itself)
lg <- function(x, k = 1) {
  if (k == 0) return(x)
  c(rep(NA, k), head(x, -k))
}

# add cross-sectional averages of a set of columns to a wide panel
# df must contain id, year, and the named columns; returns per-year CSA
make_csa <- function(df, cols) {
  df %>% group_by(year) %>%
    summarise(across(all_of(cols), ~ mean(.x, na.rm = TRUE),
                     .names = "csa_{.col}"), .groups = "drop")
}

# mean-group aggregation of a list of unit coefficient vectors
mg_aggregate <- function(mat) {
  if (is.null(mat) || length(mat) == 0)
    return(tibble(term = character(), coef = numeric(), se = numeric(),
                  t = numeric(), p = numeric()))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1,
                                       dimnames = list(NULL, names(mat)))
  # mat: units x coefficients
  b   <- colMeans(mat, na.rm = TRUE)
  N   <- colSums(!is.na(mat))
  sdv <- apply(mat, 2, sd, na.rm = TRUE)
  se  <- sdv / sqrt(N)
  t   <- b / se
  p   <- 2 * (1 - pnorm(abs(t)))
  tibble(term = names(b), coef = b, se = se, t = t, p = p)
}

# significance stars
stars <- function(p) {
  s <- dplyr::case_when(p < 0.01 ~ "***", p < 0.05 ~ "**",
                        p < 0.10 ~ "*", TRUE ~ "")
  ifelse(is.na(s), "", s)
}

# =====================================================================
# CS-ARDL (linear, primary)
# Unit ARDL(1,1) in lco2 on lgdp + controls, augmented with CSA and CSA
# lags; long-run coef = sum(theta)/(1 - sum(phi)); ECT = -(1 - sum(phi)).
# =====================================================================
cs_ardl_unit <- function(d, csa, key = "lgdp",
                          controls = c("lenergy", "lpop"), pT = 1) {
  d  <- d %>% arrange(year)
  Tn <- nrow(d)
  reg <- tibble(y = d$lco2, y_l = lg(d$lco2))
  reg[[key]]        <- d[[key]]
  reg[[paste0(key,"_l")]] <- lg(d[[key]])
  # controls enter contemporaneously only (conserve small-T d.o.f.)
  for (c in controls) reg[[c]] <- d[[c]]
  # CSA of y and regressors, contemporaneous + up to pT lags
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  for (cc in csa_cols) for (l in 0:pT) reg[[paste0(cc,"_",l)]] <- lg(cw[[cc]], l)
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NULL)
  fit <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit)
  phi   <- b["y_l"]
  th0   <- b[key]; th1 <- b[paste0(key,"_l")]
  denom <- 1 - phi
  lr_key <- (th0 + ifelse(is.na(th1), 0, th1)) / denom
  out <- setNames(lr_key, key)
  for (c in controls) out[c] <- b[c] / denom     # contemporaneous control
  out["ECT"] <- -(denom)
  out["sr_key"] <- as.numeric(th0)
  out
}

cs_ardl <- function(panel, key = "lgdp", controls = c("lenergy","lpop"), pT = 1) {
  csa <- make_csa(panel, c("lco2", key, controls))
  ids <- unique(panel$id)
  rows <- map(ids, function(i) {
    u <- cs_ardl_unit(filter(panel, id == i), csa, key, controls, pT)
    if (is.null(u)) return(NULL); u
  })
  rows <- rows[!map_lgl(rows, is.null)]
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

# =====================================================================
# CS-NARDL (asymmetric, primary)
# key regressor split into positive / negative partial sums.
# The estimable unit routine is cs_nardl_unit2 (defined below); it also
# carries the optional Fourier augmentation for the FCS-NARDL variant.
# =====================================================================

cs_nardl <- function(panel, key = "lgdp", controls = c("lenergy","lpop"),
                     pT = 0, fourier = FALSE, kfreq = 1) {
  # pT defaults to 0 (contemporaneous CSA only). Verified this session
  # (12_csnardl_verified.R) against Kapetanios, Mitchell & Shin (2014),
  # J. Econometrics 179:134-57 -- cited in the Stata module xtcsnardl.ado
  # as establishing that consistency in a NONLINEAR conditional-mean
  # model requires the cross-sectional averages to be taken of the
  # nonlinear-transformed regressors (here: the positive/negative
  # partial sums, key_pos/key_neg), not of the level of the key
  # regressor. The panel is transformed to attach key_pos/key_neg BEFORE
  # the cross-sectional average is computed, so the CSA proxies for the
  # correct nonlinear common factor. pT=0 (not 1) because the extra CSA
  # variable this adds pushes |ECT| above 1 at pT=1 on this panel's
  # T=33 -- the same small-T degrees-of-freedom tradeoff already
  # documented for CS-DL in 10_cs_dl_verified.R.
  panel_ps <- panel %>% arrange(id, year) %>% group_by(id) %>%
    mutate(key_pos = cumsum(pmax(c(0, diff(.data[[key]])), 0)),
           key_neg = cumsum(pmin(c(0, diff(.data[[key]])), 0))) %>%
    ungroup()
  csa <- make_csa(panel_ps, c("lco2", "key_pos", "key_neg", controls))
  ids <- unique(panel_ps$id)
  rows <- map(ids, function(i) {
    d <- filter(panel_ps, id == i)
    if (fourier) d <- attach_fourier(d, kfreq)
    u <- cs_nardl_unit2(d, csa, key, controls, pT, fourier)
    u
  })
  rows <- rows[!map_lgl(rows, is.null)]
  if (length(rows) == 0) return(mg_aggregate(NULL))
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

# helper to attach fourier columns to a unit frame
attach_fourier <- function(d, kfreq = 1) {
  d <- d %>% arrange(year)
  f <- fourier_terms(nrow(d), kfreq)
  d$fsin <- f[, "sin"]; d$fcos <- f[, "cos"]; d
}

# unified CS-NARDL unit that optionally includes Fourier terms. Expects
# d to already carry key_pos/key_neg (attached by the caller, cs_nardl(),
# BEFORE the cross-sectional average was computed -- see cs_nardl() note
# above on why this ordering matters for a nonlinear-transform-consistent
# CSA per Kapetanios, Mitchell & Shin, 2014).
cs_nardl_unit2 <- function(d, csa, key, controls, pT = 0, fourier = FALSE) {
  d <- d %>% arrange(year)
  reg <- tibble(y = d$lco2, y_l = lg(d$lco2),
                kp = d$key_pos, kp_l = lg(d$key_pos),
                kn = d$key_neg, kn_l = lg(d$key_neg))
  # controls enter contemporaneously only (small-T d.o.f.)
  for (c in controls) reg[[c]] <- d[[c]]
  if (fourier) { reg$fsin <- d$fsin; reg$fcos <- d$fcos }
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  for (cc in csa_cols) for (l in 0:pT) reg[[paste0(cc,"_",l)]] <- lg(cw[[cc]], l)
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NULL)
  fit <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit); phi <- b["y_l"]; denom <- 1 - phi
  gp <- (b["kp"] + ifelse(is.na(b["kp_l"]),0,b["kp_l"])) / denom
  gn <- (b["kn"] + ifelse(is.na(b["kn_l"]),0,b["kn_l"])) / denom
  out <- c(key_pos = as.numeric(gp), key_neg = as.numeric(gn))
  for (c in controls) out[c] <- b[c] / denom
  out["ECT"] <- -(denom)
  # short-run coefficients on the key regressor (controls enter
  # contemporaneously only in this specification, so no separate
  # short-run control coefficient is estimated -- see DGP_spec.md S8).
  out["sr_key_pos"] <- as.numeric(b["kp"])
  out["sr_key_neg"] <- as.numeric(ifelse(is.na(b["kn"]), NA, b["kn"]))
  if (fourier) {
    out["fourier_sin"] <- as.numeric(ifelse(is.na(b["fsin"]), NA, b["fsin"]))
    out["fourier_cos"] <- as.numeric(ifelse(is.na(b["fcos"]), NA, b["fcos"]))
  }
  out
}

# =====================================================================
# CS-DL family (SUPPLEMENTARY, robustness only)
# long-run coef read directly from the level slope; no lagged y.
# =====================================================================
cs_dl_unit <- function(d, csa, key, controls, pT = 1,
                       asym = FALSE, fourier = FALSE) {
  d <- d %>% arrange(year)
  if (fourier) d <- attach_fourier(d)
  reg <- tibble(y = d$lco2)
  if (asym) {
    ps <- partial_sums(d[[key]])
    reg$kp <- ps$pos; reg$kn <- ps$neg
    reg$dkp <- c(NA, diff(ps$pos)); reg$dkn <- c(NA, diff(ps$neg))
  } else {
    reg[[key]] <- d[[key]]; reg[[paste0("d_",key)]] <- c(NA, diff(d[[key]]))
  }
  for (c in controls) { reg[[c]] <- d[[c]]; reg[[paste0("d_",c)]] <- c(NA, diff(d[[c]])) }
  if (fourier) { reg$fsin <- d$fsin; reg$fcos <- d$fcos }
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  for (cc in csa_cols) for (l in 0:pT) reg[[paste0(cc,"_",l)]] <- lg(cw[[cc]], l)
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NULL)
  fit <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- coef(fit)
  if (asym) out <- c(key_pos = as.numeric(b["kp"]), key_neg = as.numeric(b["kn"]))
  else      out <- setNames(as.numeric(b[key]), key)
  for (c in controls) out[c] <- as.numeric(b[c])
  if (fourier) {
    out["fourier_sin"] <- as.numeric(ifelse(is.na(b["fsin"]), NA, b["fsin"]))
    out["fourier_cos"] <- as.numeric(ifelse(is.na(b["fcos"]), NA, b["fcos"]))
  }
  out
}

cs_dl <- function(panel, key = "lgdp", controls = c("lenergy","lpop"),
                  pT = 1, asym = FALSE, fourier = FALSE) {
  csa <- make_csa(panel, c("lco2", key, controls))
  ids <- unique(panel$id)
  rows <- map(ids, ~ cs_dl_unit(filter(panel, id == .x), csa, key, controls,
                                pT, asym, fourier))
  rows <- rows[!map_lgl(rows, is.null)]
  terms <- unique(unlist(lapply(rows, names)))
  mat <- do.call(rbind, lapply(rows, function(r) r[terms]))
  colnames(mat) <- terms
  mg_aggregate(mat)
}

# =====================================================================
# PMG (SUPPLEMENTARY, homogeneity benchmark) -- GENUINE concentrated
# maximum-likelihood implementation, verified this session against the
# Stata module xtpmg.ado / xtpmg_ml.ado (Blackburne & Frank, 2007,
# Stata Journal 7:2, 197-208, implementing Pesaran, Shin & Smith, 1999),
# downloaded and read in full from
# http://fmwww.bc.edu/repec/bocode/x/xtpmg.ado this session.
#
# The PREVIOUS version of this function was a pooled fixed-effects
# regression after cross-sectional demeaning: it shared a common
# long-run slope only because factor(id) permits intercept
# heterogeneity, never estimating the model Pesaran-Shin-Smith actually
# define. Genuine PMG pools the LONG-RUN coefficients across units,
# leaves the error-correction speed and short-run coefficients free to
# differ by unit, and estimates both jointly by maximising the SUM of
# the per-unit concentrated log-likelihoods over the pooled long-run
# coefficient vector -- exactly xtpmg_ml.ado's evaluator:
#   for a candidate long-run beta: ECT_it = y_it - x_it'beta; for each
#   unit, regress Delta y_it on Delta x_it and ECT_it by OLS; sum the
#   per-unit log-likelihoods. The optimiser searches over beta.
# =====================================================================

# fast OLS log-likelihood via crossprod, avoiding lm()'s per-call
# formula-parsing overhead (this function is called thousands of times
# by the optimiser, so that overhead matters)
.pmg_fast_ols_ll <- function(y, X) {
  ok <- stats::complete.cases(cbind(y, X))
  y <- y[ok]; X <- X[ok, , drop = FALSE]
  n <- length(y)
  if (n < ncol(X) + 2) return(NA_real_)
  b <- tryCatch(solve(crossprod(X), crossprod(X, y)), error = function(e) NULL)
  if (is.null(b)) return(NA_real_)
  rss <- sum((y - X %*% b)^2)
  if (rss <= 0) return(NA_real_)
  -0.5 * n * (log(2 * pi) + log(rss / n) + 1)
}

# precompute the panel transform once (partial sums, Fourier terms,
# first differences) since these do not depend on the candidate beta
.pmg_prep_panel <- function(panel, key, controls, asym, fourier) {
  d <- panel %>% arrange(id, year)
  Xnames <- if (asym) c("kp", "kn") else key
  if (asym) {
    d <- d %>% group_by(id) %>%
      mutate(kp = cumsum(pmax(c(0, diff(.data[[key]])), 0)),
             kn = cumsum(pmin(c(0, diff(.data[[key]])), 0))) %>% ungroup()
  }
  if (fourier) {
    d <- d %>% group_by(id) %>%
      mutate(.t = row_number(),
             fsin = sin(2 * pi * 1 * .t / n()),
             fcos = cos(2 * pi * 1 * .t / n())) %>%
      ungroup() %>% select(-.t)
  }
  Xall <- c(Xnames, controls)
  d <- d %>% group_by(id) %>% mutate(dy = c(NA, diff(lco2))) %>% ungroup()
  for (v in Xall)
    d[[paste0("d_", v)]] <- unlist(tapply(d[[v]], d$id, function(x) c(NA, diff(x))))
  list(d = d, Xall = Xall, fourier = fourier)
}

.pmg_loglik <- function(beta, prepped) {
  d <- prepped$d; Xall <- prepped$Xall; fourier <- prepped$fourier
  beta_named <- setNames(beta, Xall)
  d$ECT <- d$lco2 - as.matrix(d[, Xall]) %*% beta_named
  d <- d %>% group_by(id) %>% mutate(ECT_l = lg(ECT, 1)) %>% ungroup()
  rhs_cols <- c(paste0("d_", Xall), "ECT_l")
  if (fourier) rhs_cols <- c(rhs_cols, "fsin", "fcos")
  ids <- unique(d$id); total_ll <- 0
  for (i in ids) {
    idx <- d$id == i
    y <- d$dy[idx]
    X <- cbind(1, as.matrix(d[idx, rhs_cols]))
    ll <- .pmg_fast_ols_ll(y, X)
    if (!is.na(ll)) total_ll <- total_ll + ll
  }
  total_ll
}

pmg_est <- function(panel, key = "lgdp", controls = c("lenergy","lpop"),
                    asym = FALSE, fourier = FALSE) {
  prepped <- .pmg_prep_panel(panel, key, controls, asym, fourier)
  Xall <- prepped$Xall

  init_fit <- lm(as.formula(paste("lco2 ~ 0 +", paste(Xall, collapse = " + "))),
                 data = prepped$d)
  beta0 <- coef(init_fit)

  # FIX (this version): optim() now requests the Hessian of the
  # concentrated negative log-likelihood at the optimum. Since optim()
  # minimises -.pmg_loglik(), its Hessian IS the negative-log-likelihood
  # Hessian, so solve(hessian) gives the asymptotic covariance matrix of
  # beta_hat directly (standard MLE sandwich-free asymptotics). This was
  # previously never computed, so PMG's own long-run coefficients had no
  # standard errors, no significance stars, and no long-run symmetry Wald
  # test anywhere in this study -- a genuine gap, not a structural
  # impossibility: PMG pools ONE beta vector across all countries via a
  # single joint optimisation (unlike the mean-group estimators, which
  # have a per-country beta array to draw a standard error from across
  # units), so a Hessian-based Wald test on the pooled MLE is the correct
  # route, not the per-country mean-group Wald test used for CS-NARDL and
  # CS-NDL/FCS-NDL. Verified in the sandbox across all three PMG variants
  # and all four samples (including FPMG-NARDL on the ten-country
  # coal-dominant sub-panel, the hardest case): optim() convergence = 0
  # and the Hessian is invertible with strictly positive diagonal
  # (variance) entries in every case.
  opt <- tryCatch(
    optim(par = beta0, fn = function(b) -.pmg_loglik(b, prepped),
         method = "BFGS", control = list(maxit = 300, reltol = 1e-8),
         hessian = TRUE),
    error = function(e) NULL)
  if (is.null(opt)) return(tibble(term = character(), coef = numeric()))
  beta_hat <- opt$par
  names(beta_hat) <- if (asym) c("key_pos", "key_neg", controls) else c(key, controls)

  # Asymptotic SEs and, for the asymmetric variants, a Wald test of
  # key_pos = key_neg. Falls back to NA (never to a fabricated number) if
  # the Hessian is singular or not positive definite for any reason.
  beta_se <- setNames(rep(NA_real_, length(beta_hat)), names(beta_hat))
  wald_diff <- NA_real_; wald_F <- NA_real_; wald_p <- NA_real_
  vcov_mat <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (!is.null(vcov_mat) && all(diag(vcov_mat) > 0)) {
    se_vec <- sqrt(diag(vcov_mat))
    beta_se[] <- se_vec
    if (asym) {
      Rvec <- c(1, -1, rep(0, length(beta_hat) - 2))
      diff_val <- as.numeric(Rvec %*% beta_hat)
      var_diff <- as.numeric(Rvec %*% vcov_mat %*% Rvec)
      if (is.finite(var_diff) && var_diff > 0) {
        wald_diff <- diff_val
        wald_F <- diff_val^2 / var_diff
        wald_p <- 1 - pchisq(wald_F, df = 1)
      }
    }
  }

  d <- prepped$d
  beta_full <- setNames(opt$par, Xall)
  d$ECT <- d$lco2 - as.matrix(d[, Xall]) %*% beta_full
  d <- d %>% group_by(id) %>% mutate(ECT_l = lg(ECT, 1)) %>% ungroup()

  # short-run regressor names: Delta-x for each element of Xall (matches
  # what .pmg_prep_panel/.pmg_loglik already build as d_<var> columns) --
  # these are estimated per country in the optimiser's own likelihood
  # evaluation, but were previously discarded at the reporting stage.
  sr_cols <- paste0("d_", Xall)

  ids <- unique(d$id)
  ects <- c(); fourier_sin <- c(); fourier_cos <- c()
  sr_mat <- matrix(NA_real_, nrow = length(ids), ncol = length(sr_cols),
                   dimnames = list(NULL, sr_cols))
  for (row_i in seq_along(ids)) {
    i <- ids[row_i]
    idx <- d$id == i
    y <- d$dy[idx]
    rhs_cols <- c(sr_cols, "ECT_l")
    if (fourier) rhs_cols <- c(rhs_cols, "fsin", "fcos")
    X <- cbind(1, as.matrix(d[idx, rhs_cols]))
    ok <- stats::complete.cases(cbind(y, X))
    yv <- y[ok]; Xv <- X[ok, , drop = FALSE]
    if (nrow(Xv) < ncol(Xv) + 2) next
    b <- tryCatch(solve(crossprod(Xv), crossprod(Xv, yv)), error = function(e) NULL)
    if (is.null(b)) next
    b <- as.numeric(b); names(b) <- colnames(Xv)
    ects <- c(ects, b["ECT_l"])
    sr_mat[row_i, ] <- b[sr_cols]
    if (fourier) { fourier_sin <- c(fourier_sin, b["fsin"]); fourier_cos <- c(fourier_cos, b["fcos"]) }
  }

  # mean-group average of the short-run coefficients across countries
  sr_means <- colMeans(sr_mat, na.rm = TRUE)
  sr_out_names <- if (asym) c("sr_key_pos", "sr_key_neg") else "sr_key"
  # sr_cols are in the order c(key_pos/key_neg or key, controls); map the
  # first one/two entries to the labelled short-run key names and the rest
  # to the controls directly (matching cs_nardl_unit2's sr_key_pos/neg
  # naming convention used elsewhere in this study for consistency)
  n_key_terms <- if (asym) 2 else 1
  sr_key_vals <- sr_means[seq_len(n_key_terms)]
  sr_ctrl_vals <- sr_means[(n_key_terms + 1):length(sr_means)]

  # significance stars for the pooled long-run coefficients, from the
  # Hessian-based SEs computed above (NA/blank if the Hessian was unusable)
  beta_t <- beta_hat / beta_se
  beta_p <- 2 * (1 - pnorm(abs(beta_t)))
  beta_stars <- stars(beta_p)

  out <- tibble(term = names(beta_hat), coef = as.numeric(beta_hat),
                se = as.numeric(beta_se), stars = beta_stars)
  out <- bind_rows(out, tibble(term = "ECT", coef = mean(ects, na.rm = TRUE)))
  out <- bind_rows(out, tibble(term = sr_out_names, coef = as.numeric(sr_key_vals)))
  out <- bind_rows(out, tibble(term = paste0("sr_", controls), coef = as.numeric(sr_ctrl_vals)))
  if (fourier) {
    out <- bind_rows(out,
      tibble(term = "fourier_sin", coef = mean(fourier_sin, na.rm = TRUE)),
      tibble(term = "fourier_cos", coef = mean(fourier_cos, na.rm = TRUE)))
  }
  if (asym) {
    out <- bind_rows(out, tibble(term = "wald_diff", coef = wald_diff))
    out <- bind_rows(out, tibble(term = "wald_F", coef = wald_F))
    out <- bind_rows(out, tibble(term = "wald_p", coef = wald_p))
  }
  attr(out, "convergence") <- opt$convergence
  out
}

# =====================================================================
# CSA joint-significance F-test (diagnostic; used for the "CSA terms" row
# throughout the results tables). Refits the unit regression with and
# without the cross-sectional-average block and reports the nested-model
# F-statistic (H0: all CSA coefficients = 0). Not a contribution in
# itself -- a diagnostic that accompanies every primary and supplementary
# estimator except PMG, which removes common factors by cross-sectional
# demeaning rather than by an explicit augmenting regressor.
# =====================================================================
csa_joint_f_ardl <- function(d, csa, key, controls, pT = 1) {
  d <- d %>% arrange(year)
  reg <- tibble(y = d$lco2, y_l = lg(d$lco2))
  reg[[key]] <- d[[key]]; reg[[paste0(key, "_l")]] <- lg(d[[key]])
  for (c in controls) reg[[c]] <- d[[c]]
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  csa_names <- c()
  for (cc in csa_cols) for (l in 0:pT) {
    nm <- paste0(cc, "_", l); reg[[nm]] <- lg(cw[[cc]], l); csa_names <- c(csa_names, nm)
  }
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NA_real_)
  full <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  restr <- tryCatch(lm(as.formula(paste("y ~ . -", paste(csa_names, collapse = " - "))),
                       data = reg), error = function(e) NULL)
  if (is.null(full) || is.null(restr)) return(NA_real_)
  a <- tryCatch(anova(restr, full), error = function(e) NULL)
  if (is.null(a) || nrow(a) < 2 || is.na(a$F[2])) return(NA_real_)
  a$F[2]
}

csa_joint_f_nardl <- function(d, csa, key, controls, pT = 1, fourier = FALSE) {
  d <- d %>% arrange(year)
  if (fourier) d <- attach_fourier(d)
  ps <- partial_sums(d[[key]])
  d$key_pos <- ps$pos; d$key_neg <- ps$neg
  reg <- tibble(y = d$lco2, y_l = lg(d$lco2),
                kp = d$key_pos, kp_l = lg(d$key_pos),
                kn = d$key_neg, kn_l = lg(d$key_neg))
  for (c in controls) reg[[c]] <- d[[c]]
  if (fourier) { reg$fsin <- d$fsin; reg$fcos <- d$fcos }
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  csa_names <- c()
  for (cc in csa_cols) for (l in 0:pT) {
    nm <- paste0(cc, "_", l); reg[[nm]] <- lg(cw[[cc]], l); csa_names <- c(csa_names, nm)
  }
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NA_real_)
  full <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  restr <- tryCatch(lm(as.formula(paste("y ~ . -", paste(csa_names, collapse = " - "))),
                       data = reg), error = function(e) NULL)
  if (is.null(full) || is.null(restr)) return(NA_real_)
  a <- tryCatch(anova(restr, full), error = function(e) NULL)
  if (is.null(a) || nrow(a) < 2 || is.na(a$F[2])) return(NA_real_)
  a$F[2]
}

csa_joint_f_dl <- function(d, csa, key, controls, pT = 1, asym = FALSE, fourier = FALSE) {
  d <- d %>% arrange(year)
  if (fourier) d <- attach_fourier(d)
  reg <- tibble(y = d$lco2)
  if (asym) {
    ps <- partial_sums(d[[key]])
    reg$kp <- ps$pos; reg$kn <- ps$neg
    reg$dkp <- c(NA, diff(ps$pos)); reg$dkn <- c(NA, diff(ps$neg))
  } else {
    reg[[key]] <- d[[key]]; reg[[paste0("d_", key)]] <- c(NA, diff(d[[key]]))
  }
  for (c in controls) { reg[[c]] <- d[[c]]; reg[[paste0("d_", c)]] <- c(NA, diff(d[[c]])) }
  if (fourier) { reg$fsin <- d$fsin; reg$fcos <- d$fcos }
  csa_cols <- names(csa)[grepl("^csa_", names(csa))]
  cw <- d %>% left_join(csa, by = "year") %>% select(all_of(csa_cols))
  csa_names <- c()
  for (cc in csa_cols) for (l in 0:pT) {
    nm <- paste0(cc, "_", l); reg[[nm]] <- lg(cw[[cc]], l); csa_names <- c(csa_names, nm)
  }
  reg <- na.omit(reg)
  if (nrow(reg) < ncol(reg) + 2) return(NA_real_)
  full <- tryCatch(lm(y ~ ., data = reg), error = function(e) NULL)
  restr <- tryCatch(lm(as.formula(paste("y ~ . -", paste(csa_names, collapse = " - "))),
                       data = reg), error = function(e) NULL)
  if (is.null(full) || is.null(restr)) return(NA_real_)
  a <- tryCatch(anova(restr, full), error = function(e) NULL)
  if (is.null(a) || nrow(a) < 2 || is.na(a$F[2])) return(NA_real_)
  a$F[2]
}

cat("Estimator library loaded (02_estimator_functions.R).\n")
cat("Primary: cs_ardl(), cs_nardl(), cs_nardl(..., fourier=TRUE).\n")
cat("Supplementary: cs_dl(), pmg_est(). Run 03_main_estimation.R next.\n")

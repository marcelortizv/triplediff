library(data.table)
# Function to generate a test panel dataset
# Dataset only for testing purposes and not for validation of any statistical property
generate_test_panel <- function(seed = 123,
                                num_ids = 50,
                                time = 2,
                                initial.year = 2019,
                                treatment.year = 2020) {
  # Set seed for reproducibility
  set.seed(seed)

  # Number of IDs
  num_ids <- num_ids

  # Number of observations per ID
  num_obs <- time

  # Generate IDs
  id = rep(1:num_ids, each = num_obs)

  # Generate x1 and x2 for each state
  x1 <- rep(rnorm(num_ids), each = num_obs)
  x2 <- rep(rnorm(num_ids), each = num_obs)

  # Generate state variable
  state = rep(sample(c("A", "B", "C", "D"), num_ids, replace = TRUE), each = num_obs)

  # generate time variable
  year = rep(seq(initial.year, initial.year+(time-1), 1), times = num_ids)

  # Generate partition variable
  partition = rep(rbinom(num_ids, 1, 0.5), each = num_obs)

  # Generate target group
  s_2 = as.numeric(state %in% c("A", "B"))

  # Generate post treatment variable
  t_2 = as.numeric(year >= treatment.year)

  # TWFE outcome variable
  y = 2 + 5*s_2 - 2 * partition + 3 * t_2 +
    4*s_2*partition + 2*s_2*t_2 + 3*partition*t_2 +
    1*s_2*partition*t_2 + rnorm(num_ids*num_obs)

  # Create a data.table
  dt <- data.table::data.table(
    id = id,
    state = state,
    year = year,
    partition = partition,
    x1 = x1,
    x2 = x2,
    treat = s_2,
    outcome = y
  )
  return(dt)
}

# Shared data generators and fit wrappers for the aggregation, estimation and input-validation tests. Data are built here; no external files.

# Four-cohort design: cohorts S in {0, 2, 3, 4}, eligibility Q, periods 1..5, time-invariant x1, x2.
# Conditional DDD parallel trends hold; effects tau(g, t) = c(1, 4, 8)[g] + 0.5 (t - g) for eligible treated units.
fx_dgp_a <- function(n, seed, variant = c("A", "A2")) {
  variant <- match.arg(variant)
  set.seed(seed)
  S8 <- c(0, 0, 2, 2, 3, 3, 4, 4); Q8 <- c(0, 1, 0, 1, 0, 1, 0, 1)
  B8 <- rbind(c(0, 0), c(.4, -.2), c(-.3, .3), c(.5, .2), c(.2, -.4), c(-.4, .3), c(.3, .3), c(-.2, -.5))
  p0 <- if (variant == "A") c(.125, .125, .21, .09, .125, .125, .05, .15) else c(.075, .075, .30, .20, .075, .075, .05, .15)
  a <- log(p0 / p0[1])
  x1 <- rnorm(n); x2 <- rnorm(n); U <- runif(n)
  E <- exp(matrix(a, n, 8, byrow = TRUE) + outer(x1, B8[, 1]) + outer(x2, B8[, 2]))
  CP <- t(apply(E / rowSums(E), 1, cumsum))
  k <- pmin(rowSums(CP < U) + 1L, 8L); S <- S8[k]; Q <- Q8[k]
  u <- rnorm(n); eps <- matrix(rnorm(5 * n), n, 5)
  iS <- match(S, c(0, 2, 3, 4)); bS <- rbind(c(0, 0), c(.5, 0), c(0, .5), c(-.5, .5)); dS <- c(0, .2, .4, .6)
  eta <- 2 + 0.5 * x1 + 0.5 * S + Q + u
  m <- (x1 - 0.5 * x2) + bS[iS, 1] * x1 + bS[iS, 2] * x2 + Q * (0.3 * x1 - 0.3 * x2) + dS[iS] + 0.3 * Q
  tau <- function(g, t) c(`2` = 1, `3` = 4, `4` = 8)[as.character(g)] + 0.5 * (t - g)
  Y <- sapply(1:5, function(t) { eff <- numeric(n); tr <- Q == 1 & S > 0 & t >= S
    eff[tr] <- tau(S[tr], t); eta + t * m + eps[, t] + eff })
  d <- data.table::data.table(id = rep(seq_len(n), 5), time = rep(1:5, each = n), y = as.vector(Y),
                              S = rep(S, 5), Q = rep(Q, 5), x1 = rep(x1, 5), x2 = rep(x2, 5))
  data.table::setorderv(d, c("id", "time"))
  d[]
}

# Two-period ddd() fit of cohort g against comparison cohort cc on periods p_from < p_to (calendar order).
# Every unit keeps its p_from covariates in both rows. Returns the ATT, the per-unit IF (ids sorted), ids.
fx_two_period <- function(d, g, cc, p_from, p_to, xformla, method, panel = TRUE) {
  s <- data.table::copy(d[d$time %in% c(p_from, p_to) & d$S %in% c(g, cc)])
  data.table::setorderv(s, c("id", "time"))
  xv <- all.vars(xformla)
  if (length(xv) > 0) s[, (xv) := lapply(.SD, function(v) v[1]), by = id, .SDcols = xv]
  s[, G2 := ifelse(S == g, 2, 0)]
  o <- ddd(yname = "y", tname = "time", idname = "id", gname = "G2", pname = "Q", xformla = xformla,
           data = s, est_method = method, panel = panel, inffunc = TRUE)
  list(att = o$ATT, psi = o$att_inf_func, ids = sort(unique(s$id)))
}

# ATT(g, t) rebuilt from two-period fits: evaluated period t, base period b, comparison cohorts comps.
# When t < b (universal-base pre cell), the two-period fit runs in calendar order and its sign is flipped.
# Each comparison's IF is placed on its own units with scale n / n_c; several comparisons are combined by
# optimal GMM weights (colSums(Omega^-1) / sum(Omega^-1)), SE = (n * sum(Omega^-1))^(-1/2).
fx_hand_gmm <- function(d, g, t, b, comps, xformla, method, panel = TRUE) {
  ids <- sort(unique(d$id)); n <- length(ids); sgn <- if (b < t) 1 else -1
  fits <- lapply(comps, function(cc) fx_two_period(d, g, cc, min(b, t), max(b, t), xformla, method, panel))
  M <- sapply(fits, function(f) { v <- numeric(n); v[match(f$ids, ids)] <- sgn * (n / length(f$ids)) * f$psi; v })
  M <- matrix(M, nrow = n)
  atts <- sgn * vapply(fits, function(f) f$att, 0)
  if (length(comps) == 1) return(list(att = atts, se = sqrt(sum(M[, 1]^2)) / n, IF = M[, 1]))
  iO <- solve(stats::cov(M)); w <- colSums(iO) / sum(iO)
  list(att = sum(w * atts), se = sqrt(1 / (n * sum(iO))), IF = drop(M %*% w), w = w)
}

# Central difference of a scalar estimate under the unit-weight perturbation w0 * (1 + s * h).
fx_num_deriv <- function(est_fun, w0, h, s = 1e-4) {
  (est_fun(w0 * (1 + s * h)) - est_fun(w0 * (1 - s * h))) / (2 * s)
}

# Three periods, cohorts S in {0, 3}; covariate x = W (independent noise) in period 1 and Z in periods 2-3.
# fixed = TRUE records x = Z in every period. Truth: ATT(3,1) = 0, ATT(3,3) = 2.
fx_dgp_c <- function(n, seed, fixed = FALSE) {
  set.seed(seed)
  Z <- runif(n, -1, 1); W <- runif(n, -1, 1); p <- plogis(2 * Z)
  A <- rbinom(n, 1, p); Q <- rbinom(n, 1, p)
  e1 <- rnorm(n); e2 <- rnorm(n); e3 <- rnorm(n)
  data.table::data.table(id = rep(seq_len(n), 3), time = rep(1:3, each = n),
                         y = c(e1, Z + e2, 2 * Z + e3 + 2 * A * Q),
                         x = c(if (fixed) Z else W, Z, Z), S = rep(3 * A, 3), Q = rep(Q, 3),
                         W = rep(W, 3), Z = rep(Z, 3))
}

# Four periods, cohorts S in {0, 3, 4}; binary Z drives assignment, W = Z w.p. 3/4 (else -Z);
# x = W in period 1 and Z in periods 2-4. Trends Y1 - Y2 depend on W only.
fx_dgp_d2 <- function(n, seed) {
  set.seed(seed)
  Z <- ifelse(runif(n) < .5, 1, -1); W <- ifelse(runif(n) < .75, Z, -Z); U <- runif(n)
  pp <- rbind(c(.15, .15, .15, .30, .10, .15), c(.25, .25, .20, .05, .15, .10))
  P <- pp[ifelse(Z == 1, 1, 2), , drop = FALSE]
  k <- pmin(rowSums(t(apply(P, 1, cumsum)) < U) + 1L, 6L)
  S <- c(0, 0, 3, 3, 4, 4)[k]; Q <- c(0, 1, 0, 1, 0, 1)[k]; e <- matrix(rnorm(4 * n), n, 4)
  tr3 <- S == 3 & Q == 1; tr4 <- S == 4 & Q == 1
  data.table::data.table(id = rep(seq_len(n), 4), time = rep(1:4, each = n),
                         y = c(W + e[, 1], e[, 2], e[, 3] + 2 * tr3, e[, 4] + 2 * tr3 + 3 * tr4),
                         x = c(W, Z, Z, Z), S = rep(S, 4), Q = rep(Q, 4))
}

# Two-period panel from the package generator (balanced, time-invariant covariates, gname 'state' in {0, 1}).
fx_2p_base <- function() {
  set.seed(505)
  as.data.frame(gen_dgp_2periods(1000, 1)$data)
}
fx_fit_2p <- function(d, ...) {
  ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
      xformla = ~cov1 + cov2, data = d, est_method = "dr", inffunc = TRUE, ...)
}
# Multi-period panel from the package generator (cohorts 0, 2, 3; periods 1-3; column 'cluster').
fx_mp_base <- function() {
  set.seed(1234)
  as.data.frame(gen_dgp_mult_periods(1000, 1)$data)
}
fx_fit_mp <- function(d, cg = "nevertreated", bp = "varying", ...) {
  ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
      xformla = ~cov1 + cov2 + cov3 + cov4, data = d, control_group = cg, base_period = bp,
      est_method = "dr", inffunc = TRUE, ...)
}
# Two-period design with unbalanced units: 10 units lose period 2 (da); 10 other units also lose period 1 (db).
fx_2p_unbalanced <- function(d) {
  ids <- sort(unique(d$id)); set.seed(2)
  A <- sample(ids, 10); B <- sample(setdiff(ids, A), 10)
  da <- d[!(d$id %in% A & d$time == 2), ]
  list(A = A, B = B, da = da, db = da[!(da$id %in% B & da$time == 1), ])
}
# No never-treated group: cohorts S in {1, 2, 3, 4} (cohort 1 is treated in the first period), periods 1-4.
fx_dgp_nonever <- function(n, seed) {
  set.seed(seed)
  S <- sample(c(1, 2, 3, 4), n, TRUE, c(.1, .3, .3, .3)); Q <- rbinom(n, 1, .5); x1 <- rnorm(n); u <- rnorm(n)
  d <- data.table::data.table(id = rep(seq_len(n), 4), time = rep(1:4, each = n), S = rep(S, 4), Q = rep(Q, 4),
                              x1 = rep(x1, 4))
  d[, y := rep(u, 4) + time * (0.5 * x1 + 0.2 * S + 0.3 * Q) + rnorm(.N) + ifelse(Q == 1 & time >= S, 2, 0)]
  d[]
}

# Repeated cross-section and unbalanced-panel designs for the influence function tests (test-rc-inffunc.R).
# Cohorts and eligibility form J = 2K cells chosen by a softmax whose log-odds contain 0.2 (k - 1) + 0.5 q times
# pmax(x1, 0)^2, so a logit that is linear in x1, x2 is misspecified for every comparison. The outcome is linear in
# x1, x2 (correct outcome regressions). The error SD is 4 for ineligible and never-treated rows in the periods before
# the first treated period, and 0.5 elsewhere. Effect: 1 + 0.5 (t - g) for eligible units of cohort g in periods t >= g.
fx_rc_cells <- function(n, ft) {
  K <- length(ft); J <- 2 * K; jj <- seq_len(J)
  cell_k <- rep(seq_len(K), 2); cell_q <- rep(0:1, each = K)
  x1 <- pmin(pmax(rnorm(n), -2), 2); x2 <- pmin(pmax(rnorm(n), -2), 2)
  a <- 0.3 * cos(jj); b1 <- 0.4 * sin(2 * jj); b2 <- 0.2 * cos(3 * jj)
  kap <- 0.2 * (cell_k - 1) + 0.5 * cell_q
  S <- sapply(jj, function(j) a[j] + b1[j] * x1 + b2[j] * x2 + kap[j] * pmax(x1, 0)^2)
  G <- -log(-log(matrix(runif(n * J), n, J)))          # Gumbel-max draw from the softmax probabilities
  j <- max.col(S + G, ties.method = "first")
  list(x1 = x1, x2 = x2, j = j, k = cell_k[j], q = cell_q[j], first = ft[cell_k[j]])
}
fx_rc_y <- function(cl, time, ft) {
  J <- 2 * length(ft); jj <- seq_len(J)
  alpha <- 1 + 0.5 * jj; d1 <- 0.8 + 0.3 * cos(jj); d2 <- -0.5 + 0.2 * sin(jj)
  trend <- time * (0.4 * (cl$k - 1) + 0.5 * cl$q + 0.6 * cl$x1 - 0.3 * cl$x2)
  elig <- cl$first > 0 & cl$q == 1
  tau <- ifelse(elig & time >= cl$first, 1 + 0.5 * (time - cl$first), 0)
  sdv <- ifelse(!elig & time <= min(ft[ft > 0]) - 1, 4, 0.5)
  alpha[cl$j] + d1[cl$j] * cl$x1 + d2[cl$j] * cl$x2 + trend + tau + rnorm(length(time)) * sdv
}
fx_rc_state <- function(cl, T) if (T == 2) as.numeric(cl$first > 0) else cl$first
# Repeated cross-section: each row is a different unit observed in one random period.
fx_dgp_rcs <- function(n, seed, T = 2, ft = c(0, 2)) {
  set.seed(seed); cl <- fx_rc_cells(n, ft)
  time <- sample.int(T, n, replace = TRUE)
  data.frame(id = seq_len(n), time = time, y = fx_rc_y(cl, time, ft),
             state = fx_rc_state(cl, T), partition = cl$q, x1 = cl$x1, x2 = cl$x2)
}
# Panel: every unit in every period (time-invariant covariates); `drop` is the share of rows removed at random.
fx_dgp_panel <- function(n, seed, T = 2, ft = c(0, 2), drop = 0) {
  set.seed(seed); cl <- fx_rc_cells(n, ft)
  d <- do.call(rbind, lapply(seq_len(T), function(t) {
    time <- rep(t, n)
    data.frame(id = seq_len(n), time = time, y = fx_rc_y(cl, time, ft),
               state = fx_rc_state(cl, T), partition = cl$q, x1 = cl$x1, x2 = cl$x2)
  }))
  d <- d[order(d$id, d$time), ]
  if (drop > 0) d <- d[runif(nrow(d)) >= drop, ]
  rownames(d) <- NULL; d
}
# Package subgroup codes: 1 untreated and ineligible, 2 untreated and eligible, 3 treated and ineligible, 4 treated and eligible.
fx_rc_sub <- function(d) ifelse(d$state > 0 & d$partition == 1, 4, ifelse(d$state > 0, 3, ifelse(d$partition == 1, 2, 1)))

# ddd() fit on the designs above. path: "rcs2" two-period panel = FALSE, "unb2" two-period unbalanced panel,
# "pan2" two-period balanced panel, "mp_rcs" / "mp_unb" / "mp_pan" the same with several periods.
fx_rc_fit <- function(d, est = "dr", path = "rcs2", w = NULL, cg = "nevertreated", bp = "varying", xf = ~x1 + x2) {
  multi <- path %in% c("mp_rcs", "mp_unb", "mp_pan")
  suppressWarnings(ddd(yname = "y", tname = "time", idname = "id", gname = "state", pname = "partition",
                       xformla = xf, data = d, control_group = if (multi) cg else NULL,
                       base_period = if (multi) bp else NULL, est_method = est,
                       panel = !(path %in% c("rcs2", "mp_rcs")), allow_unbalanced_panel = path %in% c("unb2", "mp_unb"),
                       weightsname = w, inffunc = TRUE))
}
# Accessors for fx_rc_deriv(): the ATT and the influence function of a two-period fit (cells = NULL) or of the listed
# (g, t) cells of a multiple-period fit.
fx_rc_est <- function(cells = NULL) function(o) {
  if (is.null(cells)) return(as.numeric(o$ATT))
  vapply(seq_len(nrow(cells)), function(i) o$ATT[o$groups == cells$g[i] & o$periods == cells$t[i]], 0)
}
fx_rc_psi <- function(cells = NULL) function(o) {
  if (is.null(cells)) return(matrix(as.numeric(o$att_inf_func), ncol = 1))
  j <- vapply(seq_len(nrow(cells)), function(i) which(o$groups == cells$g[i] & o$periods == cells$t[i]), 0L)
  as.matrix(o$inf_func_mat)[, j, drop = FALSE]
}
# Gap between the numerical derivative of the estimate under unit weights w0 * (1 + s * h) and mean(h * psi), in units
# of the tolerance 1e-6 * max(1, |mean(h * psi)|): a value below 1 passes. Directions h = rnorm(n), seeds 4561 to 4563;
# base weights 1 ("ones") and 0.5 + runif(n) with seed 99 ("rand"); s = 1e-4. Rows of one id share a weight.
# fit_fun(d) must fit with weightsname = "wt".
fx_rc_deriv <- function(d, fit_fun, est = fx_rc_est(), psi = fx_rc_psi(), w0s = c("ones", "rand")) {
  d <- data.table::copy(data.table::as.data.table(d))
  ids <- sort(unique(d$id)); n <- length(ids)
  fit_w <- function(wv) { d[, wt := wv[match(id, ids)]]; fit_fun(d) }
  out <- list()
  for (w0n in w0s) {
    set.seed(99); w0 <- if (w0n == "ones") rep(1, n) else 0.5 + runif(n)
    P <- psi(fit_w(w0)); stopifnot(nrow(P) == n)
    for (k in 1:3) {
      set.seed(4560 + k); h <- rnorm(n)
      cd <- fx_num_deriv(function(wv) est(fit_w(wv)), w0, h)
      hif <- colMeans(h * P)
      out[[length(out) + 1]] <- data.frame(w0 = w0n, dir = k, cell = seq_along(cd), cd = cd, hif = hif,
                                           gap = abs(cd - hif), ratio = abs(cd - hif) / (1e-6 * pmax(1, abs(hif))))
    }
  }
  do.call(rbind, out)
}
# Trimming design: treated rows move 3 to the right in x1, and the first 6 rows of each comparison subgroup are placed
# in the treated region, where the fitted propensity score of a comparison row is above 0.995.
fx_rc_trim_design <- function(d, shift = 3, m = 6) {
  s <- fx_rc_sub(d)
  d$x1[s == 4] <- d$x1[s == 4] + shift
  for (a in 1:3) d$x1[head(which(s == a), m)] <- shift + 0.5 + 0.3 * seq_len(m) / m
  d
}
# Per comparison subgroup (3, 2, 1): comparison rows with fitted propensity score >= 0.995 (trimmed), and comparison
# rows within 2e-3 of 0.995 (a derivative step could flip them). Same logit as the package: subgroup a plus subgroup 4.
fx_rc_trim_count <- function(d, w = rep(1, nrow(d))) {
  s <- fx_rc_sub(d)
  t(sapply(c(3, 2, 1), function(a) {
    k <- s %in% c(a, 4)
    p <- suppressWarnings(fitted(glm(I(s[k] == 4) ~ x1 + x2, binomial, data = d[k, ], weights = w[k])))
    c(trimmed = sum(p >= 0.995 & s[k] == a), near = sum(abs(p - 0.995) < 2e-3 & s[k] == a))
  }))
}

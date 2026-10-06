# Helpers for the regression tests of the 0.3.0 fixes. Data are built here; no external files.

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
fx_two_period <- function(d, g, cc, p_from, p_to, xformla, method) {
  s <- data.table::copy(d[d$time %in% c(p_from, p_to) & d$S %in% c(g, cc)])
  data.table::setorderv(s, c("id", "time"))
  xv <- all.vars(xformla)
  if (length(xv) > 0) s[, (xv) := lapply(.SD, function(v) v[1]), by = id, .SDcols = xv]
  s[, G2 := ifelse(S == g, 2, 0)]
  o <- ddd(yname = "y", tname = "time", idname = "id", gname = "G2", pname = "Q", xformla = xformla,
           data = s, est_method = method, inffunc = TRUE)
  list(att = o$ATT, psi = o$att_inf_func, ids = sort(unique(s$id)))
}

# ATT(g, t) rebuilt from two-period fits: evaluated period t, base period b, comparison cohorts comps.
# When t < b (universal-base pre cell), the two-period fit runs in calendar order and its sign is flipped.
# Each comparison's IF is placed on its own units with scale n / n_c; several comparisons are combined by
# optimal GMM weights (colSums(Omega^-1) / sum(Omega^-1)), SE = (n * sum(Omega^-1))^(-1/2).
fx_hand_gmm <- function(d, g, t, b, comps, xformla, method) {
  ids <- sort(unique(d$id)); n <- length(ids); sgn <- if (b < t) 1 else -1
  fits <- lapply(comps, function(cc) fx_two_period(d, g, cc, min(b, t), max(b, t), xformla, method))
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

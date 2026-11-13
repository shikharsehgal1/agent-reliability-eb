# Bounds on lfdr(x0, k0; G) = Pr(theta <= tau | X = x0, k0 runs) over every mixing
# distribution G on a grid whose binomial marginals satisfy linear constraints
#   lower_g(x) <= int dbinom(x; k_g, p) dG(p) <= upper_g(x),   x = 0..k_g,
# for each run-budget group g.
#
# Two uses, which must not be confused:
#   (A) Identified set. lower = upper = a known marginal (the true one in a
#       population calculation, or the fitted one as a plug-in). The result is the
#       range of lfdr values that are observationally equivalent. It is NOT a
#       confidence interval: it ignores that the marginal is estimated.
#   (B) Inferential interval. lower/upper form a simultaneous confidence box for the
#       marginal (cp_box_constraints). Whenever the box contains the true marginal,
#       the true G is feasible, so the interval contains the true lfdr. Coverage
#       >= 1 - alpha follows from the box's coverage, up to grid discretisation.
#       This is the F-localization idea of Ignatiadis and Wager (2022) specialised
#       to binomial counts with a Bonferroni Clopper-Pearson box.
#
# The ratio objective is linearised with the Charnes-Cooper transform:
# y_j = t * w_j, t >= 0, normalised so int dbinom(x0; k0, p) dG = 1 after scaling.
# Restricting G to a grid gives an inner approximation of the bounds over all G on
# [0, 1]; the grid contains tau and tau + 1e-6 so the closure is approximated.

lfdr_bounds <- function(x0, k0, constraints, grid, tau, return_extremes = FALSE) {
  M <- length(grid)
  b0 <- dbinom(x0, k0, grid)
  rows <- list(c(b0, 0), c(rep(1, M), -1))
  dirs <- c("=", "=")
  rhs <- c(1, 0)
  for (con in constraints) {
    B <- kernel_matrix(0:con$k, rep(con$k, con$k + 1), grid)
    for (r in seq_len(nrow(B))) {
      rows <- c(rows, list(c(B[r, ], -con$upper[r]), c(B[r, ], -con$lower[r])))
      dirs <- c(dirs, "<=", ">=")
      rhs <- c(rhs, 0, 0)
    }
  }
  A <- do.call(rbind, rows)
  obj <- c(b0 * (grid <= tau), 0)
  # lpSolve occasionally reports a numerical failure (status 5) on badly scaled
  # problems; retry with other scaling modes before giving up.
  solve <- function(direction) {
    for (sc in c(196, 4, 0)) {
      s <- lpSolve::lp(direction, obj, A, dirs, rhs, scale = sc)
      if (s$status %in% c(0, 2)) break
    }
    s
  }
  lo <- solve("min")
  hi <- solve("max")
  out <- list(lower = if (lo$status == 0) lo$objval else NA_real_,
              upper = if (hi$status == 0) hi$objval else NA_real_,
              status = c(lo$status, hi$status))
  if (return_extremes) {
    to_w <- function(s) if (s$status == 0) s$solution[1:M] / s$solution[M + 1] else NULL
    out$w_lower <- to_w(lo); out$w_upper <- to_w(hi)
  }
  out
}

# Per-task upper end of the (B) interval: U_i = sup lfdr(x_i, k_i; G) over all G whose
# marginals lie in the confidence box. Because the binomial family has monotone
# likelihood ratio, lfdr(x, k; G) is nonincreasing in x for every G, and so is the
# supremum. For each k we therefore scan x downward from k and stop after the first
# value whose bound exceeds `stop_above`; tasks with smaller x are then also above it
# and get U_i = 1 (a valid, conservative bound). LP failures also give U_i = 1.
robust_lfdr_upper <- function(x, k, constraints, grid, tau, stop_above = 1) {
  U <- rep(1, length(x))
  for (kk in sort(unique(k))) {
    for (x0 in kk:0) {
      idx <- k == kk & x == x0
      if (!any(idx)) next
      u <- lfdr_bounds(x0, kk, constraints, grid, tau)$upper
      if (is.na(u)) u <- 1
      U[idx] <- u
      if (u > stop_above) break
    }
  }
  U
}

# (A) Equality constraints on a known marginal pmf, one entry per run budget.
# A small band (eps) absorbs floating-point error; it is reported, not hidden.
marginal_constraints <- function(marginals, eps = 1e-9) {
  lapply(names(marginals), function(kk) {
    f <- marginals[[kk]]
    list(k = as.integer(kk), lower = pmax(f - eps, 0), upper = pmin(f + eps, 1))
  })
}

# (B) Simultaneous (1 - alpha) confidence box for the marginal pmf of every run
# budget present: Clopper-Pearson intervals for each cell count, Bonferroni-adjusted
# over all cells in all groups. Valid conditional on the group sizes when k_i is
# independent of theta_i and tasks are independent.
cp_box_constraints <- function(x, k, alpha = 0.05) {
  ks <- sort(unique(k))
  n_cells <- sum(ks + 1)
  a <- alpha / n_cells
  lapply(ks, function(kk) {
    xs <- x[k == kk]
    N <- length(xs)
    n <- tabulate(xs + 1, nbins = kk + 1)
    lower <- ifelse(n == 0, 0, qbeta(a / 2, n, N - n + 1))
    upper <- ifelse(n == N, 1, qbeta(1 - a / 2, n + 1, N - n))
    list(k = kk, lower = lower, upper = upper)
  })
}

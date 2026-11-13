# Nonparametric maximum likelihood (Kiefer-Wolfowitz) estimate of G on a grid.
#
# G is restricted to a discrete grid p_1 < ... < p_M in [0, 1] and the weights
# w (on the simplex) maximise sum_i log sum_j w_j dbinom(x_i; k_i, p_j).
# The problem is concave in w. It is solved by EM with SQUAREM acceleration
# (Varadhan and Roland, 2008). Convergence is certified with the standard NPMLE
# optimality condition: with D_j(w) = (1/J) sum_i dbinom(x_i; k_i, p_j) / f_w(x_i),
# w is optimal iff max_j D_j(w) <= 1, and loglik(opt) - loglik(w) <= J * (max_j D_j - 1).
#
# The fitted *marginal* is unique; the fitted G is generally NOT. With at most K
# runs per task, only the first K moments of G enter the likelihood, so which
# optimal G the algorithm returns depends on the grid, the start, and smoothing.

# Grid on [0, 1] containing each tau exactly and a point just above it, so mass can
# sit on either side of the null boundary theta = tau.
make_grid <- function(M = 301, taus = c(0.8, 0.9, 0.95), eps = 1e-6) {
  sort(unique(c(seq(0, 1, length.out = M), taus, pmin(taus + eps, 1))))
}

# Collapse tasks to distinct (x, k) pairs; the likelihood depends only on these.
collapse_counts <- function(x, k) {
  key <- paste(x, k, sep = "/")
  first <- !duplicated(key)
  lev <- key[first]
  list(x = x[first], k = k[first],
       n = as.vector(table(factor(key, levels = lev))),
       index = match(key, lev))
}

kernel_matrix <- function(x, k, grid) {
  outer(seq_along(x), seq_along(grid), function(i, j) dbinom(x[i], k[i], grid[j]))
}

fit_npmle <- function(x, k, grid = make_grid(), method = c("barrier", "em"), ...) {
  method <- match.arg(method)
  if (method == "barrier") fit_npmle_barrier(x, k, grid, ...) else fit_npmle_em(x, k, grid, ...)
}

# Log-barrier interior-point method: maximise
#   sum_u n_u log (L w)_u + mu * sum_j log w_j   subject to sum_j w_j = 1
# by equality-constrained Newton steps, for mu decreasing to mu_final. On the central
# path max_j D_j(w) <= 1 + mu * M / J, so the KKT gap is certified by the final mu.
# For mu > 0 this objective is exactly the MAP under a Dirichlet(1 + mu) prior, so
# `smooth` > 0 stops the path early at mu = smooth. As mu -> 0 the path converges to
# one particular maximiser (the one maximising sum log w over the optimal face).
fit_npmle_barrier <- function(x, k, grid, smooth = 0, gap_target = 1e-7,
                              max_newton = 500) {
  cc <- collapse_counts(x, k)
  J <- sum(cc$n); M <- length(grid)
  L <- kernel_matrix(cc$x, cc$k, grid)
  w <- rep(1 / M, M)
  mu_final <- if (smooth > 0) smooth else gap_target * J / M
  phi <- function(w, mu) sum(cc$n * log(drop(L %*% w))) + mu * sum(log(w))
  mu <- J / M
  newton <- 0
  repeat {
    repeat {
      f <- drop(L %*% w)
      g <- drop(crossprod(L, cc$n / f)) + mu / w
      P <- crossprod(L * sqrt(cc$n) / f) + diag(mu / w^2, M)
      P <- P + diag(1e-12 * max(diag(P)), M)
      R <- tryCatch(chol(P), error = function(e) NULL)
      if (is.null(R)) break
      Pg <- backsolve(R, forwardsolve(t(R), g))
      P1 <- backsolve(R, forwardsolve(t(R), rep(1, M)))
      d <- Pg - P1 * sum(Pg) / sum(P1)
      decrement <- sum(d * g)
      newton <- newton + 1
      if (decrement / 2 < 1e-10 || newton >= max_newton) break
      neg <- d < 0
      s <- if (any(neg)) min(1, 0.99 * min(-w[neg] / d[neg])) else 1
      f0 <- phi(w, mu)
      while (s > 1e-12 && !(phi(w + s * d, mu) >= f0 + 0.25 * s * decrement)) s <- s / 2
      if (s <= 1e-12) break
      w <- w + s * d
    }
    if (mu <= mu_final || newton >= max_newton) break
    mu <- max(mu / 10, mu_final)
  }
  w <- w / sum(w)
  gap <- max(crossprod(L, cc$n / drop(L %*% w))) / J - 1
  list(grid = grid, w = w, loglik = sum(cc$n * log(drop(L %*% w))),
       kkt_gap = gap, iter = newton,
       converged = if (smooth > 0) newton < max_newton else gap < 10 * gap_target,
       smooth = smooth, J = J, method = "barrier")
}

fit_npmle_em <- function(x, k, grid, init = NULL, smooth = 0,
                         tol = 1e-6, max_iter = 20000, accelerate = TRUE) {
  cc <- collapse_counts(x, k)
  J <- sum(cc$n)
  M <- length(grid)
  L <- kernel_matrix(cc$x, cc$k, grid)
  w <- if (is.null(init)) rep(1 / M, M) else init / sum(init)

  # EM for the MAP under a symmetric Dirichlet(1 + smooth) prior on w;
  # smooth = 0 is the unpenalised NPMLE.
  em <- function(w) {
    g <- drop(crossprod(L, cc$n / drop(L %*% w)))
    (w * g + smooth) / (J + M * smooth)
  }
  objective <- function(w) {
    ll <- sum(cc$n * log(drop(L %*% w)))
    if (smooth > 0) ll + smooth * sum(log(w)) else ll
  }
  kkt_gap <- function(w) max(crossprod(L, cc$n / drop(L %*% w))) / J - 1

  obj <- objective(w)
  converged <- FALSE
  it <- 0
  while (it < max_iter) {
    it <- it + 1
    w1 <- em(w); w2 <- em(w1)
    wn <- w2
    if (accelerate) {
      r <- w1 - w; v <- w2 - w1 - r
      sv <- sqrt(sum(v^2))
      if (sv > 0) {
        alpha <- min(-sqrt(sum(r^2)) / sv, -1)
        ws <- w - 2 * alpha * r + alpha^2 * v
        # Keep every grid point alive: EM can never revive an exact zero.
        ws <- pmax(ws, 0) + 1e-15
        ws <- em(ws / sum(ws))
        if (is.finite(objective(ws)) && objective(ws) >= objective(w2)) wn <- ws
      }
    }
    obj_new <- objective(wn)
    w <- wn
    if (smooth == 0) {
      if (kkt_gap(w) < tol) { converged <- TRUE; obj <- obj_new; break }
    } else if (abs(obj_new - obj) < 1e-12 * abs(obj)) {
      converged <- TRUE; obj <- obj_new; break
    }
    obj <- obj_new
  }
  list(grid = grid, w = w, loglik = sum(cc$n * log(drop(L %*% w))),
       kkt_gap = kkt_gap(w), iter = it, converged = converged,
       smooth = smooth, J = J, method = "em")
}

# Plug-in lfdr: Pr(theta <= tau | x, k) under the fitted grid G.
lfdr_npmle <- function(fit, x, k, tau) {
  cc <- collapse_counts(x, k)
  L <- kernel_matrix(cc$x, cc$k, fit$grid)
  num <- drop(L %*% (fit$w * (fit$grid <= tau)))
  den <- drop(L %*% fit$w)
  (num / den)[cc$index]
}

# Fitted marginal pmf of X for tasks with k runs.
npmle_marginal <- function(fit, k) drop(kernel_matrix(0:k, rep(k, k + 1), fit$grid) %*% fit$w)

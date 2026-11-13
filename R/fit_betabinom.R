# Pooled parametric empirical Bayes: G = Beta(a, b) fitted by maximum marginal
# (beta-binomial) likelihood from all tasks at once.
#
# Parametrised by mean mu = a / (a + b) and precision phi = a + b. Safeguards:
#   * optimisation is box-constrained, |logit mu| <= 10 and phi in [0.05, 1e4];
#   * several starting values, best kept;
#   * fits on a box boundary are flagged. phi at its upper bound means the data show
#     no extra-binomial spread: the fitted G is nearly a point mass and the plug-in
#     lfdr becomes a near step function of x. phi at its lower bound means G is
#     pushed toward a two-point distribution on {0, 1}.

fit_betabinom <- function(x, k, mu_logit_bound = 10, log_phi_bounds = log(c(0.05, 1e4))) {
  cc <- collapse_counts(x, k)
  nll <- function(par) {
    mu <- plogis(par[1]); phi <- exp(par[2])
    a <- mu * phi; b <- (1 - mu) * phi
    -sum(cc$n * (lchoose(cc$k, cc$x) + lbeta(a + cc$x, b + cc$k - cc$x) - lbeta(a, b)))
  }
  pbar <- min(max(sum(x) / sum(k), 0.01), 0.99)
  lower <- c(-mu_logit_bound, log_phi_bounds[1])
  upper <- c(mu_logit_bound, log_phi_bounds[2])
  best <- NULL
  for (lphi in log(c(1, 10, 100))) {
    o <- optim(c(qlogis(pbar), lphi), nll, method = "L-BFGS-B", lower = lower, upper = upper)
    if (is.null(best) || o$value < best$value) best <- o
  }
  mu <- plogis(best$par[1]); phi <- exp(best$par[2])
  tol <- 1e-3
  flag <- c(phi_upper = best$par[2] >= upper[2] - tol,
            phi_lower = best$par[2] <= lower[2] + tol,
            mu_boundary = abs(best$par[1]) >= mu_logit_bound - tol)
  list(a = mu * phi, b = (1 - mu) * phi, mu = mu, phi = phi, loglik = -best$value,
       convergence = best$convergence, boundary = any(flag), flags = flag)
}

# Posterior null probability under the fitted Beta(a, b): the posterior of theta
# given x successes in k runs is Beta(a + x, b + k - x).
lfdr_betabinom <- function(fit, x, k, tau) pbeta(tau, fit$a + x, fit$b + k - x)

# Unpooled baseline: the same fixed Beta(1, 1) prior for every task, nothing estimated.
lfdr_unpooled <- function(x, k, tau, a = 1, b = 1) pbeta(tau, a + x, b + k - x)

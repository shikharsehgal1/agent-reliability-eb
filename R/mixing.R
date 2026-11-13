# Mixing distributions G on [0, 1], represented as finite mixtures of Beta
# distributions, and exact oracle quantities under such a G.

beta_mix <- function(w, a, b, label = "") {
  stopifnot(length(w) == length(a), length(a) == length(b),
            all(w > 0), all(a > 0), all(b > 0))
  list(w = w / sum(w), a = a, b = b, label = label)
}

r_mix <- function(n, G) {
  comp <- sample.int(length(G$w), n, replace = TRUE, prob = G$w)
  rbeta(n, G$a[comp], G$b[comp])
}

p_mix <- function(p, G) {
  out <- 0
  for (c in seq_along(G$w)) out <- out + G$w[c] * pbeta(p, G$a[c], G$b[c])
  out
}

d_mix <- function(p, G) {
  out <- 0
  for (c in seq_along(G$w)) out <- out + G$w[c] * dbeta(p, G$a[c], G$b[c])
  out
}

# Exact integrals of the binomial kernel against G. For one Beta(a, b) component,
#   int_0^tau dbinom(x; k, p) dBeta(p; a, b)
#     = choose(k, x) B(a + x, b + k - x) / B(a, b) * pbeta(tau, a + x, b + k - x).
# Returns num = int_[0,tau] and den = int_[0,1] for each (x, k) pair.
mix_kernel_integrals <- function(x, k, G, tau) {
  n <- max(length(x), length(k))
  x <- rep_len(x, n); k <- rep_len(k, n)
  num <- den <- numeric(n)
  for (c in seq_along(G$w)) {
    a2 <- G$a[c] + x
    b2 <- G$b[c] + k - x
    m <- G$w[c] * exp(lchoose(k, x) + lbeta(a2, b2) - lbeta(G$a[c], G$b[c]))
    den <- den + m
    num <- num + m * pbeta(tau, a2, b2)
  }
  list(num = num, den = den)
}

# lfdr(x, k; G) = Pr(theta <= tau | X = x, k runs) under the true G.
oracle_lfdr <- function(x, k, G, tau) {
  I <- mix_kernel_integrals(x, k, G, tau)
  I$num / I$den
}

# Marginal pmf of X for a task with k runs: f(x) = int dbinom(x; k, p) dG(p).
oracle_marginal <- function(k, G) mix_kernel_integrals(0:k, k, G, tau = 1)$den

# Pr(theta > tau) under G: the expected fraction of reliable tasks.
prob_reliable <- function(G, tau) 1 - p_mix(tau, G)

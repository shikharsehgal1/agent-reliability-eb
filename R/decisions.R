# Decision rules. A "discovery" (rejection of H_i: theta_i <= tau) means the
# task is declared reliable.

# Exact one-sided binomial p-value for H_i: theta_i <= tau, evaluated at the
# least-favourable null theta_i = tau: Pr(Bin(k_i, tau) >= x_i). Valid (super-uniform)
# for every theta_i <= tau because the binomial family is stochastically increasing.
pval_exact_binom <- function(x, k, tau) pbinom(x - 1, k, tau, lower.tail = FALSE)

# Benjamini-Hochberg step-up at level q. Under independence (or PRDS) of the
# null p-values, FDR <= q * m0 / m.
reject_bh <- function(p, q) {
  m <- length(p)
  o <- order(p)
  ok <- which(p[o] <= q * seq_len(m) / m)
  rej <- logical(m)
  if (length(ok)) rej[o[seq_len(max(ok))]] <- TRUE
  rej
}

# Benjamini-Yekutieli: BH at q / sum(1/i). Controls FDR under arbitrary dependence.
reject_by <- function(p, q) reject_bh(p, q / sum(1 / seq_along(p)))

# Threshold rule: declare reliable when the (estimated) posterior null probability
# is at most q. Each individual call has posterior error probability <= q *if* the
# lfdr is correct; it is not a frequentist FDR guarantee.
reject_lfdr_threshold <- function(lfdr, q) lfdr <= q

# Step-up rule: reject the largest set of smallest lfdr values whose average is <= q.
# If the lfdr values are correct, the posterior expected FDP of that set is <= q.
reject_lfdr_stepup <- function(lfdr, q) {
  o <- order(lfdr)
  ok <- which(cumsum(lfdr[o]) / seq_along(o) <= q)
  rej <- logical(length(lfdr))
  if (length(ok)) rej[o[seq_len(max(ok))]] <- TRUE
  rej
}

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

# Box-robust rule with a frequentist guarantee. Build a (1 - q/2) confidence box for
# the marginal pmf, take U_i = the largest lfdr(x_i, k_i; G) over every G consistent
# with the box, and call task i reliable when U_i <= q/2 (or use the step-up rule on
# U at level q/2). On the event that the box covers the true marginal, every call has
# true lfdr <= U_i, so E[V / R | X] = sum_called lfdr_i / R <= q/2. The box fails with
# probability <= q/2, where FDP <= 1. Hence FDR <= q under the model
# (theta_i iid ~ G, independent binomial runs, k_i independent of theta_i), for ANY
# G. No assumption about the shape of G and no consistency of an estimator is needed.
# Caveats: G is restricted to a fine grid (an inner approximation of the sup), and the
# guarantee is for FDR averaged over theta ~ G, not conditional on fixed theta's.
reject_robust <- function(x, k, q, tau, grid, rule = c("threshold", "stepup")) {
  rule <- match.arg(rule)
  con <- cp_box_constraints(x, k, alpha = q / 2)
  # For the threshold rule only U_i <= q/2 matters; for step-up, values above the
  # stopping point are set to 1, which can only lose power.
  U <- robust_lfdr_upper(x, k, con, grid, tau, stop_above = q / 2)
  if (rule == "threshold") U <= q / 2 else reject_lfdr_stepup(U, q / 2)
}

# Step-up rule: reject the largest set of smallest lfdr values whose average is <= q.
# If the lfdr values are correct, the posterior expected FDP of that set is <= q.
reject_lfdr_stepup <- function(lfdr, q) {
  o <- order(lfdr)
  ok <- which(cumsum(lfdr[o]) / seq_along(o) <= q)
  rej <- logical(length(lfdr))
  if (length(ok)) rej[o[seq_len(max(ok))]] <- TRUE
  rej
}

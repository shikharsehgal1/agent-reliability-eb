# Simulate one benchmark: J tasks, k_i runs each, X_i successes.
#
# Primary model (dependence = "none"):
#   theta_i iid ~ G,  X_i | theta_i ~ Binomial(k_i, theta_i), independent across tasks.
#
# Stress cases, which violate the primary model on purpose:
#   "within_corr": runs within a task are exchangeable but positively correlated
#                  (intra-task correlation rho), so X_i | theta_i is beta-binomial
#                  with mean k_i * theta_i, not binomial. theta_i keeps its meaning
#                  as the per-run success probability.
#   "clustered":   tasks come in clusters of near-duplicates that share one theta,
#                  so the theta_i are not independent across tasks.

# Run budgets: a single number gives a fixed budget; "het" draws k_i uniformly
# from {3, 5, 10, 20}, independently of theta_i.
het_k_values <- c(3, 5, 10, 20)

draw_k <- function(J, k_design) {
  if (identical(k_design, "het")) sample(het_k_values, J, replace = TRUE)
  else rep(as.integer(k_design), J)
}

simulate_dataset <- function(J, k_design, G, dependence = "none",
                             rho = 0.3, cluster_size = 10) {
  if (dependence == "clustered") {
    n_cl <- ceiling(J / cluster_size)
    theta <- rep(r_mix(n_cl, G), each = cluster_size)[seq_len(J)]
  } else {
    theta <- r_mix(J, G)
  }
  k <- draw_k(J, k_design)
  if (dependence == "within_corr") {
    th <- pmin(pmax(theta, 1e-10), 1 - 1e-10)
    p_run <- rbeta(J, th * (1 - rho) / rho, (1 - th) * (1 - rho) / rho)
    x <- rbinom(J, k, p_run)
  } else {
    x <- rbinom(J, k, theta)
  }
  data.frame(task = seq_len(J), theta = theta, k = k, x = x)
}

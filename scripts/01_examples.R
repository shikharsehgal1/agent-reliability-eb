# Worked examples that need no Monte Carlo:
#   1. how weak all-successes evidence is for theta > 0.9;
#   2. two (three) mixing distributions with the same count distribution but
#      different posterior null probabilities (non-identification);
#   3. population identified sets for lfdr(x | k) under the TRUE marginal, by k.
source("R/load.R")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
tau <- main_tau
grid <- make_grid(M = 601, taus = tau)
G_ex <- G_library$overlap_hi

# ---- 1. All successes ---------------------------------------------------------
ks <- c(1, 2, 3, 5, 10, 15, 20, 22, 25, 29, 30)
all_succ <- do.call(rbind, lapply(ks, function(k) {
  A <- lfdr_bounds(k, k, marginal_constraints(setNames(list(oracle_marginal(k, G_ex)), k)), grid, tau)
  data.frame(k = k, x = k,
             exact_p_value = pval_exact_binom(k, k, tau),
             post_reliable_beta11 = 1 - lfdr_unpooled(k, k, tau),
             post_reliable_oracle_overlap_hi = 1 - oracle_lfdr(k, k, G_ex, tau),
             idset_post_reliable_lower = 1 - A$upper,
             idset_post_reliable_upper = 1 - A$lower)
}))
write.csv(all_succ, "results/tables/all_successes.csv", row.names = FALSE)
cat("Smallest k with all successes giving exact p <= 0.10:", ceiling(log(0.10) / log(tau)),
    "; <= 0.05:", ceiling(log(0.05) / log(tau)), "\n")
print(signif(all_succ, 3))

# ---- 2. Non-identification at k = 3 --------------------------------------------
k <- 3
f_true <- oracle_marginal(k, G_ex)
b <- lfdr_bounds(k, k, marginal_constraints(list("3" = f_true)), grid, tau, return_extremes = TRUE)
grid_lfdr <- function(w, x, k) sum(w * dbinom(x, k, grid) * (grid <= tau)) / sum(w * dbinom(x, k, grid))
Gs <- list(true = NULL, min_lfdr = b$w_lower, max_lfdr = b$w_upper)
nonid <- do.call(rbind, lapply(names(Gs), function(nm) {
  w <- Gs[[nm]]
  f <- if (is.null(w)) f_true else drop(kernel_matrix(0:k, rep(k, k + 1), grid) %*% w)
  lf <- if (is.null(w)) oracle_lfdr(0:k, k, G_ex, tau) else sapply(0:k, function(x) grid_lfdr(w, x, k))
  data.frame(G = nm, x = 0:k, marginal_pmf = f, lfdr = lf,
             prob_reliable = if (is.null(w)) prob_reliable(G_ex, tau) else sum(w[grid > tau]),
             n_atoms = if (is.null(w)) NA else sum(w > 1e-8))
}))
write.csv(nonid, "results/tables/nonidentification_k3.csv", row.names = FALSE)
atoms <- do.call(rbind, lapply(c("min_lfdr", "max_lfdr"), function(nm) {
  w <- Gs[[nm]]; keep <- w > 1e-8
  data.frame(G = nm, atom = grid[keep], weight = w[keep])
}))
write.csv(atoms, "results/tables/nonidentification_k3_atoms.csv", row.names = FALSE)
cat("\nNon-identification, k = 3, G = overlap_hi\n"); print(nonid, digits = 4); print(atoms, digits = 4)
saveRDS(list(grid = grid, Gs = Gs, G_ex = G_ex), "results/nonidentification_k3.rds")

# ---- 3. Population identified sets under the true marginal ---------------------
pop <- do.call(rbind, lapply(names(G_library), function(gn) {
  G <- G_library[[gn]]
  do.call(rbind, lapply(c(3, 5, 10, 20, 50), function(k) {
    con <- marginal_constraints(setNames(list(oracle_marginal(k, G)), k))
    # lfdr_bounds with (x0, k0) = (0, 0) has kernel 1, so it bounds G(tau) itself.
    pi0 <- lfdr_bounds(0, 0, con, grid, tau)
    do.call(rbind, lapply(c(k, k - 1), function(x0) {
      bb <- lfdr_bounds(x0, k, con, grid, tau)
      data.frame(G = gn, k = k, x0 = x0, oracle = oracle_lfdr(x0, k, G, tau),
                 lower = bb$lower, upper = bb$upper, width = bb$upper - bb$lower,
                 status = max(bb$status),
                 prob_reliable = prob_reliable(G, tau),
                 prob_reliable_lower = 1 - pi0$upper, prob_reliable_upper = 1 - pi0$lower)
    }))
  }))
}))
write.csv(pop, "results/tables/population_identified_sets.csv", row.names = FALSE)
cat("\nPopulation identified sets (x0 = k):\n")
print(pop[pop$x0 == pop$k, c("G", "k", "oracle", "lower", "upper", "width", "status",
                              "prob_reliable_lower", "prob_reliable_upper")], digits = 3)

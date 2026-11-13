# Stability of discovery sets conditional on the tasks: theta_1..theta_J are drawn
# once and held fixed; only the run outcomes X are redrawn. Reports the mean pairwise
# Jaccard similarity of discovery sets across redraws, and each task's probability of
# being called reliable as a function of its theta (an operating-characteristic curve).
source("R/load.R")
library(parallel)
ncores <- as.integer(Sys.getenv("NCORES", max(1, detectCores() - 1)))
nredraw <- as.integer(Sys.getenv("NREDRAW", "100"))
tau <- main_tau; q <- main_q; J <- 200
grid <- make_grid()

cases <- expand.grid(G = c("overlap_hi", "separated_hi"), k = c("5", "10", "30", "het"),
                     stringsAsFactors = FALSE)

out <- mclapply(seq_len(nrow(cases)), function(ci) {
  G <- G_library[[cases$G[ci]]]
  set.seed(700000 + ci)
  theta <- r_mix(J, G)
  k <- draw_k(J, cases$k[ci])
  rej <- list()
  for (r in seq_len(nredraw)) {
    set.seed(710000 + 1000 * ci + r)
    x <- rbinom(J, k, theta)
    fit_np <- fit_npmle(x, k, grid = grid)
    fit_bb <- fit_betabinom(x, k)
    rr <- list(
      "oracle:threshold" = reject_lfdr_threshold(oracle_lfdr(x, k, G, tau), q),
      "npmle:threshold" = reject_lfdr_threshold(lfdr_npmle(fit_np, x, k, tau), q),
      "betabinom:threshold" = reject_lfdr_threshold(lfdr_betabinom(fit_bb, x, k, tau), q),
      "unpooled:threshold" = reject_lfdr_threshold(lfdr_unpooled(x, k, tau), q),
      "binom_test:BH" = reject_bh(pval_exact_binom(x, k, tau), q))
    for (m in names(rr)) rej[[m]] <- rbind(rej[[m]], rr[[m]])
  }
  summ <- do.call(rbind, lapply(names(rej), function(m) {
    R <- rej[[m]]
    data.frame(G = cases$G[ci], k = cases$k[ci], rule = m,
               mean_discoveries = mean(rowSums(R)),
               mean_pairwise_jaccard = if (any(R)) mean_pairwise_jaccard(R) else NA,
               mean_false = mean(rowSums(R[, theta <= tau, drop = FALSE])))
  }))
  curve <- do.call(rbind, lapply(names(rej), function(m)
    data.frame(G = cases$G[ci], k = cases$k[ci], rule = m, theta = theta, k_i = k,
               p_called = colMeans(rej[[m]]))))
  list(summary = summ, curve = curve)
}, mc.cores = ncores)

summary_tab <- do.call(rbind, lapply(out, `[[`, "summary"))
curves <- do.call(rbind, lapply(out, `[[`, "curve"))
write.csv(summary_tab, "results/tables/stability.csv", row.names = FALSE)
saveRDS(curves, "results/stability_curves.rds")
print(summary_tab, digits = 3, row.names = FALSE)

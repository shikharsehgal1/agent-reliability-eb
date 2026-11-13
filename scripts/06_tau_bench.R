# Secondary, illustrative analysis of public repeated-run agent data:
# tau-bench historical trajectories (Yao et al., 2024), MIT license,
# github.com/sierra-research/tau-bench at commit 59a200c6d575d595120f1cb70fea53cef0632f6b.
# Files are downloaded by scripts/fetch_tau_bench.sh into data/raw/tau_bench/.
#
# Each record has task_id, trial, and reward in {0, 1}. Trials are separate episodes
# (run.py loops over trials, each a pass over all tasks). Whether trials are
# conditionally independent given the task is an ASSUMPTION: the sampling settings
# used for these files are not recorded. We check one symptom (trial-batch effects).
# The tasks are a fixed benchmark, not a random sample from a population; the
# exchangeable model theta_i ~ G is a working assumption.
source("R/load.R")
dir <- "data/raw/tau_bench"
files <- list.files(dir, pattern = "json$", full.names = TRUE)
if (!length(files)) stop("No tau-bench files; run scripts/fetch_tau_bench.sh first.")
grid <- make_grid(M = 401)
q <- main_q

rows <- list(); batch <- list(); bounds <- list(); counts <- list()
for (f in files) {
  d <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  runs <- data.frame(task = vapply(d, function(r) r$task_id, 0),
                     trial = vapply(d, function(r) r$trial, 0),
                     reward = vapply(d, function(r) r$reward, 0))
  stopifnot(all(runs$reward %in% c(0, 1)), !anyDuplicated(runs[, c("task", "trial")]))
  name <- sub("\\.json$", "", basename(f))
  per_task <- aggregate(reward ~ task, runs, function(v) c(x = sum(v), k = length(v)))
  x <- per_task$reward[, "x"]; k <- per_task$reward[, "k"]
  counts[[name]] <- data.frame(data = name, x = 0:max(k), n_tasks = tabulate(x + 1, max(k) + 1))

  # Trial-batch check: is the success rate the same in every trial pass?
  tab <- table(runs$trial, runs$reward)
  ct <- suppressWarnings(chisq.test(tab))
  batch[[name]] <- data.frame(data = name, trials = nrow(tab),
                              min_trial_rate = min(prop.table(tab, 1)[, "1"]),
                              max_trial_rate = max(prop.table(tab, 1)[, "1"]),
                              chisq_p = ct$p.value)

  fit_np <- fit_npmle(x, k, grid = grid)
  fit_bb <- fit_betabinom(x, k)
  for (tau in c(0.8, 0.9)) {
    lf_np <- lfdr_npmle(fit_np, x, k, tau)
    lf_bb <- lfdr_betabinom(fit_bb, x, k, tau)
    lf_un <- lfdr_unpooled(x, k, tau)
    p <- pval_exact_binom(x, k, tau)
    rows[[length(rows) + 1]] <- data.frame(
      data = name, J = length(x), k = paste(unique(k), collapse = "/"), tau = tau, q = q,
      all_success_tasks = sum(x == k),
      npmle_calls = sum(lf_np <= q), betabinom_calls = sum(lf_bb <= q),
      unpooled_calls = sum(lf_un <= q), BH_calls = sum(reject_bh(p, q)),
      npmle_lfdr_all_success = lfdr_npmle(fit_np, k[1], k[1], tau),
      betabinom_lfdr_all_success = lfdr_betabinom(fit_bb, k[1], k[1], tau),
      unpooled_lfdr_all_success = lfdr_unpooled(k[1], k[1], tau),
      min_p_value = min(p),
      npmle_kkt_gap = fit_np$kkt_gap, bb_a = fit_bb$a, bb_b = fit_bb$b, bb_boundary = fit_bb$boundary)
    k0 <- k[1]
    fitted <- setNames(list(npmle_marginal(fit_np, k0)), k0)
    for (x0 in c(k0, k0 - 1)) {
      A <- lfdr_bounds(x0, k0, marginal_constraints(fitted), grid, tau)
      B <- lfdr_bounds(x0, k0, cp_box_constraints(x, k, 0.05), grid, tau)
      bounds[[length(bounds) + 1]] <- data.frame(
        data = name, tau = tau, x0 = x0, k0 = k0, n_tasks_with_x0 = sum(x == x0),
        npmle_plugin = lfdr_npmle(fit_np, x0, k0, tau),
        betabinom_plugin = lfdr_betabinom(fit_bb, x0, k0, tau),
        A_lower = A$lower, A_upper = A$upper, B_lower = B$lower, B_upper = B$upper,
        lp_status = max(c(A$status, B$status)))
    }
  }
}
out <- list(decisions = do.call(rbind, rows), batch = do.call(rbind, batch),
            bounds = do.call(rbind, bounds), counts = do.call(rbind, counts))
for (nm in names(out)) write.csv(out[[nm]], sprintf("results/tables/tau_bench_%s.csv", nm), row.names = FALSE)
for (nm in names(out)) { cat("\n==", nm, "\n"); print(out[[nm]], digits = 3, row.names = FALSE) }

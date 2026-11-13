# One simulated replicate, and the aggregation of many replicates of a scenario.

lfdr_methods <- c("oracle", "npmle", "betabinom", "unpooled")

run_replicate <- function(scen, seed, grid = make_grid(), do_bounds = TRUE,
                          bounds_alpha = 0.05) {
  set.seed(seed)
  G <- G_library[[scen$G]]
  dat <- simulate_dataset(scen$J, scen$k, G, dependence = scen$dependence)
  x <- dat$x; k <- dat$k; theta <- dat$theta

  fit_np <- fit_npmle(x, k, grid = grid)
  fit_bb <- fit_betabinom(x, k)

  metric_rows <- list()
  calib <- list(); accuracy <- list(); lfdr_main <- list()
  for (tau in analysis_taus) {
    lf <- list(oracle = oracle_lfdr(x, k, G, tau),
               npmle = lfdr_npmle(fit_np, x, k, tau),
               betabinom = lfdr_betabinom(fit_bb, x, k, tau),
               unpooled = lfdr_unpooled(x, k, tau))
    p <- pval_exact_binom(x, k, tau)
    for (q in analysis_qs) {
      rules <- list()
      for (m in lfdr_methods) rules[[paste0(m, ":threshold")]] <- reject_lfdr_threshold(lf[[m]], q)
      for (m in c("oracle", "npmle", "betabinom")) rules[[paste0(m, ":stepup")]] <- reject_lfdr_stepup(lf[[m]], q)
      rules[["binom_test:BH"]] <- reject_bh(p, q)
      rules[["binom_test:BY"]] <- reject_by(p, q)
      for (nm in names(rules)) {
        metric_rows[[length(metric_rows) + 1]] <-
          c(tau = tau, q = q, decision_metrics(rules[[nm]], theta, tau))
      }
      rule_names <- names(rules)
    }
    if (tau == main_tau) {
      lfdr_main <- lf
      for (m in lfdr_methods) {
        calib[[m]] <- calibration_counts(lf[[m]], theta, tau)
        accuracy[[m]] <- c(mae_vs_oracle = mean(abs(lf[[m]] - lf$oracle)),
                           brier = mean((lf[[m]] - (theta <= tau))^2))
      }
    }
  }
  mm <- do.call(rbind, metric_rows)
  metrics <- data.frame(rule = rep(rule_names, length.out = nrow(mm)), mm,
                        row.names = NULL, stringsAsFactors = FALSE)

  bounds <- NULL
  if (do_bounds) {
    k0 <- if (scen$k == "het") 10L else as.integer(scen$k)
    ks <- sort(unique(k))
    fitted <- setNames(lapply(ks, function(kk) npmle_marginal(fit_np, kk)), ks)
    con_A <- marginal_constraints(fitted)
    con_B <- cp_box_constraints(x, k, alpha = bounds_alpha)
    bounds <- do.call(rbind, lapply(c(k0, k0 - 1L), function(x0) {
      bA <- lfdr_bounds(x0, k0, con_A, grid, main_tau)
      bB <- lfdr_bounds(x0, k0, con_B, grid, main_tau)
      data.frame(x0 = x0, k0 = k0,
                 oracle = oracle_lfdr(x0, k0, G, main_tau),
                 npmle = lfdr_npmle(fit_np, x0, k0, main_tau),
                 betabinom = lfdr_betabinom(fit_bb, x0, k0, main_tau),
                 A_lower = bA$lower, A_upper = bA$upper,
                 B_lower = bB$lower, B_upper = bB$upper,
                 A_status = max(bA$status), B_status = max(bB$status))
    }))
  }

  list(metrics = metrics, calib = calib, accuracy = accuracy, bounds = bounds,
       diag = c(npmle_converged = fit_np$converged, npmle_iter = fit_np$iter,
                npmle_kkt_gap = fit_np$kkt_gap, bb_boundary = fit_bb$boundary,
                bb_phi_upper = unname(fit_bb$flags["phi_upper"]),
                bb_convergence = fit_bb$convergence),
       lfdr_main = lfdr_main, data = dat)
}

scenario_seed <- function(id, rep) 100000L * id + rep

run_scenario <- function(scen, nrep, grid = make_grid(), do_bounds = TRUE) {
  reps <- lapply(seq_len(nrep), function(r) {
    out <- run_replicate(scen, scenario_seed(scen$id, r), grid = grid, do_bounds = do_bounds)
    out$lfdr_main <- NULL; out$data <- NULL
    out
  })

  # Operating characteristics: mean and MCSE across replicates.
  all_m <- do.call(rbind, lapply(seq_along(reps), function(r) cbind(rep = r, reps[[r]]$metrics)))
  groups <- split(all_m, list(all_m$rule, all_m$tau, all_m$q), drop = TRUE)
  oc <- do.call(rbind, lapply(groups, function(g) {
    s <- function(v) mc_summary(g[[v]])
    data.frame(rule = g$rule[1], tau = g$tau[1], q = g$q[1],
               FDR = s("FDP")["mean"], FDR_mcse = s("FDP")["mcse"],
               power = s("power")["mean"], power_mcse = s("power")["mcse"],
               mean_R = mean(g$R), mFDR = if (sum(g$R) > 0) sum(g$V) / sum(g$R) else 0,
               P_any_false = mean(g$any_false),
               near_false_per_rep = mean(g$near_false),
               near_false_mcse = s("near_false")["mcse"],
               near_false_rate = if (sum(g$n_near) > 0) sum(g$near_false) / sum(g$n_near) else NA,
               nrep = length(unique(g$rep)), row.names = NULL)
  }))

  calib <- lapply(lfdr_methods, function(m) Reduce(`+`, lapply(reps, function(r) r$calib[[m]])))
  names(calib) <- lfdr_methods
  accuracy <- do.call(rbind, lapply(lfdr_methods, function(m) {
    a <- do.call(rbind, lapply(reps, function(r) r$accuracy[[m]]))
    data.frame(method = m, mae_vs_oracle = mean(a[, 1]), mae_mcse = sd(a[, 1]) / sqrt(nrow(a)),
               brier = mean(a[, 2]), brier_mcse = sd(a[, 2]) / sqrt(nrow(a)))
  }))
  bounds <- if (do_bounds) do.call(rbind, lapply(seq_along(reps), function(r) cbind(rep = r, reps[[r]]$bounds))) else NULL
  diag <- do.call(rbind, lapply(reps, function(r) r$diag))

  list(scenario = scen, oc = oc, calib = calib, accuracy = accuracy,
       bounds = bounds, diag = diag)
}

# Numerical sensitivity of the NPMLE plug-in lfdr to grid resolution, smoothing,
# solver, and starting values. The log-likelihood should agree across optimal fits;
# the plug-in lfdr need not, because the maximiser of the likelihood over G is not
# unique when k is small.
source("R/load.R")
library(parallel)
ncores <- as.integer(Sys.getenv("NCORES", max(1, detectCores() - 1)))
nsets <- as.integer(Sys.getenv("NSETS", "20"))
tau <- main_tau; q <- main_q

cases <- list(list(G = "overlap_hi", k = "5", J = 200),
              list(G = "near_tau", k = "10", J = 200),
              list(G = "separated_hi", k = "het", J = 200))

configs <- list(
  list(name = "barrier M=801 (reference)", M = 801, method = "barrier", smooth = 0),
  list(name = "barrier M=51",  M = 51,  method = "barrier", smooth = 0),
  list(name = "barrier M=201", M = 201, method = "barrier", smooth = 0),
  list(name = "barrier M=201 smooth=0.01", M = 201, method = "barrier", smooth = 0.01),
  list(name = "barrier M=201 smooth=0.1",  M = 201, method = "barrier", smooth = 0.1),
  list(name = "barrier M=201 smooth=1",    M = 201, method = "barrier", smooth = 1),
  list(name = "EM M=201 uniform start",    M = 201, method = "em", init = "uniform"),
  list(name = "EM M=201 random start 1",   M = 201, method = "em", init = 1),
  list(name = "EM M=201 random start 2",   M = 201, method = "em", init = 2)
)

jobs <- expand.grid(case = seq_along(cases), set = seq_len(nsets))
rows <- mclapply(seq_len(nrow(jobs)), function(j) {
  cs <- cases[[jobs$case[j]]]
  G <- G_library[[cs$G]]
  set.seed(500000 + 1000 * jobs$case[j] + jobs$set[j])
  d <- simulate_dataset(cs$J, cs$k, G)
  k0 <- if (cs$k == "het") 10 else as.integer(cs$k)
  lf_or <- oracle_lfdr(d$x, d$k, G, tau)
  ref <- NULL
  do.call(rbind, lapply(configs, function(cf) {
    grid <- make_grid(M = cf$M, taus = tau)
    t0 <- Sys.time()
    fit <- if (cf$method == "em") {
      init <- if (identical(cf$init, "uniform")) NULL else {
        set.seed(cf$init); rexp(length(grid))
      }
      fit_npmle(d$x, d$k, grid = grid, method = "em", init = init)
    } else fit_npmle(d$x, d$k, grid = grid, smooth = cf$smooth)
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    lf <- lfdr_npmle(fit, d$x, d$k, tau)
    if (is.null(ref)) ref <<- lf
    data.frame(case = paste(cs$G, "k =", cs$k, "J =", cs$J), set = jobs$set[j],
               config = cf$name, loglik = fit$loglik, kkt_gap = fit$kkt_gap,
               converged = fit$converged, iter = fit$iter, seconds = secs,
               n_reliable_calls = sum(lf <= q),
               n_reliable_calls_oracle = sum(lf_or <= q),
               lfdr_at_all_success = lfdr_npmle(fit, k0, k0, tau),
               oracle_at_all_success = oracle_lfdr(k0, k0, G, tau),
               max_abs_diff_vs_reference = max(abs(lf - ref)),
               mae_vs_oracle = mean(abs(lf - lf_or)))
  }))
}, mc.cores = ncores)
raw <- do.call(rbind, rows)
raw$loglik_minus_reference <- raw$loglik - ave(ifelse(grepl("reference", raw$config), raw$loglik, NA),
                                                raw$case, raw$set, FUN = function(v) max(v, na.rm = TRUE))
write.csv(raw, "results/tables/grid_sensitivity_raw.csv", row.names = FALSE)

agg <- do.call(rbind, lapply(split(raw, list(raw$case, raw$config), drop = TRUE), function(g) {
  data.frame(case = g$case[1], config = g$config[1], n_datasets = nrow(g),
             converged_frac = mean(g$converged), max_kkt_gap = max(g$kkt_gap),
             mean_loglik_minus_ref = mean(g$loglik_minus_reference),
             mean_seconds = mean(g$seconds),
             mean_calls = mean(g$n_reliable_calls), mean_calls_oracle = mean(g$n_reliable_calls_oracle),
             mean_lfdr_all_success = mean(g$lfdr_at_all_success),
             sd_lfdr_all_success = sd(g$lfdr_at_all_success),
             oracle_lfdr_all_success = g$oracle_at_all_success[1],
             mean_max_abs_diff_vs_ref = mean(g$max_abs_diff_vs_reference),
             mean_mae_vs_oracle = mean(g$mae_vs_oracle))
}))
agg <- agg[order(agg$case, match(agg$config, sapply(configs, `[[`, "name"))), ]
write.csv(agg, "results/tables/grid_sensitivity.csv", row.names = FALSE)
print(agg, digits = 3, row.names = FALSE)

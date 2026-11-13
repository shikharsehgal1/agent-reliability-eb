# Main Monte Carlo study: primary factorial design (independent model) and the
# separate dependence stress cases. Each scenario is cached in results/raw/ so an
# interrupted run can resume; delete results/raw/ to recompute from scratch.
#
# Environment variables: NREP (replicates per scenario, default 200),
# NCORES (parallel workers, default detectCores() - 1).
source("R/load.R")
library(parallel)
nrep <- as.integer(Sys.getenv("NREP", "200"))
ncores <- as.integer(Sys.getenv("NCORES", max(1, detectCores() - 1)))
dir.create("results/raw", recursive = TRUE, showWarnings = FALSE)

design <- rbind(primary_design(), stress_design())
grid <- make_grid()
cat(sprintf("%d scenarios x %d replicates on %d cores\n", nrow(design), nrep, ncores))

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(design)), function(i) {
  scen <- as.list(design[i, ])
  path <- sprintf("results/raw/scen_%04d_nrep%d.rds", scen$id, nrep)
  if (file.exists(path)) return(readRDS(path))
  out <- tryCatch(run_scenario(scen, nrep, grid = grid), error = function(e) {
    list(scenario = scen, error = conditionMessage(e))
  })
  if (is.null(out$error)) saveRDS(out, path)
  out
}, mc.cores = ncores, mc.preschedule = FALSE)

failed <- vapply(res, function(r) !is.null(r$error), logical(1))
if (any(failed)) {
  cat("FAILED scenarios:\n")
  for (r in res[failed]) cat(" ", r$scenario$id, r$scenario$G, r$scenario$k, r$scenario$J, ":", r$error, "\n")
}
saveRDS(list(design = design, results = res, nrep = nrep,
             session = sessionInfo(), elapsed = Sys.time() - t0),
        "results/sim_results.rds")
cat("Elapsed:", format(Sys.time() - t0), "\n")

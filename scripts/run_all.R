# Reproduce every table and figure, in order. From the project root:
#   Rscript scripts/install_deps.R   # once
#   sh scripts/fetch_tau_bench.sh    # once (optional; public data download)
#   Rscript scripts/run_all.R
# Set NREP / NCORES to change the Monte Carlo size or parallelism (see 02_simulation.R).
steps <- c("tests/run_tests.R", "scripts/01_examples.R", "scripts/02_simulation.R",
           "scripts/03_grid_sensitivity.R", "scripts/04_stability.R")
if (length(list.files("data/raw/tau_bench", pattern = "json$"))) steps <- c(steps, "scripts/06_tau_bench.R")
steps <- c(steps, "scripts/05_figures_tables.R")
for (s in steps) {
  cat("\n==========", s, "==========\n")
  t0 <- Sys.time()
  status <- system2("Rscript", s)
  cat(sprintf("---------- %s finished in %.1f min (exit %d)\n", s,
              as.numeric(Sys.time() - t0, units = "mins"), status))
  if (status != 0) stop("Step failed: ", s)
}

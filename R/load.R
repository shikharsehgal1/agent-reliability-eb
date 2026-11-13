# Load project dependencies and source all function files.
# Run scripts from the project root, e.g. `Rscript scripts/run_all.R`.

if (dir.exists(".Rlib")) .libPaths(c(normalizePath(".Rlib"), .libPaths()))

for (pkg in c("lpSolve")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is missing. Run `Rscript scripts/install_deps.R` first.")
  }
}

for (f in c("mixing.R", "simulate.R", "fit_npmle.R", "fit_betabinom.R",
            "decisions.R", "bounds.R", "metrics.R", "scenarios.R", "replicate.R")) {
  source(file.path("R", f))
}

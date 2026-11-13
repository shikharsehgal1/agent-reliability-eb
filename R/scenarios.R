# Mixing distributions and the compact factorial design.

G_library <- list(
  beta_single  = beta_mix(1, 4, 1, "Single Beta(4,1)"),
  separated_lo = beta_mix(c(0.9, 0.1), c(4, 60), c(6, 1.5), "Separated, 10% reliable comp."),
  separated_hi = beta_mix(c(0.6, 0.4), c(4, 60), c(6, 1.5), "Separated, 40% reliable comp."),
  overlap_lo   = beta_mix(c(0.9, 0.1), c(16, 45), c(4, 3), "Overlapping, 10% reliable comp."),
  overlap_hi   = beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3), "Overlapping, 40% reliable comp."),
  near_tau     = beta_mix(c(0.3, 0.7), c(2, 180), c(2, 20), "Mass concentrated near tau"),
  edges        = beta_mix(c(0.4, 0.6), c(0.4, 12), c(6, 0.35), "Mass near 0 and near 1")
)

primary_design <- function() {
  d <- expand.grid(G = names(G_library), k = c("3", "5", "10", "30", "het"),
                   J = c(50, 200, 1000), dependence = "none",
                   stringsAsFactors = FALSE)
  d$id <- seq_len(nrow(d))
  d
}

# Dependence stress cases: separate from the primary (independent) model.
stress_design <- function() {
  d <- expand.grid(G = c("overlap_hi", "separated_hi"), k = c("5", "10"), J = 200,
                   dependence = c("none", "within_corr", "clustered"),
                   stringsAsFactors = FALSE)
  d$id <- 1000 + seq_len(nrow(d))
  d
}

analysis_taus <- c(0.8, 0.9, 0.95)
analysis_qs <- c(0.05, 0.10, 0.20)
main_tau <- 0.9
main_q <- 0.10

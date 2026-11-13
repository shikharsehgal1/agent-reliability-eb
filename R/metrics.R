# Per-replicate operating characteristics of a set of reliability calls.
# FDP is V / max(R, 1); averaging it over replicates estimates the frequentist FDR.

near_width <- 0.05  # "near the boundary" means theta in (tau - near_width, tau]

decision_metrics <- function(reject, theta, tau) {
  null <- theta <= tau
  near <- null & theta > tau - near_width
  R <- sum(reject); V <- sum(reject & null); n_alt <- sum(!null)
  c(R = R, V = V, FDP = if (R > 0) V / R else 0,
    power = if (n_alt > 0) (R - V) / n_alt else NA_real_,
    any_false = as.numeric(V > 0),
    near_false = sum(reject & near), n_near = sum(near))
}

# Binned calibration counts of estimated lfdr against the realised null indicator.
calib_breaks <- seq(0, 1, by = 0.1)

calibration_counts <- function(lfdr_hat, theta, tau) {
  bin <- cut(lfdr_hat, calib_breaks, include.lowest = TRUE, labels = FALSE)
  nb <- length(calib_breaks) - 1
  cbind(n = tabulate(bin, nb),
        sum_lfdr = vapply(seq_len(nb), function(b) sum(lfdr_hat[bin == b]), 0),
        sum_null = vapply(seq_len(nb), function(b) sum(theta[bin == b] <= tau), 0))
}

# Mean and Monte Carlo standard error of a replicate-level metric.
mc_summary <- function(v) {
  v <- v[!is.na(v)]
  c(mean = mean(v), mcse = if (length(v) > 1) sd(v) / sqrt(length(v)) else NA_real_, n = length(v))
}

# Pairwise Jaccard similarity of discovery sets (rows = replicates, cols = tasks).
mean_pairwise_jaccard <- function(rej_mat) {
  n <- nrow(rej_mat)
  inter <- tcrossprod(rej_mat * 1)
  sizes <- rowSums(rej_mat)
  uni <- outer(sizes, sizes, "+") - inter
  jac <- ifelse(uni > 0, inter / uni, 1)
  mean(jac[upper.tri(jac)])
}

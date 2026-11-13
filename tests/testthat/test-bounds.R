grid <- make_grid(M = 401, taus = 0.9)

test_that("k = 1 constructed example matches the analytic bounds", {
  # With one run, the marginal fixes only the mean m1 of G. For m1 > tau,
  # sup lfdr(1) = tau (1 - m1) / ((1 - tau) m1), attained by G on {tau, 1};
  # inf lfdr(1) = 0, approached by mass at 0 and just above tau.
  for (m1 in c(0.92, 0.95, 0.99)) {
    b <- lfdr_bounds(1, 1, marginal_constraints(list("1" = c(1 - m1, m1))), grid, 0.9)
    expect_equal(b$upper, 0.9 * (1 - m1) / (0.1 * m1), tolerance = 1e-6)
    expect_lt(b$lower, 1e-4)
  }
  # For m1 <= tau, a point mass at m1 is feasible, so the sup is 1.
  b <- lfdr_bounds(1, 1, marginal_constraints(list("1" = c(0.2, 0.8))), grid, 0.9)
  expect_equal(b$upper, 1, tolerance = 1e-6)
})

test_that("identified set contains the lfdr of every G with that marginal", {
  G1 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))
  f <- oracle_marginal(3, G1)
  b <- lfdr_bounds(3, 3, marginal_constraints(list("3" = f)), grid, 0.9, return_extremes = TRUE)
  expect_true(b$lower <= oracle_lfdr(3, 3, G1, 0.9) && oracle_lfdr(3, 3, G1, 0.9) <= b$upper)
  # The extremal G's reproduce the same marginal but give the bound values.
  for (w in list(b$w_lower, b$w_upper)) {
    fw <- drop(kernel_matrix(0:3, rep(3, 4), grid) %*% w)
    expect_equal(fw, f, tolerance = 1e-7)
  }
  lf_ext <- sapply(list(b$w_lower, b$w_upper), function(w)
    sum(w * dbinom(3, 3, grid) * (grid <= 0.9)) / sum(w * dbinom(3, 3, grid)))
  expect_equal(lf_ext, c(b$lower, b$upper), tolerance = 1e-6)
  expect_gt(b$upper - b$lower, 0.1)  # genuinely not identified
})

test_that("bounds narrow as more moments are identified", {
  G1 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))
  width <- sapply(c(3, 10, 30), function(k) {
    b <- lfdr_bounds(k, k, marginal_constraints(setNames(list(oracle_marginal(k, G1)), k)), grid, 0.9)
    b$upper - b$lower
  })
  expect_true(all(diff(width) < 0))
})

test_that("confidence-box interval contains the identified set, and infeasibility is reported", {
  set.seed(1)
  G1 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))
  d <- simulate_dataset(500, "5", G1)
  fit <- fit_npmle(d$x, d$k, grid = grid)
  A <- lfdr_bounds(5, 5, marginal_constraints(list("5" = npmle_marginal(fit, 5))), grid, 0.9)
  B <- lfdr_bounds(5, 5, cp_box_constraints(d$x, d$k, 0.05), grid, 0.9)
  expect_true(B$lower <= A$lower + 1e-8 && A$upper <= B$upper + 1e-8)
  # A pmf that is not a binomial mixture (X = 0 or 2 only, never 1, would need
  # theta in {0, 1} -- feasible; but Pr(X = 1) = 1 with k = 2 forces 2p(1-p) = 1: impossible).
  bad <- lfdr_bounds(2, 2, marginal_constraints(list("2" = c(0, 1, 0))), grid, 0.9)
  expect_true(all(bad$status == 2))
})

test_that("Clopper-Pearson box covers the true marginal at roughly the nominal rate or better", {
  set.seed(2)
  G1 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))
  f <- oracle_marginal(5, G1)
  cover <- replicate(300, {
    d <- simulate_dataset(200, "5", G1)
    box <- cp_box_constraints(d$x, d$k, 0.10)[[1]]
    all(box$lower <= f & f <= box$upper)
  })
  expect_gte(mean(cover), 0.90 - 3 * sqrt(0.09 / 300))
})

test_that("early-stopped robust upper bounds equal the full computation where it matters", {
  set.seed(3)
  d <- simulate_dataset(300, "het", G_library$separated_hi)
  con <- cp_box_constraints(d$x, d$k, 0.05)
  full <- robust_lfdr_upper(d$x, d$k, con, grid, 0.9, stop_above = 1)
  early <- robust_lfdr_upper(d$x, d$k, con, grid, 0.9, stop_above = 0.05)
  keep <- early <= 0.05
  expect_equal(early[keep], full[keep])
  expect_true(all(full[!keep] > 0.05 - 1e-9))
  # The upper bound is nonincreasing in x within each run budget.
  for (kk in unique(d$k)) {
    u <- full[d$k == kk][order(d$x[d$k == kk])]
    expect_true(all(diff(u) <= 1e-7))
  }
})

test_that("box-robust rule calls nothing the box cannot certify, and is a subset of oracle calls when the box covers", {
  set.seed(4)
  G <- G_library$separated_hi
  d <- simulate_dataset(1000, "30", G)
  g <- make_grid(M = 301, taus = 0.9)
  rej <- reject_robust(d$x, d$k, 0.1, 0.9, g)
  expect_gt(sum(rej), 0)
  box <- cp_box_constraints(d$x, d$k, 0.05)[[1]]
  f <- oracle_marginal(30, G)
  if (all(box$lower <= f & f <= box$upper)) {
    expect_true(all(oracle_lfdr(d$x[rej], d$k[rej], G, 0.9) <= 0.05 + 1e-6))
  }
})

test_that("box-robust rule keeps FDR <= q in the hard near-tau, heterogeneous-k case", {
  set.seed(5)
  g <- make_grid(M = 201, taus = 0.9)
  fdp <- replicate(150, {
    d <- simulate_dataset(200, "het", G_library$near_tau)
    c(robust = decision_metrics(reject_robust(d$x, d$k, 0.1, 0.9, g), d$theta, 0.9)[["FDP"]],
      npmle = decision_metrics(lfdr_npmle(fit_npmle(d$x, d$k, grid = g), d$x, d$k, 0.9) <= 0.1,
                               d$theta, 0.9)[["FDP"]])
  })
  se <- apply(fdp, 1, sd) / sqrt(ncol(fdp))
  expect_lt(mean(fdp["robust", ]), 0.1 + 3 * se[["robust"]])
})

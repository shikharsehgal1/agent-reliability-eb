G2 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))

test_that("closed-form oracle lfdr matches numerical integration", {
  for (k in c(3, 10)) for (x in c(0, k - 1, k)) {
    num <- integrate(function(p) dbinom(x, k, p) * d_mix(p, G2), 0, 0.9)$value
    den <- integrate(function(p) dbinom(x, k, p) * d_mix(p, G2), 0, 1)$value
    expect_equal(oracle_lfdr(x, k, G2, 0.9), num / den, tolerance = 1e-6)
  }
})

test_that("oracle marginal is a probability vector", {
  for (G in G_library) for (k in c(1, 3, 10, 20)) {
    f <- oracle_marginal(k, G)
    expect_true(all(f >= 0))
    expect_equal(sum(f), 1, tolerance = 1e-10)
  }
})

test_that("lfdr lies in [0, 1] and decreases in x for fixed k", {
  for (G in G_library) {
    lf <- oracle_lfdr(0:10, 10, G, 0.9)
    expect_true(all(lf >= 0 & lf <= 1))
    expect_true(all(diff(lf) <= 1e-12))
  }
})

test_that("Beta-prior posteriors agree with the oracle formula for a Beta G", {
  fit <- list(a = 2.5, b = 0.7)
  expect_equal(lfdr_betabinom(fit, 0:5, 5, 0.9),
               oracle_lfdr(0:5, 5, beta_mix(1, 2.5, 0.7), 0.9), tolerance = 1e-10)
  expect_equal(lfdr_unpooled(3, 3, 0.9), 0.9^4)  # Beta(4, 1) posterior: Pr(theta <= .9) = .9^4
})

G2 <- beta_mix(c(0.6, 0.4), c(16, 45), c(4, 3))

test_that("NPMLE weights are on the simplex and the fit is certified", {
  set.seed(1)
  d <- simulate_dataset(300, "het", G2)
  for (method in c("barrier", "em")) {
    fit <- fit_npmle(d$x, d$k, method = method)
    expect_true(all(fit$w >= 0))
    expect_equal(sum(fit$w), 1, tolerance = 1e-12)
    # The certified gap bounds the log-likelihood suboptimality by J * gap.
    # EM is slow and may stop at its iteration cap; the barrier must converge.
    if (method == "barrier") expect_true(fit$converged)
    expect_lt(fit$kkt_gap, 1e-5)
    for (kk in c(3, 20)) expect_equal(sum(npmle_marginal(fit, kk)), 1, tolerance = 1e-10)
  }
})

test_that("barrier and EM reach the same maximised log-likelihood", {
  set.seed(2)
  d <- simulate_dataset(500, "5", G2)
  fb <- fit_npmle(d$x, d$k, method = "barrier")
  fe <- fit_npmle(d$x, d$k, method = "em")
  expect_equal(fb$loglik, fe$loglik, tolerance = 1e-4 * abs(fb$loglik))
  # ... and the same fitted marginal, which is the identified object.
  expect_equal(npmle_marginal(fb, 5), npmle_marginal(fe, 5), tolerance = 1e-3)
})

test_that("with fixed k and an interior empirical pmf, the fitted marginal is the empirical pmf", {
  set.seed(3)
  d <- simulate_dataset(5000, "3", G2)
  fit <- fit_npmle(d$x, d$k)
  emp <- tabulate(d$x + 1, 4) / nrow(d)
  expect_equal(npmle_marginal(fit, 3), emp, tolerance = 1e-4)
})

test_that("known-G recovery: the marginal is recovered; lfdr only when G avoids tau", {
  set.seed(4)
  d <- simulate_dataset(20000, "10", G2)
  fit <- fit_npmle(d$x, d$k)
  tv <- 0.5 * sum(abs(npmle_marginal(fit, 10) - oracle_marginal(10, G2)))
  expect_lt(tv, 0.02)
  # When G has little mass near tau and k is moderate, the plug-in lfdr tracks the
  # oracle. (When G has substantial mass near tau it does not, even for k = 60 and
  # J = 5000: task-averaged errors of 0.05-0.09 were observed. See progress_report.md.)
  Gs <- G_library$separated_hi
  d <- simulate_dataset(5000, "20", Gs)
  fit <- fit_npmle(d$x, d$k)
  expect_lt(mean(abs(lfdr_npmle(fit, d$x, d$k, 0.9) - oracle_lfdr(d$x, d$k, Gs, 0.9))), 0.03)
})

test_that("Beta-binomial fit recovers Beta parameters and flags the no-overdispersion boundary", {
  set.seed(5)
  d <- simulate_dataset(5000, "10", beta_mix(1, 6, 2))
  fit <- fit_betabinom(d$x, d$k)
  expect_equal(fit$mu, 0.75, tolerance = 0.02)
  expect_equal(fit$phi, 8, tolerance = 0.2 * 8)
  expect_false(fit$boundary)
  # Under-dispersed counts (variance below binomial) put the MLE of phi at +infinity.
  # (With a true point mass G this happens in roughly half of samples, not always.)
  x <- rep(c(6, 7, 8), length.out = 2000)
  fit0 <- fit_betabinom(x, rep(10, 2000))
  expect_true(fit0$flags[["phi_upper"]])
})

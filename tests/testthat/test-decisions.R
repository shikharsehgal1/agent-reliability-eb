test_that("exact binomial p-value equals Pr(Bin(k, tau) >= x)", {
  expect_equal(pval_exact_binom(3, 3, 0.9), 0.9^3)
  expect_equal(pval_exact_binom(0, 5, 0.9), 1)
  expect_equal(pval_exact_binom(4, 5, 0.9), 0.9^5 + 5 * 0.9^4 * 0.1)
})

test_that("BH and BY match p.adjust", {
  set.seed(1)
  p <- c(runif(80), rbeta(20, 0.1, 5))
  for (q in c(0.05, 0.1)) {
    expect_identical(reject_bh(p, q), p.adjust(p, "BH") <= q)
    expect_identical(reject_by(p, q), p.adjust(p, "BY") <= q)
  }
})

test_that("lfdr step-up rejects the largest prefix with mean lfdr <= q", {
  lf <- c(0.30, 0.01, 0.05, 0.20, 0.90)
  # sorted: .01 .05 .20 .30 .90; running means .01 .03 .0867 .14 -> 3 rejections at q = .1
  expect_identical(reject_lfdr_stepup(lf, 0.1), c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(reject_lfdr_threshold(lf, 0.1), c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("BH controls FDR in a small independent simulation", {
  set.seed(2)
  fdp <- replicate(400, {
    d <- simulate_dataset(200, "10", G_library$overlap_hi)
    m <- decision_metrics(reject_bh(pval_exact_binom(d$x, d$k, 0.9), 0.1), d$theta, 0.9)
    m["FDP"]
  })
  expect_lt(mean(fdp), 0.1 + 3 * sd(fdp) / sqrt(length(fdp)))
})

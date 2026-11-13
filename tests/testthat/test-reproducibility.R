test_that("simulation is reproducible from a seed", {
  set.seed(11); a <- simulate_dataset(100, "het", G_library$edges)
  set.seed(11); b <- simulate_dataset(100, "het", G_library$edges)
  expect_identical(a, b)
})

test_that("a replicate is reproducible from its scenario seed", {
  scen <- list(G = "near_tau", k = "5", J = 80, dependence = "none", id = 7)
  r1 <- run_replicate(scen, scenario_seed(7, 3))
  r2 <- run_replicate(scen, scenario_seed(7, 3))
  expect_identical(r1$metrics, r2$metrics)
  expect_identical(r1$bounds, r2$bounds)
  r3 <- run_replicate(scen, scenario_seed(7, 4))
  expect_false(identical(r1$data$x, r3$data$x))
})

test_that("stress generators keep theta as the per-run success probability", {
  set.seed(12)
  d <- simulate_dataset(20000, "10", G_library$overlap_hi, dependence = "within_corr")
  expect_equal(mean(d$x / d$k), mean(d$theta), tolerance = 0.01)
  # ... but X is overdispersed relative to Binomial(k, theta).
  expect_gt(var(d$x - d$k * d$theta), 1.5 * mean(d$k * d$theta * (1 - d$theta)))
  d2 <- simulate_dataset(100, "5", G_library$overlap_hi, dependence = "clustered")
  expect_equal(length(unique(d2$theta)), 10)
})

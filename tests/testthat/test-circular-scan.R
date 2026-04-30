test_that("circular_scan returns correct class and structure", {
  set.seed(1)
  n <- 20
  res <- circular_scan(
    cases       = rpois(n, lambda = 10),
    population  = rep(1000, n),
    region_id   = 1:n,
    x           = runif(n, 0, 10),
    y           = runif(n, 0, 10),
    nsim = 49, seed = 42
  )
  expect_s3_class(res, "circular_scan")
  expect_true(all(c("most_likely_cluster", "pvalue") %in% names(res)))
})

test_that("circular_scan detects obvious cluster", {
  set.seed(99)
  n <- 25
  res <- circular_scan(
    cases       = c(rpois(5, 50), rpois(n - 5, 5)),
    population  = rep(1000, n),
    region_id   = 1:n,
    x           = c(rep(0, 5), runif(n - 5, 5, 15)),
    y           = c(rep(0, 5), runif(n - 5, 5, 15)),
    nsim = 199, seed = 99
  )
  expect_lt(res$pvalue, 0.10)
  expect_true(all(1:3 %in% res$most_likely_cluster$region_ids))
})

test_that("circular_scan errors when cases is not aggregated", {
  # 4 rows, 2 unique regions -> not aggregated
  expect_error(
    circular_scan(
      cases      = c(5, 3, 4, 2),
      population = c(1000, 1000, 2000, 2000),
      region_id  = c(1, 1, 2, 2),
      x          = c(1, 1, 2, 2),
      y          = c(1, 1, 2, 2)
    ),
    "must already be aggregated"
  )
})

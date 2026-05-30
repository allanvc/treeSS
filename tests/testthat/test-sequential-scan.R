test_that("sequential_scan returns valid object structure (tree-spatial)", {
  set.seed(1L)
  ex <- generate_example_data(n_regions = 20L, seed = 11L)
  res <- sequential_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 3L, nsim = 29L, seed = 7L,
    verbose = FALSE
  )

  expect_s3_class(res, "sequential_scan")
  expect_true(all(c("clusters", "iterations", "regions", "tree",
                    "alpha", "nsim", "buffer_size", "n_iter",
                    "scan_type") %in% names(res)))
  expect_equal(res$scan_type, "treespatial")
  expect_true(is.data.frame(res$clusters))
  expect_true(all(c("iteration", "node_id", "n_regions", "cases",
                    "expected", "rr", "llr", "pvalue", "significant")
                  %in% names(res$clusters)))
  expect_gte(res$n_iter, 1L)
  expect_lte(res$n_iter, 3L)
  expect_equal(nrow(res$clusters), res$n_iter)
})

test_that("sequential_scan respects max_iter", {
  set.seed(2L)
  ex <- generate_example_data(n_regions = 18L, seed = 13L)
  res <- sequential_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 2L, nsim = 29L, seed = 1L, verbose = FALSE
  )
  expect_lte(nrow(res$clusters), 2L)
})

test_that("sequential_scan stops at first non-significant MLC", {
  set.seed(3L)
  # No injected clusters: scan should detect at most one (likely non-sig)
  n <- 15
  cases <- rpois(n, lambda = 5)
  pop   <- rep(1000, n)
  res <- sequential_scan(
    cases = cases, population = pop,
    region_id = seq_len(n),
    x = runif(n, 0, 10), y = runif(n, 0, 10),
    max_iter = 5L, nsim = 99L, alpha = 0.05,
    seed = 7L, verbose = FALSE
  )
  expect_s3_class(res, "sequential_scan")
  expect_equal(res$scan_type, "circular")
  # First non-significant MLC stops the loop, so the LAST recorded
  # cluster is the first non-sig one (if any), preceded by sig ones.
  if (nrow(res$clusters) > 1L) {
    # All but the last must be significant
    expect_true(all(res$clusters$significant[-nrow(res$clusters)]))
  }
})

test_that("sequential_scan with buffer_size > 0 works", {
  set.seed(4L)
  ex <- generate_example_data(n_regions = 25L, seed = 21L)
  res0 <- sequential_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 2L, nsim = 29L, buffer_size = 0L,
    seed = 11L, verbose = FALSE
  )
  res2 <- sequential_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 2L, nsim = 29L, buffer_size = 2L,
    seed = 11L, verbose = FALSE
  )
  expect_equal(res0$buffer_size, 0L)
  expect_equal(res2$buffer_size, 2L)
  # Both produce valid sequential_scan objects
  expect_s3_class(res0, "sequential_scan")
  expect_s3_class(res2, "sequential_scan")
})

test_that("sequential_scan circular mode works", {
  set.seed(5L)
  n <- 20
  cases <- rpois(n, lambda = 8)
  cases[1:4] <- rpois(4, lambda = 40)
  res <- sequential_scan(
    cases = cases, population = rep(500, n),
    region_id = seq_len(n),
    x = runif(n, 0, 10), y = runif(n, 0, 10),
    max_iter = 3L, nsim = 49L, seed = 2L, verbose = FALSE
  )
  expect_s3_class(res, "sequential_scan")
  expect_equal(res$scan_type, "circular")
})

test_that("sequential_scan tree-only mode works", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3)
  )
  cases <- c(60, 5, 4, 3)
  pop   <- c(100, 100, 100, 100)
  res <- sequential_scan(
    cases = cases, population = pop, tree = tree,
    max_iter = 2L, nsim = 49L, seed = 9L, verbose = FALSE
  )
  expect_s3_class(res, "sequential_scan")
  expect_equal(res$scan_type, "tree")
  expect_null(res$regions)
})

test_that("sequential_scan print runs without error", {
  set.seed(6L)
  ex <- generate_example_data(n_regions = 14L, seed = 31L)
  res <- sequential_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 2L, nsim = 19L, seed = 1L, verbose = FALSE
  )
  expect_output(print(res), "Sequential Scan")
})

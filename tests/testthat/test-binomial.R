test_that("circular_scan binomial runs", {
  set.seed(123)
  n <- 15
  cases <- rpois(n, lambda = 10)
  cases[1:4] <- rpois(4, lambda = 50)
  res <- circular_scan(
    cases = cases, population = rep(1000, n),
    region_id = 1:n, x = runif(n, 0, 10), y = runif(n, 0, 10),
    nsim = 49, seed = 42, model = "binomial"
  )
  expect_s3_class(res, "circular_scan")
  expect_equal(res$model, "binomial")
})

test_that("tree_scan binomial requires population", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7, 8),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3, 3)
  )
  cases <- c(50, 5, 3, 2, 4)
  pop   <- c(100, 100, 100, 100, 100)

  expect_error(
    tree_scan(tree, cases, population = NULL, nsim = 9, seed = 1,
              model = "binomial"),
    "binomial model requires 'population'"
  )

  res <- tree_scan(tree, cases, population = pop, nsim = 49, seed = 42,
                   model = "binomial")
  expect_s3_class(res, "tree_scan")
})

test_that("treespatial_scan binomial runs", {
  ex <- generate_example_data(n_regions = 16, seed = 77)
  res <- treespatial_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    nsim = 49, seed = 1, model = "binomial"
  )
  expect_s3_class(res, "treespatial_scan")
  expect_equal(res$model, "binomial")
})

test_that("n_cores > 1 path runs (treespatial)", {
  skip_on_cran()
  ex <- generate_example_data(n_regions = 16, seed = 11)
  args <- list(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    nsim = 49, seed = 1
  )
  r1 <- do.call(treespatial_scan, c(args, list(n_cores = 1L)))
  r2 <- do.call(treespatial_scan, c(args, list(n_cores = 2L)))
  expect_equal(r1$most_likely_cluster$llr, r2$most_likely_cluster$llr)
  expect_setequal(r1$most_likely_cluster$region_ids,
                  r2$most_likely_cluster$region_ids)
})

test_that("iterative_scan accepts model = 'binomial' and n_cores", {
  skip_on_cran()
  ex <- generate_example_data(n_regions = 12, seed = 33)
  res <- iterative_scan(
    cases = ex$cases, population = ex$population,
    region_id = ex$region_id, x = ex$x, y = ex$y,
    node_id = ex$node_id, tree = ex$tree,
    max_iter = 2, nsim = 29, seed = 1,
    model = "binomial", n_cores = 2L, verbose = FALSE
  )
  expect_s3_class(res, "iterative_scan")
})

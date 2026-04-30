test_that("tree_scan returns correct class and structure", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3)
  )
  cases <- c(50, 5, 3, 2)
  pop   <- c(100, 100, 100, 100)

  result <- tree_scan(tree, cases, population = pop, nsim = 49, seed = 42)

  expect_s3_class(result, "tree_scan")
  expect_true("most_likely_cluster" %in% names(result))
  expect_true("all_cuts" %in% names(result))
  expect_true("significant_cuts" %in% names(result))
  expect_equal(nrow(result$all_cuts), 7)
})

test_that("tree_scan detects obvious cluster", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3)
  )
  # Leaf 4 has far more cases
  cases <- c(200, 5, 3, 2)
  pop   <- c(100, 100, 100, 100)

  result <- tree_scan(tree, cases, population = pop, nsim = 199, seed = 77)

  # Node 4 or its parent (2) should be the most likely cluster
  expect_true(result$most_likely_cluster$node_id %in% c(2, 4))
  expect_lt(result$pvalue, 0.05)
})

test_that("tree_scan handles default population", {
  tree <- data.frame(
    node_id   = c(1, 2, 3),
    parent_id = c(NA, 1, 1)
  )
  cases <- c(10, 5)

  result <- tree_scan(tree, cases, nsim = 19, seed = 1)
  expect_s3_class(result, "tree_scan")
  expect_equal(result$total_population, 2)  # default pop = 1 per leaf
})

test_that("tree_scan errors on wrong cases length", {
  tree <- data.frame(
    node_id   = c(1, 2, 3),
    parent_id = c(NA, 1, 1)
  )
  expect_error(tree_scan(tree, c(10, 5, 3)), "must equal the number of leaf")
})

test_that("tree_scan print method works", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5),
    parent_id = c(NA, 1, 1, 2, 2)
  )
  result <- tree_scan(tree, c(30, 5, 2), population = rep(100, 3),
                      nsim = 19, seed = 1)
  expect_output(print(result), "Tree-Based Scan")
})

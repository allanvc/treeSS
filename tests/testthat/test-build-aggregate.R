test_that("aggregate_tree sums correctly", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3)
  )
  full <- aggregate_tree(
    cases     = c(10, 5, 3, 8,  2, 7, 4, 1,  6, 3, 9, 2),
    region_id = rep(1:3, each = 4),
    node_id   = rep(c(4, 5, 6, 7), times = 3),
    tree      = tree
  )

  expect_equal(nrow(full), 7)
  expect_equal(ncol(full), 3)

  leaf4 <- which(rownames(full) == "4")
  leaf5 <- which(rownames(full) == "5")
  node2 <- which(rownames(full) == "2")
  expect_equal(full[node2, ], full[leaf4, ] + full[leaf5, ])

  root <- which(rownames(full) == "1")
  per_region <- c(10 + 5 + 3 + 8,
                   2 + 7 + 4 + 1,
                   6 + 3 + 9 + 2)
  expect_equal(unname(full[root, ]), per_region)
})

test_that("aggregate_tree errors on length mismatch", {
  tree <- data.frame(node_id = c(1, 2, 3), parent_id = c(NA, 1, 1))
  expect_error(
    aggregate_tree(cases = 1:5, region_id = 1:5, node_id = rep(2, 3),
                   tree = tree),
    "same length"
  )
})

test_that("aggregate_tree accepts tree_node_id/tree_parent_id form", {
  res1 <- aggregate_tree(
    cases     = c(10, 5, 3, 8),
    region_id = rep(1:2, each = 2),
    node_id   = c(2, 3, 2, 3),
    tree      = data.frame(node_id = c(1, 2, 3),
                            parent_id = c(NA, 1, 1))
  )
  res2 <- aggregate_tree(
    cases     = c(10, 5, 3, 8),
    region_id = rep(1:2, each = 2),
    node_id   = c(2, 3, 2, 3),
    tree_node_id = c(1, 2, 3),
    tree_parent_id = c(NA, 1, 1)
  )
  expect_equal(res1, res2)
})

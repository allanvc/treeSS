test_that("treespatial_scan requires data and validates columns", {
  tree <- data.frame(node_id = c(1, 2, 3), parent_id = c(NA, 1, 1))
  dat <- data.frame(cases = 1:3, population = c(100, 100, 100),
                    region_id = 1:3, x = 1:3, y = 1:3, node_id = c(2, 2, 3))
  # `data` is required
  expect_error(treespatial_scan(), "is required")
  # unknown column name is reported
  expect_error(
    treespatial_scan(dat, cases = nope, population = population,
                     region_id = region_id, x = x, y = y, node_id = node_id,
                     tree = tree),
    "nope"
  )
})

test_that("treespatial_scan returns expected class and structure", {
  ex <- generate_example_data(n_regions = 16, seed = 1)
  dat <- data.frame(cases = ex$cases, population = ex$population,
                    region_id = ex$region_id, x = ex$x, y = ex$y,
                    node_id = ex$node_id)
  res <- treespatial_scan(dat, cases, population, region_id, x, y, node_id,
                          tree = ex$tree, nsim = 29, seed = 42)
  expect_s3_class(res, "treespatial_scan")
  expect_true(all(c("most_likely_cluster", "pvalue", "nsim") %in% names(res)))
  expect_gte(res$most_likely_cluster$llr, 0)
  expect_gte(res$pvalue, 0); expect_lte(res$pvalue, 1)
})

test_that("treespatial_scan detects injected cluster", {
  ex <- generate_example_data(n_regions = 25, cluster_regions = 1:5,
                              cluster_leaves = c(3, 4),
                              rr = 4.0, Cg = 300L, seed = 77)
  dat <- data.frame(cases = ex$cases, population = ex$population,
                    region_id = ex$region_id, x = ex$x, y = ex$y,
                    node_id = ex$node_id)
  res <- treespatial_scan(dat, cases, population, region_id, x, y, node_id,
                          tree = ex$tree, nsim = 199, seed = 77)
  detected <- res$most_likely_cluster$region_ids
  injected <- ex$true_cluster$region_ids
  overlap <- length(intersect(detected, injected)) / length(injected)
  expect_gte(overlap, 0.5)
})

test_that("treespatial_scan accepts a column name held in a variable", {
  ex <- generate_example_data(n_regions = 16, seed = 1)
  dat <- data.frame(cases = ex$cases, population = ex$population,
                    region_id = ex$region_id, x = ex$x, y = ex$y,
                    node_id = ex$node_id)
  pop_col <- "population"          # name held in a variable (programmatic use)
  res <- treespatial_scan(dat, cases, population = pop_col,
                          region_id = region_id, x = x, y = y,
                          node_id = node_id, tree = ex$tree,
                          nsim = 19, seed = 1)
  res2 <- treespatial_scan(dat, cases, population = population,
                           region_id = region_id, x = x, y = y,
                           node_id = node_id, tree = ex$tree,
                           nsim = 19, seed = 1)
  expect_s3_class(res, "treespatial_scan")
  expect_equal(res$most_likely_cluster$llr, res2$most_likely_cluster$llr)
})

test_that("treespatial_scan accepts tree_node_id/tree_parent_id form", {
  ex <- generate_example_data(n_regions = 16, seed = 1)
  dat <- data.frame(cases = ex$cases, population = ex$population,
                    region_id = ex$region_id, x = ex$x, y = ex$y,
                    node_id = ex$node_id)
  res1 <- treespatial_scan(dat, cases, population, region_id, x, y, node_id,
                           tree = ex$tree, nsim = 19, seed = 1)
  res2 <- treespatial_scan(dat, cases, population, region_id, x, y, node_id,
                           tree_node_id = ex$tree$node_id,
                           tree_parent_id = ex$tree$parent_id,
                           nsim = 19, seed = 1)
  expect_equal(res1$most_likely_cluster$llr,
               res2$most_likely_cluster$llr)
})

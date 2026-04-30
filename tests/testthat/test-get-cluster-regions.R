test_that("get_cluster_regions dispatches on iterative_scan", {
  data(london_collisions)
  data(london_tree)

  iter <- iterative_scan(
    cases       = london_collisions$cases,
    population  = london_collisions$population,
    region_id   = london_collisions$region_id,
    x           = london_collisions$x,
    y           = london_collisions$y,
    node_id     = london_collisions$node_id,
    tree        = london_tree,
    max_iter    = 3, nsim = 19, seed = 42,
    max_pop_pct = 0.25, n_cores = 1L, verbose = FALSE
  )

  expect_s3_class(iter, "iterative_scan")
  expect_gt(iter$n_iter, 0)

  # overlap = TRUE: each iteration replicates the full region universe
  cr <- get_cluster_regions(iter, overlap = TRUE)
  expect_s3_class(cr, "data.frame")
  expect_true(all(c("cluster", "node_id", "llr", "pvalue",
                     "pvalue_adjusted", "significant", "panel")
                   %in% names(cr)))
  # Each iteration covers the full set of regions
  per_panel <- table(cr$panel)
  expect_true(all(per_panel == nrow(iter$regions)))
  # Total rows = n_iter * n_regions
  expect_equal(nrow(cr), iter$n_iter * nrow(iter$regions))

  # overlap = FALSE: one row per region
  cr_no <- get_cluster_regions(iter, overlap = FALSE)
  expect_equal(nrow(cr_no), nrow(iter$regions))
  expect_false("panel" %in% names(cr_no))
})


test_that("get_cluster_regions still works on single-pass scans", {
  data(london_collisions)
  data(london_tree)

  res <- treespatial_scan(
    cases       = london_collisions$cases,
    population  = london_collisions$population,
    region_id   = london_collisions$region_id,
    x           = london_collisions$x,
    y           = london_collisions$y,
    node_id     = london_collisions$node_id,
    tree        = london_tree,
    max_pop_pct = 0.25, nsim = 19, seed = 42, n_cores = 1L
  )

  cr1 <- get_cluster_regions(res, n_clusters = 1, overlap = FALSE)
  expect_equal(nrow(cr1), nrow(res$regions))
  expect_true("cluster" %in% names(cr1))
})


test_that("get_cluster_regions on iterative_scan rejects tree-only", {
  data(london_collisions)
  data(london_tree)

  # Tree-only iterative: aggregate cases per leaf to make a tree-only call
  agg <- stats::aggregate(cases ~ node_id, data = london_collisions, FUN = sum)
  pop <- rep(sum(london_collisions$population) / nrow(agg), nrow(agg))
  iter_t <- iterative_scan(
    cases = agg$cases, population = pop, tree = london_tree,
    tree_node_id = NULL, tree_parent_id = NULL,
    node_id = NULL, region_id = NULL, x = NULL, y = NULL,
    max_iter = 2, nsim = 19, seed = 42, verbose = FALSE
  )

  # Tree-only iter has no $regions; mapping should error politely
  expect_error(get_cluster_regions(iter_t),
                regexp = "no \\$regions table")
})

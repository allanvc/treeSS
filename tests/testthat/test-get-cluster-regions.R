test_that("get_cluster_regions dispatches on sequential_scan", {
  data(london_collisions)
  data(london_tree)

  seq_res <- sequential_scan(
    london_collisions,
    cases = cases, population = population, region_id = region_id,
    x = x, y = y, node_id = node_id, tree = london_tree,
    max_iter = 3, nsim = 19, seed = 42,
    max_pop_pct = 0.25, n_cores = 1L, verbose = FALSE
  )

  expect_s3_class(seq_res, "sequential_scan")
  expect_gt(seq_res$n_iter, 0)

  cr <- get_cluster_regions(seq_res, overlap = TRUE)
  expect_s3_class(cr, "data.frame")
  expect_true(all(c("cluster", "node_id", "llr", "pvalue", "panel")
                   %in% names(cr)))
  per_panel <- table(cr$panel)
  expect_true(all(per_panel == nrow(seq_res$regions)))
  expect_equal(nrow(cr), seq_res$n_iter * nrow(seq_res$regions))

  cr_no <- get_cluster_regions(seq_res, overlap = FALSE)
  expect_equal(nrow(cr_no), nrow(seq_res$regions))
  expect_false("panel" %in% names(cr_no))
})


test_that("get_cluster_regions still works on single-pass scans", {
  data(london_collisions)
  data(london_tree)

  res <- treespatial_scan(
    london_collisions,
    cases = cases, population = population, region_id = region_id,
    x = x, y = y, node_id = node_id, tree = london_tree,
    max_pop_pct = 0.25, nsim = 19, seed = 42, n_cores = 1L
  )

  cr1 <- get_cluster_regions(res, n_clusters = 1, overlap = FALSE)
  expect_equal(nrow(cr1), nrow(res$regions))
  expect_true("cluster" %in% names(cr1))
})


test_that("get_cluster_regions on sequential_scan rejects tree-only", {
  data(london_collisions)
  data(london_tree)

  agg <- stats::aggregate(cases ~ node_id, data = london_collisions,
                          FUN = sum)
  agg$population <- rep(sum(london_collisions$population) / nrow(agg),
                        nrow(agg))
  seq_t <- sequential_scan(
    agg, cases = cases, population = population, node_id = node_id,
    tree = london_tree, max_iter = 2, nsim = 19, seed = 42, verbose = FALSE
  )

  expect_error(get_cluster_regions(seq_t),
                regexp = "no \\$regions table")
})

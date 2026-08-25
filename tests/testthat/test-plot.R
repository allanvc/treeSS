test_that("plot.treess draws centroids in base graphics and returns regions", {
  data(london_collisions); data(london_tree)
  res <- treespatial_scan(
    london_collisions,
    cases = cases, population = population, region_id = region_id,
    x = x, y = y, node_id = node_id, tree = london_tree,
    nsim = 19, seed = 42
  )
  pdf(NULL); on.exit(dev.off(), add = TRUE)

  cr1 <- plot(res)
  expect_s3_class(cr1, "data.frame")
  expect_equal(length(unique(cr1$panel)), 1L)
  expect_true(all(c("region_id", "cluster", "panel") %in% names(cr1)))

  cr2 <- plot(res, n_clusters = 2L)
  expect_true(length(unique(cr2$panel)) <= 2L)
})

test_that("plot.treess works for circular and sequential scans", {
  ex <- generate_example_data(n_regions = 25, seed = 3)
  dat <- data.frame(cases = ex$cases, population = ex$population,
                    region_id = ex$region_id, x = ex$x, y = ex$y,
                    node_id = ex$node_id)
  pdf(NULL); on.exit(dev.off(), add = TRUE)

  circ <- circular_scan(dat[!duplicated(dat$region_id), ],
                        cases = cases, population = population,
                        region_id = region_id, x = x, y = y,
                        nsim = 19, seed = 1)
  expect_s3_class(plot(circ), "data.frame")

  seq_res <- sequential_scan(dat, cases = cases, population = population,
                             region_id = region_id, x = x, y = y,
                             node_id = node_id, tree = ex$tree,
                             max_iter = 2L, nsim = 19, seed = 1,
                             verbose = FALSE)
  out <- plot(seq_res)
  expect_s3_class(out, "data.frame")
})

test_that("plot.treess with an sf map returns a ggplot and joins via key", {
  skip_if_not_installed("sf")
  skip_if_not_installed("ggplot2")
  data(london_collisions); data(london_tree); data(london_boroughs_map)
  res <- treespatial_scan(
    london_collisions,
    cases = cases, population = population, region_id = region_id,
    x = x, y = y, node_id = node_id, tree = london_tree,
    nsim = 19, seed = 42
  )
  key <- unique(london_collisions[, c("region_id", "borough")])
  p <- plot(res, map = london_boroughs_map, map_id = "NAME", key = key)
  expect_s3_class(p, "ggplot")

  # Members of the MLC must be marked exactly once in the joined layer
  built <- ggplot2::ggplot_build(p)
  mlc <- most_likely_cluster(res)
  expect_equal(sum(!is.na(p$data$cluster)), length(mlc$region_ids))

  expect_error(plot(res, map = london_boroughs_map),
               "map_id")
  expect_error(plot(res, map = data.frame(a = 1), map_id = "a"),
               "sf")
  # Without a key the integer region_id does not match borough names
  expect_warning(plot(res, map = london_boroughs_map, map_id = "NAME"),
                 "no polygon")
})

test_that("plot.tree_scan draws a dot chart of top cuts", {
  tree <- data.frame(node_id = c(1, 2, 3, 4, 5, 6, 7, 8),
                     parent_id = c(NA, 1, 1, 2, 2, 3, 3, 3))
  leaf_data <- data.frame(node_id = c(4, 5, 6, 7, 8),
                          cases = c(50, 5, 3, 2, 4))
  res <- tree_scan(leaf_data, cases = cases, node_id = node_id,
                   tree = tree, nsim = 19, seed = 1)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  cuts <- plot(res, n = 5)
  expect_equal(nrow(cuts), 5L)
  expect_true(all(diff(cuts$llr) <= 0))
})

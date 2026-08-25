# Tests for the treess parent class and accessor generics
# (most_likely_cluster, secondary_clusters, pvalue) added in 0.2.7.

make_dat <- function() {
  ex <- generate_example_data(n_regions = 16, seed = 1)
  list(dat = data.frame(cases = ex$cases, population = ex$population,
                        region_id = ex$region_id, x = ex$x, y = ex$y,
                        node_id = ex$node_id),
       tree = ex$tree)
}

make_ts <- function() {
  ex <- make_dat()
  treespatial_scan(ex$dat, cases = cases, population = population,
                   region_id = region_id, x = x, y = y,
                   node_id = node_id, tree = ex$tree,
                   nsim = 19, seed = 42)
}

test_that("all scan results inherit from treess", {
  ex <- make_dat()

  res_ts <- make_ts()
  expect_s3_class(res_ts, "treespatial_scan")
  expect_s3_class(res_ts, "treess")

  agg <- aggregate(cbind(cases = ex$dat$cases) ~ region_id,
                   data = ex$dat, FUN = sum)
  first <- ex$dat[!duplicated(ex$dat$region_id), ]
  dat_reg <- data.frame(region_id = agg$region_id, cases = agg$cases,
                        population = first$population[match(agg$region_id,
                                                            first$region_id)],
                        x = first$x[match(agg$region_id, first$region_id)],
                        y = first$y[match(agg$region_id, first$region_id)])
  res_c <- circular_scan(dat_reg, cases = cases, population = population,
                         region_id = region_id, x = x, y = y,
                         nsim = 19, seed = 42)
  expect_s3_class(res_c, "circular_scan")
  expect_s3_class(res_c, "treess")

  res_t <- tree_scan(ex$dat, cases = cases, node_id = node_id,
                     tree = ex$tree, population = population,
                     nsim = 19, seed = 42)
  expect_s3_class(res_t, "tree_scan")
  expect_s3_class(res_t, "treess")

  res_s <- sequential_scan(ex$dat, cases = cases, population = population,
                           region_id = region_id, x = x, y = y,
                           node_id = node_id, tree = ex$tree,
                           nsim = 19, seed = 42, max_iter = 2, verbose = FALSE)
  expect_s3_class(res_s, "sequential_scan")
  expect_s3_class(res_s, "treess")
})

test_that("most_likely_cluster() returns the internal MLC", {
  res <- make_ts()
  mlc <- most_likely_cluster(res)
  expect_type(mlc, "list")
  expect_identical(mlc, res$most_likely_cluster)
  expect_true(all(c("node_id", "leaf_ids", "region_ids",
                    "cases", "llr") %in% names(mlc)))
})

test_that("secondary_clusters() dispatches per class", {
  res <- make_ts()
  sc <- secondary_clusters(res)
  expect_s3_class(sc, "data.frame")
  expect_identical(sc, res$secondary_clusters)

  ex <- make_dat()
  res_t <- tree_scan(ex$dat, cases = cases, node_id = node_id,
                     tree = ex$tree, population = population,
                     nsim = 19, seed = 42)
  expect_identical(secondary_clusters(res_t), res_t$significant_cuts)
})

test_that("pvalue() returns scalar for single scans", {
  res <- make_ts()
  p <- pvalue(res)
  expect_type(p, "double")
  expect_length(p, 1L)
  expect_identical(p, res$pvalue)
  expect_true(p > 0 && p <= 1)
})

test_that("sequential_scan accessors handle per-iteration structure", {
  ex <- make_dat()
  res_s <- sequential_scan(ex$dat, cases = cases, population = population,
                           region_id = region_id, x = x, y = y,
                           node_id = node_id, tree = ex$tree,
                           nsim = 19, seed = 42, max_iter = 2, verbose = FALSE)
  mlc <- most_likely_cluster(res_s)
  if (nrow(res_s$clusters) > 0L) {
    expect_type(mlc, "list")
    expect_identical(mlc$iteration, res_s$clusters$iteration[1L])
    p <- pvalue(res_s)
    expect_length(p, nrow(res_s$clusters))
    expect_named(p)
    sc <- secondary_clusters(res_s)
    expect_s3_class(sc, "data.frame")
    expect_equal(nrow(sc), max(0L, nrow(res_s$clusters) - 1L))
  } else {
    expect_null(mlc)
  }
})

test_that("filter_clusters() is a generic with informative default", {
  expect_error(filter_clusters(list()), "must be of class")
  res <- make_ts()
  fc <- filter_clusters(res)
  expect_s3_class(fc, "data.frame")
})

test_that("tree_scan print reports unsupplied population honestly", {
  ex <- make_dat()
  res_no_pop <- tree_scan(ex$dat, cases = cases, node_id = node_id,
                          tree = ex$tree, nsim = 19, seed = 42)
  expect_false(res_no_pop$population_supplied)
  out <- paste(capture.output(print(res_no_pop)), collapse = "\n")
  expect_match(out, "not supplied")
  expect_false(grepl("Total population:", out))

  res_pop <- tree_scan(ex$dat, cases = cases, node_id = node_id,
                       tree = ex$tree, population = population,
                       nsim = 19, seed = 42)
  expect_true(res_pop$population_supplied)
  out2 <- paste(capture.output(print(res_pop)), collapse = "\n")
  expect_match(out2, "Total population:")
})

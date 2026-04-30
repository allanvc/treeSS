test_that("validate_regions catches missing columns", {
  bad_regions <- data.frame(id = 1:3, pop = c(100, 200, 300))
  expect_error(treeSS:::.validate_regions(bad_regions), "missing required columns")
})

test_that("validate_regions catches negative population", {
  regions <- data.frame(
    region_id = 1:3, population = c(100, -1, 300), x = 1:3, y = 1:3
  )
  expect_error(treeSS:::.validate_regions(regions), "non-negative")
})

test_that("validate_regions catches single region", {
  regions <- data.frame(region_id = 1, population = 100, x = 1, y = 1)
  expect_error(treeSS:::.validate_regions(regions), "At least 2")
})

test_that("validate_tree catches missing columns", {
  bad_tree <- data.frame(id = 1:3, parent = c(NA, 1, 1))
  expect_error(treeSS:::.validate_tree(bad_tree), "missing required columns")
})

test_that("validate_tree catches missing root", {
  tree <- data.frame(node_id = 1:3, parent_id = c(2, 3, 1))
  expect_error(treeSS:::.validate_tree(tree), "root node")
})

test_that("get_leaves works correctly", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5),
    parent_id = c(NA, 1, 1, 2, 2)
  )
  leaves <- treeSS:::.get_leaves(tree)
  expect_setequal(leaves, c(3, 4, 5))
})

test_that("get_descendants works correctly", {
  tree <- data.frame(
    node_id   = c(1, 2, 3, 4, 5, 6, 7),
    parent_id = c(NA, 1, 1, 2, 2, 3, 3)
  )
  desc <- treeSS:::.get_descendants(tree, 2)
  expect_setequal(desc, c(2, 4, 5))
  desc_root <- treeSS:::.get_descendants(tree, 1)
  expect_setequal(desc_root, 1:7)
})

test_that("poisson_llr returns 0 when cases <= expected", {
  expect_equal(treeSS:::.poisson_llr(5, 100, 100, 1000), 0)
})

test_that("poisson_llr returns positive for excess cases", {
  llr <- treeSS:::.poisson_llr(50, 100, 100, 1000)
  expect_gt(llr, 0)
})

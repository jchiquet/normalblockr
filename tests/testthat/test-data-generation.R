###############################################################################
## Simulation helpers (R/normal_block_data_generation.R).

test_that("generate_blocks() gives one cluster per variable, no empty cluster, and follows alpha", {
  for (seed in 1:30) {
    set.seed(seed)
    C <- normalblockr:::generate_blocks(8, 5)
    expect_true(all(rowSums(C) == 1))
    expect_true(all(colSums(C) > 0))
  }
  set.seed(1)
  C <- normalblockr:::generate_blocks(3000, 3, alpha = c(0.6, 0.3, 0.1))
  expect_equal(colMeans(C), c(0.6, 0.3, 0.1), tolerance = 0.05)
  expect_error(normalblockr:::generate_blocks(2, 3), "more clusters than variables")
})

test_that("a uniform alpha keeps the plain uniform draw", {
  set.seed(3); C <- normalblockr:::generate_blocks(20, 4)
  set.seed(3); groups <- sample(1:4, size = 20, replace = TRUE)
  expect_equal(max.col(C), groups)
})

set.seed(3)

testdata <- generate_normal_block_data(
  n = 80,
  p = 20,
  d = 1,
  q = 3
)

saveRDS(testdata, file = "tests/testthat/testdata/testdata_normal.RDS")

testdata <- generate_normal_block_data(
  n = 80,
  p = 20,
  d = 1,
  q = 3,
  kappa = runif(20, min = .2, max = .6)
)

saveRDS(testdata, file = "tests/testthat/testdata/testdata_normal_zi.RDS")

## Zero-inflated mean-block fixture (test-ZINormalBlockMeanUnknownClusters.R).
## Frozen because generate_normal_block_mean_data() draws its noise with
## MASS::mvrnorm(), whose eigen() decomposition depends on the BLAS/LAPACK:
## the same seed gives different data on different machines.
set.seed(1234)
ex  <- generate_normal_block_mean_data(n = 100, p = 20, d = 2, q = 3)
Yzi <- ex$Y; Yzi[runif(length(Yzi)) < 0.25] <- 0
saveRDS(list(Y = ex$Y, Yzi = Yzi, X = ex$X, C = ex$parameters$C),
        file = "tests/testthat/testdata/testdata_normal_mean_block_zi.RDS")

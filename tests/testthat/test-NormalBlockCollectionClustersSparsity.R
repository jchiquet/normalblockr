###############################################################################
## Collections crossing a range of q with a sparsity path at each q
## (R/NormalBlockCollectionClustersSparsity.R), for both families.
test_that("q and sparsity can be crossed in both families", {
  for (family in c("var", "mean")) {
    td <- readRDS(if (family == "var") "testdata/testdata_normal.RDS" else "testdata/testdata_normal_mean_block.RDS")
    data <- NormalBlockData$new(td$Y, td$X)
    set.seed(1)
    coll <- normal_block(data, blocks = 2:3, sparsity = TRUE, model = family,
                         control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
    prefix <- if (family == "var") "NormalBlockVar" else "NormalBlockMean"
    expect_s3_class(coll, paste0(prefix, "CollectionClustersSparsity"))
    expect_equal(coll$q_list, c(2, 3))
    expect_equal(nrow(coll$criteria), 6)
    expect_equal(coll$criteria$q, rep(2:3, each = 3))
    path <- coll$get_model(2)
    expect_s3_class(path, paste0(prefix, "CollectionSparsity"))
    lambda <- coll$get_model(3)$sparsity[2]
    expect_equal(coll$get_model(3, lambda)$sparsity, lambda)
    best <- coll$get_best_model("BIC")
    expect_equal(best$BIC, min(coll$criteria$BIC))
    expect_s3_class(coll$plot("BIC"), "ggplot")
    expect_equal(nrow(coll$sparsity), 2)
  }
})

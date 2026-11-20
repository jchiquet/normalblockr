###############################################################################
## The network of a model: latent_network() and plot_network().
test_that("plot_network() accepts force-directed layouts despite signed edge weights", {
  set.seed(3)
  ex    <- generate_normal_block_var_data(n = 100, p = 20, d = 1, q = 4)
  data  <- NormalBlockData$new(ex$Y, ex$X)
  model <- normal_block(data, blocks = ex$parameters$C, control = NB_control(verbose = FALSE))
  pdf(NULL)
  on.exit(dev.off())
  expect_true(any(model$latent_network("partial_cor") < 0))
  expect_no_error(model$plot_network(layout = igraph::layout_with_fr))
  expect_no_error(model$plot_network(layout = igraph::layout_in_circle, remove.isolated = TRUE))
  coords <- matrix(rnorm(2 * model$q), model$q, 2)
  expect_no_error(model$plot_network(layout = coords))
})

test_that("latent_network() extracts the network without plotting it", {
  td <- readRDS("testdata/testdata_normal.RDS")
  set.seed(1)
  model <- normal_block(NormalBlockData$new(td$Y, td$X), 3, sparsity = 0.05, control = NB_control(verbose = FALSE))
  net <- model$latent_network()
  expect_equal(dim(net), c(model$q, model$q))
  expect_equal(unname(diag(as.matrix(net))), rep(1, model$q))
  support <- model$latent_network("support")
  expect_equal(sum(as.matrix(support)[upper.tri(support)]), model$n_edges)
  expect_equal(as.matrix(model$latent_network("precision")), model$model_par$Omega, ignore_attr = TRUE)
  G <- model$plot_network(output = "corrplot", plot = FALSE)
  expect_equal(G, net)
})

test_that("plot_network() does not error when the inferred clustering has an empty cluster", {
  set.seed(7)
  ex    <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 3)
  model <- NormalBlockVarUnknownClusters$new(NormalBlockData$new(ex$Y, ex$X), 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 2, threshold = -1))
  ## column 3 of the memberships never wins the argmax, as happens on real
  ## data when q exceeds the number of clusters the VEM actually uses
  model$.__enclos_env__$private$C <- cbind(matrix(c(0.6, 0.4), model$p, 2, byrow = TRUE), 0)
  expect_length(tabulate(model$clustering, nbins = model$q), model$q)
  expect_no_error(model$plot_network(plot = FALSE))
})

test_that("the criteria plots name their x axis after the explored quantity", {
  set.seed(3)
  ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
  data <- NormalBlockData$new(ex$Y, ex$X)
  ctrl <- NB_control(verbose = FALSE, refine = FALSE)
  co <- normal_block(data, blocks = 2:4, control = ctrl)
  expect_equal(co$plot()$labels$x, "q")
  sp <- normal_block(data, blocks = ex$parameters$C, sparsity = TRUE, control = ctrl)
  expect_equal(sp$plot()$labels$x, "sparsity")
})

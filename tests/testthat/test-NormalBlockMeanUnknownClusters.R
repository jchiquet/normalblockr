###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
q <- ncol(C)
B <- testdata$parameters$B


test_that("normal block mean with unknown clusters - basic tests", {
  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockMeanUnknownClusters$new(data, q, 0.1)
  model$optimize()
  expect_lt(model$BIC, 11224)
})

test_that("predict() maps the cluster-level predictor back to the variables", {
  data  <- NormalBlockData$new(Y, X)
  model <- normal_block(data, blocks = q, model = "mean", control = NB_control(verbose = FALSE))
  pred  <- predict(model, X)
  expect_equal(dim(pred), dim(Y))               # n x p, not n x q
  expect_equal(unname(pred), unname(model$fitted))
})

test_that("NB_control(heuristic = TRUE) returns the moment-based fit", {
  data  <- NormalBlockData$new(Y, X)
  model <- normal_block(data, blocks = q, model = "mean",
                        control = NB_control(verbose = FALSE, heuristic = TRUE))
  expect_equal(model$inference_method, "heuristic")
  expect_true(is.na(model$loglik))
  expect_length(model$clustering, ncol(Y))
})

test_that("NB_control(fixed_tau = TRUE) leaves tau at its initial value", {
  data <- NormalBlockData$new(Y, X)
  init <- NormalBlockMeanUnknownClusters$new(data, q, control = NB_control(verbose = FALSE)
            )$.__enclos_env__$private$optim_initialize()
  model <- NormalBlockMeanUnknownClusters$new(data, q,
             control = NB_control(verbose = FALSE, clustering_init = get_clusters(init$tau)))
  model$optimize(NB_control(verbose = FALSE, fixed_tau = TRUE))
  expect_equal(model$var_par$tau, init$tau)
})

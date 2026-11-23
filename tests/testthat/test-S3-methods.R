###############################################################################
## The S3 methods of the models (print, summary, plot, logLik, BIC), for both
## families, and of the collections (print, summary, logLik, BIC).
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

data  <- NormalBlockData$new(Y, X)
set.seed(1)
model <- NormalBlockVarUnknownClusters$new(data, q)
model$optimize()

set.seed(11)
ex_small    <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 2)
data_small  <- NormalBlockData$new(ex_small$Y, ex_small$X)
fast        <- NB_control(niter = 5, threshold = -1, verbose = FALSE)
collection  <- normal_block(data_small, blocks = 2:4, control = fast)

test_that("logLik.NormalBlockVarBase() returns a well-formed logLik object", {
  ll <- logLik(model)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), model$loglik)
  expect_equal(attr(ll, "df"), model$nb_param)
  expect_equal(attr(ll, "nobs"), model$n)
})

test_that("BIC.NormalBlockVarBase() matches model$BIC", {
  expect_equal(BIC(model), model$BIC)
})

test_that("print.NormalBlockVarBase() reports the model type and criteria", {
  expect_output(print(model), model$who_am_I, fixed = TRUE)
  expect_output(print(model), "Useful S3 methods")
  capture.output(expect_invisible(print(model)))
  ## the model type ends with a period, without a space
  expect_false(any(grepl(" \\.$", capture.output(print(model)))))
  expect_false(any(grepl(" \\.$", capture.output(print(summary(model))))))
})

test_that("summary.NormalBlockVarBase() returns cluster sizes, edges and density", {
  s <- summary(model)
  expect_s3_class(s, "summary.NormalBlockBase")
  expect_equal(s$who_am_I, model$who_am_I)
  expect_equal(s$criteria, model$criteria)
  expect_equal(s$cluster_sizes, model$cluster_sizes)
  expect_equal(s$n_edges, model$n_edges)
  expect_equal(s$density, model$n_edges / choose(model$q, 2))
  expect_output(print(s), "Cluster sizes")
  expect_output(print(s), "Network")
})

test_that("plot.NormalBlockVarBase() returns the underlying ggplot object", {
  p <- plot(model)
  expect_s3_class(p, "ggplot")
})

test_that("logLik.NormalBlockCollection() and BIC.NormalBlockCollection() return one value per model", {
  ll  <- logLik(collection)
  bic <- BIC(collection)
  expect_length(ll, length(collection$models))
  expect_length(bic, length(collection$models))
  expect_equal(ll, collection$criteria$loglik)
  expect_equal(bic, collection$criteria$BIC)
})

test_that("collection$loglik raises an informative error instead of silently returning NULL", {
  expect_error(collection$loglik, "not defined for a collection")
  expect_error(collection$loglik, "logLik", fixed = TRUE)
})

test_that("print.NormalBlockCollection() reports the model type and the q range explored", {
  expect_output(print(collection), collection$who_am_I, fixed = TRUE)
  expect_output(print(collection), "model(s) explored", fixed = TRUE)
  expect_output(print(collection), "q ranging from 2 to 4")
  capture.output(expect_invisible(print(collection)))
})

test_that("summary.NormalBlockCollection() returns the full criteria table and q range", {
  s <- summary(collection)
  expect_s3_class(s, "summary.NormalBlockCollection")
  expect_equal(s$who_am_I, collection$who_am_I)
  expect_equal(s$criteria, collection$criteria)
  expect_equal(s$q_range, range(collection$criteria$q))
  expect_output(print(s), "q ranging from 2 to 4")
  expect_output(print(s), "loglik")
})

test_that("print(), summary() and plot() work on a mean-block model, its network relating the variables", {
  td <- readRDS("testdata/testdata_normal_mean_block.RDS")
  set.seed(1)
  m <- normal_block(NormalBlockData$new(td$Y, td$X), 3, sparsity = 0.1, model = "mean",
                    control = NB_control(verbose = FALSE))
  expect_output(print(m), m$who_am_I, fixed = TRUE)
  s <- summary(m)
  expect_equal(s$cluster_sizes, m$cluster_sizes)
  expect_gt(s$n_edges, 0)
  expect_equal(s$density, m$n_edges / choose(m$p, 2))
  expect_output(print(s), "density")
  expect_s3_class(plot(m), "ggplot")
})

test_that("plot() of a collection names the model selected by each criterion", {
  set.seed(3)
  ex <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
  coll <- normal_block(NormalBlockData$new(ex$Y, ex$X), blocks = 2:4,
                       control = NB_control(verbose = FALSE, refine = FALSE))
  pdf(NULL); on.exit(dev.off())
  p <- plot(coll, criteria = c("deviance", "BIC", "ICL"))
  labels <- unlist(lapply(p$layers, function(l) l$data$label))
  expect_setequal(labels, paste0(c("BIC", "ICL"), ": q = ",
                                 c(coll$get_best_model("BIC")$q, coll$get_best_model("ICL")$q)))
})

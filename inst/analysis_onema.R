library(normalblockr)

data(onema)

# select and format data to run normalblockr on it
X <- model.matrix( ~ 1, data = onema$covariates)
Y <- log(1 + onema$biomass) |> as.matrix()
data <- normalblockr::NormalBlockData$new(Y, X)

## ZI inflated normal with diagonal covariance
out <- normal_block(data, blocks = 2:10, zero_inflation = TRUE, control = NB_control(clustering_init = "ward2"))

out$plot(criteria = c("BIC", "EBIC", "deviance"))
myModel <- out$get_best_model("EBIC")

mySparseModel <- normal_block(data, blocks = myModel$q, sparsity = TRUE, zero_inflation = TRUE))
mySparseModel$plot(c("EBIC", "BIC"))
mySparseModel$get_best_model("EBIC")$plot_network()
mySparseModel$get_best_model("EBIC")$clustering

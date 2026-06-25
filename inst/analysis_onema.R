library(normalblockr)

data(onema)

# select and format data to run normalblockr on it
X <- model.matrix( ~ 1 + temperature_med, data = onema$covariates)
Y <- log(1 + onema$biomass) |> as.matrix()
data <- NormalBlockData$new(Y, X)

## ZI inflated normal with diagonal covariance
out <- normal_block(data, blocks = 2:10, zero_inflation = TRUE)

out$plot(criteria = c("BIC", "EBIC"))
myModel <- out$get_best_model("EBIC")

## Groups of species
myModel$elements_per_cluster

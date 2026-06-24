library(reshape2)

## loading data
onema_community <- read.csv("community_data.csv")
env             <- read.csv("environment.csv")
protocol        <- read.csv("fishing_protocol.csv") # link between opcods & stations

onema_community$species <- factor(onema_community$species,
                                  levels = unique(onema_community$species))
onema_community$opcod   <- factor(onema_community$opcod,
                                  levels = unique(onema_community$opcod))

# filter data for stations with environmental info
protocol        <- protocol[protocol$station %in% env$station,]
onema_community <- onema_community[onema_community$opcod %in% protocol$opcod,] # does not contain NA
species_order <- unique(onema_community$species)
opcod_order   <- unique(onema_community$opcod)
stations        <- protocol[match(opcod_order, protocol$opcod), , drop=FALSE]$station
idx             <- match(stations, env$station)
stations_env    <- env[idx, ,drop = FALSE] # does not contain NA

# select and format data to run normalblockr on it
X    <- stations_env[ , 2:ncol(stations_env)]
Y    <- setNames(apply(as.matrix(dcast(onema_community, species ~ opcod, value.var = "biomass", fill = 0))[,-1], 1, as.numeric), species_order)
Y    <- log(1 + Y)
n <- nrow(Y)
Xones <- matrix(rep(1, n), nrow = n)

## ZI inflated normal with diagonal covariance
data <- normalblockr::NormalBlockData$new(Y, Xones)
out <- normal_block(data, blocks = 2:10, zero_inflation = TRUE)

out$plot(criteria = c("BIC", "EBIC"))
myModel <- out$get_best_model("BIC")

mySparseModel <- normal_block(data, blocks = myModel$q, sparsity = TRUE, zero_inflation = TRUE)

mySparseModel$plot(c("EBIC", "BIC"))
mySparseModel$get_best_model("EBIC")$plot_network()
mySparseModel$get_best_model("EBIC")$clustering

## response is approximation error
## predictors are network (control), dynamics (control), and node set type (predictor of interest)
## something like aov(error ~ network + dynamics + nodeset, data = df)

if(interactive()) {
    if(getwd() != "/user/neilmacl/Documents/reduction/analysis/") {
        setwd("/user/neilmacl/Documents/reduction/analysis/")
    }
}

library(igraph)

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
dynamics <- c("dw", "SIS", "mutualistic", "genereg")

### this is a function
gather_data <- function(network, dynamic) {
## load optimized node sets -- these are stored by network/dynamics combination
## generate random node sets (four types)
## make into data frame
## data frame as a var for approx error, then vars for network, dynamics, and node set type
## return the data frame

### lapply this to the right combination of parameters, then do.call(rbind, ...)
## collect all data frames

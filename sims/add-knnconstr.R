library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

networks <- c("dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr")
dynamics <- c("doublewell", "SIS", "genereg", "mutualistic")
ntrials <- 100
optweights <- c(FALSE, TRUE)

## This function re-saves the .rds file with knnconstr node sets concatenated to the end
add_knnconstr <- function(network, dynamics, ntrials, optweights) {
    nodesetfile <- paste0(
        "../data/ns-",
        network, "_",
        dynamics,
        switch(optweights + 1, "", "_w"),
        ".rds")
    fullstatefile <- paste0("../data/fullstate-", network, ".rds")

    g <- readRDS(paste0("../data/", network, ".rds"))
    N <- vcount(g)
    AL <- as_adj_list(g, "all")
    k <- degree(g)
    knn <- knn(g)$knn
    
    Y <- readRDS(fullstatefile)[[dynamics]]
    y <- rowMeans(Y)
    n <- floor(log(N))

    nodesets <- readRDS(nodesetfile)
    nodesets$knnconstr <- make_dataset(
        ntrials = ntrials, ns.type = "knnconstr", ncores = ncores,
        n = n, g = g, y = y, Y = Y, optimize_weights = optweights
    )

    saveRDS(nodesets, file = nodesetfile)
}

for(network in networks) {
    for(dynamic in dynamics) {
        for(optweight in optweights) {
            add_knnconstr(network, dynamic, ntrials, optweight)
        }
    }
}

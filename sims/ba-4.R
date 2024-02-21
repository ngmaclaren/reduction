library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

g <- readRDS("../data/ba.rds")
N <- vcount(g)
Y <- readRDS("../data/fullstate-ba.rds")$doublewell
y <- rowMeans(Y)

ns <- c(1:3, 6)

ntrials <- 100

opts <- lapply(ns, function(n) {
    make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
})

saveRDS(opts, file = "../data/ba-4.rds")

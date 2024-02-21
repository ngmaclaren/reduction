library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

g <- readRDS("../data/dolphin.rds")
N <- vcount(g)
Y <- readRDS("../data/fullstate-dolphin.rds")$doublewell
y <- rowMeans(Y)

ns <- 1:12

ntrials <- 100

opts <- lapply(ns, function(n) {
    make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
})

saveRDS(opts, file = "../data/dolphin-1-to-12.rds")

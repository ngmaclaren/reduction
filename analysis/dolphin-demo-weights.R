library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(optNS)
ns <- 1:4

g <- readRDS("../data/dolphin.rds")
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
Y <- readRDS("../data/fullstate-dolphin.rds")$doublewell
y <- rowMeans(Y)

solns <- mclapply(ns, function(n) select_optimized(n, g, y, Y, optimize_weights = TRUE), mc.cores = ncores)

save.image("dolphin-demo-weights.RData")

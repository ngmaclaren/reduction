## This file will output an .RData file

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(optNS)

doublewell1D <- function(t, x, params) {
    with(params, {
        dx <- -(x - r[1]) * (x - r[2]) * (x - r[3]) + D*b*x
        return(list(c(dx)))
    })
}

ns <- 1:4

g <- readRDS("../data/dolphin.rds")
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
Y <- readRDS("../data/fullstate-dolphin.rds")$doublewell
y <- rowMeans(Y)

solns <- mclapply(ns, function(n) select_optimized(n, g, y, Y), mc.cores = ncores)

### GBB and DART ###
control <- list(times = 0:15, ncores = ncores)

## GBB
GBB.b <- mean(k^2)/mean(k)
                                        # this is the theory
GBB <- solve_in_range( 
    .doublewell$Ds, "D", doublewell1D, .doublewell$xinit.low,
    params = c(.doublewell, list(b = GBB.b)), control = control
)
                                        # and this the ground truth (simulated, numerical)
GBB.obs <- apply(Y, 1, function(x) mean(k*x)/mean(k))

## DART
eigs <- eigen(A, symmetric = TRUE)
DART.a <- eigs$vectors[, 1]/sum(eigs$vectors[, 1])
DART.b <- eigs$values[1]
                                        # theory
DART <- solve_in_range(
    .doublewell$Ds, "D", doublewell1D, .doublewell$xinit.low,
    params = c(.doublewell, list(b = DART.b)), control = control
)
                                        # ground truth
DART.obs <- apply(Y, 1, function(x) sum(DART.a*x))

save.image("dolphin-demo.RData")

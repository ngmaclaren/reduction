## Exp 1: vary n
## Exp 2: vary maxit
##
## Choose two large networks. How about Spanish and FlyBi.
## For now, double-well only
## Vary n \in {1, ... 20}
## Get the time complexity.

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)
library(sfsmisc)

## initial trial, with a small network
network <- "proximity"
dynamics <- "doublewell"

fullstatefile <- paste0("../data/fullstate-", network, ".rds")

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
AL <- as_adj_list(g, "all")

Y <- readRDS(fullstatefile)[[dynamics]]
y <- rowMeans(Y)

ns <- 1:20
                                        # This is parallelized, but I don't think this is a good idea
## runtimes <- mclapply(ns, function(n) system.time(select_optimized(n, g, y, Y))[3], mc.cores = ncores)
runtimes <- numeric(length(ns))
results <- vector("list", length(ns))

for(i in seq_along(ns)) {
    n <- ns[i]
    runtime <- system.time(
        result <- select_optimized(n, g, y, Y)
    )
    runtimes[i] <- runtime[3]
    results[[i]] <- result
}
rm(i, n, runtime)

par(mar = c(5, 5, 1, 1), pty = "s")
plot(
    ns, runtimes, xlab = "", ylab = "", axes = FALSE,
    cex = 2
)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5)
title(xlab = "n", font.lab = 3, cex.lab = 2)
title(ylab = "Wall time", cex.lab = 2)


## Now, vary maxit
maxits <- seq(10, 100, by = 10)*N
n <- 20
runtimes <- numeric(length(maxits))
results <- vector("list", length(maxits))

for(i in seq_along(maxits)) {
    maxit <- maxits[i]
    runtime <- system.time(
        result <- select_optimized(n, g, y, Y, maxit = maxit)
    )
    runtimes[i] <- runtime[3]
    results[[i]] <- result
}
rm(i, maxit, runtime)
errors <- get_error(results)

par(mar = c(5, 5, 1, 1), pty = "s")
plot(
    maxits, errors, xlab = "", ylab = "", axes = FALSE, log = "y",
    cex = 2
)
abline(v = 50*N, col = 2, lwd = 2)
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5, las = 0)
title(xlab = expression(italic(h)[max]), cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2)
legend("topright", bty = "n", legend = "50N", lwd = 2, col = 2)

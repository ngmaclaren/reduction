library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)
palette("R4")
set.seed(1234)

network <- "dolphin"
dynamics <- "doublewell"

fullstatefile <- paste0("./data/fullstate-", network, ".rds")

g <- readRDS(paste0("./data/", network, ".rds"))
N <- vcount(g)

Y <- readRDS(fullstatefile)[[dynamics]]
y <- rowMeans(Y)
n <- floor(log(N))

system.time(
    soln <- select_optimized(n, g, y, Y)
)
str(soln) # list with vertex indices, approximation error, and degree sequence; optimized node weights are null
Z <- Y[, soln$vs]
z <- rowMeans(Z)
Ds <- seq(0, 1, length.out = 100)

matplot(
    Ds, Y,
    type = "l", lty = 1, lwd = 0.25, col = "gray50",
    xlab = "D", ylab = expression(x[i]), font.lab = 3
)
lines(Ds, y, lty = 1, lwd = 8, col = 1)
matlines(Ds, Z, lty = 1, lwd = 4, col = 2)
lines(Ds, z, lty = 1, lwd = 8, col = 3)

## to make several node sets:
ntrials <- 10

system.time(
    opts <- make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
)

best <- opts[[which.min(get_error(opts))]]

rands <- make_dataset(
    ntrials = 10*ntrials, ns.type = "rand", ncores = ncores,
    n = n, g = g, y = y, Y = Y
)

error <- list(
    opt = get_error(opts),
    rand = get_error(rands)
)

xlim <- range(unlist(error))

par(mar = c(4, 8, 1, 1))
plot(NULL, xlim = xlim, ylim = c(0.5, 2.5), xlab = "Approximation error", ylab = "", log = "x", yaxt = "n")
points(error$opt, jitter(rep(2, length(error$opt)), amount = 0.1), pch = 1, col = 2)
points(error$rand, jitter(rep(1, length(error$rand)), amount = 0.1), pch = 0, col = 1)
axis(2, at = c(2, 1), tick = FALSE, labels = c("Optimized", "Random"), las = 2)

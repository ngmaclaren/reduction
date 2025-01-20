## On the cluster, need to use a larger network and draw multiple node sets per network.

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(sdn)
library(optNS)
library(sfsmisc)

network <- "dolphin"
dynamics <- "doublewell"

fullstatefile <- paste0("../data/fullstate-", network, ".rds")
nsfile <- paste0("../data/ns-", network, "_", dynamics, ".rds")

                                        # this is the original, against which we'll test
g <- readRDS(paste0("../data/", network, ".rds")) 
N <- vcount(g)
M <- ecount(g)

                                        # the test dynamics results
Y <- readRDS(fullstatefile)[[dynamics]] 
y <- rowMeans(Y)
Ds <- seq(0, 1, length.out = 100)
n <- floor(log(N))
nodesets <- readRDS(nsfile)[c("opt", "rand")]

control <- list(times = 0:15, ncores = ncores)

## First, remove 10%, 20%, and 30% of edges
toremove <- floor(M * seq(0.05, 0.5, by = 0.05)) # remove this many edges from g
removethese <- lapply(toremove, function(x) sample(seq(M), size = x)) # a list of integer vectors. These integer vectors identify edges to remove from g

glist <- lapply(removethese, function(x) delete_edges(g, x))
                                        # check
## sapply(glist, vcount)
## sapply(glist, ecount)
## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(g. in glist) plot(g., vertex.size = 5, vertex.label = "")

## Now, need ground truth sims for glist
Ylist <- lapply(
    glist, function(g.) {
        with(list(
            params = c(.doublewell, list(AL = as_adj_list(g., "all")))
        ), {
            solve_in_range(
                params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode",
                method = "adams", maxsteps = 20000
            )
        })
    }
)

## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(df in Ylist) bifplot(df, Ds, TRUE, lwd = 0.5, col = 1)

uncertain <- mapply(function(g., Y.) select_optimized(n, g., rowMeans(Y.), Y.), glist, Ylist, SIMPLIFY = FALSE) 
uncertain_error <- sapply(uncertain, function(ns) obj_fn(ns$vs, y, Y))

rand_error <- get_error(nodesets$rand)
opt_error <- get_error(nodesets$opt)
ylim <- range(c(uncertain_error, rand_error, opt_error))
optcolor <- "#3584e4"
randcolor <- "#33d17a"

dev.new(width = 8)
lyt <- matrix(c(
    1, 1, 1, 2,
    1, 1, 1, 2), byrow = TRUE, nrow = 2)
layout(lyt)
par(mar = c(5, 6, 1, 0))
plot(toremove/M, uncertain_error, log = "y",
     cex = 2, lwd = 2, ylim = ylim,
     axes = FALSE, xlab = "", ylab = "")
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5, n.axp = 1)
title(xlab = "Proportion of edges removed", cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2, line = 3.75)
legend("topright", bty = "n", pch = c(1, 0, 1), col = c(optcolor, randcolor, 1), cex = 1.5,
       pt.lwd = 2, pt.cex = 2,
       legend = c("Optimized, full info", "Random", "Optimized, partial info"))
par(mar = c(5, 0, 1, 1))
plot(
    jitter(rep(1, 100), amount = 0.1), rand_error, log = "y", ylim = ylim, cex = 2, lwd = 2, pch = 0,
    col = randcolor, axes = FALSE, xlab = "", ylab = ""
)
points(jitter(rep(1, 100), amount = 0.1), opt_error, cex = 2, pch = 1, col = optcolor, lwd = 2)




## Removing nodes shouldn't be too much harder.
tokeep <- ceiling(N * seq(0.95, 0.5, by = -0.05))
keepthese <- lapply(tokeep, function(x) sample(seq(N), size = x))
glist <- lapply(keepthese, function(x) subgraph(g, x))
ns <- sapply(glist, function(g.) floor(log(vcount(g.))))

                                        # check
## sapply(glist, vcount)
## sapply(glist, ecount)
## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(g. in glist) plot(g., vertex.size = 5, vertex.label = "")

Ylist <- lapply(
    glist, function(g.) {
        with(list(
            params = c(.doublewell, list(AL = as_adj_list(g., "all"))),
            control = list(times = 0:15, ncores = ncores)
        ), {
            solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode",
                           method = "adams", maxsteps = 20000)
        }
        )
    }
)

## dev.new(height = 6, width = 15)
## par(mfrow = c(2, 5), pty = "s")
## for(df in Ylist) bifplot(df, Ds, TRUE, lwd = 0.5, col = 1)

uncertain <- mapply(
    function(n., g., Y.) select_optimized(n., g., rowMeans(Y.), Y.),
    ns, glist, Ylist,
    SIMPLIFY = FALSE
)

uncertain_error <- sapply(uncertain, function(ns) obj_fn(ns$vs, y, Y))

rand_error <- get_error(nodesets$rand)
opt_error <- get_error(nodesets$opt)
ylim <- range(c(uncertain_error, rand_error, opt_error))
optcolor <- "#3584e4"
randcolor <- "#33d17a"
dev.new(width = 8)
lyt <- matrix(c(
    1, 1, 1, 2,
    1, 1, 1, 2), byrow = TRUE, nrow = 2)
layout(lyt)
par(mar = c(5, 6, 1, 0))
plot((N - tokeep)/N, uncertain_error, log = "y",
     cex = 2, lwd = 2, ylim = ylim,
     axes = FALSE, xlab = "", ylab = "")
box()
eaxis(1, cex.axis = 1.5)
eaxis(2, cex.axis = 1.5, n.axp = 1)
title(xlab = "Proportion of nodes removed", cex.lab = 2)
title(ylab = "Approximation error", cex.lab = 2, line = 3.75)
legend("topright", bty = "n", pch = c(1, 0, 1), col = c(optcolor, randcolor, 1), cex = 1.5,
       pt.lwd = 2, pt.cex = 2,
       legend = c("Optimized, full info", "Random", "Optimized, partial info"))
par(mar = c(5, 0, 1, 1))
plot(
    jitter(rep(1, 100), amount = 0.1), rand_error, log = "y", xlim = c(0.75, 1.25), ylim = ylim,
    cex = 2, lwd = 2, pch = 0,
    col = randcolor, axes = FALSE, xlab = "", ylab = ""
)
points(jitter(rep(1, 100), amount = 0.1), opt_error, cex = 2, pch = 1, col = optcolor, lwd = 2)

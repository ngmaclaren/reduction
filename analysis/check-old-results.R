if(getwd() != "/user/neilmacl/Documents/reduction/analysis") if(interactive()) setwd("./analysis")

library(parallel)
ncores <- detectCores()-1
library(igraph)
source("../src/functions.R")

get_error <- function(dl) sapply(dl, `[[`, "error")
get_bparam <- function(dynamics) {
    switch(
        dynamics,
        dw = doublewell_parms$Ds,
        SIS = SIS_parms$Ds,
        mutualistic = mutualistic_parms$Ds,
        genereg = genereg_parms$Ds
    )
}
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]

    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])

    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos

    text(xpos, ypos, label, ...)
}
network <- "proximity"
A.dynamics <- "dw" # "genereg"
B.dynamics <- "SIS"

load(paste0("../data/", network, ".rda"))
g <- upgrade_graph(get(network))
k <- degree(g)
N <- vcount(g)
load(paste0("../data/fullstate-", network, ".rda"))

filedir <- "../data/optimized-nodesets/"
filenames <- list.files(filedir)
opts <- lapply(
    filenames, function(filename) readRDS(paste0(filedir, filename))
)
listnames <- gsub("-", "_", filenames)
listnames <- gsub(".rds", "", listnames)
names(opts) <- listnames

## Fix the network.
## New results: evaluate all five kinds of node sets from one dynamics against those of another.
## Old results: compare the optimized from A, evaluated against B, to the fixed-degree of B.
A <- list(opt = opts[[paste(c(network, A.dynamics), collapse = "_")]])
n <- length(A$opt[[1]]$vs)
ntrials <- length(A$opt)
A.bparam <- get_bparam(A.dynamics)
A.Y <- fullstate[[A.dynamics]]
A.y <- rowMeans(fullstate[[A.dynamics]])
A.bestopt <- A$opt[which.min(get_error(A$opt))][[1]]
A$fixed <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, comps = A.bestopt$vs)
A$rand <- make_dataset(n, ntrials, A.bparam, A.y, A.Y)
A$constr <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, use_connections = TRUE)
A$quant <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, use_connections = TRUE, use_quantiles = TRUE)

B <- list(opt = opts[[paste(c(network, B.dynamics), collapse = "_")]])
B.bparam <- get_bparam(B.dynamics)
B.Y <- fullstate[[B.dynamics]]
B.y <- rowMeans(fullstate[[B.dynamics]])
B.bestopt <- B$opt[which.min(get_error(B$opt))][[1]]
B$fixed <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, comps = B.bestopt$vs)
B$rand <- make_dataset(n, ntrials, B.bparam, B.y, B.Y)
B$constr <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, use_connections = TRUE)
B$quant <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, use_connections = TRUE, use_quantiles = TRUE)


## Replicate new results
xdf <- data.frame(lapply(A, get_error))
ydf <- data.frame(
    lapply(
        A,
        function(dl) {
            sapply(dl, function(l) obj_fn(l$vs, B.y, B.Y, B.bparam))
        }
    )
)

oldcomps <- get_error(B$fixed)

xlim <- range(c(unlist(xdf), .9*min(xdf$opt)))
ylim <- range(unlist(ydf))

pdf("test.pdf")
palette("Set 1")
plot(
    NULL, xlim = xlim, ylim = ylim, log = "xy",
    xlab = A.dynamics, ylab = B.dynamics, main = paste(network, "error")
)
legend(
    "bottomright", col = c(1:5, "black"), pch = c(rep(1, 5), 4), pt.lwd = 2, bty = "n",
    legend = c("Optimized", "Fixed-degree", "Random", "Constrained", "Quantiled", "Old comps")
)
for(i in seq_along(xdf)) points(xdf[[i]], ydf[[i]], col = adjustcolor(i, .5), pch = 1)
for(i in seq_along(xdf)) points(mean(xdf[[i]]), mean(ydf[[i]]), col = i, pch = 1, cex = 2, lwd = 3)
points(rep(.9*min(xdf$opt), length(xdf$opt)), oldcomps, col = adjustcolor("black", .5), pch = 4)
points(.9*min(xdf$opt), mean(oldcomps), col = "black", pch = 4, lwd = 3, cex = 2)
placelabel(
    paste("Old metric:", round(mean(ydf$opt)/mean(get_error(B$fixed)), 2),
          "\n(smaller is better, if =1 then no difference)"),
    .02, .98, adj = c(0, 1)
)
dev.off()

## Replicate old results
## basically, the old results compared A.opt eval on B.dyn vs. B.fixed eval on B.dyn
## summary(ydf$opt)
## summary(get_error(B$fixed))

## mean(ydf$opt)/mean(get_error(B$fixed))

## this file should accept a network and two dynamics and return a plot

library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network to analyze."
    ),
    make_option(
        "--dynamicsA", type = "character", default = "dw",
        help = "The dynamics on which the node sets were optimzied."
    ),
    make_option(
        "--dynamicsB", type = "character", default = "SIS",
        help = "The dynamics on which to evaluate the node sets."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

if(interactive()) {
    ## if(getwd() != "/user/neilmacl/Documents/reduction/analysis") setwd("./analysis")
    network <- "ba"
    A.dynamics <- "dw" # "genereg"
    B.dynamics <- "SIS"
    useall <- "yes" # "no"
    outfile <- paste0(
        "../img/dynamics-comparisons/",
        paste(c(network, A.dynamics, B.dynamics, useall), collapse = "-"),
        ".pdf"
    )
} else {
    network <- args$network
    A.dynamics <- args$dynamicsA
    B.dynamics <- args$dynamicsB
    outfile <- paste0(
        "../img/dynamics-comparisons/",
        paste(c(network, A.dynamics, B.dynamics), collapse = "-"),
        ".pdf"
    )
    useall <- "no"
}

library(parallel)
ncores <- detectCores() - 1
library(igraph)#, lib.loc = "/user/neilmacl/rlocal")
library(deSolve)
library(sfsmisc) # for axis functions
source("../src/functions.R")
RNGkind("L'Ecuyer-CMRG")

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

load(paste0("../data/", network, ".rda"))
g <- upgrade_graph(get(network))
k <- degree(g)
N <- vcount(g)
## vdf <- data.frame(
##     dc = degree(g, normalized = TRUE),
##     bc = betweenness(g, directed = FALSE, normalized = TRUE),
##     cc = closeness(g, normalized = TRUE),
##     ec = eigen_centrality(g, directed = FALSE)$vector
## )
load(paste0("../data/fullstate-", network, ".rda"))

filedir <- "../data/optimized-nodesets/"
filenames <- list.files(filedir)
opts <- lapply(
    filenames, function(filename) readRDS(paste0(filedir, filename))
)
listnames <- gsub("-", "_", filenames)
listnames <- gsub(".rds", "", listnames)
names(opts) <- listnames

A <- list(opt = opts[[paste(c(network, A.dynamics), collapse = "_")]])
n <- length(A$opt[[1]]$vs)
ntrials <- length(A$opt)
A.bparam <- get_bparam(A.dynamics)
A.Y <- fullstate[[A.dynamics]]
A.y <- rowMeans(fullstate[[A.dynamics]])
A.bestopt <- A$opt[which.min(get_error(A$opt))][[1]]
A$fixed <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, comps = A.bestopt$vs)
A$rand <- make_dataset(n, ntrials, A.bparam, A.y, A.Y)
if(useall == "yes") {
    A$constr <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, use_connections = TRUE)
    A$quant <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, use_connections = TRUE, use_quantiles = TRUE)
    A$comm <- make_dataset(n, ntrials, A.bparam, A.y, A.Y, use_connections = TRUE, use_communities = TRUE)
}

B <- list(opt = opts[[paste(c(network, B.dynamics), collapse = "_")]])
B.bparam <- get_bparam(B.dynamics)
B.Y <- fullstate[[B.dynamics]]
B.y <- rowMeans(fullstate[[B.dynamics]])
B.bestopt <- B$opt[which.min(get_error(B$opt))][[1]]
B$fixed <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, comps = B.bestopt$vs)
B$rand <- make_dataset(n, ntrials, B.bparam, B.y, B.Y)
if(useall == "yes") {
    B$constr <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, use_connections = TRUE)
    B$quant <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, use_connections = TRUE, use_quantiles = TRUE)
    B$comm <- make_dataset(n, ntrials, B.bparam, B.y, B.Y, use_connections = TRUE, use_communities = TRUE)
}
                                        # demonstration node sets
## demo_VS <- lapply(vdf, function(cm) order(cm, decreasing = TRUE)[1:n])
## demo_error_A <- lapply(demo_VS, function(vs) obj_fn(vs, A.y, A.Y, A.bparam))
## demo_error_B <- lapply(demo_VS, function(vs) obj_fn(vs, B.y, B.Y, B.bparam))


xdf <- data.frame(lapply(A, get_error))
ydf <- data.frame(
    lapply(
        A,
        function(dl) {
            sapply(dl, function(l) obj_fn(l$vs, B.y, B.Y, B.bparam))
        }
    )
)

## oldcomps <- get_error(B$fixed)

xlim <- range(c(unlist(xdf), .9*min(xdf$opt)))#, unlist(demo_error_A)))
ylim <- range(unlist(ydf))#, unlist(demo_error_B))
axislabels <- c(
    dw = "Double-well", SIS = "SIS", genereg = "Gene-regulatory", mutualistic = "Mutualistic species"
)
labelsize <- 1.5

pdf(outfile)
palette("Set 1")

par(mai = c(1.25, 1.25, 0.25, 0.25))
plot(
    NULL, xlim = xlim, ylim = ylim, log = "xy", axes = FALSE,
    xlab = "", ylab = "" #, main = paste(network, "error")
)
eaxis(1, at.small = FALSE, cex.axis = labelsize*.6)
eaxis(2, at.small = FALSE, cex.axis = labelsize*.6)
title(xlab = axislabels[A.dynamics], line = 3, cex.lab = labelsize)
title(ylab = axislabels[B.dynamics], line = 4.5, cex.lab = labelsize)
box()
legend(
    "bottomright", col = seq_along(B), pch = seq_along(B), pt.lwd = 2, bty = "n",
    legend = switch(
        useall,
        yes = c("Optimized", "Fixed-degree", "Random", "Constrained", "Quantiled", "Community-based"),
        no = c("Optimized", "Fixed-degree", "Random")
    )
)
for(i in seq_along(xdf)) points(xdf[[i]], ydf[[i]], col = adjustcolor(i, .6), pch = i, cex = .9, lwd = 1.5)
for(i in seq_along(xdf)) points(mean(xdf[[i]]), mean(ydf[[i]]), col = i, pch = i, cex = 3, lwd = 3)
## points(rep(.9*min(xdf$opt), length(xdf$opt)), oldcomps, col = adjustcolor("black", .5), pch = 4)
## points(.9*min(xdf$opt), mean(oldcomps), col = "black", pch = 4, lwd = 3, cex = 2)
## placelabel(
##     paste("Old metric:", round(mean(ydf$opt)/mean(get_error(B$fixed)), 2),
##           "\n(smaller is better, if =1 then no difference)"),
##     .02, .98, adj = c(0, 1)
## )
## cols <- palette.colors(palette = "Okabe-Ito")[-1]
## for(i in seq_along(vdf)) points(demo_error_A[[i]], demo_error_B[[i]], col = cols[i], pch = 0, cex = 2, lwd = 3)
## legend("bottom", bty = "n", col = cols[seq_along(vdf)], pch = 0, pt.lwd = 2,
##        legend = c("Degree", "Betweenness", "Closeness", "Eigenvector"))

dev.off()

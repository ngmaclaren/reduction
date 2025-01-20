## This analysis focuses on the equilibrium state of each selected sentinel node.

library(igraph)
library(optNS)
library(sfsmisc)

save_plots <- FALSE # TRUE

get_zs <- function(vs, fs, dispersion = c("mad", "sd"), center = c("mean", "median")) { 
    dispersion <- get(match.arg(dispersion))
    center <- get(match.arg(center))

    m <- center(fs)
    s <- dispersion(fs)
    if(s == 0) {
        return(rep(NA, length(vs)))
    } else {
        (fs[vs] - m)/s
    }
}

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
## nets <- nets[-which(nets == "er")]
dynamics <- c(
    "doublewell", "mutualistic", "SIS", "genereg"
)

net <- "proximity" # dolphin
dispersion <- "mad" # sd
center <- "median" # mean
## dyn <- "doublewell" # doublewell mutualistic SIS genereg

fss <- readRDS(paste0("../data/fullstate-", net, ".rds"))

## for each sentinel node set I need to find the z scores separately.
## So, get the z scores for each sentinel node set separately. That's a lapply() call, down ns.
## Then, plot the collection of lines

analyze <- function(net, dyn, dispersion, center) {
    fs <- fss[[dyn]]
    nss <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))$opt

    zvals <- lapply(nss, function(ns) {
        vs <- ns$vs
        t(apply(fs, 1, function(xi) get_zs(vs, xi, dispersion, center)))
    })

    ## return(zvals)
    Ds <- switch(dyn, seq(0, 1, length.out = 100), mutualistic = seq(0, 3, length.out = 100))

    ## plot(Ds, rowMeans(fs), type = "l", lty = 1, lwd = 4, col = 1,
    matplot(Ds, fs, type = "l", lty = 1, lwd = 0.5, col = adjustcolor(8, 0.4),
            xlab = "D", ylab = "", font.lab = 3, cex.lab = 2, axes = FALSE)
    for(ns in nss) {
        matlines(Ds, fs[, ns$vs], lty = 1, lwd = 0.5, col = "#3584e4")
    }
    box()
    eaxis(1, cex.axis = 1.5)

    par(new = TRUE)
    ylim <- range(unlist(zvals), na.rm = TRUE)
    plot(NULL, xlim = range(Ds), ylim = ylim,
         xlab = "", ylab = "z", font.lab = 3, cex.lab = 2, axes = FALSE)
    for(i in seq_along(zvals)) {
        matlines(Ds, zvals[[i]], lty = 1, lwd = 0.5, col = adjustcolor(2, 0.4))
    }
    eaxis(2, cex.axis = 1.5)

    par(new = TRUE)
    disps <- apply(fs, 1, get(dispersion))
    dispcolor <- 3
    plot(Ds, disps, type = "l", lty = 1, lwd = 2, col = dispcolor, axes = FALSE, xlab = "", ylab = "")
    eaxis(4, cex.axis = 1.5, col.axis = dispcolor, col = dispcolor, small.args = list(col = dispcolor))
    ## mtext("Measure of dispersion", side = 4, line = 3, cex = 2, col = dispcolor)
}

ht <- 10
wd <- 10
if(save_plots) {
    pdf(paste0("../img/state-outliers-", dispersion, "-", net, ".pdf"), height = ht, width = wd)
} else {
    dev.new(width = wd, height = ht)
}
par(mfrow = c(2, 2), mar = c(5, 5, 1, 4), pty = "s")
for(i in seq_along(dynamics)) {
    analyze(net, dynamics[i], dispersion, center)
    if(i == 1) {
        legend("bottomright", bty = "n",
               legend = c(
                   "Node state, x scale",
                   "Sentinel node state, x scale",
                   "Sentinel node state, z scale",
                   "Measure of dispersion"
               ),
               col = c(8, "#3584e4", 2, 3), lty = 1, lwd = 1.5)
    }
    mtext(dynamics[i], line = -2, cex = 1.5)
}
if(save_plots) dev.off()

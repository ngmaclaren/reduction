library(latex2exp)
library(sfsmisc)
get_error <- optNS::get_error

dynamics <- "doublewell"

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)

allns <- lapply(networks, function(network) readRDS(paste0("../data/ns-", network, "_", dynamics, ".rds")))
allerrors <- lapply(allns, function(ns) lapply(ns, get_error))

palette("Set 1")
## palette("Tableau 10")

imgfile <- paste0("../img/compare-networks-", dynamics, ".pdf")

placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}

plotit <- function(dl, panellabel = "", xlab = "", ylab = "") { # list of errors for a given network
    re <- dl$rand; fe <- dl$fixed; oe <- dl$opt
    ht <- 2
    wd <- 7
    amnt <- 0.05
    alpha <- 1
    ptcex <- 0.75
    ylim <- c(.5, 1.5)
    xlim <- range(c(re, fe, oe))
    ypos <- 1 + c(.15, 0, -.15)
    plot(
        NULL,
        xlim = xlim, ylim = ylim, log = "x",
        axes = FALSE, xlab = xlab, ylab = ylab, las = 1
    )
    points(
        oe, jitter(rep(ypos[1], length(oe)), amount = amnt),
        col = adjustcolor(1, alpha.f = alpha), pch = 19, cex = ptcex
    )
    points(
        fe, jitter(rep(ypos[2], length(fe)), amount = amnt),
        col = adjustcolor(2, alpha.f = alpha), pch = 19, cex = ptcex
    )
    points(
        re, jitter(rep(ypos[3], length(re)), amount = amnt),
        col = adjustcolor(3, alpha.f = alpha), pch = 19, cex = ptcex
    )
}

labelsize <- 1.75
ht <- 10
wd <- 10
dev.new(height = ht, width = wd)
par(mfcol = c(length(networks)/2, 2), mar = c(2, 2, 0.25, 2)+0.75)
for(i in seq_along(networks)) {
    plotit(allerrors[[i]])
    eaxis(1, axTicks(1)[c(TRUE, FALSE)], cex.axis = labelsize, at.small = FALSE, sub10 = "10")
    placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
    if(i == 1) {
        legend(
            "topright", bty = "n", col = 1:3, pch = 19, pt.cex = 1.5, cex = 1.5,
            legend = c("Optimized", "Degree-preserving random", "Completely random")
        )
    }
}

library(latex2exp)
library(sfsmisc)
get_error <- optNS::get_error

save_plots <- FALSE # TRUE
two_panels <- TRUE # FALSE
optweights <- "no" # "yes"
use_weighted_networks <- FALSE # TRUE
use_directed_networks <- FALSE # TRUE

dynamics <- "doublewell" # genereg SIS genereg mutualistic

networks <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)
weighted <- c(
    "windsurfers", "macaques", "train_terrorists", "highschool", "drug", "residence_hall", "netsci_weighted",
    "proximity_weighted", "gap_junction_herm", "intl_trade"
)
directed <- c(
    "canton", "physician_trust", "email_company", "flamingo", "ecoli", "yeast", "usair", "jung-c", "email_uni",
    "faa"
)
if(use_weighted_networks) networks <- weighted
if(use_directed_networks) networks <- directed

allns <- lapply(
    networks,
    function(network) {
        readRDS(
            paste0("../data/ns-", network, "_", dynamics,
                   switch(optweights, no = "", yes = "_w"),
                   ".rds")
        )
    }
)
allerrors <- lapply(allns, function(ns) lapply(ns, get_error))

## palette("Set 1")
## palette("Tableau 10")

if(use_weighted_networks) {
    imgfile <- paste0("../img/compare-weighted-networks-", dynamics, ".pdf")
} else if(use_directed_networks) {
    imgfile <- paste0("../img/compare-directed-networks-", dynamics, ".pdf")
} else {
    imgfile <- paste0("../img/compare-networks-", dynamics,
                      switch(optweights, no = "", yes = "-w"),
                      ".pdf")
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

plotit <- function(dl, panellabel = "", xlab = "", ylab = "") { # list of errors for a given network
    re <- dl$rand; fe <- dl$fixed; oe <- dl$opt
    ht <- 2
    wd <- 7
    amnt <- 0.05
    alpha <- 0.85
    ptcex <- 1
    ylim <- c(.5, 1.5)
    xlim <- range(c(re, fe, oe))
    ypos <- 1 + c(.25, 0, -.25)
    plot(
        NULL,
        xlim = xlim, ylim = ylim, log = "x",
        axes = FALSE, xlab = xlab, ylab = ylab, las = 1
    )
    points(
        oe, jitter(rep(ypos[1], length(oe)), amount = amnt),
        col = adjustcolor("#3584e4", alpha), #adjustcolor(1, alpha.f = alpha),
        pch = 1, cex = ptcex
    )
    points(
        fe, jitter(rep(ypos[2], length(fe)), amount = amnt),
        col = adjustcolor("#ff7800", alpha), #adjustcolor(2, alpha.f = alpha),
        pch = 2, cex = ptcex
    )
    points(
        re, jitter(rep(ypos[3], length(re)), amount = amnt),
        col = adjustcolor("#33d17a", alpha), #adjustcolor(3, alpha.f = alpha),
        pch = 0, cex = ptcex
    )
}

if(two_panels) {
    labelsize <- 1.75
    ht <- 6
    wd <- 8

    if(save_plots) {
        pdf("../img/compare-networks-MT.pdf", height = ht, width = wd)
    } else {
        dev.new(height = ht, width = wd)
    }

    par(mfcol = c(2, 1), mar = c(4, 16, 0.5, 0.5))
    nets <- c("dolphin", "ba")
    yticks <- list(
        labels = c("Optimized", "Degree-preserving random", "Completely random"),
        pos = 1 + c(.25, 0, -.25)
    )

    for(i in seq_along(nets)) {
        plotit(allerrors[[i]])
        box()
        eaxis(1, n.axp = 1, cex.axis = labelsize)
        title(xlab = "Approximation error", cex.lab = labelsize)
        axis(2, at = yticks$pos, labels = yticks$labels, tick = TRUE, las = 2, cex.axis = 0.75*labelsize)
        placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
        ## if(i == 1) {
        ##     legend(
        ##         "topright", bty = "n", col = c("#3584e4", "#ff7800", "#33d17a"), #col = 1:3,
        ##         pch = c(1, 2, 0), pt.cex = 1.5, cex = 1.5, pt.lwd = 2,
        ##         legend = c("Optimized", "Degree-preserving random", "Completely random")
        ##     )
        ## }
    }

    if(save_plots) dev.off()
} else {
    labelsize <- 1.75
    ht <- 10
    wd <- 12

    if(save_plots) {
        pdf(imgfile, height = ht, width = wd)
    } else {
        dev.new(height = ht, width = wd)
    }

    par(mfcol = c(length(networks)/2, 2), mar = c(2, .5, .5, .5))
    for(i in seq_along(networks)) {
        plotit(allerrors[[i]])
        ##eaxis(1, axTicks(1)[c(TRUE, FALSE)], cex.axis = labelsize, at.small = FALSE, sub10 = "10")
        box()
        eaxis(1, n.axp = 1, cex.axis = labelsize)
        placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
        if(i == 1) {
            legend(
                "topright", bty = "n", col = c("#3584e4", "#ff7800", "#33d17a"), #col = 1:3,
                pch = c(1, 2, 0), pt.cex = 1.5, cex = 1.5, pt.lwd = 2,
                legend = c("Optimized", "Degree-preserving random", "Completely random")
            )
        }
    }

    if(save_plots) dev.off()
}

## plotit <- function(dl) {
##     dens <- lapply(dl[1:3], function(x) density(x, bw = "SJ"))
##     idxs <- lapply(dens, function(d) which(d$x > 0 & d$y > 10^(-4)))

##     plotx <- mapply(function(d, i) d$x[i], dens, idxs, SIMPLIFY = FALSE)
##     ploty <- mapply(function(d, i) d$y[i], dens, idxs, SIMPLIFY = FALSE)
    
##     xlim <- range(unlist(plotx))
##     ylim <- range(unlist(ploty))

##     plot(NULL, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "", ylab = "")#, log = "xy")
##     lines(plotx$opt, ploty$opt, col = "#3584e4", lwd = 3)
##     lines(plotx$fixed, ploty$fixed, col = "#ff7800", lwd = 3)
##     lines(plotx$rand, ploty$rand, col = "#33d17a", lwd = 3)
## }

## plotit <- function(dl) {
##     hists <- lapply(dl[1:3], function(x) hist(x, plot = FALSE))

##     xlim <- range(unlist(lapply(hists, `[[`, "breaks")))
##     ylim <- c(0, max(unlist(lapply(hists, `[[`, "counts"))))

##     plot(hists$rand, col = "#33d17a", xlim = xlim, ylim = ylim)
##     plot(hists$fixed, col = "#ff7800", add = TRUE)
##     plot(hists$opt, col = "#3584e4", add = TRUE)
## y}

library(sfsmisc)
get_error <- optNS::get_error

save_plots <- TRUE # FALSE
two_panels <- TRUE # FALSE
large_and_model <- FALSE # TRUE
                                        # The next three options are broken now b/c of the new file names, l.51+
optweights <- "no" # "yes"
use_weighted_networks <- FALSE # TRUE
use_directed_networks <- FALSE # TRUE

dynamics <- "doublewell" # doublewell mutualistic SIS genereg

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
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
if(large_and_model) networks <- networks[11:20] else networks <- networks[1:10]
if(two_panels) networks <- c("dolphin", "ba")

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

if(use_weighted_networks) {
    imgfile <- paste0("../img/compare-weighted-networks-", dynamics, ".pdf")
} else if(use_directed_networks) {
    imgfile <- paste0("../img/compare-directed-networks-", dynamics, ".pdf")
} else if(large_and_model) {
    imgfile <- paste0("../img/compare-large_and_model-networks-", dynamics, ".pdf")
} else {
    imgfile <- paste0("../img/compare-small_empirical-networks-", dynamics,
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
        col = adjustcolor("#3584e4", alpha), 
        pch = 1, cex = ptcex
    )
    points(
        fe, jitter(rep(ypos[2], length(fe)), amount = amnt),
        col = adjustcolor("#ff7800", alpha), 
        pch = 2, cex = ptcex
    )
    points(
        re, jitter(rep(ypos[3], length(re)), amount = amnt),
        col = adjustcolor("#33d17a", alpha), 
        pch = 0, cex = ptcex
    )
}

if(two_panels) {
    labelsize <- 1.75
    ht <- 6
    wd <- 8

    if(save_plots) {
        pdf("../img/compare-networks-MT-v2.pdf", height = ht, width = wd)
    } else {
        dev.new(height = ht, width = wd)
    }

    par(mfcol = c(2, 1), mar = c(4, 16, 0.5, 0.5))
    ## nets <- c("dolphin", "ba")
    yticks <- list(
        labels = c("Optimized", "Degree-preserving random", "Completely random"),
        pos = 1 + c(.25, 0, -.25)
    )

    for(i in seq_along(networks)) {
        plotit(allerrors[[i]])
        box()
        eaxis(1, n.axp = 1, cex.axis = labelsize)
        title(xlab = "Approximation error", cex.lab = labelsize)
        axis(2, at = yticks$pos, labels = yticks$labels, tick = TRUE, las = 2, cex.axis = 0.75*labelsize)
        placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
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

    if(large_and_model) IDadj <- 10 else IDadj <- 0

    par(mfcol = c(length(networks)/2, 2), mar = c(2, 1.25, .5, 1.25))
    for(i in seq_along(networks)) {
        plotit(allerrors[[i]])
        box()
        eaxis(1, n.axp = 1, cex.axis = labelsize)
        placelabel(paste0("(", letters[i + IDadj], ")"), 0.01, 0.95, adj = c(0, 1), cex = labelsize)
        if(i == 1) {
            legend(
                "topright", bty = "n", col = c("#3584e4", "#ff7800", "#33d17a"), 
                pch = c(1, 2, 0), pt.cex = 1.5, cex = 1.5, pt.lwd = 2,
                legend = c("Optimized", "Degree-preserving random", "Completely random")
            )
        }
    }

    if(save_plots) dev.off()
}


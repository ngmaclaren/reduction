library(sfsmisc)
library(optNS)

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
networklabels <- c(
    "Dolphin", "Proximity", "Metabolic", "Road", "Email",
    "FlyBi", "Reactome", "Route views", "Spanish words", "FOLDOC",
    "Tree of life", "English words", "Enron", "Marker Cafe", "Propser",
    "ER", "BA", "HK", "GKK", "LFR"
)
dynamics <- c(
    "doublewell", "mutualistic", "SIS", "genereg"
)
dynamicslabels <- c(
    "double-well", "mutualistic", "SIS", "gene"
)

## make into function here
plotit <- function(net, dyn) {
    ## net <- "dolphin"
    ## dyn <- "doublewell"
                                        # optimized node set
    ns <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))$opt
    best <- ns[[which.min(get_error(ns))]]
                                        # full state
    Y <- readRDS(paste0("../data/fullstate-", net, ".rds"))[[dyn]]
                                        # the black line
    y <- rowMeans(Y)
                                        # the purple line
    z <- rowMeans(Y[, best$vs])
                                        # x axis
    Ds <- switch(dyn, seq(0, 1, length.out = 100), mutualistic = seq(0, 3, length.out = 100))
                                        # plot
    plot(NULL, xlab = "", ylab = "", xlim = range(Ds), ylim = range(c(z, y)), axes = FALSE, font.lab = 3)
    lines(Ds, y, lwd = 2, col = "black")
    lines(Ds, z, lwd = 2, col = "#9141ac")
    box()
    eaxis(1)
    eaxis(2)
    mtext(networklabels[which(networks == net)], side = 1, line = -1, adj = 0.98)
    mtext(dynamicslabels[which(dynamics == dyn)], side = 1, line = -2.5, adj = 0.98)
}

save_plots <- TRUE # FALSE
wd <- 6.5
ht <- (wd/4)*5

if(save_plots) {
    pdf("../img/track-critical-states.pdf", height = ht, width = wd)
} else {
    dev.new(width = wd, height = ht)
}
par(mfrow = c(5, 4), mai = c(0.25, 0.25, 0.1, 0.1), pty = "s")
for(i in seq_along(networks)) {
## for(i in 1:5) {
    net <- networks[i]
    for(j in seq_along(dynamics)) {
        dyn <- dynamics[j]
        plotit(net, dyn)
    }
}
if(save_plots) dev.off()

library(igraph)

save_plots <- TRUE # FALSE

networks <- c("dolphin", "proximity", "celegans", "euroroad", "email", "lfr")
graphlist <- lapply(
    paste0("../data/", networks, ".rds"),
    readRDS
)
names(graphlist) <- networks

parts <- lapply(graphlist, cluster_louvain)
lengths(parts)
mbrs <- lapply(parts, membership)

max_fromsame <- function(ns, part) {
    mbr <- part[ns$vs]
    max(sapply(mbr, function(x) sum(mbr == x)))
}

opts <- lapply(networks, function(network) readRDS(paste0("../data/ns-", network, "_doublewell.rds"))$opt)
rands <- lapply(networks, function(network) readRDS(paste0("../data/ns-", network, "_doublewell.rds"))$rand)
names(opts) <- names(rands) <- networks

maxFromSame.opt <- mapply(function(opt, mbr) sapply(opt, max_fromsame, mbr), opts, mbrs, SIMPLIFY = FALSE)
maxFromSame.rand <- mapply(function(rand, mbr) sapply(rand, max_fromsame, mbr), rands, mbrs, SIMPLIFY = FALSE)

ht <- 8
wd <- 5
labelsize <- 1.5
## palette("Dark 2")
if(save_plots) {
    pdf("../img/community-histograms.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(length(opts), 1), mar = c(4.5, 4.5, 1.5, 0.5))
for(i in seq_along(opts)) {
    nbins <- length(opts[[i]][[1]]$vs)
    ylim <- c(0, max(c(table(maxFromSame.rand[[i]]), table(maxFromSame.opt[[i]]))))
    barplot(
        matrix(
            c(tabulate(maxFromSame.opt[[i]], nbins = nbins),
              tabulate(maxFromSame.rand[[i]], nbins = nbins)),
            byrow = TRUE, nrow = 2),
        col = c("#3584e4", "#33d17a"), beside = TRUE,
        ylim = ylim, cex.axis = labelsize,
        names.arg = seq(nbins), cex.names = labelsize, cex.lab = labelsize,
        ylab = "Frequency",
        xlab = switch(i, "", "", "", "", "Max. nodes from same community")
    )
    mtext(paste0("(", letters[i], ")"), line = 0.1, adj = 0.01, cex = labelsize)
    if(i == 1) legend(
                   "topright", bty = "n", pch = 15, col = c("#3584e4", "#33d17a"),
                   legend = c("Optimized", "Random"),
                   cex = 0.75*labelsize, pt.cex = 2
               )
}
if(save_plots) dev.off()

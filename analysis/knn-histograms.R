library(igraph)

save_plots <- TRUE # FALSE

networks <- c("dolphin", "proximity", "celegans", "euroroad", "email")
graphlist <- lapply(
    paste0("../data/", networks, ".rds"),
    readRDS
)
names(graphlist) <- networks

knns <- lapply(graphlist, function(g) knn(g)$knn)

opts <- lapply(networks, function(network) readRDS(paste0("../data/ns-", network, "_doublewell.rds"))$opt)
rands <- lapply(networks, function(network) readRDS(paste0("../data/ns-", network, "_doublewell.rds"))$rand)

ht <- 8
wd <- 5
labelsize <- 1.5
palette("Dark 2")
if(save_plots) {
    pdf("../img/knn-histograms.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(length(opts), 1), mar = c(4.5, 4.5, 1.5, 0.5))
for(i in seq_along(opts)) {
    pd.opt <- knns[[i]][as.numeric(optNS::get_vs(opts[[i]]))]
    ## pd.opt <- knns[[i]][colMeans(optNS::get_vs(opts[[i]]))]
    pd.rand <- knns[[i]][as.numeric(optNS::get_vs(rands[[i]]))]
    ## pd.rand <- knns[[i]][colMeans(optNS::get_vs(rands[[i]]))]
    hist.all <- hist(c(pd.opt, pd.rand), breaks = 20, plot = FALSE)
    hist.opt <- hist(pd.opt, breaks = hist.all$breaks, plot = FALSE)
    hist.rand <- hist(pd.rand, breaks = hist.all$breaks, plot = FALSE)
    xlim <- range(c(hist.opt$breaks, hist.rand$breaks))
    ylim <- range(c(hist.opt$counts, hist.rand$counts))
    plot(hist.rand, col = adjustcolor(3, 0.75), border = 3, density = 10, angle = 45,
         xlim = xlim, ylim = ylim, main = "", 
         xlab = switch(i, "", "", "", "", "Node average nearest neighbor degree"),
         cex.axis = labelsize, cex.lab = labelsize
         )
    plot(hist.opt, col = adjustcolor(1, 0.75), border = 1, density = 10, angle = 135,
         xlim = xlim, ylim = ylim, axes = FALSE, add = TRUE, main = "")
    abline(v = mean(pd.opt), lwd = 3, col = 1)
    abline(v = mean(pd.rand), lwd = 3, col = 3)
    if(i == 1) legend("topleft", bty = "n", col = c(1, 3), lwd = 2, legend = c("Optimized", "Random"))
    mtext(paste0("(", letters[i], ")"), cex = labelsize, adj = 0.01, line = 0.1)
}
if(save_plots) dev.off()

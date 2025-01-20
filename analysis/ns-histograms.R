library(igraph)
library(sfsmisc)
library(optNS)

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
nodesets <- lapply(nets, function(network) {
    readRDS(paste0("../data/ns-", network, "_doublewell.rds"))[c("opt", "rand")]
})
networks <- lapply(nets, function(network) {
    readRDS(paste0("../data/", network, ".rds"))
})
names(networks) <- names(nodesets) <- nets

res <- list(rand = numeric(), opt = numeric(), mark = numeric())
dev.new(height = 12, width = 15)
par(mfrow = c(4, 5), mar = c(3, 3, 0, 0))
for(i in seq_along(nets)) {
    randdata <- as.numeric(get_ks(nodesets[[i]]$rand))
    optdata <- as.numeric(get_ks(nodesets[[i]]$opt))
    rand <- hist(randdata, breaks = 20, plot = FALSE)
    opt <- hist(optdata, breaks = rand$breaks, plot = FALSE)
    g <- networks[[i]]
    k <- degree(g)
    mark <- quantile(k, 0.99) # 0.99 0.95
    ## print(nets[i])
    ## print("rand")
    ## print(sum(randdata > mark)/length(randdata))
    res$rand[i] <- sum(randdata > mark)/length(randdata)
    ## print("opt")
    ## print(sum(optdata > mark)/length(optdata))
    res$opt[i] <- sum(optdata > mark)/length(optdata)
    res$mark[i] <- mark
    ylim <- range(c(rand$counts, opt$counts))
    plot(rand, main = "", xlab = "Degree", col = adjustcolor("gray50", .5), ylim = ylim)
    plot(opt, add = TRUE, col = adjustcolor(2, .5))
    if(i == 1) legend("topright", bty = "n", pch = 15, col = c(adjustcolor(2, .5), adjustcolor("gray50", .5)),
                      legend = c("Optimized", "Random"), pt.cex = 2)
}

dev.new()
par(mar = c(5, 6, 1, 1), pty = "s")
plot(res$mark, # the 95th percentile of k
     res$rand, # the prop of rand nodes larger than mark
     log = "x",
     pch = 0, col = 1,
     xlab = "95th % of degree distribution",
     ylab = "",
     cex.lab = 2, cex = 2,
     ylim = range(c(res$rand, res$opt)), axes = FALSE)
points(res$mark, res$opt, pch = 1, cex = 2, col = 1)
arrows(res$mark, y0 = res$rand, y1 = res$opt, lty = 1, lwd = 1, col = 1, length = 0.1)
box()
eaxis(1, cex.axis = 2, n.axp = 1)
eaxis(2, cex.axis = 2)
title(ylab = "Prop. selected nodes larger than 95th %", cex.lab = 2, line = 4.5)
legend("bottomright", bty = "n", pch = c(0, 1), pt.cex = 1.5, cex = 1.5, col = 1,
       legend = c("Random", "Optimized"))

sum(res$rand > res$opt)/length(res$rand)

## look at the actual fraction of nodes in each type of node set that belongs to the top 1% or 5%

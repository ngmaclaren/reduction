library(igraph)
library(optNS)
library(sfsmisc)

save_plots <- TRUE # FALSE

nets <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)
nets <- nets[-which(nets == "er")]
dynamics <- c(
    "doublewell", "mutualistic", "SIS", "genereg"
)
networks <- lapply(nets, function(network) {
    readRDS(paste0("../data/", network, ".rds"))
})
names(networks) <- nets

countit <- function(net, dyn, threshold = 0.99) {
    nss <- readRDS(paste0("../data/ns-", net, "_", dyn, ".rds"))#[c("opt", "rand")]
    opt <- list(ns = nss$opt)
    rand <- list(ns = nss$rand)
    g <- networks[[net]]
    N <- vcount(g)
    k <- degree(g)
    kbar <- mean(k)
    mark <- quantile(k, threshold)

    opt$ks <- as.numeric(get_ks(opt$ns))
    rand$ks <- as.numeric(get_ks(rand$ns))

    opt$count <- sum(opt$ks >= mark) # how many selected nodes were above the mark
    rand$count <- sum(rand$ks >= mark)

    opt$prop <- opt$count/length(opt$ks) # what proportion was that count of all selected nodes
    rand$prop <- rand$count/length(rand$ks)

    return(data.frame(
        network = net, kbar = kbar, dynamics = dyn, threshold = threshold,
        opt.count = opt$count, rand.count = rand$count,
        opt.prop = opt$prop, rand.prop = rand$prop
    ))
}

## countit("dolphin", "doublewell", threshold = 0.99)
## countit("dolphin", "doublewell", threshold = 0.95)

conds <- expand.grid(nets, dynamics, stringsAsFactors = FALSE)
df99 <- do.call(rbind, apply(conds, 1, function(row) countit(row[1], row[2], threshold = 0.99), simplify = FALSE))
df95 <- do.call(rbind, apply(conds, 1, function(row) countit(row[1], row[2], threshold = 0.95), simplify = FALSE))
df <- rbind(df99, df95)
df$dynamics <- factor(df$dynamics, levels = dynamics)

round(sum(df99$rand.count > df99$opt.count)/nrow(df99), 4) # 0.9473684
round(sum(df95$rand.count > df95$opt.count)/nrow(df95), 4) # 0.6710526

plotit <- function(bigdf, threshold = 0.95, optcolor = "#3584e4", randcolor = "#33d17a",
                   show.legend = TRUE, show.title = FALSE, show.smooth = FALSE, show.guides = TRUE) {
    df <- bigdf[bigdf$threshold == threshold, ]
    ylim <- range(c(df$opt.prop, df$rand.prop))
    ylim[1] <- ylim[1] - 0.02*ylim[2]
    ylim[2] <- ylim[2] + 0.02*ylim[2]
    par(mar = c(5, 6, 1, 1), pty = "s")
    plot(
        NULL, log = "x", yaxs = "i",
        xlim = range(df$kbar), ylim = ylim,
        xlab = "Average degree", ylab = "", axes = FALSE, cex.lab = 2
    )
    points(df$kbar, df$rand.prop, pch = as.numeric(df$dynamics), col = randcolor, cex = 1.5, lwd = 1.5)
    points(df$kbar, df$opt.prop, pch = as.numeric(df$dynamics), col = optcolor, cex = 1.5, lwd = 1.5)
    if(show.guides) abline(h = 1 - threshold, lwd = 2, col = "#c01c28", lty = 2)
    if(show.smooth) {
        lines(lowess(df$kbar, df$rand.prop), col = randcolor, lwd = 2)
        lines(lowess(df$kbar, df$opt.prop), col = optcolor, lwd = 2)
    }
    box()
    eaxis(1, cex.axis = 1.75, n.axp = 1)
    eaxis(2, cex.axis = 1.75)
    title(ylab = "Fraction of nodes", cex.lab = 2, line = 4.5)
    if(show.title) title(main = paste0("Threshold: ", df$threshold[1]))
    if(show.legend) {
        legend("topright", bty = "n", legend = c("Optimized", "Random"),
               col = c(optcolor, randcolor), lwd = 2, cex = 1.75, pt.cex = 1.75)
        legend("topleft", bty = "n", legend = dynamics, pch = 1:4, col = "gray50", cex = 1.75, pt.cex = 1.75)
    }
}    

ht <- 7
wd <- 14
if(save_plots) {
    pdf("../img/ns-count-extreme.pdf", height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(1, 2))
plotit(df, 0.99, show.legend = TRUE, show.smooth = FALSE)
plotit(df, 0.95, show.legend = FALSE, show.smooth = FALSE)
if(save_plots) dev.off()

## File doesn't take long, but loading the node features data frame is the slowest part.
## Cycle through conditions at the bottom to save many files.

library(sfsmisc)
library(optNS)

networks <- c(
    "dolphin", "proximity", "celegans", "euroroad", "email",
    "drosophila", "reactome", "route_views", "spanish", "foldoc",
    "tree_of_life", "word_assoc", "enron", "marker_cafe", "prosper",
    "er", "ba", "hk", "gkk", "lfr"
)

nodesets <- lapply(networks, function(network) {
    readRDS(paste0("../data/ns-", network, "_doublewell.rds"))[c("opt", "rand")]
})
nodefeatures <- lapply(networks, function(network) {
    subset(read.csv(paste0("../data/nodefeatures-", network, ".csv")), dynamics == "doublewell")
})
names(nodefeatures) <- names(nodesets) <- networks

## For each network,
### For each variable,
#### Every time a node appears in a node set, record that node's variable value.
#### So, if v_i has k_i=5 and v_i appears four times, it should contribute 5 four times to the k_i vector.

plotit_surv <- function(net, var, log = "",
                   optcolor = "#3584e4", randcolor = "#33d17a", graycolor = "#7f7f7f", labelsize = 2,
                   show.legend = "no", show.xlabel = "no", show.ylabel = "no") {
    nf <- nodefeatures[[net]]
    ns <- nodesets[[net]]

    opt <- list()
    rand <- list()
    legendtext <- c("Optimized", "Random")
    
    opt$vs <- as.numeric(get_vs(ns$opt))
    rand$vs <- as.numeric(get_vs(ns$rand))

    opt$var <- sort(nf[match(opt$vs, nf[, "v"]), var])
    rand$var <- sort(nf[match(rand$vs, nf[, "v"]), var])
    ## if(var == "lcl") xlim <- c(0, 1) else xlim <- range(c(opt$var, rand$var))
    xlim <- range(c(opt$var, rand$var))
    if(grepl("x", log)) xlim <- range(c(opt$var[opt$var > 0], rand$var[rand$var > 0]))
    if(var == "lcl" & net == "er") xlim <- c(0, max(c(opt$var, rand$var)))
    if(var == "lcl" & net == "prosper") xlim <- c(0, 0.04)
    if(var == "kcore" & net %in% c("ba", "hk")) xlim <- c(1, 3)
    
    opt$frac <- sapply(opt$var, function(x) sum(opt$var > x)/length(opt$var))
    rand$frac <- sapply(rand$var, function(x) sum(rand$var > x)/length(rand$var))
    opt$cprob <- sapply(opt$var, function(x) sum(opt$var <= x)/length(opt$var))
    rand$cprob <- sapply(rand$var, function(x) sum(rand$var <= x)/length(rand$var))
    
    all_y <- c(opt$frac, rand$frac, opt$cprob, rand$cprob)
    ylim <- range(c(all_y, 1))

    ## if(grepl("y", log)) ylim <- with(list(y = c(opt$frac, rand$frac)), range(y[y > 0]))
    if(grepl("y", log)) ylim <- range(c(all_y[all_y > 0], 1))
    ## if(ylim[1] > 0.1 & ylim[2] < 1.0) ylim <- c(0.1, 1)

    xlab <- switch(show.xlabel, no = "", yes = var)
    ylab <- switch(show.ylabel, no = "", yes = expression(p(k) > k))

    if(grepl("y", log) & ylim[1] > 10^(-1.1)) ylim[1] <- 10^(-1.1)
    
    plot(NULL, xlim = xlim, ylim = ylim, log = log, axes = FALSE, xaxs = "i", yaxs = "i",
         main = "", ylab = ylab, xlab = xlab, cex.axis = labelsize, cex.lab = labelsize)
    ## print(net)
    ## print(par("usr"))
    ## print(c(xlim, ylim))
    lines(rand$var, rand$frac, type = "s", col = randcolor, lwd = 2)
    lines(opt$var, opt$frac, type = "s", col = optcolor, lwd = 2)
    if(var %in% c("k", "bc", "knn")) { # , "cc"
        lines(rand$var, rand$cprob, type = "s", col = randcolor, lwd = 2, lty = 2)
        lines(opt$var, opt$cprob, type = "s", col = optcolor, lwd = 2, lty = 2)
        legendtext <- c(legendtext, "Cumulative prob.")
    }
    eaxis(1, cex.axis = labelsize, n.axp = 1)
    eaxis(2, cex.axis = labelsize, n.axp = 1)
    box()
    if(show.legend == "yes") {
        legend(
            "bottomleft", cex = .75*labelsize, lwd = 2, lty = c(1, 1, 2), bty = "n",
            col = c(optcolor, randcolor, graycolor),
            legend = legendtext
        )
    }
    if(var == "kcore" & net %in% c("ba", "hk")) {
        lines(seq(1, 3, length.out = 1000), c(rep(1, 500), rep(0, 500)), type = "s",
              col = randcolor, lwd = 2, lty = 1)
        lines(seq(1, 3, length.out = 1000), c(rep(1, 500), rep(0, 500)), type = "s",
              col = optcolor, lwd = 2, lty = 1)
        ## abline(v = mean(rand$var), col = randcolor)
        ## abline(v = mean(opt$var), col = optcolor)
    }
}

save_plots <- TRUE # FALSE
var <- "kcore" # k knn lcl cc bc kcore
log <- switch(var, "", k = "xy", knn = "xy", lcl = "y", bc = "xy")
ht <- 2*(4*6.5/5)
wd <- 2*(6.5)
if(save_plots) {
    pdf(paste0("../img/nodefeatures-", var, ".pdf"), height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mfrow = c(4, 5), mai = c(.3, .5, .125, .2))
for(i in seq_along(networks)) {
    net <- networks[i]
    if(i == 1) show.legend <- "yes" else show.legend <- "no"
    ## plotit(net, var, show.legend = show.legend)
    plotit_surv(net, var, log = log, show.legend = show.legend)
    mtext(paste0("(", letters[i], ")"), line = -3, adj = 0.025, cex = 2)
}
if(save_plots) dev.off()

## plotit <- function(net, var, 
##                    optcolor = "#3584e4", randcolor = "#33d17a", labelsize = 2,
##                    show.legend = "no", show.xlabel = "no", show.ylabel = "no") {
##     nf <- nodefeatures[[net]]
##     ns <- nodesets[[net]]
##     opt_vs <- as.numeric(get_vs(ns$opt))
##     rand_vs <- as.numeric(get_vs(ns$rand))

##     opt <- nf[match(opt_vs, nf[, "v"]), var]
##     rand <- nf[match(rand_vs, nf[, "v"]), var]

##     hist_opt_rand <- hist(c(opt, rand), breaks = 20, plot = FALSE)
##     hist_opt <- hist(opt, breaks = hist_opt_rand$breaks, plot = FALSE)
##     hist_rand <- hist(rand, breaks = hist_opt_rand$breaks, plot = FALSE)
##     ylim <- c(0, max(c(hist_opt$counts, hist_rand$counts)))

##     plot(hist_rand, col = adjustcolor(randcolor, .5), ylim = ylim, main = "",
##          cex.axis = labelsize, cex.lab = labelsize)
##     plot(hist_opt, col = adjustcolor(optcolor, .5), add = TRUE)
##     if(show.legend == "yes") {
##         legend(
##             "topright", cex = .75*labelsize, pch = 22, pt.bg = adjustcolor(c(optcolor, randcolor), .5), pt.cex = 2,
##             legend = c("Optimized", "Random"), bty = "n"
##         )
##     }
## }

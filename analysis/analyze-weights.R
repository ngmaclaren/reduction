## I want to see the typical weight assigned to a node of a certain degree.
## If the y-axis goes from zero to one, it will be apparent that balancing of some kind is happening.
## The x-axis should cover the whole range of k_i chosen for 1:maxn.

library(optparse)
optionlist <- list(
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files. Default is %default. "
    ),
    make_option(
        c("-w", "--use-weights"), action = "store_true", default = TRUE,#FALSE,
        help = "Use the simulation results with optimized node weights when generating analysis results. Default is %default. "
    ),
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is 'dw'. Options: 'dw', 'SIS', 'genereg', 'mutualistic', and 'wilsoncowan'."
    ),
    make_option(
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "The network to use in generating analysis results. Ignored if only one option is available. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist), # args = c("-s", "-w", "-g", "dolphin"),
    convert_hyphens_to_underscores = TRUE
)

palette("Tableau 10")

use_weights <- args$use_weights
network <- args$network
dynamics <- args$dynamics
saveplots <- args$save_plots # ensure set locally, rather than by what was used before
infile <- switch(
    use_weights + 1,
    paste0("../data/degree-sequences-", network, "-", dynamics, ".RData"),
    paste0("../data/degree-sequences-", network, "-", dynamics, "-weighted.RData")
)
outfile <- paste0("../img/ks-and-ws-", network, "-", dynamics, ".pdf")

load(infile)
save_plots <- saveplots
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}

labelsize <- 1.75
ticksize <- 1.75

dfs <- lapply(seq(maxn), function(n) {
    df <- data.frame(
        v = unlist(lapply(dl[[n]], `[[`, "vs")),
        k = unlist(lapply(dl[[n]], `[[`, "ks")),
        w = unlist(lapply(dl[[n]], `[[`, "ws")),
        n = n,
        trial = rep(1:ntrials, each = n)
    )
    df$order <- unlist(lapply(split(df, df$trial), function(x) order(x$k)))
    df
})
df <- do.call(rbind, dfs)

                                        # simple test
palette("Tableau 10")#palette(hcl.colors(8))
pdf(paste0("../img/new-ksws-tst-", network, ".pdf"), height = 2*3, width = 4*3)
par(mfrow = c(2, 4))
par(mar = c(4, 4, 0, 0) + 0.5)
for(n in 1:8) {
    plotdf <- df[df$n == n, ]
    plot(
        plotdf$k, plotdf$w, ylim = c(0, 1), xlim = range(k),
        cex = 3, pch = 16, col = adjustcolor(1, alpha.f = .25), #plotdf$order,
        xlab = "k", ylab = "w", font.lab = 3, cex.lab = labelsize, cex.axis = ticksize
    )
    usr <- par("usr")
    par(usr = c(0, 1, 0, 1))
    ## text(0.025, 1 - 0.025, LETTERS[n], font = 2, cex = labelsize)
    placelabel(paste0("(", letters[n], ")\nn = ", n), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
    par(usr = usr)
    ## if(n == 1) legend(
    ##                "topright", pch = 19, col = 1:8,
    ##                legend = c("Smallest node, by degree", rep("", 6), "Largest node, by degree")
    ##            )
}
dev.off()

## pdf(outfile, height = 9, width = 12)
## par(mfrow = c(3, 4), mar = c(4, 4, 0, 0) + .5)

## xlim <- range(unlist(df$k))
## ylim <- c(0, 1)

## for(n in seq(maxn)) {
##     plotdf <- df[df$n == n, ]
##     aggweights <- aggregate(w ~ k, data = plotdf, FUN = mean)
##     aggcounts <- aggregate(w ~ k, data = plotdf, FUN = length)
##     colnames(aggcounts) <- c("k", "c")

##     plot(
##         jitter(plotdf$k, amount = .25), plotdf$w,
##         pch = 19, cex = .5, col = adjustcolor(4, alpha.f = .5),
##         xlim = xlim, ylim = ylim, xlab = "Degree", ylab = "Weight",
##         cex.axis = ticksize, cex.lab = labelsize
##     )
##     points(
##         aggweights$k, aggweights$w, pch = 19, cex = 2*log10(aggcounts$c), col = 3
##     )
##     mtext(paste("n", n, sep = "="), line = -2, adj = 0.98, font = 1, cex = .75*labelsize)
## }

## dev.off()

##                                         # v1: θ is order; v2: θ is k
## ver <- "order"
## n <- 6
## test <- df[df$n == n, ]
## ks <- sort(unique(test$k))
## maxk <- max(test$k)
## test$order <- rep(seq(n), times = ntrials)
## test$x <- switch(ver, k = test$w*cos(2*pi*(test$k/maxk)), order = test$w*cos(2*pi*(test$order/n)))
## ##test$x <- sqrt(test$k)*cos(2*pi*test$w)
## ##test$x <- test$k*cos(2*pi*test$w)
## test$y <- switch(ver, k = test$w*sin(2*pi*(test$k/maxk)), order = test$w*sin(2*pi*(test$order/n)))
## ##test$y <- sqrt(test$k)*sin(2*pi*test$w)
## ##test$y <- test$k*sin(2*pi*test$w)
## palette(hcl.colors(max(test$k)))
## maxr <- max(abs(test[, c("x", "y")]))
## pdf(paste0("../img/ks-and-ws-dolphin-dw-polar-", ver, "-", n, ".pdf"))
## plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(-maxr, maxr), ylim = c(-maxr, maxr))
## lines(test$x, test$y, col = adjustcolor("gray60", alpha.f = .5), lwd = .5, lty = 1)
## points(test$x, test$y, col = adjustcolor(test$k, alpha.f = .5), pch = 19, cex = 2)
## legend("topleft", col = ks, title = "k", pch = 19, pt.cex = 1.5, legend = ks, bty = "n")
## abline(v = 0, col = "black", lty = 1, lwd = .5)
## abline(h = 0, col = "black", lty = 1, lwd = .5)
## dev.off()

## ks <- sort(unique(test$k))
## guides <- lapply(ks, function(k) {
##     list(
##         x = k*cos(2*pi*seq(0, 1, length.out = 1000)),
##         y = k*sin(2*pi*seq(0, 1, length.out = 1000))
##     )
## })

## pdf("../img/ks-and-ws-dolphin-dw-polar.pdf")
## plot(NULL, xlim = range(test$x), ylim = range(test$y), axes = FALSE)
## ##plot(test$x, test$y, pch = 19, col = adjustcolor(test$k, alpha.f = .5), axes = FALSE)
## abline(v = 0, col = "black")
## abline(h = 0, col = "black")
## for(k in ks) lines(guides[[k]]$x, guides[[k]]$y, lty = 1, lwd = .5, col = k)
## points(test$x, test$y, pch = 1, col = test$k)
## dev.off()

## test.r <- reshape(test, timevar = "k", idvar = c("v", "n", "trial"), direction = "wide")

## cor(as.matrix(test.r[, 4:ncol(test.r)]), use = "pairwise.complete.obs")

## test.r <- test.r[, order(colnames(test.r))]
## colnames(test.r) <- c("n", "trial", "v", paste0("k", 1:9))

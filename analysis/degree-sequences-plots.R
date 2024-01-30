library(optparse)
optionlist <- list(
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files. Default is %default. "
    ),
    make_option(
        c("-w", "--use-weights"), action = "store_true", default = FALSE,
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
    OptionParser(option_list = optionlist), #args = c("-s"),#, "-w", "-g", "dolphin"),
    convert_hyphens_to_underscores = TRUE
)

library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)

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
outfile_ds <- switch(
    use_weights + 1,
    paste0("../img/degree-sequences-", network, "-", dynamics, ".pdf"),
    paste0("../img/degree-sequences-", network, "-", dynamics, "-weighted.pdf")
)
outfile_kld <- switch(
    use_weights + 1,
    paste0("../img/kld-fig-", network, "-", dynamics, ".pdf"),
    paste0("../img/kld-fig-", network, "-", dynamics, "-weighted.pdf")
)

load(infile)
save_plots <- saveplots # hack-ish. Fix?
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]

    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])

    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos

    text(xpos, ypos, label, ...)
}

                                        # Plotting 1
ht <- 10
wd <- 7
labelsize <- 2
ticksize <- 2
if(save_plots) {
    pdf(outfile_ds, height = ht, width = wd)
} else {
    dev.new(width = wd, height = ht)
}
par(mfrow = c(5, 1))#c(length(histdat), 1))
for(i in 1:5) {#seq_along(histdat)) {
    par(mar = c(4, 4, 2, 0)+0.5)
    if(i == 1) {
        barplot(histdat[[i]], names.arg = kposs, col = "gray60",
                cex.axis = ticksize, cex.names = ticksize, cex.lab = labelsize,
                xlab = "Degree", ylab = "Frequency") # /max(histdat[[i]])
    } else if(i %in% 2:4){
        barplot(histdat[[i]], names.arg = kposs, col = 1:nrow(histdat[[i]]), # /max(histdat[[i]])
                beside = FALSE, cex.axis = ticksize, cex.names = ticksize, cex.lab = labelsize,
                xlab = "Degree", ylab = "Frequency")
    } else if(i == 5) {
        lnNi <- floor(log(N)) + 1
        barplot(histdat[[lnNi]], names.arg = kposs, col = 1:nrow(histdat[[lnNi]]),
                beside = FALSE, cex.axis = ticksize, cex.names = ticksize, cex.lab = labelsize,
                xlab = "Degree", ylab = "Frequency")
    }
}
if(save_plots) dev.off()

                                        # Plotting 2
ht <- 7
wd <- 14
labelsize <- 2
ticksize <- 1.75
if(save_plots) {
    pdf(outfile_kld, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))
                                        # KLD
plot(
    1:maxn, klds, type = "p", col = 1, pch = 16, cex = 2, lwd = 2,
    xlab = "", ylab = "Kullback-Leibler divergence", yaxt = "n",
    ylim = range(c(klds, as.numeric(rklds))), log = "y", 
    cex.axis = ticksize, cex.lab = labelsize
)
eaxis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize, at.small = FALSE, sub10 = "10", las = 0) # 
title(xlab = "n", font.lab = 3, cex.lab = labelsize)
points(1:maxn, colMeans(rklds), col = 2, pch = 15, cex = 2)
segments(
    x0 = 1:maxn,
    y0 = apply(rklds, 2, function(x) quantile(x, probs = 0.025)),
    y1 = apply(rklds, 2, function(x) quantile(x, probs = 0.975)),
    lty = 1, lwd = 2, col = 2
)
##mtext("A", line = -2.2, adj = 0.02, cex = labelsize, font = 2)
placelabel(paste0("(", letters[1], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
legend(
    "topright", legend = c("Optimized", "Random"), bty = "n",
    col = c(1, 2), pch = c(16, 15), pt.cex = 2, cex = labelsize
)
                                        # Error
plot(
    1:length(dl), sapply(errs, mean), type = "p", col = 1, pch = 16, cex = 2, lwd = 2, log = "y",
    xlab = "", ylab = "Approximation error", yaxt = "n",
    ylim = range(c(unlist(errs), as.numeric(rerrs))), 
    cex.axis = ticksize, cex.lab = labelsize
)
title(xlab = "n", font.lab = 3, cex.lab = labelsize)
## axis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize)
eaxis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize, at.small = FALSE, sub10 = "10", las = 0) # 
segments(
    x0 = 1:maxn,
    y0 = sapply(errs, function(x) quantile(x, probs = 0.025)),
    y1 = sapply(errs, function(x) quantile(x, probs = 0.975)),
    lty = 1, lwd = 2, col = 1
)
points(1:maxn, colMeans(rerrs), col = 2, pch = 15, cex = 2)
segments(
    x0 = 1:maxn,
    y0 = apply(rerrs, 2, function(x) quantile(x, probs = 0.025)),
    y1 = apply(rerrs, 2, function(x) quantile(x, probs = 0.975)),
    lty = 1, lwd = 2, col = 2
)
##mtext("B", line = -2.2, adj = 0.02, cex = labelsize, font = 2)
placelabel(paste0("(", letters[2], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
##axis(4, at = pretty(range(sapply(errs, mean))))
##mtext("Error", side = 4, font = 1, line = 3)
##abline(h = 0, col = 10, lwd = .75)
if(save_plots) dev.off()

## if(use_weights) {
##                                         # additional analysis
##                                         # scratch
##     ## dl is 1:maxn (12)
##     ## want ks and ws from each
##     palette("Tableau 10")
##     labelsize <- 1.75
##     ticksize <- 1.75
##     ks <- lapply(seq(maxn), function(n) unlist(lapply(dl[[n]], `[[`, "ks")))
##     ws <- lapply(seq(maxn), function(n) unlist(lapply(dl[[n]], `[[`, "ws")))
##     colors <- lapply(seq(maxn), function(n) rep(seq(n), times = ntrials))
##     pdf(paste0("./img/ks-and-ws.pdf"), height = 9, width = 12)
##     par(mfrow = c(3, 4))
##     par(mar = c(4, 4, 1, 1) + .5)
##     xlim <- range(unlist(ks))
##     ylim <- range(unlist(ws))
##     for(n in seq_along(ks)) {
##         plot(
##             NULL, type = "o", xlim = xlim, ylim = ylim, xlab = "Degree", ylab = "Weight",
##             cex.axis = ticksize, cex.lab = labelsize
##         )
##         for(i in seq(ntrials)) {
##             points(dl[[n]][[i]]$ks, dl[[n]][[i]]$ws, col = 1, pch = 19, lwd = 1, cex = .75)
##             lines(dl[[n]][[i]]$ks, dl[[n]][[i]]$ws, col = 1, lwd = 1, lty = 1)
##         }
##         mtext(paste("n", n, sep = "="), line = -2, adj = 0.98, font = 1, cex = labelsize)
##     }
##     dev.off()   
## }

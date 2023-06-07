library(optparse)
optionlist <- list(
    make_option(
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is %default. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

net <- args$network
dyn <- args$dynamics
infile <- paste0("../data/transfer-learning-", net, "-", dyn, ".RData")
outfile <- paste0("../img/transfer-learning-", net, "-", dyn, ".pdf")

load(infile)

sapply(soln_errors, median)
lapply(comp_errors, quantile, probs = c(0.025, 0.5, 0.975))

palette("Tableau 10")
## three panels
ticksize <- 1.75
labelsize <- 1.75
ypos <- 1:4
ylim <- c(0.5, 4.5)
pdf(outfile, width = 14, height = 7)
panels <- layout(matrix(c(1, 1, 2, 3), ncol = 2))
##layout.show(panels)
                                        # Bifurcation diagram
plot(NULL, xlim = range(bp), ylim = range(c(gt, a1, a2)), xlab = "Bifurcation parameter", ylab = "x",
     cex.axis = ticksize, cex.lab = labelsize)
lines(bp, gt, lty = 1, lwd = 3, col = 1)
lines(bp, a1, lty = 1, lwd = 3, col = 2)
lines(bp, a2, lty = 1, lwd = 3, col = 3)
legend("bottomright", bty = "n", legend = c("(1, 3, 5)", "(1, 2, 5)", "(1, 4, 5)"), title = "r =",
       col = 1:3, lty = 1, lwd = 3, cex = 1.25)
                                        # low
low_errors <- sapply(solns$low, `[[`, "error")
plot(
    NULL, xlab = "Error", ylab = "", log = "x", cex.axis = ticksize, cex.lab = labelsize,
    xlim = range(c(low_errors, rand_errors$low)), ylim = ylim, axes = FALSE
)
points(
    low_errors, jitter(rep(4, length(low_errors)), amount = 0.1), col = 1, pch = 19
)
points(
    as.numeric(comp_errors$low), jitter(rep(3, length(as.numeric(comp_errors$low))), amount = 0.1),
    col = 2, pch = 19
)
points(
    rand_errors$low, jitter(rep(2, length(rand_errors$low)), amount = 0.1),
    col = 3, pch = 19
)
points(
    soln_errors$low, jitter(rep(1, length(soln_errors$low)), amount = 0.1),
    col = 4, pch = 19
)
box()
axis(1, cex.axis = ticksize)
legend(
    "topright", bty = "n", pch = 19, cex = 1.25, col = 1:4,
    legend = c("Optimized on r = (1, 2, 5)", "Random, fixed-degree", "Random",
               "Optimized on r = (1, 3, 5)")
)
                                        # high
high_errors <- sapply(solns$high, `[[`, "error")
plot(
    NULL, xlab = "Error", ylab = "", log = "x", cex.axis = ticksize, cex.lab = labelsize,
    xlim = range(c(high_errors, rand_errors$high)), ylim = ylim, axes = FALSE
)
points(
    high_errors, jitter(rep(4, length(high_errors)), amount = 0.1), col = 1, pch = 19
)
points(
    as.numeric(comp_errors$high), jitter(rep(3, length(as.numeric(comp_errors$high))), amount = 0.1),
    col = 2, pch = 19
)
points(
    rand_errors$high, jitter(rep(2, length(rand_errors$high)), amount = 0.1),
    col = 3, pch = 19
)
points(
    soln_errors$high, jitter(rep(1, length(soln_errors$high)), amount = 0.1),
    col = 4, pch = 19
)
box()
axis(1, cex.axis = ticksize)
legend(
    "topright", bty = "n", pch = 19, cex = 1.25, col = 1:4,
    legend = c("Optimized on r = (1, 4, 5)", "Random, fixed-degree", "Random",
               "Optimized on r = (1, 3, 5)")
)
dev.off()

## ## better to do another strip plot here?

## ## running the transfer-learning.R script for the dolphin network with ntrials = 100
## ## try to get the above into a set of histograms somehow

## plotthis <- function(cond, ensemble_col = 4, soln_col = 1) { # orig low high
##     palette("Tableau 10")
##     soln_e <- soln_errors[[cond]]
##     comp_e <- as.numeric(comp_errors[[cond]])

##     xlim <- range(c(soln_e, comp_e))

##     soln <- hist(soln_e, breaks = 20, plot = FALSE)
##     comp <- hist(comp_e, breaks = 20, plot = FALSE)

##     ylim <- range(c(soln$counts, comp$counts))

##     plot(comp, xlim = xlim, ylim = ylim, col = ensemble_col,# freq = FALSE,
##          xlab = "Error", ylab = "Probability mass")
##     plot(soln, col = soln_col, add = TRUE)
## }

## pdf(outfile, width = 21)
## par(mfrow = c(1, 3))

## labelsize <- 1.75
## ticksize <- 1.75
## plotthis("orig")
## mtext("A", line = -2, at = 0.01, cex = labelsize, font = 2)
## plotthis("low")
## mtext("B", line = -2, at = 0.01, cex = labelsize, font = 2)
## plotthis("high")
## mtext("C", line = -2, at = 0.01, cex = labelsize, font = 2)

## dev.off()

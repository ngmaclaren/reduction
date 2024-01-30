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

library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)

load(infile)
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}


sapply(soln_errors, median)
lapply(comp_errors, quantile, probs = c(0.025, 0.5, 0.975))

palette("Tableau 10")

                                        # for the three self-versions
                                        # Also use the dark color for the main plot
                                        # Use Tableau 10 [1] for the main one
colorvecs <- list(# dark to light
    ##low = c("#6A3D9A", "#9F73BF", "#CAB2D6"),
    low = c("#F28E2B", "#F5A85C", "#F5C390"),
    ##high = c("#FF7F00", "#FCA036", "#FDBF6F")
    high = c("#E15759", "#E98081", "#EEAEAF")
)

## three panels
ticksize <- 1.75
labelsize <- 1.75
ypos <- 1:4
ylim <- c(0.5, 4.5)
pdf(outfile, width = 14, height = 7)
panels <- layout(matrix(c(1, 1, 2, 3), ncol = 2))
##layout.show(panels)

legend_A <- switch(
    dyn,
    dw = c(
        TeX(r"($r_2 = 2$)", italic = TRUE),
        TeX(r"($r_2 = 3$)", italic = TRUE),
        TeX(r"($r_2 = 4$)", italic = TRUE)
    ),
    mutualistic = c(
        TeX(r"($K = 3$)", italic = TRUE),
        TeX(r"($K = 5$)", italic = TRUE),
        TeX(R"($K = 7$)", italic = TRUE)
    ),
    genereg = c(
        TeX(r"($f = 2, h = 2$)", italic = TRUE),
        TeX(r"($f = 1, h = 2$)", italic = TRUE),
        TeX(r"($f = 1, h = 4$)", italic = TRUE)
    ),
    SIS = c(
        TeX(r"($\mu = 0.5$)", italic = TRUE),
        TeX(r"($\mu = 1$)", italic = TRUE),
        TeX(r"($\mu = 2$)", italic = TRUE)
    )
)
    
legend_B <- switch(
    dyn,
    dw = c(
        TeX(r"(Optimized, $r_2 = 2$)"), "Degree-preserving random", "Completely random", TeX(r"(Optimized, $r_2 = 3$)")
    ),
    mutualistic = c(
        TeX(r"(Optimized, $K = 3$)"), "Degree-preserving random", "Completely random", TeX(r"(Optimized, $K = 5$)")
    ),
    genereg = c(
        TeX(r"(Optimized, $f = 2, h = 2$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $f = 1, h = 2$)")
    ),
    SIS = c(
        TeX(r"(Optimized, $\mu = 0.5$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $\mu = 1$)")
    )
)

legend_C <- switch(
    dyn,
    dw = c(
        TeX(r"(Optimized, $r_2 = 4$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $r_2 = 3$)")
    ),
    mutualistic = c(
        TeX(r"(Optimized, $K = 7$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $K = 5$)")
    ),
    genereg = c(
        TeX(r"(Optimized, $f = 1, h = 4$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $f = 1, h = 2$)")
    ),
    SIS = c(
        TeX(r"(Optimized, $\mu = 2$)"), "Degree-preserving random", "Completely random",
        TeX(r"(Optimized, $\mu = 1$)")
    )
)

                                        # Bifurcation diagram
par(mar = c(5, 5, 1, 1))
plot(
    NULL, xlim = range(bp), ylim = range(c(gt, a1, a2)),
    xlab = TeX(r"($D$)", italic = TRUE), ylab = TeX(r"($\bar{x}$)", italic = TRUE),
    cex.axis = ticksize, cex.lab = labelsize
)
lines(bp, gt, lty = 1, lwd = 8, col = 1)
lines(bp, a1, lty = 1, lwd = 8, col = colorvecs$low[1])
lines(bp, a2, lty = 1, lwd = 8, col = colorvecs$high[1])
legend("bottomright", bty = "n",#  title = "r =",
       ## legend = c("(1, 3, 5)", "(1, 2, 5)", "(1, 4, 5)"),
       legend = legend_A,
       col = c(colorvecs$low[1], 1, colorvecs$high[1]), lty = 1, lwd = 4, cex = 1.5)
##mtext("A", line = -2.1, adj = 0.01, font = 2, cex = labelsize)
placelabel(paste0("(", letters[1], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
                                        # low
par(mar = c(5, 0, 1, 2))
low_errors <- low$oe# sapply(solns$low, `[[`, "error")
plot(
    NULL, xlab = "Approximation error", ylab = "", log = "x", cex.axis = ticksize, cex.lab = labelsize,
    xlim = range(c(low_errors, rand_errors$low)), ylim = ylim, axes = FALSE
)
points(
    low_errors, jitter(rep(4, length(low_errors)), amount = 0.1), col = colorvecs$low[1], pch = 19
)
points(
    as.numeric(comp_errors$low), jitter(rep(3, length(as.numeric(comp_errors$low))), amount = 0.1),
    col = colorvecs$low[2], pch = 19
)
points(
    rand_errors$low, jitter(rep(2, length(rand_errors$low)), amount = 0.1),
    col = colorvecs$low[3], pch = 19
)
points(
    soln_errors$low, jitter(rep(1, length(soln_errors$low)), amount = 0.1),
    col = 1, pch = 19
)
box()
## axis(1, cex.axis = ticksize)
eaxis(1, axTicks(1)[c(TRUE, FALSE)], at.small = FALSE, cex.axis = ticksize)
legend(
    "topright", bty = "n", pch = 19, cex = 1.5, col = c(colorvecs$low, 1),
    legend = legend_B
)
##mtext("B", line = -2.1, adj = 0.01, font = 2, cex = labelsize)
placelabel(paste0("(", letters[2], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
                                        # high
high_errors <- high$oe# sapply(solns$high, `[[`, "error")
plot(
    NULL, xlab = "Approximation error", ylab = "", log = "x", cex.axis = ticksize, cex.lab = labelsize,
    xlim = range(c(high_errors, rand_errors$high)), ylim = ylim, axes = FALSE
)
points(
    high_errors, jitter(rep(4, length(high_errors)), amount = 0.1), col = colorvecs$high[1], pch = 19
)
points(
    as.numeric(comp_errors$high), jitter(rep(3, length(as.numeric(comp_errors$high))), amount = 0.1),
    col = colorvecs$high[2], pch = 19
)
points(
    rand_errors$high, jitter(rep(2, length(rand_errors$high)), amount = 0.1),
    col = colorvecs$high[3], pch = 19
)
points(
    soln_errors$high, jitter(rep(1, length(soln_errors$high)), amount = 0.1),
    col = 1, pch = 19
)
box()
##axis(1, cex.axis = ticksize)
eaxis(1, axTicks(1)[c(TRUE, FALSE)], at.small = FALSE, cex.axis = ticksize)
legend(
    "topright", bty = "n", pch = 19, cex = 1.5, col = c(colorvecs$high, 1),
    legend = legend_C
)
##mtext("C", line = -2.1, adj = 0.01, font = 2, cex = labelsize)
placelabel(paste0("(", letters[3], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
dev.off()

library(optparse)
optionlist <- list(
    make_option(
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is %default. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
    ),
    make_option(
        c("-n", "--ntrials"), type = "integer", default = 3,
        help = "The number of independent simulations on each network [default %default]. To be more efficient, set to an even multiple of the number of usable cores. In the code, this defaults to (total number of available cores) - 1. The default is set based on many personal computers, which have 4 CPUs."
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

## running the transfer-learning.R script for the dolphin network with ntrials = 100
## try to get the above into a set of histograms somehow

plotthis <- function(cond, ensemble_col = 4, soln_col = 1) { # orig low high
    palette("Tableau 10")
    soln_e <- soln_errors[[cond]]
    comp_e <- as.numeric(comp_errors[[cond]])

    xlim <- range(c(soln_e, comp_e))

    soln <- hist(soln_e, breaks = 20, plot = FALSE)
    comp <- hist(comp_e, breaks = 20, plot = FALSE)

    ylim <- range(c(soln$counts, comp$counts))

    plot(comp, xlim = xlim, ylim = ylim, col = ensemble_col,# freq = FALSE,
         xlab = "Error", ylab = "Probability mass")
    plot(soln, col = soln_col, add = TRUE)
}

pdf(outfile, width = 21)
par(mfrow = c(1, 3))

labelsize <- 1.75
ticksize <- 1.75
plotthis("orig")
mtext("A", line = -2, at = 0.01, cex = labelsize, font = 2)
plotthis("low")
mtext("B", line = -2, at = 0.01, cex = labelsize, font = 2)
plotthis("high")
mtext("C", line = -2, at = 0.01, cex = labelsize, font = 2)

dev.off()

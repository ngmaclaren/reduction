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
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "The network to use in generating analysis results. Ignored if only one option is available. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist), #args = c("-s", "-w"),
    convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/") # lib.loc only necessary on the cluster

                                        # Options
palette("Tableau 10")
use_weights <- args$use_weights
network <- args$network
saveplots <- args$save_plots # ensure set locally, rather than by what was used before
infile <- switch(
    use_weights + 1,
    "../data/dolphin-demo.RData",
    "../data/dolphin-demo-weighted.RData"
)
outfile <- switch(
    use_weights + 1,
    "../img/dolphin-demo.pdf", # false
    "../img/dolphin-demo-weighted.pdf" # true
)

load(infile)
save_plots <- saveplots # hack-ish. Fix?

                                        # colors
colors <- list(
    nodestates = adjustcolor("gray60", alpha.f = .5),
    systemstate = "gray30",
    approximation = 1,
    selectednodes = 5,
    DART = 3,
    GBB = 4,
    exact = 2
)
ht <- 10
wd <- 12
if(save_plots) {
    pdf(outfile, height = ht, width = wd)
} else dev.new(width = wd, height = ht)
par(mfrow = c(2, 2), mar = c(5, 5, 1, 1)) # c(1, length(ns))
labelsize <- 2
ticksize <- 2
for(i in seq_along(solns)) {
    with(doublewell_parms, {
        matplot(Ds, Y, type = "l", lty = 1, lwd = .5, col = colors$nodestates,
                xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
                cex.lab = labelsize, cex.axis = ticksize)
        lines(Ds, y, lty = 1, lwd = 6, col = colors$systemstate)
        if(ns[i] == 1) {
            lines(Ds, DART[, 1], lty = 1, lwd = 8, col = colors$DART)
            lines(Ds, GBB[, 1], lty = 1, lwd = 8, col = colors$GBB)
            lines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 8, col = colors$approximation)
        } else {
            matlines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 4, col = colors$selectednodes)
            lines(
                Ds,
                switch(
                    use_weights + 1,
                    rowMeans(Y[, solns[[i]]$vs]), # false
                    apply(Y[, solns[[i]]$vs], 1, weighted.mean, solns[[i]]$ws)# true
                ),
                lty = 1, lwd = 8, col = colors$approximation
            )
        }
        if(ns[i] == 1 || ns[i] == 2) {
            lines(
                Ds,
                switch(
                    use_weights + 1,
                    rowMeans(as.matrix(Y[, exact[[i]]$vs])), # false
                    apply(as.matrix(Y[, exact[[i]]$vs]), 1, weighted.mean, exact[[i]]$ws) # true
                ),
                lty = 2, lwd = 5, col = colors$exact
            )
        }   
    })
    mtext(LETTERS[i], line = -2.2, adj = 0.02, font = 2, cex = labelsize)
    if(i == 2) { # length(ns)
        legend(
            "bottomright", bty = "n", cex = 1.5, lwd = 3, lty = c(rep(1, 6), 2),
            col = unlist(colors),
            legend = c(
                "Node states",
                "System state",
                "Approximation",
                "Selected nodes",
                "DART",
                "GBB",
                "Exact solultion"
            )
        )
    }
}
if(save_plots) dev.off()

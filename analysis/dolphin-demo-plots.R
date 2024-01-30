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
    OptionParser(option_list = optionlist), # args = c("-s"),# "-w"),
    convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/") # lib.loc only necessary on the cluster

                                        # Options
palette("Tableau 10")
useweights <- args$use_weights
network <- args$network
saveplots <- args$save_plots # ensure set locally, rather than by what was used before
infile <- switch(
    useweights + 1,
    "../data/dolphin-demo.RData",
    "../data/dolphin-demo-weighted.RData"
)
outfile <- switch(
    useweights + 1,
    "../img/dolphin-demo.pdf", # false
    "../img/dolphin-demo-weighted.pdf" # true
)

load(infile)
save_plots <- saveplots # hack-ish. Fix?
use_weights <- useweights
placelabel <- function(label, x, y, ...) {
    xlim <- par("usr")[1:2]
    ylim <- par("usr")[3:4]
    xpos <- xlim[1] + x*(xlim[2] - xlim[1])
    ypos <- ylim[1] + y*(ylim[2] - ylim[1])
    if(par("xlog")) xpos <- 10^xpos
    if(par("ylog")) ypos <- 10^ypos
    text(xpos, ypos, label, ...)
}

if(use_weights) {
    noweights <- new.env()
    with(noweights, load(gsub("-weighted", "", infile)))
}

                                        # colors
colors <- list(
    nodestates = adjustcolor("gray60", alpha.f = .5),
    systemstate = "gray30",
    approximation = 1,
    selectednodes = 2,
    GBB = 3,
    DART = 4## ,
    ## exact = 2
)
legendtext <- switch(
    use_weights + 1,
    c(
        "Node states", "System state", "Approximation", "Selected nodes", "GBB", "DART"## ,
        ## "Exact solultion"
    ),
    c("Node states", "System state", "Approximation", "Selected nodes"##, "Exact solution"
      )
)
legendlines <- switch(use_weights + 1, c(rep(1, 6), 2), c(rep(1, 4), 2))
legendcolors <- switch(
    use_weights + 1,
    unlist(colors),
    unlist(colors[-which(names(colors) %in% c("GBB", "DART"))])
)
ht <- switch(use_weights + 1, 12, 8)
wd <- 9
if(save_plots) {
    pdf(outfile, height = ht, width = wd)
} else dev.new(width = wd, height = ht)
par(mfrow = switch(use_weights + 1, c(3, 2), c(2, 2)),
    ##mai = c(.66, .66, .1, .1)
    mar = c(5, 5, 1, 1)
    )
labelsize <- 2
ticksize <- 2
for(i in seq_along(solns)) {
    with(doublewell_parms, {
        matplot(Ds, Y, type = "l", lty = 1, lwd = .5, col = colors$nodestates,
                xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
                cex.lab = labelsize, cex.axis = ticksize)
        lines(Ds, y, lty = 1, lwd = 6, col = colors$systemstate)
        if(ns[i] == 1) {
            ## lines(Ds, DART[, 1], lty = 1, lwd = 8, col = colors$DART)
            ## lines(Ds, GBB[, 1], lty = 1, lwd = 8, col = colors$GBB)
            lines(Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 8, col = colors$approximation)
            if(use_weights) {
                with(noweights, {
                    lines(doublewell_parms$Ds, Y[, solns[[i]]$vs], lty = 1, lwd = 8,
                          col = adjustcolor(colors$approximation, alpha.f = .5))
                })
            }
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
            if(use_weights) {
                with(noweights, {
                    lines(doublewell_parms$Ds, rowMeans(Y[, solns[[i]]$vs]), lty = 1, lwd = 8,
                          col = adjustcolor(colors$approximation, alpha.f = .5))
                })
            }
        }
        ## if(ns[i] == 1 || ns[i] == 2) {
        ##     lines(
        ##         Ds,
        ##         switch(
        ##             use_weights + 1,
        ##             rowMeans(as.matrix(Y[, exact[[i]]$vs])), # false
        ##             apply(as.matrix(Y[, exact[[i]]$vs]), 1, weighted.mean, exact[[i]]$ws) # true
        ##         ),
        ##         lty = 2, lwd = 5, col = colors$exact
        ##     )
        ## }   
    })
    ## mtext(LETTERS[i], line = -2.2, adj = 0.02, font = 2, cex = labelsize)
    placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
    if(i == 1) { # length(ns)
        legend(
            "bottomright", bty = "n", cex = switch(use_weights + 1, 1.5, 1.3),
            lwd = 3, lty = legendlines, col = legendcolors, legend = legendtext
        )
    }
}
if(!use_weights) {
    with(doublewell_parms, {
                                        # GBB
        matplot(Ds, Y, type = "l", lty = 1, lwd = .5, col = colors$nodestates,
                xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
                cex.lab = labelsize, cex.axis = ticksize)
        lines(Ds, GBB.obs, lty = 1, lwd = 8, col = "gray30")
        lines(Ds, GBB[, 1], lty = 1, lwd = 8, col = colors$GBB)
        ##mtext(LETTERS[6], line = -2.2, adj = 0.02, font = 2, cex = labelsize)
        placelabel(paste0("(", letters[5], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
                                        # DART
        matplot(Ds, Y, type = "l", lty = 1, lwd = .5, col = colors$nodestates,
                xlab = TeX(r"(D)", italic = TRUE), ylab = TeX(r"($x^*$)", italic = TRUE),
                cex.lab = labelsize, cex.axis = ticksize)
        lines(Ds, DART.obs, lty = 1, lwd = 8, col = "gray30")
        lines(Ds, DART[, 1], lty = 1, lwd = 8, col = colors$DART)
        ##mtext(LETTERS[5], line = -2.2, adj = 0.02, font = 2, cex = labelsize)
        placelabel(paste0("(", letters[6], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
    })
}
if(save_plots) dev.off()


                                        # Calculate ε for D <= 0.6
Dsr <- doublewell_parms$Ds[which(doublewell_parms$Ds < 0.6)]
idx <- seq_along(Dsr)
GBBr <- GBB[idx, ]
GBB.obsr <- GBB.obs[idx]
DARTr <- DART[idx, ]
DART.obsr <- DART.obs[idx]
Yr <- Y[idx, ]
yr <- y[idx]

                                        # ε for the full range (0, 1)
calc_obj(Y[, solns[[1]]$vs], y)#, Y)
calc_obj(DART, DART.obs)#, Y)
calc_obj(GBB, GBB.obs)#, Y)

                                        # ε for the reduced range, D < 0.6
calc_obj(Yr[, solns[[1]]$vs], yr)#, Yr)
calc_obj(DARTr, DART.obsr)#, Yr)
calc_obj(GBBr, GBB.obsr)#, Yr)

                                        # ε for the full range, n ∈ {2, 3, 4}
sapply(2:4, function(i) calc_obj(rowMeans(Y[, solns[[i]]$vs]), y))
sapply(2:4, function(i) calc_obj(rowMeans(Yr[, solns[[i]]$vs]), yr))

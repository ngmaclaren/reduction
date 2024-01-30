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
    OptionParser(option_list = optionlist), #args = c("-s"),#, "-w"),
    convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)
                                        # Options
palette("Tableau 10")
use_weights <- args$use_weights
network <- args$network
dynamics <- args$dynamics
saveplots <- args$save_plots # ensure set locally, rather than by what was used before
infile <- switch(
    use_weights + 1,
    paste0("../data/compare-networks-", dynamics, ".RData"),
    paste0("../data/compare-networks-", dynamics, "-weighted.RData")
)
outfile <- switch(
    use_weights + 1,
    paste0("../img/compare-networks-", dynamics, ".pdf"), # false
    paste0("../img/compare-networks-", dynamics, "-weighted.pdf") # true
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

                                        # Plotting
plotit <- function(dl, panellabel = "", xlab = "", ylab = "") {
    re <- dl$re; fe <- dl$fe; oe <- dl$oe
    ht <- 2
    wd <- 7
    amnt <- 0.05
    alpha <- 1
    ptcex <- 0.75
    ##par(mar = rep(0, 4))
    ylim <- c(.5, 1.5)
    xlim <- range(c(re, fe, oe))
    ypos <- 1 + c(.15, 0, -.15)
    plot(
        NULL, #c(median(oe), median(fe), median(re)), ypos, type = "p", pch = 1, lwd = 3, cex = 3,
        xlim = xlim, ylim = ylim, log = "x", # c("dodgerblue", "goldenrod", "firebrick")
        axes = FALSE, xlab = xlab, ylab = ylab, las = 1
    )
    points(
        oe, jitter(rep(ypos[1], ntrials), amount = amnt),
        col = adjustcolor(1, alpha.f = alpha), pch = 19, cex = ptcex
    )
    points(
        fe, jitter(rep(ypos[2], ntrials), amount = amnt),
        col = adjustcolor(2, alpha.f = alpha), pch = 19, cex = ptcex
    )
    points(
        re, jitter(rep(ypos[3], ntrials), amount = amnt),
        col = adjustcolor(3, alpha.f = alpha), pch = 19, cex = ptcex
    )
}

ticksize <- 1.75
labelsize <- 1.75
ht <- 10
wd <- 10
if(save_plots) {
    pdf(outfile, height = ht, width = wd)
} else {
    dev.new(height = 12)
}
par(mfcol = c(length(nets)/2, 2))
par(mar = c(2, 2, .25, 2)+0.75)
for(i in seq_len(length(nets))) {
    plotit(dl[[i]])
    ##axis(1, cex.axis = ticksize)
    eaxis(1, axTicks(1)[c(TRUE, FALSE)], cex.axis = ticksize, at.small = FALSE, sub10 = "10")
    ##mtext(LETTERS[i], cex = labelsize, line = -1, adj = 0.02, font = 2)
    placelabel(paste0("(", letters[i], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
    if(i == 1) {
        legend(
            "topright", bty = "n", col = 1:3, pch = 19, pt.cex = 1.5, cex = 1.5,
            legend = c("Optimized", "Degree-preserving random", "Completely random")
        )
    }
}
if(save_plots) dev.off()

tbl <- data.frame(
    fe = sapply(dl, function(x) sum(x$fe <= max(x$oe))/length(x$fe)),
    re = sapply(dl, function(x) sum(x$re <= max(x$oe))/length(x$re)),
    row.names = nets
)
    
print(tbl)

sum(dl[[1]]$re < max(dl[[1]]$oe))
sum(dl[[1]]$fe < max(dl[[1]]$oe))

summary(dl[[1]]$fe)
summary(dl[[1]]$oe)

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
    OptionParser(option_list = optionlist), # args = c("-s", "-w", "-g", "dolphin"),
    convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")

                                        # Options
palette("Tableau 10")
use_weights <- args$use_weights
network <- args$network
dynamics <- args$dynamics
saveplots <- args$save_plots # ensure set locally, rather than by what was used before
infile <- switch(
    use_weights + 1,
    paste0("../data/compare-dynamics-", network, ".RData"),
    paste0("../data/compare-dynamics-", network, "-weighted.RData")
)
outfile <- switch(
    use_weights + 1,
    paste0("../img/compare-dynamics-", network, ".pdf"), # false
    paste0("../img/compare-dynamics-", network, "-weighted.pdf") # true
)

load(infile)
save_plots <- saveplots # hack-ish. Fix?

colorbar <- function(colors, rng, title = "") {
    ##crp <- colorRampPalette(colors)
    z <- matrix(1:100, nrow = 1)
    x <- 1
    y <- seq(rng[1], rng[2], length.out = 100)
    par(mar = c(3, 0, 6, 3.5) + .5)
    image(x, y, z, col = colors,#crp(100),
          axes = FALSE, xlab = "", ylab = "")
    axis(4, at = pretty(rng), labels = pretty(rng), cex.axis = ticksize)
    mtext(title, 4, line = 2.5, cex = labelsize)
    box()
}

save_plots <- TRUE # FALSE

##maxval <- max(abs(range(compvals)))
maxval <- max(compvals)
##zlim <- c(-maxval, maxval)
zlim <- c(0, 5)# maxval)

zvals <- as.data.frame(cbind(expand.grid(5:1, 1:5), as.numeric(compvals)))
colnames(zvals) <- c("opt", "eval", "ratio")

plotvals <- ifelse(compvals > 5, 5, compvals)

ht <- 7
wd <- 8
labelsize <- 1.75
ticksize <- 1.75
if(save_plots) {
    pdf(outfile, height = ht, width = wd)
} else {
    dev.new(height = ht, width = wd)
}
par(mar = c(0, 4, 4, 0)+0.5)
hclpal <- "Spectral"# "Geyser"
nf <- layout(matrix(c(1, 2), nrow = 1, ncol = 2), c(4.1, 0.9), 5)
tm <- function(m) t(m)[, nrow(m): 1]
image(
    x = 1:5, y = 1:5, z = tm(plotvals), # log10(R/comps)
    col = hcl.colors(100, hclpal, rev = TRUE), zlim = zlim,
    ## xlab = "Evaluated on",
    ylab = "Optimized on", 
    cex.axis = ticksize, cex.lab = labelsize, axes = FALSE
)
axis(3, at = 1:5, labels = 1:5, cex.axis = ticksize)
axis(2, at = 1:5, labels = 5:1, cex.axis = ticksize)
mtext("Evaluated on", line = 3, cex = labelsize)
text(
    zvals$eval, zvals$opt, labels = format(round(zvals$ratio, 2), nsmall = 2),
    adj = .5, col = "black", cex = labelsize, font = 2
)
## labels <- c(
##     "Double-well", "SIS", "Michaelis-\nMertens", "Mutualistic\nspecies", "Wilson-\nCowan"
## )
colorbar(
    hcl.colors(100, hclpal, rev = TRUE),
    zlim, # range(as.numeric(compvals)), # log10(R/comps)
    "Average relative error" 
)
if(save_plots) dev.off()

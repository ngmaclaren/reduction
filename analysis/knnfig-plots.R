library(optparse)
optionlist <- list(
    make_option(
        c("-g", "--network"), type = "character", default = "email",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'. To run analysis for networks other than %default, run the associated simulation code first."
    ),
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is %default. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
    ),
    make_option(
        c("-w", "--use-weights"), action = "store_true", default = FALSE,
        help = "Use the simulation results with optimized node weights when generating analysis results. Default is %default. "
    ),
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files. Default is %default. "
    ),
    make_option(
        c("-r", "--random-seed"), type = "character", default = "none",
        help = "If used, the random seed used to build the file. Shown in the .Rdata file name."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist), # args = c("-s"),
    convert_hyphens_to_underscores = TRUE
)

library(igraph)
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
palette("Tableau 10")

saveplots <- args$save_plots
use_weights <- args$use_weights
network <- args$network
dynamics <- args$dynamics
if(args$random_seed == "none") {
    seedname <- ""
} else {
    seed <- as.integer(args$random_seed)
    seedname <- paste0("-", seed)
}

inputfile <- paste0(
    "../data/knnfig-", network, "-", dynamics, seedname,
    switch(use_weights + 1, "", "-weighted"),
    ".RData"
)

## load("../data/knnfig-email-dw.RData")
load(inputfile)
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

outputfile <- paste0(
    "../img/knnfig-", network, "-", dynamics, seedname,
    switch(use_weights + 1, "", "-weighted"),
    ".pdf"
)

                                        # assign colors    
colorref <- data.frame(
    k = sort(unique(dat$k)),
    color = 1:length(unique(dat$k))
)

                                        # plot
if(save_plots) {
    ## pdf("../img/knn-fig-noweights.pdf", height = 7, width = 14)
    pdf(outputfile, height = 7, width = 14)
} else {
    dev.new(height = 7, width = 14)
}
ticksize <- 1.75
labelsize <- 1.75
par(mar = c(4, 4, 1, 1)+0.5, mfrow = c(1, 2))
plot(
    jitter(dat$k, amount = 0.2), dat$error, pch = 1, cex = 1.5, lwd = 3,
    ##ylim = round(range(dat$error), digits = 3),
    col = adjustcolor(colorref$color[match(dat$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize, xlab = "Degree", ylab = "Approximation error"
)
## points(as.numeric(k[ref$vs]), rep(ref$error, n),
##        pch = 3, cex = 3, lwd = 3, col = "gray30")
abline(h = ref$error, lwd = 3, lty = 2, col = "gray30")
vquants <- data.frame(
    v = as.numeric(V(g)),
    k = degree(g),
    bc = betweenness(g, directed = FALSE, normalized = TRUE),
    cc = closeness(g, mode = "all", normalized = TRUE),
    tr = transitivity(g, "local"),
    knn = knn(g)$knn
)
##mtext("A", line = -2, adj = 0.02, cex = labelsize, font = 2)
placelabel(paste0("(", letters[1], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
datr <- aggregate(error ~ v, data = dat, FUN = mean)
## df <- merge(vquants, datr, by = "v")
df <- merge(dat, vquants, by = c("v", "k"))
plot(
    error ~ knn, data = df, pch = 1, cex = 1.5, lwd = 3,
    ## ylim = round(range(dat$error), digits = 3),
    col = adjustcolor(colorref$color[match(df$k, colorref$k)], alpha.f = 1),
    cex.axis = ticksize, cex.lab = labelsize,
    xlab = "Avgerage nearest neighbor degree", ylab = "Approximation error"
)
abline(h = ref$error, lwd = 3, lty = 2, col = "gray30")
## points(
##     df$knn[df$v %in% ref$vs], df$error[df$v %in% ref$vs],
##     col = "gray30", pch = 3, cex = 3, lwd = 3
## )
##mtext("B", line = -2, adj = 0.02, cex = labelsize, font = 2)
placelabel(paste0("(", letters[2], ")"), 0.01, 0.99, adj = c(0, 1), cex = labelsize)
if(save_plots) dev.off()


##                                         # weighted
## load("../data/knnfig-email-dw-weighted.RData")
## save_plots <- saveplots
## colorref <- data.frame(
##     k = sort(unique(dat$k)),
##     color = 1:length(unique(dat$k))
## )

## if(save_plots) {
##     pdf("../img/knn-fig-withweights.pdf", height = 7, width = 14)
## } else {
##     dev.new(height = 7, width = 14)
## }
## ticksize <- 1.75
## labelsize <- 1.75
## par(mar = c(4, 4, 1, 1)+0.5, mfrow = c(1, 2))
## plot(
##     jitter(dat$k, amount = 0.2), dat$error, pch = 1, cex = 1.5, lwd = 3,
##     col = adjustcolor(colorref$color[match(dat$k, colorref$k)], alpha.f = 1),
##     cex.axis = ticksize, cex.lab = labelsize, xlab = "Degree", ylab = "Error"
## )
## points(as.numeric(k[ref$vs]), rep(ref$error, n),
##        pch = 3, cex = 3, lwd = 3, col = "gray30")
## mtext("A", line = -2, adj = 0.02, cex = labelsize, font = 2)
## vquants <- data.frame(
##     v = as.numeric(V(g)),
##     k = degree(g),
##     bc = betweenness(g, directed = FALSE, normalized = TRUE),
##     cc = closeness(g, mode = "all", normalized = TRUE),
##     tr = transitivity(g, "local"),
##     knn = knn(g)$knn
## )
## datr <- aggregate(error ~ v, data = dat, FUN = mean)
## df <- merge(vquants, datr, by = "v")
## plot(
##     error ~ knn, data = df, pch = 1, cex = 1.5, lwd = 3,
##     col = adjustcolor(colorref$color[match(df$k, colorref$k)], alpha.f = 1),
##     cex.axis = ticksize, cex.lab = labelsize,
##     xlab = "Avg. nearest neighbor degree", ylab = "Avg. error"
## )
## points(
##     df$knn[df$v %in% ref$vs], df$error[df$v %in% ref$vs],
##     col = "gray30", pch = 3, cex = 3, lwd = 3
## )
## mtext("B", line = -2, adj = 0.02, cex = labelsize, font = 2)
## if(save_plots) dev.off()

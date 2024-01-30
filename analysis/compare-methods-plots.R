if(interactive()) {
    setwd("/user/neilmacl/Documents/reduction/analysis/")
    dynamics <- "dw"
} else {
    library(optparse)
    optionlist <- list(
        make_option(
            c("-d", "--dynamics"), type = "character", default = "dw",
            help = "The dynamics to simulate on each network. Default is 'dw'. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
        )
    )
    args <- parse_args(
        OptionParser(option_list = optionlist),
        convert_hyphens_to_underscores = TRUE
    )
    dynamics <- args$dynamics
}
                                        # Libraries and functions
library(parallel)
ncores <- detectCores()-1
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)

                                        # Options
palette("Tableau 10")

infile <- paste0("../data/compare-methods-", dynamics, ".RData")
outfile <- paste0("../img/compare-methods-", dynamics, "test.pdf")
nets <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)
envs <- paste(nets, "env", sep = "_")

load(infile)

pdf(outfile, height = 15, width = 10)
par(mar = c(4, 10, 1, 1), mfcol = c(5, 2))
labels <- c("Optimized", "Random", "Fixed-degree", "Quantiles", "Testing") # , "Discarding"
##for(env in lapply(envs, get)) {
for(i in seq_along(nets)) {
    env <- get(envs[i])
    net <- nets[i]
    with(env, {
        source("../src/functions.R", local = TRUE)
        testing <- make_dataset( # to update to new definition
            n, ntrials, bparam, y, Y, use_quantiles = TRUE
        )
        errors <- lapply(
            list(opts, rands, fixing, respecting, testing), # , discarding
            function(dl) sapply(dl, `[[`, "error")
        )
        colors <- seq_along(errors)
        names(errors) <- c("opts", "rands", "fixing", "respecting", "testing") # "discarding", 
        xlim <- range(unlist(errors))
        ylim <- range(colors) + c(-0.5, 0.5)
        plot(
            NULL, xlim = xlim, ylim = rev(ylim), log = "x",
            xlab = expression(epsilon), ylab = "", axes = FALSE, cex.lab = 1.75
        )
        box()
        axis(1, cex.axis = 1.75)
        axis(2, at = colors, labels = labels, las = 2, cex.axis = 1.75)
        for(i in seq_along(errors)) {
            points(
                errors[[i]], jitter(rep(colors[i], ntrials), amount = 0.25), col = colors[i], pch = 16
            )
        }
        mtext(net, adj = 0, line = -1)
        ##print(lapply(errors, summary))
    })
}
dev.off()

## ## do the t.tests
## for(env in lapply(envs, get)) {
##     with(env, {
##         print(
##             t.test(
##                 colMeans(sapply(rands, `[[`, "ks")),
##                 colMeans(sapply(opts, `[[`, "ks"))
##             )
##         )
##         print(
##             t.test(
##                 apply(sapply(rands, `[[`, "ks"), 2, sd),
##                 apply(sapply(opts, `[[`, "ks"), 2, sd)
##             )
##         )
##         print(
##             ##knn <- knn(g)$knn
##             t.test(
##                 apply(sapply(rands, `[[`, "vs"), 2, function(col) mean(knn[col])),
##                 apply(sapply(opts, `[[`, "vs"), 2, function(col) mean(knn[col]))
##             )
##         )
##     })
## }

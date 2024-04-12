library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

if(interactive()) {
    args$network <- "dolphin"
}

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(optNS)

network <- args$network

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
Y <- readRDS(paste0("../data/fullstate-", network, ".rds"))$doublewell
y <- rowMeans(Y)

ns <- 1:12

ntrials <- 100

opts <- lapply(ns, function(n) {
    make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y
    )
})

saveRDS(opts, file = paste0("../data/", network, "-1-to-12.rds"))

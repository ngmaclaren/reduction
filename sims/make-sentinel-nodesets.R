if(interactive()) setwd("/user/neilmacl/Documents/reduction/sims/")

                                        # Options
library(optparse)
optionlist <- list(
    make_option(
        "--nnodes", type = "integer", default = 0,
        help = "The number of nodes in the sentinel node set. Default is %default: if %default, program will reset to floor(log(N))."
    ),
    make_option(
        "--ntrials", type = "integer", default = 3,
        help = "The number of node sets to make. Will make ntrials optimized (sentinel) node sets and ntrials random node sets."
    ),
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to simulate dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        "--dynamics", type = "character", default = "dw",
        help = "The dynamics to simulate on the network. Default is %default. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist), convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
library(parallel)
ncores <- detectCores() - 1
library(igraph)
library(deSolve)
source("../src/functions.R")

                                        # Parameters
if(interactive()) {
    args$ntrials <- ncores
    args$network <- "proximity"
}

ntrials <- args$ntrials # the number of node sets to optimize. Default is 3

network <- args$network
load(paste0("../data/", network, ".rda")) # these should be .rds files. 
load(paste0("../data/fullstate-", network, ".rda")) # so should these.
g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)

if(args$nnodes == 0) { # the number of nodes in the node set. Default is floor(log(N))
    n <- floor(log(N))
} else {
    n <- args$nnodes 
}

dynamics <- args$dynamics
bparam <- switch(
    dynamics,
    dw = doublewell_parms$Ds,
    SIS = SIS_parms$Ds,
    genereg = genereg_parms$Ds,
    mutualistic = mutualistic_parms$Cs
)
Y <- fullstate[[dynamics]]
y <- rowMeans(Y)

outfile <- paste0("../data/optimized-nodesets/", paste(c(network, dynamics), collapse = "-"), ".rds") # .RData

                                        # Main call
opts <- make_dataset(
    n = n,
    ntrials = ntrials,
    bparam = bparam,
    y = y,
    Y = Y,
    optimize = TRUE
)

## rands <- make_dataset(
##     n = n, ntrials = ntrials, bparam = bparam, y = y, Y = Y
## )

## save.image(outfile) # <- change this to saveRDS(opts)
saveRDS(opts, outfile)

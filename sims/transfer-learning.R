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

library(parallel)
ncores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
library(igraph)

network <- args$network
dynamics <- args$dynamics
ntrials <- args$ntrials
outputfile <- paste0("../data/transfer-learning-", network, "-", dynamics, ".RData")

source("../src/functions.R")
load(paste0("../data/", network, ".rda"))
load(paste0("../data/fullstate-", network, ".rda"))
load(paste0("../data/fullstate-alt-", network, ".rda"))

g <- get(network)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
n <- floor(log(N))

bp <- switch(
    dynamics,
    dw = seq(0.001, 1, length.out = lout),
    SIS = seq(0.001, 1, length.out = lout),
    genereg = seq(0.001, 1, length.out = lout),
    mutualistic = seq(0.1, 2, length.out = lout)
)
alt1 <- switch(
    dynamics,
    dw = "dw1",
    SIS = "sis1",
    genereg = "gr1",
    mutualistic = "ms1"
)
alt2 <- switch(
    dynamics,
    dw = "dw2",
    SIS = "sis2",
    genereg = "gr2",
    mutualistic = "ms2"
)

GT <- fullstate[[dynamics]] # full ground truth
gt <- rowMeans(fullstate[[dynamics]]) # mean state

A1 <- fullstate_alt[[alt1]]
a1 <- rowMeans(A1)

A2 <- fullstate_alt[[alt2]]
a2 <- rowMeans(A2)

solns <- list(
    orig = mclapply(seq(ntrials), function(x) experiment(n, gt, GT, bp), mc.cores = ncores),
    low = mclapply(seq(ntrials), function(x) experiment(n, a1, A1, bp), mc.cores = ncores),
    high = mclapply(seq(ntrials), function(x) experiment(n, a2, A2, bp), mc.cores = ncores)
)
        
comps <- list(
    orig = mclapply(
        solns$orig,
        function(opt) make_dataset(n, ntrials, bp, gt, GT, comps = opt$vs),
        mc.cores = ncores
    ),
    low = mclapply(
        solns$low,
        function(opt) make_dataset(n, ntrials, bp, a1, A1, comps = opt$vs),
        mc.cores = ncores
    ),
    high = mclapply(
        solns$high,
        function(opt) make_dataset(n, ntrials, bp, a2, A2, comps = opt$vs),
        mc.cores = ncores
    )
)

soln_errors <- list(
    orig = sapply(solns$orig, `[[`, "error"),
    low = sapply(solns$orig, function(soln) obj_fn(soln$vs, a1, A1, bp)),
    high = sapply(solns$orig, function(soln) obj_fn(soln$vs, a2, A2, bp))
)

comp_errors <- lapply(comps, function(complist) {
    do.call(cbind, lapply(complist, function(comp) sapply(comp, `[[`, "error")))
})

save.image(outputfile)

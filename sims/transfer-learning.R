## Switch these to use generate_data() instead.
## Check transfer-learning-plots.R for the needed variables
## solns, soln_errors, comp_errors, rand_errors
## [x (I think)] For this sim and compare-networks.R, use the min error soln instead of the median.


library(optparse)
optionlist <- list(
    make_option(
        c("-v", "--verbose"), action = "store_true", default = FALSE,
        help = "Show optimization details. Setting this flag sets `trace = TRUE` in `optim()` and `mc.silent = FALSE` in `mclapply()`."
    ),
    make_option(
        c("-w", "--optimize-weights"), action = "store_true", default = FALSE,
        help = "Use quadratic programming to optimize weights of nodes selected by simulated annealing. Substantially increases run time (approximately 4--5x)."
    ),
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

verbose <- args$verbose
optimize_weights <- args$optimize_weights
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

                                        # Change this to take the sequences from the standard objects
bp <- switch(
    dynamics,
    dw = doublewell_parms$Ds, # seq(0.001, 1, length.out = lout),
    SIS = SIS_parms$Ds, # seq(0.001, 1, length.out = lout),
    genereg = genereg_parms$Ds, # seq(0.001, 1, length.out = lout),
    mutualistic = mutualistic_parms$Ds # seq(0.1, 2, length.out = lout)
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

orig <- generate_data(
    g, GT, gt, n, ntrials = ntrials, optimize_weights = optimize_weights, verbose = verbose
)
low <- generate_data(
    g, A1, a1, n, ntrials = ntrials, optimize_weights = optimize_weights, verbose = verbose
)
high <- generate_data(
    g, A2, a2, n, ntrials = ntrials, optimize_weights = optimize_weights, verbose = verbose
)

soln_errors <- list(
    orig = orig$oe, # sapply(solns$orig, `[[`, "error"),
    ## low = sapply(solns$orig, function(soln) obj_fn(soln$vs, a1, A1, bp)),
    low = apply(orig$ovs, 1, function(vs) obj_fn(vs, a1, A1, bp)),
    ## high = sapply(solns$orig, function(soln) obj_fn(soln$vs, a2, A2, bp))
    high = apply(orig$ovs, 1, function(vs) obj_fn(vs, a2, A2, bp))
)

comp_errors <- list(
    orig = orig$fe,
    low = low$fe,
    high = high$fe
)

## comp_errors <- lapply(comps, function(complist) {
##     do.call(cbind, lapply(complist, function(comp) sapply(comp, `[[`, "error")))
## })

rand_errors <- list(
    orig = orig$re,
    low = low$re,
    high = high$re
)

## rand_errors <- lapply(rand, function(comp) sapply(comp, `[[`, "error"))

save.image(outputfile)

## old code ##
## solns <- list(
##     orig = mclapply(seq(ntrials), function(x) experiment(n, gt, GT, bp), mc.cores = ncores),
##     low = mclapply(seq(ntrials), function(x) experiment(n, a1, A1, bp), mc.cores = ncores),
##     high = mclapply(seq(ntrials), function(x) experiment(n, a2, A2, bp), mc.cores = ncores)
## )
## comps <- list(
##     orig = mclapply(
##         solns$orig,
##         function(opt) make_dataset(n, ntrials, bp, gt, GT, comps = opt$vs),
##         mc.cores = ncores
##     ),
##     low = mclapply(
##         solns$low,
##         function(opt) make_dataset(n, ntrials, bp, a1, A1, comps = opt$vs),
##         mc.cores = ncores
##     ),
##     high = mclapply(
##         solns$high,
##         function(opt) make_dataset(n, ntrials, bp, a2, A2, comps = opt$vs),
##         mc.cores = ncores
##     )
## )

## rand <- list(
##     orig = make_dataset(n, ntrials^2, bp, gt, GT),
##     low = make_dataset(n, ntrials^2, bp, a1, A1),
##     high = make_dataset(n, ntrials^2, bp, a2, A2)
## )

## bad code ##
## solns <- lapply(list(orig, low, high), function(dl) {
##     idx <- which.min(dl$oe)
##     list(
##         vs = dl$ovs[[idx]],
##         error = dl$oe[[idx]],
##         ks = dl$oks[[idx]]
##     )
## })

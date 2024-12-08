## This is a file for making all the different kinds of sentinel nodesets I could want.
## The output file format is ./[network]-[dynamics].rds
## The output object is a list, the elements of which are "ns.type"
## Each element is a list, the output of make_dataset(), with elements "vs", "error", "ks", and "ws", with "ws" potentially NULL

## optparse
library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        "--dynamics", type = "character", default = "doublewell",
        help = "The dynamics on the network. Default is %default. Option: 'doublewell', 'SIS', 'genereg', 'mutualistic'."
    ),
    make_option(
        "--ntrials", type = "integer", default = 3,
        help = "The number of independent node sets of each type. Default is %default."
    ),
    make_option(
        "--optimize-weights", type = "character", default = "false",
        help = "Whether or not to optimize the node weights. Default is %default, set to 'true' or 'false'."
    ),
    make_option(
        "--whichhalf", type = "character", default = "first",
        help = "Default is %default. Use either first or last. Only for separating file names---has no other effect."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

if(interactive()) {
    args$network <- "prosper"
    args$dynamics <- "doublewell"
    args$ntrials <- 3
    args$optimize_weights <- FALSE # TRUE
}

## packages
library(parallel)
ncores <- detectCores()-1
library(igraph)
library(ROI)
library(ROI.plugin.qpoases)
library(optNS)

## organize variables
network <- args$network
dynamics <- args$dynamics
ntrials <- as.numeric(args$ntrials)
optimize_weights <- as.logical(toupper(args$optimize_weights))
whichhalf <- args$whichhalf

if(optimize_weights) {
    cond <- paste(c(network, dynamics, "w"), collapse = "_")
} else {
    cond <- paste(c(network, dynamics), collapse = "_")
}
print(cond)
fullstatefile <- paste0("../data/fullstate-", network, ".rds")
outfile <- paste0("../data/ns-", cond, "-", whichhalf, ".rds")

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
AL <- as_adj_list(g, "all")

Y <- readRDS(fullstatefile)[[dynamics]]
y <- rowMeans(Y)
n <- floor(log(N))

## main call
runtime <- system.time(
    opts <- make_dataset(
        ntrials = ntrials, ns.type = "opt", ncores = ncores,
        n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
    )
)


                                        # For timing the SA alg
timingDF <- data.frame(
    network = network, N = N, dynamics = dynamics, whichhalf = whichhalf, optweight = optimize_weights, runtime = as.numeric(runtime[3])
)

write.csv(timingDF, paste0("../shell/output/ns-times/", cond, "-", whichhalf, ".csv"), row.names = FALSE)

best <- opts[[which.min(get_error(opts))]]

fixeds <- make_dataset(
    ntrials = ntrials, ns.type = "fixed", ncores = ncores,
    n = n, g = g, comps = best$vs, y = y, Y = Y, optimize_weights = optimize_weights
)

rands <- make_dataset(
    ntrials = ntrials, ns.type = "rand", ncores = ncores,
    n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
)

constrs <- make_dataset(
    ntrials = ntrials, ns.type = "constr", ncores = ncores,
    n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
)

quants <- make_dataset(
    ntrials = ntrials, ns.type = "quant", ncores = ncores,
    n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
)

knnconstr <- make_dataset(
    ntrials = ntrials, ns.type = "knnconstr", ncores = ncores,
    n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
)

comms <- make_dataset(
    ntrials = ntrials, ns.type = "comm", ncores = ncores,
    n = n, g = g, y = y, Y = Y, optimize_weights = optimize_weights
)

allns <- list(
    opt = opts, fixed = fixeds, rand = rands, constr = constrs, quant = quants, knnconstr = knnconstr, comm = comms
)

## store
saveRDS(allns, file = outfile)

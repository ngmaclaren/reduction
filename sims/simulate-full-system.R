library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        "--ncparam", type = "numeric", default = 100,
        help = "The number of evenly spaced control parameter values to sample from the range. In this study, the only control parameter is D and the range is [0, 1] for doublewell, SIS, and genereg but [0, 3] for mutualistic."
    )
)
args <- parse_args(
    OptionParser(option_list = optionlist),
    convert_hyphens_to_underscores = TRUE
)

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(deSolve)
library(sdn)

if(interactive()) {
    args$network <- "celegans"
    args$ncparam <- ncores
}

network <- args$network
ncparam <- args$ncparam
print(network)
print(ncparam)
## fullstatefile <- paste0("../data/fullstate-", network, "-L", ncparam, ".rds")
fullstatefile <- paste0("../data/fullstate-", network, ".rds")

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
## if(is_weighted(g)) {
##     A <- as_adj(g, "both", attr = "weight", sparse = FALSE)
## } else {
##     A <- as_adj(g, "both", sparse = FALSE)
## }
AL <- as_adj_list(g, "all")

times <- 0:15 # for all? was 0:20 for SIS?

print("solving doublewell...")
Y.doublewell <- with(
    list(
        params = c(.doublewell, list(AL = AL)),
        control = list(times = times, ncores = ncores)
    ), {
        params$Ds <- seq(0, 1, length.out = ncparam)
        solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
    }
)

print("solving SIS...")
Y.SIS <- with(
    list(
        params = c(.SIS, list(AL = AL)),
        control = list(times = times, ncores = ncores)
    ), {
        params$Ds <- seq(0, 1, length.out = ncparam)
        solve_in_range(params$Ds, "D", SIS, rep(params$xinit.low, N), params, control, "ode")
    }
)

print("solving gene regulatory...")
Y.genereg <- with(
    list(
        params = c(.genereg, list(AL = AL)),
        control = list(times = times, ncores = ncores)
    ), {
        params$Ds <- seq(0, 1, length.out = ncparam)
        solve_in_range(params$Ds, "D", genereg, rep(params$xinit.high, N), params, control, "ode")
    }
)

print("solving mutualistic species...")
Y.mutualistic <- with(
    list(
        params = c(.mutualistic, list(AL = AL)),
        control = list(times = times, ncores = ncores)
    ),  {
        params$Ds <- seq(0, 3, length.out = ncparam)
        solve_in_range(params$Ds, "D", mutualistic, rep(params$xinit.low, N), params, control, "ode")
    }
)

fullstate <- list(
    doublewell = Y.doublewell, # a change!
    SIS = Y.SIS,
    genereg = Y.genereg,
    mutualistic = Y.mutualistic
)

attr(fullstate, "graph") <- g

objs <- sum(sapply(ls(), function(nomen) object.size(get(nomen))))
print(paste(round(objs/1024^2, 2), "Mb"))
saveRDS(fullstate, file = fullstatefile)

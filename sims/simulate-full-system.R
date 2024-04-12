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
    args$network <- ""
}

library(parallel)
ncores <- detectCores()-1
library(igraph)
library(deSolve)
library(localsolver)

network <- args$network
print(network)
fullstatefile <- paste0("../data/fullstate-", network, ".rds")

g <- readRDS(paste0("../data/", network, ".rds"))
N <- vcount(g)
if(is_weighted(g)) {
    A <- as_adj(g, "both", attr = "weight", sparse = FALSE)
} else {
    A <- as_adj(g, "both", sparse = FALSE)
}

times <- 0:15 # for all? was 0:20 for SIS?

print("solving doublewell...")
Y.doublewell <- with(
    list(
        params = c(.doublewell, list(A = A)),
        control = list(times = times, ncores = ncores)
    ), {
        solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
    }
)

print("solving SIS...")
Y.SIS <- with(
    list(
        params = c(.SIS, list(A = A)),
        control = list(times = times, ncores = ncores)
    ), {
        solve_in_range(params$Ds, "D", SIS, rep(params$xinit.low, N), params, control, "ode")
    }
)

print("solving gene regulatory...")
Y.genereg <- with(
    list(
        params = c(.genereg, list(A = A)),
        control = list(times = times, ncores = ncores)
    ), {
        solve_in_range(params$Ds, "D", genereg, rep(params$xinit.high, N), params, control, "ode")
    }
)

print("solving mutualistic species...")
Y.mutualistic <- with(
    list(
        params = c(.mutualistic, list(A = A)),
        control = list(times = times, ncores = ncores)
    ),  {
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
saveRDS(fullstate, file = fullstatefile)

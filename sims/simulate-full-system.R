library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options include 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', and 'lfr', among others."
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
    args$network <- "dolphin"
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
AL <- as_adj_list(g, "all") # this method will not work directly with weighted networks

times <- 0:15

                                        # Record total time to stdout, for comp time analyis
t_start <- Sys.time()

print("solving doublewell...")
system.time(
    Y.doublewell <- with(
        list(
            params = c(.doublewell, list(AL = AL)),
            control = list(times = times, ncores = ncores)
        ), {
            params$Ds <- seq(0, 1, length.out = ncparam)
            solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode")
        }
    )
)

print("solving SIS...")
system.time(
    Y.SIS <- with(
        list(
            params = c(.SIS, list(AL = AL)),
            control = list(times = times, ncores = ncores)
        ), {
            params$Ds <- seq(0, 1, length.out = ncparam)
            params$xinit.low <- 0.01 # changed in sdn to 0.001; for this study we use 0.01
            solve_in_range(params$Ds, "D", SIS, rep(params$xinit.low, N), params, control, "ode")
        }
    )
)

print("solving gene regulatory...")
system.time(
    Y.genereg <- with(
        list(
            params = c(.genereg, list(AL = AL)),
            control = list(times = times, ncores = ncores)
        ), {
            params$Ds <- seq(0, 1, length.out = ncparam)
            solve_in_range(params$Ds, "D", genereg, rep(params$xinit.high, N), params, control, "ode")
        }
    )
)

print("solving mutualistic species...")
system.time(
    Y.mutualistic <- with(
        list(
            params = c(.mutualistic, list(AL = AL)),
            control = list(times = times, ncores = ncores)
        ),  {
            params$Ds <- seq(0, 3, length.out = ncparam)
            solve_in_range(params$Ds, "D", mutualistic, rep(params$xinit.low, N), params, control, "ode")
        }
    )
)

t_stop <- Sys.time()

fullstate <- list(
    doublewell = Y.doublewell,
    SIS = Y.SIS,
    genereg = Y.genereg,
    mutualistic = Y.mutualistic
)

attr(fullstate, "graph") <- g

objs <- sum(sapply(ls(), function(nomen) object.size(get(nomen))))

print("Memory used")
print(paste(round(objs/1024^2, 2), "Mb"))
print("")
print("Total time elapsed")
print(t_stop - t_start)

saveRDS(fullstate, file = fullstatefile)

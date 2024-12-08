library(optparse)
optionlist <- list(
    make_option(
        "--network", type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options include 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', and 'lfr', among others."
    ),
    make_option(
        "--ncparam", type = "numeric", default = 100,
        help = "The number of evenly spaced control parameter values to sample from the range. Default is %default. In this study, the only control parameter is D and the range is [0, 1] for doublewell, SIS, and genereg but [0, 3] for mutualistic."
    ),
    make_option(
        "--splitsims", type = "logical", default = FALSE,
        help = "Use args to specify which model to simulate. Default is %default."
    ),
    make_option(
        "--dynamics", type = "character", default = "doublewell",
        help = "The dynamics on the network. Default is %default. Option: 'doublewell', 'SIS', 'genereg', 'mutualistic'."
    ),
    make_option(
        "--split50", type = "logical", default = FALSE,
        help = "If splitting sims, also split the ncparam sims between two nodes? Only implemented for mutualistic species: all other simulations are finishing within the 72 hr limit. Default is %default."
    ),
    make_option(
        "--which50", type = "character", default = "first",
        help = "Default is %default, which simulates on indices 1:50; other option is 'last', which simulates on indices 51:100. Options '--splitsims' and '--split50' must also be used or the option has no effect."
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
    args$dynamics <- "mutualistic"
    args$splitsims <- TRUE # FALSE
    args$split50 <- TRUE
    args$which50 <- "last"
    args$ncparam <- ncores
}

network <- args$network
ncparam <- args$ncparam
splitsims <- args$splitsims
split50 <- args$split50
which50 <- args$which50
dynamics <- args$dynamics
print(network)
print(ncparam)
## fullstatefile <- paste0("../data/fullstate-", network, "-L", ncparam, ".rds")
if(isFALSE(splitsims)) {
    fullstatefile <- paste0("../data/fullstate-", network, ".rds")
} else {
                                        # Will need to then combine later to support follow-on workflow
    fullstatefile <- paste0("../data/fullstate-", network, "-", dynamics, ".rds") 
}

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

if(isFALSE(splitsims) | dynamics == "doublewell") {
    print("solving doublewell...")
    system.time(
        Y.doublewell <- with(
            list(
                params = c(.doublewell, list(AL = AL)),
                control = list(times = times, ncores = ncores)
            ), {
                params$Ds <- seq(0, 1, length.out = ncparam)
                solve_in_range(params$Ds, "D", doublewell, rep(params$xinit.low, N), params, control, "ode", method = "adams", maxsteps = 20000)
            }
        )
    )
}

if(isFALSE(splitsims) | dynamics == "SIS") {
    print("solving SIS...")
    system.time(
        Y.SIS <- with(
            list(
                params = c(.SIS, list(AL = AL)),
                control = list(times = times, ncores = ncores)
            ), {
                params$Ds <- seq(0, 1, length.out = ncparam)
                params$xinit.low <- 0.01 # changed in sdn to 0.001; for this study we use 0.01
                solve_in_range(params$Ds, "D", SIS, rep(params$xinit.low, N), params, control, "ode", method = "adams", maxsteps = 20000)
            }
        )
    )
}

if(isFALSE(splitsims) | dynamics == "genereg") {
    print("solving gene regulatory...")
    system.time(
        Y.genereg <- with(
            list(
                params = c(.genereg, list(AL = AL)),
                control = list(times = times, ncores = ncores)
            ), {
                params$Ds <- seq(0, 1, length.out = ncparam)
                solve_in_range(params$Ds, "D", genereg, rep(params$xinit.high, N), params, control, "ode", method = "adams", maxsteps = 20000)
            }
        )
    )
}

if(isFALSE(splitsims) | dynamics == "mutualistic") {
    print("solving mutualistic species...")
    if(split50) {
        if(which50 == "first") {
            newstart <- 1
            newstop <- floor(ncparam/2)
        } else if(which50 == "last") {
            newstart <- floor(ncparam/2) + 1
            newstop <- ncparam
        }

        system.time(
            Y.mutualistic <- with(
                list(
                    params = c(.mutualistic, list(AL = AL)),
                    control = list(times = times, ncores = ncores)
                ), {
                    params$Ds <- seq(0, 3, length.out = ncparam)[seq(newstart, newstop, by = 1)]
                    solve_in_range(params$Ds, "D", mutualistic, rep(params$xinit.low, N), params, control, "ode", method = "adams", maxsteps = 100000)
                }
            )
        )

        fullstatefile <- gsub(".rds", paste0("-", which50, ".rds"), fullstatefile)
    } else {
        system.time(
            Y.mutualistic <- with(
                list(
                    params = c(.mutualistic, list(AL = AL)),
                    control = list(times = times, ncores = ncores)
                ),  {
                    params$Ds <- seq(0, 3, length.out = ncparam)
                    solve_in_range(params$Ds, "D", mutualistic, rep(params$xinit.low, N), params, control, "ode", method = "adams", maxsteps = 100000)
                }
            )
        )
    }
}

t_stop <- Sys.time()

## fullstate <- list(
##     doublewell = Y.doublewell,
##     SIS = Y.SIS,
##     genereg = Y.genereg,
##     mutualistic = Y.mutualistic
## )
fullstate <- list()
if("Y.doublewell" %in% ls()) fullstate$doublewell <- Y.doublewell
if("Y.SIS" %in% ls()) fullstate$SIS <- Y.SIS
if("Y.genereg" %in% ls()) fullstate$genereg <- Y.genereg
if("Y.mutualistic" %in% ls()) fullstate$mutualistic <- Y.mutualistic

print("dims")
sapply(fullstate, dim)
print("range")
sapply(fullstate, function(Y) range(apply(Y, 1, max)))

attr(fullstate, "graph") <- g

objs <- sum(sapply(ls(), function(nomen) object.size(get(nomen))))

print("Memory used")
print(paste(round(objs/1024^2, 2), "Mb"))
print("")
print("Total time elapsed")
print(t_stop - t_start)

saveRDS(fullstate, file = fullstatefile)

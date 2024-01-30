library(optparse)
optionlist <- list(
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is 'dw'. Options: 'dw', 'SIS', 'genereg', and 'mutualistic'."
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
library(igraph, lib.loc = "/user/neilmacl/rlocal")
library(deSolve)
source("../src/functions.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(54321)

dynamics <- args$dynamics
ntrials <- args$ntrials
outputfile <- paste0("../data/compare-methods-", dynamics, ".RData")

nets <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "ba", "hk", "gkk", "lfr"
)
envs <- paste(nets, "env", sep = "_")

for(i in seq_along(nets)) {
    net <- nets[i]
    env <- envs[i]
    assign(env, new.env())
    with(get(env), {
        load(paste0("../data/", net, ".rda"), envir = get(env))
        load(paste0("../data/fullstate-", net, ".rda"), envir = get(env))
        source("../src/functions.R", local = TRUE)
        g <- get(net, envir = get(env))
        N <- vcount(g)
        A <- as_adj(g, "both", sparse = FALSE)
        k <- degree(g)
        knn <- knn(g)$knn
        Y <- fullstate[[dynamics]]#$dw
        y <- rowMeans(Y)
        n <- floor(log(N))
        bparam <- doublewell_parms$Ds

        ## optimized
        opts <- make_dataset(n, ntrials, bparam, y, Y, optimize = TRUE)
        
        ## random
        rands <- make_dataset(n, ntrials, bparam, y, Y)

        ## random, fixing degree
        idx <- which.min(sapply(opts, `[[`, "error"))
        ref <- opts[[idx]]
        fixing <- make_dataset(n, ntrials, bparam, y, Y, comps = ref$vs)

        ## random, discarding top and bottom 10%
        ##discarding <- make_dataset(n, ntrials, bparam, y, Y, use_connections = TRUE)

        ## random, discarding top and bottom 10%, respecting quantiles
        respecting <- make_dataset(
            n, ntrials, bparam, y, Y, use_quantiles = TRUE # , use_connections = TRUE
        )
    })
}

save.image(outputfile)

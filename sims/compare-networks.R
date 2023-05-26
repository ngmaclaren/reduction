                                        # Options
library(optparse)
optionlist <- list(
    make_option(
        c("-v", "--verbose"), action = "store_true", default = FALSE,
        help = "Show optimization details. Setting this flag sets `trace = TRUE` in `optim()` and `mc.silent = FALSE` in `mclapply()`."
    ),
    make_option(
        c("-s", "--save-plots"), action = "store_true", default = FALSE,
        help = "Store output figures as PDFs. If tabular output is created, saves to a text file. If FALSE [default], will display plot in an R graphics device (and silently fail to do so if graphics devices are not available.) Ignored for simulation files."
    ),
    make_option(
        c("-w", "--optimize-weights"), action = "store_true", default = FALSE,
        help = "Use quadratic programming to optimize weights of nodes selected by simulated annealing. Substantially increases run time (approximately 4--5x)."
    ),
    make_option(
        c("-d", "--dynamics"), type = "character", default = "dw",
        help = "The dynamics to simulate on each network. Default is 'dw'. Options: 'dw', 'SIS', 'genereg', 'mutualistic', and 'wilsoncowan'."
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

                                        # Libraries, functions
library(parallel)
ncores <- detectCores() - 1
library(igraph)
library(deSolve)
source("../src/functions.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

                                        # Extract options
verbose <- args$verbose
save_plots <- args$save_plots
optimize_weights <- args$optimize_weights
dynamics <- args$dynamics #commandArgs(trailingOnly = TRUE)[1]
ntrials <- args$ntrials#100
outputfile <- switch(
    optimize_weights + 1,
    paste0("../data/compare-networks-", dynamics, ".RData"), #FALSE
    paste0("../data/compare-networks-", dynamics, "-weighted.RData") #TRUE
)

                                        # Data setup
nets <- c(
    "dolphin", "celegans", "proximity", "euroroad", "email", "er", "gkk", "ba", "hk", "lfr"
)
envs <- paste(nets, "env", sep = "_")

for(i in seq_along(nets)) {
    net <- nets[i] # a string
    env <- envs[i] # a string
    assign(env, new.env())
    with(get(env), {
        load(paste0("../data/", net, ".rda"), envir = get(env))
        load(paste0("../data/fullstate-", net, ".rda"), envir = get(env))
        source("../src/functions.R", local = TRUE)
        g <- get(net, envir = get(env))
        N <- vcount(g)
        A <- as_adj(g, "both", sparse = FALSE)
        k <- degree(g)
        Y <- fullstate[[dynamics]]#$dw
        y <- rowMeans(Y)
        n <- floor(log(N))
    })
}

                                        # Simulations
dl <- lapply(envs, function(env) {
    with(get(env), {
        generate_data(
            g, Y, y, n, ntrials = ntrials, optimize_weights = optimize_weights, verbose = verbose)
    })
})

save.image(outputfile)

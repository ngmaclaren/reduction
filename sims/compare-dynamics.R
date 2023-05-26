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
        c("-g", "--network"), type = "character", default = "dolphin",
        help = "The network on which to compare each dynamics. Default is %default. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
    ),
    make_option(
        c("-n", "--ntrials"), type = "integer", default = 3,
        help = "The number of independent simulations on each network [default %default]. To be more efficient, set to an even multiple of the number of usable cores. In the code, this defaults to (total number of available cores) - 1. The default is set based on many personal computers, which have 4 CPUs."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist), # args = c("-w", "-n", "25", "-v"),
    convert_hyphens_to_underscores = TRUE
)

                                        # Libraries and functions
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
net <- args$network
ntrials <- args$ntrials # 100. But is 25 sufficient?
outputfile <- switch(
    optimize_weights + 1,
    paste0("../data/compare-dynamics-", net, ".RData"), #FALSE
    paste0("../data/compare-dynamics-", net, "-weighted.RData") #TRUE
)

                                        # Load data and set parameters
load(paste0("../data/", net, ".rda"))
load(paste0("../data/fullstate-", net, ".rda"))
g <- get(net)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
n <- floor(log(N))

bparms <- list(
    dw = doublewell_parms$Ds,
    SIS = SIS_parms$Ds,
    genereg = genereg_parms$Ds,
    mutualistic = mutualistic_parms$Cs,
    wilsoncowan = wilsoncowan_parms$c5s
)

                                        # Optimize a set of nodes for each network
optimized <- list(
    dw = mclapply(seq_len(ntrials), function(x) {
        experiment(
            n, rowMeans(fullstate$dw), fullstate$dw, doublewell_parms$Ds,
            optimize_weights = optimize_weights
        )
    }, mc.cores = ncores, mc.silent = isFALSE(verbose)),
    SIS = mclapply(seq_len(ntrials), function(x) {
        experiment(
            n, rowMeans(fullstate$SIS), fullstate$SIS, SIS_parms$Ds,
            optimize_weights = optimize_weights
        )
    }, mc.cores = ncores, mc.silent = isFALSE(verbose)),
    genereg = mclapply(seq_len(ntrials), function(x) {
        experiment(
            n, rowMeans(fullstate$genereg), fullstate$genereg, genereg_parms$Ds,
            optimize_weights = optimize_weights
        )
    }, mc.cores = ncores, mc.silent = isFALSE(verbose)),
    mutualistic = mclapply(seq_len(ntrials), function(x) {
        experiment(
            n, rowMeans(fullstate$mutualistic), fullstate$mutualistic, mutualistic_parms$Cs,
            optimize_weights = optimize_weights
        )
    }, mc.cores = ncores, mc.silent = isFALSE(verbose)),
    wilsoncowan = mclapply(seq_len(ntrials), function(x) {
        experiment(
            n, rowMeans(fullstate$wilsoncowan), fullstate$wilsoncowan, wilsoncowan_parms$c5s,
            optimize_weights = optimize_weights
        )
    }, mc.cores = ncores, mc.silent = isFALSE(verbose))
)

                                        # Generate ntrials comps for each of the ntrials optimized
                                        # node sets
comps <- vector("list", length(optimized))
for(i in seq_along(optimized)) {
    comps[[i]] <- lapply(optimized[[i]], function(opt) {
        make_dataset(
            n, ntrials, bparms[[i]], rowMeans(fullstate[[i]]), fullstate[[i]],
            comps = opt$vs, optimize_weights = optimize_weights, verbose = verbose
        )
    })
}

compmats <- lapply(comps, function(comp) {
    ## each column is 25 random fixed degree with potentialy unique degree re other columns
    do.call(cbind, lapply(comp, function(optset) sapply(optset, `[[`, "error")))
})
compmeans <- sapply(compmats, mean)
names(compmeans) <- names(optimized) # these are the average errors for the fixed degree random. They align with the columns of R. This means that every element of R[, 1] is divided by colmeans[1]
compmeanmat <- matrix(rep(compmeans, times = length(optimized)), nrow = length(optimized), byrow = TRUE)


                                        # Calculate the error for each combination
R <- matrix(NA, nrow = length(optimized), ncol = length(fullstate))
for(i in seq_along(optimized)) {
    for(j in seq_along(fullstate)) {
        R[i, j] <- mean(
            sapply(optimized[[i]], function(opt) {
                obj_fn(
                    opt$vs, rowMeans(fullstate[[j]]), fullstate[[j]], bparms[[j]],
                    optimize_weights = optimize_weights
                )
            })
        )
    }
}
colnames(R) <- rownames(R) <- names(optimized)

##comps <- matrix(rep(diag(R), times = length(fullstate)), ncol = length(fullstate), byrow = TRUE)
compvals <- R/compmeanmat

save.image(outputfile)

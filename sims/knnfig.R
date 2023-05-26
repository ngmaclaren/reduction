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
        c("-g", "--network"), type = "character", default = "email",
        help = "The network on which to compare each dynamics. Default is %default to have a larger sample of nodes at each degree. Options: 'dolphin', 'celegans', 'proximity', 'euroroad', 'email', 'er', 'gkk', 'ba', 'hk', 'lfr'."
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

                                        # Libraries and functions
library(parallel)
ncores <- detectCores() - 1
library(igraph)
library(deSolve)
source("../src/functions.R")

                                        # Extract options
verbose <- args$verbose
save_plots <- args$save_plots
optimize_weights <- args$optimize_weights
dynamics <- args$dynamics
net <- args$network
ntrials <- args$ntrials

load(paste0("../data/", net, ".rda"))
load(paste0("../data/fullstate-", net, ".rda"))
outputfile <- switch(
    optimize_weights + 1,
    paste0("../data/knnfig-", net, "-", dynamics, ".RData"), # FALSE
    paste0("../data/knnfig-", net, "-", dynamics, "-weighted.RData")
)

g <- get(net)
N <- vcount(g)
A <- as_adj(g, "both", sparse = FALSE)
k <- degree(g)
bparam <- switch(
    dynamics,
    dw = doublewell_parms$Ds,
    SIS = SIS_parms$Ds,
    genereg = genereg_parms$Ds,
    mutualistic = mutualistic_parms$Cs,
    wilsoncowan = wilsoncowan_parms$c5s
)
Y <- fullstate[[dynamics]]
y <- rowMeans(Y)
n <- floor(log(N))

ref <- experiment(n, y, Y, bparam, optimize_weights = optimize_weights, trace = verbose)

dat <- mclapply(
    ref$vs, function(i) {
        comps <- as.numeric(V(g)[which(k == k[i])])
        lapply(comps, function(j) {
            newset <- c(ref$vs[-which(ref$vs == i)], j)
            obj <- obj_fn(newset, y, Y, bparam, optimize_weights = optimize_weights)
            c(v = j, error = obj, k = as.numeric(k[j]))
        })
    }, mc.cores = ncores, mc.silent = isFALSE(verbose)
)
dat <- lapply(dat, function(x) do.call(rbind, x))
dat <- as.data.frame(do.call(rbind, dat))

save.image(outputfile)

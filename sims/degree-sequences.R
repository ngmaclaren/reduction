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

library(philentropy, lib.loc = "/user/neilmacl/rlocal/")
KLD <- function(x, y) {
    x <- x/sum(x) # should be refcounts
    y <- y/sum(y) # and the degree distribution to compare; need to remake the data
    KL(rbind(x, y))
}

                                        # Extract options
verbose <- args$verbose
save_plots <- args$save_plots
optimize_weights <- args$optimize_weights
dynamics <- args$dynamics
net <- args$network
ntrials <- args$ntrials
outputfile <- switch(
    optimize_weights + 1,
    paste0("../data/degree-sequences-", net, "-", dynamics, ".RData"),          # FALSE
    paste0("../data/degree-sequences-", net, "-", dynamics, "-weighted.RData")  # TRUE
)

load(paste0("../data/", net, ".rda"))
load(paste0("../data/fullstate-", net, ".rda"))
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

maxn <- 12
                                        # Generate `ntrials` optimized node sets for each `n` from 1:12
dl <- lapply(1:maxn, function(x) {
    make_dataset(x, ntrials, bparam, y, Y, optimize = TRUE, optimize_weights = optimize_weights,
                 verbose = verbose)
})

                                        # Tally k_i for the original distribution. Could have been
                                        # written `table(k)`.
kposs <- 1:max(k)
refcounts <- sapply(kposs, function(x) sum(k == x))

                                        # Tallies the number of nodes at each degree, summed over the
                                        # trials, for each `n`. For the dolphin network, the max degree
                                        # is 12. The final `histdat` is a list of matrices. The first
                                        # element is the original degree distribution. Elements 2:9 are
                                        # for the optimized node sets with n = 1:8. Within each
                                        # element, the first row is for the first node in S, the second
                                        # row is for the second node in S, and so on; the column
                                        # corresponds to the degree of the node; and the entry
                                        # corresponds to the count. The purpose of this data structure
                                        # is to support `barplot()` in a downstream file.
tally_row <- function(pos, n) {
    sapply(kposs, function(x) sum(sapply(dl[[n]], function(y) y$ks[pos]) == x))
}
histdat <- lapply(1:8, function(n) lapply(1:n, function(pos) tally_row(pos, n)))
for(i in seq_along(histdat)) {
    histdat[[i]] <- do.call(rbind, histdat[[i]])
}
histdat <- c(list(t(as.matrix(refcounts))), histdat)

                                        # Now, to show larger n, do the KLD figure
tally_all <- function(n, dl) {
    sapply(kposs, function(x) sum(sapply(dl[[n]], function(y) y$ks) == x))
}

klddat <- lapply(1:length(dl), tally_all, dl)
klds <- sapply(klddat, function(x) KLD(refcounts, x))
errs <- lapply(dl, function(ds) sapply(ds, `[[`, "error"))

                                        # Need to do the random version
make_rds <- function(maxn) {
    dl <- lapply(1:maxn, function(x) {
        make_dataset(x, ntrials, bparam, y, Y, optimize = FALSE,
                     optimize_weights = optimize_weights, verbose = verbose)
    })
    klddat <- lapply(1:length(dl), tally_all, dl)
    klds <- sapply(klddat, function(x) KLD(refcounts, x))
    errs <- lapply(dl, function(run) sapply(run, `[[`, "error"))
    list(klds = klds, errors = sapply(errs, mean))
}
rds <- mclapply(
    seq(ntrials),
    function(trial) make_rds(maxn),
    mc.cores = ncores, mc.silent = isFALSE(verbose)
)
rklds <- do.call(
    rbind,
    lapply(rds, `[[`, "klds")
)
rerrs <- do.call(
    rbind,
    lapply(rds, `[[`, "errors") # function(x) do.call(cbind, x$errors)
)

save.image(outputfile)
